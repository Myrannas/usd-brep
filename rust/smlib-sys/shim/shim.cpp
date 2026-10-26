// C shim over SMLib / SM_API for the smlib-sys Rust crate. See shim.h.

#include "shim.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <new>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <SmApiBrep.h>
#include <SmApiCurves.h>
#include <SmApiGeneral.h>
#include <SmApiImportExport.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmApiSurfaces.h>
#include <SmApiTrimmedSurfaces.h>

#include <SmAttribute.h>
#include <SmAxis2Placement.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmCircle.h>
#include <SmCone.h>
#include <SmCopyBrepMap.h>
#include <SmCurveClass.h>
#include <SmCylinder.h>
#include <SmDatabaseIO.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmEllipse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmFilletCorner.h>
#include <SmFilletExecutive.h>
#include <SmFilletSolver.h>
#include <SmFilletStandardSolver.h>
#include <SmGraphicsExtern.h>
#include <SmLine.h>
#include <SmLoop.h>
#include <SmLoopuse.h>
#include <SmMerge.h>
#include <SmOffsetExecutive.h>
#include <SmOffsetGeometryCreation.h>
#include <SmPlane.h>
#include <SmPoly.h>
#include <SmPrimitiveCreation.h>
#include <SmRotationalSweepGeometry.h>
#include <SmShell.h>
#include <SmSphere.h>
#include <SmSurfOfExtrusion.h>
#include <SmSurfOfRevolution.h>
#include <SmSurfTypes.h>
#include <SmTess.h>
#include <SmTopologyIntersector.h>
#include <SmTopologySweep.h>
#include <SmTorus.h>
#include <SmVertex.h>
#include <SmVertexuse.h>

struct sm_mesh {
    std::vector<float> positions;
    std::vector<float> normals;
    std::vector<float> surface_uvs;
    std::vector<uint32_t> indices;
    std::vector<uint32_t> face_ids;
    std::vector<sm_face_info> faces;
    std::vector<sm_edge_info> edges;
    std::vector<float> edge_points;
    size_t polygons_split = 0;
    size_t corners_missing_uv = 0;
};

namespace {

SmBrep* brep_of(sm_brep* brep) { return reinterpret_cast<SmBrep*>(brep); }
const SmBrep* brep_of(const sm_brep* brep) { return reinterpret_cast<const SmBrep*>(brep); }
sm_brep* handle_of(SmBrep* brep) { return reinterpret_cast<sm_brep*>(brep); }

// Runs f, turning C++ exceptions into status codes so none cross the C boundary.
template <typename F>
sm_status guarded(F&& f) {
    try {
        return f();
    } catch (const std::bad_alloc&) {
        return SM_SHIM_ERR_OUT_OF_MEMORY;
    } catch (...) {
        return SM_SHIM_ERR_EXCEPTION;
    }
}

void put3(double out[3], const SmVector3d& v) {
    out[0] = v.x;
    out[1] = v.y;
    out[2] = v.z;
}

// Labels ride on SMLib's attribute system (SmAttribute.h) under one
// application attribute id (applications use ids above 10000). The kernel
// copies a COPY attribute onto both pieces of a split face or edge, onto
// copies, and onto a coincident face from the other brep of a boolean.
const ULONG LABEL_ATTRIBUTE_ID = 0x534d4c42;  // "SMLB"

// A label as a generic attribute (longs: feature, index; characters: role),
// so it also reads back from files as one. Merging is the only behaviour it
// changes: when two labelled faces or edges merge, the survivor takes the
// older label (the lower feature id, then the lower index), whichever of the
// two the kernel keeps.
class LabelAttribute : public SmGenericAttribute {
public:
    LabelAttribute(const SmTArray<long>& longs, const SmTArray<char>& role)
        : SmGenericAttribute(LABEL_ATTRIBUTE_ID, SM_AB_COPY, longs, SmTArray<double>(), role) {}
    LabelAttribute(const LabelAttribute& other) : SmGenericAttribute(other) {}

    SmAttribute* MakeCopy(const SmContext& context) const override { return new (context) LabelAttribute(*this); }

    // SmAObject::Notify(SM_NO_MERGE_IN_BREP) calls this on the survivor's
    // attributes with the deleted object as the second argument, and only
    // copies the deleted object's attributes across when the survivor has
    // none with the same id.
    void Merge(SmAObject* survivor, SmAObject* deleted, SmAObject* result) override {
        const SmAttribute* other =
            deleted && deleted != survivor ? deleted->FindAttribute(LABEL_ATTRIBUTE_ID) : nullptr;
        if (other && result == survivor && survivor && survivor->FindAttribute(LABEL_ATTRIBUTE_ID) == this &&
            other->GetNumLongElements() >= 2 && m_vLongElements.GetSize() >= 2 &&
            std::make_pair(other->GetLongElementsAddress()[0], other->GetLongElementsAddress()[1]) <
                std::make_pair(m_vLongElements[0], m_vLongElements[1])) {
            m_vLongElements.RemoveAll();
            for (ULONG i = 0; i < other->GetNumLongElements(); ++i) {
                m_vLongElements.Add(other->GetLongElementsAddress()[i]);
            }
            m_vCharacterElements.RemoveAll();
            for (ULONG i = 0; i < other->GetNumCharacterElements(); ++i) {
                m_vCharacterElements.Add(other->GetCharacterElementsAddress()[i]);
            }
        }
        SmGenericAttribute::Merge(survivor, deleted, result);
    }
};

sm_label read_label(const SmAObject* object) {
    sm_label label;
    std::memset(&label, 0, sizeof label);
    const SmAttribute* attribute = object ? object->FindAttribute(LABEL_ATTRIBUTE_ID) : nullptr;
    if (!attribute || attribute->GetNumLongElements() < 2 || attribute->GetNumCharacterElements() == 0) return label;
    label.present = 1;
    label.feature = static_cast<uint32_t>(attribute->GetLongElementsAddress()[0]);
    label.index = static_cast<uint32_t>(attribute->GetLongElementsAddress()[1]);
    const size_t length = std::min<size_t>(attribute->GetNumCharacterElements(), SM_LABEL_ROLE_CAPACITY - 1);
    std::memcpy(label.role, attribute->GetCharacterElementsAddress(), length);
    return label;
}

// The role's length, or 0 when it's empty or not NUL-terminated in its buffer.
size_t role_length(const sm_label& label) {
    const void* end = std::memchr(label.role, '\0', SM_LABEL_ROLE_CAPACITY);
    return end ? static_cast<size_t>(static_cast<const char*>(end) - label.role) : 0;
}

void write_label(SmAObject& object, const sm_label& label) {
    if (SmAttribute* old = object.FindAttribute(LABEL_ATTRIBUTE_ID)) object.RemoveAttribute(old);
    if (!label.present) return;
    SmTArray<long> longs;
    longs.Add(static_cast<long>(label.feature));
    longs.Add(static_cast<long>(label.index));
    SmTArray<char> role;
    const size_t length = role_length(label);
    for (size_t i = 0; i < length; ++i) role.Add(label.role[i]);
    object.AddAttribute(new (*object.GetContext()) LabelAttribute(longs, role));
}

void describe_placement(sm_face_info& info, const SmAxis2Placement& placement) {
    put3(info.origin, placement.GetOrigin());
    put3(info.axis, placement.GetZAxis());
    put3(info.ref_dir, placement.GetXAxis());
}

sm_face_info describe_face(const SmFace* face) {
    sm_face_info info;
    std::memset(&info, 0, sizeof info);
    info.label = read_label(face);
    const SmExtent2d domain = face->GetUVDomain();
    info.uv_min[0] = domain.GetUMin();
    info.uv_min[1] = domain.GetVMin();
    info.uv_max[0] = domain.GetUMax();
    info.uv_max[1] = domain.GetVMax();

    const SmSurface* surface = face->GetSurface();
    if (!surface) {
        info.kind = SM_SURFACE_OTHER;
        return info;
    }
    const SM_TYPE type = surface->GetType();
    info.kernel_type = static_cast<int32_t>(type);
    // Exact type ids, most derived first: SmPlane derives from
    // SmBSplineSurface, and SmCylinder from SmCone.
    if (type == SmPlane_TYPE) {
        info.kind = SM_SURFACE_PLANE;
        describe_placement(info, static_cast<const SmPlane*>(surface)->GetPosition());
    } else if (type == SmCylinder_TYPE || type == SmCone_TYPE) {
        // SmApiCreateCylinder makes an SmCone with a zero semi-angle.
        const SmCone* cone = static_cast<const SmCone*>(surface);
        const bool cylinder = type == SmCylinder_TYPE || cone->IsCylinder();
        info.kind = cylinder ? SM_SURFACE_CYLINDER : SM_SURFACE_CONE;
        describe_placement(info, cone->GetPosition());
        info.radius = cone->GetRadius(0.0);
        // GetRadius(h) = base + tan(semi-angle) * h.
        info.half_angle_rad = cylinder ? 0.0 : std::atan(cone->GetRadius(1.0) - cone->GetRadius(0.0));
    } else if (type == SmSphere_TYPE) {
        const SmSphere* sphere = static_cast<const SmSphere*>(surface);
        info.kind = SM_SURFACE_SPHERE;
        describe_placement(info, sphere->GetPosition());
        info.radius = sphere->GetRadius();
    } else if (type == SmTorus_TYPE) {
        const SmTorus* torus = static_cast<const SmTorus*>(surface);
        info.kind = SM_SURFACE_TORUS;
        describe_placement(info, torus->GetPosition());
        info.radius = torus->GetMajorRadius();
        info.minor_radius = torus->GetMinorRadius();
    } else if (type == SmBSplineSurface_TYPE) {
        const SmBSplineSurface* nurbs = static_cast<const SmBSplineSurface*>(surface);
        info.kind = SM_SURFACE_NURBS;
        info.degree_u = static_cast<uint32_t>(nurbs->GetDegree(SM_SP_U));
        info.degree_v = static_cast<uint32_t>(nurbs->GetDegree(SM_SP_V));
        info.poles_u = static_cast<uint32_t>(nurbs->GetNumberControlPoints(SM_SP_U));
        info.poles_v = static_cast<uint32_t>(nurbs->GetNumberControlPoints(SM_SP_V));
        info.rational = nurbs->IsRational() ? 1 : 0;
    } else {
        info.kind = SM_SURFACE_OTHER;
    }
    return info;
}

struct PolyBrepDeleter {
    void operator()(SmPolyBrep* poly) const {
        poly->SetOKBackPtrs(FALSE);
        delete poly;
    }
};

using BooleanFn = SmApiStatus (*)(SmBrep*, SmBrep*, SmBrep*&);

sm_status run_boolean(BooleanFn op, sm_brep* a, sm_brep* b, sm_brep** out) {
    if (!a || !b || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* result = nullptr;
        const SmApiStatus status = op(brep_of(a), brep_of(b), result);
        if (status != SM_SUCCESS) return status;
        if (!result) return SM_ERR;
        *out = handle_of(result);
        return SM_SUCCESS;
    });
}

// One output vertex per (poly vertex, aux data, B-rep face): faces never share
// vertices, and the two sides of a seam get their own (u, v).
using VertexKey = std::tuple<const SmPolyVertex*, const SmPolyVertAuxData*, uint32_t>;

struct MeshBuilder {
    sm_mesh& mesh;
    std::map<VertexKey, uint32_t> vertices;

    uint32_t vertex(const SmPolyEdge* corner, uint32_t face_id) {
        const SmPolyVertex* pv = corner->GetStartPolyVertex();
        const SmPolyVertAuxData* aux = corner->GetAuxData();
        const VertexKey key(pv, aux, face_id);
        const auto found = vertices.find(key);
        if (found != vertices.end()) return found->second;

        const uint32_t index = static_cast<uint32_t>(mesh.positions.size() / 3);
        vertices.emplace(key, index);
        const SmPoint3d p = pv->GetPoint();
        mesh.positions.insert(mesh.positions.end(),
                              {static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z)});
        if (aux) {
            const SmVector3d& n = aux->m_vNormal;
            mesh.normals.insert(mesh.normals.end(),
                                {static_cast<float>(n.x), static_cast<float>(n.y), static_cast<float>(n.z)});
            mesh.surface_uvs.insert(mesh.surface_uvs.end(),
                                    {static_cast<float>(aux->m_vUV.x), static_cast<float>(aux->m_vUV.y)});
        } else {
            ++mesh.corners_missing_uv;
            mesh.normals.insert(mesh.normals.end(), {0.0f, 0.0f, 0.0f});
            mesh.surface_uvs.insert(mesh.surface_uvs.end(), {NAN, NAN});
        }
        return index;
    }
};

sm_status build_triangles(SmPolyBrep& poly, const std::unordered_map<const SmFace*, uint32_t>& face_index,
                          sm_mesh& mesh) {
    SmTArray<SmPolyFace*> polygons;
    poly.GetPolyFaces(polygons);
    MeshBuilder builder{mesh, {}};
    SmTArray<SmPolyLoop*> loops;
    SmTArray<SmPolyEdge*> corners;
    std::vector<uint32_t> ring;
    for (ULONG i = 0; i < polygons.GetSize(); ++i) {
        const SmPolyFace* polygon = polygons[i];
        const auto found = face_index.find(polygon->GetOriginalFace());
        if (found == face_index.end()) return SM_SHIM_ERR_UNMAPPED_FACE;
        const uint32_t face_id = found->second;

        loops.RemoveAll();
        polygon->GetPolyLoops(loops);
        if (loops.GetSize() != 1) return SM_SHIM_ERR_UNSUPPORTED_POLYGON;
        corners.RemoveAll();
        loops[0]->GetPolyEdges(corners);
        if (corners.GetSize() < 3) return SM_SHIM_ERR_UNSUPPORTED_POLYGON;

        ring.clear();
        for (ULONG j = 0; j < corners.GetSize(); ++j) ring.push_back(builder.vertex(corners[j], face_id));
        if (ring.size() > 3) ++mesh.polygons_split;
        // Poly loops run counter-clockwise about the face's outward normal.
        for (size_t j = 1; j + 1 < ring.size(); ++j) {
            mesh.indices.insert(mesh.indices.end(), {ring[0], ring[j], ring[j + 1]});
            mesh.face_ids.push_back(face_id);
        }
    }
    return SM_SUCCESS;
}

using FaceIndex = std::unordered_map<const SmFace*, uint32_t>;

FaceIndex index_faces(const SmBrep& brep) {
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    FaceIndex face_index;
    for (ULONG i = 0; i < faces.GetSize(); ++i) face_index.emplace(faces[i], static_cast<uint32_t>(i));
    return face_index;
}

// The faces an edge joins: the first two, and how many faces use it.
struct EdgeFaces {
    SmFace* faces[2] = {nullptr, nullptr};
    uint32_t use_count = 0;
};

EdgeFaces faces_of(const SmEdge* edge) {
    EdgeFaces result;
    SmTArray<SmEdgeuse*> uses;
    edge->GetEdgeuses(uses);
    // Each use of the edge by a face comes as an edgeuse and its mate, one
    // per side of the face; count each pair once. A seam is used twice by
    // the same face.
    std::vector<const SmEdgeuse*> counted;
    for (ULONG j = 0; j < uses.GetSize(); ++j) {
        const SmEdgeuse* use = uses[j];
        if (std::find(counted.begin(), counted.end(), use->GetMate()) != counted.end()) continue;
        counted.push_back(use);
        if (counted.size() <= 2) result.faces[counted.size() - 1] = use->GetFace();
    }
    result.use_count = static_cast<uint32_t>(counted.size());
    return result;
}

int32_t face_id(const FaceIndex& face_index, const SmFace* face) {
    const auto found = face_index.find(face);
    return found == face_index.end() ? -1 : static_cast<int32_t>(found->second);
}

sm_status build_edges(SmBrep& brep, double angle_tol_deg, const FaceIndex& face_index, sm_mesh& mesh) {
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    std::unordered_map<const SmEdge*, size_t> edge_index;
    mesh.edges.resize(edges.GetSize());
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        edge_index.emplace(edges[i], i);
        sm_edge_info& info = mesh.edges[i];
        std::memset(&info, 0, sizeof info);
        const EdgeFaces faces = faces_of(edges[i]);
        info.face_a = faces.faces[0] ? face_id(face_index, faces.faces[0]) : -1;
        info.face_b = faces.faces[1] ? face_id(face_index, faces.faces[1]) : -1;
        info.use_count = faces.use_count;
        info.label = read_label(edges[i]);
    }

    SmTArray<SmPoint3d> points;
    SmTArray<ULONG> counts;
    SmTArray<SmEdge*> sampled;
    const SmApiStatus status = SmApiTessellateBoundaries(&brep, angle_tol_deg, points, counts, sampled);
    if (status != SM_SUCCESS) return status;

    // Copy each sampled polyline to its edge, in the brep's edge order.
    std::vector<std::pair<ULONG, ULONG>> ranges(edges.GetSize(), {0, 0});
    ULONG start = 0;
    for (ULONG i = 0; i < counts.GetSize() && i < sampled.GetSize(); ++i) {
        const auto found = edge_index.find(sampled[i]);
        if (found != edge_index.end()) ranges[found->second] = {start, counts[i]};
        start += counts[i];
    }
    for (size_t i = 0; i < ranges.size(); ++i) {
        mesh.edges[i].first_point = static_cast<uint32_t>(mesh.edge_points.size() / 3);
        mesh.edges[i].point_count = static_cast<uint32_t>(ranges[i].second);
        for (ULONG k = 0; k < ranges[i].second; ++k) {
            const SmPoint3d& p = points[ranges[i].first + k];
            mesh.edge_points.insert(mesh.edge_points.end(),
                                    {static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z)});
        }
    }
    return SM_SUCCESS;
}

// The face's outward unit normal at the point on it nearest `near`, as
// usd_brep's Face.outward_normal_at_uv orients it: the surface normal,
// flipped when the face's upward faceuse isn't in the infinite region.
bool outward_normal(const SmFace* face, const SmPoint3d& near, double out[3]) {
    const SmSurface* surface = face->GetSurface();
    if (!surface) return false;
    SmPoint3d closest;
    SmTArray<SmPoint2d> uvs;
    double distance = 0.0;
    if (SmApiSurfaceClosestPoint(surface, near, closest, uvs, distance) != SM_SUCCESS || uvs.GetSize() == 0) {
        return false;
    }
    SmVector3d normal;
    if (SmApiEvaluateSurfaceNormal(surface, SmVector2d(uvs[0].x, uvs[0].y), normal) != SM_SUCCESS) return false;
    const SmBrep* brep = face->GetBrep();
    const SmFaceuse* upward = face->GetUpwardFaceuse();
    if (brep && upward && upward->GetShell() && brep->GetInfiniteRegion() != upward->GetShell()->GetRegion()) {
        normal *= -1.0;
    }
    put3(out, normal);
    return true;
}

sm_status describe_edge(const SmEdge* edge, const FaceIndex& face_index, sm_edge_geometry& info) {
    std::memset(&info, 0, sizeof info);
    double t0 = 0.0;
    double t1 = 0.0;
    SmApiStatus status = SmApiEdgeParameterRange(edge, t0, t1);
    if (status != SM_SUCCESS) return status;
    // A fillet can leave an edge's range an ulp past its curve's domain, which
    // SmApiEvaluateCurve rejects, so samples stay inside both.
    const SmExtent1d domain = edge->GetCurve() ? edge->GetCurve()->GetNaturalInterval() : SmExtent1d(t0, t1);
    for (int k = 0; k < SM_EDGE_SAMPLES; ++k) {
        SmVector3d point;
        SmVector3d* point_out = &point;
        SmVector3d* no_derivative = nullptr;
        SmVector3d* no_second = nullptr;
        const double t =
            std::min(std::max(t0 + (t1 - t0) * k / (SM_EDGE_SAMPLES - 1), domain.GetMin()), domain.GetMax());
        status = SmApiEvaluateCurve(edge->GetCurve(), t, point_out, no_derivative, no_second);
        if (status != SM_SUCCESS) return status;
        put3(info.points[k], point);
    }

    const EdgeFaces faces = faces_of(edge);
    info.use_count = faces.use_count;
    info.label = read_label(edge);
    info.label_a = read_label(faces.faces[0]);
    info.label_b = read_label(faces.faces[1]);
    const double* mid = info.points[SM_EDGE_SAMPLES / 2];
    const SmPoint3d midpoint(mid[0], mid[1], mid[2]);
    int32_t* ids[2] = {&info.face_a, &info.face_b};
    double* normals[2] = {info.normal_a, info.normal_b};
    for (int side = 0; side < 2; ++side) {
        const SmFace* face = faces.faces[side];
        *ids[side] = face ? face_id(face_index, face) : -1;
        if (face && !outward_normal(face, midpoint, normals[side])) {
            normals[side][0] = normals[side][1] = normals[side][2] = 0.0;
        }
    }
    return SM_SUCCESS;
}

bool is_manifold_solid(const SmBrep& brep) {
    SmBoolean manifold = FALSE;
    return SmApiBrepIsManifoldSolid(&brep, manifold) == SM_SUCCESS && manifold;
}

// Runs a fillet operation and holds the kernel to its word: some report
// success yet leave a manifold solid open (see sm_brep_chamfer_edges).
template <typename F>
sm_status keeping_solid(SmBrep& brep, F&& fillet) {
    const bool was_solid = is_manifold_solid(brep);
    const sm_status status = fillet();
    if (status != SM_SUCCESS) return status;
    if (was_solid && !is_manifold_solid(brep)) return SM_SHIM_ERR_BROKE_SOLID;
    return SM_SUCCESS;
}

// Resolves edge indices into the brep's current edge list. Repeats are
// rejected: the kernel would fillet the edge twice.
sm_status edges_at(SmBrep& brep, const uint32_t* indices, size_t count, SmTArray<SmEdge*>& out) {
    SmTArray<SmEdge*> all;
    brep.GetEdges(all);
    std::vector<bool> taken(all.GetSize(), false);
    for (size_t i = 0; i < count; ++i) {
        const uint32_t index = indices[i];
        if (index >= all.GetSize() || taken[index]) return SM_SHIM_ERR_BAD_EDGE;
        taken[index] = true;
        out.Add(all[index]);
    }
    return SM_SUCCESS;
}

// SmApiChamferFillet's set-up, loading solvers for the picked edges only.
// The executive deletes the solvers; the generator and the self-intersection
// handler only need to outlive it.
//
// Unlike SmApiChamferFillet, this leaves global merge off. The executive
// already merges the chamfer in by topology insertion, and its header says
// only one merge mode should be on; with both, a chamfer on a circular edge
// can come back as an open shell. SmApiFilletEdgesPerEdge also ends up with
// global merge off, since SetFilletParameters resets it. Corners that need
// global merge turn it on themselves (SmFilletNx1Corner).
sm_status chamfer(SmBrep& target, SmTArray<SmEdge*>& picked, double width) {
    const SmContext* context = SmApiGetOrCreateContext();
    SmLinearCrossSectionFSG linear;
    SmMakeSurfaceBlendSIH self_intersections(1.0);
    SmFilletExecutive executive(*context, &target);
    executive.SetSelfIntersectionHandler(&self_intersections);
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        SmEdgeuse* use = picked[i]->GetBlendEdgeuse();
        if (!use) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
        SmFilletSolver* solver = new (*context) SmConstantDistanceFS(
            *context, target.GetTolerance(), 30.0 * SM_PI / 180.0, 2.0 * SM_PI / 180.0, width, use);
        solver->SetFilletSurfaceGenerator(&linear);
        const SmStatus loaded = executive.LoadFilletSolver(solver);
        if (loaded != SM_SUCCESS) return static_cast<sm_status>(loaded);
    }
    SmStatus result = executive.CreateFilletCorners();
    if (result != SM_SUCCESS) return static_cast<sm_status>(result);
    SmTArray<SmFilletCorner*>& corners = executive.GetFilletCorners();
    for (ULONG i = 0; i < corners.GetSize(); ++i) corners[i]->SetBevel(TRUE);
    return static_cast<sm_status>(executive.DoFilleting());
}

SmFace* face_at(const SmBrep& brep, size_t index) {
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    return index < faces.GetSize() ? faces[static_cast<ULONG>(index)] : nullptr;
}

// Planarity tolerance, in model units: how far a surface's control points may
// stray from its plane.
const double PLANAR_TOLERANCE = 1e-7;

// The normal of a planar face pointing out of the solid, if it's planar. A
// plane's normal is the same everywhere, so it's evaluated mid-domain and then
// oriented as outward_normal orients it. SmSurface::IsPlanar bounds the
// surface with a box no thinner than about 1e-5, so a NURBS surface whose
// control points share a plane (a tapered extrusion's walls) is checked
// directly (plane_of).
bool plane_of(const SmSurface* surface, SmPoint3d& point, SmVector3d& normal);

bool outward_plane(const SmFace* face, double normal[3]) {
    const SmSurface* surface = face->GetSurface();
    SmPoint3d on_plane;
    SmVector3d plane_normal;
    if (!surface || (!surface->IsPlanar(PLANAR_TOLERANCE) && !plane_of(surface, on_plane, plane_normal))) return false;
    const SmExtent2d domain = face->GetUVDomain();
    const SmVector2d middle(0.5 * (domain.GetUMin() + domain.GetUMax()), 0.5 * (domain.GetVMin() + domain.GetVMax()));
    SmVector3d direction;
    if (SmApiEvaluateSurfaceNormal(surface, middle, direction) != SM_SUCCESS) return false;
    const SmBrep* brep = face->GetBrep();
    const SmFaceuse* upward = face->GetUpwardFaceuse();
    if (brep && upward && upward->GetShell() && brep->GetInfiniteRegion() != upward->GetShell()->GetRegion()) {
        direction *= -1.0;
    }
    put3(normal, direction);
    return true;
}

sm_status describe_face_geometry(const SmFace* face, sm_face_geometry& out) {
    std::memset(&out, 0, sizeof out);
    out.info = describe_face(face);
    SmExtent3d box;
    const SmStatus status = face->CalculateBoundingBox(box);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    put3(out.bounds_min, box.GetMin());
    put3(out.bounds_max, box.GetMax());
    out.planar = outward_plane(face, out.normal) ? 1 : 0;
    if (!out.planar) out.normal[0] = out.normal[1] = out.normal[2] = 0.0;
    return SM_SUCCESS;
}

// As SmApiFaceComputeCentroid, but without its 1e-4 accuracy floor.
sm_status face_centroid(const SmFace* face, double accuracy, double& area, SmPoint3d& centroid) {
    SmExtent3d box;
    SmStatus status = face->CalculateBoundingBox(box);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    const SmPoint3d origin = box.GetMid();
    const SmVector3d size = box.GetSize();
    const double extent = std::max(size.z, std::max(size.x, size.y));
    double signed_area = 0.0;
    double unused_volume = 0.0;
    SmTArray<SmVector3d> moments;
    status = face->ComputePreciseProperties(SM_OT_SAME, accuracy, origin, extent * extent / 100.0, 0.0, signed_area,
                                            unused_volume, moments, SM_PPF_AREA | SM_PPF_AREA_FIRST_MOMENTS);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    if (std::fabs(signed_area) < SM_EFF_ZERO || moments.GetSize() < 1) return SM_ERR_INVALID_INPUT;
    area = std::fabs(signed_area);
    centroid = origin + moments[0] / signed_area;
    return SM_SUCCESS;
}

bool unit_vector(const sm_vec3& v) {
    const double length = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return std::isfinite(length) && std::fabs(length - 1.0) <= 1e-9;
}

SmPoint3d midpoint_of(const SmEdge* edge) {
    SmPoint3d point(0.0, 0.0, 0.0);
    const SmExtent1d interval = edge->GetInterval();
    edge->GetCurve()->EvaluatePoint(0.5 * (interval.GetMin() + interval.GetMax()), point);
    return point;
}

// Which edge of the swept face each side face came from. The sweep works on
// copies of the face's edge curves, so the link is geometric: a side's edge
// in the start cap lies on the source edge it was swept from.
int32_t source_edge(const SmFace* side, const SmFace* start, const SmTArray<SmEdge*>& source_edges,
                    const std::unordered_map<const SmEdge*, uint32_t>& edge_index) {
    SmTArray<SmEdge*> side_edges;
    side->GetEdges(side_edges);
    SmTArray<SmEdge*> start_edges;
    if (start) start->GetEdges(start_edges);
    int32_t best = -1;
    double best_distance = SM_BIG_DOUBLE;
    for (ULONG i = 0; i < side_edges.GetSize(); ++i) {
        ULONG unused;
        if (start && !start_edges.FindElement(side_edges[i], unused)) continue;
        const SmPoint3d mid = midpoint_of(side_edges[i]);
        for (ULONG j = 0; j < source_edges.GetSize(); ++j) {
            SmPoint3d closest;
            double parameter = 0.0;
            double distance = SM_BIG_DOUBLE;
            if (SmApiCurveClosestPoint(source_edges[j]->GetCurve(), mid, closest, parameter, distance) != SM_SUCCESS) {
                continue;
            }
            if (distance < best_distance) {
                best_distance = distance;
                const auto found = edge_index.find(source_edges[j]);
                best = found == edge_index.end() ? -1 : static_cast<int32_t>(found->second);
            }
        }
    }
    // The copy lies on the source curve to well within the sweep's tolerance.
    return best_distance <= 1e3 * SM_ZONE_TOL_3D ? best : -1;
}

// Puts `analytic` under the face in place of its surface, which it must
// match in shape and parameters. SmFace::SetSurface drops the edges' UV trim
// curves, so this keeps copies and puts them back; the face keeps its domain.
void replace_surface(SmFace* face, std::unique_ptr<SmSurface> analytic) {
    const SmContext* context = face->GetBrep()->GetContext();
    // Each edgeuse pair holds its trim curve on one of the two.
    std::vector<std::pair<SmEdgeuse*, std::unique_ptr<SmBSplineCurve>>> trims;
    SmTArray<SmEdgeuse*> uses;
    face->GetUpwardEdgeuses(uses);
    for (ULONG j = 0; j < uses.GetSize(); ++j) {
        SmEdgeuse* owner = nullptr;
        const SmBSplineCurve* trim = uses[j]->GetUVTrimCurve(&owner);
        const bool kept = std::any_of(trims.begin(), trims.end(), [&](const auto& t) { return t.first == owner; });
        if (!trim || !owner || kept) continue;
        SmCurve* copy = nullptr;
        if (trim->Copy(*context, copy) != SM_SUCCESS || !copy) {
            delete copy;
            continue;
        }
        trims.emplace_back(owner, static_cast<SmBSplineCurve*>(copy));
    }
    const SmExtent2d domain = face->GetUVDomain();
    face->SetSurface(analytic.release(), TRUE);
    face->SetUVDomain(domain);
    for (auto& [owner, trim] : trims) owner->SetUVTrimCurve(trim.release(), 0.0, TRUE);
}

// CreateCurveSweepFromFaces makes every side a NURBS surface, so lines sweep
// into planar NURBS and arcs and circles into cylindrical ones;
// CreateTaperExtrude's ruled sides are planar or conical NURBS. Gives each
// side that is a plane, a cylinder or a cone that analytic surface instead
// (SmBSplineSurface::CopyAndAddAnalytics, as SmPrimitiveCreation makes its
// planes). SMLib's analytic surfaces are B-splines too: the copy is checked
// against the NURBS and keeps its parameters, so the face keeps its domain
// (replace_surface). Any other side keeps its NURBS.
void make_sides_analytic(SmBrep& brep, const SmTArray<SmFace*>& sides) {
    const SmContext* context = brep.GetContext();
    SmTemporaryChangeValue<SmBoolean> editing(brep.m_bEditingEnabled, TRUE);
    for (ULONG i = 0; i < sides.GetSize(); ++i) {
        SmFace* side = sides[i];
        const SmSurface* surface = side->GetSurface();
        if (!surface || surface->GetType() != SmBSplineSurface_TYPE) continue;
        SmSurface* raw = nullptr;
        if (surface->CopyAndAddAnalytics(*context, raw) != SM_SUCCESS || !raw) {
            delete raw;
            continue;
        }
        std::unique_ptr<SmSurface> analytic(raw);
        const SM_TYPE type = analytic->GetType();
        const bool round = type == SmCylinder_TYPE || type == SmCone_TYPE;
        if ((type != SmPlane_TYPE && !round) || !(analytic->GetNaturalUVDomain() == surface->GetNaturalUVDomain())) {
            continue;
        }
        replace_surface(side, std::move(analytic));
    }
}

sm_status sweep_face(const SmBrep& brep, SmFace* face, double distance, SmBrep*& prism,
                     std::vector<sm_swept_face>& roles) {
    double normal[3];
    if (!outward_plane(face, normal)) return SM_SHIM_ERR_NOT_PLANAR;
    SmPoint3d inside;
    if (SmApiFaceInternalPoint(face, inside) != SM_SUCCESS) return SM_SHIM_ERR_NOT_PLANAR;
    const SmVector3d direction(normal[0], normal[1], normal[2]);

    const SmContext* context = SmApiGetOrCreateContext();
    SmBSplineCurve* raw_path = nullptr;
    SmStatus status = SmBSplineCurve::CreateLineSegment(*context, 3, inside, inside + direction * distance, raw_path);
    if (status != SM_SUCCESS || !raw_path) return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
    const std::unique_ptr<SmBSplineCurve> path(raw_path);

    // Translate the face's own loops along the path and cap both ends, in
    // place rather than moved to the path's start.
    SmSweepOptions options;
    options.bCapEndsArg = TRUE;
    options.bTranslationalSweep = TRUE;
    options.bMoveToProfile = TRUE;
    std::unique_ptr<SmBrep> result(new (*context) SmBrep());
    SmPrimitiveCreation creator(result->GetInfiniteRegion());
    SmTArray<SmFace*> profile;
    profile.Add(face);
    SmTArray<SmFace*> starts, sides, ends;
    status = creator.CreateCurveSweepFromFaces(profile, path.get(), nullptr, nullptr, SM_ZONE_TOL_3D, &options, starts,
                                               sides, ends);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    if (!is_manifold_solid(*result) || starts.GetSize() != 1 || ends.GetSize() != 1) return SM_SHIM_ERR_BAD_SWEEP;
    make_sides_analytic(*result, sides);

    SmTArray<SmEdge*> all_edges;
    brep.GetEdges(all_edges);
    std::unordered_map<const SmEdge*, uint32_t> edge_index;
    for (ULONG i = 0; i < all_edges.GetSize(); ++i) edge_index.emplace(all_edges[i], i);
    SmTArray<SmEdge*> face_edges;
    face->GetEdges(face_edges);

    SmTArray<SmFace*> faces;
    result->GetFaces(faces);
    roles.assign(faces.GetSize(), sm_swept_face{SM_SWEPT_SIDE, -1});
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        ULONG unused;
        if (faces[i] == starts[0]) {
            roles[i].role = SM_SWEPT_START;
        } else if (faces[i] == ends[0]) {
            roles[i].role = SM_SWEPT_END;
        } else if (sides.FindElement(faces[i], unused)) {
            roles[i].source_edge = source_edge(faces[i], starts[0], face_edges, edge_index);
        } else {
            return SM_SHIM_ERR_BAD_SWEEP;
        }
    }
    prism = result.release();
    return SM_SUCCESS;
}

// A curve of sm_brep_create_planar_faces, made in the plane z = 0, where
// SM_API makes its lines, arcs and circles.
sm_status plane_curve(const sm_plane_curve& curve, SmCurve*& out) {
    out = nullptr;
    const SmPoint3d a(curve.a[0], curve.a[1], 0.0);
    if (!std::isfinite(a.x) || !std::isfinite(a.y)) return SM_ERR_INVALID_INPUT;
    const bool round = curve.kind == SM_PLANE_ARC || curve.kind == SM_PLANE_CIRCLE;
    if (round && !(std::isfinite(curve.radius) && curve.radius > 0.0)) return SM_ERR_INVALID_INPUT;
    SmApiStatus status = SM_ERR_INVALID_INPUT;
    if (curve.kind == SM_PLANE_LINE) {
        const SmPoint3d b(curve.b[0], curve.b[1], 0.0);
        if (!std::isfinite(b.x) || !std::isfinite(b.y)) return SM_ERR_INVALID_INPUT;
        SmLine* line = nullptr;
        status = SmApiCreateLineSegment(a, b, line);
        out = line;
    } else if (curve.kind == SM_PLANE_ARC) {
        const double sweep = std::fabs(curve.sweep_deg);
        if (!std::isfinite(curve.start_deg) || !(sweep > 0.0 && sweep <= 360.0)) return SM_ERR_INVALID_INPUT;
        // SmApiCreateArc runs counter-clockwise between angles in [-360, 360];
        // the faces don't mind which way each curve runs.
        double from = std::fmod(std::min(curve.start_deg, curve.start_deg + curve.sweep_deg), 360.0);
        if (from < 0.0) from += 360.0;
        if (from + sweep > 360.0) from -= 360.0;
        SmBSplineCurve* arc = nullptr;
        status = SmApiCreateArc(a, curve.radius, from, from + sweep, arc);
        out = arc;
    } else if (curve.kind == SM_PLANE_CIRCLE) {
        SmBSplineCurve* circle = nullptr;
        status = SmApiCreateCircle(a, curve.radius, circle);
        out = circle;
    }
    if (status != SM_SUCCESS) {
        delete out;
        out = nullptr;
        return static_cast<sm_status>(status);
    }
    return out ? SM_SUCCESS : SM_ERR;
}

// The input curve an edge of a new planar face lies on: the one it holds, or
// failing that the one nearest its midpoint.
long curve_of(const SmEdge* edge, const SmTArray<SmCurve*>& curves) {
    for (ULONG i = 0; i < curves.GetSize(); ++i) {
        if (edge->GetCurve() == curves[i]) return static_cast<long>(i);
    }
    const SmPoint3d mid = midpoint_of(edge);
    long best = -1;
    double best_distance = SM_BIG_DOUBLE;
    for (ULONG i = 0; i < curves.GetSize(); ++i) {
        SmPoint3d closest;
        double parameter = 0.0;
        double distance = SM_BIG_DOUBLE;
        if (SmApiCurveClosestPoint(curves[i], mid, closest, parameter, distance) == SM_SUCCESS &&
            distance < best_distance) {
            best_distance = distance;
            best = static_cast<long>(i);
        }
    }
    return best_distance <= 1e3 * SM_ZONE_TOL_3D ? best : -1;
}

// Planar faces from curves in the plane z = 0, each facing +z, with each
// curve's edge labelled.
sm_status planar_faces(const sm_plane_curve* input, size_t count, SmBrep*& out) {
    SmTArray<SmCurve*> curves;
    // The curves are the caller's to free until the faces take them.
    struct Free {
        SmTArray<SmCurve*>& curves;
        bool owned = true;
        ~Free() {
            for (ULONG i = 0; owned && i < curves.GetSize(); ++i) delete curves[i];
        }
    } free_curves{curves};
    for (size_t i = 0; i < count; ++i) {
        if (input[i].label.present && role_length(input[i].label) == 0) return SM_SHIM_ERR_BAD_LABEL;
        SmCurve* curve = nullptr;
        const sm_status status = plane_curve(input[i], curve);
        if (status != SM_SUCCESS) return status;
        curves.Add(curve);
    }
    SmBrep* raw = nullptr;
    const SmApiStatus status = SmApiCreatePlanarFaces(curves, raw);
    if (status != SM_SUCCESS || !raw) return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
    free_curves.owned = false;
    std::unique_ptr<SmBrep> brep(raw);

    SmTArray<SmEdge*> edges;
    brep->GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        const long curve = curve_of(edges[i], curves);
        if (curve < 0) return SM_ERR;
        write_label(*edges[i], input[curve].label);
    }
    // The kernel picks each face's side from the curves; turn those facing
    // -z over. SwapUV mirrors the surface's parameters and so its normal.
    SmTArray<SmFace*> faces;
    brep->GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        double normal[3];
        if (!outward_plane(faces[i], normal)) return SM_SHIM_ERR_NOT_PLANAR;
        if (normal[2] < 0.0) {
            const SmStatus swapped = faces[i]->SwapUV();
            if (swapped != SM_SUCCESS) return static_cast<sm_status>(swapped);
            if (!outward_plane(faces[i], normal) || normal[2] < 0.0) return SM_ERR;
        }
    }
    out = brep.release();
    return SM_SUCCESS;
}

// The largest distance from a vertex to its edges or faces, or from an edge
// to its faces (SmBrep::ValidateAndUpdateTolerances, validating only).
sm_status largest_gap(SmBrep& brep, double& out) {
    SmBoolean violated = FALSE;
    double trim = 0.0;
    double vertex_edge = 0.0;
    double vertex_face = 0.0;
    double edge_face = 0.0;
    const SmStatus status = brep.ValidateAndUpdateTolerances(TRUE, violated, trim, vertex_edge, vertex_face, NULL,
                                                             SM_LEVEL_0, NULL, NULL, NULL, NULL, NULL, NULL,
                                                             &edge_face, NULL);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    out = std::max({vertex_edge, vertex_face, edge_face});
    return SM_SUCCESS;
}

// The volume each face contributes about `origin`, as SmApiBrepComputeVolume
// sums it (SmBrep::ComputePreciseProperties, from the faces' trim curves),
// but integrating the volume alone. The integrator's step control follows
// the face's area whatever is asked for, so each face's volume comes out
// exactly as it does with every property requested. A face with solid on
// both sides contributes from both; one with the same region on both sides
// contributes nothing, as there.
sm_status face_volumes(const SmTArray<SmFace*>& faces, const SmPoint3d& origin, double estimated_area, double accuracy,
                       std::vector<double>& out) {
    out.assign(faces.GetSize(), 0.0);
    SmTArray<SmVector3d> unused_moments;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        SmFaceuse* uses[2] = {faces[i]->GetUpwardFaceuse(), faces[i]->GetUpwardFaceuse()->GetMate()};
        if (uses[0]->GetShell()->GetRegion() == uses[1]->GetShell()->GetRegion()) continue;
        for (SmFaceuse* use : uses) {
            if (use->GetShell()->GetRegion()->IsVoid()) continue;
            double unused_area = 0.0;
            double volume = 0.0;
            const SmStatus status = faces[i]->ComputePreciseProperties(use->GetOrientation(), accuracy, origin,
                                                                       estimated_area, 0.0, unused_area, volume,
                                                                       unused_moments, SM_PPF_VOLUME);
            if (status != SM_SUCCESS) return static_cast<sm_status>(status);
            out[i] += volume;
        }
    }
    return SM_SUCCESS;
}

// Half the integral of (x - x0) dy - (y - y0) dx along the edge, start to
// end: summed around a face's loops, the z component of its vector area.
// Gauss-Legendre on equal pieces of the edge's range, exact for lines and
// far tighter than the check needs for arcs and splines.
sm_status edge_swept_area_z(const SmEdge* edge, const SmPoint3d& origin, double& out) {
    out = 0.0;
    if (!edge->GetCurve()) return SM_ERR;
    double t0 = 0.0;
    double t1 = 0.0;
    SmApiStatus status = SmApiEdgeParameterRange(edge, t0, t1);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    // As in describe_edge, a fillet can leave the range an ulp past the
    // curve's domain.
    const SmExtent1d domain = edge->GetCurve()->GetNaturalInterval();
    t0 = std::min(std::max(t0, domain.GetMin()), domain.GetMax());
    t1 = std::min(std::max(t1, domain.GetMin()), domain.GetMax());
    static const double nodes[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                                    0.8611363115940526};
    static const double weights[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                                      0.3478548451374538};
    const int pieces = 8;
    const double half = 0.5 * (t1 - t0) / pieces;
    for (int piece = 0; piece < pieces; ++piece) {
        const double middle = t0 + (2 * piece + 1) * half;
        for (int k = 0; k < 4; ++k) {
            SmVector3d point;
            SmVector3d tangent;
            SmVector3d* point_out = &point;
            SmVector3d* tangent_out = &tangent;
            SmVector3d* no_second = nullptr;
            status = SmApiEvaluateCurve(edge->GetCurve(), middle + nodes[k] * half, point_out, tangent_out, no_second);
            if (status != SM_SUCCESS) return static_cast<sm_status>(status);
            out += weights[k] * half * 0.5 * ((point.x - origin.x) * tangent.y - (point.y - origin.y) * tangent.x);
        }
    }
    return SM_SUCCESS;
}

// The z component of each face's vector area from its edges, the normal
// pointing out of the solid: around each solid side's faceuse, whose
// edgeuses run clockwise seen from outside.
sm_status face_areas_z(const SmBrep& brep, const SmTArray<SmFace*>& faces, const SmPoint3d& origin,
                       std::vector<double>& out) {
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    std::unordered_map<const SmEdge*, double> swept;
    swept.reserve(edges.GetSize());
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        double area = 0.0;
        const sm_status status = edge_swept_area_z(edges[i], origin, area);
        if (status != SM_SUCCESS) return status;
        swept.emplace(edges[i], area);
    }
    out.assign(faces.GetSize(), 0.0);
    SmTArray<SmEdgeuse*> edgeuses;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        SmFaceuse* uses[2] = {faces[i]->GetUpwardFaceuse(), faces[i]->GetUpwardFaceuse()->GetMate()};
        if (uses[0]->GetShell()->GetRegion() == uses[1]->GetShell()->GetRegion()) continue;
        for (SmFaceuse* use : uses) {
            if (use->GetShell()->GetRegion()->IsVoid()) continue;
            edgeuses.RemoveAll();
            use->GetEdgeuses(edgeuses);
            for (ULONG e = 0; e < edgeuses.GetSize(); ++e) {
                const auto found = swept.find(edgeuses[e]->GetEdge());
                if (found == swept.end()) return SM_ERR;
                out[i] += edgeuses[e]->GetOrientation() == SM_OT_OPPOSITE ? found->second : -found->second;
            }
        }
    }
    return SM_SUCCESS;
}

// What the brep's tessellation encloses: coarse, by angle alone, so an exact
// solid's mesh is within about 1% of its volume.
sm_status tessellated_volume(SmBrep& brep, double& enclosed) {
    enclosed = 0.0;
    SmTessellationParams quality;
    quality.dChordHeightTol = 0.0;
    quality.dCurveAngleTolDeg = 15.0;
    quality.dSurfaceAngleTolDeg = 15.0;
    quality.dMaxEdgeLength = 0.0;
    quality.dMaxAspectRatio = 0.0;
    SmPolyBrep* poly = nullptr;
    const SmApiStatus status = SmApiTessellate(&brep, poly, quality);
    if (status != SM_SUCCESS || !poly) return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
    const std::unique_ptr<SmPolyBrep, PolyBrepDeleter> owned(poly);
    SmTArray<SmPolyFace*> polygons;
    poly->GetPolyFaces(polygons);
    SmTArray<SmPolyLoop*> loops;
    SmTArray<SmPolyEdge*> corners;
    for (ULONG i = 0; i < polygons.GetSize(); ++i) {
        loops.RemoveAll();
        polygons[i]->GetPolyLoops(loops);
        for (ULONG l = 0; l < loops.GetSize(); ++l) {
            corners.RemoveAll();
            loops[l]->GetPolyEdges(corners);
            if (corners.GetSize() < 3) continue;
            const SmPoint3d first = corners[0]->GetStartPolyVertex()->GetPoint();
            for (ULONG k = 1; k + 1 < corners.GetSize(); ++k) {
                const SmPoint3d b = corners[k]->GetStartPolyVertex()->GetPoint();
                const SmPoint3d c = corners[k + 1]->GetStartPolyVertex()->GetPoint();
                enclosed += SmVector3d(first.x, first.y, first.z).Dot(SmVector3d(b.x, b.y, b.z) * SmVector3d(c.x, c.y, c.z)) / 6.0;
            }
        }
    }
    return SM_SUCCESS;
}

// Whether the brep's own volume (SmApiBrepComputeVolume, from its faces' trim
// curves) can be trusted. The offset executive can return a solid that
// tessellates right and passes the gap checks but whose faces integrate
// wrongly; mass properties, and booleans after them, then go wrong. Seen with
// fillets sliding along drafted walls: the fillets' integrals over their
// trims disagree with their own edges, and the volume comes out 20-30% short.
//
// SMLib integrates volume as the sum over faces of the integral of
// (z - z0) n_z dA, so raising the origin by h changes a face's share by -h
// times its projected area, the integral of n_z dA; for a face that
// integrates right, that is also what its edges enclose (Stokes). So the
// volume is integrated about two origins a box height apart, and each face's
// change is held to its edges. That costs two integrations, against one and
// a tessellation. When a face disagrees, or the volume isn't positive, the
// old test decides: the volume must agree to 2% with what the brep's
// tessellation encloses, above a coarse mesh's chord error and well inside
// the failures seen.
sm_status volume_agrees(SmBrep& brep, bool& agrees) {
    agrees = false;
    // SmApiBrepComputeVolume's origin, area estimate and accuracy.
    SmExtent3d box;
    SmStatus kernel_status = brep.CalculateBoundingBox(box);
    if (kernel_status != SM_SUCCESS) return static_cast<sm_status>(kernel_status);
    const SmVector3d size = box.GetSize();
    const double extent = std::max(size.z, std::max(size.x, size.y));
    const double estimated_area = extent * extent / 100.0;
    const double accuracy = 1e-4;
    const SmPoint3d origin = box.GetMid();
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    std::vector<double> about;
    sm_status status = face_volumes(faces, origin, estimated_area, accuracy, about);
    if (status != SM_SUCCESS) return status;
    double volume = 0.0;
    for (double share : about) volume += share;

    // Anything the cheaper test can't measure goes to the tessellation too.
    const double height = size.z;
    std::vector<double> raised;
    std::vector<double> projected;
    if (volume > 0.0 && height > 0.0 &&
        face_volumes(faces, origin + SmVector3d(0.0, 0.0, height), estimated_area, accuracy, raised) ==
            SM_SUCCESS &&
        face_areas_z(brep, faces, origin, projected) == SM_SUCCESS) {
        double disagreement = 0.0;
        double scale = 0.0;
        for (ULONG i = 0; i < faces.GetSize(); ++i) {
            disagreement += std::fabs(about[i] - raised[i] - height * projected[i]);
            scale += std::fabs(about[i]) + std::fabs(raised[i]);
        }
        if (disagreement <= 1e-3 * scale) {
            agrees = true;
            return SM_SUCCESS;
        }
    }

    double enclosed = 0.0;
    status = tessellated_volume(brep, enclosed);
    if (status != SM_SUCCESS) return status;
    agrees = std::fabs(volume - enclosed) <= 2e-2 * std::fabs(enclosed);
    return SM_SUCCESS;
}

// Points along an edge's curve, ends included, within the curve's domain.
bool edge_points(const SmEdge* edge, int count, std::vector<SmPoint3d>& out) {
    double t0 = 0.0;
    double t1 = 0.0;
    if (!edge->GetCurve() || SmApiEdgeParameterRange(edge, t0, t1) != SM_SUCCESS) return false;
    const SmExtent1d domain = edge->GetCurve()->GetNaturalInterval();
    for (int k = 0; k < count; ++k) {
        SmVector3d point;
        SmVector3d* point_out = &point;
        SmVector3d* no_derivative = nullptr;
        SmVector3d* no_second = nullptr;
        const double t = std::min(std::max(t0 + (t1 - t0) * k / (count - 1), domain.GetMin()), domain.GetMax());
        if (SmApiEvaluateCurve(edge->GetCurve(), t, point_out, no_derivative, no_second) != SM_SUCCESS) return false;
        out.emplace_back(point.x, point.y, point.z);
    }
    return true;
}

// Whether every planar face's loops keep apart: no two cross, and every
// inner loop lies inside the outer one. A local operation can move a face's
// outer edge across features standing on it (a boss near a lid's edge when
// the wall under that edge is pushed in), which needs new topology; the
// executive leaves the loops crossing, or the boss's rim outside the lid,
// and the solid still passes for manifold, with a volume a few hundredths
// off. Each loop is a polyline in the face's plane through points along its
// edges.
bool loops_stay_apart(const SmBrep& brep) {
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        SmTArray<SmLoop*> loops;
        faces[i]->GetLoops(loops);
        SmPoint3d origin;
        SmVector3d normal;
        if (loops.GetSize() < 2 || !faces[i]->GetSurface() || !plane_of(faces[i]->GetSurface(), origin, normal)) {
            continue;
        }
        SmVector3d u = std::fabs(normal.x) < 0.9 ? SmVector3d(1.0, 0.0, 0.0) : SmVector3d(0.0, 1.0, 0.0);
        u = u - normal * u.Dot(normal);
        u = u / u.Length();
        const SmVector3d v = normal * u;  // cross product
        struct Polyline {
            std::vector<std::array<double, 4>> segments;  // (x0, y0, x1, y1)
            double box[4] = {HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
            bool outer = false;
        };
        std::vector<Polyline> polylines(loops.GetSize());
        for (ULONG l = 0; l < loops.GetSize(); ++l) {
            polylines[l].outer = loops[l]->IsOuterLoop();
            SmTArray<SmEdgeuse*> uses;
            loops[l]->GetEdgeuses(uses);
            for (ULONG e = 0; e < uses.GetSize(); ++e) {
                std::vector<SmPoint3d> points;
                if (!uses[e]->GetEdge() || !edge_points(uses[e]->GetEdge(), 24, points)) return true;
                for (size_t k = 0; k + 1 < points.size(); ++k) {
                    const SmVector3d a = points[k] - origin;
                    const SmVector3d b = points[k + 1] - origin;
                    std::array<double, 4> segment = {a.Dot(u), a.Dot(v), b.Dot(u), b.Dot(v)};
                    polylines[l].segments.push_back(segment);
                    double* box = polylines[l].box;
                    box[0] = std::min({box[0], segment[0], segment[2]});
                    box[1] = std::min({box[1], segment[1], segment[3]});
                    box[2] = std::max({box[2], segment[0], segment[2]});
                    box[3] = std::max({box[3], segment[1], segment[3]});
                }
            }
        }
        auto side = [](double ax, double ay, double bx, double by, double px, double py) {
            return (bx - ax) * (py - ay) - (by - ay) * (px - ax);
        };
        auto crosses = [&](const std::array<double, 4>& p, const std::array<double, 4>& q) {
            const double d1 = side(p[0], p[1], p[2], p[3], q[0], q[1]);
            const double d2 = side(p[0], p[1], p[2], p[3], q[2], q[3]);
            const double d3 = side(q[0], q[1], q[2], q[3], p[0], p[1]);
            const double d4 = side(q[0], q[1], q[2], q[3], p[2], p[3]);
            return ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0));
        };
        for (size_t a = 0; a < polylines.size(); ++a) {
            for (size_t b = a + 1; b < polylines.size(); ++b) {
                const double* p = polylines[a].box;
                const double* q = polylines[b].box;
                if (p[0] > q[2] || q[0] > p[2] || p[1] > q[3] || q[1] > p[3]) continue;
                for (const auto& s1 : polylines[a].segments) {
                    for (const auto& s2 : polylines[b].segments) {
                        if (crosses(s1, s2)) return false;
                    }
                }
            }
        }
        // Even-odd: a ray from an inner loop's first point along +x crosses
        // the outer loop an odd number of times when it starts inside.
        const auto outer = std::find_if(polylines.begin(), polylines.end(), [](const Polyline& p) { return p.outer; });
        if (outer == polylines.end()) continue;
        for (const Polyline& inner : polylines) {
            if (inner.outer || inner.segments.empty()) continue;
            const double x = inner.segments[0][0];
            const double y = inner.segments[0][1];
            bool inside = false;
            for (const auto& s : outer->segments) {
                if ((s[1] > y) != (s[3] > y) && x < s[0] + (y - s[1]) * (s[2] - s[0]) / (s[3] - s[1])) inside = !inside;
            }
            if (!inside) return false;
        }
    }
    return true;
}

// Runs a local operation and holds it to what it promises: the same faces and
// edges, every planar face's loops apart, still a manifold solid, every
// vertex and edge on its faces to the brep's tolerance, and a volume its own
// tessellation agrees with. The offset
// executive can report success while failing any of them.
template <typename F>
sm_status keeping_topology(SmBrep& brep, F&& operation) {
    const bool was_solid = is_manifold_solid(brep);
    const ULONG faces = brep.GetNumFaces();
    const ULONG edges = brep.GetNumEdges();
    const sm_status status = operation();
    if (status != SM_SUCCESS) return status;
    if (was_solid && !is_manifold_solid(brep)) return SM_SHIM_ERR_BROKE_SOLID;
    if (brep.GetNumFaces() != faces || brep.GetNumEdges() != edges) return SM_SHIM_ERR_CHANGED_TOPOLOGY;
    if (!loops_stay_apart(brep)) return SM_SHIM_ERR_CHANGED_TOPOLOGY;
    double gap = 0.0;
    sm_status measured = largest_gap(brep, gap);
    if (measured != SM_SUCCESS) return measured;
    if (gap > brep.GetTolerance()) return SM_SHIM_ERR_BROKE_SOLID;
    bool agrees = false;
    measured = volume_agrees(brep, agrees);
    if (measured != SM_SUCCESS) return measured;
    return agrees ? SM_SUCCESS : SM_SHIM_ERR_BROKE_SOLID;
}

// The plane a surface lies on, if it is planar: SmPlane, or a NURBS surface
// whose control points all lie on one plane (such as a tapered extrusion's
// walls), which makes it exactly flat. SmSurface::IsPlanar only bounds the
// surface with a box, about 1e-5 thick at best, so it just proposes the plane.
bool plane_of(const SmSurface* surface, SmPoint3d& point, SmVector3d& normal) {
    if (surface->GetType() == SmPlane_TYPE) {
        const SmAxis2Placement& placement = static_cast<const SmPlane*>(surface)->GetPosition();
        point = placement.GetOrigin();
        normal = placement.GetZAxis();
        return true;
    }
    if (surface->GetType() != SmBSplineSurface_TYPE || !surface->IsPlanar(1e-3, &point, &normal)) return false;
    ULONG u_count = 0;
    ULONG v_count = 0;
    SmTArray<SmPoint3d> poles;
    SmTArray<double> weights;
    if (static_cast<const SmBSplineSurface*>(surface)->GetControlPointNet(u_count, v_count, poles, weights) !=
        SM_SUCCESS) {
        return false;
    }
    for (ULONG i = 0; i < weights.GetSize(); ++i) {
        if (!(weights[i] > 0.0)) return false;
    }
    for (ULONG i = 0; i < poles.GetSize(); ++i) {
        if (std::fabs((poles[i] - point).Dot(normal)) > 1e-9 * (1.0 + (poles[i] - point).Length())) return false;
    }
    return true;
}

bool is_cylinder(const SmSurface* surface) {
    const SM_TYPE type = surface->GetType();
    return type == SmCylinder_TYPE || (type == SmCone_TYPE && static_cast<const SmCone*>(surface)->IsCylinder());
}

// Whether translating the surface along `direction` leaves it where it is: a
// plane containing the direction, or a cylinder along it.
bool invariant_under(const SmSurface* surface, const SmVector3d& direction) {
    SmPoint3d point;
    SmVector3d normal;
    if (plane_of(surface, point, normal)) return std::fabs(normal.Dot(direction)) < 1e-9;
    if (is_cylinder(surface)) {
        const SmVector3d axis = static_cast<const SmCone*>(surface)->GetPosition().GetZAxis();
        return 1.0 - std::fabs(axis.Dot(direction)) < 1e-12;
    }
    return false;
}

// Whether a face smoothly joined to a moving face is a blend that can ride
// along with it: a round (cylinder, sphere or torus, as the fillet executive
// makes them), or any other curved face that SmFace::IdentifyFillet takes
// for a fillet, such as the surface of revolution on a hole's rim.
bool is_blend(const SmFace* face, double tolerance) {
    const SmSurface* surface = face->GetSurface();
    const SM_TYPE type = surface->GetType();
    if (is_cylinder(surface) || type == SmSphere_TYPE || type == SmTorus_TYPE) return true;
    SmPoint3d point;
    SmVector3d normal;
    if (plane_of(surface, point, normal)) return false;
    SmBoolean fillet = FALSE;
    ULONG section = 0;
    ULONG solver = 0;
    SmSurfParamType rails = SM_SP_NEITHER;
    double radius = 0.0;
    double distance = 0.0;
    SmTArray<SmEdge*> low;
    SmTArray<SmEdge*> high;
    return face->IdentifyFillet(tolerance, 1.0, TRUE, TRUE, fillet, section, solver, rails, radius, distance, low,
                                high) == SM_SUCCESS &&
           fillet;
}

// Whether the two faces meet smoothly along the edge: their outward normals
// agree at its midpoint.
bool smooth_edge(const SmEdge* edge, const SmFace* a, const SmFace* b) {
    double t0 = 0.0;
    double t1 = 0.0;
    if (SmApiEdgeParameterRange(edge, t0, t1) != SM_SUCCESS) return false;
    SmVector3d mid;
    SmVector3d* point_out = &mid;
    SmVector3d* no_derivative = nullptr;
    SmVector3d* no_second = nullptr;
    if (SmApiEvaluateCurve(edge->GetCurve(), 0.5 * (t0 + t1), point_out, no_derivative, no_second) != SM_SUCCESS) {
        return false;
    }
    double na[3];
    double nb[3];
    const SmPoint3d at(mid.x, mid.y, mid.z);
    if (!outward_normal(a, at, na) || !outward_normal(b, at, nb)) return false;
    return na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2] > std::cos(0.5 * SM_PI / 180.0);
}

// The faces smoothly joined to `face` across its edges, each once.
std::vector<SmFace*> smooth_neighbours(SmFace* face) {
    std::vector<SmFace*> out;
    SmTArray<SmEdge*> edges;
    face->GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        const EdgeFaces sides = faces_of(edges[i]);
        if (sides.use_count != 2) continue;
        SmFace* other = sides.faces[0] == face ? sides.faces[1] : sides.faces[0];
        if (!other || other == face || std::find(out.begin(), out.end(), other) != out.end()) continue;
        if (smooth_edge(edges[i], face, other)) out.push_back(other);
    }
    return out;
}

// How an offset changes one face's surface.
struct Motion {
    enum Kind {
        kKept,    // left as it is
        kShift,   // translated by `shift`: picked planes and cones, and the blends on them
        kOffset,  // replaced by its offset: other picked curved faces
    };
    Kind kind = kKept;
    SmVector3d shift = SmVector3d(0.0, 0.0, 0.0);
    // A rigid move's turn, for the faces sm_brep_move_faces turns: a kShift
    // face's points then go to turn * p + shift. Offsets never turn.
    bool turned = false;
    double turn[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
};

using Motions = std::unordered_map<const SmFace*, Motion>;

Motion motion_of(const Motions& motions, const SmFace* face) {
    const auto found = motions.find(face);
    return found == motions.end() ? Motion() : found->second;
}

// A direction after the motion's turn.
SmVector3d turned_vector(const Motion& motion, const SmVector3d& v) {
    if (!motion.turned) return v;
    const double(&r)[3][3] = motion.turn;
    return SmVector3d(r[0][0] * v.x + r[0][1] * v.y + r[0][2] * v.z, r[1][0] * v.x + r[1][1] * v.y + r[1][2] * v.z,
                      r[2][0] * v.x + r[2][1] * v.y + r[2][2] * v.z);
}

// Where the motion takes a point.
SmPoint3d moved_point(const Motion& motion, const SmPoint3d& p) {
    if (!motion.turned) return p + motion.shift;
    const SmVector3d turned = turned_vector(motion, SmVector3d(p.x, p.y, p.z));
    return SmPoint3d(turned.x, turned.y, turned.z) + motion.shift;
}

// The motion as a placement for SmSurface::Transform and SmCurve::Transform,
// which take a point p to origin + p.x X + p.y Y + p.z Z.
SmAxis2Placement placement_of(const Motion& motion) {
    SmAxis2Placement placement;
    if (!motion.turned) {
        placement.Translate(motion.shift);
        return placement;
    }
    const SmVector3d x = turned_vector(motion, SmVector3d(1.0, 0.0, 0.0));
    const SmVector3d y = turned_vector(motion, SmVector3d(0.0, 1.0, 0.0));
    return SmAxis2Placement(SmPoint3d(motion.shift.x, motion.shift.y, motion.shift.z), x, y);
}

// Whether two rigid motions are the same, to 1e-12.
bool same_motion(const Motion& a, const Motion& b) {
    if (a.kind == Motion::kOffset || b.kind == Motion::kOffset) return false;
    if ((a.shift - b.shift).Length() > 1e-12) return false;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            if (std::fabs(a.turn[r][c] - b.turn[r][c]) > 1e-12) return false;
        }
    }
    return true;
}

// One linear condition on a blend's translation v: normal . v == value.
struct Row {
    SmVector3d normal;
    double value;
};

// What a blend's translation v must satisfy to stay tangent to a neighbour
// whose surface moves by `shift`: v - shift must leave the neighbour's
// surface in place. A plane allows any move within it, a cylinder any move
// along its axis, anything else none.
void tangency_rows(const SmSurface* surface, const SmVector3d& shift, std::vector<Row>& rows) {
    SmPoint3d point;
    SmVector3d normal;
    if (plane_of(surface, point, normal)) {
        rows.push_back({normal, normal.Dot(shift)});
        return;
    }
    if (is_cylinder(surface)) {
        const SmVector3d axis = static_cast<const SmCone*>(surface)->GetPosition().GetZAxis();
        SmVector3d across = std::fabs(axis.x) < 0.9 ? SmVector3d(1.0, 0.0, 0.0) : SmVector3d(0.0, 1.0, 0.0);
        across = across - axis * across.Dot(axis);
        across = across / across.Length();
        const SmVector3d other = axis * across;  // cross product
        rows.push_back({across, across.Dot(shift)});
        rows.push_back({other, other.Dot(shift)});
        return;
    }
    for (const SmVector3d& e : {SmVector3d(1.0, 0.0, 0.0), SmVector3d(0.0, 1.0, 0.0), SmVector3d(0.0, 0.0, 1.0)}) {
        rows.push_back({e, e.Dot(shift)});
    }
}

// The shortest v meeting every row, or false when the rows conflict.
bool solve_rows(const std::vector<Row>& rows, double tolerance, SmVector3d& v) {
    std::vector<Row> basis;  // orthonormal normals, with values scaled to match
    for (const Row& row : rows) {
        SmVector3d a = row.normal;
        double b = row.value;
        for (const Row& e : basis) {
            const double k = a.Dot(e.normal);
            a = a - e.normal * k;
            b -= k * e.value;
        }
        const double length = a.Length();
        if (length < 1e-9) {
            if (std::fabs(b) > tolerance) return false;
            continue;
        }
        basis.push_back({a / length, b / length});
    }
    v = SmVector3d(0.0, 0.0, 0.0);
    for (const Row& e : basis) v = v + e.normal * e.value;
    return true;
}

// The face's outward unit normal at a point inside it.
sm_status outward_normal_inside(const SmFace* face, SmPoint3d& inside, SmVector3d& normal) {
    const SmApiStatus status = SmApiFaceInternalPoint(face, inside);
    if (status != SM_SUCCESS) return status;
    double n[3];
    if (!outward_normal(face, inside, n)) return SM_ERR;
    normal = SmVector3d(n[0], n[1], n[2]);
    return SM_SUCCESS;
}

sm_status carry_blends(const SmTArray<SmFace*>& picked, double tolerance, double slack, Motions& motions);

// Works out how every face changes when the picked faces move `distance`
// along their outward normals:
//
// - A picked plane is translated by distance x its normal, and a picked cone
//   along its axis (its offset is that translation); any other picked curved
//   face is replaced by its offset.
// - Blends joined to a moving face through smooth edges (fillets, the
//   corners between them) are translated by the shortest move that keeps
//   them tangent to every neighbour, layer by layer outward from the picked
//   faces: a fillet between a moved top and a kept wall slides along the
//   wall. A blend that the move leaves in place, such as a vertical fillet
//   under a raised top, stretches like a wall instead.
//
// Every other face keeps its surface and is extended or trimmed along it.
//
// Refused (SM_SHIM_ERR_UNSUPPORTED_FACE): blends on a picked curved face,
// whose new shape a translation can't give, and blends no translation keeps
// tangent to all their neighbours.
sm_status plan_offset(const SmTArray<SmFace*>& picked, double distance, double tolerance, Motions& motions) {
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        SmPoint3d inside;
        SmVector3d normal;
        const sm_status status = outward_normal_inside(picked[i], inside, normal);
        if (status != SM_SUCCESS) return status;
        Motion motion;
        SmPoint3d on_plane;
        SmVector3d plane_normal;
        const SmSurface* surface = picked[i]->GetSurface();
        if (plane_of(surface, on_plane, plane_normal)) {
            motion.kind = Motion::kShift;
            motion.shift = normal * distance;
        } else if (surface->GetType() == SmCone_TYPE && !static_cast<const SmCone*>(surface)->IsCylinder()) {
            // A cone's offset is the cone slid along its axis, by whatever
            // moves each point `distance` along the normal: d / (n . axis),
            // the same everywhere on it. Its apex lands on the new apex.
            const SmVector3d axis = static_cast<const SmCone*>(surface)->GetPosition().GetZAxis();
            motion.kind = Motion::kShift;
            motion.shift = axis * (distance / normal.Dot(axis));
        } else {
            motion.kind = Motion::kOffset;
        }
        motions[picked[i]] = motion;
    }
    return carry_blends(picked, tolerance, tolerance * (1.0 + std::fabs(distance)), motions);
}

// The blends joined to the picked faces through smooth edges, placed as
// plan_offset says, given the picked faces' motions; `slack` is how far a
// blend may miss a neighbour. A blend on a turned face can't follow it by
// sliding, and is refused, as on a picked curved face.
sm_status carry_blends(const SmTArray<SmFace*>& picked, double tolerance, double slack, Motions& motions) {
    std::vector<SmFace*> blends;
    std::vector<SmFace*> layer;
    for (ULONG i = 0; i < picked.GetSize(); ++i) layer.push_back(picked[i]);
    std::unordered_map<const SmFace*, bool> seen;
    for (SmFace* face : layer) seen[face] = true;
    while (!layer.empty()) {
        std::vector<SmFace*> next;
        for (SmFace* face : layer) {
            const Motion moving = motion_of(motions, face);
            for (SmFace* other : smooth_neighbours(face)) {
                if (seen[other] || !is_blend(other, tolerance)) continue;
                if (moving.kind == Motion::kShift && !moving.turned && moving.shift.Length() > 0.0 &&
                    invariant_under(other->GetSurface(), moving.shift / moving.shift.Length())) {
                    continue;
                }
                if (moving.kind == Motion::kOffset || moving.turned) return SM_SHIM_ERR_UNSUPPORTED_FACE;
                seen[other] = true;
                next.push_back(other);
            }
        }
        // Each blend from what its already-placed neighbours require; blends
        // further out aren't placed yet and are checked below.
        Motions placed;
        for (SmFace* blend : next) {
            std::vector<Row> rows;
            for (SmFace* other : smooth_neighbours(blend)) {
                const bool pending = std::find(next.begin(), next.end(), other) != next.end() ||
                                     (is_blend(other, tolerance) && !seen[other] && motions.count(other) == 0);
                if (pending) continue;
                const Motion motion = motion_of(motions, other);
                if (motion.kind == Motion::kOffset || motion.turned) return SM_SHIM_ERR_UNSUPPORTED_FACE;
                tangency_rows(other->GetSurface(), motion.shift, rows);
            }
            Motion motion;
            motion.kind = Motion::kShift;
            if (!solve_rows(rows, slack, motion.shift)) return SM_SHIM_ERR_UNSUPPORTED_FACE;
            placed[blend] = motion;
        }
        for (const auto& entry : placed) motions[entry.first] = entry.second;
        blends.insert(blends.end(), next.begin(), next.end());
        layer = next;
    }

    // Every blend must now be tangent to all its neighbours.
    for (SmFace* blend : blends) {
        const SmVector3d v = motions[blend].shift;
        for (SmFace* other : smooth_neighbours(blend)) {
            const Motion motion = motion_of(motions, other);
            if (motion.kind == Motion::kOffset || motion.turned) return SM_SHIM_ERR_UNSUPPORTED_FACE;
            std::vector<Row> rows;
            tangency_rows(other->GetSurface(), motion.shift, rows);
            for (const Row& row : rows) {
                if (std::fabs(row.normal.Dot(v) - row.value) > slack) return SM_SHIM_ERR_UNSUPPORTED_FACE;
            }
        }
    }

    return SM_SUCCESS;
}

// A copy of the face's surface translated by `shift`.
//
// A cone slid along its axis keeps its old height band, which no longer
// reaches the faces it meets: the executive would then trim it on its G1
// extension, a NURBS surface, rather than the cone. So a cone's band is
// widened by the slide on both sides.
sm_status translated_surface(const SmFace* face, const SmVector3d& shift, SmSurface*& out) {
    const SmSurface* surface = face->GetSurface();
    if (!surface) return SM_SHIM_ERR_UNSUPPORTED_FACE;
    SmSurface* moved = nullptr;
    SmStatus status = surface->Copy(*SmApiGetOrCreateContext(), moved);
    if (status != SM_SUCCESS) return status;
    SmAxis2Placement placement;
    placement.Translate(shift);
    status = moved->Transform(placement);
    if (status == SM_SUCCESS && moved->GetType() == SmCone_TYPE && !static_cast<SmCone*>(moved)->IsCylinder()) {
        SmCone* cone = static_cast<SmCone*>(moved);
        const double slide = std::fabs(shift.Dot(cone->GetPosition().GetZAxis()));
        SmExtent2d band = moved->GetSTEPUVDomain();
        double low = band.GetVMin() - 2.0 * slide;
        double high = band.GetVMax() + 2.0 * slide;
        // Never past the apex.
        SmBoolean cylinder = FALSE;
        SmPoint3d apex;
        double apex_height = 0.0;
        if (cone->GetApex(cylinder, apex, apex_height) == SM_SUCCESS) {
            if (apex_height <= band.GetVMin()) low = std::max(low, apex_height);
            if (apex_height >= band.GetVMax()) high = std::min(high, apex_height);
        }
        band.SetMinMax(band.GetUMin(), low, band.GetUMax(), high);
        status = cone->AdjustSTEPUVDomain(band);
    }
    if (status != SM_SUCCESS) {
        delete moved;
        return status;
    }
    out = moved;
    return SM_SUCCESS;
}

// A copy of the face's surface moved rigidly: translated_surface for a
// translation, and for a turn the surface turned and moved, a cone's band
// widened by `reach` on both sides (as a slid cone's is), never past its
// apex.
sm_status moved_surface(const SmFace* face, const Motion& motion, double reach, SmSurface*& out) {
    if (!motion.turned) return translated_surface(face, motion.shift, out);
    const SmSurface* surface = face->GetSurface();
    if (!surface) return SM_SHIM_ERR_UNSUPPORTED_FACE;
    SmSurface* moved = nullptr;
    SmStatus status = surface->Copy(*SmApiGetOrCreateContext(), moved);
    if (status != SM_SUCCESS) return status;
    status = moved->Transform(placement_of(motion));
    if (status == SM_SUCCESS && moved->GetType() == SmCone_TYPE && !static_cast<SmCone*>(moved)->IsCylinder()) {
        SmCone* cone = static_cast<SmCone*>(moved);
        SmExtent2d band = moved->GetSTEPUVDomain();
        double low = band.GetVMin() - 2.0 * reach;
        double high = band.GetVMax() + 2.0 * reach;
        SmBoolean cylinder = FALSE;
        SmPoint3d apex;
        double apex_height = 0.0;
        if (cone->GetApex(cylinder, apex, apex_height) == SM_SUCCESS) {
            if (apex_height <= band.GetVMin()) low = std::max(low, apex_height);
            if (apex_height >= band.GetVMax()) high = std::min(high, apex_height);
        }
        band.SetMinMax(band.GetUMin(), low, band.GetUMax(), high);
        status = cone->AdjustSTEPUVDomain(band);
    }
    if (status != SM_SUCCESS) {
        delete moved;
        return status;
    }
    out = moved;
    return SM_SUCCESS;
}

// The face's surface offset `distance` outward (SmSurface::CreateOffsetSurface,
// as the solid offset makes them: exact for analytic surfaces, so a cylinder,
// sphere or torus changes radius), checked to pass through the face's inside
// point moved along its outward normal.
sm_status offset_surface(const SmFace* face, double distance, double tolerance, SmSurface*& out) {
    SmPoint3d inside;
    SmVector3d normal;
    sm_status status = outward_normal_inside(face, inside, normal);
    if (status != SM_SUCCESS) return status;
    const SmPoint3d target = inside + normal * distance;
    const SmSurface* surface = face->GetSurface();
    const SmContext& context = *SmApiGetOrCreateContext();
    // The surface's own normal may point in or out; try outward first.
    SmPoint3d closest;
    SmTArray<SmPoint2d> uvs;
    double gap = 0.0;
    SmVector3d own;
    if (SmApiSurfaceClosestPoint(surface, inside, closest, uvs, gap) != SM_SUCCESS || uvs.GetSize() == 0 ||
        SmApiEvaluateSurfaceNormal(surface, SmVector2d(uvs[0].x, uvs[0].y), own) != SM_SUCCESS) {
        return SM_ERR;
    }
    const double sign = own.Dot(normal) >= 0.0 ? 1.0 : -1.0;
    SmSurface* offset = nullptr;
    status = static_cast<sm_status>(surface->CreateOffsetSurface(context, sign * distance, tolerance / 10.0, offset));
    if (status != SM_SUCCESS || !offset) {
        delete offset;
        return status != SM_SUCCESS ? status : SM_SHIM_ERR_UNSUPPORTED_FACE;
    }
    uvs.RemoveAll();
    if (SmApiSurfaceClosestPoint(offset, target, closest, uvs, gap) != SM_SUCCESS || gap > tolerance) {
        delete offset;
        return SM_SHIM_ERR_UNSUPPORTED_FACE;
    }
    out = offset;
    return SM_SUCCESS;
}

// Whether the face's surface, after its own motion, is also where `shift`
// would put it: the motion is `shift`, or the difference moves the surface
// within itself.
bool lands_like(const SmFace* face, const Motion& motion, const SmVector3d& shift) {
    if (motion.kind == Motion::kOffset || motion.turned) return false;
    const SmVector3d difference = shift - motion.shift;
    const double length = difference.Length();
    return length <= 1e-12 || invariant_under(face->GetSurface(), difference / length);
}

// The same for a rigid motion that may turn: a turned one lands alike only
// as the same motion.
bool lands_like(const SmFace* face, const Motion& motion, const Motion& candidate) {
    if (motion.kind == Motion::kOffset || candidate.kind == Motion::kOffset) return false;
    if (!motion.turned && !candidate.turned) return lands_like(face, motion, candidate.shift);
    return same_motion(motion, candidate);
}

// SmBrep::LocalOperation's geometry, except for edges that move rigidly.
//
// An edge moves rigidly by v when every face at both of its ends lands where
// v would put it: it moves by v, or v moves it within its own surface (a
// wall under a raised fillet, the stretched fillet beside a pushed wall).
// The new edge is then exactly the old one moved by v, ends and all, so
// nothing needs intersecting. Each face's copy of the edge (its old trim
// curve lifted onto its new surface) is right already where the face itself
// moved by v; elsewhere it is replaced with the moved edge here.
//
// Mating such edges instead re-intersects surfaces whose intersection is
// known. That loses precision, fails outright where two equal fillets meet
// at a corner (their cylinders touch tangentially at its tip), and on tangent
// edges between a moved face and a kept one lifts whichever side's old trim
// curve the executive meets first, which is the old edge half the time.
class RigidMoveGeometry : public SmLocalOperationGeometryCreation {
public:
    RigidMoveGeometry(double distance, double tolerance, SmTArray<SmFace*>& faces, SmTArray<SmSurface*>& surfaces,
                      const Motions& motions)
        : SmLocalOperationGeometryCreation(distance, tolerance, faces, surfaces), m_motions(motions) {}

    SmStatus EdgeuseOffset(const SmContext& context, SmOffsetExecutive* executive, const SmEdgeuse* original,
                           const SmEdgeuse* other_original, SmEdgeuse* offset, SmEdgeuse* radial_offset,
                           SmSurface*& new_surface, SmTArray<SmCurve*>& trims, SmTArray<SmOrientType>& orients,
                           SmSSIData& self_intersections, SmMarkType mark) override {
        if (!m_prepared) {
            m_prepared = true;
            move_rigid_edges(context, executive, original->GetFace()->GetBrep(), mark);
        }
        if (m_moved.count(original->GetEdge())) {
            new_surface = NULL;
            self_intersections.m_v3DCurves.ReSet();
            self_intersections.m_vUVCurves1.ReSet();
            self_intersections.m_vUVCurves2.ReSet();
            self_intersections.m_vCurveTypes.ReSet();
            self_intersections.m_vDeviations.ReSet();
            return SM_SUCCESS;
        }
        return SmLocalOperationGeometryCreation::EdgeuseOffset(context, executive, original, other_original, offset,
                                                               radial_offset, new_surface, trims, orients,
                                                               self_intersections, mark);
    }

private:
    // The motion the edge moves rigidly by, tried as each of its faces'
    // motions.
    bool rigid_motion(const SmEdge* edge, Motion& out) const {
        std::vector<const SmFace*> around;
        const SmVertex* ends[2] = {edge->GetStartVertex(), edge->GetOtherVertex(edge->GetStartVertex())};
        for (const SmVertex* end : ends) {
            if (!end) continue;
            SmTArray<SmFace*> faces;
            end->GetFaces(faces);
            for (ULONG i = 0; i < faces.GetSize(); ++i) around.push_back(faces[i]);
        }
        const EdgeFaces sides = faces_of(edge);
        for (const SmFace* side : sides.faces) {
            if (!side) continue;
            const Motion candidate = motion_of(m_motions, side);
            if (candidate.kind == Motion::kOffset) return false;
            bool rigid = true;
            for (const SmFace* face : around) {
                rigid = rigid && lands_like(face, motion_of(m_motions, face), candidate);
            }
            if (rigid) {
                out = candidate;
                return true;
            }
        }
        return false;
    }

    // `curve` (owned) as the copy's edge, with its trim curve on the copy's
    // (extended) surface, or false if the curve won't drop onto it whole. A
    // curve on a closed surface's seam drops as a pair, one on each side
    // (u = 0 and u = 2πr on a cylinder); each of a seam's edgeuses takes the
    // one on its original's side.
    bool dropped_copy(const SmContext& context, SmOffsetExecutive* executive, const SmEdgeuse* original,
                      SmEdgeuse* copy, SmBSplineCurve* curve, SmCurve*& out, SmBSplineCurve*& uv) const {
        out = nullptr;
        uv = nullptr;
        if (!curve) return false;
        SmSurface* extended = nullptr;
        SmExtent2d domain;
        SmTArray<SmBSplineCurve*> drops;
        double drop = 0.0;
        double deviation = 0.0;
        // The original edge may sit looser on its faces than the brep's
        // tolerance (fillet mitres do): it moves no less precise than it was.
        const double tolerance = std::max(copy->GetEdge()->GetTolerance(), original->GetEdge()->GetTolerance());
        const EdgeFaces sides = faces_of(original->GetEdge());
        const bool seam = sides.faces[0] && sides.faces[0] == sides.faces[1];
        const SmBSplineCurve* reference = original->GetUVTrimCurve();
        ULONG pick = 0;
        if (executive->GetCreateExtendedSurface(copy->GetFace()->GetSurface(), domain, extended) != SM_SUCCESS ||
            !extended ||
            extended->DropAndTrimCurve(context, domain, *curve, curve->GetNaturalInterval(), tolerance / 3.0, drop,
                                       deviation, drops, TRUE, TRUE) != SM_SUCCESS ||
            drop > tolerance || drops.GetSize() == 0 ||
            (seam ? !reference || drops.GetSize() > 2 || !on_side_of(*extended, domain, *reference, drops, pick)
                  : drops.GetSize() != 1)) {
            delete curve;
            for (ULONG i = 0; i < drops.GetSize(); ++i) delete drops[i];
            return false;
        }
        out = curve;
        uv = drops[pick];
        for (ULONG i = 0; i < drops.GetSize(); ++i) {
            if (i != pick) delete drops[i];
        }
        return true;
    }

    // Which of `drops` (one curve's UV curves over the same interval: one, or
    // a seam's pair) is on the side of `surface`'s seam that `reference` is
    // on: the one whose middle is nearest, moved by whole periods where the
    // surface is closed over `domain` and it still lies a period away.
    static bool on_side_of(const SmSurface& surface, const SmExtent2d& domain, const SmBSplineCurve& reference,
                           SmTArray<SmBSplineCurve*>& drops, ULONG& pick) {
        const SmPoint3d want = middle(reference);
        double nearest = SM_BIG_DOUBLE;
        for (ULONG i = 0; i < drops.GetSize(); ++i) {
            const SmExtent1d interval = drops[i]->GetNaturalInterval();
            if (!SM_ARE_SAME(interval.GetMin(), drops[0]->GetNaturalInterval().GetMin()) ||
                !SM_ARE_SAME(interval.GetMax(), drops[0]->GetNaturalInterval().GetMax())) {
                return false;
            }
            const double distance = (middle(*drops[i]) - want).Length();
            if (distance < nearest) {
                nearest = distance;
                pick = i;
            }
        }
        const SmPoint3d have = middle(*drops[pick]);
        SmVector3d shift(0.0, 0.0, 0.0);
        if (surface.IsClosed(domain, SM_SP_U)) {
            shift.x = std::round((want.x - have.x) / domain.XLength()) * domain.XLength();
        }
        if (surface.IsClosed(domain, SM_SP_V)) {
            shift.y = std::round((want.y - have.y) / domain.YLength()) * domain.YLength();
        }
        if (shift.x == 0.0 && shift.y == 0.0) return true;
        SmAxis2Placement placement;
        placement.Translate(shift);
        return drops[pick]->Transform(placement) == SM_SUCCESS;
    }

    static SmPoint3d middle(const SmBSplineCurve& curve) {
        SmPoint3d point;
        curve.EvaluatePoint(curve.GetNaturalInterval().Evaluate(0.5), point);
        return point;
    }

    // Puts a new curve on every copy of `edge` that `replace` picks, from
    // `make` (given the original edgeuse; the curve runs as the edge's own);
    // all or nothing. A seam's copies are both uses of one edge of the copied
    // face: the first puts its curve on that edge, the other only its own UV
    // curve, which runs as the edge's new curve (ReplaceEdgeCurves reverses
    // neither, since both run as the edge's own).
    template <typename Replace, typename Make>
    bool replace_copies(const SmContext& context, SmOffsetExecutive* executive, SmEdge* edge, Replace&& replace,
                        Make&& make) {
        SmTArray<SmEdgeuse*> uses;
        edge->GetEdgeuses(uses);
        std::vector<std::pair<SmEdgeuse*, SmEdgeuse*>> sides;  // (original, its copy)
        for (ULONG i = 0; i < uses.GetSize(); ++i) {
            SmEdgeuse* copy = executive->GetOffsetEdgeuse(uses[i]);
            if (copy) sides.emplace_back(uses[i], copy);
        }
        if (sides.size() < 2) return false;
        std::vector<SmCurve*> curves(sides.size(), nullptr);
        std::vector<SmBSplineCurve*> uvs(sides.size(), nullptr);
        bool ready = true;
        for (size_t i = 0; i < sides.size() && ready; ++i) {
            if (!replace(sides[i].first->GetFace())) continue;
            ready = dropped_copy(context, executive, sides[i].first, sides[i].second, make(sides[i].first), curves[i],
                                 uvs[i]);
        }
        // Whether each copy's edge took its new curve, once one copy has.
        std::unordered_map<const SmEdge*, bool> placed;
        for (size_t i = 0; i < sides.size() && ready; ++i) {
            if (!curves[i]) continue;
            SmEdgeuse* copy = sides[i].second;
            const auto shared = placed.find(copy->GetEdge());
            if (shared != placed.end()) {
                // ReplaceEdgeCurves would double the edge's tolerance again.
                if (shared->second) {
                    copy->SetUVTrimCurve(uvs[i], 0.0, TRUE);
                    uvs[i] = nullptr;
                }
                continue;
            }
            SmBoolean delete_curve = FALSE;
            SmBoolean delete_uv = FALSE;
            const SmStatus status = ReplaceEdgeCurves(copy, curves[i], uvs[i], sides[i].first, delete_curve, delete_uv);
            placed.emplace(copy->GetEdge(), !delete_curve);
            if (delete_curve) delete curves[i];
            if (delete_uv) delete uvs[i];
            curves[i] = nullptr;
            uvs[i] = nullptr;
            ready = status == SM_SUCCESS;
        }
        for (size_t i = 0; i < sides.size(); ++i) {
            delete curves[i];
            delete uvs[i];
        }
        return ready;
    }

    // Puts the original edge moved by `motion` on every copy of `edge` whose
    // face didn't itself move by `motion`.
    bool move_copies(const SmContext& context, SmOffsetExecutive* executive, SmEdge* edge, const Motion& motion) {
        const SmAxis2Placement placement = placement_of(motion);
        return replace_copies(
            context, executive, edge,
            [&](const SmFace* face) {
                const Motion own = motion_of(m_motions, face);
                if (!own.turned && !motion.turned) return (own.shift - motion.shift).Length() > 1e-12;
                return !same_motion(own, motion);
            },
            [&](const SmEdgeuse* original) {
                SmBSplineCurve* moved = original->GetEdge()->CreateTrimmedNURBSCurve(context);
                if (moved && moved->Transform(placement) != SM_SUCCESS) {
                    delete moved;
                    moved = nullptr;
                }
                return moved;
            });
    }

    // Where a straight edge goes when it stays on its own line and only its
    // ends move along it: a wall's edge with a fillet under a raised top, or
    // beside a pushed wall, or a hole's seam under a moved floor. Both its
    // faces land as no move would put them (kept, or moved within their own
    // surfaces), and both ends are placed (aim_vertices) on the line, in the
    // same order. Mating such an edge lifts the old trim curve of whichever
    // face (or side of the seam) the executive meets first, which is off in
    // the last digits.
    bool stretched_line(const SmContext& context, const SmEdge* edge, SmPoint3d& from, SmPoint3d& to) const {
        SmVector3d direction;
        if (!straight(edge, direction)) return false;
        const EdgeFaces sides = faces_of(edge);
        if (!sides.faces[0] || !sides.faces[1]) return false;
        const SmVector3d none(0.0, 0.0, 0.0);
        for (const SmFace* side : sides.faces) {
            if (!lands_like(side, motion_of(m_motions, side), none)) return false;
        }
        const SmVertex* ends[2] = {edge->GetStartVertex(), edge->GetOtherVertex(edge->GetStartVertex())};
        if (!ends[0] || !ends[1] || ends[0] == ends[1] || !m_placed.count(ends[0]) || !m_placed.count(ends[1])) {
            return false;
        }
        // The edge's curve as the edge runs, against where its ends were.
        const std::unique_ptr<SmBSplineCurve> curve(edge->CreateTrimmedNURBSCurve(context));
        if (!curve) return false;
        SmPoint3d start;
        SmPoint3d end;
        curve->GetEnds(start, end);
        const SmPoint3d old0 = old_point(ends[0]);
        const SmPoint3d old1 = old_point(ends[1]);
        const bool forward = (start - old0).Length() + (end - old1).Length() <= (start - old1).Length() + (end - old0).Length();
        from = (forward ? ends[0] : ends[1])->GetPoint();
        to = (forward ? ends[1] : ends[0])->GetPoint();
        const SmVector3d before = end - start;
        const SmVector3d after = to - from;
        if (after.Dot(before) <= 0.0 || after.Length() <= edge->GetTolerance()) return false;
        // Both new ends on the old line.
        for (const SmPoint3d& point : {from, to}) {
            const SmVector3d offset = point - start;
            if ((offset - direction * offset.Dot(direction)).Length() > edge->GetTolerance()) return false;
        }
        return true;
    }

    // Puts the stretched line (see stretched_line) on both copies of `edge`.
    bool stretch_copies(const SmContext& context, SmOffsetExecutive* executive, SmEdge* edge) {
        SmPoint3d from;
        SmPoint3d to;
        if (!stretched_line(context, edge, from, to)) return false;
        return replace_copies(
            context, executive, edge, [](const SmFace*) { return true; },
            [&](const SmEdgeuse*) {
                SmBSplineCurve* line = nullptr;
                if (SmBSplineCurve::CreateLineSegment(context, 3, from, to, line) != SM_SUCCESS) {
                    delete line;
                    line = nullptr;
                }
                return line;
            });
    }

    SmPoint3d old_point(const SmVertex* vertex) const {
        for (const auto& aimed : m_aimed) {
            if (aimed.first == vertex) return aimed.second;
        }
        return vertex->GetPoint();
    }

    // Before the executive mates any edge, moves every rigid one and marks its
    // edgeuses as mated. Its neighbours then find its moved ends whatever
    // order the executive mates them in: a tangent neighbour that
    // isn't marked yet is trimmed against its old trim curve instead, which
    // puts it back at the old corner.
    void move_rigid_edges(const SmContext& context, SmOffsetExecutive* executive, SmBrep* brep, SmMarkType mark) {
        if (!brep) return;
        aim_vertices(*brep);
        SmTArray<SmEdge*> edges;
        brep->GetEdges(edges);
        for (ULONG i = 0; i < edges.GetSize(); ++i) {
            Motion motion;
            const bool rigid = rigid_motion(edges[i], motion);
            if (rigid ? !move_copies(context, executive, edges[i], motion) : !stretch_copies(context, executive, edges[i])) {
                continue;
            }
            m_moved.insert(edges[i]);
            SmTArray<SmEdgeuse*> uses;
            edges[i]->GetEdgeuses(uses);
            for (ULONG j = 0; j < uses.GetSize(); ++j) uses[j]->Mark(mark);
        }
    }

    // Where an edge's end goes is chosen among candidate points by nearness to
    // the original vertex (TrimExtendedEdge), which after a long move can be
    // the wrong one: a stretched fillet's end 10 mm from a pushed wall picked
    // the far side of the moved fillet's cylinder. So, once the executive has
    // built its faces (it reads the original brep's vertices only as such
    // references from then on), each vertex that every face around it would
    // carry the same way is moved to where it will land, until restore().
    void aim_vertices(SmBrep& brep) {
        SmTArray<SmVertex*> vertices;
        brep.GetVertices(vertices);
        for (ULONG i = 0; i < vertices.GetSize(); ++i) {
            SmPoint3d aimed;
            if (!predict_vertex(vertices[i], aimed)) continue;
            m_placed.insert(vertices[i]);
            if ((aimed - vertices[i]->GetPoint()).Length() <= 1e-12) continue;
            m_aimed.emplace_back(vertices[i], vertices[i]->GetPoint());
            vertices[i]->SetPoint(aimed);
        }
    }

    // Where a vertex will land. When every face around it lands as one shift
    // would put it, that shift; otherwise the point the moved planes and
    // straight edges around it agree on: each planar face's plane after its
    // motion, and each straight edge whose two faces both land as its moved
    // copy would (a tangent line under a sliding fillet, a wall's edge left
    // in place), as a line. False when they don't pin one point.
    bool predict_vertex(const SmVertex* vertex, SmPoint3d& out) const {
        const SmPoint3d old = vertex->GetPoint();
        SmTArray<SmFace*> around;
        vertex->GetFaces(around);
        for (ULONG j = 0; j < around.GetSize(); ++j) {
            if (motion_of(m_motions, around[j]).kind == Motion::kOffset) return false;
        }
        for (ULONG j = 0; j < around.GetSize(); ++j) {
            const Motion candidate = motion_of(m_motions, around[j]);
            bool lands = true;
            for (ULONG k = 0; k < around.GetSize() && lands; ++k) {
                lands = lands_like(around[k], motion_of(m_motions, around[k]), candidate);
            }
            if (lands) {
                out = moved_point(candidate, old);
                return true;
            }
        }
        std::vector<Row> rows;
        for (ULONG j = 0; j < around.GetSize(); ++j) {
            SmPoint3d point;
            SmVector3d normal;
            if (!plane_of(around[j]->GetSurface(), point, normal)) continue;
            const Motion motion = motion_of(m_motions, around[j]);
            if (!motion.turned) {
                rows.push_back({normal, normal.Dot(point + motion.shift)});
                continue;
            }
            const SmVector3d turned = turned_vector(motion, normal);
            rows.push_back({turned, turned.Dot(moved_point(motion, point))});
        }
        SmTArray<SmEdge*> edges;
        vertex->GetEdges(edges);
        for (ULONG j = 0; j < edges.GetSize(); ++j) {
            SmVector3d direction;
            if (!straight(edges[j], direction)) continue;
            const EdgeFaces sides = faces_of(edges[j]);
            if (!sides.faces[0] || !sides.faces[1]) continue;
            for (const SmFace* side : sides.faces) {
                const Motion candidate = motion_of(m_motions, side);
                if (!lands_like(sides.faces[0], motion_of(m_motions, sides.faces[0]), candidate) ||
                    !lands_like(sides.faces[1], motion_of(m_motions, sides.faces[1]), candidate)) {
                    continue;
                }
                const SmVector3d along = turned_vector(candidate, direction);
                SmVector3d across = std::fabs(along.x) < 0.9 ? SmVector3d(1.0, 0.0, 0.0) : SmVector3d(0.0, 1.0, 0.0);
                across = across - along * across.Dot(along);
                across = across / across.Length();
                const SmVector3d other = along * across;  // cross product
                const SmPoint3d through = moved_point(candidate, old);
                rows.push_back({across, across.Dot(through)});
                rows.push_back({other, other.Dot(through)});
                break;
            }
        }
        // Least squares: the normal equations, solved by Cramer's rule.
        double m[3][3] = {{0.0}};
        double b[3] = {0.0};
        for (const Row& row : rows) {
            const double n[3] = {row.normal.x, row.normal.y, row.normal.z};
            for (int r = 0; r < 3; ++r) {
                b[r] += n[r] * row.value;
                for (int c = 0; c < 3; ++c) m[r][c] += n[r] * n[c];
            }
        }
        auto det = [](const double a[3][3]) {
            return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                   a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
        };
        const double d = det(m);
        if (std::fabs(d) < 1e-6) return false;
        double p[3];
        for (int c = 0; c < 3; ++c) {
            double a[3][3];
            for (int r = 0; r < 3; ++r) {
                for (int k = 0; k < 3; ++k) a[r][k] = k == c ? b[r] : m[r][k];
            }
            p[c] = det(a) / d;
        }
        const SmPoint3d solved(p[0], p[1], p[2]);
        for (const Row& row : rows) {
            if (std::fabs(row.normal.Dot(solved) - row.value) > 1e-6 * (1.0 + std::fabs(row.value))) return false;
        }
        out = solved;
        return true;
    }

    // Whether the edge is a straight segment, and its direction.
    static bool straight(const SmEdge* edge, SmVector3d& direction) {
        double t0 = 0.0;
        double t1 = 0.0;
        if (SmApiEdgeParameterRange(edge, t0, t1) != SM_SUCCESS) return false;
        SmVector3d points[3];
        for (int k = 0; k < 3; ++k) {
            SmVector3d* point_out = &points[k];
            SmVector3d* no_derivative = nullptr;
            SmVector3d* no_second = nullptr;
            if (SmApiEvaluateCurve(edge->GetCurve(), t0 + (t1 - t0) * k / 2.0, point_out, no_derivative, no_second) !=
                SM_SUCCESS) {
                return false;
            }
        }
        const SmVector3d span = points[2] - points[0];
        const double length = span.Length();
        if (length <= 0.0) return false;
        direction = span / length;
        const SmVector3d off = points[1] - points[0];
        return (off - direction * off.Dot(direction)).Length() <= 1e-9 * (1.0 + length);
    }

public:
    // Puts the original brep's vertices back.
    void restore() {
        for (auto& aimed : m_aimed) aimed.first->SetPoint(aimed.second);
        m_aimed.clear();
        m_placed.clear();
    }

private:
    bool m_prepared = false;
    std::unordered_set<const SmEdge*> m_moved;
    std::vector<std::pair<SmVertex*, SmPoint3d>> m_aimed;
    std::unordered_set<const SmVertex*> m_placed;  // every vertex aim_vertices placed, moved or not
    const Motions& m_motions;
};

// SmBrep::LocalOperation with RigidMoveGeometry: gives each face in `moved`
// its surface in `surfaces`, which it takes, and every other face a copy of
// its own, then rebuilds the brep from the offset executive's result. The
// brep is untouched unless this succeeds.
// Whether every edge of the face lies on the surface, sampled, and none is a
// seam (used twice by the face), whose trim curves rebuilding by projection
// would put on one side of the seam.
bool edges_lie_on(const SmFace* face, const SmSurface* surface, double tolerance) {
    SmTArray<SmEdge*> edges;
    face->GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        const EdgeFaces sides = faces_of(edges[i]);
        if (sides.faces[0] == sides.faces[1]) return false;
        double t0 = 0.0;
        double t1 = 0.0;
        if (SmApiEdgeParameterRange(edges[i], t0, t1) != SM_SUCCESS) return false;
        for (int k = 0; k <= 4; ++k) {
            SmVector3d point;
            SmVector3d* point_out = &point;
            SmVector3d* no_derivative = nullptr;
            SmVector3d* no_second = nullptr;
            if (SmApiEvaluateCurve(edges[i]->GetCurve(), t0 + (t1 - t0) * k / 4.0, point_out, no_derivative,
                                   no_second) != SM_SUCCESS) {
                return false;
            }
            SmPoint3d closest;
            SmTArray<SmPoint2d> uvs;
            double gap = 0.0;
            if (SmApiSurfaceClosestPoint(surface, SmPoint3d(point.x, point.y, point.z), closest, uvs, gap) !=
                    SM_SUCCESS ||
                gap > tolerance) {
                return false;
            }
        }
    }
    return true;
}

// Whether the surface is one of the kinds plan_offset recognises: SMLib's
// IsAnalytic also counts any surface of extrusion or revolution.
bool canonical(const SmSurface* surface) {
    const SM_TYPE type = surface->GetType();
    return type == SmPlane_TYPE || type == SmCylinder_TYPE || type == SmCone_TYPE || type == SmSphere_TYPE ||
           type == SmTorus_TYPE;
}

// After the executive, a face it left on an extension of one of the
// `analytic` surfaces gets that surface back. The executive only keeps the
// original where the new trims fit its domain to SM_EFF_ZERO, and an offset
// cone's can spill past it by a hair; it stretches a cylinder along its axis
// by making it a surface of extrusion. The trims are rebuilt from the edges,
// so only where every edge lies on the analytic surface and none is a seam;
// elsewhere the face keeps the extension, which has the same shape.
void restore_analytic_surfaces(SmBrep& brep, const std::vector<const SmSurface*>& analytic) {
    if (analytic.empty()) return;
    const SmContext& context = *brep.GetContext();
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const SmSurface* current = faces[i]->GetSurface();
        if (!current || canonical(current)) continue;
        SmPoint3d inside;
        SmVector3d normal;
        if (outward_normal_inside(faces[i], inside, normal) != SM_SUCCESS) continue;
        for (const SmSurface* candidate : analytic) {
            SmPoint3d closest;
            SmTArray<SmPoint2d> uvs;
            double gap = 0.0;
            if (SmApiSurfaceClosestPoint(candidate, inside, closest, uvs, gap) != SM_SUCCESS || uvs.GetSize() == 0 ||
                gap > brep.GetTolerance()) {
                continue;
            }
            SmVector3d candidate_normal;
            SmVector3d current_normal;
            SmTArray<SmPoint2d> current_uvs;
            double current_gap = 0.0;
            if (SmApiEvaluateSurfaceNormal(candidate, SmVector2d(uvs[0].x, uvs[0].y), candidate_normal) != SM_SUCCESS ||
                SmApiSurfaceClosestPoint(current, inside, closest, current_uvs, current_gap) != SM_SUCCESS ||
                current_uvs.GetSize() == 0 ||
                SmApiEvaluateSurfaceNormal(current, SmVector2d(current_uvs[0].x, current_uvs[0].y), current_normal) !=
                    SM_SUCCESS ||
                candidate_normal.Dot(current_normal) < 0.999) {
                continue;
            }
            if (!edges_lie_on(faces[i], candidate, brep.GetTolerance())) continue;
            SmSurface* copy = nullptr;
            if (candidate->Copy(context, copy) != SM_SUCCESS || !copy) continue;
            if (brep.ReplaceSurface(faces[i], copy, FALSE, NULL, TRUE) != SM_SUCCESS) continue;
            break;
        }
    }
}

sm_status rigid_local_operation(SmBrep& brep, const SmTArray<SmFace*>& moved, SmTArray<SmSurface*>& surfaces,
                                const Motions& motions, double reach) {
    // Copies of the analytic offsets, to restore afterwards.
    std::vector<std::unique_ptr<SmSurface>> owned;
    std::vector<const SmSurface*> analytic;
    for (ULONG i = 0; i < moved.GetSize(); ++i) {
        if (motion_of(motions, moved[i]).kind != Motion::kOffset || !surfaces[i]->IsAnalytic()) continue;
        SmSurface* copy = nullptr;
        if (surfaces[i]->Copy(*brep.GetContext(), copy) == SM_SUCCESS && copy) {
            owned.emplace_back(copy);
            analytic.push_back(copy);
        }
    }
    // And of the kept cylinders, longer by twice the reach at each end, for
    // the ones the executive stretches along their axis (the fillet beside a
    // pushed wall): left as surfaces of extrusion, the next offset wouldn't
    // take them for blends.
    {
        SmTArray<SmFace*> kept;
        brep.GetFaces(kept);
        for (ULONG i = 0; i < kept.GetSize(); ++i) {
            ULONG index = 0;
            const SmSurface* surface = kept[i]->GetSurface();
            if (moved.FindElement(kept[i], index) || !surface || !is_cylinder(surface)) continue;
            SmSurface* copy = nullptr;
            if (surface->Copy(*brep.GetContext(), copy) != SM_SUCCESS || !copy) continue;
            std::unique_ptr<SmSurface> longer(copy);
            SmExtent2d band = copy->GetSTEPUVDomain();
            band.SetMinMax(band.GetUMin(), band.GetVMin() - 2.0 * reach, band.GetUMax(), band.GetVMax() + 2.0 * reach);
            if (static_cast<SmCone*>(copy)->AdjustSTEPUVDomain(band) != SM_SUCCESS) continue;
            analytic.push_back(copy);
            owned.push_back(std::move(longer));
        }
    }
    SmTemporaryChangeValue<SmBoolean> editing(brep.m_bEditingEnabled, TRUE);
    const SmContext* context = brep.GetContext();
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    SmTArray<SmSurface*> all_surfaces;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        ULONG index = 0;
        if (moved.FindElement(faces[i], index)) {
            all_surfaces.Add(surfaces[index]);
            continue;
        }
        SmSurface* copy = nullptr;
        const SmStatus status = faces[i]->GetSurface()->CopyAndAddAnalytics(*context, copy);
        if (status != SM_SUCCESS || !copy) {
            for (ULONG j = 0; j < all_surfaces.GetSize(); ++j) {
                if (!surfaces.FindElement(all_surfaces[j], index)) delete all_surfaces[j];
            }
            for (ULONG j = 0; j < surfaces.GetSize(); ++j) delete surfaces[j];
            return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
        }
        all_surfaces.Add(copy);
    }

    SmExtent3d box;
    brep.CalculateBoundingBox(box);
    SmPoint3d centre;
    double radius = 0.0;
    box.ComputeSphereBound(centre, radius);
    // The distance only sizes the surface extensions, at five times it:
    // LocalOperation's 2.5 x the bounding radius, or far enough for the move.
    RigidMoveGeometry geometry(std::max(radius / 2.0, reach), brep.GetTolerance() / 2.0, faces, all_surfaces,
                               motions);
    SmOffsetExecutive executive(*context, SM_OO_LOCAL_OPERATION, geometry, &brep);
    SmTArray<SmFace*> no_shell_faces;
    executive.SetShellFaces(no_shell_faces);
    executive.SetExtendConvexEdges(TRUE);
    SmBrep* result = nullptr;
    // The executive deletes its working brep, on failure too. On failure the
    // surfaces it hadn't yet given a face leak, as in LocalOperation.
    const SmStatus status = executive.DoExtendedOffset(result);
    geometry.restore();
    if (status != SM_SUCCESS || !result) {
        delete result;
        return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
    }
    const std::unique_ptr<SmBrep> rebuilt(result);
    SmTArray<SmFace*> old_faces;
    brep.GetFaces(old_faces);
    SmStatus moved_status = brep.RemoveFaces(old_faces);
    if (moved_status != SM_SUCCESS) return static_cast<sm_status>(moved_status);
    SmTArray<SmFace*> new_faces;
    rebuilt->GetFaces(new_faces);
    moved_status = rebuilt->CopyFaces(new_faces, &brep);
    if (moved_status != SM_SUCCESS) return static_cast<sm_status>(moved_status);
    restore_analytic_surfaces(brep, analytic);
    return SM_SUCCESS;
}

// Moves the picked faces `distance` along their outward normals (see
// plan_offset), re-intersecting every other neighbour along its own surface.
sm_status offset_faces(SmBrep& brep, const SmTArray<SmFace*>& picked, double distance) {
    Motions motions;
    sm_status result = plan_offset(picked, distance, brep.GetTolerance(), motions);
    if (result != SM_SUCCESS) return result;

    // In the brep's face order, not the map's (by address): the analytic
    // surfaces restored afterwards are tried in this order.
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    SmTArray<SmFace*> changed;
    SmTArray<SmSurface*> surfaces;
    double reach = std::fabs(distance);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        SmFace* face = faces[i];
        const Motion motion = motion_of(motions, face);
        SmSurface* surface = nullptr;
        if (motion.kind == Motion::kShift) {
            result = translated_surface(face, motion.shift, surface);
            reach = std::max(reach, motion.shift.Length());
        } else if (motion.kind == Motion::kOffset) {
            result = offset_surface(face, distance, brep.GetTolerance(), surface);
        } else {
            continue;
        }
        if (result != SM_SUCCESS) {
            for (ULONG j = 0; j < surfaces.GetSize(); ++j) delete surfaces[j];
            return result;
        }
        changed.Add(face);
        surfaces.Add(surface);
    }
    return rigid_local_operation(brep, changed, surfaces, motions, reach);
}

// Moving faces rigidly (sm_brep_move_faces).
//
// The picked faces all take one rigid motion: their surfaces are moved, and
// turned when it turns, and every other neighbour is extended or trimmed
// along its own surface to meet them, as an offset's are, through the same
// local operation. A translation carries the blends on them as an offset
// does (carry_blends): a round between a moved face and a kept one slides to
// stay tangent to both, and one the move leaves in place stretches. A blend
// on a turned face can't follow it by sliding and is refused.

// The motion every face takes when the picked ones take `rigid`, whose
// translation is `translation`: `rigid` for them, and for the blends
// carry_blends places their slides. A picked plane takes only the part of
// the translation along its normal, after the turn: moved within itself a
// plane is where it was, and its parametrisation shifting under the face's
// trims would carry them off its edges (a lid moved up and sideways
// measured 7e-5 more than moved up). A curved face keeps the whole
// translation: a cylinder held back along its axis would have to be
// extended, which the executive does with a NURBS surface (a boss lifted
// with its top came out 6e-6 over).
sm_status plan_move(const SmTArray<SmFace*>& picked, const Motion& rigid, const SmVector3d& translation, double tolerance,
                    Motions& motions) {
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        Motion motion = rigid;
        SmPoint3d point;
        SmVector3d normal;
        if (plane_of(picked[i]->GetSurface(), point, normal) && normal.Length() > 0.0) {
            const SmVector3d n = turned_vector(rigid, normal / normal.Length());
            motion.shift = motion.shift - (translation - n * translation.Dot(n));
        }
        motions[picked[i]] = motion;
    }
    return carry_blends(picked, tolerance, tolerance * (1.0 + rigid.shift.Length()), motions);
}

// How far the motion takes any point of the face: the most it takes a corner
// of the face's bounding box, which bounds it (the distance moved is convex
// in the point).
double displacement(const SmFace* face, const Motion& motion) {
    if (!motion.turned) return motion.shift.Length();
    SmExtent3d box;
    if (face->CalculateBoundingBox(box) != SM_SUCCESS) return motion.shift.Length();
    const SmPoint3d low = box.GetMin();
    const SmPoint3d high = box.GetMax();
    double most = 0.0;
    for (int corner = 0; corner < 8; ++corner) {
        const SmPoint3d p((corner & 1) ? high.x : low.x, (corner & 2) ? high.y : low.y, (corner & 4) ? high.z : low.z);
        most = std::max(most, (moved_point(motion, p) - p).Length());
    }
    return most;
}

// Moves every face plan_move moved, re-intersecting every other neighbour
// along its own surface.
sm_status move_planned(SmBrep& brep, const Motions& motions) {
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    double reach = 0.0;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const Motion motion = motion_of(motions, faces[i]);
        if (motion.kind == Motion::kShift) reach = std::max(reach, displacement(faces[i], motion));
    }
    SmTArray<SmFace*> changed;
    SmTArray<SmSurface*> surfaces;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const Motion motion = motion_of(motions, faces[i]);
        if (motion.kind != Motion::kShift) continue;
        SmSurface* surface = nullptr;
        const sm_status result = moved_surface(faces[i], motion, reach, surface);
        if (result != SM_SUCCESS) {
            for (ULONG j = 0; j < surfaces.GetSize(); ++j) delete surfaces[j];
            return result;
        }
        changed.Add(faces[i]);
        surfaces.Add(surface);
    }
    return rigid_local_operation(brep, changed, surfaces, motions, reach);
}

// What each face of the brep is to a move planned as `motions`
// (sm_move_role): moved, carried (a blend that slides with them), followed
// (sharing a vertex with a face that moves, so its boundary changes) or kept.
std::vector<uint8_t> move_roles(const SmBrep& brep, const SmTArray<SmFace*>& picked, const Motions& motions) {
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    std::unordered_set<const SmFace*> moving;
    std::unordered_set<const SmFace*> touched;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        if (motion_of(motions, faces[i]).kind == Motion::kKept) continue;
        moving.insert(faces[i]);
        SmTArray<SmVertex*> vertices;
        faces[i]->GetVertices(vertices);
        for (ULONG v = 0; v < vertices.GetSize(); ++v) {
            SmTArray<SmFace*> around;
            vertices[v]->GetFaces(around);
            for (ULONG a = 0; a < around.GetSize(); ++a) touched.insert(around[a]);
        }
    }
    std::vector<uint8_t> roles(faces.GetSize(), SM_MOVE_KEPT);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        ULONG index = 0;
        if (picked.FindElement(faces[i], index)) {
            roles[i] = SM_MOVE_MOVED;
        } else if (moving.count(faces[i])) {
            roles[i] = SM_MOVE_CARRIED;
        } else if (touched.count(faces[i])) {
            roles[i] = SM_MOVE_FOLLOWED;
        }
    }
    return roles;
}

// Moving a face through what stands in its way.
//
// A face pushed in sweeps the space between its old plane and its new one.
// Its neighbours trim back with it, and a feature standing on one of them
// near the moving edge (a boss near a lid's edge when the wall below is
// pushed in) is then crossed by the neighbour's new boundary: new topology,
// which a local operation can't make. So is a feature under a blend that
// moves with the face (a pocket in a wall that a corner round slides over),
// and one on the face itself that the move swallows a face of (a pocket
// pushed past its floor, a chamfer pushed past its width, a boss a pulled
// face rises past). A curved face moved either way crosses features beside
// it the same way: a boss on a lid near a hole that widens, a hole beside a
// boss that grows, a pin on a post that thins. Such features are set aside
// and the face is offset without them. A boss comes back as its material
// behind the face's new plane, or off what a curved face swept, joined back;
// or not at all when it's under it. Holes, filled for the offset, come back
// by taking what the offset swept away from the part with its holes (or
// adding what it grew), which cuts a hole the face runs through along its new
// surface.

// A feature the offset runs into: faces standing on `base`, attached to it
// along inner loops, that the move crosses or swallows a face of.
struct Obstacle {
    SmFace* base = nullptr;
    // Whether `base` is one of the faces moved.
    bool on_picked = false;
    std::vector<SmFace*> faces;
    // +1: on the base's outer side (a boss); -1: inside it (a hole); 0:
    // neither, or the base isn't planar.
    int side = 0;
    // Its bounds in the base's frame (u, v, out of the solid), its box.
    SmPoint3d origin;
    SmVector3d u, v, w;
    double low[3] = {0.0, 0.0, 0.0};
    double high[3] = {0.0, 0.0, 0.0};
};

// Two unit axes along a plane with unit normal `w`, making (u, v, w)
// right-handed.
void plane_axes(const SmVector3d& w, SmVector3d& u, SmVector3d& v) {
    u = std::fabs(w.x) < 0.9 ? SmVector3d(1.0, 0.0, 0.0) : SmVector3d(0.0, 1.0, 0.0);
    u = u - w * u.Dot(w);
    u = u / u.Length();
    v = w * u;  // cross product
}

// A picked plane pushed in: a point on it, its outward normal, and how far;
// and the box each blend moving with it sweeps through (see add_blends).
struct Push {
    SmPoint3d point;
    SmVector3d normal;
    double distance;
    // A box in the frame (u, v, normal) at `point`.
    struct Box {
        double low[3];
        double high[3];
    };
    SmVector3d u, v;
    std::vector<Box> blends;
    // Whether `p` lies in the space the face sweeps: behind its old plane and
    // beyond its new one.
    bool sweeps(const SmPoint3d& p, double tolerance) const {
        const double behind_old = (p - point).Dot(normal);
        const double beyond_new = behind_old + distance;
        return behind_old <= tolerance && beyond_new > tolerance;
    }
    // Whether `p` lies in a box a blend moving with the face sweeps.
    bool blend_sweeps(const SmPoint3d& p, double tolerance) const {
        const SmVector3d d = p - point;
        const double c[3] = {d.Dot(u), d.Dot(v), d.Dot(normal)};
        for (const Box& box : blends) {
            bool inside = true;
            for (int k = 0; k < 3; ++k) {
                inside = inside && c[k] >= box.low[k] - tolerance && c[k] <= box.high[k] + tolerance;
            }
            if (inside) return true;
        }
        return false;
    }
};

std::vector<Push> pushes_of(const SmTArray<SmFace*>& picked, double distance) {
    std::vector<Push> pushes;
    if (!(distance < 0.0)) return pushes;
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        SmPoint3d on_plane;
        SmVector3d plane_normal;
        SmPoint3d inside;
        SmVector3d normal;
        if (!picked[i]->GetSurface() || !plane_of(picked[i]->GetSurface(), on_plane, plane_normal) ||
            outward_normal_inside(picked[i], inside, normal) != SM_SUCCESS) {
            continue;
        }
        Push push{inside, normal, -distance, {}, {}, {}};
        plane_axes(normal, push.u, push.v);
        pushes.push_back(push);
    }
    return pushes;
}

// A picked curved face moved `distance` along its normal, in or out: the
// space between its surface and its offset, near the face.
struct Shell {
    const SmFace* face;
    double distance;
    // The face's bounding box grown by the distance, past which nothing is
    // in the way.
    SmPoint3d low, high;
    // Whether `p` lies between the surface and its offset: its height off
    // the surface, along the face's outward normal at the nearest point,
    // runs from 0 to `distance`.
    bool sweeps(const SmPoint3d& p, double tolerance) const {
        if (p.x < low.x || p.y < low.y || p.z < low.z || p.x > high.x || p.y > high.y || p.z > high.z) return false;
        SmPoint3d closest;
        SmTArray<SmPoint2d> uvs;
        double gap = 0.0;
        double n[3];
        if (SmApiSurfaceClosestPoint(face->GetSurface(), p, closest, uvs, gap) != SM_SUCCESS ||
            !outward_normal(face, closest, n)) {
            return false;
        }
        const double height = (p - closest).Dot(SmVector3d(n[0], n[1], n[2]));
        return height >= std::min(distance, 0.0) - tolerance && height <= std::max(distance, 0.0) + tolerance;
    }
};

std::vector<Shell> shells_of(const SmTArray<SmFace*>& picked, double distance) {
    std::vector<Shell> shells;
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        SmPoint3d on_plane;
        SmVector3d plane_normal;
        SmExtent3d box;
        if (!picked[i]->GetSurface() || plane_of(picked[i]->GetSurface(), on_plane, plane_normal) ||
            picked[i]->CalculateBoundingBox(box) != SM_SUCCESS) {
            continue;
        }
        const double grow = std::fabs(distance) + 1.0;
        const SmVector3d by(grow, grow, grow);
        shells.push_back({picked[i], distance, box.GetMin() - by, box.GetMax() + by});
    }
    return shells;
}

// A frame on a face's plane: origin, two axes in it and the outward normal.
bool frame_of(const SmFace* face, SmPoint3d& origin, SmVector3d& u, SmVector3d& v, SmVector3d& w) {
    SmPoint3d on_plane;
    SmVector3d plane_normal;
    if (!face->GetSurface() || !plane_of(face->GetSurface(), on_plane, plane_normal) ||
        outward_normal_inside(face, origin, w) != SM_SUCCESS) {
        return false;
    }
    plane_axes(w, u, v);
    return true;
}

// Points that bound a face: along its edges, and its bounding box's corners
// (for faces that bulge between their edges).
void face_points(const SmFace* face, std::vector<SmPoint3d>& out) {
    SmTArray<SmEdge*> edges;
    face->GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) edge_points(edges[i], 24, out);
}

void box_corners(const SmFace* face, std::vector<SmPoint3d>& out) {
    SmExtent3d box;
    if (face->CalculateBoundingBox(box) != SM_SUCCESS) return;
    const SmPoint3d low = box.GetMin();
    const SmPoint3d high = box.GetMax();
    for (int k = 0; k < 8; ++k) {
        out.emplace_back(k & 1 ? high.x : low.x, k & 2 ? high.y : low.y, k & 4 ? high.z : low.z);
    }
}

// Each push's blends: the faces moving with the picked ones (plan_offset),
// each boxed in the push's frame from where it is to where it goes.
void add_blends(const SmBrep& brep, const SmTArray<SmFace*>& picked, const Motions& motions,
                std::vector<Push>& pushes) {
    if (pushes.empty()) return;
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        ULONG index = 0;
        const Motion motion = motion_of(motions, faces[i]);
        if (picked.FindElement(faces[i], index) || motion.kind != Motion::kShift) continue;
        std::vector<SmPoint3d> points;
        face_points(faces[i], points);
        box_corners(faces[i], points);
        for (Push& push : pushes) {
            Push::Box box = {{HUGE_VAL, HUGE_VAL, HUGE_VAL}, {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL}};
            const double shift[3] = {motion.shift.Dot(push.u), motion.shift.Dot(push.v), motion.shift.Dot(push.normal)};
            for (const SmPoint3d& p : points) {
                const SmVector3d d = p - push.point;
                const double c[3] = {d.Dot(push.u), d.Dot(push.v), d.Dot(push.normal)};
                for (int k = 0; k < 3; ++k) {
                    box.low[k] = std::min(box.low[k], c[k] + std::min(shift[k], 0.0));
                    box.high[k] = std::max(box.high[k], c[k] + std::max(shift[k], 0.0));
                }
            }
            push.blends.push_back(box);
        }
    }
}

// Whether moving an obstacle's base `distance` along its normal swallows one
// of the obstacle's faces whole (or, with `all`, every one): the face lies in
// the space the base sweeps, so a local operation would have to remove it.
bool swallows(const Obstacle& obstacle, double distance, double tolerance, bool all = false) {
    const double low = std::min(distance, 0.0) - tolerance;
    const double high = std::max(distance, 0.0) + tolerance;
    for (SmFace* face : obstacle.faces) {
        std::vector<SmPoint3d> points;
        face_points(face, points);
        SmPoint3d interior;
        if (face->CalculateAnInternalPoint(interior) == SM_SUCCESS) points.push_back(interior);
        bool inside = !points.empty();
        for (const SmPoint3d& p : points) {
            const double height = (p - obstacle.origin).Dot(obstacle.w);
            inside = inside && height >= low && height <= high;
        }
        if (inside != all) return inside;
    }
    return all;
}

// The features the faces' move runs into (see Obstacle), in the brep's face
// order. A planar face pushed in sweeps space, as do the blends moving with
// it, and a curved face moved either way; `motions` (plan_offset) says which
// faces those are, and the faces beside them are the ones whose boundaries
// move. On a picked plane itself, pushed in or pulled out, only a feature the
// move swallows a face of is in the way: a hole pushed in, a boss pulled out.
void find_obstacles(const SmBrep& brep, const SmTArray<SmFace*>& picked, double distance, const Motions& motions,
                    std::vector<Obstacle>& out) {
    std::vector<Push> pushes = pushes_of(picked, distance);
    add_blends(brep, picked, motions, pushes);
    const std::vector<Shell> shells = shells_of(picked, distance);
    const double tolerance = brep.GetTolerance();
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    auto moving = [&](const SmFace* face) {
        ULONG index = 0;
        return picked.FindElement(const_cast<SmFace*>(face), index) || motion_of(motions, face).kind != Motion::kKept;
    };
    std::unordered_set<const SmFace*> beside;
    for (ULONG i = 0; i < faces.GetSize() && !(pushes.empty() && shells.empty()); ++i) {
        if (!moving(faces[i])) continue;
        SmTArray<SmEdge*> edges;
        faces[i]->GetEdges(edges);
        for (ULONG e = 0; e < edges.GetSize(); ++e) {
            for (const SmFace* other : faces_of(edges[e]).faces) {
                if (other && !moving(other)) beside.insert(other);
            }
        }
    }
    std::unordered_set<const SmFace*> planes;
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        SmPoint3d on_plane;
        SmVector3d plane_normal;
        if (picked[i]->GetSurface() && plane_of(picked[i]->GetSurface(), on_plane, plane_normal)) {
            planes.insert(picked[i]);
        }
    }
    const int swallowed_side = distance < 0.0 ? -1 : 1;
    std::unordered_set<const SmFace*> taken;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        SmFace* base = faces[i];
        const bool on_picked = planes.count(base) > 0;
        if (!beside.count(base) && !on_picked) continue;
        SmTArray<SmLoop*> loops;
        base->GetLoops(loops);
        for (ULONG l = 0; l < loops.GetSize(); ++l) {
            if (loops[l]->IsOuterLoop()) continue;
            SmTArray<SmEdgeuse*> uses;
            loops[l]->GetEdgeuses(uses);
            // A loop along a moving face is the edge of what moves, such as
            // a boss's rim on the lid when its side moves.
            bool along_moving = false;
            for (ULONG e = 0; e < uses.GetSize() && !along_moving; ++e) {
                if (!uses[e]->GetEdge()) continue;
                for (const SmFace* other : faces_of(uses[e]->GetEdge()).faces) {
                    along_moving = along_moving || (other && other != base && moving(other));
                }
            }
            if (along_moving) continue;
            bool reached = on_picked;
            for (ULONG e = 0; e < uses.GetSize() && !reached; ++e) {
                std::vector<SmPoint3d> points;
                if (!uses[e]->GetEdge()) continue;
                edge_points(uses[e]->GetEdge(), 24, points);
                for (const SmPoint3d& p : points) {
                    for (const Push& push : pushes) {
                        reached = reached || push.sweeps(p, tolerance) || push.blend_sweeps(p, tolerance);
                    }
                    for (const Shell& shell : shells) reached = reached || shell.sweeps(p, tolerance);
                }
            }
            if (!reached) continue;
            // The faces standing on the loop: everything reachable from it
            // without crossing the base. Reaching a moving face means the
            // feature is part of what moves, which can't be set aside.
            Obstacle obstacle;
            obstacle.base = base;
            obstacle.on_picked = on_picked;
            std::vector<SmFace*> queue;
            std::unordered_set<const SmFace*> seen;
            bool escapes = false;
            for (ULONG e = 0; e < uses.GetSize(); ++e) {
                if (!uses[e]->GetEdge()) continue;
                for (SmFace* other : faces_of(uses[e]->GetEdge()).faces) {
                    if (other && other != base && seen.insert(other).second) queue.push_back(other);
                }
            }
            for (size_t q = 0; q < queue.size() && !escapes; ++q) {
                if (moving(queue[q]) || queue.size() > faces.GetSize() / 2 + 1) {
                    escapes = true;
                    break;
                }
                SmTArray<SmEdge*> edges;
                queue[q]->GetEdges(edges);
                for (ULONG e = 0; e < edges.GetSize(); ++e) {
                    for (SmFace* other : faces_of(edges[e]).faces) {
                        if (other && other != base && seen.insert(other).second) queue.push_back(other);
                    }
                }
            }
            // In the brep's order, so what follows is the same every run.
            for (ULONG f = 0; f < faces.GetSize(); ++f) {
                if (seen.count(faces[f])) obstacle.faces.push_back(faces[f]);
            }
            if (obstacle.faces.empty() || taken.count(obstacle.faces[0])) continue;
            if (!escapes && frame_of(base, obstacle.origin, obstacle.u, obstacle.v, obstacle.w)) {
                std::vector<SmPoint3d> points;
                for (SmFace* face : obstacle.faces) {
                    face_points(face, points);
                    box_corners(face, points);
                }
                double low[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
                double high[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
                bool out_side = true;
                bool in_side = true;
                std::vector<SmPoint3d> edge_samples;
                for (SmFace* face : obstacle.faces) face_points(face, edge_samples);
                for (const SmPoint3d& p : edge_samples) {
                    const double height = (p - obstacle.origin).Dot(obstacle.w);
                    out_side = out_side && height >= -10.0 * tolerance;
                    in_side = in_side && height <= 10.0 * tolerance;
                }
                for (const SmPoint3d& p : points) {
                    const SmVector3d d = p - obstacle.origin;
                    const double c[3] = {d.Dot(obstacle.u), d.Dot(obstacle.v), d.Dot(obstacle.w)};
                    for (int k = 0; k < 3; ++k) {
                        low[k] = std::min(low[k], c[k]);
                        high[k] = std::max(high[k], c[k]);
                    }
                }
                obstacle.side = out_side ? 1 : in_side ? -1 : 0;
                // A margin keeps the box's walls off the feature's own faces.
                const double margin =
                    0.5 + 0.02 * std::max({high[0] - low[0], high[1] - low[1], high[2] - low[2], 0.0});
                for (int k = 0; k < 2; ++k) {
                    obstacle.low[k] = low[k] - margin;
                    obstacle.high[k] = high[k] + margin;
                }
                obstacle.low[2] = obstacle.side > 0 ? 0.0 : low[2] - margin;
                obstacle.high[2] = obstacle.side > 0 ? high[2] + margin : 0.0;
            }
            // A feature on a picked face that the move leaves whole (a boss
            // on a face pushed in grows; a hole on one pulled out deepens)
            // isn't in the way, and neither is one that reaches round to
            // other faces, such as a hole right through.
            if (on_picked &&
                (escapes || obstacle.side != swallowed_side || !swallows(obstacle, distance, 10.0 * tolerance))) {
                continue;
            }
            for (SmFace* face : obstacle.faces) taken.insert(face);
            out.push_back(obstacle);
        }
    }
}

// A face beside a picked one that the picked face, pushed in, swallows whole:
// a chamfer along a lid's edge pushed past its width, a vent's chamfer when
// the vent widens past it, a post's when it thins. The picked face then meets
// the plane beyond it directly, as if the edge had never been chamfered,
// which is new topology. So the swallowed face is plugged, for the offset: a
// box in the plane's frame, from deeper than the swallowed face reaches to
// just outside the plane, covering it where it is and where the move takes
// it, and the picked face is extended through the box along its own surface,
// for which it must run along the plane's normal (a lid beside a chamfered
// wall, a vent's wall under its chamfered lid). A picked plane extends as
// the box, which ends on it, merges into it; a picked curved face as the
// box's inner face, where it bounds the space the face looks into, is
// pushed out past the plane. The picked face is then offset as usual, and
// the part cut by the result: what the offset takes from the plugged part.
struct Plug {
    // A frame on the plane beyond, w its outward normal, and the box in it.
    // For a picked plane, u is its outward normal and the box ends on it.
    SmPoint3d origin;
    SmVector3d u, v, w;
    double low[3] = {0.0, 0.0, 0.0};
    double high[3] = {0.0, 0.0, 0.0};
    // Whether its inner face is pushed out: for a picked curved face.
    bool opens = true;
};

// The faces across the face's edges, each once.
std::vector<SmFace*> neighbours_of(const SmFace* face) {
    std::vector<SmFace*> out;
    SmTArray<SmEdge*> edges;
    face->GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        for (SmFace* other : faces_of(edges[i]).faces) {
            if (other && other != face && std::find(out.begin(), out.end(), other) == out.end()) out.push_back(other);
        }
    }
    return out;
}

// Points that bound a face: along its edges, and one inside it.
void bounding_points(const SmFace* face, std::vector<SmPoint3d>& out) {
    face_points(face, out);
    SmPoint3d interior;
    if (face->CalculateAnInternalPoint(interior) == SM_SUCCESS) out.push_back(interior);
}

// How far `p` stands off the face's surface along its outward normal at the
// nearest point, and that normal.
bool height_off(const SmFace* face, const SmPoint3d& p, double& height, SmVector3d& normal) {
    SmPoint3d closest;
    SmTArray<SmPoint2d> uvs;
    double gap = 0.0;
    double n[3];
    if (SmApiSurfaceClosestPoint(face->GetSurface(), p, closest, uvs, gap) != SM_SUCCESS ||
        !outward_normal(face, closest, n)) {
        return false;
    }
    normal = SmVector3d(n[0], n[1], n[2]);
    height = (p - closest).Dot(normal);
    return true;
}

// The plugs a push needs (see Plug), in the brep's face order; none unless
// every swallowed face can be plugged, so that the offset then fails as it
// would have. Faces of `obstacles` are set aside instead.
void find_plugs(const SmBrep& brep, const SmTArray<SmFace*>& picked, double distance, const Motions& motions,
                const std::vector<Obstacle>& obstacles, std::vector<Plug>& out) {
    if (!(distance < 0.0)) return;
    const double tolerance = 10.0 * brep.GetTolerance();
    auto moving = [&](const SmFace* face) {
        ULONG index = 0;
        return picked.FindElement(const_cast<SmFace*>(face), index) || motion_of(motions, face).kind != Motion::kKept;
    };
    std::unordered_set<const SmFace*> set_aside;
    for (const Obstacle& obstacle : obstacles) set_aside.insert(obstacle.faces.begin(), obstacle.faces.end());
    // The kept faces beside a picked face that lie wholly between its
    // surface and its offset.
    std::unordered_map<const SmFace*, SmFace*> swallowed;  // and the picked face
    std::vector<SmFace*> order;
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        for (SmFace* face : neighbours_of(picked[i])) {
            if (moving(face) || swallowed.count(face) || set_aside.count(face)) continue;
            std::vector<SmPoint3d> points;
            bounding_points(face, points);
            bool inside = !points.empty();
            for (const SmPoint3d& p : points) {
                double height = 0.0;
                SmVector3d normal;
                inside = inside && height_off(picked[i], p, height, normal) && height >= distance - tolerance &&
                         height <= tolerance;
                if (!inside) break;
            }
            if (!inside) continue;
            swallowed[face] = picked[i];
            order.push_back(face);
        }
    }
    if (swallowed.empty()) return;
    std::unordered_set<const SmVertex*> on_moving;
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        if (!moving(faces[i])) continue;
        SmTArray<SmVertex*> vertices;
        faces[i]->GetVertices(vertices);
        for (ULONG k = 0; k < vertices.GetSize(); ++k) on_moving.insert(vertices[k]);
    }
    std::vector<Plug> plugs;
    for (SmFace* face : order) {
        SmFace* by = swallowed[face];
        // Beyond it: the one kept, unswallowed face it meets that touches no
        // moving face, not even at a corner, a plane the picked face runs
        // along. (The faces it meets at its ends, such as the next chamfer
        // along a lid's edges, touch the picked face at a corner.)
        SmFace* beyond = nullptr;
        int count = 0;
        bool alone = true;
        for (SmFace* other : neighbours_of(face)) {
            if (moving(other)) {
                alone = alone && other == by;
                continue;
            }
            if (swallowed.count(other)) continue;
            SmTArray<SmVertex*> vertices;
            other->GetVertices(vertices);
            bool touches = false;
            for (ULONG k = 0; k < vertices.GetSize(); ++k) touches = touches || on_moving.count(vertices[k]) > 0;
            if (touches) continue;
            beyond = other;
            ++count;
        }
        Plug plug;
        if (!alone || count != 1 || !frame_of(beyond, plug.origin, plug.u, plug.v, plug.w) ||
            !invariant_under(by->GetSurface(), plug.w)) {
            return;
        }
        SmPoint3d on_picked;
        SmVector3d picked_normal;
        if (plane_of(by->GetSurface(), on_picked, picked_normal)) {
            if (outward_normal_inside(by, on_picked, picked_normal) != SM_SUCCESS) return;
            plug.u = picked_normal;
            plug.v = plug.w * plug.u;  // cross product
            plug.opens = false;
        }
        // Where the face is and where the move takes it: each point carried
        // along the picked face's normal, as the move carries the face.
        std::vector<SmPoint3d> points;
        bounding_points(face, points);
        const size_t here = points.size();
        for (size_t k = 0; k < here; ++k) {
            double height = 0.0;
            SmVector3d normal;
            if (!height_off(by, points[k], height, normal)) return;
            points.push_back(points[k] + normal * distance);
        }
        double low[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
        double high[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
        for (const SmPoint3d& p : points) {
            const SmVector3d d = p - plug.origin;
            const double c[3] = {d.Dot(plug.u), d.Dot(plug.v), d.Dot(plug.w)};
            for (int k = 0; k < 3; ++k) {
                low[k] = std::min(low[k], c[k]);
                high[k] = std::max(high[k], c[k]);
            }
        }
        for (int k = 0; k < 3; ++k) {
            plug.low[k] = low[k];
            plug.high[k] = high[k];
        }
        if (!plug.opens) plug.high[0] = (on_picked - plug.origin).Dot(plug.u);
        plugs.push_back(plug);
    }
    // One margin for every plug, and every plug as deep as the deepest, so
    // that where plugs meet (along a chamfered lid's corners) they cover
    // each other's ends, and one push takes all their inner faces out.
    double margin = 0.0;
    double deepest = 0.0;
    for (const Plug& plug : plugs) {
        const double size = std::max({plug.high[0] - plug.low[0], plug.high[1] - plug.low[1], -plug.low[2]});
        margin = std::max(margin, 0.5 + 0.02 * size);
        deepest = std::min(deepest, plug.low[2]);
    }
    for (Plug& plug : plugs) {
        for (int k = 0; k < 2; ++k) {
            plug.low[k] -= margin;
            plug.high[k] += k == 0 && !plug.opens ? 0.0 : margin;
        }
        plug.low[2] = deepest - margin;
        plug.high[2] = margin;
    }
    // A plug may only cover the swallowed faces and what meets them: else it
    // fills more than the swallowed faces leave.
    std::unordered_set<const SmFace*> allowed;
    for (const auto& entry : swallowed) {
        allowed.insert(entry.first);
        SmTArray<SmVertex*> vertices;
        entry.first->GetVertices(vertices);
        for (ULONG i = 0; i < vertices.GetSize(); ++i) {
            SmTArray<SmFace*> faces;
            vertices[i]->GetFaces(faces);
            for (ULONG f = 0; f < faces.GetSize(); ++f) allowed.insert(faces[f]);
        }
    }
    for (const Plug& plug : plugs) {
        for (ULONG i = 0; i < faces.GetSize(); ++i) {
            if (allowed.count(faces[i])) continue;
            std::vector<SmPoint3d> points;
            bounding_points(faces[i], points);
            for (const SmPoint3d& p : points) {
                const SmVector3d d = p - plug.origin;
                const double c[3] = {d.Dot(plug.u), d.Dot(plug.v), d.Dot(plug.w)};
                bool inside = true;
                for (int k = 0; k < 3; ++k) {
                    inside = inside && c[k] > plug.low[k] + tolerance && c[k] < plug.high[k] - tolerance;
                }
                if (inside) return;
            }
        }
    }
    out = std::move(plugs);
}

// Whether any face but the obstacle's own and its base reaches into its box:
// then setting the box's contents aside would take more than the feature.
bool box_holds_only(const SmBrep& brep, const Obstacle& obstacle) {
    const double tolerance = 10.0 * brep.GetTolerance();
    auto inside = [&](const SmPoint3d& p) {
        const SmVector3d d = p - obstacle.origin;
        const double c[3] = {d.Dot(obstacle.u), d.Dot(obstacle.v), d.Dot(obstacle.w)};
        for (int k = 0; k < 3; ++k) {
            if (c[k] <= obstacle.low[k] + tolerance || c[k] >= obstacle.high[k] - tolerance) return false;
        }
        return true;
    };
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        if (faces[i] == obstacle.base ||
            std::find(obstacle.faces.begin(), obstacle.faces.end(), faces[i]) != obstacle.faces.end()) {
            continue;
        }
        std::vector<SmPoint3d> points;
        face_points(faces[i], points);
        SmPoint3d interior;
        if (faces[i]->CalculateAnInternalPoint(interior) == SM_SUCCESS) points.push_back(interior);
        for (const SmPoint3d& p : points) {
            if (inside(p)) return false;
        }
    }
    return true;
}

// A box from `low` to `high` in the frame (origin; u, v, and w = u x v).
sm_status frame_box(const SmPoint3d& origin, const SmVector3d& u, const SmVector3d& v, const double low[3],
                    const double high[3], std::unique_ptr<SmBrep>& out) {
    SmBrep* box = nullptr;
    SmApiStatus status = SmApiCreateBox(SmVector3d(low[0], low[1], low[2]), high[0] - low[0], high[1] - low[1],
                                        high[2] - low[2], box);
    if (status != SM_SUCCESS || !box) return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
    out.reset(box);
    return static_cast<sm_status>(SmApiTransform(box, SmAxis2Placement(origin, u, v)));
}

sm_status copy_of(const SmBrep& brep, std::unique_ptr<SmBrep>& out) {
    SmBrep* copy = nullptr;
    const SmApiStatus status = SmApiBrepCopy(&brep, copy);
    if (status != SM_SUCCESS || !copy) return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
    out.reset(copy);
    return SM_SUCCESS;
}

// `a` op `b` into `a`. The kernel modifies `a` in place and deletes `b` on
// success; on failure both stay allocated, and are freed here.
template <typename Op>
sm_status boolean_into(Op op, std::unique_ptr<SmBrep>& a, std::unique_ptr<SmBrep> b) {
    SmBrep* result = nullptr;
    const SmApiStatus status = op(a.get(), b.get(), result);
    if (status != SM_SUCCESS || !result) return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
    b.release();
    if (result != a.get()) a.reset(result);
    return SM_SUCCESS;
}

// Replaces the brep's faces with work's, in place, so the caller's handle
// stays good. MergeBrep rather than CopyFaces, which classifies regions again
// and can turn the solid inside out (a vent cut along a moved face).
sm_status take_faces(SmBrep& brep, SmBrep& work) {
    SmTemporaryChangeValue<SmBoolean> editing(brep.m_bEditingEnabled, TRUE);
    SmTArray<SmFace*> old_faces;
    brep.GetFaces(old_faces);
    const sm_status status = static_cast<sm_status>(brep.RemoveFaces(old_faces));
    if (status != SM_SUCCESS) return status;
    return static_cast<sm_status>(brep.MergeBrep(work, NULL, NULL, FALSE));
}

// Marks a picked face through the booleans that set obstacles aside.
const ULONG PICK_ATTRIBUTE_ID = 0x534d504b;  // "SMPK"

void mark_picked(SmFace* face, long index) {
    SmTArray<long> longs;
    longs.Add(index);
    face->AddAttribute(new (*face->GetContext())
                           SmGenericAttribute(PICK_ATTRIBUTE_ID, SM_AB_COPY, longs, SmTArray<double>(), SmTArray<char>()));
}

void unmark_all(SmBrep& brep) {
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        while (SmAttribute* mark = faces[i]->FindAttribute(PICK_ATTRIBUTE_ID)) faces[i]->RemoveAttribute(mark);
    }
}

// The picked faces again, by their marks; false if one is gone or split.
bool find_picked(const SmBrep& brep, ULONG count, SmTArray<SmFace*>& out) {
    std::vector<SmFace*> found(count, nullptr);
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const SmAttribute* mark = faces[i]->FindAttribute(PICK_ATTRIBUTE_ID);
        if (!mark || mark->GetNumLongElements() < 1) continue;
        const long index = mark->GetLongElementsAddress()[0];
        if (index < 0 || static_cast<ULONG>(index) >= count || found[index]) return false;
        found[index] = faces[i];
    }
    for (SmFace* face : found) {
        if (!face) return false;
        out.Add(face);
    }
    return true;
}

// A solid made by booleans held to what an offset is held to, but for its
// faces and edges: closed, every planar face's loops apart, no gap past the
// tolerance, and a volume its tessellation agrees with.
sm_status check_solid(SmBrep& brep) {
    if (!is_manifold_solid(brep)) return SM_SHIM_ERR_BROKE_SOLID;
    if (!loops_stay_apart(brep)) return SM_SHIM_ERR_BROKE_SOLID;
    double gap = 0.0;
    sm_status status = largest_gap(brep, gap);
    if (status != SM_SUCCESS) return status;
    if (gap > brep.GetTolerance()) return SM_SHIM_ERR_BROKE_SOLID;
    bool agrees = false;
    status = volume_agrees(brep, agrees);
    if (status != SM_SUCCESS) return status;
    return agrees ? SM_SUCCESS : SM_SHIM_ERR_BROKE_SOLID;
}

ULONG shell_count(const SmBrep& brep) {
    SmTArray<SmShell*> shells;
    brep.GetShells(shells);
    return shells.GetSize();
}

// The largest height of a brep's bounding box corners in an obstacle's
// frame: how far out of its base it reaches, at most.
double reach_out(const SmBrep& brep, const Obstacle& obstacle) {
    SmExtent3d bounds;
    if (brep.CalculateBoundingBox(bounds) != SM_SUCCESS) return HUGE_VAL;
    const SmPoint3d low = bounds.GetMin();
    const SmPoint3d high = bounds.GetMax();
    double furthest = -HUGE_VAL;
    for (int k = 0; k < 8; ++k) {
        const SmPoint3d corner(k & 1 ? high.x : low.x, k & 2 ? high.y : low.y, k & 4 ? high.z : low.z);
        furthest = std::max(furthest, (corner - obstacle.origin).Dot(obstacle.w));
    }
    return furthest;
}

// Whether two bounding boxes overlap, or come within `margin` of it.
bool boxes_meet(const SmExtent3d& a, const SmExtent3d& b, double margin) {
    return a.GetMin().x <= b.GetMax().x + margin && b.GetMin().x <= a.GetMax().x + margin &&
           a.GetMin().y <= b.GetMax().y + margin && b.GetMin().y <= a.GetMax().y + margin &&
           a.GetMin().z <= b.GetMax().z + margin && b.GetMin().z <= a.GetMax().z + margin;
}

// Cuts from a boss what stands over the region an offset swept, where that
// meets the boss's base: each face of `swept` on the base's plane, swept out
// of it past the boss, is taken from the boss. A boss beside a curved face
// pushed in is cut along the face's new surface so.
sm_status cut_over(const SmBrep& swept, const Obstacle& obstacle, std::unique_ptr<SmBrep>& boss) {
    const double tolerance = 10.0 * swept.GetTolerance();
    const double height = reach_out(*boss, obstacle) + 1.0;
    SmTArray<SmFace*> faces;
    swept.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        double normal[3];
        SmPoint3d inside;
        if (!outward_plane(faces[i], normal) || SmApiFaceInternalPoint(faces[i], inside) != SM_SUCCESS) continue;
        const SmVector3d n(normal[0], normal[1], normal[2]);
        if (n.Dot(obstacle.w) < 1.0 - 1e-9 || std::fabs((inside - obstacle.origin).Dot(obstacle.w)) > tolerance) {
            continue;
        }
        SmBrep* raw = nullptr;
        std::vector<sm_swept_face> roles;
        sm_status status = sweep_face(swept, faces[i], height, raw, roles);
        if (status != SM_SUCCESS) return status;
        std::unique_ptr<SmBrep> over(raw);
        SmExtent3d a;
        SmExtent3d b;
        if (over->CalculateBoundingBox(a) != SM_SUCCESS || boss->CalculateBoundingBox(b) != SM_SUCCESS) return SM_ERR;
        if (!boxes_meet(a, b, tolerance)) continue;
        if ((status = boolean_into(SmApiBooleanDifference, boss, std::move(over))) != SM_SUCCESS) return status;
        if (boss->GetNumFaces() == 0) return SM_SUCCESS;
    }
    return SM_SUCCESS;
}

// Where an offset of one face can sweep, as boxes in a frame (origin; u, v,
// and w = u x v): for a plane, in its frame, the slab between its old and
// new places across the whole part; for a curved face, in the world's axes,
// its box grown by the distance; and each blend moving with it, from where it
// is to where it goes. What the part has alike before and after away from
// the move, such as a vent beside a pushed side, is outside them.
struct Region {
    struct Box {
        double low[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
        double high[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    };
    SmPoint3d origin;
    SmVector3d axes[3];
    // The picked face's first: a plane's slab when `slab`, which spans the
    // part along u and v.
    std::vector<Box> boxes;
    bool slab = false;
    double distance = 0.0;

    void at(const SmPoint3d& p, double c[3]) const {
        for (int k = 0; k < 3; ++k) c[k] = (p - origin).Dot(axes[k]);
    }
    // Grows `box` to take in `p` and everything within `grow` of it, along
    // the first `count` axes.
    void add(Box& box, const SmPoint3d& p, double grow, int count = 3) const {
        double c[3];
        at(p, c);
        for (int k = 0; k < count; ++k) {
            box.low[k] = std::min(box.low[k], c[k] - grow);
            box.high[k] = std::max(box.high[k], c[k] + grow);
        }
    }
    // Whether `c` lies inside a box but the one at `except`, `by` clear of its
    // walls.
    bool holds(const double c[3], size_t except, double by) const {
        for (size_t b = 0; b < boxes.size(); ++b) {
            bool inside = b != except;
            for (int k = 0; k < 3; ++k) inside = inside && c[k] > boxes[b].low[k] + by && c[k] < boxes[b].high[k] - by;
            if (inside) return true;
        }
        return false;
    }
};

// The region an offset of `moved` by `distance` sweeps, short of a slab's
// span; false when that isn't one picked face.
bool region_of(const SmBrep& brep, const SmTArray<SmFace*>& moved, double distance, Region& out) {
    Motions motions;
    if (moved.GetSize() != 1 || plan_offset(moved, distance, brep.GetTolerance(), motions) != SM_SUCCESS) return false;
    out.slab = frame_of(moved[0], out.origin, out.axes[0], out.axes[1], out.axes[2]);
    if (!out.slab) {
        out.origin = SmPoint3d(0.0, 0.0, 0.0);
        out.axes[0] = SmVector3d(1.0, 0.0, 0.0);
        out.axes[1] = SmVector3d(0.0, 1.0, 0.0);
        out.axes[2] = SmVector3d(0.0, 0.0, 1.0);
    }
    out.distance = distance;
    out.boxes.emplace_back();
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const Motion motion = motion_of(motions, faces[i]);
        if (motion.kind == Motion::kKept) continue;
        if (faces[i] != moved[0]) out.boxes.emplace_back();
        Region::Box& box = faces[i] == moved[0] ? out.boxes.front() : out.boxes.back();
        std::vector<SmPoint3d> points;
        face_points(faces[i], points);
        box_corners(faces[i], points);
        for (const SmPoint3d& p : points) {
            if (motion.kind == Motion::kShift) {
                out.add(box, p, 0.0);
                out.add(box, p + motion.shift, 0.0);
            } else {
                out.add(box, p, std::fabs(distance));
            }
        }
    }
    return out.boxes.front().low[0] <= out.boxes.front().high[0];
}

// A face's extents in a region's frame; along which axes a wall across
// would cut it cleanly, the face running straight through it (a plane square
// to the wall, a cylinder along its normal); and whether it's a surface a
// boolean is slow to match where it coincides (canonical).
struct Extent {
    double low[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
    double high[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    bool clean[3] = {false, false, false};
    bool slow = false;
};

void add_extents(const Region& region, const SmBrep& brep, std::vector<Extent>& out) {
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const SmSurface* surface = faces[i]->GetSurface();
        Extent extent;
        std::vector<SmPoint3d> points;
        face_points(faces[i], points);
        box_corners(faces[i], points);
        for (const SmPoint3d& p : points) {
            double c[3];
            region.at(p, c);
            for (int k = 0; k < 3; ++k) {
                extent.low[k] = std::min(extent.low[k], c[k]);
                extent.high[k] = std::max(extent.high[k], c[k]);
            }
        }
        for (int k = 0; k < 3; ++k) extent.clean[k] = surface && invariant_under(surface, region.axes[k]);
        extent.slow = !surface || !canonical(surface);
        out.push_back(extent);
    }
}

// Whether a face with `extent` would meet the wall across `axis` at `wall`
// of the box at `index` only inside another box, where the wall isn't the
// region's.
bool met_inside(const Region& region, size_t index, int axis, double wall, const Extent& extent, double tolerance) {
    const Region::Box& box = region.boxes[index];
    for (size_t b = 0; b < region.boxes.size(); ++b) {
        const Region::Box& other = region.boxes[b];
        bool holds = b != index && wall > other.low[axis] + tolerance && wall < other.high[axis] - tolerance;
        for (int j = 0; j < 3; ++j) {
            holds = holds && (j == axis || (std::max(extent.low[j], box.low[j]) > other.low[j] + tolerance &&
                                            std::min(extent.high[j], box.high[j]) < other.high[j] - tolerance));
        }
        if (holds) return true;
    }
    return false;
}

// Sets the region's walls `margin` out from its boxes, then moves each wall
// off each face that it would cut but not cleanly (see Extent), or only
// touch: in between the face and what the box must hold, if there's room,
// or else out past the face. Cut any other way (a vent's wall from the side,
// its B-spline chamfer), what's left inside meets the whole face in the other
// operand, partly coincident, and the kernel can crash on that. False if the
// walls don't settle.
bool clear_walls(Region& region, const std::vector<Extent>& extents, double margin, double tolerance) {
    const std::vector<Region::Box> held = region.boxes;
    for (Region::Box& box : region.boxes) {
        for (int k = 0; k < 3; ++k) {
            box.low[k] -= margin;
            box.high[k] += margin;
        }
    }
    // The least room between a face and what a box holds to put a wall in.
    const double room = std::max(0.1, 100.0 * tolerance);
    for (int round = 0; round < 64; ++round) {
        bool moved = false;
        for (size_t b = 0; b < region.boxes.size(); ++b) {
            Region::Box& box = region.boxes[b];
            for (int k = 0; k < 3; ++k) {
                for (int side = 0; side < 2; ++side) {
                    double& wall = side ? box.high[k] : box.low[k];
                    for (const Extent& extent : extents) {
                        bool across = extent.low[k] <= wall + tolerance && extent.high[k] >= wall - tolerance;
                        for (int j = 0; j < 3; ++j) {
                            across = across && (j == k || (extent.low[j] <= box.high[j] + tolerance &&
                                                           extent.high[j] >= box.low[j] - tolerance));
                        }
                        const bool through = extent.low[k] < wall - tolerance && extent.high[k] > wall + tolerance;
                        if (!across || (through && extent.clean[k]) ||
                            met_inside(region, b, k, wall, extent, tolerance)) {
                            continue;
                        }
                        const double out = side ? 1.0 : -1.0;
                        const double hold = side ? held[b].high[k] : held[b].low[k];
                        const double near = side ? extent.low[k] : extent.high[k];
                        const double far = side ? extent.high[k] : extent.low[k];
                        wall = (near - hold) * out > room ? 0.5 * (hold + near) : far + out * margin;
                        moved = true;
                    }
                }
            }
        }
        if (!moved) return true;
    }
    return false;
}

// Whether the region leaves out a face that's slow to match (see Extent):
// otherwise taking the swept region inside it saves nothing.
bool leaves_out_slow(const Region& region, const std::vector<Extent>& extents) {
    for (const Extent& extent : extents) {
        bool apart = extent.slow;
        for (const Region::Box& box : region.boxes) {
            bool meets = true;
            for (int k = 0; k < 3; ++k) meets = meets && extent.low[k] <= box.high[k] && extent.high[k] >= box.low[k];
            apart = apart && !meets;
        }
        if (apart) return true;
    }
    return false;
}

// Whether a face of `brep` lies on the region's boundary: then the region
// cut through it.
bool reaches_walls(const SmBrep& brep, const Region& region) {
    const double tolerance = 10.0 * brep.GetTolerance();
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        double normal[3];
        if (!outward_plane(faces[i], normal)) continue;
        SmPoint3d inside;
        if (SmApiFaceInternalPoint(faces[i], inside) != SM_SUCCESS) return true;
        const SmVector3d n(normal[0], normal[1], normal[2]);
        double c[3];
        region.at(inside, c);
        for (size_t b = 0; b < region.boxes.size(); ++b) {
            const Region::Box& box = region.boxes[b];
            for (int k = 0; k < 3; ++k) {
                if (std::fabs(std::fabs(n.Dot(region.axes[k])) - 1.0) > 1e-9 ||
                    (std::fabs(c[k] - box.low[k]) > tolerance && std::fabs(c[k] - box.high[k]) > tolerance)) {
                    continue;
                }
                bool on_wall = true;
                for (int j = 0; j < 3; ++j) {
                    on_wall =
                        on_wall && (j == k || (c[j] >= box.low[j] - tolerance && c[j] <= box.high[j] + tolerance));
                }
                if (on_wall && !region.holds(c, b, tolerance)) return true;
            }
        }
    }
    return false;
}

// What an offset swept through, from the part `before` and `after` it: what
// it took, pushed in, or added, pulled out; into `before`. The two share
// most faces, and where they share one a boolean is slow to match (a vent's
// B-spline chamfer beside a pushed side), so when the region leaves such a
// face out, only what lies in the region is taken. That's the same unless
// the region cut through what was swept, and then the whole is taken after
// all.
sm_status swept_between(std::unique_ptr<SmBrep>& before, const SmBrep& after, bool pushed, const Region* given) {
    if (given) {
        Region region = *given;
        if (region.slab) {
            for (const SmBrep* brep : {static_cast<const SmBrep*>(before.get()), &after}) {
                SmExtent3d bounds;
                if (brep->CalculateBoundingBox(bounds) != SM_SUCCESS) return SM_ERR;
                const SmPoint3d low = bounds.GetMin();
                const SmPoint3d high = bounds.GetMax();
                for (int k = 0; k < 8; ++k) {
                    const SmPoint3d corner(k & 1 ? high.x : low.x, k & 2 ? high.y : low.y, k & 4 ? high.z : low.z);
                    region.add(region.boxes.front(), corner, 0.0, 2);
                }
            }
        }
        // A margin keeps the walls off what was swept.
        const double margin = 1.0 + 0.01 * std::fabs(region.distance);
        std::vector<Extent> extents;
        add_extents(region, *before, extents);
        add_extents(region, after, extents);
        const double tolerance = 10.0 * std::max(before->GetTolerance(), after.GetTolerance());
        bool good = clear_walls(region, extents, margin, tolerance) && leaves_out_slow(region, extents);
        std::unique_ptr<SmBrep> shape;
        for (size_t b = 0; b < region.boxes.size() && good; ++b) {
            std::unique_ptr<SmBrep> box;
            const Region::Box& at = region.boxes[b];
            good = frame_box(region.origin, region.axes[0], region.axes[1], at.low, at.high, box) == SM_SUCCESS;
            if (good && !shape) {
                shape = std::move(box);
            } else if (good) {
                good = boolean_into(SmApiBooleanUnion, shape, std::move(box)) == SM_SUCCESS;
            }
        }
        std::unique_ptr<SmBrep> inside;
        std::unique_ptr<SmBrep> other;
        if (good && copy_of(pushed ? *before : after, inside) == SM_SUCCESS &&
            boolean_into(SmApiBooleanIntersection, inside, std::move(shape)) == SM_SUCCESS &&
            copy_of(pushed ? after : *before, other) == SM_SUCCESS &&
            boolean_into(SmApiBooleanDifference, inside, std::move(other)) == SM_SUCCESS &&
            !reaches_walls(*inside, region)) {
            before = std::move(inside);
            return SM_SUCCESS;
        }
    }
    std::unique_ptr<SmBrep> other;
    sm_status status = copy_of(after, other);
    if (status != SM_SUCCESS) return status;
    if (pushed) return boolean_into(SmApiBooleanDifference, before, std::move(other));
    status = boolean_into(SmApiBooleanDifference, other, std::move(before));
    if (status != SM_SUCCESS) return status;
    before = std::move(other);
    return SM_SUCCESS;
}

// The plugs' inner faces in `work`: on a plug's inner plane, meeting a
// picked face.
sm_status inner_faces(const SmBrep& work, ULONG picked_count, const std::vector<Plug>& plugs,
                      SmTArray<SmFace*>& inner) {
    SmTArray<SmFace*> picked;
    if (!find_picked(work, picked_count, picked)) return SM_SHIM_ERR_OBSTRUCTED;
    const double tolerance = 10.0 * work.GetTolerance();
    SmTArray<SmFace*> faces;
    work.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        double normal[3];
        SmPoint3d inside;
        if (!outward_plane(faces[i], normal) || SmApiFaceInternalPoint(faces[i], inside) != SM_SUCCESS) continue;
        const SmVector3d n(normal[0], normal[1], normal[2]);
        for (const Plug& plug : plugs) {
            if (!plug.opens) continue;
            const SmVector3d d = inside - plug.origin;
            const double c[3] = {d.Dot(plug.u), d.Dot(plug.v), d.Dot(plug.w)};
            if (std::fabs(n.Dot(plug.w)) < 1.0 - 1e-9 || std::fabs(c[2] - plug.low[2]) > tolerance ||
                c[0] < plug.low[0] || c[0] > plug.high[0] || c[1] < plug.low[1] || c[1] > plug.high[1]) {
                continue;
            }
            bool meets = false;
            for (SmFace* other : neighbours_of(faces[i])) {
                ULONG index = 0;
                meets = meets || picked.FindElement(other, index);
            }
            if (meets) inner.Add(faces[i]);
            break;
        }
    }
    return inner.GetSize() > 0 ? SM_SUCCESS : SM_SHIM_ERR_OBSTRUCTED;
}

// The plug an inner face lies on.
const Plug* plug_under(const SmFace* face, const std::vector<Plug>& plugs, double tolerance) {
    SmPoint3d inside;
    if (SmApiFaceInternalPoint(face, inside) != SM_SUCCESS) return nullptr;
    for (const Plug& plug : plugs) {
        if (plug.opens && std::fabs((inside - plug.origin).Dot(plug.w) - plug.low[2]) <= tolerance) return &plug;
    }
    return nullptr;
}

bool same_surface(const SmFace* a, const SmFace* b, double tolerance);

// Edges that already join two faces on one surface when an offset starts,
// as an imprint leaves them: the push-through's merges keep them, merging
// only the faces its own boxes and cuts split. Split edges pass the mark to
// both pieces, and copies keep it.
const ULONG KEPT_EDGE_ATTRIBUTE_ID = 0x534d4b45;  // "SMKE"

void mark_kept_edges(SmBrep& brep) {
    const double tolerance = brep.GetTolerance();
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        const EdgeFaces sides = faces_of(edges[i]);
        if (sides.use_count != 2 || !sides.faces[0] || !sides.faces[1] || sides.faces[0] == sides.faces[1] ||
            !same_surface(sides.faces[0], sides.faces[1], tolerance)) {
            continue;
        }
        edges[i]->AddAttribute(new (*edges[i]->GetContext()) SmGenericAttribute(
            KEPT_EDGE_ATTRIBUTE_ID, SM_AB_COPY, SmTArray<long>(), SmTArray<double>(), SmTArray<char>()));
    }
}

void unmark_kept_edges(SmBrep& brep) {
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        while (SmAttribute* mark = edges[i]->FindAttribute(KEPT_EDGE_ATTRIBUTE_ID)) edges[i]->RemoveAttribute(mark);
    }
}

// RemoveTopologicalEdgesAndVertices, keeping the marked edges (see
// mark_kept_edges). With none marked it's the plain merge, as before.
sm_status merge_unmarked(SmBrep& brep) {
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    SmTArray<SmEdge*> kept;
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        if (edges[i]->FindAttribute(KEPT_EDGE_ATTRIBUTE_ID)) kept.Add(edges[i]);
    }
    if (kept.GetSize() == 0) return static_cast<sm_status>(brep.RemoveTopologicalEdgesAndVertices());
    return static_cast<sm_status>(brep.RemoveTopologicalEdgesAndVertices(&kept));
}

// Unmarks a brep's kept edges when it goes out of scope, whatever the
// push-through returns.
struct KeptEdges {
    SmBrep& brep;
    explicit KeptEdges(SmBrep& marked) : brep(marked) { mark_kept_edges(brep); }
    ~KeptEdges() { unmark_kept_edges(brep); }
    KeptEdges(const KeptEdges&) = delete;
    KeptEdges& operator=(const KeptEdges&) = delete;
};

// Fills `work` with the plugs (see Plug). A plug over a blind hole's mouth
// seals the hole in, a void, whose faces nothing here orients right; so each
// inner face over a void is then bored through, along the plug's normal,
// clear of its edges, which opens the void again and leaves the inner face a
// ring whose bore stretches as it is pushed out.
sm_status fill_plugs(std::unique_ptr<SmBrep>& work, ULONG picked_count, const std::vector<Plug>& plugs) {
    const ULONG shells = shell_count(*work);
    sm_status status = SM_SUCCESS;
    for (const Plug& plug : plugs) {
        std::unique_ptr<SmBrep> box;
        if ((status = frame_box(plug.origin, plug.u, plug.v, plug.low, plug.high, box)) != SM_SUCCESS ||
            (status = boolean_into(SmApiBooleanUnion, work, std::move(box))) != SM_SUCCESS) {
            return status;
        }
    }
    if ((status = merge_unmarked(*work)) != SM_SUCCESS) return status;
    if (shell_count(*work) == shells) return SM_SUCCESS;
    SmTArray<SmFace*> inner;
    if ((status = inner_faces(*work, picked_count, plugs, inner)) != SM_SUCCESS) return status;
    const double tolerance = 10.0 * work->GetTolerance();
    std::vector<std::unique_ptr<SmBrep>> bores;
    for (ULONG i = 0; i < inner.GetSize(); ++i) {
        const Plug* plug = plug_under(inner[i], plugs, tolerance);
        SmPoint3d centre;
        if (!plug || SmApiFaceInternalPoint(inner[i], centre) != SM_SUCCESS) return SM_SHIM_ERR_OBSTRUCTED;
        std::vector<SmPoint3d> points;
        face_points(inner[i], points);
        double clear = HUGE_VAL;
        for (const SmPoint3d& p : points) {
            const SmVector3d d = p - centre;
            clear = std::min(clear, (d - plug->w * d.Dot(plug->w)).Length());
        }
        const double radius = 0.5 * clear;
        // Into the void a little, and out of the plug's outer face.
        const double past = 0.05;
        if (!(radius > tolerance) || points.empty()) return SM_SHIM_ERR_OBSTRUCTED;
        const SmVector3d d = centre - plug->origin;
        SmBrep* raw = nullptr;
        SmApiStatus made = SmApiCreateCylinder(SmVector3d(d.Dot(plug->u), d.Dot(plug->v), plug->low[2] - past),
                                               radius, plug->high[2] - plug->low[2] + 2.0 * past, raw);
        if (made != SM_SUCCESS || !raw) return made != SM_SUCCESS ? static_cast<sm_status>(made) : SM_ERR;
        std::unique_ptr<SmBrep> bore(raw);
        made = SmApiTransform(raw, SmAxis2Placement(plug->origin, plug->u, plug->v));
        if (made != SM_SUCCESS) return static_cast<sm_status>(made);
        bores.push_back(std::move(bore));
    }
    for (auto& bore : bores) {
        if ((status = boolean_into(SmApiBooleanDifference, work, std::move(bore))) != SM_SUCCESS) return status;
    }
    if ((status = merge_unmarked(*work)) != SM_SUCCESS) return status;
    return shell_count(*work) == shells ? SM_SUCCESS : SM_SHIM_ERR_OBSTRUCTED;
}

// Pushes the plugs' inner faces in `work` out past the planes beyond, so
// the picked faces reach through them.
sm_status open_plugs(SmBrep& work, ULONG picked_count, const std::vector<Plug>& plugs) {
    if (std::none_of(plugs.begin(), plugs.end(), [](const Plug& plug) { return plug.opens; })) return SM_SUCCESS;
    SmTArray<SmFace*> inner;
    sm_status status = inner_faces(work, picked_count, plugs, inner);
    if (status != SM_SUCCESS) return status;
    const double push = plugs[0].low[2] - plugs[0].high[2] / 2.0;
    return keeping_topology(work, [&] { return offset_faces(work, inner, push); });
}

// Offsets the picked faces through `obstacles` (see find_obstacles) and the
// faces they swallow (see Plug), in place.
sm_status push_through(SmBrep& brep, const SmTArray<SmFace*>& picked, double distance,
                       const std::vector<Obstacle>& obstacles, const std::vector<Plug>& plugs) {
    for (const Obstacle& obstacle : obstacles) {
        if (obstacle.side == 0 || !box_holds_only(brep, obstacle)) return SM_SHIM_ERR_OBSTRUCTED;
    }
    const std::vector<Push> pushes = pushes_of(picked, distance);
    const bool pushed = distance < 0.0;
    const bool curved = !shells_of(picked, distance).empty();
    const KeptEdges kept(brep);
    for (ULONG i = 0; i < picked.GetSize(); ++i) mark_picked(picked[i], static_cast<long>(i));
    std::unique_ptr<SmBrep> work;
    sm_status status = copy_of(brep, work);
    unmark_all(brep);
    if (status != SM_SUCCESS) return status;

    // Set each obstacle aside: a boss's material, taken out; then a hole's
    // space, filled in, keeping the part as it was with its holes.
    std::vector<std::pair<const Obstacle*, std::unique_ptr<SmBrep>>> bosses;
    for (const Obstacle& obstacle : obstacles) {
        if (obstacle.side < 0) continue;
        std::unique_ptr<SmBrep> box;
        std::unique_ptr<SmBrep> second;
        std::unique_ptr<SmBrep> part;
        if ((status = frame_box(obstacle.origin, obstacle.u, obstacle.v, obstacle.low, obstacle.high, box)) !=
                SM_SUCCESS ||
            (status = copy_of(*box, second)) != SM_SUCCESS || (status = copy_of(brep, part)) != SM_SUCCESS ||
            (status = boolean_into(SmApiBooleanIntersection, part, std::move(box))) != SM_SUCCESS ||
            (status = boolean_into(SmApiBooleanDifference, work, std::move(second))) != SM_SUCCESS) {
            return status;
        }
        bosses.emplace_back(&obstacle, std::move(part));
    }
    // A hole on a picked face that the push swallows whole has nothing to
    // come back as; any other, and the part around a plug, needs the part as
    // it is to come back from.
    std::unique_ptr<SmBrep> holed;
    const double tolerance = 10.0 * brep.GetTolerance();
    if (!plugs.empty() && (status = copy_of(*work, holed)) != SM_SUCCESS) return status;
    for (const Obstacle& obstacle : obstacles) {
        if (obstacle.side > 0) continue;
        const bool gone = obstacle.on_picked && swallows(obstacle, distance, tolerance, true);
        if (!gone && !holed && (status = copy_of(*work, holed)) != SM_SUCCESS) return status;
    }
    for (const Obstacle& obstacle : obstacles) {
        if (obstacle.side > 0) continue;
        std::unique_ptr<SmBrep> box;
        if ((status = frame_box(obstacle.origin, obstacle.u, obstacle.v, obstacle.low, obstacle.high, box)) !=
                SM_SUCCESS ||
            (status = boolean_into(SmApiBooleanUnion, work, std::move(box))) != SM_SUCCESS) {
            return status;
        }
    }
    if ((status = merge_unmarked(*work)) != SM_SUCCESS) return status;
    if (!plugs.empty() && (status = fill_plugs(work, picked.GetSize(), plugs)) != SM_SUCCESS) return status;
    // What the offset sweeps through, from the filled part before and after:
    // what it takes, pushed in, or adds, pulled out. Holes come back by it,
    // and bosses beside a curved face pushed in are cut by it; around a plug,
    // it takes the swallowed faces, and what the plug opened only takes
    // what's empty in the part.
    const bool cut_bosses = curved && pushed && !bosses.empty();
    std::unique_ptr<SmBrep> swept;
    if ((holed || cut_bosses) && (status = copy_of(*work, swept)) != SM_SUCCESS) return status;
    if (!plugs.empty() && (status = open_plugs(*work, picked.GetSize(), plugs)) != SM_SUCCESS) return status;
    SmTArray<SmFace*> moved;
    if (!find_picked(*work, picked.GetSize(), moved)) return SM_SHIM_ERR_OBSTRUCTED;
    // The swept region is taken only around the move (see swept_between),
    // but not around a plug: what opening it took is no move region_of
    // knows of, and with the swallowed face plugged the whole is quick.
    Region region;
    const bool bounded = swept && plugs.empty() && region_of(*work, moved, distance, region);
    status = keeping_topology(*work, [&] { return offset_faces(*work, moved, distance); });
    if (status != SM_SUCCESS) return status;
    if (swept) {
        if ((status = swept_between(swept, *work, pushed, bounded ? &region : nullptr)) != SM_SUCCESS) return status;
    }

    // The holes back: the part with its holes, less what the offset took (or
    // with what it added). That cuts what the moved face runs through along
    // its new surface, as a cut of the part with a slab would; cutting each
    // hole's space back out of the offset part instead can leave a broken
    // solid, as where a vent's chamfer meets the lid.
    if (holed) {
        std::unique_ptr<SmBrep> sweep;
        if (cut_bosses) {
            if ((status = copy_of(*swept, sweep)) != SM_SUCCESS) return status;
        } else {
            sweep = std::move(swept);
        }
        status = boolean_into(pushed ? SmApiBooleanDifference : SmApiBooleanUnion, holed, std::move(sweep));
        if (status != SM_SUCCESS) return status;
        work = std::move(holed);
    }

    // Each boss back as far as the moved faces leave it: behind a plane
    // pushed in, off what a curved face pushed in swept, and not at all under
    // a plane pulled out past it.
    for (auto& [obstacle, boss] : bosses) {
        bool gone = obstacle->on_picked && distance > 0.0 && reach_out(*boss, *obstacle) <= distance + tolerance;
        if (!gone && cut_bosses) {
            if ((status = cut_over(*swept, *obstacle, boss)) != SM_SUCCESS) return status;
            gone = boss->GetNumFaces() == 0;
        }
        for (const Push& push : pushes) {
            if (gone || curved) break;
            SmExtent3d bounds;
            if (boss->CalculateBoundingBox(bounds) != SM_SUCCESS) return SM_ERR;
            const SmPoint3d plane = push.point - push.normal * push.distance;
            double nearest = HUGE_VAL;
            double furthest = -HUGE_VAL;
            const SmPoint3d low = bounds.GetMin();
            const SmPoint3d high = bounds.GetMax();
            for (int k = 0; k < 8; ++k) {
                const SmPoint3d corner(k & 1 ? high.x : low.x, k & 2 ? high.y : low.y, k & 4 ? high.z : low.z);
                const double beyond = (corner - plane).Dot(push.normal);
                nearest = std::min(nearest, beyond);
                furthest = std::max(furthest, beyond);
            }
            if (nearest >= -brep.GetTolerance()) {
                gone = true;
                break;
            }
            if (furthest <= brep.GetTolerance()) continue;
            // Keep what's behind the plane: a box on it reaching back past
            // the boss, wide enough to cover it.
            const double reach = (high - low).Length() + 1.0;
            const SmVector3d centre = (SmVector3d(low.x, low.y, low.z) + SmVector3d(high.x, high.y, high.z)) * 0.5;
            const SmVector3d offset = SmPoint3d(centre.x, centre.y, centre.z) - plane;
            const double box_low[3] = {offset.Dot(push.u) - reach, offset.Dot(push.v) - reach, -2.0 * reach};
            const double box_high[3] = {offset.Dot(push.u) + reach, offset.Dot(push.v) + reach, 0.0};
            std::unique_ptr<SmBrep> behind;
            if ((status = frame_box(plane, push.u, push.v, box_low, box_high, behind)) != SM_SUCCESS ||
                (status = boolean_into(SmApiBooleanIntersection, boss, std::move(behind))) != SM_SUCCESS) {
                return status;
            }
        }
        if (gone) continue;
        if ((status = boolean_into(SmApiBooleanUnion, work, std::move(boss))) != SM_SUCCESS) return status;
    }
    if ((status = merge_unmarked(*work)) != SM_SUCCESS) return status;
    unmark_all(*work);
    if ((status = check_solid(*work)) != SM_SUCCESS) return status;
    // A boss whose face pulled back from under it (behind a round that moved
    // with the pushed face) would come back on its own.
    if (shell_count(*work) != shell_count(brep)) return SM_SHIM_ERR_OBSTRUCTED;

    return take_faces(brep, *work);
}

// Steps. A picked face beside an unpicked face on its own surface, as the two
// pieces of a face an imprint splits are, can't stretch that neighbour to
// meet it: the neighbour's surface runs alongside the face's new place and
// never reaches it. The move needs a new face there, a step wall swept along
// the edge between them, and the offset executive never makes faces.
//
// For planar faces all moving one way, the step is what sweeping the faces
// along the move into prisms and joining them (pulled out) or cutting them
// (pushed in) makes, as long as every other neighbour stands along the move:
// then each prism side beside one lies on its surface and merges into it,
// stretched or trimmed as the offset would. The prisms must add or take
// exactly their own volume; anything else means the move runs into
// something, and it's refused rather than cut through.

// How much of the prisms' own volume a step may add or take beyond them:
// well past the 1e-8 volumes are found to, and far short of any
// obstruction worth the name.
const double STEP_VOLUME_TOLERANCE = 1e-6;

// Whether two faces lie on one surface: the same surface, the same plane, or
// the same cylinder, sphere or torus.
bool same_surface(const SmFace* a, const SmFace* b, double tolerance) {
    const SmSurface* first = a->GetSurface();
    const SmSurface* second = b->GetSurface();
    if (!first || !second) return false;
    if (first == second) return true;
    SmPoint3d point_a, point_b;
    SmVector3d normal_a, normal_b;
    const bool planar_a = plane_of(first, point_a, normal_a);
    const bool planar_b = plane_of(second, point_b, normal_b);
    if (planar_a || planar_b) {
        return planar_a && planar_b && 1.0 - std::fabs(normal_a.Dot(normal_b)) < 1e-12 &&
               std::fabs((point_b - point_a).Dot(normal_a)) <= tolerance;
    }
    const sm_face_info info_a = describe_face(a);
    const sm_face_info info_b = describe_face(b);
    if (info_a.kind != info_b.kind || std::fabs(info_a.radius - info_b.radius) > tolerance ||
        std::fabs(info_a.minor_radius - info_b.minor_radius) > tolerance) {
        return false;
    }
    const SmVector3d axis(info_a.axis[0], info_a.axis[1], info_a.axis[2]);
    const SmVector3d other_axis(info_b.axis[0], info_b.axis[1], info_b.axis[2]);
    const SmVector3d apart(info_b.origin[0] - info_a.origin[0], info_b.origin[1] - info_a.origin[1],
                           info_b.origin[2] - info_a.origin[2]);
    switch (info_a.kind) {
        case SM_SURFACE_CYLINDER:
            return 1.0 - std::fabs(axis.Dot(other_axis)) < 1e-12 &&
                   (apart - axis * apart.Dot(axis)).Length() <= tolerance;
        case SM_SURFACE_SPHERE:
            return apart.Length() <= tolerance;
        case SM_SURFACE_TORUS:
            return 1.0 - std::fabs(axis.Dot(other_axis)) < 1e-12 && apart.Length() <= tolerance;
        default:
            return false;
    }
}

// Whether a picked face meets an unpicked face on its own surface.
bool needs_steps(const SmTArray<SmFace*>& picked, double tolerance) {
    std::unordered_set<const SmFace*> moving;
    for (ULONG i = 0; i < picked.GetSize(); ++i) moving.insert(picked[i]);
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        for (const SmFace* other : neighbours_of(picked[i])) {
            if (!moving.count(other) && same_surface(picked[i], other, tolerance)) return true;
        }
    }
    return false;
}

sm_status brep_volume(SmBrep& brep, double& volume) {
    return static_cast<sm_status>(SmApiBrepComputeVolume(&brep, 1e-8, volume));
}

// How much of the prisms a step would add or take that isn't theirs: what of
// them lies inside the part when pulled out, or outside it when pushed in.
// Booleans with the small prisms first, so only what's left, usually nothing,
// is measured, not the whole part.
sm_status step_overlap(const SmBrep& brep, const SmBrep& tool, bool pulled, double& overlap) {
    overlap = 0.0;
    std::unique_ptr<SmBrep> part;
    std::unique_ptr<SmBrep> prisms;
    sm_status status = SM_SUCCESS;
    if ((status = copy_of(brep, part)) != SM_SUCCESS || (status = copy_of(tool, prisms)) != SM_SUCCESS) return status;
    status = pulled ? boolean_into(SmApiBooleanIntersection, prisms, std::move(part))
                    : boolean_into(SmApiBooleanDifference, prisms, std::move(part));
    if (status != SM_SUCCESS || prisms->GetNumFaces() == 0) return status;
    return brep_volume(*prisms, overlap);
}

// The move with step walls (see above). SM_SHIM_ERR_UNSUPPORTED_STEP when a
// picked face isn't planar, they don't all face one way, or a neighbour that
// isn't on a picked face's own surface doesn't stand along the move;
// SM_SHIM_ERR_OBSTRUCTED when the prisms reach into the part pulled out, or
// out of it pushed in, so wouldn't add or take exactly their own volume.
// The moved faces keep their labels, and each wall takes the label of the
// edge it stands on, if that has one, so callers can tell which it is.
sm_status offset_with_steps(SmBrep& brep, const SmTArray<SmFace*>& picked, double distance) {
    const double tolerance = brep.GetTolerance();
    std::unordered_set<const SmFace*> moving;
    for (ULONG i = 0; i < picked.GetSize(); ++i) moving.insert(picked[i]);
    SmVector3d move;
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        double outward[3];
        if (!outward_plane(picked[i], outward)) return SM_SHIM_ERR_UNSUPPORTED_STEP;
        const SmVector3d normal(outward[0], outward[1], outward[2]);
        if (i == 0) {
            move = normal;
        } else if (normal.Dot(move) < 1.0 - 1e-12) {
            return SM_SHIM_ERR_UNSUPPORTED_STEP;
        }
        for (const SmFace* other : neighbours_of(picked[i])) {
            if (moving.count(other) || same_surface(picked[i], other, tolerance)) continue;
            if (!other->GetSurface() || !invariant_under(other->GetSurface(), move)) {
                return SM_SHIM_ERR_UNSUPPORTED_STEP;
            }
        }
    }

    std::unique_ptr<SmBrep> tool;
    sm_status status = SM_SUCCESS;
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    for (ULONG i = 0; i < picked.GetSize(); ++i) {
        SmBrep* swept = nullptr;
        std::vector<sm_swept_face> roles;
        if ((status = sweep_face(brep, picked[i], distance, swept, roles)) != SM_SUCCESS) return status;
        std::unique_ptr<SmBrep> prism(swept);
        // The prism's far end is the face where it moved to, and a side
        // swept from an edge the face shares with another piece of its own
        // surface is a wall.
        SmTArray<SmFace*> faces;
        prism->GetFaces(faces);
        for (ULONG f = 0; f < faces.GetSize(); ++f) {
            if (roles[f].role == SM_SWEPT_END) write_label(*faces[f], read_label(picked[i]));
            const int32_t source = roles[f].source_edge;
            if (roles[f].role != SM_SWEPT_SIDE || source < 0 || static_cast<ULONG>(source) >= edges.GetSize()) continue;
            const EdgeFaces sides = faces_of(edges[static_cast<ULONG>(source)]);
            const SmFace* other = sides.faces[0] == picked[i] ? sides.faces[1] : sides.faces[0];
            if (other && !moving.count(other) && same_surface(picked[i], other, tolerance)) {
                const sm_label label = read_label(edges[static_cast<ULONG>(source)]);
                if (label.present) write_label(*faces[f], label);
            }
        }
        if (!tool) {
            tool = std::move(prism);
        } else if ((status = boolean_into(SmApiBooleanUnion, tool, std::move(prism))) != SM_SUCCESS) {
            return status;
        }
    }
    double swept = 0.0;
    double overlap = 0.0;
    if ((status = brep_volume(*tool, swept)) != SM_SUCCESS ||
        (status = step_overlap(brep, *tool, distance > 0.0, overlap)) != SM_SUCCESS) {
        return status;
    }
    if (std::fabs(overlap) > STEP_VOLUME_TOLERANCE * swept) return SM_SHIM_ERR_OBSTRUCTED;
    std::unique_ptr<SmBrep> work;
    if ((status = copy_of(brep, work)) != SM_SUCCESS) return status;
    status = distance > 0.0 ? boolean_into(SmApiBooleanUnion, work, std::move(tool))
                            : boolean_into(SmApiBooleanDifference, work, std::move(tool));
    if (status != SM_SUCCESS) return status;
    if ((status = check_solid(*work)) != SM_SUCCESS) return status;
    if (shell_count(*work) != shell_count(brep)) return SM_SHIM_ERR_OBSTRUCTED;
    return take_faces(brep, *work);
}

// The offset with whatever it runs into cut through, or as before when
// nothing is in its way; with step walls where a face moves away from
// another piece of its own surface.
sm_status offset_through(SmBrep& brep, const SmTArray<SmFace*>& picked, double distance) {
    if (needs_steps(picked, brep.GetTolerance())) return offset_with_steps(brep, picked, distance);
    Motions motions;
    std::vector<Obstacle> obstacles;
    std::vector<Plug> plugs;
    if (plan_offset(picked, distance, brep.GetTolerance(), motions) == SM_SUCCESS) {
        find_obstacles(brep, picked, distance, motions, obstacles);
        find_plugs(brep, picked, distance, motions, obstacles, plugs);
    }
    if (obstacles.empty() && plugs.empty()) {
        return keeping_topology(brep, [&] { return offset_faces(brep, picked, distance); });
    }
    return push_through(brep, picked, distance, obstacles, plugs);
}

}  // namespace

extern "C" {

sm_status sm_context_create(void) {
    return guarded([] {
        SmApiCreateContext();
        return SmApiGetOrCreateContext() ? SM_SUCCESS : SM_ERR;
    });
}

sm_status sm_brep_create_box(sm_vec3 origin, double length, double width, double height, sm_brep** out) {
    if (!out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* brep = nullptr;
        const SmApiStatus status = SmApiCreateBox(SmVector3d(origin.x, origin.y, origin.z), length, width, height, brep);
        if (status != SM_SUCCESS) return status;
        *out = handle_of(brep);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_create_cylinder(sm_vec3 base_center, double radius, double height, sm_brep** out) {
    if (!out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* brep = nullptr;
        const SmApiStatus status =
            SmApiCreateCylinder(SmVector3d(base_center.x, base_center.y, base_center.z), radius, height, brep);
        if (status != SM_SUCCESS) return status;
        *out = handle_of(brep);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_create_sphere(sm_vec3 center, double radius, sm_brep** out) {
    if (!out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* brep = nullptr;
        const SmApiStatus status = SmApiCreateSphere(SmVector3d(center.x, center.y, center.z), radius, brep);
        if (status != SM_SUCCESS) return status;
        *out = handle_of(brep);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_copy(const sm_brep* brep, sm_brep** out) {
    if (!brep || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* copy = nullptr;
        const SmApiStatus status = SmApiBrepCopy(brep_of(brep), copy);
        if (status != SM_SUCCESS) return status;
        *out = handle_of(copy);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_union(sm_brep* a, sm_brep* b, sm_brep** out) { return run_boolean(SmApiBooleanUnion, a, b, out); }

sm_status sm_brep_difference(sm_brep* a, sm_brep* b, sm_brep** out) {
    return run_boolean(SmApiBooleanDifference, a, b, out);
}

sm_status sm_brep_intersection(sm_brep* a, sm_brep* b, sm_brep** out) {
    return run_boolean(SmApiBooleanIntersection, a, b, out);
}

sm_status sm_brep_is_manifold_solid(const sm_brep* brep, int32_t* out) {
    if (!brep || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBoolean manifold = FALSE;
        const SmApiStatus status = SmApiBrepIsManifoldSolid(brep_of(brep), manifold);
        if (status != SM_SUCCESS) return status;
        *out = manifold ? 1 : 0;
        return SM_SUCCESS;
    });
}

sm_status sm_brep_volume(const sm_brep* brep, double relative_accuracy, double* out) {
    if (!brep || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        double volume = 0.0;
        const SmApiStatus status = SmApiBrepComputeVolume(brep_of(brep), relative_accuracy, volume);
        if (status != SM_SUCCESS) return status;
        *out = volume;
        return SM_SUCCESS;
    });
}

sm_status sm_brep_counts(const sm_brep* brep, size_t* faces, size_t* edges) {
    if (!brep || !faces || !edges) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        *faces = brep_of(brep)->GetNumFaces();
        *edges = brep_of(brep)->GetNumEdges();
        return SM_SUCCESS;
    });
}

void sm_brep_free(sm_brep* brep) {
    try {
        delete brep_of(brep);
    } catch (...) {
    }
}

sm_status sm_brep_edge_geometry(const sm_brep* brep, sm_edge_geometry* out, size_t capacity, size_t* count) {
    if (!brep || !count || (capacity > 0 && !out)) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        const SmBrep* source = brep_of(brep);
        SmTArray<SmEdge*> edges;
        source->GetEdges(edges);
        const FaceIndex face_index = index_faces(*source);
        const size_t written = std::min(capacity, static_cast<size_t>(edges.GetSize()));
        for (size_t i = 0; i < written; ++i) {
            const sm_status status = describe_edge(edges[static_cast<ULONG>(i)], face_index, out[i]);
            if (status != SM_SUCCESS) return status;
        }
        *count = edges.GetSize();
        return SM_SUCCESS;
    });
}

sm_status sm_brep_edge_length(const sm_brep* brep, size_t edge, double relative_accuracy, double* out) {
    if (!brep || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!std::isfinite(relative_accuracy)) return SM_ERR_INVALID_INPUT;
    const double accuracy = std::min(1e-1, std::max(1e-10, relative_accuracy));
    return guarded([&] {
        SmTArray<SmEdge*> edges;
        brep_of(brep)->GetEdges(edges);
        if (edge >= edges.GetSize()) return SM_SHIM_ERR_BAD_EDGE;
        double length = 0.0;
        const SmApiStatus status = SmApiEdgeComputeLength(edges[static_cast<ULONG>(edge)], accuracy, length);
        if (status != SM_SUCCESS) return static_cast<sm_status>(status);
        *out = length;
        return SM_SUCCESS;
    });
}

sm_status sm_brep_fillet_edges(sm_brep* brep, const uint32_t* edges, const double* radii, const int32_t* sections,
                               size_t count, uint32_t continuity, double thumbweight) {
    if (!brep || (count > 0 && (!edges || !radii || !sections))) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (continuity < 1 || continuity > 3 || !(thumbweight > 0.0) || !std::isfinite(thumbweight)) {
        return SM_SHIM_ERR_BAD_FILLET;
    }
    return guarded([&] {
        SmBrep* target = brep_of(brep);
        SmTArray<SmEdge*> picked;
        sm_status status = edges_at(*target, edges, count, picked);
        if (status != SM_SUCCESS) return status;
        SmTArray<double> radius_values;
        SmTArray<SmFilletSurfaceGeneratorType> section_types;
        for (size_t i = 0; i < count; ++i) {
            if (!(radii[i] > 0.0) || !std::isfinite(radii[i])) return SM_SHIM_ERR_BAD_FILLET;
            radius_values.Add(radii[i]);
            switch (sections[i]) {
                case SM_FILLET_LINEAR: section_types.Add(SM_FSG_LINEAR); break;
                case SM_FILLET_CIRCULAR: section_types.Add(SM_FSG_CIRCULAR); break;
                case SM_FILLET_BLEND: section_types.Add(SM_FSG_BLEND_CURVE); break;
                default: return SM_SHIM_ERR_BAD_FILLET;
            }
        }
        return keeping_solid(*target, [&] {
            return static_cast<sm_status>(
                SmApiFilletEdgesPerEdge(target, picked, radius_values, section_types, continuity, thumbweight));
        });
    });
}

sm_status sm_brep_chamfer_edges(sm_brep* brep, const uint32_t* edges, size_t count, double width) {
    if (!brep || (count > 0 && !edges)) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!(width > 0.0) || !std::isfinite(width)) return SM_SHIM_ERR_BAD_FILLET;
    return guarded([&] {
        SmBrep* target = brep_of(brep);
        SmTArray<SmEdge*> picked;
        const sm_status status = edges_at(*target, edges, count, picked);
        if (status != SM_SUCCESS) return status;
        if (picked.GetSize() == 0) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
        return keeping_solid(*target, [&] { return chamfer(*target, picked, width); });
    });
}

sm_status sm_brep_face_geometry(const sm_brep* brep, sm_face_geometry* out, size_t capacity, size_t* count) {
    if (!brep || !count || (capacity > 0 && !out)) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmTArray<SmFace*> faces;
        brep_of(brep)->GetFaces(faces);
        const size_t written = std::min(capacity, static_cast<size_t>(faces.GetSize()));
        for (size_t i = 0; i < written; ++i) {
            const sm_status status = describe_face_geometry(faces[static_cast<ULONG>(i)], out[i]);
            if (status != SM_SUCCESS) return status;
        }
        *count = faces.GetSize();
        return SM_SUCCESS;
    });
}

sm_status sm_brep_face_centroid(const sm_brep* brep, size_t face, double relative_accuracy, double* area,
                                double centroid[3]) {
    if (!brep || !area || !centroid) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!std::isfinite(relative_accuracy)) return SM_ERR_INVALID_INPUT;
    const double accuracy = std::min(1e-1, std::max(1e-10, relative_accuracy));
    return guarded([&] {
        const SmFace* target = face_at(*brep_of(brep), face);
        if (!target) return SM_SHIM_ERR_BAD_FACE;
        double value = 0.0;
        SmPoint3d point;
        const sm_status status = face_centroid(target, accuracy, value, point);
        if (status != SM_SUCCESS) return status;
        *area = value;
        put3(centroid, point);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_bounding_box(const sm_brep* brep, double min[3], double max[3]) {
    if (!brep || !min || !max) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmPoint3d low;
        SmPoint3d high;
        const SmApiStatus status = SmApiBrepBoundingBox(brep_of(brep), FALSE, low, high);
        if (status != SM_SUCCESS) return status;
        put3(min, low);
        put3(max, high);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_set_face_label(sm_brep* brep, size_t face, const sm_label* label) {
    if (!brep || !label) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (label->present && role_length(*label) == 0) return SM_SHIM_ERR_BAD_LABEL;
    return guarded([&] {
        SmFace* target = face_at(*brep_of(brep), face);
        if (!target) return SM_SHIM_ERR_BAD_FACE;
        write_label(*target, *label);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_set_edge_label(sm_brep* brep, size_t edge, const sm_label* label) {
    if (!brep || !label) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (label->present && role_length(*label) == 0) return SM_SHIM_ERR_BAD_LABEL;
    return guarded([&] {
        SmTArray<SmEdge*> edges;
        brep_of(brep)->GetEdges(edges);
        if (edge >= edges.GetSize()) return SM_SHIM_ERR_BAD_EDGE;
        write_label(*edges[static_cast<ULONG>(edge)], *label);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_place(sm_brep* brep, sm_vec3 origin, sm_vec3 x_axis, sm_vec3 y_axis) {
    if (!brep) return SM_SHIM_ERR_NULL_ARGUMENT;
    const double dot = x_axis.x * y_axis.x + x_axis.y * y_axis.y + x_axis.z * y_axis.z;
    if (!std::isfinite(origin.x) || !std::isfinite(origin.y) || !std::isfinite(origin.z) || !unit_vector(x_axis) ||
        !unit_vector(y_axis) || std::fabs(dot) > 1e-9) {
        return SM_SHIM_ERR_BAD_FRAME;
    }
    return guarded([&] {
        const SmAxis2Placement frame(SmPoint3d(origin.x, origin.y, origin.z), SmVector3d(x_axis.x, x_axis.y, x_axis.z),
                                     SmVector3d(y_axis.x, y_axis.y, y_axis.z));
        return static_cast<sm_status>(SmApiTransform(brep_of(brep), frame));
    });
}

sm_status sm_brep_transform(sm_brep* brep, const double rotation[9], sm_vec3 translation) {
    if (!brep || !rotation) return SM_SHIM_ERR_NULL_ARGUMENT;
    // R's columns are where the x, y and z axes go.
    const sm_vec3 x = {rotation[0], rotation[3], rotation[6]};
    const sm_vec3 y = {rotation[1], rotation[4], rotation[7]};
    const sm_vec3 z = {rotation[2], rotation[5], rotation[8]};
    const double dot = x.x * y.x + x.y * y.y + x.z * y.z;
    const sm_vec3 cross = {x.y * y.z - x.z * y.y, x.z * y.x - x.x * y.z, x.x * y.y - x.y * y.x};
    const double off = std::max({std::fabs(cross.x - z.x), std::fabs(cross.y - z.y), std::fabs(cross.z - z.z)});
    if (!std::isfinite(translation.x) || !std::isfinite(translation.y) || !std::isfinite(translation.z) ||
        !unit_vector(x) || !unit_vector(y) || std::fabs(dot) > 1e-9 || !(off <= 1e-9)) {
        return SM_SHIM_ERR_BAD_FRAME;
    }
    return guarded([&] {
        const SmAxis2Placement motion(SmPoint3d(translation.x, translation.y, translation.z),
                                      SmVector3d(x.x, x.y, x.z), SmVector3d(y.x, y.y, y.z));
        return static_cast<sm_status>(SmApiTransform(brep_of(brep), motion));
    });
}

sm_status sm_brep_scale(sm_brep* brep, sm_vec3 point, sm_vec3 factors) {
    if (!brep) return SM_SHIM_ERR_NULL_ARGUMENT;
    const auto positive = [](double f) { return std::isfinite(f) && f > SM_EFF_ZERO; };
    if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z) || !positive(factors.x) ||
        !positive(factors.y) || !positive(factors.z)) {
        return SM_SHIM_ERR_BAD_SCALE;
    }
    return guarded([&] {
        SmBrep& solid = *brep_of(brep);
        SmTArray<SmVertex*> vertices;
        SmTArray<SmEdge*> edges;
        SmTArray<SmFace*> faces;
        solid.GetVertices(vertices);
        solid.GetEdges(edges);
        solid.GetFaces(faces);
        std::vector<double> vertex_tolerances, edge_tolerances, face_tolerances;
        for (ULONG i = 0; i < vertices.GetSize(); ++i) vertex_tolerances.push_back(vertices[i]->GetTolerance().val);
        for (ULONG i = 0; i < edges.GetSize(); ++i) edge_tolerances.push_back(edges[i]->GetTolerance().val);
        for (ULONG i = 0; i < faces.GetSize(); ++i) face_tolerances.push_back(faces[i]->GetTolerance().val);
        const double brep_tolerance = solid.GetTolerance().val;
        // SmApiScaleByPt's composition: scale about the origin, then move the
        // origin's image so that the point stays put.
        SmAxis2Placement move;
        move.Translate(SmVector3d(point.x * (1.0 - factors.x), point.y * (1.0 - factors.y),
                                  point.z * (1.0 - factors.z)));
        const SmVector3d scale(factors.x, factors.y, factors.z);
        // Not the kernel's tolerance scaling, which resets every tolerance to
        // a band picked from the brep's new size (1e-5 to 1e-4 for a ×2 crate),
        // shrinking any a boolean grew, and moves the context's model size.
        const SmStatus status = solid.Transform(move, &scale, FALSE);
        if (status != SM_SUCCESS) return static_cast<sm_status>(status);
        const double by = std::max({factors.x, factors.y, factors.z});
        for (ULONG i = 0; i < vertices.GetSize(); ++i) {
            vertices[i]->SetTolerance(SmZoneTol3d(vertex_tolerances[i] * by), FALSE, FALSE);
        }
        for (ULONG i = 0; i < edges.GetSize(); ++i) {
            edges[i]->SetTolerance(SmZoneTol3d(edge_tolerances[i] * by), FALSE, FALSE);
        }
        for (ULONG i = 0; i < faces.GetSize(); ++i) {
            faces[i]->SetTolerance(SmZoneTol3d(face_tolerances[i] * by), FALSE, FALSE);
        }
        solid.SetTolerance(SmZoneTol3d(brep_tolerance * by), FALSE, FALSE);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_tolerance(const sm_brep* brep, double* out) {
    if (!brep || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        *out = brep_of(brep)->GetTolerance().val;
        return SM_SUCCESS;
    });
}

sm_status sm_brep_sweep_face(const sm_brep* brep, size_t face, double distance, sm_brep** out, sm_swept_face** faces,
                             size_t* face_count) {
    if (!brep || !out || !faces || !face_count) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!std::isfinite(distance) || distance == 0.0) return SM_SHIM_ERR_BAD_SWEEP;
    return guarded([&] {
        const SmBrep& source = *brep_of(brep);
        SmFace* target = face_at(source, face);
        if (!target) return SM_SHIM_ERR_BAD_FACE;
        SmBrep* prism = nullptr;
        std::vector<sm_swept_face> roles;
        const sm_status status = sweep_face(source, target, distance, prism, roles);
        if (status != SM_SUCCESS) return status;
        std::unique_ptr<SmBrep> owned(prism);
        std::unique_ptr<sm_swept_face[]> array(new sm_swept_face[roles.size()]);
        std::copy(roles.begin(), roles.end(), array.get());
        *out = handle_of(owned.release());
        *faces = array.release();
        *face_count = roles.size();
        return SM_SUCCESS;
    });
}

void sm_swept_faces_free(sm_swept_face* faces) { delete[] faces; }

sm_status sm_brep_merge_faces(sm_brep* brep) {
    if (!brep) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] { return static_cast<sm_status>(brep_of(brep)->RemoveTopologicalEdgesAndVertices()); });
}

sm_status sm_brep_combine(sm_brep* brep, const sm_brep* other) {
    if (!brep || !other) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (brep_of(brep) == brep_of(other)) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    // FALSE: keep every surface's type, as copied.
    return guarded(
        [&] { return static_cast<sm_status>(brep_of(brep)->MergeBrep(*brep_of(other), NULL, NULL, FALSE)); });
}

sm_status sm_brep_create_tapered_prism(const sm_vec3* points, size_t count, double height, double draft_deg,
                                       sm_brep** out) {
    if (!points || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (count < 3) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    return guarded([&] {
        SmTArray<SmCurve*> sides;
        struct Free {
            SmTArray<SmCurve*>& curves;
            ~Free() {
                for (ULONG i = 0; i < curves.GetSize(); ++i) delete curves[i];
            }
        } free_sides{sides};
        for (size_t i = 0; i < count; ++i) {
            const sm_vec3& a = points[i];
            const sm_vec3& b = points[(i + 1) % count];
            SmLine* side = nullptr;
            const SmApiStatus status = SmApiCreateLineSegment(SmPoint3d(a.x, a.y, a.z), SmPoint3d(b.x, b.y, b.z), side);
            if (status != SM_SUCCESS) return status;
            sides.Add(side);
        }
        SmBrep* brep = nullptr;
        const SmApiStatus status = SmApiCreateTaperExtrude(sides, height, draft_deg, 3, FALSE, brep);
        if (status != SM_SUCCESS) return status;
        *out = handle_of(brep);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_create_planar_faces(const sm_plane_curve* curves, size_t count, sm_vec3 origin, sm_vec3 x_axis,
                                      sm_vec3 y_axis, sm_brep** out) {
    if (!curves || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    const double dot = x_axis.x * y_axis.x + x_axis.y * y_axis.y + x_axis.z * y_axis.z;
    if (!std::isfinite(origin.x) || !std::isfinite(origin.y) || !std::isfinite(origin.z) || !unit_vector(x_axis) ||
        !unit_vector(y_axis) || std::fabs(dot) > 1e-9) {
        return SM_SHIM_ERR_BAD_FRAME;
    }
    if (count == 0) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    return guarded([&] {
        SmBrep* raw = nullptr;
        const sm_status status = planar_faces(curves, count, raw);
        if (status != SM_SUCCESS) return status;
        std::unique_ptr<SmBrep> brep(raw);
        const SmAxis2Placement frame(SmPoint3d(origin.x, origin.y, origin.z), SmVector3d(x_axis.x, x_axis.y, x_axis.z),
                                     SmVector3d(y_axis.x, y_axis.y, y_axis.z));
        const SmApiStatus placed = SmApiTransform(brep.get(), frame);
        if (placed != SM_SUCCESS) return static_cast<sm_status>(placed);
        *out = handle_of(brep.release());
        return SM_SUCCESS;
    });
}

sm_status sm_brep_create_cone(sm_vec3 base_center, double base_radius, double top_radius, double height,
                              sm_brep** out) {
    if (!out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* brep = nullptr;
        const SmApiStatus status = SmApiCreateCone(SmVector3d(base_center.x, base_center.y, base_center.z),
                                                   base_radius, top_radius, height, brep);
        if (status != SM_SUCCESS) return status;
        *out = handle_of(brep);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_create_torus(sm_vec3 center, double major_radius, double minor_radius, sm_brep** out) {
    if (!out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* brep = nullptr;
        const SmApiStatus status =
            SmApiCreateTorus(SmVector3d(center.x, center.y, center.z), major_radius, minor_radius, brep);
        if (status != SM_SUCCESS) return status;
        *out = handle_of(brep);
        return SM_SUCCESS;
    });
}

// The faces at `faces`, each once; false if one is out of range or repeated.
bool picked_faces(const SmBrep& brep, const uint32_t* faces, size_t count, SmTArray<SmFace*>& picked) {
    SmTArray<SmFace*> all;
    brep.GetFaces(all);
    std::vector<bool> taken(all.GetSize(), false);
    for (size_t i = 0; i < count; ++i) {
        if (faces[i] >= all.GetSize() || taken[faces[i]]) return false;
        taken[faces[i]] = true;
        picked.Add(all[faces[i]]);
    }
    return true;
}

sm_status sm_brep_offset_faces(sm_brep* brep, const uint32_t* faces, size_t count, double distance) {
    if (!brep || (count > 0 && !faces)) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!std::isfinite(distance) || count == 0) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    return guarded([&] {
        SmBrep* target = brep_of(brep);
        SmTArray<SmFace*> picked;
        if (!picked_faces(*target, faces, count, picked)) return SM_SHIM_ERR_BAD_FACE;
        if (distance == 0.0) return SM_SUCCESS;
        return offset_through(*target, picked, distance);
    });
}

sm_status sm_brep_offset_faces_unchecked(sm_brep* brep, const uint32_t* faces, size_t count, double distance) {
    if (!brep || (count > 0 && !faces)) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!std::isfinite(distance) || count == 0) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    return guarded([&] {
        SmBrep* target = brep_of(brep);
        SmTArray<SmFace*> picked;
        if (!picked_faces(*target, faces, count, picked)) return SM_SHIM_ERR_BAD_FACE;
        if (distance == 0.0) return SM_SUCCESS;
        return offset_faces(*target, picked, distance);
    });
}

sm_status sm_brep_offset_obstacles(const sm_brep* brep, const uint32_t* faces, size_t count, double distance,
                                   uint32_t* out_faces, uint32_t* out_obstacles, size_t capacity,
                                   size_t* out_count) {
    if (!brep || !out_count || (count > 0 && !faces) || (capacity > 0 && (!out_faces || !out_obstacles))) {
        return SM_SHIM_ERR_NULL_ARGUMENT;
    }
    if (!std::isfinite(distance) || count == 0) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    return guarded([&] {
        const SmBrep& source = *brep_of(brep);
        SmTArray<SmFace*> picked;
        if (!picked_faces(source, faces, count, picked)) return SM_SHIM_ERR_BAD_FACE;
        std::vector<Obstacle> obstacles;
        Motions motions;
        if (distance != 0.0 && plan_offset(picked, distance, source.GetTolerance(), motions) == SM_SUCCESS) {
            find_obstacles(source, picked, distance, motions, obstacles);
        }
        const FaceIndex index = index_faces(source);
        size_t written = 0;
        for (size_t o = 0; o < obstacles.size(); ++o) {
            for (const SmFace* face : obstacles[o].faces) {
                if (written < capacity) {
                    out_faces[written] = static_cast<uint32_t>(face_id(index, face));
                    out_obstacles[written] = static_cast<uint32_t>(o);
                }
                ++written;
            }
        }
        *out_count = written;
        return SM_SUCCESS;
    });
}

sm_status sm_brep_move_faces(sm_brep* brep, const uint32_t* faces, size_t count, sm_vec3 translation,
                             sm_vec3 axis_point, sm_vec3 axis_direction, double angle_deg, uint8_t* out_roles,
                             size_t role_count) {
    if (!brep || !out_roles || (count > 0 && !faces)) return SM_SHIM_ERR_NULL_ARGUMENT;
    const auto finite = [](sm_vec3 v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); };
    if (count == 0 || !finite(translation) || !std::isfinite(angle_deg)) {
        return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    }
    const SmVector3d axis(axis_direction.x, axis_direction.y, axis_direction.z);
    if (angle_deg != 0.0 && (!finite(axis_point) || !finite(axis_direction) || !(axis.Length() > 0.0))) {
        return SM_SHIM_ERR_BAD_AXIS;
    }
    return guarded([&] {
        SmBrep* target = brep_of(brep);
        SmTArray<SmFace*> picked;
        if (!picked_faces(*target, faces, count, picked)) return SM_SHIM_ERR_BAD_FACE;
        if (role_count != target->GetNumFaces()) return SM_SHIM_ERR_BAD_FACE;
        Motion rigid;
        rigid.kind = Motion::kShift;
        rigid.shift = SmVector3d(translation.x, translation.y, translation.z);
        if (angle_deg != 0.0) {
            // Rodrigues: a turn by the angle about the axis through the point,
            // anticlockwise looking down the axis from its tip, then the
            // translation: p goes to R (p - c) + c + t.
            const SmVector3d k = axis / axis.Length();
            const double theta = angle_deg * SM_PI / 180.0;
            const double c = std::cos(theta);
            const double s = std::sin(theta);
            const double v = 1.0 - c;
            const double r[3][3] = {{c + k.x * k.x * v, k.x * k.y * v - k.z * s, k.x * k.z * v + k.y * s},
                                    {k.y * k.x * v + k.z * s, c + k.y * k.y * v, k.y * k.z * v - k.x * s},
                                    {k.z * k.x * v - k.y * s, k.z * k.y * v + k.x * s, c + k.z * k.z * v}};
            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < 3; ++col) rigid.turn[row][col] = r[row][col];
            }
            rigid.turned = true;
            const SmPoint3d centre(axis_point.x, axis_point.y, axis_point.z);
            const SmVector3d turned = turned_vector(rigid, SmVector3d(centre.x, centre.y, centre.z));
            rigid.shift = rigid.shift + SmVector3d(centre.x - turned.x, centre.y - turned.y, centre.z - turned.z);
        }
        if (!rigid.turned && rigid.shift.Length() == 0.0) {
            // Nothing moves, and nothing changes.
            SmTArray<SmFace*> all;
            target->GetFaces(all);
            for (ULONG i = 0; i < all.GetSize(); ++i) {
                ULONG index = 0;
                out_roles[i] = picked.FindElement(all[i], index) ? SM_MOVE_MOVED : SM_MOVE_KEPT;
            }
            return SM_SUCCESS;
        }
        // A face beside another piece of its own surface (an imprinted
        // panel) can't stretch that piece to meet it.
        if (needs_steps(picked, target->GetTolerance())) return SM_SHIM_ERR_UNSUPPORTED_STEP;
        Motions motions;
        const SmVector3d moved_by(translation.x, translation.y, translation.z);
        sm_status status = plan_move(picked, rigid, moved_by, target->GetTolerance(), motions);
        if (status != SM_SUCCESS) return status;
        const std::vector<uint8_t> roles = move_roles(*target, picked, motions);
        status = keeping_topology(*target, [&] { return move_planned(*target, motions); });
        if (status != SM_SUCCESS) return status;
        std::copy(roles.begin(), roles.end(), out_roles);
        return SM_SUCCESS;
    });
}

sm_status sm_brep_tessellate(sm_brep* brep, const sm_tessellation_params* params, sm_mesh** out) {
    if (!brep || !params || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* source = brep_of(brep);
        SmTessellationParams quality;
        quality.dChordHeightTol = params->chord_height_tol;
        quality.dCurveAngleTolDeg = params->curve_angle_tol_deg;
        quality.dSurfaceAngleTolDeg = params->surface_angle_tol_deg;
        quality.dMaxEdgeLength = params->max_edge_length;
        quality.dMaxAspectRatio = params->max_aspect_ratio;

        SmPolyBrep* poly = nullptr;
        const SmApiStatus status = SmApiTessellate(source, poly, quality);
        if (status != SM_SUCCESS) return status;
        if (!poly) return SM_ERR;
        // Freed as SmApiTessellate frees a rejected mesh: its back pointers
        // label the source brep's faces, which it doesn't own.
        const std::unique_ptr<SmPolyBrep, PolyBrepDeleter> owned_poly(poly);

        // Faces are identified by their index in the brep's own list.
        SmTArray<SmFace*> faces;
        source->GetFaces(faces);
        std::unordered_map<const SmFace*, uint32_t> face_index;
        auto mesh = std::make_unique<sm_mesh>();
        mesh->faces.reserve(faces.GetSize());
        for (ULONG i = 0; i < faces.GetSize(); ++i) {
            face_index.emplace(faces[i], static_cast<uint32_t>(i));
            mesh->faces.push_back(describe_face(faces[i]));
        }

        sm_status result = build_triangles(*poly, face_index, *mesh);
        if (result != SM_SUCCESS) return result;
        result = build_edges(*source, params->edge_angle_tol_deg, face_index, *mesh);
        if (result != SM_SUCCESS) return result;
        *out = mesh.release();
        return SM_SUCCESS;
    });
}

void sm_mesh_view_get(const sm_mesh* mesh, sm_mesh_view* out) {
    if (!out) return;
    std::memset(out, 0, sizeof *out);
    if (!mesh) return;
    out->positions = mesh->positions.data();
    out->normals = mesh->normals.data();
    out->surface_uvs = mesh->surface_uvs.data();
    out->vertex_count = mesh->positions.size() / 3;
    out->indices = mesh->indices.data();
    out->face_ids = mesh->face_ids.data();
    out->triangle_count = mesh->face_ids.size();
    out->faces = mesh->faces.data();
    out->face_count = mesh->faces.size();
    out->edges = mesh->edges.data();
    out->edge_count = mesh->edges.size();
    out->edge_points = mesh->edge_points.data();
    out->edge_point_count = mesh->edge_points.size() / 3;
    out->polygons_split = mesh->polygons_split;
    out->corners_missing_uv = mesh->corners_missing_uv;
}

void sm_mesh_free(sm_mesh* mesh) { delete mesh; }

}  // extern "C"

// Non-uniform scales that keep surfaces and curves analytic where the scale
// leaves them so (sm_brep_scale_keeping).
//
// SmBrep::Transform turns every surface and curve into a NURBS before a
// non-uniform scale (TurnToNURBS), since SMLib's analytic surfaces and curves
// refuse one. They are B-splines too, though: each keeps a NURBS of itself,
// and its own Transform scales that NURBS and derives its analytic description
// from it again (SmPlane::Transform), refusing only a scale that isn't
// uniform. So this scales each face's and edge's NURBS exactly, which keeps
// its knots and weights and so its parameters, and makes it analytic again
// where the scaled NURBS still is: a plane always, a cylinder or cone along a
// scale axis whose other two factors are equal, a circle in a plane the scale
// keeps round. A circle the scale stretches becomes an ellipse, and with
// extrusions a cylinder stretched across its axis becomes an ellipse swept
// along a line (SmSurfOfExtrusion). A candidate is taken only when it has the
// scaled NURBS's domain and lies on it at the same parameters, so every edge
// keeps its interval and every UV trim curve stays right. Replacing a face's
// surface or an edge's curve drops the trims, so they are put back after.
//
// A reflection (sm_brep_mirror) is the same scale with a factor of -1 and a
// turn after it. It keeps every kind, but it turns each surface's (u, v)
// round, so the normal Su x Sv now points into the solid; the faces' sides
// and the order of faces round each edge flip as SmBrep::Mirror flips them.

namespace {

// A scale about the origin by factors, then the rigid motion move.
struct Scaling {
    SmVector3d factors;
    SmAxis2Placement move;

    SmPoint3d point(const SmPoint3d& p) const {
        SmPoint3d out;
        move.TransformPoint(p.Multiply(factors), out);
        return out;
    }
    SmVector3d vector(const SmVector3d& v) const {
        SmVector3d out;
        move.TransformVector(v.Multiply(factors), out);
        return out;
    }
    double largest() const { return std::max({std::fabs(factors.x), std::fabs(factors.y), std::fabs(factors.z)}); }
    bool reflects() const { return factors.x * factors.y * factors.z < 0.0; }
    // Equal factors but for their signs: angles and shapes stay.
    bool similar() const {
        const double by = largest();
        return SM_ARE_SAME(std::fabs(factors.x), by) && SM_ARE_SAME(std::fabs(factors.y), by) &&
               SM_ARE_SAME(std::fabs(factors.z), by);
    }
    bool mirror() const { return reflects() && similar() && SM_ARE_SAME(largest(), 1.0); }

    // For a mirror, the rigid motion that moves everything in the plane
    // through `on` with normal `normal` exactly as the mirror does: the
    // mirror after a reflection across that plane, which leaves the plane
    // where it is. Planar surfaces and curves take it with their own
    // Transform, which keeps their kind, their NURBS and their parameters.
    SmAxis2Placement turn_for(const SmPoint3d& on, const SmVector3d& normal) const {
        SmVector3d n = normal;
        n.Unitize();
        const auto across = [&](const SmVector3d& v) { return vector(v - n * (2.0 * v.Dot(n))); };
        // SmPoint3d is SmVector3d.
        const SmPoint3d origin = point(n * (2.0 * on.Dot(n)));
        return SmAxis2Placement(origin, across(SmVector3d(1.0, 0.0, 0.0)), across(SmVector3d(0.0, 1.0, 0.0)));
    }
};

// The largest distance between two surfaces at the same (u, v), over a grid
// on the domain.
double surface_gap(const SmSurface& a, const SmSurface& b, const SmExtent2d& domain) {
    double gap = 0.0;
    for (int i = 0; i <= 4; ++i) {
        for (int j = 0; j <= 4; ++j) {
            const SmPoint2d uv = domain.Evaluate(i / 4.0, j / 4.0);
            SmPoint3d p, q;
            if (a.EvaluatePoint(uv, p) != SM_SUCCESS || b.EvaluatePoint(uv, q) != SM_SUCCESS) return HUGE_VAL;
            gap = std::max(gap, (p - q).Length());
        }
    }
    return gap;
}

// The largest distance between two curves at the same parameter.
double curve_gap(const SmCurve& a, const SmCurve& b, const SmExtent1d& interval) {
    double gap = 0.0;
    for (int i = 0; i <= 8; ++i) {
        const double t = interval.Evaluate(i / 8.0);
        SmPoint3d p, q;
        if (a.EvaluatePoint(t, p) != SM_SUCCESS || b.EvaluatePoint(t, q) != SM_SUCCESS) return HUGE_VAL;
        gap = std::max(gap, (p - q).Length());
    }
    return gap;
}

// Takes the candidate if it stands in for the scaled NURBS exactly; deletes it
// otherwise.
template <typename T>
T* if_same_surface(T* candidate, const SmBSplineSurface& scaled, double tolerance) {
    const SmExtent2d domain = scaled.GetNaturalUVDomain();
    if (candidate && candidate->GetNaturalUVDomain() == domain &&
        surface_gap(*candidate, scaled, domain) <= tolerance) {
        return candidate;
    }
    delete candidate;
    return nullptr;
}

template <typename T>
T* if_same_curve(T* candidate, const SmBSplineCurve& scaled, double tolerance) {
    const SmExtent1d interval = scaled.GetNaturalInterval();
    if (candidate && candidate->GetNaturalInterval() == interval &&
        curve_gap(*candidate, scaled, interval) <= tolerance) {
        return candidate;
    }
    delete candidate;
    return nullptr;
}

// Whether a curve runs along the ellipse centred on `centre` with unit axes
// x and y and radii rx and ry, from the angle step.min to step.max in
// degrees, as SmEllipse describes one.
bool runs_on_ellipse(const SmBSplineCurve& curve, const SmPoint3d& centre, const SmVector3d& x, const SmVector3d& y,
                     double rx, double ry, const SmExtent1d& step, double tolerance) {
    const auto at = [&](double degrees) {
        const double angle = degrees * SM_PI / 180.0;
        return centre + x * (rx * std::cos(angle)) + y * (ry * std::sin(angle));
    };
    const SmExtent1d interval = curve.GetNaturalInterval();
    SmPoint3d start, end;
    if (curve.EvaluatePoint(interval.GetMin(), start) != SM_SUCCESS ||
        curve.EvaluatePoint(interval.GetMax(), end) != SM_SUCCESS || (start - at(step.GetMin())).Length() > tolerance ||
        (end - at(step.GetMax())).Length() > tolerance) {
        return false;
    }
    const SmVector3d z = x * y;
    for (int i = 1; i < 8; ++i) {
        SmPoint3d p;
        if (curve.EvaluatePoint(interval.Evaluate(i / 8.0), p) != SM_SUCCESS) return false;
        const SmVector3d d = p - centre;
        const double u = d.Dot(x) / rx;
        const double v = d.Dot(y) / ry;
        if (std::fabs(std::sqrt(u * u + v * v) - 1.0) * std::max(rx, ry) > tolerance ||
            std::fabs(d.Dot(z)) > tolerance) {
            return false;
        }
    }
    return true;
}

// The ellipse the scale takes an ellipse or circle to, over the scaled NURBS
// of it, which it keeps as its own, as SmCircle::IsNurbCurveCircle keeps the
// NURBS it recognises: when the images of its axes stay square to each other
// and the NURBS runs along it. NULL otherwise.
SmEllipse* scaled_ellipse(const SmContext& context, const SmEllipse& ellipse, const SmBSplineCurve& scaled,
                          const Scaling& scaling, double tolerance) {
    if (ellipse.GetInsideOut()) return nullptr;
    const SmAxis2Placement& position = ellipse.GetPosition();
    const SmVector3d x = scaling.vector(position.GetXAxisRef() * ellipse.GetXRadius());
    const SmVector3d y = scaling.vector(position.GetYAxisRef() * ellipse.GetYRadius());
    const double rx = x.Length();
    const double ry = y.Length();
    if (rx <= 0.0 || ry <= 0.0 || std::fabs(x.Dot(y)) > 1e-12 * rx * ry) return nullptr;
    const SmPoint3d centre = scaling.point(position.GetOriginRef());
    const SmExtent1d step = ellipse.GetSTEPInterval();
    if (!runs_on_ellipse(scaled, centre, x / rx, y / ry, rx, ry, step, tolerance)) return nullptr;
    return new (context) SmEllipse(centre, x / rx, y / ry, step, rx, ry, 3, &context, &scaled);
}

// A cylinder stretched across its axis as an ellipse swept along a line, with
// the scaled NURBS's parameters (as SmSurfOfExtrusion::IsNurbSurfaceSurfOfExtrusion
// builds one, but from the ellipse the cylinder's circles scale to rather than
// from a NURBS of it); NULL when it isn't one.
SmSurfOfExtrusion* ellipse_extrusion(const SmContext& context, const SmCone& cylinder, const SmBSplineSurface& scaled,
                                     const Scaling& scaling, double tolerance) {
    if (!cylinder.IsCylinder()) return nullptr;
    const ULONG degree_u = scaled.GetDegree(SM_SP_U);
    const ULONG degree_v = scaled.GetDegree(SM_SP_V);
    if ((degree_u == 1) == (degree_v == 1)) return nullptr;
    const SmExtent2d domain = scaled.GetNaturalUVDomain();
    // The generator is the iso curve at the start of the linear direction.
    const bool linear_u = degree_u == 1;
    const SmSurfParamType across = linear_u ? SM_SP_U : SM_SP_V;
    const double start = linear_u ? domain.GetMin().x : domain.GetMin().y;
    const double end = linear_u ? domain.GetMax().x : domain.GetMax().y;
    SmPoint3d origin, line_end;
    scaled.EvaluatePoint(domain.GetMin(), origin);
    scaled.EvaluatePoint(domain.Evaluate(linear_u ? 1.0 : 0.0, linear_u ? 0.0 : 1.0), line_end);
    SmVector3d sweep = (line_end - origin) / (end - start);

    // The cylinder's circle there, exactly, and the ellipse it scales to.
    SmBSplineCurve* raw = nullptr;
    if (cylinder.CreateIsoParametricCurve(context, across, start, 0.0, raw) != SM_SUCCESS || !raw) return nullptr;
    const std::unique_ptr<SmBSplineCurve> iso(raw);
    SmCircle* circle_raw = nullptr;
    if (!SmCircle::IsNurbCurveCircle(context, iso.get(), circle_raw)) return nullptr;
    const std::unique_ptr<SmCircle> circle(if_same_curve(circle_raw, *iso, tolerance));
    if (!circle) return nullptr;
    raw = nullptr;
    if (scaled.CreateIsoParametricCurve(context, across, start, 0.0, raw) != SM_SUCCESS || !raw) return nullptr;
    const std::unique_ptr<SmBSplineCurve> scaled_iso(raw);
    std::unique_ptr<SmEllipse> ellipse(scaled_ellipse(context, *circle, *scaled_iso, scaling, tolerance));
    if (!ellipse) return nullptr;

    // As the recogniser does, so that the analytic and NURBS domains agree.
    if (!SM_IS_ZERO(start)) {
        SmAxis2Placement shift;
        shift.Translate(-start * sweep);
        ellipse->Transform(shift);
        origin -= start * sweep;
    }
    const SmExtent1d profile = ellipse->GetSTEPInterval();
    const SmExtent2d analytic(SmPoint2d(profile.GetMin(), start), SmPoint2d(profile.GetMax(), end));
    SmSurfOfExtrusion* extrusion =
        new (context) SmSurfOfExtrusion(ellipse.release(), origin, sweep, analytic, linear_u ? TRUE : FALSE);
    if (extrusion->Reparameterize(domain) != SM_SUCCESS) {
        delete extrusion;
        return nullptr;
    }
    return if_same_surface(extrusion, scaled, tolerance);
}

// The torus a reflection (with equal factors otherwise) takes a torus to,
// from its own description, since SmTorus::IsNurbSurfaceTorus can't find one
// in a NURBS: centre and axes reflected, radii scaled, and the minor circle
// running the other way round its new axis (SmTorus keeps that as its
// generating circle's inside-out flag). NULL for any other scale, or when it
// doesn't lie on the scaled NURBS at the same parameters.
SmTorus* reflected_torus(const SmContext& context, const SmTorus& original, const SmBSplineSurface& scaled,
                         const Scaling& scaling, double tolerance) {
    const double by = scaling.largest();
    if (!scaling.reflects() || !scaling.similar()) return nullptr;
    const SmEllipse* circle = dynamic_cast<const SmEllipse*>(original.GetGenCurve());
    if (!circle) return nullptr;
    SmAxis2Placement frame;
    double major = 0.0, minor = 0.0;
    SmBoolean swap = FALSE;
    if (original.GetCanonical(frame, major, minor, &swap) != SM_SUCCESS) return nullptr;
    SmVector3d x = scaling.vector(frame.GetXAxisRef());
    SmVector3d y = scaling.vector(frame.GetYAxisRef());
    x.Unitize();
    y.Unitize();
    SmTorus* torus = new (context)
        SmTorus(scaling.point(frame.GetOriginRef()), x, y, original.GetSTEPUVDomain(), major * by, minor * by, swap,
                circle->GetInsideOut() ? FALSE : TRUE, nullptr, &context);
    if (torus->Reparameterize(scaled.GetNaturalUVDomain()) != SM_SUCCESS) {
        delete torus;
        return nullptr;
    }
    return if_same_surface(torus, scaled, tolerance);
}

SmCurve* analytic_curve(const SmContext& context, const SmEllipse& original, const SmBSplineCurve& scaled,
                        const Scaling& scaling, double tolerance, bool& stretched);

// A copy of a curve, reflected (or scaled) exactly: a circle or ellipse that
// stays one (under a mirror, moved by the turn that moves its plane), and
// anything else as its NURBS.
SmBSplineCurve* scaled_copy(const SmContext& context, const SmBSplineCurve& curve, const Scaling& scaling,
                            double tolerance) {
    if (scaling.mirror() && curve.IsKindOf(SmEllipse_TYPE)) {
        SmCurve* copy = nullptr;
        if (curve.Copy(context, copy) != SM_SUCCESS || !copy) return nullptr;
        const SmAxis2Placement& position = static_cast<const SmEllipse&>(curve).GetPosition();
        if (copy->Transform(scaling.turn_for(position.GetOriginRef(), position.GetZAxis())) != SM_SUCCESS) {
            delete copy;
            return nullptr;
        }
        return static_cast<SmBSplineCurve*>(copy);
    }
    std::unique_ptr<SmBSplineCurve> scaled(new (context) SmBSplineCurve(curve));
    if (scaled->Transform(scaling.move, &scaling.factors) != SM_SUCCESS) return nullptr;
    if (curve.IsKindOf(SmEllipse_TYPE)) {
        bool stretched = false;
        SmCurve* analytic =
            analytic_curve(context, static_cast<const SmEllipse&>(curve), *scaled, scaling, tolerance, stretched);
        if (analytic) return static_cast<SmBSplineCurve*>(analytic);
    }
    return scaled.release();
}

// The surface of revolution a reflection (with equal factors otherwise)
// takes one to, from its own description: its generating curve reflected and
// revolved the other way round the reflected axis, which keeps the angle
// running as before. SMLib makes a fillet on a circular edge as one. NULL for
// any other scale, or when it doesn't lie on the scaled NURBS at the same
// parameters.
SmSurfOfRevolution* reflected_revolution(const SmContext& context, const SmSurfOfRevolution& original,
                                         const SmBSplineSurface& scaled, const Scaling& scaling, double tolerance) {
    if (!scaling.reflects() || !scaling.similar()) return nullptr;
    const SmBSplineCurve* generator = original.GetGenCurve();
    if (!generator) return nullptr;
    SmBSplineCurve* reflected = scaled_copy(context, *generator, scaling, tolerance);
    if (!reflected) return nullptr;
    const SmAxis2Placement& position = original.GetPosition();
    const SmVector3d axis = scaling.vector(position.GetXAxisRef()) * scaling.vector(position.GetYAxisRef());
    SmSurfOfRevolution* surface = nullptr;
    // It owns the generating curve from here, even when it fails.
    if (SmSurfOfRevolution::CreateCanonical(context, reflected, scaling.point(position.GetOriginRef()), axis, surface,
                                            original.GetSwapUV()) != SM_SUCCESS ||
        !surface) {
        delete surface;
        return nullptr;
    }
    if (surface->AdjustSTEPUVDomain(original.GetSTEPUVDomain()) != SM_SUCCESS ||
        surface->Reparameterize(scaled.GetNaturalUVDomain()) != SM_SUCCESS) {
        delete surface;
        return nullptr;
    }
    return if_same_surface(surface, scaled, tolerance);
}

// The analytic surface the scaled NURBS of `original` still is, of
// `original`'s own kind, or with extrusions an ellipse swept along a line;
// NULL when it is neither.
SmSurface* analytic_surface(const SmContext& context, const SmSurface& original, const SmBSplineSurface& scaled,
                            const Scaling& scaling, bool extrusions, double tolerance, bool& extruded) {
    extruded = false;
    const SM_TYPE type = original.GetType();
    if (type == SmPlane_TYPE) {
        // SmPlane::Transform's ladder of tolerances.
        for (double scale : {1.0, 10.0, 100.0, 1000.0}) {
            SmPlane* plane = nullptr;
            if (SmPlane::IsNurbSurfacePlane(context, &scaled, plane, scale)) {
                return if_same_surface(plane, scaled, tolerance);
            }
        }
        return nullptr;
    }
    if (type == SmCone_TYPE || type == SmCylinder_TYPE) {
        SmCone* cone = nullptr;
        if (SmCone::IsNurbSurfaceCone(context, &scaled, cone)) {
            if (SmCone* same = if_same_surface(cone, scaled, tolerance)) return same;
        }
        if (!extrusions) return nullptr;
        SmSurfOfExtrusion* extrusion =
            ellipse_extrusion(context, static_cast<const SmCone&>(original), scaled, scaling, tolerance);
        extruded = extrusion != nullptr;
        return extrusion;
    }
    if (type == SmSphere_TYPE) {
        SmSphere* sphere = nullptr;
        if (SmSphere::IsNurbSurfaceSphere(context, &scaled, sphere)) return if_same_surface(sphere, scaled, tolerance);
        return nullptr;
    }
    if (type == SmSurfOfRevolution_TYPE) {
        return reflected_revolution(context, static_cast<const SmSurfOfRevolution&>(original), scaled, scaling,
                                    tolerance);
    }
    if (type == SmTorus_TYPE) {
        if (SmTorus* torus = reflected_torus(context, static_cast<const SmTorus&>(original), scaled, scaling, tolerance)) {
            return torus;
        }
        SmTorus* torus = nullptr;
        if (SmTorus::IsNurbSurfaceTorus(context, &scaled, torus)) return if_same_surface(torus, scaled, tolerance);
    }
    return nullptr;
}

// The analytic curve the scaled NURBS of an ellipse or circle still is: a
// circle, or the ellipse its axes scale to. NULL when it is neither.
SmCurve* analytic_curve(const SmContext& context, const SmEllipse& original, const SmBSplineCurve& scaled,
                        const Scaling& scaling, double tolerance, bool& stretched) {
    stretched = false;
    SmCircle* circle = nullptr;
    if (SmCircle::IsNurbCurveCircle(context, &scaled, circle)) {
        if (SmCircle* same = if_same_curve(circle, scaled, tolerance)) return same;
    }
    SmEllipse* ellipse =
        if_same_curve(scaled_ellipse(context, original, scaled, scaling, tolerance), scaled, tolerance);
    stretched = ellipse != nullptr;
    return ellipse;
}

// A reflection reverses the order of the faces round an edge, so each pair
// of edgeuses (one face's two sides) swaps places, as in SmBrep::Mirror.
void reverse_radial_order(SmEdge& edge) {
    SmTArray<SmEdgeuse*> uses;
    edge.GetEdgeuses(uses);
    for (ULONG i = 0; i < uses.GetSize(); ++i) edge.Remove(uses[i]);
    for (ULONG i = 0; i + 1 < uses.GetSize(); i += 2) {
        edge.PostInsert(uses[i + 1]);
        edge.PostInsert(uses[i]);
    }
}

// Where an object sits in its owner's list: first, or after another.
struct ListPlace {
    bool first;
    SmTopology* after;
};

ListPlace place_of(const SmOwningTopology& owner, const SmTopology* object) {
    return {owner.GetList() == object || owner.GetSize() == 1, object->GetLast()};
}

void put_at(SmOwningTopology& owner, SmTopology* object, const ListPlace& place) {
    if (place.first) {
        owner.PreInsert(object);
    } else {
        owner.InsertAfter(object, place.after);
    }
}

// A reflection turns a face's normal to its other side, so its two faceuses
// trade shells, as in SmBrep::Mirror, each taking the other's place in its
// shell's list so that a second reflection gives back the same order.
void swap_sides(SmFace& face) {
    SmFaceuse *first = nullptr, *second = nullptr;
    face.GetFaceuses(first, second);
    if (!first || !second) return;
    SmShell* a = first->GetShell();
    SmShell* b = second->GetShell();
    if (a == b) return;
    const ListPlace in_a = place_of(*a, first);
    const ListPlace in_b = place_of(*b, second);
    a->Remove(first);
    b->Remove(second);
    put_at(*a, second, in_a);
    put_at(*b, first, in_b);
}

sm_status scale_keeping(SmBrep& brep, const Scaling& scaling, bool extrusions, sm_scale_report& report) {
    const SmContext& context = *brep.GetContext();
    // Analytic candidates must lie on the scaled NURBS far inside the
    // tolerance the scale leaves the brep, so they change nothing.
    const double by = scaling.largest();
    const double tolerance = 1e-3 * brep.GetTolerance().val * by;
    // Every analytic surface and curve scales itself uniformly, and keeps
    // its kind; SMLib can't recognise a torus in a NURBS.
    const bool uniform = SM_ARE_SAME(scaling.factors.x, scaling.factors.y) &&
                         SM_ARE_SAME(scaling.factors.x, scaling.factors.z);
    SmTArray<SmVertex*> vertices;
    SmTArray<SmEdge*> edges;
    SmTArray<SmFace*> faces;
    brep.GetVertices(vertices);
    brep.GetEdges(edges);
    brep.GetFaces(faces);
    // Only what this knows how to scale exactly, each surface and curve once.
    SmTArray<SmSurface*> surfaces;
    SmTArray<SmCurve*> curves;
    brep.GetSurfaces(surfaces);
    brep.GetCurves(curves);
    if (surfaces.GetSize() != faces.GetSize() || curves.GetSize() != edges.GetSize()) return SM_SHIM_ERR_BAD_SCALE;
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        const SmCurve* curve = edges[i]->GetCurve();
        if (!curve || !curve->IsKindOf(SmBSplineCurve_TYPE)) return SM_SHIM_ERR_BAD_SCALE;
    }
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const SmSurface* surface = faces[i]->GetSurface();
        if (!surface || !surface->IsKindOf(SmBSplineSurface_TYPE)) return SM_SHIM_ERR_BAD_SCALE;
    }
    std::vector<double> vertex_tolerances, edge_tolerances, face_tolerances;
    for (ULONG i = 0; i < vertices.GetSize(); ++i) vertex_tolerances.push_back(vertices[i]->GetTolerance().val);
    for (ULONG i = 0; i < edges.GetSize(); ++i) edge_tolerances.push_back(edges[i]->GetTolerance().val);
    for (ULONG i = 0; i < faces.GetSize(); ++i) face_tolerances.push_back(faces[i]->GetTolerance().val);
    const double brep_tolerance = brep.GetTolerance().val;
    std::vector<SmBoolean> rectangular;
    for (ULONG i = 0; i < faces.GetSize(); ++i) rectangular.push_back(faces[i]->GetRectangularTrim());

    // Every edgeuse pair's UV trim curve, held by one of the two.
    std::vector<std::pair<SmEdgeuse*, std::unique_ptr<SmBSplineCurve>>> trims;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        SmTArray<SmEdgeuse*> uses;
        faces[i]->GetUpwardEdgeuses(uses);
        for (ULONG j = 0; j < uses.GetSize(); ++j) {
            SmEdgeuse* owner = nullptr;
            const SmBSplineCurve* trim = uses[j]->GetUVTrimCurve(&owner);
            const bool kept = std::any_of(trims.begin(), trims.end(), [&](const auto& t) { return t.first == owner; });
            if (!trim || !owner || kept) continue;
            SmCurve* copy = nullptr;
            if (trim->Copy(context, copy) != SM_SUCCESS || !copy) {
                delete copy;
                return SM_ERR;
            }
            trims.emplace_back(owner, static_cast<SmBSplineCurve*>(copy));
        }
    }

    brep.Notify(SM_NO_PRE_EDIT, &brep, NULL, NULL);
    SmTemporaryChangeValue<SmBoolean> editing(brep.m_bEditingEnabled, TRUE);
    for (ULONG i = 0; i < vertices.GetSize(); ++i) {
        vertices[i]->Notify(SM_NO_PRE_EDIT, vertices[i], &brep, NULL);
        vertices[i]->SetPoint(scaling.point(vertices[i]->GetPoint()));
    }

    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        SmEdge* edge = edges[i];
        SmCurve* curve = edge->GetCurve();
        edge->Notify(SM_NO_PRE_EDIT, edge, &brep, NULL);
        const SM_TYPE type = curve->GetType();
        // Lines and NURBS scale as they are; SmLine::Transform keeps it a line.
        if (uniform || type == SmLine_TYPE || type == SmBSplineCurve_TYPE) {
            SER(curve->Transform(scaling.move, &scaling.factors));
            ++(type == SmBSplineCurve_TYPE ? report.curves_were_nurbs : report.curves_kept);
            continue;
        }
        if (scaling.mirror() && curve->IsKindOf(SmEllipse_TYPE)) {
            const SmAxis2Placement& position = static_cast<SmEllipse*>(curve)->GetPosition();
            SER(curve->Transform(scaling.turn_for(position.GetOriginRef(), position.GetZAxis())));
            ++report.curves_kept;
            continue;
        }
        std::unique_ptr<SmBSplineCurve> scaled(new (context) SmBSplineCurve(*static_cast<SmBSplineCurve*>(curve)));
        SER(scaled->Transform(scaling.move, &scaling.factors));
        SmCurve* replacement = nullptr;
        if (curve->IsKindOf(SmEllipse_TYPE)) {
            bool stretched = false;
            replacement =
                analytic_curve(context, *static_cast<SmEllipse*>(curve), *scaled, scaling, tolerance, stretched);
            if (replacement) ++(stretched ? report.curves_stretched : report.curves_kept);
        }
        if (!replacement) {
            replacement = scaled.release();
            ++report.curves_nurbs;
        }
        // As TurnToNURBS replaces an edge's curve.
        edge->SetCurve(NULL, TRUE);
        replacement->SetOwner(NULL);
        edge->SetCurve(replacement, FALSE);
        replacement->SetOwner(edge);
    }

    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        SmFace* face = faces[i];
        SmSurface* surface = face->GetSurface();
        face->Notify(SM_NO_PRE_EDIT, face, &brep, NULL);
        if (uniform || surface->GetType() == SmBSplineSurface_TYPE) {
            SER(surface->Transform(scaling.move, &scaling.factors));
            ++(surface->GetType() == SmBSplineSurface_TYPE ? report.faces_were_nurbs : report.faces_kept);
            continue;
        }
        if (scaling.mirror() && surface->GetType() == SmPlane_TYPE) {
            const SmAxis2Placement& position = static_cast<SmPlane*>(surface)->GetPosition();
            SER(surface->Transform(scaling.turn_for(position.GetOriginRef(), position.GetZAxis())));
            ++report.faces_kept;
            continue;
        }
        std::unique_ptr<SmBSplineSurface> scaled(
            new (context) SmBSplineSurface(*static_cast<SmBSplineSurface*>(surface)));
        SER(scaled->Transform(scaling.move, &scaling.factors));
        bool extruded = false;
        SmSurface* replacement = analytic_surface(context, *surface, *scaled, scaling, extrusions, tolerance, extruded);
        if (replacement) {
            ++(extruded ? report.faces_extruded : report.faces_kept);
        } else {
            replacement = scaled.release();
            ++report.faces_nurbs;
        }
        const SmExtent2d domain = face->GetUVDomain();
        face->SetSurface(replacement, TRUE);
        face->SetUVDomain(domain);
    }
    if (scaling.reflects()) {
        for (ULONG i = 0; i < edges.GetSize(); ++i) reverse_radial_order(*edges[i]);
        for (ULONG i = 0; i < faces.GetSize(); ++i) swap_sides(*faces[i]);
    }
    for (auto& [owner, trim] : trims) owner->SetUVTrimCurve(trim.release(), 0.0, TRUE);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        if (rectangular[i] == TRUE) faces[i]->SetIsRectangularTrim(TRUE);
    }

    for (ULONG i = 0; i < vertices.GetSize(); ++i) {
        vertices[i]->SetTolerance(SmZoneTol3d(vertex_tolerances[i] * by), FALSE, FALSE);
    }
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        edges[i]->SetTolerance(SmZoneTol3d(edge_tolerances[i] * by), FALSE, FALSE);
    }
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        faces[i]->SetTolerance(SmZoneTol3d(face_tolerances[i] * by), FALSE, FALSE);
    }
    brep.SetTolerance(SmZoneTol3d(brep_tolerance * by), FALSE, FALSE);
    SmTArray<SmRegion*> regions;
    brep.GetRegions(regions);
    if (scaling.reflects() && regions.GetSize() > 1) {
        SmTemporaryChangeValue<SmBoolean> quiet(brep.m_bEditingEnabled, FALSE);
        SER(brep.FindAndSetInfiniteRegion());
    }
    return SM_SUCCESS;
}

// A unit vector square to n.
SmVector3d square_to(const SmVector3d& n) {
    const SmVector3d a = std::fabs(n.x) < 0.9 ? SmVector3d(1.0, 0.0, 0.0) : SmVector3d(0.0, 1.0, 0.0);
    SmVector3d x = a - n * a.Dot(n);
    x.Unitize();
    return x;
}

}  // namespace

extern "C" {

sm_status sm_brep_scale_keeping(sm_brep* brep, sm_vec3 point, sm_vec3 factors, int32_t extrusions,
                                sm_scale_report* report) {
    if (!brep) return SM_SHIM_ERR_NULL_ARGUMENT;
    const auto positive = [](double f) { return std::isfinite(f) && f > SM_EFF_ZERO; };
    if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z) || !positive(factors.x) ||
        !positive(factors.y) || !positive(factors.z)) {
        return SM_SHIM_ERR_BAD_SCALE;
    }
    return guarded([&] {
        Scaling scaling;
        scaling.factors = SmVector3d(factors.x, factors.y, factors.z);
        scaling.move.Translate(SmVector3d(point.x * (1.0 - factors.x), point.y * (1.0 - factors.y),
                                          point.z * (1.0 - factors.z)));
        sm_scale_report counts;
        std::memset(&counts, 0, sizeof counts);
        const sm_status status = scale_keeping(*brep_of(brep), scaling, extrusions != 0, counts);
        if (report) *report = counts;
        return status;
    });
}

sm_status sm_brep_mirror(sm_brep* brep, sm_vec3 point, sm_vec3 normal, sm_scale_report* report) {
    if (!brep) return SM_SHIM_ERR_NULL_ARGUMENT;
    const double length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z) || !std::isfinite(length) ||
        length < SM_EFF_ZERO) {
        return SM_SHIM_ERR_BAD_FRAME;
    }
    return guarded([&] {
        // p goes to R p + 2 (point . n) n, with R = I - 2 n n^T: a scale by -1
        // along z, then the turn taking x and y where R takes them (and z to
        // -R z), then the move.
        const SmVector3d n(normal.x / length, normal.y / length, normal.z / length);
        const SmVector3d x = SmVector3d(1.0, 0.0, 0.0) - n * (2.0 * n.x);
        const SmVector3d y = SmVector3d(0.0, 1.0, 0.0) - n * (2.0 * n.y);
        const double offset = 2.0 * (point.x * n.x + point.y * n.y + point.z * n.z);
        Scaling scaling;
        scaling.factors = SmVector3d(1.0, 1.0, -1.0);
        scaling.move = SmAxis2Placement(SmPoint3d(0.0, 0.0, 0.0) + n * offset, x, y);
        sm_scale_report counts;
        std::memset(&counts, 0, sizeof counts);
        const sm_status status = scale_keeping(*brep_of(brep), scaling, false, counts);
        if (report) *report = counts;
        return status;
    });
}

sm_status sm_brep_mirror_stock(sm_brep* brep, sm_vec3 point, sm_vec3 normal) {
    if (!brep) return SM_SHIM_ERR_NULL_ARGUMENT;
    const double length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z) || !std::isfinite(length) ||
        length < SM_EFF_ZERO) {
        return SM_SHIM_ERR_BAD_FRAME;
    }
    return guarded([&] {
        const SmVector3d n(normal.x / length, normal.y / length, normal.z / length);
        const SmVector3d x = square_to(n);
        const SmAxis2Placement plane(SmPoint3d(point.x, point.y, point.z), x, n * x);
        return static_cast<sm_status>(brep_of(brep)->Mirror(plane));
    });
}

sm_status sm_brep_write_file(const sm_brep* brep, const char* path) {
    if (!brep || !path) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] { return static_cast<sm_status>(SmApiWriteBrepToFile(brep_of(brep), path, FALSE)); });
}

sm_status sm_brep_read_file(const char* path, sm_brep** out) {
    if (!path || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* brep = nullptr;
        const SmApiStatus status = SmApiReadBrepFromFile(path, FALSE, FALSE, brep);
        if (status != SM_SUCCESS) {
            delete brep;
            return static_cast<sm_status>(status);
        }
        *out = handle_of(brep);
        return SM_SUCCESS;
    });
}

}  // extern "C"

// Face-by-face tessellation (sm_tessellation).
//
// SmApiTessellate runs SmBrep::ConvertToPolyBrep, which drives SmTess in three
// phases on a copy of the brep. Phase 1 splits every edge for the whole solid
// and builds each face's surface cache and boundary polygon from those
// splits. Phase 2 then triangulates one face inside its boundary, and the
// polygon output turns that into 3D polygons; both read only that face's
// cache, boundary, surface, trim curves and vertices, and leave the copy's
// topology alone (SmTess.cpp). sm_tessellation drives the phases itself: phase
// 1 for every face, as ConvertToPolyBrep runs it, then phase 2 and the output
// only for the faces asked for, each into a poly brep of its own. So a face
// comes out the same whichever other faces are triangulated, and its key, a
// digest of everything phase 2 and the output read for it, says when an
// earlier result for it can stand in.

struct sm_tessellation {
    SmBrep* source;
    SmCurveTessDriver curves;
    SmSurfaceTessDriver surfaces;
    SmTess tess;
    double max_edge_length;
    sm_tessellation_params params;
    // The copy's faces, in the source's face order.
    SmTArray<SmFace*> faces;
    // Each face's key digest, and whether it couldn't be keyed.
    std::vector<std::array<uint64_t, 2>> keys;
    std::vector<bool> unkeyed;
    bool keyed = false;
    bool meshed = false;

    sm_tessellation(SmBrep& brep, const sm_tessellation_params& p, double edge_length)
        : source(&brep),
          curves(0, p.chord_height_tol, p.curve_angle_tol_deg, edge_length, 0.0001, FALSE),
          surfaces(p.chord_height_tol, p.surface_angle_tol_deg, edge_length, 0.0, 0.001, p.max_aspect_ratio),
          tess(*brep.GetContext(), curves, surfaces),
          max_edge_length(edge_length),
          params(p) {
        tess.SetAdvancingFront(FALSE);
    }
};

namespace {

// ConvertToPolyBrep tessellates with the context's boolean flag set.
struct DoingBoolean {
    SmTemporaryChangeValue<SmBoolean> change;
    explicit DoingBoolean(const SmBrep& brep)
        : change(SM_CONST_CAST(SmContext*, brep.GetContext())->GetDoingBooleanRef(), TRUE) {}
};

// SmApiTessellate's max edge length: a negative value is relative to the
// brep's bounding-box diagonal.
double max_edge_length(SmBrep& brep, double requested) {
    if (requested >= 0.0) return requested;
    SmExtent3d box;
    brep.CalculateBoundingBox(box);
    double size = box.GetSize().Length();
    if (size >= 1.0e15 || size <= 0.0) size = 1.0;
    return -0.025 * size * requested;
}

// Collects a face's key. SMLib writes curves and surfaces in full through
// SmDatabaseIO (WriteToDB); this binary writer keeps the exact bits of every
// value instead of writing a file. Anything written straight to the stream
// is kept too.
class KeyWriter final : public SmDatabaseIO {
public:
    KeyWriter() : SmDatabaseIO(SM_BINARY) { SetOutStreamPtr(&stream_); }

    SmStatus WriteLongs(const ULONG* values, ULONG count) override { return put(values, count * sizeof *values); }
    SmStatus WriteDoubles(const double* values, ULONG count) override { return put(values, count * sizeof *values); }
    SmStatus WriteCharacters(const char* values, ULONG count) override { return put(values, count); }
    SmStatus WriteBooleans(const SmBoolean* values, ULONG count) override {
        return put(values, count * sizeof *values);
    }
    SmStatus WriteLong(long value) override { return put(&value, sizeof value); }
    SmStatus WriteBoolean(SmBoolean value) override { return put(&value, sizeof value); }
    SmStatus WriteDouble(double value) override { return put(&value, sizeof value); }
    SmStatus WriteShort(short value) override { return put(&value, sizeof value); }
    SmStatus WriteChar(char value) override { return put(&value, sizeof value); }

    void point(const SmPoint3d& p) {
        WriteDouble(p.x);
        WriteDouble(p.y);
        WriteDouble(p.z);
    }

    void words(const uint64_t* values, size_t count) { put(values, count * sizeof *values); }

    // Writes a curve or surface the first time it's seen and its index after
    // that. False when SMLib can't write it.
    template <typename Geometry>
    bool geometry(const Geometry* object) {
        if (!object) {
            WriteLong(-1);
            return true;
        }
        const auto found = seen_.find(object);
        if (found != seen_.end()) {
            WriteLong(found->second);
            return true;
        }
        WriteLong(static_cast<long>(seen_.size()));
        seen_.emplace(object, static_cast<long>(seen_.size()));
        WriteLong(static_cast<long>(object->GetType()));
        return object->WriteToDB(*this, SM_CURRENT_DATABASE_VERSION) == SM_SUCCESS;
    }

    std::string finish() {
        bytes_ += stream_.str();
        return std::move(bytes_);
    }

private:
    SmStatus put(const void* data, size_t size) {
        bytes_.append(static_cast<const char*>(data), size);
        return SM_SUCCESS;
    }

    std::ostringstream stream_;
    std::string bytes_;
    std::map<const void*, long> seen_;
};

// MurmurHash3_x64_128 (Austin Appleby, public domain) with seed 0: a face
// key's bytes run to tens of kilobytes, of which callers keep this digest.
void digest(const std::string& bytes, uint64_t out[2]) {
    const auto rotl = [](uint64_t x, int r) { return (x << r) | (x >> (64 - r)); };
    const auto fmix = [](uint64_t k) {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdULL;
        k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53ULL;
        k ^= k >> 33;
        return k;
    };
    const uint64_t c1 = 0x87c37b91114253d5ULL;
    const uint64_t c2 = 0x4cf5ad432745937fULL;
    uint64_t h1 = 0;
    uint64_t h2 = 0;
    const size_t blocks = bytes.size() / 16;
    for (size_t i = 0; i < blocks; ++i) {
        uint64_t k1 = 0;
        uint64_t k2 = 0;
        std::memcpy(&k1, bytes.data() + i * 16, 8);
        std::memcpy(&k2, bytes.data() + i * 16 + 8, 8);
        h1 ^= rotl(k1 * c1, 31) * c2;
        h1 = (rotl(h1, 27) + h2) * 5 + 0x52dce729;
        h2 ^= rotl(k2 * c2, 33) * c1;
        h2 = (rotl(h2, 31) + h1) * 5 + 0x38495ab5;
    }
    // The tail, zero-padded: a zero word leaves the hash as the reference's
    // switch over the tail length does.
    unsigned char tail[16] = {};
    std::memcpy(tail, bytes.data() + blocks * 16, bytes.size() % 16);
    uint64_t k1 = 0;
    uint64_t k2 = 0;
    std::memcpy(&k1, tail, 8);
    std::memcpy(&k2, tail + 8, 8);
    h2 ^= rotl(k2 * c2, 33) * c1;
    h1 ^= rotl(k1 * c1, 31) * c2;
    h1 ^= bytes.size();
    h2 ^= bytes.size();
    h1 += h2;
    h2 += h1;
    h1 = fmix(h1);
    h2 = fmix(h2);
    h1 += h2;
    h2 += h1;
    out[0] = h1;
    out[1] = h2;
}

// The key of one face of the tessellation's copy after phase 1: the
// parameters, the face's surface, domain, flags and tolerances, every
// edgeuse of its loops (UV trim curve, edge curve and interval, start vertex),
// and the boundary polygon phase 1 built, vertex by vertex. Empty when SMLib
// can't write the face's geometry, so it can't be keyed.
std::string face_key(const sm_tessellation& t, SmFace* face) {
    KeyWriter key;
    key.WriteDouble(t.params.chord_height_tol);
    key.WriteDouble(t.params.curve_angle_tol_deg);
    key.WriteDouble(t.params.surface_angle_tol_deg);
    key.WriteDouble(t.max_edge_length);
    key.WriteDouble(t.params.max_aspect_ratio);

    if (!key.geometry(face->GetSurface())) return std::string();
    const SmExtent2d domain = face->GetUVDomain();
    key.WriteDouble(domain.GetUMin());
    key.WriteDouble(domain.GetVMin());
    key.WriteDouble(domain.GetUMax());
    key.WriteDouble(domain.GetVMax());
    key.WriteDouble(face->GetTolerance().val);
    key.WriteBoolean(face->GetRectangularTrim());
    const SmFaceuse* upward = face->GetUpwardFaceuse();
    const SmBrep* brep = face->GetBrep();
    if (!upward || !upward->GetShell() || !brep) return std::string();
    key.WriteBoolean(upward->GetShell()->GetRegion() != brep->GetInfiniteRegion());
    key.WriteDouble(brep->GetTolerance().val);

    std::map<const SmEdgeuse*, long> edgeuses;
    SmTArray<SmLoopuse*> loops;
    face->GetUpwardLoopuses(loops);
    SmTArray<SmEdgeuse*> uses;
    for (ULONG i = 0; i < loops.GetSize(); ++i) {
        key.WriteLong(-2);
        uses.RemoveAll();
        loops[i]->GetEdgeuses(uses);
        for (ULONG j = 0; j < uses.GetSize(); ++j) {
            const SmEdgeuse* use = uses[j];
            edgeuses.emplace(use, static_cast<long>(edgeuses.size()));
            key.WriteLong(static_cast<long>(use->GetOrientation()));
            if (!key.geometry(use->GetUVTrimCurvePointer())) return std::string();
            const SmEdge* edge = use->GetEdge();
            if (!edge) return std::string();
            if (!key.geometry(edge->GetCurve())) return std::string();
            key.WriteDouble(edge->GetInterval().GetMin());
            key.WriteDouble(edge->GetInterval().GetMax());
            key.WriteDouble(edge->GetTolerance().val);
            const SmVertexuse* start = use->GetVertexuse();
            const SmVertex* vertex = start ? start->GetVertex() : nullptr;
            if (!vertex) return std::string();
            key.point(vertex->GetPoint());
            key.WriteDouble(vertex->GetTolerance().val);
        }
    }

    const SmPolyBrep* boundary = t.tess.GetTessBrepOfFace(face);
    if (!boundary) return std::string();
    key.WriteDouble(boundary->GetTolerance().val);
    std::map<const SmPolyVertex*, long> vertices;
    SmTArray<SmPolyFace*> polygons;
    SM_CONST_CAST(SmPolyBrep*, boundary)->GetPolyFaces(polygons);
    SmTArray<SmPolyLoop*> rings;
    SmTArray<SmPolyEdge*> sides;
    for (ULONG i = 0; i < polygons.GetSize(); ++i) {
        key.WriteLong(-3);
        rings.RemoveAll();
        polygons[i]->GetPolyLoops(rings);
        for (ULONG j = 0; j < rings.GetSize(); ++j) {
            key.WriteLong(-4);
            sides.RemoveAll();
            rings[j]->GetPolyEdges(sides);
            for (ULONG k = 0; k < sides.GetSize(); ++k) {
                const SmPolyEdge* side = sides[k];
                const auto found = edgeuses.find(side->GetOriginalEdgeuse());
                key.WriteLong(found == edgeuses.end() ? -1 : found->second);
                const SmPolyVertex* vertex = side->GetStartPolyVertex();
                const auto seen = vertices.find(vertex);
                if (seen != vertices.end()) {
                    key.WriteLong(seen->second);
                    continue;
                }
                key.WriteLong(static_cast<long>(vertices.size()));
                vertices.emplace(vertex, static_cast<long>(vertices.size()));
                if (!vertex) continue;
                key.point(vertex->GetPoint());
                const auto owner = edgeuses.find(vertex->GetOriginalEdgeuse());
                key.WriteLong(owner == edgeuses.end() ? -1 : owner->second);
                key.WriteBoolean(vertex->GetOriginalVertex() != nullptr);
            }
        }
    }
    return key.finish();
}

// Adds the polygons of one face's own poly brep to the mesh, as
// build_triangles does for a whole one.
sm_status add_face_triangles(SmPolyBrep& poly, const SmFace* face, uint32_t face_id, sm_mesh& mesh) {
    SmTArray<SmPolyFace*> polygons;
    poly.GetPolyFaces(polygons);
    MeshBuilder builder{mesh, {}};
    SmTArray<SmPolyLoop*> loops;
    SmTArray<SmPolyEdge*> corners;
    std::vector<uint32_t> ring;
    for (ULONG i = 0; i < polygons.GetSize(); ++i) {
        const SmPolyFace* polygon = polygons[i];
        if (polygon->GetOriginalFace() != face) return SM_SHIM_ERR_UNMAPPED_FACE;
        loops.RemoveAll();
        polygon->GetPolyLoops(loops);
        if (loops.GetSize() != 1) return SM_SHIM_ERR_UNSUPPORTED_POLYGON;
        corners.RemoveAll();
        loops[0]->GetPolyEdges(corners);
        if (corners.GetSize() < 3) return SM_SHIM_ERR_UNSUPPORTED_POLYGON;

        ring.clear();
        for (ULONG j = 0; j < corners.GetSize(); ++j) ring.push_back(builder.vertex(corners[j], face_id));
        if (ring.size() > 3) ++mesh.polygons_split;
        for (size_t j = 1; j + 1 < ring.size(); ++j) {
            mesh.indices.insert(mesh.indices.end(), {ring[0], ring[j], ring[j + 1]});
            mesh.face_ids.push_back(face_id);
        }
    }
    return SM_SUCCESS;
}

// Phase 2 and the polygon output for one face, into a poly brep of its own,
// then the face's triangles into the mesh.
sm_status triangulate_face(sm_tessellation& t, ULONG index, sm_mesh& mesh) {
    SmFace* face = t.faces[index];
    // SmApiTessellate fails when any face does.
    if (t.tess.Phase2CreateFacePolygons(face) != SM_SUCCESS) return SM_ERR;
    SmPolygonOutputCallback output;
    output.SetOutputType(SM_PO_CREATE_POLYBREP);
    SmPolyShell* shell = nullptr;
    const SmStatus status = t.tess.OutputFacePolygons(face, output, &shell);
    t.tess.Phase3DeleteFacePolygons(face);
    SmPolyBrep* poly = t.tess.Get3DPolyBrep();
    t.tess.Set3DPolyBrep(nullptr);
    const std::unique_ptr<SmPolyBrep, PolyBrepDeleter> owned_poly(poly);
    if (status != SM_SUCCESS) return status;
    if (!poly) return SM_SUCCESS;
    // The copy's vertices point at their output vertices, which go with this
    // poly brep; the next face makes its own.
    SmTArray<SmPolyVertex*> outputs;
    poly->GetPolyVertices(outputs);
    for (ULONG i = 0; i < outputs.GetSize(); ++i) {
        SmVertex* vertex = outputs[i]->GetOriginalVertex();
        if (vertex && vertex->GetPolyVertex() == outputs[i]) vertex->SetPolyVertex(nullptr);
    }
    return add_face_triangles(*poly, face, static_cast<uint32_t>(index), mesh);
}

// The mesh of a tessellation whose wanted faces are set up (NULL: every
// face): their triangles, and the brep's faces and edges.
sm_status mesh_faces(sm_tessellation& t, const uint8_t* wanted, sm_mesh** out) {
    SmTArray<SmFace*> faces;
    t.source->GetFaces(faces);
    auto mesh = std::make_unique<sm_mesh>();
    FaceIndex face_index;
    mesh->faces.reserve(faces.GetSize());
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        face_index.emplace(faces[i], static_cast<uint32_t>(i));
        mesh->faces.push_back(describe_face(faces[i]));
    }
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        if (wanted && !wanted[i]) {
            if (i < t.faces.GetSize()) t.tess.Phase3DeleteFacePolygons(t.faces[i]);
            continue;
        }
        const sm_status result = triangulate_face(t, i, *mesh);
        if (result != SM_SUCCESS) return result;
    }
    const sm_status result = build_edges(*t.source, t.params.edge_angle_tol_deg, face_index, *mesh);
    if (result != SM_SUCCESS) return result;
    *out = mesh.release();
    return SM_SUCCESS;
}

// Copies the tessellation's brep, with the copy's faces in t.faces in the
// source's order. The caller owns the copy until phase 1 takes it.
sm_status copy_into(sm_tessellation& t, SmBrep*& copy) {
    SmCopyBrepMap copies;
    copy = new (*t.source->GetContext()) SmBrep(*t.source, &copies);
    SmTArray<SmFace*> originals;
    t.source->GetFaces(originals);
    copy->GetFaces(t.faces);
    bool mapped = t.faces.GetSize() == originals.GetSize();
    for (ULONG i = 0; mapped && i < t.faces.GetSize(); ++i) {
        mapped = copies.GetAt(t.faces[i]) == originals[i];
    }
    if (!mapped) {
        delete copy;
        copy = nullptr;
        t.faces.RemoveAll();
        return SM_SHIM_ERR_UNMAPPED_FACE;
    }
    return SM_SUCCESS;
}

}  // namespace

extern "C" {

sm_status sm_brep_tessellation_begin(sm_brep* brep, const sm_tessellation_params* params, sm_tessellation** out) {
    if (!brep || !params || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* source = brep_of(brep);
        const DoingBoolean doing_boolean(*source);
        auto t = std::make_unique<sm_tessellation>(*source, *params, max_edge_length(*source, params->max_edge_length));
        SmBrep* copy = nullptr;
        const sm_status copied = copy_into(*t, copy);
        if (copied != SM_SUCCESS) return copied;

        // From here the tessellator owns the copy.
        SmBoolean failed = FALSE;
        SmTArray<SmStatus> boundaries;
        const SmStatus status = t->tess.Phase1SetupBrep(copy, failed, NULL, &boundaries);
        if (status != SM_SUCCESS) return status;
        for (ULONG i = 0; i < boundaries.GetSize(); ++i) {
            if (boundaries[i] != SM_SUCCESS) return SM_ERR;
        }
        *out = t.release();
        return SM_SUCCESS;
    });
}

size_t sm_tessellation_face_count(const sm_tessellation* tessellation) {
    return tessellation ? tessellation->faces.GetSize() : 0;
}

sm_status sm_tessellation_face_key(sm_tessellation* tessellation, size_t face, uint64_t key[2], int32_t* keyed) {
    if (!tessellation || !key || !keyed) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (face >= tessellation->faces.GetSize()) return SM_SHIM_ERR_BAD_FACE;
    return guarded([&] {
        sm_tessellation& t = *tessellation;
        if (!t.keyed) {
            if (t.meshed) return SM_SHIM_ERR_SPENT_TESSELLATION;
            const DoingBoolean doing_boolean(*t.source);
            std::vector<std::array<uint64_t, 2>> keys(t.faces.GetSize(), {0, 0});
            t.unkeyed.assign(t.faces.GetSize(), false);
            for (ULONG i = 0; i < t.faces.GetSize(); ++i) {
                const std::string bytes = face_key(t, t.faces[i]);
                if (bytes.empty()) {
                    t.unkeyed[i] = true;
                } else {
                    digest(bytes, keys[i].data());
                }
            }
            t.keys = std::move(keys);
            t.keyed = true;
        }
        key[0] = t.keys[face][0];
        key[1] = t.keys[face][1];
        *keyed = t.unkeyed[face] ? 0 : 1;
        return SM_SUCCESS;
    });
}

sm_status sm_tessellation_mesh(sm_tessellation* tessellation, const uint8_t* wanted, size_t wanted_count,
                               sm_mesh** out) {
    if (!tessellation || !out) return SM_SHIM_ERR_NULL_ARGUMENT;
    sm_tessellation& t = *tessellation;
    if (wanted && wanted_count != t.faces.GetSize()) return SM_SHIM_ERR_BAD_FACE;
    if (t.meshed) return SM_SHIM_ERR_SPENT_TESSELLATION;
    t.meshed = true;
    return guarded([&] {
        const DoingBoolean doing_boolean(*t.source);
        if (t.source->GetNumFaces() != t.faces.GetSize()) return SM_SHIM_ERR_UNMAPPED_FACE;
        return mesh_faces(t, wanted, out);
    });
}

void sm_tessellation_free(sm_tessellation* tessellation) {
    if (!tessellation) return;
    const DoingBoolean doing_boolean(*tessellation->source);
    delete tessellation;
}

}  // extern "C"

// Tessellation of some faces only (sm_brep_tessellate_some).
//
// SmTess::Phase1SetupFaces, a fork patch, sets up some faces exactly as
// phase 1 for the whole brep would, provided no face's boundary needs
// refining: it splits only their edges, with the surface caches of the faces
// using those edges. It checks the faces it sets up; a caller knows the rest
// don't refine by their keys, taken before tessellating. A face's key covers
// everything phase 1 reads for it without refinement, so a face keyed as one
// an earlier sm_brep_tessellate_some set up doesn't refine either, and keeps
// that face's triangles.

namespace {

// Each curve's and surface's digest, written once as SMLib writes it.
class Digests {
public:
    // Puts the object's digest into key, or a mark for none. False when SMLib
    // can't write it.
    template <typename Geometry>
    bool put(KeyWriter& key, const Geometry* object) {
        if (!object) {
            key.WriteLong(-1);
            return true;
        }
        auto found = seen_.find(object);
        if (found == seen_.end()) {
            KeyWriter bytes;
            bytes.WriteLong(static_cast<long>(object->GetType()));
            Seen seen;
            seen.written = object->WriteToDB(bytes, SM_CURRENT_DATABASE_VERSION) == SM_SUCCESS;
            if (seen.written) digest(bytes.finish(), seen.digest.data());
            found = seen_.emplace(object, seen).first;
        }
        if (!found->second.written) return false;
        key.WriteLong(1);
        key.words(found->second.digest.data(), found->second.digest.size());
        return true;
    }

private:
    struct Seen {
        bool written = false;
        std::array<uint64_t, 2> digest{};
    };
    std::unordered_map<const void*, Seen> seen_;
};

// Where a face is in the brep's face list, which orders phase 1's passes.
using FaceOrder = std::unordered_map<const SmFace*, ULONG>;

// What phase 1 reads of a face for an edge it shares: its cache, which
// depends only on its surface, and what its trim curve's creation reads.
bool put_neighbour(KeyWriter& key, Digests& digests, const SmFace& face) {
    if (!digests.put(key, face.GetSurface())) return false;
    const SmExtent2d domain = face.GetUVDomain();
    key.WriteDouble(domain.GetUMin());
    key.WriteDouble(domain.GetVMin());
    key.WriteDouble(domain.GetUMax());
    key.WriteDouble(domain.GetVMax());
    key.WriteDouble(face.GetTolerance().val);
    return true;
}

// A face's key before tessellating (sm_brep_tessellation_keys). Empty when
// SMLib can't write its geometry, or a neighbour's.
std::string input_key(const sm_tessellation_params& params, double edge_length, const SmBrep& brep,
                      const SmFace& face, const FaceOrder& order, Digests& digests) {
    KeyWriter key;
    // Unlike any face_key, which starts with the chord height.
    key.WriteLong(-100);
    key.WriteDouble(params.chord_height_tol);
    key.WriteDouble(params.curve_angle_tol_deg);
    key.WriteDouble(params.surface_angle_tol_deg);
    key.WriteDouble(edge_length);
    key.WriteDouble(params.max_aspect_ratio);
    key.WriteDouble(brep.GetTolerance().val);

    if (!put_neighbour(key, digests, face)) return std::string();
    const SmFaceuse* upward = face.GetUpwardFaceuse();
    if (!upward || !upward->GetShell()) return std::string();
    key.WriteBoolean(upward->GetShell()->GetRegion() != brep.GetInfiniteRegion());
    const auto here = order.find(&face);
    if (here == order.end()) return std::string();

    SmTArray<SmLoopuse*> loops;
    face.GetUpwardLoopuses(loops);
    SmTArray<SmEdgeuse*> uses;
    SmTArray<SmEdgeuse*> sharers;
    for (ULONG i = 0; i < loops.GetSize(); ++i) {
        key.WriteLong(-2);
        uses.RemoveAll();
        loops[i]->GetEdgeuses(uses);
        for (ULONG j = 0; j < uses.GetSize(); ++j) {
            const SmEdgeuse* use = uses[j];
            key.WriteLong(-3);
            key.WriteLong(static_cast<long>(use->GetOrientation()));
            if (!digests.put(key, use->GetUVTrimCurvePointer())) return std::string();
            SmEdge* edge = use->GetEdge();
            if (!edge) return std::string();
            if (!digests.put(key, edge->GetCurve())) return std::string();
            key.WriteDouble(edge->GetInterval().GetMin());
            key.WriteDouble(edge->GetInterval().GetMax());
            key.WriteDouble(edge->GetTolerance().val);
            const SmVertexuse* start = use->GetVertexuse();
            const SmVertex* vertex = start ? start->GetVertex() : nullptr;
            if (!vertex) return std::string();
            key.point(vertex->GetPoint());
            key.WriteDouble(vertex->GetTolerance().val);

            // Every use of the edge, in the order phase 1 visits them, and the
            // faces sharing it, before or after this one.
            sharers.RemoveAll();
            edge->GetEdgeuses(sharers);
            for (ULONG k = 0; k < sharers.GetSize(); ++k) {
                const SmEdgeuse* sharer = sharers[k];
                key.WriteLong(-4);
                key.WriteLong(static_cast<long>(sharer->GetOrientation()));
                if (!digests.put(key, sharer->GetUVTrimCurvePointer())) return std::string();
                const SmFace* other = sharer->GetFace();
                if (other == &face) {
                    key.WriteLong(0);
                    continue;
                }
                const auto there = other ? order.find(other) : order.end();
                if (there == order.end()) return std::string();
                key.WriteLong(there->second < here->second ? 1 : 2);
                if (!put_neighbour(key, digests, *other)) return std::string();
            }
        }
    }
    return key.finish();
}

}  // namespace

// Surface caches kept between sm_brep_tessellate_some calls, for
// SmTess::Phase1SetupFaces: a cache depends only on its surface and the
// parameters, so each is kept by a digest of both.
struct sm_surface_caches final : SmTessSrfCacheStore {
    sm_surface_caches() = default;
    sm_surface_caches(const sm_surface_caches&) = delete;
    sm_surface_caches& operator=(const sm_surface_caches&) = delete;

    ~sm_surface_caches() override {
        for (auto& [key, caches] : kept) {
            for (const Kept& cache : caches) delete cache.cache;
        }
    }

    // Starts a call with these parameters.
    void begin(const sm_tessellation_params& p, double length) {
        ++call;
        params = p;
        edge_length = length;
        digests = Digests();
        missed[call % 3].clear();
    }

    // Ends a call, dropping the caches neither it nor the call before kept.
    void end() {
        for (auto entry = kept.begin(); entry != kept.end();) {
            auto& caches = entry->second;
            for (auto cache = caches.begin(); cache != caches.end();) {
                if (cache->call + 1 >= call) {
                    ++cache;
                    continue;
                }
                delete cache->cache;
                cache = caches.erase(cache);
            }
            entry = caches.empty() ? kept.erase(entry) : std::next(entry);
        }
    }

    size_t count() const {
        size_t n = 0;
        for (const auto& [key, caches] : kept) n += caches.size();
        return n;
    }

    SmTessSrfCache* Take(SmFace* face) override {
        std::array<uint64_t, 2> digest;
        if (!key(*face, digest)) return nullptr;
        const auto found = kept.find(digest);
        if (found == kept.end()) {
            missed[call % 3].insert(digest);
            return nullptr;
        }
        SmTessSrfCache* cache = found->second.back().cache;
        found->second.pop_back();
        if (found->second.empty()) kept.erase(found);
        return cache;
    }

    void Keep(SmFace* face, SmTessSrfCache* cache) override {
        std::array<uint64_t, 2> digest;
        if (!key(*face, digest)) {
            delete cache;
            return;
        }
        kept[digest].push_back({cache, call});
    }

    // A face being tessellated built its cache: keeping a copy pays when one
    // of the last two calls built a cache for the same surface too, as for a
    // face whose boundary moves in a drag while its surface stays. Surfaces
    // that change every call (the faces of the feature dragged) aren't kept.
    SmBoolean WantsCopy(SmFace* face) override {
        std::array<uint64_t, 2> digest;
        if (!key(*face, digest)) return FALSE;
        return missed[(call + 1) % 3].count(digest) > 0 || missed[(call + 2) % 3].count(digest) > 0;
    }

private:
    struct Kept {
        SmTessSrfCache* cache;
        uint64_t call;
    };

    // What fixes a face's cache: its surface, and the parameters the
    // tessellation gives its surface driver. False when SMLib can't write the
    // surface.
    bool key(const SmFace& face, std::array<uint64_t, 2>& out) {
        KeyWriter bytes;
        bytes.WriteDouble(params.chord_height_tol);
        bytes.WriteDouble(params.surface_angle_tol_deg);
        bytes.WriteDouble(edge_length);
        bytes.WriteDouble(params.max_aspect_ratio);
        if (!digests.put(bytes, face.GetSurface())) return false;
        digest(bytes.finish(), out.data());
        return true;
    }

    std::map<std::array<uint64_t, 2>, std::vector<Kept>> kept;
    // The surfaces each of the last three calls found no cache for, by call % 3.
    std::array<std::set<std::array<uint64_t, 2>>, 3> missed;
    uint64_t call = 0;
    sm_tessellation_params params{};
    double edge_length = 0.0;
    Digests digests;
};

extern "C" {

sm_status sm_brep_tessellation_keys(sm_brep* brep, const sm_tessellation_params* params, uint64_t* keys,
                                    int32_t* keyed, size_t face_count) {
    if (!brep || !params || (face_count > 0 && (!keys || !keyed))) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* source = brep_of(brep);
        SmTArray<SmFace*> faces;
        source->GetFaces(faces);
        if (face_count != faces.GetSize()) return SM_SHIM_ERR_BAD_FACE;
        const double edge_length = max_edge_length(*source, params->max_edge_length);
        FaceOrder order;
        for (ULONG i = 0; i < faces.GetSize(); ++i) order.emplace(faces[i], i);
        Digests digests;
        for (ULONG i = 0; i < faces.GetSize(); ++i) {
            const std::string bytes = input_key(*params, edge_length, *source, *faces[i], order, digests);
            keyed[i] = bytes.empty() ? 0 : 1;
            keys[2 * i] = 0;
            keys[2 * i + 1] = 0;
            if (!bytes.empty()) digest(bytes, keys + 2 * i);
        }
        return SM_SUCCESS;
    });
}

sm_status sm_brep_tessellate_some(sm_brep* brep, const sm_tessellation_params* params, const uint8_t* wanted,
                                  size_t wanted_count, sm_surface_caches* caches, sm_mesh** out) {
    if (!brep || !params || !out || (wanted_count > 0 && !wanted)) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* source = brep_of(brep);
        if (wanted_count != source->GetNumFaces()) return SM_SHIM_ERR_BAD_FACE;
        const DoingBoolean doing_boolean(*source);
        auto t = std::make_unique<sm_tessellation>(*source, *params, max_edge_length(*source, params->max_edge_length));
        if (std::any_of(wanted, wanted + wanted_count, [](uint8_t w) { return w != 0; })) {
            if (caches) caches->begin(*params, t->max_edge_length);
            // Ends the call however it goes, once phase 1 has given caches back.
            const std::unique_ptr<sm_surface_caches, void (*)(sm_surface_caches*)> ending(
                caches, [](sm_surface_caches* c) { c->end(); });
            SmBrep* copy = nullptr;
            const sm_status copied = copy_into(*t, copy);
            if (copied != SM_SUCCESS) return copied;
            SmTArray<SmFace*> some;
            for (ULONG i = 0; i < t->faces.GetSize(); ++i) {
                if (wanted[i]) some.Add(t->faces[i]);
            }
            // From here the tessellator owns the copy.
            SmBoolean same = FALSE;
            const SmStatus status = t->tess.Phase1SetupFaces(copy, some, same, caches);
            if (status != SM_SUCCESS) return status;
            if (!same) {
                *out = nullptr;
                return SM_SUCCESS;
            }
        }
        return mesh_faces(*t, wanted, out);
    });
}

sm_surface_caches* sm_surface_caches_new(void) { return new (std::nothrow) sm_surface_caches(); }

size_t sm_surface_caches_count(const sm_surface_caches* caches) { return caches ? caches->count() : 0; }

void sm_surface_caches_free(sm_surface_caches* caches) { delete caches; }

}  // extern "C"

// Imprinting. SmApiBooleanWithOptions can imprint, but it deletes the tool
// or keeps it without saying which of its faces each new edge lies on, and
// it can't classify faces and imprint without the rest of the boolean. So
// this drives SmMerge as it does, keeping the tool until the edges' mates on
// it have been read.

namespace {

// The labels of the faces an edge of the tool lies on, at most two, each
// once, in the order the edge lists its faces.
void tool_face_labels(const SmEdge* edge, sm_label out[2]) {
    std::memset(out, 0, 2 * sizeof(sm_label));
    SmTArray<SmFace*> faces;
    edge->GetFaces(faces);
    size_t written = 0;
    for (ULONG i = 0; i < faces.GetSize() && written < 2; ++i) {
        const sm_label label = read_label(faces[i]);
        if (!label.present) continue;
        if (written == 1 && std::memcmp(&out[0], &label, sizeof label) == 0) continue;
        out[written++] = label;
    }
}

// Takes the classifying merge's marks off a brep: SM_AI_BOOLEAN_SAVE and
// SM_AI_BOOLEAN_DELETE on faces, SM_AI_BOOLEAN_INTERSECT on edges and
// vertices.
void unmark_merge(SmBrep& brep) {
    const auto unmark = [](SmAObject* object, ULONG id) {
        while (SmAttribute* mark = object->FindAttribute(id)) object->RemoveAttribute(mark);
    };
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        unmark(faces[i], SM_AI_BOOLEAN_SAVE);
        unmark(faces[i], SM_AI_BOOLEAN_DELETE);
    }
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) unmark(edges[i], SM_AI_BOOLEAN_INTERSECT);
    SmTArray<SmVertex*> vertices;
    brep.GetVertices(vertices);
    for (ULONG i = 0; i < vertices.GetSize(); ++i) unmark(vertices[i], SM_AI_BOOLEAN_INTERSECT);
}

// Marks every edge of a brep before an imprint, so the edges it adds are
// the ones without the mark; split edges pass it to both pieces.
const ULONG OLD_EDGE_ATTRIBUTE_ID = 0x534d4f45;  // "SMOE"

void mark_edges(SmBrep& brep) {
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        edges[i]->AddAttribute(new (*edges[i]->GetContext()) SmGenericAttribute(
            OLD_EDGE_ATTRIBUTE_ID, SM_AB_COPY, SmTArray<long>(), SmTArray<double>(), SmTArray<char>()));
    }
}

void unmark_edges(SmBrep& brep) {
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        while (SmAttribute* mark = edges[i]->FindAttribute(OLD_EDGE_ATTRIBUTE_ID)) edges[i]->RemoveAttribute(mark);
    }
}

}  // namespace

extern "C" {

sm_status sm_brep_imprint(sm_brep* brep, sm_brep* tool, int32_t classify, sm_brep** out, sm_imprint_edge** edges,
                          size_t* edge_count, uint8_t** inside, size_t* face_count) {
    if (!brep || !tool || !out || !edges || !edge_count || !inside || !face_count) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (brep == tool) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    return guarded([&] {
        SmBrep* target = brep_of(brep);
        SmBrep* other = brep_of(tool);
        std::vector<sm_imprint_edge> along;
        std::vector<uint8_t> in;
        mark_edges(*target);
        {
            // Imprinting stops before the boolean deletes any face or merges
            // faces on one surface; classifying goes on to sort the faces as
            // an intersection would keep them, marking them instead.
            SmMergeOptions options(FALSE, classify ? FALSE : TRUE, FALSE, classify ? TRUE : FALSE, TRUE);
            SmMerge merge(*SmApiGetOrCreateContext(), target, other, 0, 0, &options);
            const SmBooleanOperationType operation = classify ? SM_BO_INTERSECTION : SM_BO_IMPRINT;
            SmBrep* result = nullptr;
            const SmStatus status = target->IsManifoldSolid() && other->IsManifoldSolid()
                                        ? merge.ManifoldBoolean(operation, result)
                                        : merge.NonManifoldBoolean(operation, result);
            if (status != SM_SUCCESS) return static_cast<sm_status>(status);
            if (result != target) return static_cast<sm_status>(SM_ERR);

            SmTArray<SmEdge*> common;
            SmTArray<SmEdge*> mates;
            merge.GetTopologyIntersector().GetCommonEdges(common, mates);
            SmTArray<SmEdge*> all;
            target->GetEdges(all);
            std::unordered_map<const SmEdge*, size_t> index;
            for (ULONG i = 0; i < all.GetSize(); ++i) index.emplace(all[i], i);
            SmTArray<SmEdge*> tool_edges;
            other->GetEdges(tool_edges);
            std::unordered_set<const SmEdge*> on_tool;
            for (ULONG i = 0; i < tool_edges.GetSize(); ++i) on_tool.insert(tool_edges[i]);
            for (ULONG i = 0; i < common.GetSize(); ++i) {
                const auto found = index.find(common[i]);
                if (found == index.end() || !on_tool.count(mates[i])) continue;
                sm_imprint_edge edge;
                std::memset(&edge, 0, sizeof edge);
                edge.edge = found->second;
                edge.is_new = common[i]->FindAttribute(OLD_EDGE_ATTRIBUTE_ID) == NULL ? 1 : 0;
                tool_face_labels(mates[i], edge.tool_faces);
                along.push_back(edge);
            }
            std::sort(along.begin(), along.end(),
                      [](const sm_imprint_edge& a, const sm_imprint_edge& b) { return a.edge < b.edge; });

            if (classify) {
                SmTArray<SmFace*> faces;
                target->GetFaces(faces);
                for (ULONG i = 0; i < faces.GetSize(); ++i) {
                    in.push_back(faces[i]->FindAttribute(SM_AI_BOOLEAN_SAVE) != NULL ? 1 : 0);
                }
                unmark_merge(*target);
            }
            unmark_edges(*target);
        }
        std::unique_ptr<sm_imprint_edge[]> edge_array(new sm_imprint_edge[along.size()]);
        std::copy(along.begin(), along.end(), edge_array.get());
        std::unique_ptr<uint8_t[]> inside_array(in.empty() ? nullptr : new uint8_t[in.size()]);
        std::copy(in.begin(), in.end(), inside_array.get());
        delete other;
        *out = brep;
        *edges = edge_array.release();
        *edge_count = along.size();
        *inside = inside_array.release();
        *face_count = in.size();
        return SM_SUCCESS;
    });
}

void sm_imprint_free(sm_imprint_edge* edges, uint8_t* inside) {
    delete[] edges;
    delete[] inside;
}

sm_status sm_brep_merge_faces_keeping(sm_brep* brep, const uint32_t* keep, size_t count) {
    if (!brep || (count > 0 && !keep)) return SM_SHIM_ERR_NULL_ARGUMENT;
    return guarded([&] {
        SmBrep* target = brep_of(brep);
        SmTArray<SmEdge*> kept;
        if (edges_at(*target, keep, count, kept) != SM_SUCCESS) return SM_SHIM_ERR_BAD_EDGE;
        return static_cast<sm_status>(target->RemoveTopologicalEdgesAndVertices(&kept));
    });
}

}  // extern "C"

// Solids swept from faces, beyond sm_brep_sweep_face: revolves and tapered
// sweeps.

namespace {

// A new brep holding a copy of `face` alone (SmBrep::CopyFaces, without
// merging): *copy is that face, and source_edge gives each of its edges the
// index in the source's edge list of the edge it was copied from (-1 when
// none is found). The copy's edges are new curves, so they're traced back by
// where they lie: an edge's midpoint lies within the edge it was copied from.
sm_status face_alone(const SmBrep& brep, const SmFace* face, std::unique_ptr<SmBrep>& alone, SmFace*& copy,
                     std::unordered_map<const SmEdge*, int32_t>& source_edge) {
    const SmContext* context = SmApiGetOrCreateContext();
    alone.reset(new (*context) SmBrep());
    SmTArray<SmFace*> faces;
    faces.Add(SM_CONST_CAST(SmFace*, face));
    SmTArray<SmFace*> copies;
    const SmStatus status = SM_CONST_CAST(SmBrep*, &brep)->CopyFaces(faces, alone.get(), &copies, FALSE);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    if (copies.GetSize() != 1 || !copies[0]) return SM_ERR;
    copy = copies[0];

    SmTArray<SmEdge*> all;
    brep.GetEdges(all);
    std::unordered_map<const SmEdge*, int32_t> index;
    for (ULONG i = 0; i < all.GetSize(); ++i) index.emplace(all[i], static_cast<int32_t>(i));
    SmTArray<SmEdge*> originals;
    face->GetEdges(originals);
    SmTArray<SmEdge*> edges;
    copy->GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        const SmPoint3d mid = midpoint_of(edges[i]);
        int32_t best = -1;
        double best_distance = SM_BIG_DOUBLE;
        for (ULONG j = 0; j < originals.GetSize(); ++j) {
            SmPoint3d closest;
            double parameter = 0.0;
            double distance = SM_BIG_DOUBLE;
            const SmExtent1d interval = originals[j]->GetInterval();
            if (SmApiCurveClosestPoint(originals[j]->GetCurve(), mid, closest, parameter, distance) != SM_SUCCESS ||
                !interval.ContainsValue(parameter, 1e-6 * interval.GetLength()) || distance >= best_distance) {
                continue;
            }
            best_distance = distance;
            best = index[originals[j]];
        }
        source_edge[edges[i]] = best_distance <= 1e3 * SM_ZONE_TOL_3D ? best : -1;
    }
    return SM_SUCCESS;
}

// Takes every label off the brep's faces and edges.
void remove_labels(SmBrep& brep) {
    const sm_label none{};
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) write_label(*faces[i], none);
    SmTArray<SmEdge*> edges;
    brep.GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize(); ++i) write_label(*edges[i], none);
}

// Whether the line through `base` along `axis` runs through the inside of the
// face (SmFace::CurveOnClassify, as SmPrimitiveCreation's rotational sweep
// checks it). Along its boundary or through a vertex doesn't count.
sm_status axis_crosses(const SmFace* face, const SmPoint3d& base, const SmVector3d& axis, bool& crosses) {
    crosses = false;
    SmExtent3d box;
    SmStatus status = face->CalculateBoundingBox(box);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    const double reach = 2.0 * (box.GetSize().Length() + (box.GetMid() - base).Length());
    const SmContext* context = SmApiGetOrCreateContext();
    SmBSplineCurve* raw = nullptr;
    status = SmBSplineCurve::CreateLineSegment(*context, 3, base - axis * reach, base + axis * reach, raw);
    if (status != SM_SUCCESS || !raw) return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
    const std::unique_ptr<SmBSplineCurve> line(raw);
    SmCurveClassification classes(line.get(), line->GetNaturalInterval(), NULL, face->GetTolerance());
    status = face->CurveOnClassify(TRUE, classes);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    for (ULONG i = 0; i < classes.GetSize(); ++i) {
        if (classes[i].m_vMid.GetPointClass() == SM_PC_FACE) crosses = true;
    }
    return SM_SUCCESS;
}

// The rotational sweep makes its planes, cones, cylinders and spheres
// analytic (SmBSplineSurface::CopyAndAddAnalytics) but leaves tori as
// surfaces of revolution, since CopyAndAddAnalytics doesn't look for them.
// Gives each side that is one its torus (SmTorus::IsNurbSurfaceTorus), when
// the torus has the surface's parameters (the same domain, and the same
// point at every corner of a grid over it) and those points lie on the torus
// its centre, axis and radii describe, which is what callers read: the
// recogniser can evaluate as the surface does while describing another
// torus (before patch c7639c2). Otherwise the side keeps its surface of
// revolution.
void make_tori(SmBrep& brep) {
    const SmContext* context = brep.GetContext();
    SmTemporaryChangeValue<SmBoolean> editing(brep.m_bEditingEnabled, TRUE);
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const SmSurface* surface = faces[i]->GetSurface();
        if (!surface || (surface->GetType() != SmSurfOfRevolution_TYPE && surface->GetType() != SmBSplineSurface_TYPE)) {
            continue;
        }
        SmTorus* raw = nullptr;
        if (!SmTorus::IsNurbSurfaceTorus(*context, static_cast<const SmBSplineSurface*>(surface), raw) || !raw) {
            delete raw;
            continue;
        }
        std::unique_ptr<SmSurface> torus(raw);
        const SmExtent2d domain = surface->GetNaturalUVDomain();
        if (!(torus->GetNaturalUVDomain() == domain)) continue;
        SmExtent3d box;
        if (faces[i]->CalculateBoundingBox(box) != SM_SUCCESS) continue;
        const double tolerance = 1e-9 * (1.0 + box.GetSize().Length());
        const SmAxis2Placement& position = raw->GetPosition();
        const SmVector3d axis = position.GetZAxis();
        bool agrees = true;
        for (int a = 0; a <= 4 && agrees; ++a) {
            for (int b = 0; b <= 4 && agrees; ++b) {
                const SmPoint2d uv = domain.Evaluate(a / 4.0, b / 4.0);
                SmPoint3d on_surface, on_torus;
                agrees = surface->EvaluatePoint(uv, on_surface) == SM_SUCCESS &&
                         torus->EvaluatePoint(uv, on_torus) == SM_SUCCESS &&
                         (on_surface - on_torus).Length() <= tolerance;
                const SmVector3d from_centre = on_surface - position.GetOrigin();
                const double height = from_centre.Dot(axis);
                const double out = (from_centre - axis * height).Length() - raw->GetMajorRadius();
                agrees = agrees && std::fabs(std::hypot(out, height) - raw->GetMinorRadius()) <= tolerance;
            }
        }
        if (agrees) replace_surface(faces[i], std::move(torus));
    }
}

// What a face of a swept solid came from, kept while booleans put the solid
// together (a face with holes): a label with this feature id, never a
// caller's, whose role is the sm_swept_role and whose index the source
// edge's.
const uint32_t SWEPT_ROLE = UINT32_MAX;

void mark_role(SmFace& face, int32_t role, int32_t source_edge) {
    sm_label label{};
    label.present = 1;
    label.feature = SWEPT_ROLE;
    label.index = static_cast<uint32_t>(source_edge);
    std::snprintf(label.role, sizeof label.role, "%d", static_cast<int>(role));
    write_label(face, label);
}

void mark_side(SmFace& side, int32_t source_edge) { mark_role(side, SM_SWEPT_SIDE, source_edge); }

// The roles mark_role left on the solid's faces, in its face order, and the
// solid without them or any other label.
void take_roles(SmBrep& solid, std::vector<sm_swept_face>& roles) {
    SmTArray<SmFace*> faces;
    solid.GetFaces(faces);
    roles.assign(faces.GetSize(), sm_swept_face{SM_SWEPT_SIDE, -1});
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const sm_label label = read_label(faces[i]);
        if (!label.present || label.feature != SWEPT_ROLE) continue;
        roles[i].role = static_cast<int32_t>(std::atoi(label.role));
        roles[i].source_edge = roles[i].role == SM_SWEPT_SIDE ? static_cast<int32_t>(label.index) : -1;
    }
    remove_labels(solid);
}

// The frame a revolve works in: the axis along z through the origin, the
// face in the xz plane on the side of +x. The sweep and the analytic
// surfaces are made there whatever plane the face is in, and the solid is
// moved back after.
SmAxis2Placement revolve_frame(const SmPoint3d& base, const SmVector3d& axis, const SmPoint3d& inside) {
    SmVector3d across = inside - base;
    across -= axis * across.Dot(axis);
    across.Unitize(FALSE);
    return SmAxis2Placement(base, across, axis * across);
}

// The face copied alone and moved from `frame` to the revolve's own.
sm_status profile_in(const SmBrep& brep, const SmFace* face, const SmAxis2Placement& frame,
                     std::unique_ptr<SmBrep>& alone, SmFace*& profile,
                     std::unordered_map<const SmEdge*, int32_t>& source_edge) {
    const sm_status status = face_alone(brep, face, alone, profile, source_edge);
    if (status != SM_SUCCESS) return status;
    SmAxis2Placement to_frame;
    frame.Invert(to_frame);
    return static_cast<sm_status>(SmApiTransform(alone.get(), to_frame));
}

// A turn short of a whole one: the face swept into the solid, itself the
// start and its copy the end (SmTopologySweep's maps), each side its edge's.
sm_status turn_part(const SmBrep& brep, const SmFace* face, const SmAxis2Placement& frame, double angle,
                    std::unique_ptr<SmBrep>& solid) {
    SmFace* profile = nullptr;
    std::unordered_map<const SmEdge*, int32_t> source_edge;
    sm_status status = profile_in(brep, face, frame, solid, profile, source_edge);
    if (status != SM_SUCCESS) return status;
    SmRotationalSweepGeometry geometry(SmPoint3d(0.0, 0.0, 0.0), SmVector3d(0.0, 0.0, 1.0), angle,
                                       solid->GetTolerance());
    SmTopologySweep sweep(geometry);
    sweep.SetRepetitions(1);
    SmTArray<SmFace*> swept;
    swept.Add(profile);
    {
        SmTemporaryChangeValue<SmBoolean> editing(solid->m_bEditingEnabled, TRUE);
        const SmStatus done = sweep.DoSweep(solid.get(), solid->GetInfiniteRegion(), &swept);
        if (done != SM_SUCCESS) return static_cast<sm_status>(done);
    }
    if (!is_manifold_solid(*solid)) return SM_SHIM_ERR_BAD_SWEEP;
    make_tori(*solid);
    remove_labels(*solid);
    mark_role(*profile, SM_SWEPT_START, -1);
    if (SmFace* copy = sweep.GetToSame(profile)) mark_role(*copy, SM_SWEPT_END, -1);
    for (const auto& [edge, index] : source_edge) {
        if (SmFace* side = sweep.GetToHigher(SM_CONST_CAST(SmEdge*, edge))) mark_side(*side, index);
    }
    return SM_SUCCESS;
}

// A whole turn of one loop of the face: its edges swept round into closed
// sides, stitched shut (SmTopologySweep, then SmBrep::StitchAndOrient). The
// face itself would meet its own copy, so it isn't swept. Its loops are
// copied in one order every time, so `loop` picks the same one in each copy.
sm_status turn_loop(const SmBrep& brep, const SmFace* face, const SmAxis2Placement& frame, ULONG loop,
                    std::unique_ptr<SmBrep>& solid, bool& outer) {
    SmFace* profile = nullptr;
    std::unordered_map<const SmEdge*, int32_t> source_edge;
    sm_status status = profile_in(brep, face, frame, solid, profile, source_edge);
    if (status != SM_SUCCESS) return status;
    SmTArray<SmLoop*> loops;
    profile->GetLoops(loops);
    SmTArray<SmLoop*> outers;
    profile->GetOuterLoops(outers);
    if (loop >= loops.GetSize()) return SM_ERR;
    ULONG unused;
    outer = outers.FindElement(loops[loop], unused) == TRUE;
    SmTArray<SmEdge*> edges;
    loops[loop]->GetEdges(edges);

    SmRotationalSweepGeometry geometry(SmPoint3d(0.0, 0.0, 0.0), SmVector3d(0.0, 0.0, 1.0), 360.0,
                                       solid->GetTolerance());
    SmTopologySweep sweep(geometry);
    sweep.SetRepetitions(1);
    sweep.SetDoStitching(FALSE);
    {
        SmTemporaryChangeValue<SmBoolean> editing(solid->m_bEditingEnabled, TRUE);
        SmStatus done = sweep.DoSweep(solid.get(), solid->GetInfiniteRegion(), nullptr, &edges);
        if (done != SM_SUCCESS) return static_cast<sm_status>(done);
        SmTArray<SmFace*> start;
        start.Add(profile);
        SmTopologyIntersector intersector(*solid->GetContext());
        done = intersector.DeleteFaces(start);
        if (done != SM_SUCCESS) return static_cast<sm_status>(done);
        done = solid->StitchAndOrient();
        if (done != SM_SUCCESS) return static_cast<sm_status>(done);
    }
    if (!is_manifold_solid(*solid)) return SM_SHIM_ERR_BAD_SWEEP;
    make_tori(*solid);
    remove_labels(*solid);
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        if (SmFace* side = sweep.GetToHigher(edges[i])) mark_side(*side, source_edge[edges[i]]);
    }
    return SM_SUCCESS;
}

// Cuts each of `holes` out of the solid (SmApiBooleanDifference).
sm_status cut_holes(std::unique_ptr<SmBrep>& solid, std::vector<std::unique_ptr<SmBrep>>& holes) {
    for (auto& hole : holes) {
        SmBrep* a = solid.release();
        SmBrep* result = nullptr;
        const SmApiStatus status = SmApiBooleanDifference(a, hole.get(), result);
        if (status != SM_SUCCESS || result != a) {
            delete a;
            return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
        }
        hole.release();  // the boolean deleted it
        solid.reset(result);
    }
    return SM_SUCCESS;
}

// A whole turn: each loop turned on its own, and the holes' rings cut out of
// the outer loop's solid, leaving closed voids.
sm_status turn_whole(const SmBrep& brep, const SmFace* face, const SmAxis2Placement& frame,
                     std::unique_ptr<SmBrep>& solid) {
    SmTArray<SmLoop*> loops;
    face->GetLoops(loops);
    std::vector<std::unique_ptr<SmBrep>> holes;
    for (ULONG i = 0; i < loops.GetSize(); ++i) {
        std::unique_ptr<SmBrep> turned;
        bool outer = false;
        const sm_status status = turn_loop(brep, face, frame, i, turned, outer);
        if (status != SM_SUCCESS) return status;
        if (outer) {
            if (solid) return SM_ERR;
            solid = std::move(turned);
        } else {
            holes.push_back(std::move(turned));
        }
    }
    if (!solid) return SM_ERR;
    return cut_holes(solid, holes);
}

sm_status revolve_face(const SmBrep& brep, const SmFace* face, SmPoint3d base, SmVector3d axis, double angle,
                       SmBrep*& out, std::vector<sm_swept_face>& roles) {
    double normal[3];
    SmPoint3d inside;
    if (!outward_plane(face, normal) || SmApiFaceInternalPoint(face, inside) != SM_SUCCESS) {
        return SM_SHIM_ERR_NOT_PLANAR;
    }
    const SmVector3d n(normal[0], normal[1], normal[2]);
    if (std::fabs(axis.Dot(n)) > 1e-9 || std::fabs((base - inside).Dot(n)) > brep.GetTolerance()) {
        return SM_SHIM_ERR_BAD_AXIS;
    }
    bool crosses = false;
    sm_status status = axis_crosses(face, base, axis, crosses);
    if (status != SM_SUCCESS) return status;
    if (crosses) return SM_SHIM_ERR_AXIS_CROSSES_FACE;
    // A turn the other way is a turn this way about the axis reversed.
    if (angle < 0.0) {
        axis *= -1.0;
        angle = -angle;
    }
    const SmAxis2Placement frame = revolve_frame(base, axis, inside);
    SmBoolean whole = FALSE;
    {
        SmRotationalSweepGeometry geometry(base, axis, angle, brep.GetTolerance());
        if (geometry.IsSweepClosed(1, whole) != SM_SUCCESS) return SM_SHIM_ERR_BAD_SWEEP;
    }
    std::unique_ptr<SmBrep> solid;
    status = whole ? turn_whole(brep, face, frame, solid) : turn_part(brep, face, frame, angle, solid);
    if (status != SM_SUCCESS) return status;

    take_roles(*solid, roles);
    const SmApiStatus back = SmApiTransform(solid.get(), frame);
    if (back != SM_SUCCESS) return static_cast<sm_status>(back);
    out = solid.release();
    return SM_SUCCESS;
}

// The unit normal of `surface` at the point on it nearest `near`.
bool normal_near(const SmSurface* surface, const SmPoint3d& near, SmVector3d& normal) {
    SmPoint3d closest;
    SmTArray<SmPoint2d> uvs;
    double distance = 0.0;
    return SmApiSurfaceClosestPoint(surface, near, closest, uvs, distance) == SM_SUCCESS && uvs.GetSize() > 0 &&
           SmApiEvaluateSurfaceNormal(surface, SmVector2d(uvs[0].x, uvs[0].y), normal) == SM_SUCCESS;
}

// Puts `surface` under the face in place of its own, which it must match in
// shape though not in parameters: reversed if need be to face the way the old
// one did, as SmBrep::LocalOperation does, and with the face's UV trims made
// again on it (SmFace::CreateUVTrimCurves), over its whole domain.
sm_status rebuild_on(SmFace* face, std::unique_ptr<SmSurface> surface) {
    SmPoint3d inside;
    SmVector3d before, after;
    if (SmApiFaceInternalPoint(face, inside) != SM_SUCCESS || !normal_near(face->GetSurface(), inside, before) ||
        !normal_near(surface.get(), inside, after)) {
        return SM_ERR;
    }
    SmStatus status = SM_SUCCESS;
    if (before.Dot(after) < 0.0) {
        status = surface->Reverse(SM_SP_U);
        if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    }
    status = face->RemoveUVTrimCurves();
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    const SmExtent2d domain = surface->GetNaturalUVDomain();
    face->SetSurface(surface.release(), TRUE);
    face->SetUVDomain(domain);
    double mean_gap, max_gap, mean_vertex_gap, max_vertex_gap, mean_uv_gap, max_uv_gap;
    status = face->CreateUVTrimCurves(FALSE, nullptr, nullptr, mean_gap, max_gap, mean_vertex_gap, max_vertex_gap,
                                      mean_uv_gap, max_uv_gap);
    return static_cast<sm_status>(status);
}

// Puts each planar NURBS face of the brep on a plane of its own: SMLib's
// recogniser (SmPlane::IsNurbSurfacePlane) only takes patches whose
// parameters run square across them, which a tapered wall's trapezoid
// doesn't. The plane is made from the face's own edges
// (SmSurface::CreatePlaneWith3DCurves).
sm_status make_planes(SmBrep& brep) {
    const SmContext* context = brep.GetContext();
    SmTemporaryChangeValue<SmBoolean> editing(brep.m_bEditingEnabled, TRUE);
    SmTArray<SmFace*> faces;
    brep.GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        SmFace* face = faces[i];
        const SmSurface* surface = face->GetSurface();
        SmPoint3d on_plane;
        SmVector3d plane_normal;
        if (!surface || surface->GetType() != SmBSplineSurface_TYPE || !plane_of(surface, on_plane, plane_normal)) {
            continue;
        }
        SmTArray<ULONG> counts;
        SmTArray<SmCurve*> curves;
        SmTArray<SmOrientType> orients;
        SmTArray<SmPoint3d> points;
        SmStatus status = face->CreateCurvesFromFace(*context, SM_OT_SAME, counts, &curves, nullptr, orients, points);
        std::vector<std::unique_ptr<SmCurve>> owned;
        for (ULONG j = 0; j < curves.GetSize(); ++j) owned.emplace_back(curves[j]);
        if (status != SM_SUCCESS) return static_cast<sm_status>(status);
        SmSurface* raw = nullptr;
        SmTArray<ULONG> loop_counts;
        SmTArray<SmCurve*> ordered;
        SmTArray<SmOrientType> ordered_orients;
        status = SmSurface::CreatePlaneWith3DCurves(*context, curves, brep.GetTolerance(), raw, loop_counts, ordered,
                                                    ordered_orients);
        std::unique_ptr<SmSurface> plane(raw);
        if (status != SM_SUCCESS || !plane) return status != SM_SUCCESS ? static_cast<sm_status>(status) : SM_ERR;
        const sm_status rebuilt = rebuild_on(face, std::move(plane));
        if (rebuilt != SM_SUCCESS) return rebuilt;
    }
    return SM_SUCCESS;
}

// Copies of the curves of each edge loop of `face` (SmFace::CreateCurvesFromFace),
// loop by loop, and which loop is the outer one.
sm_status loop_curves(const SmFace* face, std::vector<std::vector<std::unique_ptr<SmCurve>>>& loops, size_t& outer) {
    const SmContext* context = SmApiGetOrCreateContext();
    SmTArray<ULONG> counts;
    SmTArray<SmCurve*> curves;
    SmTArray<SmOrientType> orients;
    SmTArray<SmPoint3d> points;
    const SmStatus status = face->CreateCurvesFromFace(*context, SM_OT_SAME, counts, &curves, nullptr, orients, points);
    std::vector<std::unique_ptr<SmCurve>> owned;
    for (ULONG i = 0; i < curves.GetSize(); ++i) owned.emplace_back(curves[i]);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    size_t next = 0;
    for (ULONG i = 0; i < counts.GetSize(); ++i) {
        if (next + counts[i] > owned.size()) return SM_ERR;
        loops.emplace_back();
        for (ULONG j = 0; j < counts[i]; ++j) loops.back().push_back(std::move(owned[next++]));
    }
    // CreateCurvesFromFace takes the edge loops in the faceuse's order.
    const SmFaceuse* use = face->GetUpwardFaceuse()->GetOrientation() == SM_OT_SAME ? face->GetUpwardFaceuse()
                                                                                     : face->GetUpwardFaceuse()->GetMate();
    SmTArray<SmLoopuse*> loopuses;
    use->GetLoopuses(loopuses);
    SmTArray<SmLoop*> outers;
    face->GetOuterLoops(outers);
    outer = loops.size();
    size_t k = 0;
    for (ULONG i = 0; i < loopuses.GetSize(); ++i) {
        if (loopuses[i]->IsVertexLoopuse()) continue;
        ULONG unused;
        if (outers.FindElement(loopuses[i]->GetLoop(), unused)) outer = k;
        ++k;
    }
    return outer < loops.size() && k == loops.size() ? SM_SUCCESS : SM_ERR;
}

// One loop's curves swept `height` along `direction` with sharp corners and
// walls leaning `draft` degrees from it: towards the loop's inside when
// positive, out when negative (SmPrimitiveCreation::CreateTaperExtrude,
// which runs the loop anticlockwise about `direction` and offsets it to the
// left for a positive draft). Both ends are capped.
sm_status taper_loop(std::vector<std::unique_ptr<SmCurve>>& loop, const SmVector3d& direction, double height,
                     double draft, std::unique_ptr<SmBrep>& solid) {
    const SmContext* context = SmApiGetOrCreateContext();
    solid.reset(new (*context) SmBrep());
    SmTArray<SmCurve*> curves;
    for (auto& curve : loop) curves.Add(curve.get());
    SmVector3d normal = direction;
    SmPrimitiveCreation creator(solid->GetInfiniteRegion());
    const SmStatus status = creator.CreateTaperExtrude(curves, height, draft, 3, normal, SM_OC_LINEAR_EXTENSION);
    // It takes the curves on success; on failure they are still ours, but
    // for any it deleted as duplicates.
    for (auto& curve : loop) {
        ULONG unused;
        if (status == SM_SUCCESS || !curves.FindElement(curve.get(), unused)) curve.release();
    }
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    return is_manifold_solid(*solid) ? SM_SUCCESS : SM_SHIM_ERR_BAD_SWEEP;
}

// Marks the faces of a loop swept by taper_loop from `face`: the start, on
// the face's plane, the end, `height` along `direction`, and each side by the
// edge of `face` its edge on the start lies along (source_edge). A side
// standing on a corner of the loop, with no edge on the start, has none.
// False when there's no start or no end: the sides met before the end, and
// CreateTaperExtrude stopped there.
bool mark_taper(SmBrep& solid, const SmBrep& brep, const SmFace* face, const SmVector3d& direction, double height,
                bool caps, SmFace*& end) {
    SmTArray<SmFace*> faces;
    solid.GetFaces(faces);
    SmPoint3d on_face;
    SmApiFaceInternalPoint(face, on_face);
    const double tolerance = std::max(1e3 * SM_ZONE_TOL_3D, 1e-9 * height);
    SmFace* start = nullptr;
    end = nullptr;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        double normal[3];
        SmPoint3d inside;
        if (!outward_plane(faces[i], normal) || SmApiFaceInternalPoint(faces[i], inside) != SM_SUCCESS ||
            std::fabs(std::fabs(SmVector3d(normal[0], normal[1], normal[2]).Dot(direction)) - 1.0) > 1e-9) {
            continue;
        }
        const double along = (inside - on_face).Dot(direction);
        if (std::fabs(along) <= tolerance) {
            start = faces[i];
            if (caps) mark_role(*faces[i], SM_SWEPT_START, -1);
        } else if (std::fabs(along - height) <= tolerance) {
            end = faces[i];
            if (caps) mark_role(*faces[i], SM_SWEPT_END, -1);
        }
    }
    SmTArray<SmEdge*> all_edges;
    brep.GetEdges(all_edges);
    std::unordered_map<const SmEdge*, uint32_t> edge_index;
    for (ULONG i = 0; i < all_edges.GetSize(); ++i) edge_index.emplace(all_edges[i], i);
    SmTArray<SmEdge*> face_edges;
    face->GetEdges(face_edges);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        if (faces[i] == start || faces[i] == end) continue;
        mark_side(*faces[i], start ? source_edge(faces[i], start, face_edges, edge_index) : -1);
    }
    return start && end;
}

// Whether every edge of `inner` lies inside `outer`, a face in its plane,
// and off its boundary (SmFace::CurveOnClassify).
sm_status lies_inside(const SmFace* inner, const SmFace* outer, bool& inside) {
    inside = true;
    SmTArray<SmEdge*> edges;
    inner->GetEdges(edges);
    for (ULONG i = 0; i < edges.GetSize() && inside; ++i) {
        SmCurveClassification classes(edges[i]->GetCurve(), edges[i]->GetInterval(), NULL, outer->GetTolerance());
        const SmStatus status = outer->CurveOnClassify(TRUE, classes);
        if (status != SM_SUCCESS) return static_cast<sm_status>(status);
        for (ULONG j = 0; j < classes.GetSize(); ++j) {
            if (classes[j].m_vMid.GetPointClass() != SM_PC_FACE) inside = false;
        }
    }
    return SM_SUCCESS;
}

// The face swept `distance` along its outward normal (into the solid when
// negative) with its walls leaning `draft` degrees: each loop swept by
// taper_loop, the outer one leaning in and the holes out so the face shrinks
// as a whole (grows for a negative draft), and the holes' solids cut from the
// outer one's.
sm_status taper_face(const SmBrep& brep, const SmFace* face, double distance, double draft, SmBrep*& out,
                     std::vector<sm_swept_face>& roles) {
    double normal[3];
    if (!outward_plane(face, normal)) return SM_SHIM_ERR_NOT_PLANAR;
    const SmVector3d direction = SmVector3d(normal[0], normal[1], normal[2]) * (distance < 0.0 ? -1.0 : 1.0);
    const double height = std::fabs(distance);
    std::vector<std::vector<std::unique_ptr<SmCurve>>> loops;
    size_t outer = 0;
    sm_status status = loop_curves(face, loops, outer);
    if (status != SM_SUCCESS) return status;
    std::unique_ptr<SmBrep> solid;
    std::vector<std::unique_ptr<SmBrep>> holes;
    SmFace* end = nullptr;
    std::vector<SmFace*> hole_ends;
    for (size_t i = 0; i < loops.size(); ++i) {
        std::unique_ptr<SmBrep> swept;
        status = taper_loop(loops[i], direction, height, i == outer ? draft : -draft, swept);
        if (status != SM_SUCCESS) return status;
        SmTArray<SmFace*> faces;
        swept->GetFaces(faces);
        make_sides_analytic(*swept, faces);
        status = make_planes(*swept);
        if (status != SM_SUCCESS) return status;
        status = static_cast<sm_status>(swept->RemoveTopologicalEdgesAndVertices());
        if (status != SM_SUCCESS) return status;
        SmFace* swept_end = nullptr;
        if (!mark_taper(*swept, brep, face, direction, height, i == outer, swept_end)) {
            return SM_SHIM_ERR_DRAFT_CLOSES;
        }
        (i == outer ? end : hole_ends.emplace_back()) = swept_end;
        (i == outer ? solid : holes.emplace_back()) = std::move(swept);
    }
    // A hole leaning out mustn't reach the outer loop's sides by the end.
    for (const SmFace* hole_end : hole_ends) {
        bool inside = false;
        status = lies_inside(hole_end, end, inside);
        if (status != SM_SUCCESS) return status;
        if (!inside) return SM_SHIM_ERR_DRAFT_CLOSES;
    }
    status = cut_holes(solid, holes);
    if (status != SM_SUCCESS) return status;
    if (!is_manifold_solid(*solid)) return SM_SHIM_ERR_BAD_SWEEP;
    take_roles(*solid, roles);
    out = solid.release();
    return SM_SUCCESS;
}

// One curve of a pipe's path, in the frame's own (x, y): where it starts and
// ends and its unit tangents there, along the way it runs.
struct PipeSegment {
    sm_plane_curve curve;
    bool round = false;  // an arc or a circle
    double start[2] = {0.0, 0.0};
    double end[2] = {0.0, 0.0};
    double start_tangent[2] = {0.0, 0.0};
    double end_tangent[2] = {0.0, 0.0};
    double sweep = 0.0;  // an arc's turn in degrees, signed; a circle's 360
};

PipeSegment pipe_segment(const sm_plane_curve& curve) {
    PipeSegment segment;
    segment.curve = curve;
    if (curve.kind == SM_PLANE_LINE) {
        const double dx = curve.b[0] - curve.a[0], dy = curve.b[1] - curve.a[1];
        const double length = std::hypot(dx, dy);
        segment.start[0] = curve.a[0];
        segment.start[1] = curve.a[1];
        segment.end[0] = curve.b[0];
        segment.end[1] = curve.b[1];
        segment.start_tangent[0] = segment.end_tangent[0] = dx / length;
        segment.start_tangent[1] = segment.end_tangent[1] = dy / length;
        return segment;
    }
    segment.round = true;
    const double from = curve.kind == SM_PLANE_CIRCLE ? 0.0 : curve.start_deg;
    segment.sweep = curve.kind == SM_PLANE_CIRCLE ? 360.0 : curve.sweep_deg;
    const double to = from + segment.sweep;
    const double turn = segment.sweep < 0.0 ? -1.0 : 1.0;
    const double a = from * SM_PI / 180.0, b = to * SM_PI / 180.0;
    segment.start[0] = curve.a[0] + curve.radius * std::cos(a);
    segment.start[1] = curve.a[1] + curve.radius * std::sin(a);
    segment.end[0] = curve.a[0] + curve.radius * std::cos(b);
    segment.end[1] = curve.a[1] + curve.radius * std::sin(b);
    segment.start_tangent[0] = -turn * std::sin(a);
    segment.start_tangent[1] = turn * std::cos(a);
    segment.end_tangent[0] = -turn * std::sin(b);
    segment.end_tangent[1] = turn * std::cos(b);
    return segment;
}

// The sharpest turn a pipe's path may take where two curves meet: 179
// degrees, in radians.
const double MAX_PIPE_TURN = 179.0 * SM_PI / 180.0;

// Where two curves of a pipe's path meet: the angle the path turns through
// there, in radians, and the plane halfway between them through the corner
// (its unit normal, in the frame's (x, y)), where they're mitred unless they
// meet tangentially.
struct PipeJoint {
    double turn = 0.0;
    double corner[2] = {0.0, 0.0};
    double mitre[2] = {0.0, 0.0};
    bool mitred() const { return turn > 1e-9; }
};

PipeJoint pipe_joint(const PipeSegment& before, const PipeSegment& after) {
    PipeJoint joint;
    const double* t1 = before.end_tangent;
    const double* t2 = after.start_tangent;
    joint.turn = std::atan2(std::fabs(t1[0] * t2[1] - t1[1] * t2[0]), t1[0] * t2[0] + t1[1] * t2[1]);
    joint.corner[0] = after.start[0];
    joint.corner[1] = after.start[1];
    const double mx = t1[0] + t2[0], my = t1[1] + t2[1];
    const double length = std::hypot(mx, my);
    joint.mitre[0] = length > 0.0 ? mx / length : t1[0];
    joint.mitre[1] = length > 0.0 ? my / length : t1[1];
    return joint;
}

// The labels a pipe's pieces carry while they're joined: SWEPT_ROLE's
// feature id, the sm_pipe_role as the role and the path curve as the index.
void mark_pipe(SmFace& face, int32_t role, int32_t segment) {
    sm_label label{};
    label.present = 1;
    label.feature = SWEPT_ROLE;
    label.index = static_cast<uint32_t>(segment);
    std::snprintf(label.role, sizeof label.role, "pipe %d", static_cast<int>(role));
    write_label(face, label);
}

// One curve's piece of a solid pipe of `radius`: a disc square to the curve
// at its start, swept along a line (sweep_face) or turned about an arc's
// centre (revolve_face), `before` and `after` longer at each end (in model
// units along the curve), its wall marked with `wall` and the curve's index,
// and its caps at the path's ends when `first` or `last`.
sm_status pipe_piece(const PipeSegment& segment, int32_t index, const SmAxis2Placement& frame, double radius,
                     int32_t wall, double before, double after, bool first, bool last,
                     std::unique_ptr<SmBrep>& piece) {
    const SmVector3d x = frame.GetXAxis(), y = frame.GetYAxis(), normal = frame.GetZAxis();
    const auto world = [&](double u, double v) { return frame.GetOrigin() + x * u + y * v; };
    // Where the piece starts, extended back along the curve.
    SmPoint3d origin;
    SmVector3d tangent;
    double angle = 0.0;
    SmPoint3d centre;
    if (!segment.round) {
        tangent = x * segment.start_tangent[0] + y * segment.start_tangent[1];
        origin = world(segment.start[0], segment.start[1]) - tangent * before;
    } else {
        const double r = segment.curve.radius;
        const double turn = segment.sweep < 0.0 ? -1.0 : 1.0;
        const double from = (segment.curve.kind == SM_PLANE_CIRCLE ? 0.0 : segment.curve.start_deg) -
                            turn * before / r * 180.0 / SM_PI;
        angle = segment.sweep + turn * (before + after) / r * 180.0 / SM_PI;
        if (std::fabs(angle) > 360.0) return SM_SHIM_ERR_BAD_PATH;
        const double a = from * SM_PI / 180.0;
        centre = world(segment.curve.a[0], segment.curve.a[1]);
        origin = world(segment.curve.a[0] + r * std::cos(a), segment.curve.a[1] + r * std::sin(a));
        tangent = (x * -std::sin(a) + y * std::cos(a)) * turn;
    }
    // The disc, facing along the curve.
    sm_plane_curve circle{};
    circle.kind = SM_PLANE_CIRCLE;
    circle.radius = radius;
    SmBrep* raw = nullptr;
    sm_status status = planar_faces(&circle, 1, raw);
    if (status != SM_SUCCESS) return status;
    std::unique_ptr<SmBrep> disc(raw);
    status = static_cast<sm_status>(SmApiTransform(disc.get(), SmAxis2Placement(origin, normal, tangent * normal)));
    if (status != SM_SUCCESS) return status;
    SmTArray<SmFace*> disc_faces;
    disc->GetFaces(disc_faces);
    if (disc_faces.GetSize() != 1) return SM_ERR;

    SmBrep* swept = nullptr;
    std::vector<sm_swept_face> roles;
    if (!segment.round) {
        const double length = std::hypot(segment.end[0] - segment.start[0], segment.end[1] - segment.start[1]);
        status = sweep_face(*disc, disc_faces[0], length + before + after, swept, roles);
    } else {
        status = revolve_face(*disc, disc_faces[0], centre, normal, angle, swept, roles);
    }
    if (status != SM_SUCCESS) return status;
    piece.reset(swept);
    SmTArray<SmFace*> faces;
    piece->GetFaces(faces);
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        if (roles[i].role == SM_SWEPT_SIDE) {
            mark_pipe(*faces[i], wall, index);
        } else if (roles[i].role == SM_SWEPT_START && first) {
            mark_pipe(*faces[i], SM_PIPE_START, -1);
        } else if (roles[i].role == SM_SWEPT_END && last) {
            mark_pipe(*faces[i], SM_PIPE_END, -1);
        }
    }
    return SM_SUCCESS;
}

// Cuts away the part of the piece beyond the plane through `corner` with
// unit normal `away` (SmApiBooleanDifference with a box that size covers).
sm_status cut_beyond(std::unique_ptr<SmBrep>& piece, const SmPoint3d& corner, const SmVector3d& away, double size) {
    SmVector3d across = away * SmVector3d(0.0, 0.0, 1.0);
    if (across.Length() < 1e-6) across = away * SmVector3d(1.0, 0.0, 0.0);
    across.Unitize(FALSE);
    SmBrep* box = nullptr;
    SmApiStatus status = SmApiCreateBox(SmVector3d(-size, -size, 0.0), 2.0 * size, 2.0 * size, size, box);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    std::unique_ptr<SmBrep> beyond(box);
    status = SmApiTransform(beyond.get(), SmAxis2Placement(corner, across, away * across));
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    std::vector<std::unique_ptr<SmBrep>> cut;
    cut.push_back(std::move(beyond));
    return cut_holes(piece, cut);
}

// A solid pipe of `radius` along the path: a piece for each curve, mitred to
// its neighbours where they meet at an angle (each made longer at that end,
// by enough for its cap to clear the plane halfway between them, and cut
// there) and met face to face where they meet tangentially, all joined
// (SmApiBooleanUnion). An open path's ends are made `margin` longer, and
// capped when `caps`. Walls are marked `wall`.
sm_status solid_pipe(const std::vector<PipeSegment>& segments, const std::vector<PipeJoint>& joints, bool closed,
                     const SmAxis2Placement& frame, double radius, int32_t wall, double margin, bool caps,
                     double size, std::unique_ptr<SmBrep>& pipe) {
    const size_t count = segments.size();
    const auto extension = [&](const PipeJoint& joint) {
        return joint.mitred() ? radius * std::tan(joint.turn / 2.0) + radius : 0.0;
    };
    const SmVector3d x = frame.GetXAxis(), y = frame.GetYAxis();
    for (size_t i = 0; i < count; ++i) {
        const PipeJoint* behind = i > 0 ? &joints[i - 1] : (closed && count > 1 ? &joints.back() : nullptr);
        const PipeJoint* ahead = i < joints.size() && (i + 1 < count || (closed && count > 1)) ? &joints[i] : nullptr;
        const double before = behind ? extension(*behind) : (closed ? 0.0 : margin);
        const double after = ahead ? extension(*ahead) : (closed ? 0.0 : margin);
        std::unique_ptr<SmBrep> piece;
        sm_status status = pipe_piece(segments[i], static_cast<int32_t>(i), frame, radius, wall, before, after,
                                      caps && i == 0 && !closed, caps && i + 1 == count && !closed, piece);
        if (status != SM_SUCCESS) return status;
        for (const auto& [joint, keep_before] : {std::make_pair(behind, false), std::make_pair(ahead, true)}) {
            if (!joint || !joint->mitred()) continue;
            const SmPoint3d corner = frame.GetOrigin() + x * joint->corner[0] + y * joint->corner[1];
            const SmVector3d mitre = x * joint->mitre[0] + y * joint->mitre[1];
            status = cut_beyond(piece, corner, keep_before ? mitre : mitre * -1.0, size);
            if (status != SM_SUCCESS) return status;
        }
        if (!pipe) {
            pipe = std::move(piece);
            continue;
        }
        SmBrep* a = pipe.release();
        SmBrep* joined = nullptr;
        const SmApiStatus union_status = SmApiBooleanUnion(a, piece.get(), joined);
        if (union_status != SM_SUCCESS || joined != a) {
            delete a;
            return union_status != SM_SUCCESS ? static_cast<sm_status>(union_status) : SM_ERR;
        }
        piece.release();  // the union deleted it
        pipe.reset(joined);
    }
    return is_manifold_solid(*pipe) ? SM_SUCCESS : SM_SHIM_ERR_BAD_SWEEP;
}

// A pipe along the path: a solid pipe of `radius`, less a solid pipe of
// `inner` (when above 0) running past the path's ends, so no ring meets
// another face to face. Only an open path's pipe can be hollow
// (SM_SHIM_ERR_BAD_PATH): a closed one's bore would be a sealed cavity, some
// of whose walls SMLib's tessellation winds the wrong way. A bend meeting its neighbour at an angle can't
// be mitred with a plane (a cylinder's and a torus's sections there differ),
// so SM_SHIM_ERR_ARC_CORNER. The faces' marks say what each is.
sm_status create_pipe(const sm_plane_curve* path, size_t count, const SmAxis2Placement& frame, double radius,
                      double inner, SmBrep*& out, std::vector<sm_pipe_face>& roles) {
    std::vector<PipeSegment> segments;
    for (size_t i = 0; i < count; ++i) segments.push_back(pipe_segment(path[i]));
    const PipeSegment& first = segments.front();
    const PipeSegment& last = segments.back();
    const bool closed = std::hypot(last.end[0] - first.start[0], last.end[1] - first.start[1]) <= 1e-6;
    if (closed && inner > 0.0) return SM_SHIM_ERR_BAD_PATH;
    // joints[i] is where curve i meets the next; the last only when closed.
    std::vector<PipeJoint> joints;
    for (size_t i = 0; i + 1 < count; ++i) joints.push_back(pipe_joint(segments[i], segments[i + 1]));
    if (closed && count > 1) joints.push_back(pipe_joint(last, first));
    for (size_t i = 0; i < joints.size(); ++i) {
        const bool round = segments[i].round || segments[(i + 1) % count].round;
        if (joints[i].turn > MAX_PIPE_TURN) return SM_SHIM_ERR_BAD_PATH;
        if (round && joints[i].mitred()) return SM_SHIM_ERR_ARC_CORNER;
    }
    // Big enough for a box to cover any piece beyond a mitre.
    double size = 8.0 * radius;
    for (const PipeSegment& segment : segments) {
        for (const double* p : {segment.start, segment.end}) size += std::hypot(p[0], p[1]);
        if (segment.round) size += 2.0 * segment.curve.radius;
    }
    for (const PipeJoint& joint : joints) size += 2.0 * radius * (std::tan(joint.turn / 2.0) + 1.0);

    std::unique_ptr<SmBrep> pipe;
    sm_status status = solid_pipe(segments, joints, closed, frame, radius, SM_PIPE_OUTER, 0.0, true, size, pipe);
    if (status != SM_SUCCESS) return status;
    if (inner > 0.0) {
        std::vector<std::unique_ptr<SmBrep>> bore(1);
        status = solid_pipe(segments, joints, closed, frame, inner, SM_PIPE_INNER, radius, false, size, bore[0]);
        if (status == SM_SUCCESS) status = cut_holes(pipe, bore);
        if (status != SM_SUCCESS) return status;
        if (!is_manifold_solid(*pipe)) return SM_SHIM_ERR_BAD_SWEEP;
    }
    SmTArray<SmFace*> faces;
    pipe->GetFaces(faces);
    roles.assign(faces.GetSize(), sm_pipe_face{SM_PIPE_OUTER, -1});
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const sm_label label = read_label(faces[i]);
        int role = 0;
        if (!label.present || label.feature != SWEPT_ROLE || std::sscanf(label.role, "pipe %d", &role) != 1) continue;
        roles[i].role = role;
        roles[i].segment = role == SM_PIPE_OUTER || role == SM_PIPE_INNER ? static_cast<int32_t>(label.index) : -1;
    }
    remove_labels(*pipe);
    out = pipe.release();
    return SM_SUCCESS;
}

// How far `point` is from the surface `face` lies on.
bool off_surface(const SmFace* face, const SmPoint3d& point, double& distance) {
    SmPoint3d closest;
    SmTArray<SmPoint2d> uvs;
    return SmApiSurfaceClosestPoint(face->GetSurface(), point, closest, uvs, distance) == SM_SUCCESS;
}

// The brep hollowed `thickness` in from its faces, the faces at `open` left
// open (SmApiShellBrepFull: SmOffsetExecutive::ShellBrep insets the solid,
// those faces left where they are, and cuts the inset from it). Each face is
// marked with its index first, in place of its label, so each face of the
// result says which face it came from: one on that face's own surface stays
// (SM_SHELL_OUTER, or SM_SHELL_RIM for what's left of an opened face, round
// its opening), one `thickness` in is its offset (SM_SHELL_INNER), and one
// with no mark is new (SM_SHELL_NEW). Faces that stay get their labels back,
// as do the edges only they use; inner faces and their edges have none.
sm_status shell_brep(const SmBrep& brep, const std::vector<bool>& opened, double thickness, SmBrep*& out,
                     std::vector<sm_shell_face>& roles) {
    if (!is_manifold_solid(brep)) return SM_SHIM_ERR_BAD_SHELL;
    const SmContext* context = SmApiGetOrCreateContext();
    SmTArray<SmFace*> sources;
    brep.GetFaces(sources);
    std::unordered_map<const SmTopology*, int32_t> source_index;
    for (ULONG i = 0; i < sources.GetSize(); ++i) source_index.emplace(sources[i], static_cast<int32_t>(i));

    SmCopyBrepMap to_original;
    std::unique_ptr<SmBrep> work(new (*context) SmBrep(brep, &to_original, nullptr, TRUE, FALSE));
    std::vector<sm_label> labels(sources.GetSize());
    SmTArray<SmFace*> faces;
    work->GetFaces(faces);
    SmTArray<const SmFace*> open_faces;
    for (ULONG i = 0; i < faces.GetSize(); ++i) {
        const auto found = source_index.find(to_original.GetAt(faces[i]));
        if (found == source_index.end()) return SM_ERR;
        labels[found->second] = read_label(faces[i]);
        sm_label mark{};
        mark.present = 1;
        mark.feature = SWEPT_ROLE;
        mark.index = static_cast<uint32_t>(found->second);
        std::memcpy(mark.role, "shell", 6);
        write_label(*faces[i], mark);
        if (opened[found->second]) open_faces.Add(faces[i]);
    }
    SmBrep* raw = nullptr;
    const SmApiStatus status = SmApiShellBrepFull(work.get(), -thickness, TRUE, TRUE, FALSE, open_faces, raw);
    std::unique_ptr<SmBrep> shell(raw);
    if (status != SM_SUCCESS) return static_cast<sm_status>(status);
    if (!shell || !is_manifold_solid(*shell)) return SM_SHIM_ERR_BAD_SHELL;
    // The offset can come back with faces whose own integrals disagree with
    // their edges (a box with its top's rim rounded, opened at the top).
    bool agrees = false;
    const sm_status checked = volume_agrees(*shell, agrees);
    if (checked != SM_SUCCESS) return checked;
    if (!agrees) return SM_SHIM_ERR_BAD_SHELL;

    SmTArray<SmFace*> result;
    shell->GetFaces(result);
    roles.assign(result.GetSize(), sm_shell_face{SM_SHELL_NEW, -1});
    const double tolerance = std::max(1e3 * SM_ZONE_TOL_3D, 1e-6 * thickness);
    std::vector<bool> stays(result.GetSize(), false);
    for (ULONG i = 0; i < result.GetSize(); ++i) {
        const sm_label mark = read_label(result[i]);
        if (!mark.present || mark.feature != SWEPT_ROLE || std::strcmp(mark.role, "shell") != 0 ||
            mark.index >= sources.GetSize()) {
            write_label(*result[i], sm_label{});
            continue;
        }
        const int32_t source = static_cast<int32_t>(mark.index);
        SmPoint3d inside;
        double distance = 0.0;
        if (SmApiFaceInternalPoint(result[i], inside) != SM_SUCCESS ||
            !off_surface(sources[static_cast<ULONG>(source)], inside, distance)) {
            return SM_ERR;
        }
        roles[i].source_face = source;
        if (distance <= tolerance) {
            roles[i].role = opened[source] ? SM_SHELL_RIM : SM_SHELL_OUTER;
            stays[i] = true;
            write_label(*result[i], labels[source]);
        } else {
            roles[i].role = std::fabs(distance - thickness) <= tolerance ? SM_SHELL_INNER : SM_SHELL_NEW;
            write_label(*result[i], sm_label{});
        }
    }
    // An edge keeps its label only where faces that stay are all that use it.
    SmTArray<SmEdge*> edges;
    shell->GetEdges(edges);
    std::unordered_map<const SmFace*, bool> face_stays;
    for (ULONG i = 0; i < result.GetSize(); ++i) face_stays[result[i]] = stays[i];
    for (ULONG i = 0; i < edges.GetSize(); ++i) {
        SmTArray<SmFace*> users;
        edges[i]->GetFaces(users);
        bool kept = users.GetSize() > 0;
        for (ULONG j = 0; j < users.GetSize(); ++j) kept = kept && face_stays[users[j]];
        if (!kept) write_label(*edges[i], sm_label{});
    }
    out = shell.release();
    return SM_SUCCESS;
}

}  // namespace

extern "C" {

sm_status sm_brep_revolve_face(const sm_brep* brep, size_t face, sm_vec3 axis_point, sm_vec3 axis_direction,
                               double angle_deg, sm_brep** out, sm_swept_face** faces, size_t* face_count) {
    if (!brep || !out || !faces || !face_count) return SM_SHIM_ERR_NULL_ARGUMENT;
    const double length = std::sqrt(axis_direction.x * axis_direction.x + axis_direction.y * axis_direction.y +
                                    axis_direction.z * axis_direction.z);
    if (!std::isfinite(axis_point.x) || !std::isfinite(axis_point.y) || !std::isfinite(axis_point.z) ||
        !std::isfinite(length) || length == 0.0) {
        return SM_SHIM_ERR_BAD_AXIS;
    }
    if (!std::isfinite(angle_deg) || angle_deg == 0.0 || std::fabs(angle_deg) > 360.0) return SM_SHIM_ERR_BAD_SWEEP;
    return guarded([&] {
        const SmBrep& source = *brep_of(brep);
        const SmFace* target = face_at(source, face);
        if (!target) return SM_SHIM_ERR_BAD_FACE;
        const SmVector3d axis(axis_direction.x / length, axis_direction.y / length, axis_direction.z / length);
        SmBrep* solid = nullptr;
        std::vector<sm_swept_face> roles;
        const sm_status status = revolve_face(source, target, SmPoint3d(axis_point.x, axis_point.y, axis_point.z),
                                              axis, angle_deg, solid, roles);
        if (status != SM_SUCCESS) return status;
        std::unique_ptr<SmBrep> owned(solid);
        std::unique_ptr<sm_swept_face[]> array(new sm_swept_face[roles.size()]);
        std::copy(roles.begin(), roles.end(), array.get());
        *out = handle_of(owned.release());
        *faces = array.release();
        *face_count = roles.size();
        return SM_SUCCESS;
    });
}

sm_status sm_brep_sweep_face_tapered(const sm_brep* brep, size_t face, double distance, double draft_deg,
                                     sm_brep** out, sm_swept_face** faces, size_t* face_count) {
    if (!brep || !out || !faces || !face_count) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!std::isfinite(distance) || distance == 0.0) return SM_SHIM_ERR_BAD_SWEEP;
    if (!std::isfinite(draft_deg) || std::fabs(draft_deg) >= 90.0) return SM_SHIM_ERR_BAD_DRAFT;
    return guarded([&] {
        const SmBrep& source = *brep_of(brep);
        SmFace* target = face_at(source, face);
        if (!target) return SM_SHIM_ERR_BAD_FACE;
        SmBrep* solid = nullptr;
        std::vector<sm_swept_face> roles;
        const sm_status status = draft_deg == 0.0 ? sweep_face(source, target, distance, solid, roles)
                                                  : taper_face(source, target, distance, draft_deg, solid, roles);
        if (status != SM_SUCCESS) return status;
        std::unique_ptr<SmBrep> owned(solid);
        std::unique_ptr<sm_swept_face[]> array(new sm_swept_face[roles.size()]);
        std::copy(roles.begin(), roles.end(), array.get());
        *out = handle_of(owned.release());
        *faces = array.release();
        *face_count = roles.size();
        return SM_SUCCESS;
    });
}

sm_status sm_brep_create_pipe(const sm_plane_curve* path, size_t count, sm_vec3 origin, sm_vec3 x_axis,
                              sm_vec3 y_axis, double radius, double inner_radius, sm_brep** out, sm_pipe_face** faces,
                              size_t* face_count) {
    if (!path || !out || !faces || !face_count) return SM_SHIM_ERR_NULL_ARGUMENT;
    const double dot = x_axis.x * y_axis.x + x_axis.y * y_axis.y + x_axis.z * y_axis.z;
    if (!std::isfinite(origin.x) || !std::isfinite(origin.y) || !std::isfinite(origin.z) || !unit_vector(x_axis) ||
        !unit_vector(y_axis) || std::fabs(dot) > 1e-9) {
        return SM_SHIM_ERR_BAD_FRAME;
    }
    if (count == 0 || !std::isfinite(radius) || radius <= 0.0 || !std::isfinite(inner_radius) || inner_radius < 0.0 ||
        inner_radius >= radius) {
        return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    }
    return guarded([&] {
        const SmAxis2Placement frame(SmPoint3d(origin.x, origin.y, origin.z), SmVector3d(x_axis.x, x_axis.y, x_axis.z),
                                     SmVector3d(y_axis.x, y_axis.y, y_axis.z));
        SmBrep* pipe = nullptr;
        std::vector<sm_pipe_face> roles;
        const sm_status status = create_pipe(path, count, frame, radius, inner_radius, pipe, roles);
        if (status != SM_SUCCESS) return status;
        std::unique_ptr<SmBrep> owned(pipe);
        std::unique_ptr<sm_pipe_face[]> array(new sm_pipe_face[roles.size()]);
        std::copy(roles.begin(), roles.end(), array.get());
        *out = handle_of(owned.release());
        *faces = array.release();
        *face_count = roles.size();
        return SM_SUCCESS;
    });
}

void sm_pipe_faces_free(sm_pipe_face* faces) { delete[] faces; }

sm_status sm_brep_shell(const sm_brep* brep, double thickness, const uint32_t* open, size_t open_count,
                        sm_brep** out, sm_shell_face** faces, size_t* face_count) {
    if (!brep || !out || !faces || !face_count || (open_count > 0 && !open)) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!std::isfinite(thickness) || thickness <= 0.0) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    return guarded([&] {
        const SmBrep& source = *brep_of(brep);
        size_t count = 0, unused = 0;
        sm_brep_counts(brep, &count, &unused);
        std::vector<bool> opened(count, false);
        for (size_t i = 0; i < open_count; ++i) {
            if (open[i] >= count || opened[open[i]]) return SM_SHIM_ERR_BAD_FACE;
            opened[open[i]] = true;
        }
        // A closed hollow is a sealed cavity, which SMLib's tessellation can
        // wind the wrong way.
        if (open_count == 0 || open_count == count) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
        SmBrep* shell = nullptr;
        std::vector<sm_shell_face> roles;
        const sm_status status = shell_brep(source, opened, thickness, shell, roles);
        if (status != SM_SUCCESS) return status;
        std::unique_ptr<SmBrep> owned(shell);
        std::unique_ptr<sm_shell_face[]> array(new sm_shell_face[roles.size()]);
        std::copy(roles.begin(), roles.end(), array.get());
        *out = handle_of(owned.release());
        *faces = array.release();
        *face_count = roles.size();
        return SM_SUCCESS;
    });
}

void sm_shell_faces_free(sm_shell_face* faces) { delete[] faces; }

}  // extern "C"
