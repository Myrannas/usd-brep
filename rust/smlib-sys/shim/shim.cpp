// C shim over SMLib / SM_API for the smlib-sys Rust crate. See shim.h.

#include "shim.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <new>
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
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmApiSurfaces.h>
#include <SmApiTrimmedSurfaces.h>

#include <SmAttribute.h>
#include <SmAxis2Placement.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmCone.h>
#include <SmCopyBrepMap.h>
#include <SmCylinder.h>
#include <SmDatabaseIO.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmFilletCorner.h>
#include <SmFilletExecutive.h>
#include <SmFilletSolver.h>
#include <SmFilletStandardSolver.h>
#include <SmGraphicsExtern.h>
#include <SmLine.h>
#include <SmLoopuse.h>
#include <SmOffsetExecutive.h>
#include <SmOffsetGeometryCreation.h>
#include <SmPlane.h>
#include <SmPoly.h>
#include <SmPrimitiveCreation.h>
#include <SmSphere.h>
#include <SmSurfOfRevolution.h>
#include <SmSurfTypes.h>
#include <SmTess.h>
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

// Runs a local operation and holds it to what it promises: the same faces and
// edges, still a manifold solid, every vertex and edge on its faces to the
// brep's tolerance, and a volume its own tessellation agrees with. The offset
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
};

using Motions = std::unordered_map<const SmFace*, Motion>;

Motion motion_of(const Motions& motions, const SmFace* face) {
    const auto found = motions.find(face);
    return found == motions.end() ? Motion() : found->second;
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

    const double slack = tolerance * (1.0 + std::fabs(distance));
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
                if (moving.kind == Motion::kShift && moving.shift.Length() > 0.0 &&
                    invariant_under(other->GetSurface(), moving.shift / moving.shift.Length())) {
                    continue;
                }
                if (moving.kind == Motion::kOffset) return SM_SHIM_ERR_UNSUPPORTED_FACE;
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
                if (motion.kind == Motion::kOffset) return SM_SHIM_ERR_UNSUPPORTED_FACE;
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
            if (motion.kind == Motion::kOffset) return SM_SHIM_ERR_UNSUPPORTED_FACE;
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
    if (motion.kind == Motion::kOffset) return false;
    const SmVector3d difference = shift - motion.shift;
    const double length = difference.Length();
    return length <= 1e-12 || invariant_under(face->GetSurface(), difference / length);
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
    // The v the edge moves rigidly by, tried as each of its faces' shifts.
    bool rigid_shift(const SmEdge* edge, SmVector3d& out) const {
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
                rigid = rigid && lands_like(face, motion_of(m_motions, face), candidate.shift);
            }
            if (rigid) {
                out = candidate.shift;
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

    // Puts the original edge moved by `shift` on every copy of `edge` whose
    // face didn't itself move by `shift`.
    bool move_copies(const SmContext& context, SmOffsetExecutive* executive, SmEdge* edge, const SmVector3d& shift) {
        SmAxis2Placement placement;
        placement.Translate(shift);
        return replace_copies(
            context, executive, edge,
            [&](const SmFace* face) { return (motion_of(m_motions, face).shift - shift).Length() > 1e-12; },
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
            SmVector3d shift;
            const bool rigid = rigid_shift(edges[i], shift);
            if (rigid ? !move_copies(context, executive, edges[i], shift) : !stretch_copies(context, executive, edges[i])) {
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
                lands = lands_like(around[k], motion_of(m_motions, around[k]), candidate.shift);
            }
            if (lands) {
                out = old + candidate.shift;
                return true;
            }
        }
        std::vector<Row> rows;
        for (ULONG j = 0; j < around.GetSize(); ++j) {
            SmPoint3d point;
            SmVector3d normal;
            if (!plane_of(around[j]->GetSurface(), point, normal)) continue;
            const SmVector3d shift = motion_of(m_motions, around[j]).shift;
            rows.push_back({normal, normal.Dot(point + shift)});
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
                if (!lands_like(sides.faces[0], motion_of(m_motions, sides.faces[0]), candidate.shift) ||
                    !lands_like(sides.faces[1], motion_of(m_motions, sides.faces[1]), candidate.shift)) {
                    continue;
                }
                SmVector3d across = std::fabs(direction.x) < 0.9 ? SmVector3d(1.0, 0.0, 0.0) : SmVector3d(0.0, 1.0, 0.0);
                across = across - direction * across.Dot(direction);
                across = across / across.Length();
                const SmVector3d other = direction * across;  // cross product
                const SmPoint3d through = old + candidate.shift;
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

sm_status sm_brep_offset_faces(sm_brep* brep, const uint32_t* faces, size_t count, double distance) {
    if (!brep || (count > 0 && !faces)) return SM_SHIM_ERR_NULL_ARGUMENT;
    if (!std::isfinite(distance) || count == 0) return static_cast<sm_status>(SM_ERR_INVALID_INPUT);
    return guarded([&] {
        SmBrep* target = brep_of(brep);
        SmTArray<SmFace*> all;
        target->GetFaces(all);
        SmTArray<SmFace*> picked;
        std::vector<bool> taken(all.GetSize(), false);
        for (size_t i = 0; i < count; ++i) {
            if (faces[i] >= all.GetSize() || taken[faces[i]]) return SM_SHIM_ERR_BAD_FACE;
            taken[faces[i]] = true;
            picked.Add(all[faces[i]]);
        }
        if (distance == 0.0) return SM_SUCCESS;
        return keeping_topology(*target, [&] { return offset_faces(*target, picked, distance); });
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
        if (found == kept.end()) return nullptr;
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
