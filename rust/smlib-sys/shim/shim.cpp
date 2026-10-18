// C shim over SMLib / SM_API for the smlib-sys Rust crate. See shim.h.

#include "shim.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <new>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <SmApiBrep.h>
#include <SmApiCurves.h>
#include <SmApiGeneral.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmApiSurfaces.h>

#include <SmAxis2Placement.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmCone.h>
#include <SmCylinder.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmFaceuse.h>
#include <SmFilletCorner.h>
#include <SmFilletExecutive.h>
#include <SmFilletSolver.h>
#include <SmFilletStandardSolver.h>
#include <SmPlane.h>
#include <SmPoly.h>
#include <SmSphere.h>
#include <SmSurfOfRevolution.h>
#include <SmSurfTypes.h>
#include <SmTorus.h>

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

void describe_placement(sm_face_info& info, const SmAxis2Placement& placement) {
    put3(info.origin, placement.GetOrigin());
    put3(info.axis, placement.GetZAxis());
    put3(info.ref_dir, placement.GetXAxis());
}

sm_face_info describe_face(const SmFace* face) {
    sm_face_info info;
    std::memset(&info, 0, sizeof info);
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
    for (int k = 0; k < SM_EDGE_SAMPLES; ++k) {
        SmVector3d point;
        SmVector3d* point_out = &point;
        SmVector3d* no_derivative = nullptr;
        SmVector3d* no_second = nullptr;
        const double t = t0 + (t1 - t0) * k / (SM_EDGE_SAMPLES - 1);
        status = SmApiEvaluateCurve(edge->GetCurve(), t, point_out, no_derivative, no_second);
        if (status != SM_SUCCESS) return status;
        put3(info.points[k], point);
    }
    status = SmApiEdgeComputeLength(edge, 1e-9, info.length);
    if (status != SM_SUCCESS) return status;

    const EdgeFaces faces = faces_of(edge);
    info.use_count = faces.use_count;
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
sm_status chamfer(SmBrep& target, SmTArray<SmEdge*>& picked, double width) {
    const SmContext* context = SmApiGetOrCreateContext();
    SmLinearCrossSectionFSG linear;
    SmMakeSurfaceBlendSIH self_intersections(1.0);
    SmFilletExecutive executive(*context, &target);
    executive.SetSelfIntersectionHandler(&self_intersections);
    executive.SetDoGlobalMerge(TRUE);
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
