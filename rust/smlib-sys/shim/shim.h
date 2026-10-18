/*
 * C interface over SMLib / SM_API for the smlib-sys Rust crate.
 *
 * Every function that can fail returns the kernel's status code
 * (SM_SUCCESS = 1000, see source/SMLib/inc/SmMessages.h) or one of the
 * SM_SHIM_ERR_* codes below. Out-pointers are written only on success.
 *
 * Threading: SMLib keeps process-global state and is not thread-safe. Call
 * sm_context_create() once at startup, then make every call from one thread.
 *
 * bindgen runs over this header only; never over the kernel's own headers.
 */
#ifndef SMLIB_SHIM_H
#define SMLIB_SHIM_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* SmApiStatus is a long. */
typedef long sm_status;

static const sm_status SM_STATUS_SUCCESS = 1000;

/* Shim-level failures, outside the kernel's own range. */
static const sm_status SM_SHIM_ERR_NULL_ARGUMENT = 2001;       /* a required pointer argument was NULL */
static const sm_status SM_SHIM_ERR_EXCEPTION = 2002;           /* the kernel threw a C++ exception */
static const sm_status SM_SHIM_ERR_OUT_OF_MEMORY = 2003;       /* std::bad_alloc while building output */
static const sm_status SM_SHIM_ERR_UNMAPPED_FACE = 2004;       /* a mesh polygon's source face isn't in the brep */
static const sm_status SM_SHIM_ERR_UNSUPPORTED_POLYGON = 2005; /* a polygon has holes or fewer than 3 sides */

/* Opaque handles. sm_brep is an SmBrep owned by the caller; sm_mesh is
 * tessellation output copied into flat arrays, independent of the brep. */
typedef struct sm_brep sm_brep;
typedef struct sm_mesh sm_mesh;

typedef struct sm_vec3 {
    double x;
    double y;
    double z;
} sm_vec3;

/* Mirrors SmTessellationParams (SmApiTessellationParams.h). 0 disables a
 * control. Chord height and max edge length are in model units. */
typedef struct sm_tessellation_params {
    double chord_height_tol;
    double curve_angle_tol_deg;
    double surface_angle_tol_deg;
    double max_edge_length;
    double max_aspect_ratio;
    /* Angle tolerance for the edge polylines (SmApiTessellateBoundaries). */
    double edge_angle_tol_deg;
} sm_tessellation_params;

/* The kind of surface a B-rep face lies on. */
enum sm_surface_kind {
    SM_SURFACE_PLANE = 0,
    SM_SURFACE_CYLINDER = 1,
    SM_SURFACE_CONE = 2,
    SM_SURFACE_SPHERE = 3,
    SM_SURFACE_TORUS = 4,
    SM_SURFACE_NURBS = 5,
    /* Any other surface class, such as an offset or swept surface. */
    SM_SURFACE_OTHER = 6
};

/* One B-rep face, in the order of the brep's own face list. Which fields
 * apply depends on kind; unused fields are zero.
 *
 * - origin, axis, ref_dir: the surface's placement. For a plane, origin is
 *   on the plane, axis is its normal and ref_dir its u direction. For
 *   cylinders, cones, spheres and tori, origin is the centre (for cones and
 *   cylinders, the axis point at v = 0), axis is the axis of revolution and
 *   ref_dir the direction of angle zero.
 * - radius: cylinder, sphere, cone radius at origin, torus major radius.
 * - minor_radius: torus minor radius.
 * - half_angle_rad: cone semi-angle, signed: positive widens along axis.
 * - degree_u/v, poles_u/v, rational: NURBS degree, control-point counts and
 *   whether it has weights. Planes are NURBS-backed in SMLib but report
 *   SM_SURFACE_PLANE.
 * - uv_min, uv_max: the face's domain in its surface's (u, v).
 * - kernel_type: SMLib's surface type id (SmSurfTypes.h), for diagnostics. */
typedef struct sm_face_info {
    int32_t kind;
    int32_t kernel_type;
    double origin[3];
    double axis[3];
    double ref_dir[3];
    double radius;
    double minor_radius;
    double half_angle_rad;
    uint32_t degree_u;
    uint32_t degree_v;
    uint32_t poles_u;
    uint32_t poles_v;
    int32_t rational;
    double uv_min[2];
    double uv_max[2];
} sm_face_info;

/* One B-rep edge, in the order of the brep's own edge list. Its polyline is
 * edge_points[first_point .. first_point + point_count]. face_a and face_b
 * index the face list (a seam joins a face to itself); -1 when the edge has
 * fewer uses. use_count is how many faces use it (2 for a manifold solid). */
typedef struct sm_edge_info {
    uint32_t first_point;
    uint32_t point_count;
    int32_t face_a;
    int32_t face_b;
    uint32_t use_count;
} sm_edge_info;

/* Borrowed view of an sm_mesh's arrays, valid until sm_mesh_free.
 *
 * Vertices are split per B-rep face (and either side of a seam), so every
 * vertex belongs to one face. Triangles wind counter-clockwise seen from
 * outside the solid. */
typedef struct sm_mesh_view {
    const float *positions;   /* 3 per vertex */
    const float *normals;     /* 3 per vertex: unit surface normal from the kernel */
    const float *surface_uvs; /* 2 per vertex: (u, v) on the face's surface */
    size_t vertex_count;
    const uint32_t *indices;  /* 3 per triangle */
    const uint32_t *face_ids; /* 1 per triangle: index into faces */
    size_t triangle_count;
    const sm_face_info *faces;
    size_t face_count;
    const sm_edge_info *edges;
    size_t edge_count;
    const float *edge_points; /* 3 per point */
    size_t edge_point_count;
    /* Kernel polygons with more than three sides, fan-triangulated here. */
    size_t polygons_split;
    /* Polygon corners the kernel gave no (u, v) and normal for. */
    size_t corners_missing_uv;
} sm_mesh_view;

/* Creates SMLib's process-global context (SmApiCreateContext). Call once,
 * before any other function and before starting other threads. */
sm_status sm_context_create(void);

/* SmApiCreateBox: origin is the minimum corner; length, width and height run
 * along x, y and z. */
sm_status sm_brep_create_box(sm_vec3 origin, double length, double width, double height, sm_brep **out);

/* SmApiCreateCylinder: base_center is the centre of the base; the axis is +z. */
sm_status sm_brep_create_cylinder(sm_vec3 base_center, double radius, double height, sm_brep **out);

/* SmApiCreateSphere. */
sm_status sm_brep_create_sphere(sm_vec3 center, double radius, sm_brep **out);

/* SmApiBrepCopy: a deep copy. */
sm_status sm_brep_copy(const sm_brep *brep, sm_brep **out);

/* Booleans (SmApiBooleanUnion / Difference / Intersection). On success *out
 * is a (modified in place) and b has been deleted: the caller owns *out and
 * must not touch a or b again. On failure both a and b are still allocated
 * but may have been modified: the caller should free both with
 * sm_brep_free. */
sm_status sm_brep_union(sm_brep *a, sm_brep *b, sm_brep **out);
sm_status sm_brep_difference(sm_brep *a, sm_brep *b, sm_brep **out);
sm_status sm_brep_intersection(sm_brep *a, sm_brep *b, sm_brep **out);

/* SmApiBrepIsManifoldSolid. */
sm_status sm_brep_is_manifold_solid(const sm_brep *brep, int32_t *out);

/* SmApiBrepComputeVolume; relative_accuracy is clamped to [1e-8, 1e-1]. */
sm_status sm_brep_volume(const sm_brep *brep, double relative_accuracy, double *out);

/* Face and edge counts of the brep's own lists. */
sm_status sm_brep_counts(const sm_brep *brep, size_t *faces, size_t *edges);

/* Deletes a brep. NULL is ignored. */
void sm_brep_free(sm_brep *brep);

/* SmApiTessellate plus SmApiTessellateBoundaries, copied into an sm_mesh
 * with every triangle and vertex linked back to its B-rep face. */
sm_status sm_brep_tessellate(sm_brep *brep, const sm_tessellation_params *params, sm_mesh **out);

/* Fills *out with views of the mesh's arrays. */
void sm_mesh_view_get(const sm_mesh *mesh, sm_mesh_view *out);

/* Deletes a mesh. NULL is ignored. */
void sm_mesh_free(sm_mesh *mesh);

#ifdef __cplusplus
}
#endif

#endif /* SMLIB_SHIM_H */
