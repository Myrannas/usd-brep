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
static const sm_status SM_SHIM_ERR_BAD_EDGE = 2006;            /* an edge index is out of range or repeated */
static const sm_status SM_SHIM_ERR_BAD_FILLET = 2007;          /* a radius, section or continuity is invalid */
static const sm_status SM_SHIM_ERR_BROKE_SOLID = 2008;         /* the kernel reported success but left a manifold solid non-manifold */
static const sm_status SM_SHIM_ERR_BAD_FACE = 2020;            /* a face index is out of range */
static const sm_status SM_SHIM_ERR_BAD_LABEL = 2021;           /* a label's role is empty or not NUL-terminated */
static const sm_status SM_SHIM_ERR_NOT_PLANAR = 2022;          /* the operation needs a planar face */
static const sm_status SM_SHIM_ERR_BAD_FRAME = 2023;           /* placement axes aren't unit and perpendicular */
static const sm_status SM_SHIM_ERR_BAD_SWEEP = 2024;           /* a zero sweep distance, or a sweep that isn't a solid */
static const sm_status SM_SHIM_ERR_UNSUPPORTED_FACE = 2025;    /* the face, or a blend on it, can't be moved this way */
static const sm_status SM_SHIM_ERR_CHANGED_TOPOLOGY = 2026;    /* the move would make a face vanish (or appear) */
static const sm_status SM_SHIM_ERR_SPENT_TESSELLATION = 2030; /* the tessellation's mesh was already made */

/* Opaque handles. sm_brep is an SmBrep owned by the caller; sm_mesh is
 * tessellation output copied into flat arrays, independent of the brep. */
typedef struct sm_brep sm_brep;
typedef struct sm_mesh sm_mesh;
/* Surface caches kept from one sm_brep_tessellate_some to the next. */
typedef struct sm_surface_caches sm_surface_caches;

/* A tessellation in progress (sm_brep_tessellation_begin), owned by the
 * caller. */
typedef struct sm_tessellation sm_tessellation;

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

/* Room for a label's role, including its terminating NUL. */
#define SM_LABEL_ROLE_CAPACITY 64

/* A caller's label on a face or edge: the feature that made it, the role it
 * plays in that feature, and an index telling apart entities one feature
 * makes with the same role. The kernel carries labels as attributes
 * (SmAttribute.h): a face or edge that splits passes its label to both
 * pieces, a copy of the brep copies them, and when two labelled faces or
 * edges merge into one, the result keeps the label with the lower feature id
 * (the older one). Entities a fillet or chamfer creates have none.
 *
 * present is 0 for no label; the other fields are then zero. role is UTF-8,
 * NUL-terminated, and not empty. */
typedef struct sm_label {
    int32_t present;
    uint32_t feature;
    uint32_t index;
    char role[SM_LABEL_ROLE_CAPACITY];
} sm_label;

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
 * - kernel_type: SMLib's surface type id (SmSurfTypes.h), for diagnostics.
 * - label: the face's label, if any. */
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
    sm_label label;
} sm_face_info;

/* One B-rep face's geometry, in the order of the brep's own face list, for
 * finding and placing features on faces. Cheap to compute; for a face's area
 * and centroid, see sm_brep_face_centroid.
 *
 * - info: as sm_brep_tessellate reports it, with the label.
 * - bounds_min, bounds_max: an axis-aligned box around the trimmed face. It
 *   follows the face's boundary, so the pieces of a split face have their own
 *   boxes, but may be loose around curves.
 * - planar: nonzero when the face lies on a plane, whatever its surface
 *   class (a planar NURBS surface counts).
 * - normal: for planar faces, the unit normal pointing out of the solid;
 *   zero otherwise. */
typedef struct sm_face_geometry {
    sm_face_info info;
    double bounds_min[3];
    double bounds_max[3];
    int32_t planar;
    double normal[3];
} sm_face_geometry;

/* One B-rep edge, in the order of the brep's own edge list. Its polyline is
 * edge_points[first_point .. first_point + point_count]. face_a and face_b
 * index the face list (a seam joins a face to itself); -1 when the edge has
 * fewer uses. use_count is how many faces use it (2 for a manifold solid).
 * label is the edge's own label, if any. */
typedef struct sm_edge_info {
    uint32_t first_point;
    uint32_t point_count;
    int32_t face_a;
    int32_t face_b;
    uint32_t use_count;
    sm_label label;
} sm_edge_info;

/* How many points sm_edge_geometry samples along each edge. */
#define SM_EDGE_SAMPLES 5

/* One B-rep edge's geometry, in the order of the brep's own edge list, for
 * picking edges by where they are rather than by index.
 *
 * - points: SM_EDGE_SAMPLES points on the edge, evenly spaced in its curve
 *   parameter from one end to the other. points[SM_EDGE_SAMPLES / 2] is the
 *   parametric midpoint.
 * - length: the edge's arc length, to a relative accuracy of 1e-7.
 * - face_a, face_b, use_count: as in sm_edge_info.
 * - normal_a, normal_b: the outward unit normals of face_a and face_b at the
 *   midpoint; zero when the face is absent or the normal can't be found.
 * - label: the edge's own label; label_a, label_b: those of face_a and
 *   face_b. */
typedef struct sm_edge_geometry {
    double points[SM_EDGE_SAMPLES][3];
    double length;
    int32_t face_a;
    int32_t face_b;
    uint32_t use_count;
    double normal_a[3];
    double normal_b[3];
    sm_label label;
    sm_label label_a;
    sm_label label_b;
} sm_edge_geometry;

/* What each face of a swept prism came from (sm_brep_sweep_face). */
enum sm_swept_role {
    /* The copy of the swept face where the sweep starts. */
    SM_SWEPT_START = 0,
    /* The copy where it ends, the swept distance away. */
    SM_SWEPT_END = 1,
    /* A face swept from one edge of the swept face. */
    SM_SWEPT_SIDE = 2
};

/* One face of a swept prism, in the order of the prism's face list. For
 * sides, source_edge is the edge of the swept face it was swept from, as an
 * index into the source brep's edge list, or -1 when the shim couldn't tell. */
typedef struct sm_swept_face {
    int32_t role;
    int32_t source_edge;
} sm_swept_face;

/* Fillet cross-sections (SmFilletSurfaceGeneratorType). */
enum sm_fillet_section {
    /* A flat face between the two lines where a rolling ball of the given
     * radius touches the faces: a chamfer set back by the radius on each face
     * of a 90-degree edge. */
    SM_FILLET_LINEAR = 1,
    /* A circular arc: the round a ball of the given radius leaves. */
    SM_FILLET_CIRCULAR = 2,
    /* A blend curve between the same lines, with the requested continuity. */
    SM_FILLET_BLEND = 3
};

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

/* SmApiCreateCone: a frustum standing on base_center along +z, with the given
 * base and top radii (either may be 0 for a point). */
sm_status sm_brep_create_cone(sm_vec3 base_center, double base_radius, double top_radius, double height,
                              sm_brep **out);

/* SmApiCreateTorus: a ring around the +z axis through center. */
sm_status sm_brep_create_torus(sm_vec3 center, double major_radius, double minor_radius, sm_brep **out);

/* SmApiCreateSphere. */
sm_status sm_brep_create_sphere(sm_vec3 center, double radius, sm_brep **out);

/* SmApiCreateTaperExtrude: extrudes the closed polygon through `points` (in
 * order, in one plane) by `height` along the polygon's normal, its walls
 * leaning in by `draft_deg` (0 to under 90) from that direction, with both
 * ends capped and sharp corners. */
sm_status sm_brep_create_tapered_prism(const sm_vec3 *points, size_t count, double height, double draft_deg,
                                       sm_brep **out);

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

/* Writes the geometry of the brep's first `capacity` edges to out, and the
 * brep's edge count to *count. Pass capacity = *count from sm_brep_counts to
 * get every edge. Edge indices below hold until the brep is next changed. */
sm_status sm_brep_edge_geometry(const sm_brep *brep, sm_edge_geometry *out, size_t capacity, size_t *count);

/* SmApiFilletEdgesPerEdge: fillets the edges at the given indices of the
 * brep's edge list in place. Edge edges[i] gets radii[i] and the
 * sm_fillet_section sections[i]. continuity (1, 2 or 3 for G1, G2 or G3) and
 * thumbweight (1.0 is the kernel default; larger is fuller) shape blend
 * sections; corners where three filleted edges meet are only G1. Returns
 * SM_SHIM_ERR_BROKE_SOLID when the brep was a manifold solid and the kernel
 * reports success without leaving one. On any failure the brep may be partly
 * modified: the caller should discard it. */
sm_status sm_brep_fillet_edges(sm_brep *brep, const uint32_t *edges, const double *radii, const int32_t *sections,
                               size_t count, uint32_t continuity, double thumbweight);

/* Chamfers the edges at the given indices of the brep's edge list in place,
 * as SmApiChamferFillet does for every edge: a constant-distance solver with
 * a flat cross-section, beveling corners where chamfered edges meet. width is
 * the distance across the chamfer face, so a 90-degree edge is set back by
 * width / sqrt(2) on each face. Unlike SmApiChamferFillet it leaves the
 * executive's global merge off, which breaks chamfers on circular edges.
 * Failures are as for sm_brep_fillet_edges. */
sm_status sm_brep_chamfer_edges(sm_brep *brep, const uint32_t *edges, size_t count, double width);

/* Writes the geometry of the brep's first `capacity` faces to out, and the
 * brep's face count to *count, as sm_brep_edge_geometry does for edges. Face
 * indices hold until the brep is next changed. */
sm_status sm_brep_face_geometry(const sm_brep *brep, sm_face_geometry *out, size_t capacity, size_t *count);

/* The area and area centroid of the trimmed face at an index of the brep's
 * face list, integrated to relative_accuracy (clamped to [1e-10, 1e-1];
 * SmApiFaceComputeCentroid stops at 1e-4). On a planar face the centroid lies
 * in the face's plane, though not necessarily inside the face. */
sm_status sm_brep_face_centroid(const sm_brep *brep, size_t face, double relative_accuracy, double *area,
                                double centroid[3]);

/* Writes the brep's axis-aligned bounding box. It encloses the brep but may
 * be loose around curved faces. */
sm_status sm_brep_bounding_box(const sm_brep *brep, double min[3], double max[3]);

/* Labels the face or edge at an index of the brep's list, replacing any label
 * it has; a label with present = 0 removes it. Labels don't change the brep's
 * shape or lists, so indices stay valid. */
sm_status sm_brep_set_face_label(sm_brep *brep, size_t face, const sm_label *label);
sm_status sm_brep_set_edge_label(sm_brep *brep, size_t edge, const sm_label *label);

/* Moves the brep rigidly into a frame (SmApiTransform): the point (x, y, z)
 * goes to origin + x * x_axis + y * y_axis + z * (x_axis × y_axis). The axes
 * must be unit length and perpendicular to within 1e-9. */
sm_status sm_brep_place(sm_brep *brep, sm_vec3 origin, sm_vec3 x_axis, sm_vec3 y_axis);

/* Sweeps the planar face at an index of the brep's face list along its
 * outward normal by distance (negative sweeps into the solid), making a new
 * closed prism (SmApiCreateCurveSweepFromFaces, translational, capped). The
 * brep is unchanged. On success *out is the prism and *faces an array of
 * *face_count entries, one per prism face in the prism's face order, which the
 * caller frees with sm_swept_faces_free. */
sm_status sm_brep_sweep_face(const sm_brep *brep, size_t face, double distance, sm_brep **out,
                             sm_swept_face **faces, size_t *face_count);

/* Frees the array sm_brep_sweep_face returns. NULL is ignored. */
void sm_swept_faces_free(sm_swept_face *faces);

/* Merges faces that lie on one surface into one face, removing the edges
 * between them, and edges that lie on one curve into one edge
 * (SmBrep::RemoveTopologicalEdgesAndVertices). The shape is unchanged. Labels
 * follow sm_label's merge rule. */
sm_status sm_brep_merge_faces(sm_brep *brep);

/* Moves the faces at the given indices of the brep's face list `distance`
 * along their outward normals (negative moves them inward), in place, as a
 * Parasolid face offset does: every other neighbour is extended or trimmed
 * along its own surface to meet them. Picked planes are translated, picked
 * cones slid along their axis, and other picked curved faces replaced by their
 * offsets (SmSurface::CreateOffsetSurface). Blends on a picked plane or cone
 * (fillets, corner blends) are translated by the shortest move that keeps
 * them tangent to all their neighbours. Built on SmBrep::LocalOperation; see
 * shim.cpp.
 *
 * The brep keeps its faces and edges, so a move that would make a face vanish
 * fails, with SM_SHIM_ERR_CHANGED_TOPOLOGY or SM_SHIM_ERR_BROKE_SOLID (the
 * kernel returned an open solid, one whose vertices or edges are off their
 * faces by more than the brep's tolerance, or one whose volume its own
 * tessellation disagrees with). SM_SHIM_ERR_UNSUPPORTED_FACE: a picked curved
 * face other than a cone has blends on it, no translation keeps a blend
 * tangent to all its neighbours, or a surface has no offset.
 * SM_SHIM_ERR_BAD_FACE: an index is out of range or repeated. The brep is
 * unchanged when the kernel itself fails; after a check fails it may be
 * modified, so the caller should discard it. */
sm_status sm_brep_offset_faces(sm_brep *brep, const uint32_t *faces, size_t count, double distance);

/* SmApiTessellate plus SmApiTessellateBoundaries, copied into an sm_mesh
 * with every triangle and vertex linked back to its B-rep face. */
sm_status sm_brep_tessellate(sm_brep *brep, const sm_tessellation_params *params, sm_mesh **out);

/* Fills *out with views of the mesh's arrays. */
void sm_mesh_view_get(const sm_mesh *mesh, sm_mesh_view *out);

/* Deletes a mesh. NULL is ignored. */
void sm_mesh_free(sm_mesh *mesh);

/* sm_brep_tessellate one face at a time. sm_brep_tessellation_begin splits
 * every edge as SmApiTessellate does for the whole brep, on a copy, then
 * sm_tessellation_mesh triangulates only the faces asked for. A face's
 * triangles and vertices come out the same whichever other faces are
 * triangulated, so meshes put together from different tessellations of the
 * same faces stay watertight. The brep must outlive the tessellation and stay
 * unchanged until sm_tessellation_mesh has returned. */
sm_status sm_brep_tessellation_begin(sm_brep *brep, const sm_tessellation_params *params, sm_tessellation **out);

/* The number of faces: the brep's. */
size_t sm_tessellation_face_count(const sm_tessellation *tessellation);

/* The face's key: a 128-bit digest (MurmurHash3) of everything that fixes
 * its triangles, from the tessellation parameters, the face's surface, trim
 * and tolerances, and its boundary as split. Faces with the same key are
 * triangulated the same, vertex for vertex. *keyed is 0 for a face that can't
 * be keyed (SMLib can't write its geometry). Keys can't be had after
 * sm_tessellation_mesh (SM_SHIM_ERR_SPENT_TESSELLATION). */
sm_status sm_tessellation_face_key(sm_tessellation *tessellation, size_t face, uint64_t key[2], int32_t *keyed);

/* The mesh, with triangles and vertices only for the faces `wanted` marks:
 * one byte per face, non-zero to triangulate (NULL triangulates every face).
 * Faces follow each other in face order, each with its own run of vertices,
 * numbered in the order its triangles first use them; each face's triangles
 * and vertices are sm_brep_tessellate's for it. faces and edges cover the
 * whole brep.
 * Once per tessellation (SM_SHIM_ERR_SPENT_TESSELLATION);
 * SM_SHIM_ERR_BAD_FACE when wanted_count isn't the face count. */
sm_status sm_tessellation_mesh(sm_tessellation *tessellation, const uint8_t *wanted, size_t wanted_count,
                               sm_mesh **out);

/* Deletes a tessellation. NULL is ignored. */
void sm_tessellation_free(sm_tessellation *tessellation);

/* Tessellation of some faces only, with SMLib's edge splitting run for just
 * those faces and their neighbours (SmTess::Phase1SetupFaces, a fork patch).
 *
 * sm_brep_tessellation_keys writes each face's key before tessellating: a
 * 128-bit digest of everything the tessellation reads for the face when no
 * face's boundary needs refining. That is the parameters, the brep's
 * tolerance, the face's surface, domain, tolerance and loops (each edgeuse's
 * UV trim curve, edge curve, interval and tolerance, and start vertex), and
 * for each of its edges, every edgeuse's trim curve and the surface, domain,
 * tolerance and order of the other faces using it. keys holds two words per
 * face and keyed one flag per face, 0 for a face that can't be keyed (SMLib
 * can't write its geometry). SM_SHIM_ERR_BAD_FACE when face_count isn't the
 * brep's. These keys are distinct from sm_tessellation_face_key's.
 *
 * sm_brep_tessellate_some triangulates the faces `wanted` marks (one byte per
 * face, non-zero to triangulate), laid out as sm_tessellation_mesh lays them
 * out, each face's triangles and vertices being sm_brep_tessellate's for it,
 * provided every face not wanted has the key of a face an earlier
 * sm_brep_tessellate_some triangulated. *out is NULL when it can't vouch for
 * that: when a wanted face's boundary needs refining, or set-up of one of
 * them fails or can't be matched; sm_brep_tessellation_begin then
 * tessellates the brep. SM_SHIM_ERR_BAD_FACE when wanted_count isn't the face
 * count. The brep isn't changed.
 *
 * With caches (or NULL), a face whose surface is identical to one a cache
 * was kept for takes that cache rather than building its own, and the
 * caches of faces used only to split shared edges are kept. A wanted face
 * takes a copy, since triangulating changes a cache, and the cache stays
 * kept; it keeps a copy of a cache it builds when one of the store's last
 * two calls built one for the same surface. A cache depends only on its
 * surface and the parameters, so the mesh comes out the same. */
sm_status sm_brep_tessellation_keys(sm_brep *brep, const sm_tessellation_params *params, uint64_t *keys,
                                    int32_t *keyed, size_t face_count);
sm_status sm_brep_tessellate_some(sm_brep *brep, const sm_tessellation_params *params, const uint8_t *wanted,
                                  size_t wanted_count, sm_surface_caches *caches, sm_mesh **out);

/* An empty store of surface caches for sm_brep_tessellate_some, which keeps
 * the caches given back by its last two calls that triangulated any face, by
 * surface and parameters, and remembers which surfaces they found no cache
 * for.
 * NULL when out of memory. Use and free it on the kernel thread. */
sm_surface_caches *sm_surface_caches_new(void);

/* The number of caches kept. */
size_t sm_surface_caches_count(const sm_surface_caches *caches);

/* Deletes the store and its caches. NULL is ignored. */
void sm_surface_caches_free(sm_surface_caches *caches);

#ifdef __cplusplus
}
#endif

#endif /* SMLIB_SHIM_H */
