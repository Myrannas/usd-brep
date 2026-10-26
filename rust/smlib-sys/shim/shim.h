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
static const sm_status SM_SHIM_ERR_OBSTRUCTED = 2027;          /* the move runs into features it can't cut through */
static const sm_status SM_SHIM_ERR_UNSUPPORTED_STEP = 2028;    /* a moved face needs a step wall that can't be built */
static const sm_status SM_SHIM_ERR_SPENT_TESSELLATION = 2030; /* the tessellation's mesh was already made */
static const sm_status SM_SHIM_ERR_BAD_AXIS = 2031;            /* an axis isn't finite, has no direction or leaves the face's plane */
static const sm_status SM_SHIM_ERR_AXIS_CROSSES_FACE = 2032;   /* a revolve's axis runs through the inside of its face */
static const sm_status SM_SHIM_ERR_BAD_DRAFT = 2033;           /* a draft angle isn't finite and under 90 degrees either way */
static const sm_status SM_SHIM_ERR_DRAFT_CLOSES = 2034;        /* a draft's sides meet before the end of the sweep */
static const sm_status SM_SHIM_ERR_BAD_PATH = 2035;            /* a pipe's path turns back on itself, or a bend overlaps itself */
static const sm_status SM_SHIM_ERR_ARC_CORNER = 2036;          /* an arc of a pipe's path meets a curve beside it at an angle */
static const sm_status SM_SHIM_ERR_BAD_SHELL = 2037;           /* a shell's input isn't a manifold solid, or its result isn't one */
static const sm_status SM_SHIM_ERR_BAD_SCALE = 2038;           /* a scale factor isn't finite and positive, or a surface or curve has no exact image under it */
static const sm_status SM_SHIM_ERR_NOT_APART = 2039;           /* solids to combine touch, overlap or lie one inside the other */

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
#define SM_LABEL_ROLE_CAPACITY 256

/* A caller's label on a face or edge: the feature that made it, the role it
 * plays in that feature, and an index telling apart entities one feature
 * makes with the same role. The kernel carries labels as attributes
 * (SmAttribute.h): a face or edge that splits passes its label to both
 * pieces, a copy of the brep copies them, and when two labelled faces or
 * edges merge into one, the result keeps the label with the lower feature id
 * (the older one), or of one feature the lower index. Entities a fillet or
 * chamfer creates have none.
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
 * - face_a, face_b, use_count: as in sm_edge_info.
 * - normal_a, normal_b: the outward unit normals of face_a and face_b at the
 *   midpoint; zero when the face is absent or the normal can't be found.
 * - label: the edge's own label; label_a, label_b: those of face_a and
 *   face_b. */
typedef struct sm_edge_geometry {
    double points[SM_EDGE_SAMPLES][3];
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

/* What a curve of sm_brep_create_planar_faces is. */
enum sm_plane_curve_kind {
    /* The segment from a to b (SmApiCreateLineSegment). */
    SM_PLANE_LINE = 0,
    /* The arc round centre a of radius, from start_deg (from x towards y) by
     * sweep_deg: counter-clockwise when positive (SmApiCreateArc). */
    SM_PLANE_ARC = 1,
    /* The whole circle round centre a of radius, a loop on its own
     * (SmApiCreateCircle). */
    SM_PLANE_CIRCLE = 2
};

/* A line, arc or circle in a plane's own (x, y), and the label for its edge.
 * Fields a kind doesn't use are ignored. */
typedef struct sm_plane_curve {
    int32_t kind;
    double a[2];
    double b[2];
    double radius;
    double start_deg;
    double sweep_deg;
    sm_label label;
} sm_plane_curve;

/* Planar faces bounded by closed loops of curves in the plane through origin
 * spanned by x_axis and y_axis (unit and perpendicular, as for
 * sm_brep_place): the curves are made with SmApiCreateLineSegment,
 * SmApiCreateArc and SmApiCreateCircle in the plane's (x, y), and the faces
 * with SmApiCreatePlanarFaces. Curves may come in any order and direction;
 * every end must meet exactly one other curve's end, and loops mustn't cross
 * or touch. Outer loops bound faces, loops inside them holes, and loops inside
 * those faces again. *out is a sheet, not a solid: each face's normal is
 * x_axis × y_axis, and each curve's edge has the curve's label (none when
 * present is 0). SM_ERR_INVALID_INPUT for a bad curve; the kernel's status
 * when the loops don't make faces, though it doesn't notice every crossing. */
sm_status sm_brep_create_planar_faces(const sm_plane_curve *curves, size_t count, sm_vec3 origin, sm_vec3 x_axis,
                                      sm_vec3 y_axis, sm_brep **out);

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

/* The arc length of the edge at an index of the brep's edge list
 * (SmApiEdgeComputeLength), integrated to relative_accuracy (clamped to
 * [1e-10, 1e-1]). It costs far more than sampling the edge, so
 * sm_brep_edge_geometry leaves it out. */
sm_status sm_brep_edge_length(const sm_brep *brep, size_t edge, double relative_accuracy, double *out);

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

/* Moves the brep rigidly (SmApiTransform): the point p goes to R p + t, with
 * R given row by row in rotation. R must be a rotation to within 1e-9: unit
 * columns, perpendicular, the third the cross product of the first two;
 * SM_SHIM_ERR_BAD_FRAME otherwise. Every surface and curve keeps its type
 * and parameterisation, and every face and edge its label. */
sm_status sm_brep_transform(sm_brep *brep, const double rotation[9], sm_vec3 translation);

/* The brep's tolerance (SmBrep::GetTolerance), in model units. */
sm_status sm_brep_tolerance(const sm_brep *brep, double *out);

/* What sm_brep_scale_keeping or sm_brep_mirror made of the brep's faces and
 * edges. */
typedef struct sm_scale_report {
    uint32_t faces_kept;        /* analytic, and still of the same kind */
    uint32_t faces_extruded;    /* cylinders stretched across their axis, now ellipses swept along a line */
    uint32_t faces_nurbs;       /* analytic, now NURBS */
    uint32_t faces_were_nurbs;  /* NURBS before, and still */
    uint32_t curves_kept;       /* lines, and circles or ellipses still of the same kind */
    uint32_t curves_stretched;  /* circles now ellipses */
    uint32_t curves_nurbs;      /* analytic, now NURBS */
    uint32_t curves_were_nurbs; /* NURBS before, and still */
} sm_scale_report;

/* Scales the brep about a point: the point p goes to
 * point + factors * (p - point), axis by axis. The factors must be finite and
 * positive; SM_SHIM_ERR_BAD_SCALE otherwise. Rather than SmBrep::Transform,
 * which turns every surface and curve into a NURBS before a non-uniform
 * scale (SmBrep::TurnToNURBS), this works out each face's surface and each
 * edge's curve anew, keeping it analytic wherever the scale leaves it so:
 * everything under equal factors (an offset surface's distance scaled with
 * it), and otherwise planes and lines, a cylinder or cone along a scale axis
 * whose other two factors are equal, circles in a plane the scale keeps
 * round, and surfaces of extrusion. A circle the scale stretches becomes an
 * ellipse when the images of its axes stay square; with extrusions, a
 * cylinder stretched across its axis becomes an ellipse swept along a line
 * (an SmSurfOfExtrusion). Every other B-spline surface or curve becomes its
 * NURBS, scaled exactly. Each new surface and curve keeps the old one's
 * parameters and lies where the scale takes it to within a thousandth of the
 * brep's tolerance, so faces keep their domains and edges their intervals and
 * UV trim curves; every face and edge keeps its label. The brep's tolerance,
 * and every face's, edge's and vertex's, scale by the largest factor, so every
 * gap stays within the tolerance it was within; the kernel's context and
 * other breps are untouched. SM_SHIM_ERR_BAD_SCALE, with the brep unchanged,
 * when a surface or curve has no exact image this can make (an offset surface
 * under unequal factors, a curve that isn't a B-spline) or two faces or edges
 * share one. report may be NULL; it is written only on success. */
sm_status sm_brep_scale_keeping(sm_brep *brep, sm_vec3 point, sm_vec3 factors, int32_t extrusions,
                                sm_scale_report *report);

/* Reflects the brep across the plane through point with the given normal
 * (any finite length but zero): the point p goes to p - 2 ((p - point) . n) n
 * for the unit normal n. Every surface and curve keeps its kind and its
 * parameters. Planes, circles and ellipses move by the rigid motion that
 * moves their own plane as the reflection does (the reflection after one
 * across that plane), through their own Transform, as a rigid move would.
 * Everything else takes sm_brep_scale_keeping's path with a factor of -1
 * along n: its NURBS is reflected exactly, and cylinders, cones and spheres
 * are found in it again, while tori, surfaces of revolution (SMLib's fillets
 * on circular edges) and of extrusion, which SMLib can't find in a NURBS, are
 * made from their own description, and an offset surface from its reflected
 * base with its distance negated. Rather than SmBrep::Mirror, which turns
 * every curve into a NURBS and rebuilds every UV trim curve. Each face's
 * normal Su x Sv then points into the solid, so the faces' two sides swap
 * shells and the order of faces round each edge reverses, as SmBrep::Mirror
 * does: the result is a valid solid facing out. Faces keep their (u, v)
 * domains, edges their intervals and UV trim curves, and every face and edge
 * its label; tolerances are unchanged. SM_SHIM_ERR_BAD_FRAME for a point or
 * normal that isn't finite, or a zero normal; SM_SHIM_ERR_BAD_SCALE as
 * sm_brep_scale_keeping, and for a face without two sides or an edge with an
 * odd number of uses, with the brep unchanged; SM_SHIM_ERR_BROKE_SOLID if
 * turning its topology round fails after the geometry has changed. report
 * may be NULL; it is written only on success. */
sm_status sm_brep_mirror(sm_brep *brep, sm_vec3 point, sm_vec3 normal, sm_scale_report *report);

/* Writes the brep to an SMLib binary file (SmApiWriteBrepToFile), and reads
 * one back (SmApiReadBrepFromFile, keeping the UV trim curves it holds). */
sm_status sm_brep_write_file(const sm_brep *brep, const char *path);
sm_status sm_brep_read_file(const char *path, sm_brep **out);

/* Sweeps the planar face at an index of the brep's face list along its
 * outward normal by distance (negative sweeps into the solid), making a new
 * closed prism (SmApiCreateCurveSweepFromFaces, translational, capped). The
 * brep is unchanged. On success *out is the prism and *faces an array of
 * *face_count entries, one per prism face in the prism's face order, which the
 * caller frees with sm_swept_faces_free. */
sm_status sm_brep_sweep_face(const sm_brep *brep, size_t face, double distance, sm_brep **out,
                             sm_swept_face **faces, size_t *face_count);

/* Frees the array sm_brep_sweep_face, sm_brep_revolve_face or
 * sm_brep_sweep_face_tapered returns. NULL is ignored. */
void sm_swept_faces_free(sm_swept_face *faces);

/* Revolves the planar face at an index of the brep's face list about the
 * axis through axis_point along axis_direction (any length but zero), by
 * angle_deg: anticlockwise about the axis when positive, up to a whole turn
 * either way. The axis must lie in the face's plane (SM_SHIM_ERR_BAD_AXIS)
 * and may run along the face's boundary or touch it, but not through its
 * inside (SM_SHIM_ERR_AXIS_CROSSES_FACE). The face is swept on a copy of it
 * alone (SmTopologySweep with SmRotationalSweepGeometry); the brep is
 * unchanged. *out is the closed solid, without labels, and *faces says what
 * each of its faces came from, as for sm_brep_sweep_face: the face where the
 * turn starts (SM_SWEPT_START, the face itself) and ends (SM_SWEPT_END),
 * which a whole turn has neither of, and a side for each edge of the face,
 * which the sweep's own maps tie to it. Free *faces with
 * sm_swept_faces_free. */
sm_status sm_brep_revolve_face(const sm_brep *brep, size_t face, sm_vec3 axis_point, sm_vec3 axis_direction,
                               double angle_deg, sm_brep **out, sm_swept_face **faces, size_t *face_count);

/* sm_brep_sweep_face with the sides leaning draft_deg (under 90 either way)
 * from the sweep's direction: positive leans them into the face, so the end
 * is the face shrunk all round (its outer loop in, its holes out), negative
 * leans them out. Each loop is swept by SmPrimitiveCreation::CreateTaperExtrude
 * with sharp corners, the holes' solids cut from the outer loop's. Sides
 * swept from lines are planes and from arcs or circles cones; *faces says
 * what each face came from as sm_brep_sweep_face's does, a side's source edge
 * found from where its edge on the face lies. A side standing on a corner,
 * which an outward draft can add where a line meets an arc, has source edge
 * -1. A draft of 0 is sm_brep_sweep_face. SM_SHIM_ERR_BAD_DRAFT for an angle
 * that isn't finite and under 90 degrees; SM_SHIM_ERR_DRAFT_CLOSES when the
 * sides meet before the end, so the face closes up or a hole swallows the
 * end; the kernel's status, or SM_SHIM_ERR_BAD_SWEEP, when it builds no
 * solid (an outward draft where a line meets an arc at an angle). Free
 * *faces with sm_swept_faces_free. */
sm_status sm_brep_sweep_face_tapered(const sm_brep *brep, size_t face, double distance, double draft_deg,
                                     sm_brep **out, sm_swept_face **faces, size_t *face_count);

/* What a face of a pipe is (sm_brep_create_pipe). */
enum sm_pipe_role {
    /* The outer wall, swept round one curve of the path. */
    SM_PIPE_OUTER = 0,
    /* The inner wall of a hollow pipe, swept round one curve of the path. */
    SM_PIPE_INNER = 1,
    /* The cap at the path's start. */
    SM_PIPE_START = 2,
    /* The cap at its end. */
    SM_PIPE_END = 3
};

/* One face of a pipe, in the order of its face list: its sm_pipe_role, and
 * for walls the index in the path of the curve it was swept along, or -1
 * when the shim couldn't tell (and for caps). */
typedef struct sm_pipe_face {
    int32_t role;
    int32_t segment;
} sm_pipe_face;

/* A pipe: a circle of `radius` swept along a path of lines and arcs in the
 * plane through origin spanned by x_axis and y_axis (unit and perpendicular),
 * as sm_brep_create_planar_faces takes its curves but in order, each
 * starting where the last ended and running its own way (an arc's sweep_deg
 * may be negative); a circle is a closed path on its own. With inner_radius
 * above 0 (and below radius) an open path's pipe is hollow. Each curve gets a piece: a disc
 * square to it at its start, swept along a line (as
 * sm_brep_sweep_face) or turned about an arc's centre (as
 * sm_brep_revolve_face), so walls along lines are cylinders and along arcs
 * tori. Where two curves meet tangentially the pieces meet face to face;
 * where they meet at an angle both are made longer and cut at the plane
 * halfway between them, a mitre: only between lines, since a plane can't
 * mitre a cylinder against a torus. The pieces are joined
 * (SmApiBooleanUnion); a hollow pipe is the solid one less a solid pipe of
 * inner_radius. An open path is capped at both ends, a closed one isn't.
 * *faces says what each face is (sm_pipe_face); free it with
 * sm_pipe_faces_free. SM_ERR_INVALID_INPUT for an empty path or bad radii;
 * SM_SHIM_ERR_BAD_PATH when the path turns back on itself (by more than 179
 * degrees), a bend would overlap itself, or a closed path's pipe is hollow
 * (its bore would be sealed inside); SM_SHIM_ERR_ARC_CORNER when an arc
 * meets a curve beside it at an angle; the kernel's status, or
 * SM_SHIM_ERR_BAD_SWEEP, when it builds no solid. */
sm_status sm_brep_create_pipe(const sm_plane_curve *path, size_t count, sm_vec3 origin, sm_vec3 x_axis,
                              sm_vec3 y_axis, double radius, double inner_radius, sm_brep **out, sm_pipe_face **faces,
                              size_t *face_count);

/* Frees the array sm_brep_create_pipe returns. NULL is ignored. */
void sm_pipe_faces_free(sm_pipe_face *faces);

/* What a face of a shelled solid is (sm_brep_shell), by the face of the
 * solid it came from. */
enum sm_shell_role {
    /* A face that stays where it was, outside. */
    SM_SHELL_OUTER = 0,
    /* What's left of a face opened, round the opening. */
    SM_SHELL_RIM = 1,
    /* A face's offset, the wall's inside. */
    SM_SHELL_INNER = 2,
    /* A face from no face of the solid. */
    SM_SHELL_NEW = 3
};

/* One face of a shelled solid, in the order of its face list: its
 * sm_shell_role, and the index in the solid's face list of the face it came
 * from (-1 for SM_SHELL_NEW). */
typedef struct sm_shell_face {
    int32_t role;
    int32_t source_face;
} sm_shell_face;

/* The manifold solid hollowed to a wall `thickness` thick inward from its
 * faces, the faces at the given indices of its face list left open
 * (SmApiShellBrepFull with extended corners and self-intersections trimmed:
 * SmOffsetExecutive::ShellBrep insets the solid, the open faces left where
 * they are, and cuts the inset from it). The brep is unchanged. *faces says
 * what each face of the result is (sm_shell_face): faces that stay keep
 * their labels, and so do edges only they use; inner faces have none.
 * SM_ERR_INVALID_INPUT for a thickness that isn't finite and positive, or
 * for opening no face (a closed hollow is a sealed cavity, which SMLib's
 * tessellation can wind the wrong way) or every face;
 * SM_SHIM_ERR_BAD_FACE for an index out of range or repeated;
 * SM_SHIM_ERR_BAD_SHELL when the brep or the result isn't a manifold solid,
 * or the result's faces' integrals disagree with their edges; the kernel's
 * status when it can't offset the faces. Free *faces with
 * sm_shell_faces_free. */
sm_status sm_brep_shell(const sm_brep *brep, double thickness, const uint32_t *open, size_t open_count,
                        sm_brep **out, sm_shell_face **faces, size_t *face_count);

/* Frees the array sm_brep_shell returns. NULL is ignored. */
void sm_shell_faces_free(sm_shell_face *faces);

/* Merges faces that lie on one surface into one face, removing the edges
 * between them, and edges that lie on one curve into one edge
 * (SmBrep::RemoveTopologicalEdgesAndVertices). The shape is unchanged. Labels
 * follow sm_label's merge rule. */
sm_status sm_brep_merge_faces(sm_brep *brep);

/* How two solids lie (SmApiBrepRelationship): relationship is
 * SM_RELATIONSHIP_APART when their surfaces are more than both breps'
 * tolerances apart and neither lies in the other's material, distance their
 * closest approach (SmTopologySolver::BrepBrepSolve); SM_RELATIONSHIP_CONTAINS
 * or SM_RELATIONSHIP_INSIDE when a lies all round b or b round a, distance
 * the gap between their walls; otherwise their surfaces meet, distance 0, and
 * SM_RELATIONSHIP_OVERLAPPING when a sample of either's vertices and face
 * points lies inside the other, SM_RELATIONSHIP_TOUCHING when none does (so
 * an overlap every sample misses reads as touching). A solid in a cavity of
 * the other is apart. */
static const int32_t SM_RELATIONSHIP_APART = 0;
static const int32_t SM_RELATIONSHIP_TOUCHING = 1;
static const int32_t SM_RELATIONSHIP_OVERLAPPING = 2;
static const int32_t SM_RELATIONSHIP_CONTAINS = 3;
static const int32_t SM_RELATIONSHIP_INSIDE = 4;
sm_status sm_brep_relationship(const sm_brep *a, const sm_brep *b, int32_t *relationship, double *distance);

/* Adds a copy of other's solids to the brep as lumps of their own, labels
 * and all, without intersecting them with the brep's (SmBrep::MergeBrep,
 * keeping surface types). For solids that neither touch nor overlap, such as
 * a pattern's copies, it's their union without a boolean's work. Unless their
 * bounding boxes are apart by more than both tolerances, it asks whether they
 * are apart as sm_brep_relationship does (their surfaces' closest approach,
 * searched no further than both tolerances, then whether a point of each lies
 * outside the other), and refuses solids that aren't with
 * SM_SHIM_ERR_NOT_APART, the brep unchanged, since merging them would make a
 * brep that isn't a valid solid. A solid in a cavity of the other isn't apart
 * here. The closest approach is a lower bound on round faces (see
 * Brep::relationship in the Rust crate), so copies of a cylinder closer than
 * about half its radius can be refused though apart. other is unchanged; it
 * may not be brep itself. */
sm_status sm_brep_combine(sm_brep *brep, const sm_brep *other);

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
 * A planar face pushed in can run into features standing on the faces
 * beside it (a boss or a hole on a lid near the edge above a pushed wall):
 * their loops would cross their face's new boundary. Each is set aside, the
 * faces are offset without it, and it comes back cut by the pushed face's
 * new plane: a boss joined back as far as the plane leaves it, a hole cut
 * back. SM_SHIM_ERR_OBSTRUCTED when one can't be set aside: it doesn't stand
 * on a planar face entirely on one side of it, or something else stands in
 * its way (see sm_brep_offset_obstacles). A face pushed in past a chamfer
 * beside it (a lid past the chamfer along its edge, a hole widened or a post
 * thinned past the chamfer on its rim) swallows the chamfer and meets the
 * plane beyond directly, as if the edge had never been chamfered, where it
 * runs along that plane's normal and nothing else stands near the chamfer.
 * Pushing through merges back only the faces its own cuts split: an edge
 * that already joined two faces on one surface (an imprint's) stays.
 * A picked face beside an unpicked face on its own surface (two pieces of
 * one face, as an imprint splits it) steps away from it instead: the picked
 * faces are swept along the move into prisms, which are joined (pulled out)
 * or cut (pushed in), so a wall stands along each edge to such a face and
 * every other neighbour stretches or is trimmed along its surface as
 * before. The moved faces keep their labels, and each wall takes the label
 * of the edge it stands on, if that has one. That
 * needs every picked face planar, all facing one way, and every other
 * neighbour standing along the move (a plane containing it, or a cylinder
 * along it): SM_SHIM_ERR_UNSUPPORTED_STEP otherwise. SM_SHIM_ERR_OBSTRUCTED
 * when the prisms don't add or take exactly their own volume, as when a
 * face pushed in runs into a pocket under it. A step never pushes through
 * anything in its way.
 *
 * Otherwise the brep keeps its faces and edges, so a move that would make a
 * face vanish fails, with SM_SHIM_ERR_CHANGED_TOPOLOGY or
 * SM_SHIM_ERR_BROKE_SOLID (the kernel
 * returned an open solid, a face whose loops cross or whose inner loop left
 * its outer one, one whose vertices or edges are off their faces by more
 * than the brep's tolerance, or one whose volume its own tessellation
 * disagrees with). SM_SHIM_ERR_UNSUPPORTED_FACE: a picked curved
 * face other than a cone has blends on it, no translation keeps a blend
 * tangent to all its neighbours, or a surface has no offset.
 * SM_SHIM_ERR_BAD_FACE: an index is out of range or repeated. The brep is
 * unchanged when the kernel itself fails; after a check fails it may be
 * modified, so the caller should discard it. */
sm_status sm_brep_offset_faces(sm_brep *brep, const uint32_t *faces, size_t count, double distance);

/* sm_brep_offset_faces's local operation alone: the kernel's own result, not
 * checked and not pushed through anything in the way, for diagnosing the
 * kernel and for testing what callers do with a broken solid. */
sm_status sm_brep_offset_faces_unchecked(sm_brep *brep, const uint32_t *faces, size_t count, double distance);

/* What an offset of the faces by `distance` runs into (see
 * sm_brep_offset_faces): the features standing on faces beside the moving
 * ones, along inner loops, that reach into the space a picked plane pushed
 * in sweeps. Writes up to `capacity` pairs of a face's index and its
 * obstacle's number (from 0, in the brep's face order) and sets `count` to
 * the number of pairs; none when nothing is in the way. */
sm_status sm_brep_offset_obstacles(const sm_brep *brep, const uint32_t *faces, size_t count, double distance,
                                   uint32_t *out_faces, uint32_t *out_obstacles, size_t capacity,
                                   size_t *out_count);

/* What a face is to sm_brep_move_faces. */
enum sm_move_role {
    SM_MOVE_KEPT = 0,     /* untouched */
    SM_MOVE_MOVED = 1,    /* picked: moved rigidly */
    SM_MOVE_CARRIED = 2,  /* a blend on a picked face, slid to stay tangent to its neighbours */
    SM_MOVE_FOLLOWED = 3, /* a neighbour extended or trimmed along its own surface to meet them */
};

/* Moves the faces at the given indices of the brep's face list rigidly, in
 * place: turned by `angle_deg` about the axis through `axis_point` along
 * `axis_direction` (anticlockwise looking down the axis from its tip; no turn
 * when the angle is 0, and the axis is then ignored), then translated by
 * `translation`. Every other neighbour is extended or trimmed along its own
 * surface to meet them, as sm_brep_offset_faces does, on the same local
 * operation. A translation carries the blends on the picked faces
 * (fillets, corner blends), each by the shortest slide that keeps it tangent
 * to all its neighbours; a blend the move leaves in place stretches. A turn
 * can't carry a blend: SM_SHIM_ERR_UNSUPPORTED_FACE when one is on a turned
 * face (and when no slide keeps a blend tangent).
 *
 * The brep keeps its faces and edges, and their order and labels, so a move
 * that would make a face vanish or appear fails, with
 * SM_SHIM_ERR_CHANGED_TOPOLOGY or SM_SHIM_ERR_BROKE_SOLID, and so does a
 * result that isn't a manifold solid, has gaps or a volume that can't be
 * trusted (see sm_brep_offset_faces): its neighbours can't reach the moved
 * faces. Nothing is pushed through features in the way.
 * SM_SHIM_ERR_UNSUPPORTED_STEP when a picked face is beside an unpicked
 * piece of its own surface (an imprint's), which can't stretch to meet it.
 * SM_SHIM_ERR_BAD_AXIS for a turn about an axis that isn't finite or has no
 * direction; SM_SHIM_ERR_BAD_FACE for an index out of range or repeated, or
 * a role count that isn't the face count. On success writes each face's
 * sm_move_role, in face order, to out_roles. The brep is unchanged when the
 * kernel itself fails; after a check fails it may be modified, so the
 * caller should discard it. */
sm_status sm_brep_move_faces(sm_brep *brep, const uint32_t *faces, size_t count, sm_vec3 translation,
                             sm_vec3 axis_point, sm_vec3 axis_direction, double angle_deg, uint8_t *out_roles,
                             size_t role_count);

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

/* An edge of an imprinted brep that lies along the tool's boundary: its
 * index in the brep's edge list, whether the imprint made it (is_new 1) or
 * it, or the edge it was split from, was there before (0), and the labels of
 * the tool's faces it lies on (sm_label.present is 0 for none). One label
 * where it lies inside a tool face; two where it runs along one of the
 * tool's edges; none where those faces have none. */
typedef struct sm_imprint_edge {
    size_t edge;
    int32_t is_new;
    sm_label tool_faces[2];
} sm_imprint_edge;

/* Imprints `tool` on `brep`: splits brep's faces and edges along the curves
 * where tool's boundary meets them, which become edges and vertices of
 * brep, and changes nothing else (SmMerge::ManifoldBoolean with
 * SM_BO_IMPRINT, keeping the tool). There is no post-process merge, so
 * pieces of a face on one surface stay apart. The shape is unchanged. Split
 * faces and edges keep their labels on every piece (see sm_label); the new
 * edges have none.
 *
 * *edges lists, by edge index, every edge of the result that lies along the
 * tool, new or already there, with the tool faces it lies on (see
 * sm_imprint_edge).
 *
 * With classify non-zero, the imprint goes on to sort the faces as an
 * intersection would, without deleting any (SM_BO_INTERSECTION with
 * bImprintAndClassifyFaces): *inside has one byte per face of the result,
 * non-zero for a face inside the tool, and *face_count is the face count.
 * Otherwise *inside is NULL and *face_count 0.
 *
 * Ownership as for the booleans: on success *out is brep (modified in
 * place) and tool has been deleted; on failure both are still allocated
 * but may have been modified, and the caller frees both. Free *edges and
 * *inside together with sm_imprint_free. */
sm_status sm_brep_imprint(sm_brep *brep, sm_brep *tool, int32_t classify, sm_brep **out, sm_imprint_edge **edges,
                          size_t *edge_count, uint8_t **inside, size_t *face_count);

/* Frees the arrays sm_brep_imprint returns. NULLs are ignored. */
void sm_imprint_free(sm_imprint_edge *edges, uint8_t *inside);

/* sm_brep_merge_faces, except that the edges at the given indices of the
 * brep's edge list stay, whether or not the faces either side lie on one
 * surface (SmBrep::RemoveTopologicalEdgesAndVertices with pOptKeepEdges).
 * SM_SHIM_ERR_BAD_EDGE when an index is out of range or repeated. */
sm_status sm_brep_merge_faces_keeping(sm_brep *brep, const uint32_t *keep, size_t count);

#ifdef __cplusplus
}
#endif

#endif /* SMLIB_SHIM_H */
