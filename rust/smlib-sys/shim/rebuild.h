/* Rebuilding faces for the smlib-sys shim (rebuild.cpp): faces taken out of
 * a solid and surfaces built in the gaps they leave from the rim's edges,
 * sewn back into a solid in one call, so an open body never leaves the
 * shim. Also the measures that check how the surfaces meet: the normal
 * curvature across each edge, and surfaces carried on past their edges.
 *
 * Kept apart from shim.h while it's a spike: its declarations are written
 * out by hand in smlib-sys's `rebuild` module rather than by bindgen. */

#ifndef SMLIB_SHIM_REBUILD_H
#define SMLIB_SHIM_REBUILD_H

#include "shim.h"

#ifdef __cplusplus
extern "C" {
#endif

static const sm_status SM_SHIM_ERR_NOT_RIM = 2101;        /* a side's edge isn't between a face taken out and a kept one */
static const sm_status SM_SHIM_ERR_OPEN_LOOP = 2102;      /* a surface's sides don't run end to end round a loop */
static const sm_status SM_SHIM_ERR_BAD_TOOL = 2103;       /* a tool, a side count or an order the tool can't take */
static const sm_status SM_SHIM_ERR_SURFACE_FAILED = 2104; /* the kernel couldn't make a surface; the report says which */
static const sm_status SM_SHIM_ERR_LEFT_OPEN = 2105;      /* the stitch left edges open; the report says how far apart */
static const sm_status SM_SHIM_ERR_NOT_SOLID = 2106;      /* the stitched result isn't a manifold solid with a volume */

/* The kernel call that makes a surface. */
enum {
    /* SmBSplineSurface::CreateNSidedPatch through 3 to 8 sides in order
     * round the loop, a Tangent side's cross-derivative curve sampled from
     * the kept face beyond it (SmTangentField, fitted with matching
     * parameters to tangent_tol_deg). It makes a patch a side, which meet
     * one another in the middle. Sides next to one another on the same
     * kept face, meeting it alike, are one side: their curves joined. */
    SM_REBUILD_FILL = 0,
    /* SmBSplineSurface::CreateCoonsPatch, bicubic, through 4 sides in
     * order, joined as a fill's are: a Tangent side's cross-derivative as a
     * fill's, a Position side's a cubic from the side before it to the side
     * after, leaving each corner at the twist the field beside it has
     * there, so the patch's corners agree. */
    SM_REBUILD_COONS = 1,
    /* SmBrep::CreateCornerBlend from the rim's edgeuses, 3 or 4 sides, any
     * order: a bicubic Coons patch (with a point for a fourth side) whose
     * cross-tangents SmTangentField fits to tangent_tol_deg. Every side
     * Tangent but a 4-sided one's fourth, which may be Position. */
    SM_REBUILD_CORNER = 2,
    /* SmBrep::CreateEdgeEdgeBlend between sides[0] and sides[1], `order`
     * 0 ruled, 1 G1, 2 G2, 3 G3, both ends alike. */
    SM_REBUILD_EDGE_BLEND = 3,
    /* SmBSplineSurface::CreateSkinnedSurface through sides[0] and
     * sides[1], cubic, a Tangent end taking its cross-derivatives from a
     * derivative surface of the kept face beyond it
     * (SmDerivSurfDefinition, towards the other side), G1 at most. */
    SM_REBUILD_SKIN = 4
};

/* How a side meets the face beyond it. */
enum { SM_MEETS_POSITION = 0, SM_MEETS_TANGENT = 1, SM_MEETS_CURVATURE = 2 };

typedef struct sm_rebuild_side {
    uint32_t edge; /* a rim edge, by its index in the brep given, before any face is taken out */
    int32_t meets;
} sm_rebuild_side;

typedef struct sm_rebuild_surface {
    int32_t tool;
    const sm_rebuild_side *sides;
    size_t side_count;
    int32_t order;          /* an edge blend's: 0 to 3 */
    double tangent_tol_deg; /* how near a fill, Coons patch or corner blend fits its cross-tangents; <= 0 for 0.01 */
    double fit_tolerance;   /* how near, in mm, a rational edge's non-rational fit stays to it, for a fill, Coons
                             * patch, corner blend or skin; <= 0 for the larger of the brep's tolerance and 1e-7 of
                             * its diagonal */
    sm_label label;         /* for every face the surface makes; present 0 for none */
} sm_rebuild_surface;

/* What a rebuild did, filled in whether it succeeds or not. */
typedef struct sm_rebuild_report {
    int32_t solid;            /* the result is a closed manifold solid */
    uint32_t stitched;        /* edges the stitch glued */
    uint32_t open_edges;      /* edges it left open */
    double vertex_gap;        /* the widest gap between vertices it glued */
    double edge_gap;          /* the widest gap between edges it glued */
    double open_vertex_gap;   /* the narrowest gap between vertices it left apart, -1 for none */
    double open_edge_gap;     /* the narrowest gap between edges it left apart, -1 for none */
    uint32_t faces_made;      /* faces the surfaces made */
    int64_t failed_surface;   /* the surface the kernel couldn't make, -1 for none */
    sm_status surface_status; /* what it said */
    double surface_seconds;   /* making the surfaces */
    double stitch_seconds;    /* taking the faces out, adding the surfaces and stitching */
} sm_rebuild_report;

/* Takes the faces at the given indices out of a copy of brep, makes each
 * surface from the rim's edges and the kept faces beyond them, adds those
 * as faces labelled with their labels (a face a surface, split only where
 * it isn't C0), and stitches the lot into a solid (SmStitch::DoStitching at
 * stitch_tolerance, splitting an edge where another's vertex lies on it, as
 * a fill's patches' do on the rim; or with stitch_tolerance <= 0
 * SmStitch::StitchIntoSolid, which widens its tolerance up to 1/200 of the
 * part's diagonal and glues what it finds there, reporting only the widest
 * gap it glued), then orients it. The surfaces are made before any face
 * goes, from the kept faces' edgeuses. *out is the new solid; brep is
 * unchanged. On any failure *out is NULL and no open body is returned:
 * SM_SHIM_ERR_LEFT_OPEN or SM_SHIM_ERR_NOT_SOLID with the report saying
 * what the stitch did, SM_SHIM_ERR_SURFACE_FAILED with the surface and its
 * status. report may not be NULL. */
sm_status sm_brep_rebuild_faces(const sm_brep *brep, const uint32_t *faces, size_t face_count,
                                const sm_rebuild_surface *surfaces, size_t surface_count, double stitch_tolerance,
                                sm_brep **out, sm_rebuild_report *report);

#define SM_EDGE_CURVATURE_SAMPLES 12

/* The normal curvature across an edge on each of its faces at
 * SM_EDGE_CURVATURE_SAMPLES points along it (the edge_creases points), in
 * 1/length: each face's curvature in the direction square to the edge,
 * signed against that face's outward normal (positive bending away from
 * it, as a convex round does). valid 0 for an edge without two faces, a
 * seam, or where the kernel can't evaluate. */
typedef struct sm_edge_curvature {
    int32_t valid;
    double a[SM_EDGE_CURVATURE_SAMPLES];
    double b[SM_EDGE_CURVATURE_SAMPLES];
} sm_edge_curvature;

sm_status sm_brep_edge_curvatures(const sm_brep *brep, sm_edge_curvature *out, size_t capacity, size_t *count);

/* How SmSurface::CreateExtendedSurface carries a surface on. */
enum {
    SM_EXTEND_LINEAR = 0,  /* SM_CT_G1: straight along its tangent */
    SM_EXTEND_SMOOTH = 1,  /* SM_CT_G1_G2: with its second derivative */
    SM_EXTEND_NATURAL = 2, /* SM_CT_CINFINITY: as the surface runs on */
    SM_EXTEND_RULED = 3    /* SM_CT_G1R */
};
/* Which of its sides. */
enum { SM_SIDE_UMIN = 0, SM_SIDE_VMIN = 1, SM_SIDE_UMAX = 2, SM_SIDE_VMAX = 3, SM_SIDE_ALL = 4 };

/* A copy of the surface carried on past one side, or all four, by
 * distance (SmSurface::CreateExtendedSurface). */
sm_status sm_surface_extend(const sm_surface *surface, int32_t side, double distance, int32_t shape,
                            sm_surface **out);

#ifdef __cplusplus
}
#endif

#endif /* SMLIB_SHIM_REBUILD_H */
