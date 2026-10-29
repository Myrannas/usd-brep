//! Rebuilding faces: faces taken out of a solid, and surfaces built in the
//! gaps from the rim's edges and the kept faces beyond them, sewn back into
//! a solid in one kernel call, so an open body never reaches Rust
//! ([`Brep::rebuild_faces`]). Also what checks how they meet: the normal
//! curvature across every edge ([`Brep::edge_curvatures`], the curvature
//! twin of [`Brep::edge_creases`]) and surfaces carried on past their sides
//! ([`FreeSurface::extended`]).

use std::fmt;
use std::ptr;

use smlib_sys as sys;
use sys::rebuild as rb;

use crate::brep::Brep;
use crate::error::{Error, Status};
use crate::freeform::FreeSurface;
use crate::label::{Label, NO_LABEL};

/// How a side meets the face beyond it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Meets {
    /// Along the edge, at any angle (G0).
    Position,
    /// With no crease (G1).
    Tangent,
    /// As curved as the face beyond (G2), where the tool can.
    Curvature,
}

impl Meets {
    fn raw(self) -> i32 {
        match self {
            Meets::Position => rb::SM_MEETS_POSITION,
            Meets::Tangent => rb::SM_MEETS_TANGENT,
            Meets::Curvature => rb::SM_MEETS_CURVATURE,
        }
    }
}

/// A side of a surface: a rim edge, by its index in the brep's edge list
/// before any face is taken out, and how the surface meets the kept face
/// beyond it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Side {
    /// The rim edge, between a face taken out and a kept one.
    pub edge: usize,
    /// How the surface meets the kept face.
    pub meets: Meets,
}

/// The kernel call that makes a surface.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tool {
    /// `SmBSplineSurface::CreateNSidedPatch` through 3 to 8 sides in order
    /// round the loop, a Tangent side's cross-derivatives sampled from the
    /// kept face (`SmTangentField`), fitted with matching parameters: a
    /// patch a side, meeting in the middle. Sides next to one another on one
    /// kept face, meeting it alike, are one side.
    Fill,
    /// `SmBSplineSurface::CreateCoonsPatch`, bicubic, through 4 sides in
    /// order, joined as a fill's: Tangent sides' cross-derivatives as a
    /// fill's, a Position side's a cubic between the sides beside it, taking
    /// the corners' twists from their fields so the corners agree.
    Coons,
    /// `SmBrep::CreateCornerBlend` from the rim's edgeuses, 3 or 4 sides in
    /// any order: a bicubic Coons patch, a three-sided one with a point for
    /// its fourth side.
    Corner,
    /// `SmBrep::CreateEdgeEdgeBlend` between the first two sides: 0 ruled,
    /// 1 G1, 2 G2, 3 G3, both ends alike.
    EdgeBlend {
        /// The blend's order.
        order: u32,
    },
    /// `SmBSplineSurface::CreateSkinnedSurface` through the first two
    /// sides, a Tangent end taking its cross-derivatives from a derivative
    /// surface of the kept face (G1 at most).
    Skin,
}

/// A surface to build in a gap.
#[derive(Clone, Debug, PartialEq)]
pub struct RebuildSurface {
    /// How it's made.
    pub tool: Tool,
    /// Its sides.
    pub sides: Vec<Side>,
    /// How near, in degrees, a fill, Coons patch or corner blend fits the
    /// kept faces' tangent planes; 0 for the shim's 0.01.
    pub tangent_tol_deg: f64,
    /// How near, in mm, a rational edge's non-rational fit stays to it, for
    /// a fill, Coons patch or corner blend; 0 for the brep's tolerance.
    pub fit_tolerance: f64,
    /// The label of every face it makes.
    pub label: Option<Label>,
}

/// What a rebuild did, whether it built or not.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RebuildReport {
    /// The result is a closed manifold solid.
    pub solid: bool,
    /// Edges the stitch glued.
    pub stitched: u32,
    /// Edges it left open.
    pub open_edges: u32,
    /// The widest gap between vertices it glued.
    pub vertex_gap: f64,
    /// The widest gap between edges it glued.
    pub edge_gap: f64,
    /// The narrowest gap between vertices it left apart.
    pub open_vertex_gap: Option<f64>,
    /// The narrowest gap between edges it left apart.
    pub open_edge_gap: Option<f64>,
    /// Faces the surfaces made.
    pub faces_made: u32,
    /// The surface the kernel couldn't make, and what it said.
    pub failed_surface: Option<(usize, Status)>,
    /// Seconds making the surfaces.
    pub surface_seconds: f64,
    /// Seconds taking the faces out, adding the surfaces and stitching.
    pub stitch_seconds: f64,
}

impl RebuildReport {
    fn from_raw(raw: &rb::sm_rebuild_report) -> RebuildReport {
        let gap = |g: f64| (g >= 0.0).then_some(g);
        RebuildReport {
            solid: raw.solid != 0,
            stitched: raw.stitched,
            open_edges: raw.open_edges,
            vertex_gap: raw.vertex_gap,
            edge_gap: raw.edge_gap,
            open_vertex_gap: gap(raw.open_vertex_gap),
            open_edge_gap: gap(raw.open_edge_gap),
            faces_made: raw.faces_made,
            failed_surface: usize::try_from(raw.failed_surface)
                .ok()
                .map(|s| (s, Status::from_raw(raw.surface_status))),
            surface_seconds: raw.surface_seconds,
            stitch_seconds: raw.stitch_seconds,
        }
    }
}

/// A rebuilt solid.
pub struct Rebuilt {
    /// The solid, with the faces taken out and the surfaces sewn in.
    pub brep: Brep,
    /// What the rebuild did.
    pub report: RebuildReport,
}

/// Why a rebuild didn't build, with what it did get to.
#[derive(Clone, Debug, PartialEq)]
pub struct RebuildError {
    /// The failure.
    pub error: Error,
    /// What the rebuild did, boxed so a failure stays small.
    pub report: Box<RebuildReport>,
}

impl fmt::Display for RebuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.error)
    }
}

impl std::error::Error for RebuildError {}

/// The normal curvature across an edge on each of its two faces at the
/// [`Brep::edge_creases`] points, in 1/length, each against its face's
/// outward normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeCurvature {
    /// The first face's.
    pub a: [f64; rb::SM_EDGE_CURVATURE_SAMPLES],
    /// The second face's.
    pub b: [f64; rb::SM_EDGE_CURVATURE_SAMPLES],
}

impl EdgeCurvature {
    /// The largest difference between the two sides' curvatures, in
    /// 1/length.
    pub fn worst_difference(&self) -> f64 {
        self.a.iter().zip(&self.b).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)
    }

    /// The largest difference relative to the more curved side there, 0
    /// where both are flatter than `flat` (1/length).
    pub fn worst_ratio(&self, flat: f64) -> f64 {
        self.a
            .iter()
            .zip(&self.b)
            .map(|(a, b)| {
                let most = a.abs().max(b.abs());
                if most < flat { 0.0 } else { (a - b).abs() / most }
            })
            .fold(0.0, f64::max)
    }
}

/// How [`FreeSurface::extended`] carries a surface on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExtendShape {
    /// Straight along its tangent (`SM_CT_G1`).
    Linear,
    /// With its second derivative (`SM_CT_G1_G2`).
    Smooth,
    /// As the surface runs on (`SM_CT_CINFINITY`).
    Natural,
    /// `SM_CT_G1R`.
    Ruled,
}

/// Which side of a surface's (u, v) domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceSide {
    /// u at its least.
    UMin,
    /// v at its least.
    VMin,
    /// u at its most.
    UMax,
    /// v at its most.
    VMax,
    /// All four.
    All,
}

impl Brep {
    /// Takes the faces at `faces` out of a copy of the solid, builds each
    /// surface from its sides (rim edges, by index in this brep) and the
    /// kept faces beyond them, and sews them in: `SmStitch::DoStitching` at
    /// `stitch_tolerance`, its report saying the narrowest gap it left open,
    /// or with 0 `SmStitch::StitchIntoSolid`, which widens its tolerance up
    /// to 1/200 of the part's diagonal and glues whatever it finds there,
    /// saying only the widest gap it glued. The new solid, or why there
    /// isn't one; this one is unchanged.
    pub fn rebuild_faces(
        &self,
        faces: &[usize],
        surfaces: &[RebuildSurface],
        stitch_tolerance: f64,
    ) -> Result<Rebuilt, RebuildError> {
        const OP: &str = "sm_brep_rebuild_faces";
        let refuse = |reason| RebuildError {
            error: Error::InvalidArgument { operation: OP, reason },
            report: Box::default(),
        };
        let faces: Vec<u32> = faces
            .iter()
            .map(|&f| u32::try_from(f).map_err(|_| refuse("a face index past u32")))
            .collect::<Result<_, _>>()?;
        let sides: Vec<Vec<rb::sm_rebuild_side>> = surfaces
            .iter()
            .map(|s| {
                s.sides
                    .iter()
                    .map(|side| {
                        Ok(rb::sm_rebuild_side {
                            edge: u32::try_from(side.edge).map_err(|_| refuse("an edge index past u32"))?,
                            meets: side.meets.raw(),
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<_, _>>()?;
        let mut raw = Vec::with_capacity(surfaces.len());
        for (surface, sides) in surfaces.iter().zip(&sides) {
            let (tool, order) = match surface.tool {
                Tool::Fill => (rb::SM_REBUILD_FILL, 0),
                Tool::Coons => (rb::SM_REBUILD_COONS, 0),
                Tool::Corner => (rb::SM_REBUILD_CORNER, 0),
                Tool::EdgeBlend { order } => (rb::SM_REBUILD_EDGE_BLEND, order as i32),
                Tool::Skin => (rb::SM_REBUILD_SKIN, 0),
            };
            let label = match &surface.label {
                Some(label) => label.to_raw(OP).map_err(|error| RebuildError {
                    error,
                    report: Box::default(),
                })?,
                None => NO_LABEL,
            };
            raw.push(rb::sm_rebuild_surface {
                tool,
                sides: sides.as_ptr(),
                side_count: sides.len(),
                order,
                tangent_tol_deg: surface.tangent_tol_deg,
                fit_tolerance: surface.fit_tolerance,
                label,
            });
        }
        let mut out = ptr::null_mut();
        let mut report = rb::sm_rebuild_report::default();
        // SAFETY: self is live on the kernel thread; faces, raw and the
        // sides each points into outlive the call.
        let status = unsafe {
            rb::sm_brep_rebuild_faces(
                self.as_raw(),
                faces.as_ptr(),
                faces.len(),
                raw.as_ptr(),
                raw.len(),
                stitch_tolerance,
                &mut out,
                &mut report,
            )
        };
        let report = RebuildReport::from_raw(&report);
        if let Err(error) = Error::check(OP, status) {
            return Err(RebuildError { error, report: Box::new(report) });
        }
        match Brep::from_raw(OP, out) {
            Ok(brep) => Ok(Rebuilt { brep, report }),
            Err(error) => Err(RebuildError { error, report: Box::new(report) }),
        }
    }

    /// The normal curvature across each edge on its two faces, in the
    /// order of [`edges`](Brep::edges), at the twelve points
    /// [`edge_creases`](Brep::edge_creases) reads; `None` for an edge
    /// without two faces, a seam, or where the kernel can't evaluate.
    pub fn edge_curvatures(&self) -> Result<Vec<Option<EdgeCurvature>>, Error> {
        let mut raw = vec![rb::sm_edge_curvature::default(); self.edge_count()];
        let mut count = 0;
        // SAFETY: self is live on the kernel thread; raw has room for
        // raw.len() edges.
        let status = unsafe {
            rb::sm_brep_edge_curvatures(self.as_raw(), raw.as_mut_ptr(), raw.len(), &mut count)
        };
        Error::check("sm_brep_edge_curvatures", status)?;
        raw.truncate(count);
        Ok(raw
            .into_iter()
            .map(|e| (e.valid != 0).then_some(EdgeCurvature { a: e.a, b: e.b }))
            .collect())
    }
}

impl FreeSurface {
    /// A copy carried on past `side` by `distance`, shaped as `shape`
    /// (`SmSurface::CreateExtendedSurface`).
    pub fn extended(&self, side: SurfaceSide, distance: f64, shape: ExtendShape) -> Result<FreeSurface, Error> {
        const OP: &str = "SmSurface::CreateExtendedSurface";
        let side = match side {
            SurfaceSide::UMin => rb::SM_SIDE_UMIN,
            SurfaceSide::VMin => rb::SM_SIDE_VMIN,
            SurfaceSide::UMax => rb::SM_SIDE_UMAX,
            SurfaceSide::VMax => rb::SM_SIDE_VMAX,
            SurfaceSide::All => rb::SM_SIDE_ALL,
        };
        let shape = match shape {
            ExtendShape::Linear => rb::SM_EXTEND_LINEAR,
            ExtendShape::Smooth => rb::SM_EXTEND_SMOOTH,
            ExtendShape::Natural => rb::SM_EXTEND_NATURAL,
            ExtendShape::Ruled => rb::SM_EXTEND_RULED,
        };
        let mut out = ptr::null_mut();
        // SAFETY: self is live on the kernel thread.
        let status = unsafe { rb::sm_surface_extend(self.as_raw(), side, distance, shape, &mut out) };
        Error::check(OP, status)?;
        FreeSurface::wrap(OP, out)
    }
}
