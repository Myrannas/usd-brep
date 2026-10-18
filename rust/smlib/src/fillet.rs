//! Fillet and chamfer settings for [`Brep::fillet`](crate::Brep::fillet).

use smlib_sys as sys;

use crate::edges::{EdgeQuery, Tolerance};

/// The shape a fillet leaves across the edge (`SmFilletSurfaceGeneratorType`).
///
/// Each runs between the two lines where a ball of the fillet's radius,
/// rolled along the edge, touches the faces. On a 90° edge those lines are
/// set back by the radius on each face.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum CrossSection {
    /// A circular arc: the ordinary round.
    #[default]
    Circular,
    /// A straight line: a chamfer set back by the radius on each face of a
    /// 90° edge. For a chamfer given by its width, see
    /// [`Brep::chamfer`](crate::Brep::chamfer).
    Linear,
    /// A blend curve meeting the faces with [`FilletOptions::continuity`].
    Blend,
}

impl CrossSection {
    pub(crate) fn to_raw(self) -> i32 {
        (match self {
            CrossSection::Circular => sys::SM_FILLET_CIRCULAR,
            CrossSection::Linear => sys::SM_FILLET_LINEAR,
            CrossSection::Blend => sys::SM_FILLET_BLEND,
        }) as i32
    }
}

/// How smoothly a [`CrossSection::Blend`] meets the faces.
///
/// Corners where three filleted edges meet are only G1 whatever this says.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Continuity {
    /// Tangent: no crease, but reflections kink.
    #[default]
    G1 = 1,
    /// Curvature: smooth reflections.
    G2 = 2,
    /// Third order: smooth highlight lines.
    G3 = 3,
}

/// One group of edges to fillet: the edges an [`EdgeQuery`] picks, with a
/// radius and a cross-section.
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeBlend {
    /// Which edges.
    pub edges: EdgeQuery,
    /// The rolling ball's radius, in model units.
    pub radius: f64,
    /// The shape across the edge.
    pub section: CrossSection,
}

impl EdgeBlend {
    /// A circular fillet of `radius` on the picked edges.
    pub fn round(edges: EdgeQuery, radius: f64) -> EdgeBlend {
        EdgeBlend {
            edges,
            radius,
            section: CrossSection::Circular,
        }
    }

    /// A flat cut on the picked edges, set back by `radius` on each face of a
    /// 90° edge.
    pub fn linear(edges: EdgeQuery, radius: f64) -> EdgeBlend {
        EdgeBlend {
            edges,
            radius,
            section: CrossSection::Linear,
        }
    }

    /// A blend curve on the picked edges, with the continuity in
    /// [`FilletOptions`].
    pub fn blend(edges: EdgeQuery, radius: f64) -> EdgeBlend {
        EdgeBlend {
            edges,
            radius,
            section: CrossSection::Blend,
        }
    }
}

/// Settings shared by every edge in one [`Brep::fillet`](crate::Brep::fillet)
/// call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilletOptions {
    /// How smoothly blend sections meet the faces.
    pub continuity: Continuity,
    /// How full blend sections are: 1 is the kernel's default, larger is
    /// fuller.
    pub thumbweight: f64,
    /// How the edge queries judge "on".
    pub tolerance: Tolerance,
}

impl Default for FilletOptions {
    fn default() -> Self {
        FilletOptions {
            continuity: Continuity::G1,
            thumbweight: 1.0,
            tolerance: Tolerance::default(),
        }
    }
}
