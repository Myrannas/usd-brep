//! What [`Brep::imprint`](crate::Brep::imprint) gives back.

use smlib_sys as sys;

use crate::brep::Brep;
use crate::label::Label;

/// A brep with a tool imprinted on it, by
/// [`Brep::imprint`](crate::Brep::imprint) or
/// [`Brep::imprint_classified`](crate::Brep::imprint_classified).
#[derive(Debug)]
pub struct Imprint {
    /// The brep, its faces and edges split where the tool's boundary meets
    /// them.
    pub brep: Brep,
    /// Every edge of `brep` that lies along the tool's boundary, new or
    /// already there, in edge order.
    pub edges: Vec<ImprintEdge>,
    /// Whether each face of `brep` lies inside the tool, in face order, for
    /// [`Brep::imprint_classified`](crate::Brep::imprint_classified); empty
    /// for [`Brep::imprint`](crate::Brep::imprint).
    pub inside: Vec<bool>,
}

/// An edge of an [`Imprint`] that lies along the tool's boundary.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ImprintEdge {
    /// Its index in the imprinted brep's [`edges`](crate::Brep::edges).
    pub edge: usize,
    /// Whether the imprint made it, rather than finding it, or the edge it
    /// was split from, already there along the tool.
    pub new: bool,
    /// The labels of the tool's faces it lies on: one where it lies inside a
    /// tool face, two where it runs along one of the tool's edges, and none
    /// where those faces have no label. Label the tool's faces before
    /// imprinting to tell which made each edge.
    pub tool_faces: Vec<Label>,
}

impl ImprintEdge {
    pub(crate) fn from_raw(raw: &sys::sm_imprint_edge) -> ImprintEdge {
        ImprintEdge {
            edge: raw.edge,
            new: raw.is_new != 0,
            tool_faces: raw.tool_faces.iter().filter_map(Label::from_raw).collect(),
        }
    }
}

/// A curve in a face's (u, v), for [`Brep::imprint_uv`](crate::Brep::imprint_uv):
/// a clamped B-spline of `degree`, each control point's u, v and weight
/// (not multiplied through), and `points.len() + degree + 1` knots.
#[derive(Clone, Debug, PartialEq)]
pub struct UvCurve {
    /// The face's index in [`faces`](crate::Brep::faces), before the call.
    pub face: usize,
    /// Its degree, at least 1.
    pub degree: u32,
    /// Its knots, clamped.
    pub knots: Vec<f64>,
    /// Its control points: u, v and weight.
    pub points: Vec<[f64; 3]>,
}
