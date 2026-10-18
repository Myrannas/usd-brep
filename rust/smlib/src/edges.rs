//! Picking B-rep edges by where they are, never by index.
//!
//! Edge indices change with every operation, so the fillet and chamfer calls
//! take an [`EdgeQuery`] instead: a geometric test that each edge's
//! [`EdgeGeometry`] either passes or fails.

use std::ops::Not;

use smlib_sys as sys;

use crate::label::{Label, LabelMatch};

/// How many points [`EdgeGeometry::samples`] holds.
pub const EDGE_SAMPLES: usize = sys::SM_EDGE_SAMPLES as usize;

/// Where a B-rep edge is: what an [`EdgeQuery`] tests.
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeGeometry {
    /// The edge's position in the brep's edge list, which is also its index
    /// in [`Mesh::edges`](crate::Mesh::edges). It holds until the brep next
    /// changes shape.
    pub index: usize,
    /// Points on the edge, evenly spaced in its curve parameter from one end
    /// to the other.
    pub samples: [[f64; 3]; EDGE_SAMPLES],
    /// The outward unit normals of the two faces the edge joins, at its
    /// midpoint. A seam joins a face to itself, so both come from one face.
    /// `None` when the edge has fewer than two faces or the kernel couldn't
    /// find the normal.
    pub normals: [Option<[f64; 3]>; 2],
    /// How many faces use the edge: 2 on a manifold solid.
    pub use_count: u32,
    /// The two faces the edge joins, as indices into
    /// [`Brep::faces`](crate::Brep::faces). A seam joins a face to itself.
    pub faces: [Option<u32>; 2],
    /// The edge's own label.
    pub label: Option<Label>,
    /// The labels of the two faces, in the order of `faces`.
    pub face_labels: [Option<Label>; 2],
}

impl EdgeGeometry {
    pub(crate) fn from_raw(index: usize, raw: &sys::sm_edge_geometry) -> EdgeGeometry {
        let normal = |n: [f64; 3]| (dot(n, n) > 0.5).then_some(n);
        let has = |face: i32| face >= 0;
        EdgeGeometry {
            index,
            faces: [raw.face_a, raw.face_b].map(|face| u32::try_from(face).ok()),
            label: Label::from_raw(&raw.label),
            face_labels: [Label::from_raw(&raw.label_a), Label::from_raw(&raw.label_b)],
            samples: raw.points,
            normals: [
                if has(raw.face_a) {
                    normal(raw.normal_a)
                } else {
                    None
                },
                if has(raw.face_b) {
                    normal(raw.normal_b)
                } else {
                    None
                },
            ],
            use_count: raw.use_count,
        }
    }

    /// One end of the edge.
    pub fn start(&self) -> [f64; 3] {
        self.samples[0]
    }

    /// The other end.
    pub fn end(&self) -> [f64; 3] {
        self.samples[EDGE_SAMPLES - 1]
    }

    /// The point halfway along the edge's curve parameter.
    pub fn midpoint(&self) -> [f64; 3] {
        self.samples[EDGE_SAMPLES / 2]
    }

    /// The angle between the two faces' outward normals at the midpoint, in
    /// degrees: 0 where the faces meet smoothly (and on a seam), 90 on the
    /// edge of a box. `None` without both normals.
    pub fn bend_deg(&self) -> Option<f64> {
        let [Some(a), Some(b)] = self.normals else {
            return None;
        };
        Some(angle_deg(a, b))
    }
}

/// How close counts as on, for an [`EdgeQuery`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerance {
    /// Distance, in model units.
    pub length: f64,
    /// Angle between directions, in degrees.
    pub angle_deg: f64,
}

impl Default for Tolerance {
    /// 1e-3 model units (a micron in millimetres) and 1°.
    fn default() -> Self {
        Tolerance {
            length: 1e-3,
            angle_deg: 1.0,
        }
    }
}

/// A geometric test that picks B-rep edges.
///
/// Directions (`normal`, `axis`, `a`, `b`) needn't be unit length but must
/// be non-zero. Combine queries with [`and`](EdgeQuery::and),
/// [`or`](EdgeQuery::or) and `!`.
///
/// ```
/// use smlib::EdgeQuery;
///
/// // The four edges around the top of a box whose top face is at z = 20.
/// let top = EdgeQuery::MidpointOnPlane { point: [0.0, 0.0, 20.0], normal: [0.0, 0.0, 1.0] };
/// // The rim of a hole of radius 10 through that face at the origin.
/// let rim = EdgeQuery::OnCircle { center: [0.0, 0.0, 20.0], axis: [0.0, 0.0, 1.0], radius: 10.0 };
/// // The top edges that aren't the rim.
/// let outer = top.and(!rim);
/// ```
#[derive(Clone, Debug, PartialEq)]
pub enum EdgeQuery {
    /// Every edge.
    All,
    /// Edges whose midpoint lies on the plane through `point` with normal
    /// `normal`.
    MidpointOnPlane {
        /// A point on the plane.
        point: [f64; 3],
        /// The plane's normal.
        normal: [f64; 3],
    },
    /// Edges lying in the plane: every sample point is on it.
    InPlane {
        /// A point on the plane.
        point: [f64; 3],
        /// The plane's normal.
        normal: [f64; 3],
    },
    /// Edges on a circle, or an arc of it: every sample point is `radius`
    /// from `center` in the plane through `center` square to `axis`.
    OnCircle {
        /// The circle's centre.
        center: [f64; 3],
        /// The circle's axis.
        axis: [f64; 3],
        /// The circle's radius.
        radius: f64,
    },
    /// Edges between two faces whose outward normals at the edge's midpoint
    /// point along `a` and `b`, either way round.
    BetweenFaces {
        /// One face's outward normal.
        a: [f64; 3],
        /// The other's.
        b: [f64; 3],
    },
    /// Edges where the faces' outward normals differ by at least
    /// `min_angle_deg` at the midpoint: creases, not seams or tangent joins.
    Sharp {
        /// The smallest angle between the normals, in degrees.
        min_angle_deg: f64,
    },
    /// Straight edges running along `direction`, either way.
    ParallelTo {
        /// The direction.
        direction: [f64; 3],
    },
    /// Edges between a face whose [`Label`] passes `a` and one whose label
    /// passes `b`, either way round.
    Between {
        /// One face's label.
        a: LabelMatch,
        /// The other's.
        b: LabelMatch,
    },
    /// The edges at these positions in [`Brep::edges`](crate::Brep::edges),
    /// for edges a caller has already found there. Positions hold only until
    /// the brep next changes shape.
    Indices(Vec<usize>),
    /// Edges every query picks.
    And(Vec<EdgeQuery>),
    /// Edges any query picks.
    Or(Vec<EdgeQuery>),
    /// Edges the query doesn't pick.
    Not(Box<EdgeQuery>),
}

impl EdgeQuery {
    /// Edges both `self` and `other` pick.
    pub fn and(self, other: EdgeQuery) -> EdgeQuery {
        match self {
            EdgeQuery::And(mut all) => {
                all.push(other);
                EdgeQuery::And(all)
            }
            first => EdgeQuery::And(vec![first, other]),
        }
    }

    /// Edges either `self` or `other` picks.
    pub fn or(self, other: EdgeQuery) -> EdgeQuery {
        match self {
            EdgeQuery::Or(mut any) => {
                any.push(other);
                EdgeQuery::Or(any)
            }
            first => EdgeQuery::Or(vec![first, other]),
        }
    }

    /// Whether the query picks this edge.
    pub fn matches(&self, edge: &EdgeGeometry, tolerance: &Tolerance) -> bool {
        let near = |distance: f64| distance.abs() <= tolerance.length;
        let along = |n: [f64; 3], d: [f64; 3]| angle_deg(n, d) <= tolerance.angle_deg;
        match self {
            EdgeQuery::All => true,
            EdgeQuery::MidpointOnPlane { point, normal } => {
                near(plane_distance(edge.midpoint(), *point, *normal))
            }
            EdgeQuery::InPlane { point, normal } => edge
                .samples
                .iter()
                .all(|&p| near(plane_distance(p, *point, *normal))),
            EdgeQuery::OnCircle {
                center,
                axis,
                radius,
            } => edge.samples.iter().all(|&p| {
                let offset = sub(p, *center);
                let axial = dot(offset, unit(*axis));
                let radial = dot(offset, offset) - axial * axial;
                near(axial) && near(radial.max(0.0).sqrt() - radius)
            }),
            EdgeQuery::BetweenFaces { a, b } => match edge.normals {
                [Some(n), Some(m)] => {
                    (along(n, *a) && along(m, *b)) || (along(n, *b) && along(m, *a))
                }
                _ => false,
            },
            EdgeQuery::Sharp { min_angle_deg } => {
                edge.bend_deg().is_some_and(|bend| bend >= *min_angle_deg)
            }
            EdgeQuery::ParallelTo { direction } => {
                let chord = sub(edge.end(), edge.start());
                let length = dot(chord, chord).sqrt();
                let straight = edge.samples.iter().all(|&p| {
                    let offset = sub(p, edge.start());
                    let along = dot(offset, chord) / length;
                    near((dot(offset, offset) - along * along).max(0.0).sqrt())
                });
                length > tolerance.length
                    && straight
                    && angle_deg(chord, *direction).min(180.0 - angle_deg(chord, *direction))
                        <= tolerance.angle_deg
            }
            EdgeQuery::Between { a, b } => match &edge.face_labels {
                [Some(m), Some(n)] => {
                    (a.matches(m) && b.matches(n)) || (a.matches(n) && b.matches(m))
                }
                _ => false,
            },
            EdgeQuery::Indices(indices) => indices.contains(&edge.index),
            EdgeQuery::And(all) => all.iter().all(|q| q.matches(edge, tolerance)),
            EdgeQuery::Or(any) => any.iter().any(|q| q.matches(edge, tolerance)),
            EdgeQuery::Not(query) => !query.matches(edge, tolerance),
        }
    }

    /// Why the query can't be evaluated, if it can't.
    pub(crate) fn problem(&self) -> Option<&'static str> {
        let finite = |v: &[f64]| v.iter().all(|c| c.is_finite());
        let direction = |d: [f64; 3]| finite(&d) && dot(d, d) > 0.0;
        match self {
            EdgeQuery::All => None,
            EdgeQuery::MidpointOnPlane { point, normal } | EdgeQuery::InPlane { point, normal } => {
                (!finite(point) || !direction(*normal))
                    .then_some("plane needs a finite point and a non-zero normal")
            }
            EdgeQuery::OnCircle {
                center,
                axis,
                radius,
            } => (!finite(center) || !direction(*axis) || !(radius.is_finite() && *radius > 0.0))
                .then_some("circle needs a finite centre, a non-zero axis and a positive radius"),
            EdgeQuery::BetweenFaces { a, b } => (!direction(*a) || !direction(*b))
                .then_some("face normals must be finite and non-zero"),
            EdgeQuery::Sharp { min_angle_deg } => {
                (!min_angle_deg.is_finite()).then_some("sharp-edge angle must be finite")
            }
            EdgeQuery::ParallelTo { direction: along } => {
                (!direction(*along)).then_some("direction must be finite and non-zero")
            }
            EdgeQuery::Between { .. } | EdgeQuery::Indices(_) => None,
            EdgeQuery::And(queries) | EdgeQuery::Or(queries) => {
                queries.iter().find_map(EdgeQuery::problem)
            }
            EdgeQuery::Not(query) => query.problem(),
        }
    }
}

impl Not for EdgeQuery {
    type Output = EdgeQuery;

    fn not(self) -> EdgeQuery {
        match self {
            EdgeQuery::Not(query) => *query,
            query => EdgeQuery::Not(Box::new(query)),
        }
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let length = dot(v, v).sqrt();
    [v[0] / length, v[1] / length, v[2] / length]
}

fn plane_distance(p: [f64; 3], point: [f64; 3], normal: [f64; 3]) -> f64 {
    dot(sub(p, point), unit(normal))
}

fn angle_deg(a: [f64; 3], b: [f64; 3]) -> f64 {
    dot(unit(a), unit(b)).clamp(-1.0, 1.0).acos().to_degrees()
}
