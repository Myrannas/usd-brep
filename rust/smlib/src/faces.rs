//! Faces as places to put things: their geometry, frames to place solids in,
//! prisms swept from them, and queries that pick them.

use std::ops::Not;

use smlib_sys as sys;

use crate::brep::Brep;
use crate::edges::Tolerance;
use crate::error::Error;
use crate::label::{Label, LabelMatch, NO_LABEL};
use crate::mesh::{Face, Surface};

/// Where a B-rep face is, for finding faces and placing solids on them. For
/// its area and centroid, which take longer to find, see
/// [`Brep::face_centroid`](crate::Brep::face_centroid).
#[derive(Clone, Debug, PartialEq)]
pub struct FaceGeometry {
    /// The face's position in the brep's face list, which is also its index in
    /// [`Mesh::faces`](crate::Mesh::faces). It holds until the brep next
    /// changes shape.
    pub index: usize,
    /// The surface the face lies on.
    pub surface: Surface,
    /// The face's domain in its surface's (u, v): minimum corner.
    pub uv_min: [f64; 2],
    /// The face's domain in its surface's (u, v): maximum corner.
    pub uv_max: [f64; 2],
    /// The face's [`Label`], if it has one.
    pub label: Option<Label>,
    /// An axis-aligned box around the trimmed face, as its minimum and maximum
    /// corners. It follows the face's boundary, so each piece of a split face
    /// has its own, but it may be loose around curves.
    pub bounds: [[f64; 3]; 2],
    /// For a planar face, its unit normal pointing out of the solid; `None`
    /// for curved faces. Planar NURBS surfaces count as planar.
    pub normal: Option<[f64; 3]>,
}

impl FaceGeometry {
    pub(crate) fn from_raw(index: usize, raw: &sys::sm_face_geometry) -> FaceGeometry {
        let Face {
            surface,
            uv_min,
            uv_max,
            label,
        } = Face::from_raw(&raw.info);
        FaceGeometry {
            index,
            surface,
            uv_min,
            uv_max,
            label,
            bounds: [raw.bounds_min, raw.bounds_max],
            normal: (raw.planar != 0).then_some(raw.normal),
        }
    }

    /// The middle of [`bounds`](FaceGeometry::bounds).
    pub fn bounds_centre(&self) -> [f64; 3] {
        let [min, max] = self.bounds;
        [0, 1, 2].map(|i| (min[i] + max[i]) / 2.0)
    }
}

/// A face's area and area centroid, from
/// [`Brep::face_centroid`](crate::Brep::face_centroid).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FaceCentroid {
    /// The trimmed face's area.
    pub area: f64,
    /// The trimmed face's area centroid. On a planar face it lies in the
    /// face's plane, though not necessarily inside the face.
    pub centroid: [f64; 3],
}

/// A right-handed orthonormal frame to place a solid in with
/// [`Brep::place`](crate::Brep::place).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    /// Where the solid's origin goes.
    pub origin: [f64; 3],
    /// Where its x axis points: a unit vector.
    pub x_axis: [f64; 3],
    /// Where its y axis points: a unit vector perpendicular to `x_axis`.
    pub y_axis: [f64; 3],
}

impl Frame {
    /// A frame at `origin` whose x axis points along `x` and whose y axis is
    /// `y` made perpendicular to it. `None` when either is zero or not finite,
    /// or they're parallel.
    pub fn new(origin: [f64; 3], x: [f64; 3], y: [f64; 3]) -> Option<Frame> {
        let x_axis = normalised(x)?;
        let y_axis = normalised(sub(y, scale(x_axis, dot(y, x_axis))))?;
        origin.iter().all(|c| c.is_finite()).then_some(Frame {
            origin,
            x_axis,
            y_axis,
        })
    }

    /// The frame's z axis, `x_axis × y_axis`.
    pub fn z_axis(&self) -> [f64; 3] {
        cross(self.x_axis, self.y_axis)
    }

    /// Where a point given in the frame's coordinates lies.
    pub fn to_world(&self, [x, y, z]: [f64; 3]) -> [f64; 3] {
        let mut point = self.origin;
        for (axis, amount) in [(self.x_axis, x), (self.y_axis, y), (self.z_axis(), z)] {
            point = add(point, scale(axis, amount));
        }
        point
    }
}

/// What [`Brep::scale_axes_keeping`](crate::Brep::scale_axes_keeping) or
/// [`Brep::mirror`](crate::Brep::mirror) made of a brep's faces and edges.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ScaleReport {
    /// Analytic faces still of the same kind.
    pub faces_kept: u32,
    /// Cylinders stretched across their axis, now ellipses swept along a
    /// line.
    pub faces_extruded: u32,
    /// Analytic faces now NURBS.
    pub faces_nurbs: u32,
    /// Faces that were NURBS already.
    pub faces_were_nurbs: u32,
    /// Lines, and circles and ellipses still of the same kind.
    pub curves_kept: u32,
    /// Circles now ellipses.
    pub curves_stretched: u32,
    /// Analytic curves now NURBS.
    pub curves_nurbs: u32,
    /// Curves that were NURBS already.
    pub curves_were_nurbs: u32,
}

/// How two solids lie, from
/// [`Brep::relationship`](crate::Brep::relationship).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Relationship {
    /// Their surfaces are more than both breps' tolerances apart, and
    /// neither lies in the other's material. A solid in a cavity of the
    /// other is apart.
    Apart {
        /// Their surfaces' closest approach.
        gap: f64,
    },
    /// Their surfaces meet, and no sample of either's vertices and face
    /// points lies inside the other.
    Touching,
    /// Their surfaces meet, and a sample of one's points lies inside the
    /// other.
    Overlapping,
    /// The other lies in this one's material.
    Contains {
        /// The gap between their walls.
        gap: f64,
    },
    /// This lies in the other's material.
    Inside {
        /// The gap between their walls.
        gap: f64,
    },
}

/// A closed solid swept from a planar face by
/// [`Brep::sweep_face`](crate::Brep::sweep_face) (a prism) or
/// [`Brep::revolve_face`](crate::Brep::revolve_face).
#[derive(Debug)]
pub struct Sweep {
    /// The solid, apart from the brep it was swept from.
    pub prism: Brep,
    /// What each of the solid's faces came from, in its face order.
    pub faces: Vec<SweptFace>,
}

impl Sweep {
    /// Takes over a solid and the roles of its faces from the shim.
    ///
    /// # Safety
    ///
    /// `faces` holds `count` entries, or is null, and is freed here.
    pub(crate) unsafe fn from_raw(
        operation: &'static str,
        prism: *mut sys::sm_brep,
        faces: *mut sys::sm_swept_face,
        count: usize,
    ) -> Result<Sweep, Error> {
        // SAFETY: as the caller promises; we copy the entries, then free them
        // once.
        let roles = unsafe {
            let roles = if faces.is_null() {
                Vec::new()
            } else {
                std::slice::from_raw_parts(faces, count)
                    .iter()
                    .map(SweptFace::from_raw)
                    .collect()
            };
            sys::sm_swept_faces_free(faces);
            roles
        };
        Ok(Sweep {
            prism: Brep::from_raw(operation, prism)?,
            faces: roles,
        })
    }
}

/// What one face of a [`Sweep`] came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SweptFace {
    /// The copy of the swept face where the sweep starts, on the face itself.
    Start,
    /// The copy of the swept face the sweep distance (or angle) away.
    End,
    /// The face swept from one edge of the swept face.
    Side {
        /// That edge, as its index in the source brep's
        /// [`edges`](crate::Brep::edges), or `None` when it couldn't be told.
        source_edge: Option<usize>,
    },
}

impl SweptFace {
    pub(crate) fn from_raw(raw: &sys::sm_swept_face) -> SweptFace {
        match raw.role as sys::sm_swept_role {
            sys::SM_SWEPT_START => SweptFace::Start,
            sys::SM_SWEPT_END => SweptFace::End,
            _ => SweptFace::Side {
                source_edge: usize::try_from(raw.source_edge).ok(),
            },
        }
    }
}

/// A line, arc or circle in a plane, for
/// [`Brep::planar_faces`](crate::Brep::planar_faces). Points are in the
/// plane's own (x, y); angles are in degrees from x towards y.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlaneCurve {
    /// The segment from `start` to `end` (`SmApiCreateLineSegment`).
    Line {
        /// Where it starts.
        start: [f64; 2],
        /// Where it ends.
        end: [f64; 2],
    },
    /// An arc from `start_deg` by `sweep_deg`, counter-clockwise when
    /// positive (`SmApiCreateArc`).
    Arc {
        /// The arc's centre.
        centre: [f64; 2],
        /// Its radius.
        radius: f64,
        /// The angle it starts at.
        start_deg: f64,
        /// How far it turns: not 0, and at most a whole turn either way.
        sweep_deg: f64,
    },
    /// A whole circle, a loop on its own (`SmApiCreateCircle`).
    Circle {
        /// The circle's centre.
        centre: [f64; 2],
        /// Its radius.
        radius: f64,
    },
}

impl PlaneCurve {
    /// Where the curve starts; a circle starts and ends on its x side.
    pub fn start(&self) -> [f64; 2] {
        match *self {
            PlaneCurve::Line { start, .. } => start,
            PlaneCurve::Arc {
                centre,
                radius,
                start_deg,
                ..
            } => on_circle(centre, radius, start_deg),
            PlaneCurve::Circle { centre, radius } => on_circle(centre, radius, 0.0),
        }
    }

    /// Where the curve ends.
    pub fn end(&self) -> [f64; 2] {
        match *self {
            PlaneCurve::Line { end, .. } => end,
            PlaneCurve::Arc {
                centre,
                radius,
                start_deg,
                sweep_deg,
            } => on_circle(centre, radius, start_deg + sweep_deg),
            PlaneCurve::Circle { centre, radius } => on_circle(centre, radius, 0.0),
        }
    }

    /// The unit direction the curve runs in where it starts; a circle runs
    /// anticlockwise.
    pub(crate) fn start_tangent(&self) -> [f64; 2] {
        match *self {
            PlaneCurve::Line { start, end } => unit([end[0] - start[0], end[1] - start[1]]),
            PlaneCurve::Arc {
                start_deg,
                sweep_deg,
                ..
            } => along_circle(start_deg, sweep_deg),
            PlaneCurve::Circle { .. } => along_circle(0.0, 1.0),
        }
    }

    /// The unit direction the curve runs in where it ends.
    pub(crate) fn end_tangent(&self) -> [f64; 2] {
        match *self {
            PlaneCurve::Arc {
                start_deg,
                sweep_deg,
                ..
            } => along_circle(start_deg + sweep_deg, sweep_deg),
            _ => self.start_tangent(),
        }
    }

    /// Why the shim would refuse the curve, if it would.
    pub(crate) fn problem(&self) -> Option<&'static str> {
        let finite = |p: [f64; 2]| p.iter().all(|c| c.is_finite());
        let radius_ok = |r: f64| r.is_finite() && r > 0.0;
        match *self {
            PlaneCurve::Line { start, end } => (!(finite(start) && finite(end)) || start == end)
                .then_some("a line needs two different, finite ends"),
            PlaneCurve::Arc {
                centre,
                radius,
                start_deg,
                sweep_deg,
            } => (!(finite(centre)
                && radius_ok(radius)
                && start_deg.is_finite()
                && sweep_deg != 0.0
                && sweep_deg.abs() <= 360.0))
                .then_some(
                    "an arc needs a finite centre, a positive radius and a sweep of up to a turn",
                ),
            PlaneCurve::Circle { centre, radius } => (!(finite(centre) && radius_ok(radius)))
                .then_some("a circle needs a finite centre and a positive radius"),
        }
    }

    pub(crate) fn to_raw(
        self,
        label: &Option<Label>,
        operation: &'static str,
    ) -> Result<sys::sm_plane_curve, Error> {
        let mut raw = sys::sm_plane_curve {
            kind: 0,
            a: [0.0; 2],
            b: [0.0; 2],
            radius: 0.0,
            start_deg: 0.0,
            sweep_deg: 0.0,
            label: label
                .as_ref()
                .map_or(Ok(NO_LABEL), |l| l.to_raw(operation))?,
        };
        match self {
            PlaneCurve::Line { start, end } => {
                raw.kind = sys::SM_PLANE_LINE as i32;
                raw.a = start;
                raw.b = end;
            }
            PlaneCurve::Arc {
                centre,
                radius,
                start_deg,
                sweep_deg,
            } => {
                raw.kind = sys::SM_PLANE_ARC as i32;
                raw.a = centre;
                raw.radius = radius;
                raw.start_deg = start_deg;
                raw.sweep_deg = sweep_deg;
            }
            PlaneCurve::Circle { centre, radius } => {
                raw.kind = sys::SM_PLANE_CIRCLE as i32;
                raw.a = centre;
                raw.radius = radius;
            }
        }
        Ok(raw)
    }
}

/// One curve of a planar face's boundary, for
/// [`Brep::planar_faces`](crate::Brep::planar_faces), and the label its edge
/// gets.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryCurve {
    /// The curve.
    pub curve: PlaneCurve,
    /// The label for the curve's edge, if any.
    pub label: Option<Label>,
}

impl BoundaryCurve {
    /// A curve whose edge gets `label`.
    pub fn labelled(curve: PlaneCurve, label: Label) -> BoundaryCurve {
        BoundaryCurve {
            curve,
            label: Some(label),
        }
    }
}

fn unit([x, y]: [f64; 2]) -> [f64; 2] {
    let length = x.hypot(y);
    [x / length, y / length]
}

/// The direction a circle runs in at `angle_deg`, anticlockwise when `turn`
/// is positive.
fn along_circle(angle_deg: f64, turn: f64) -> [f64; 2] {
    let (sin, cos) = angle_deg.to_radians().sin_cos();
    let sign = turn.signum();
    [-sin * sign, cos * sign]
}

fn on_circle(centre: [f64; 2], radius: f64, angle_deg: f64) -> [f64; 2] {
    let (sin, cos) = angle_deg.to_radians().sin_cos();
    [centre[0] + radius * cos, centre[1] + radius * sin]
}

/// The kind of surface a face lies on, for [`FaceQuery::Kind`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SurfaceKind {
    /// A plane.
    Plane,
    /// A circular cylinder.
    Cylinder,
    /// A circular cone.
    Cone,
    /// A sphere.
    Sphere,
    /// A torus.
    Torus,
    /// A NURBS surface.
    Nurbs,
    /// Any other surface class.
    Other,
}

impl SurfaceKind {
    /// The kind of `surface`.
    pub fn of(surface: &Surface) -> SurfaceKind {
        match surface {
            Surface::Plane { .. } => SurfaceKind::Plane,
            Surface::Cylinder { .. } => SurfaceKind::Cylinder,
            Surface::Cone { .. } => SurfaceKind::Cone,
            Surface::Sphere { .. } => SurfaceKind::Sphere,
            Surface::Torus { .. } => SurfaceKind::Torus,
            Surface::Nurbs { .. } => SurfaceKind::Nurbs,
            Surface::Other { .. } => SurfaceKind::Other,
        }
    }
}

/// A test that picks faces by where they are or what they're labelled, for
/// [`Brep::find_faces`](crate::Brep::find_faces). Face indices change with
/// every operation; a query finds the faces again.
///
/// Directions needn't be unit length but must be non-zero. Combine queries
/// with [`and`](FaceQuery::and), [`or`](FaceQuery::or) and `!`.
///
/// ```
/// use smlib::FaceQuery;
///
/// // The top face of a box whose top is at z = 20.
/// let top = FaceQuery::OnPlane { point: [0.0, 0.0, 20.0], normal: [0.0, 0.0, 1.0] };
/// ```
#[derive(Clone, Debug, PartialEq)]
pub enum FaceQuery {
    /// Every face.
    All,
    /// Planar faces lying on the plane through `point` with normal `normal`,
    /// facing out along `normal`.
    OnPlane {
        /// A point on the plane.
        point: [f64; 3],
        /// The plane's normal, pointing out of the solid.
        normal: [f64; 3],
    },
    /// Planar faces whose outward normal is along `direction`.
    Facing {
        /// The outward direction.
        direction: [f64; 3],
    },
    /// Faces on a surface of this kind.
    Kind(SurfaceKind),
    /// Faces whose label passes.
    Labelled(LabelMatch),
    /// Faces every query picks.
    And(Vec<FaceQuery>),
    /// Faces any query picks.
    Or(Vec<FaceQuery>),
    /// Faces the query doesn't pick.
    Not(Box<FaceQuery>),
}

impl FaceQuery {
    /// Faces both `self` and `other` pick.
    pub fn and(self, other: FaceQuery) -> FaceQuery {
        match self {
            FaceQuery::And(mut all) => {
                all.push(other);
                FaceQuery::And(all)
            }
            first => FaceQuery::And(vec![first, other]),
        }
    }

    /// Faces either `self` or `other` picks.
    pub fn or(self, other: FaceQuery) -> FaceQuery {
        match self {
            FaceQuery::Or(mut any) => {
                any.push(other);
                FaceQuery::Or(any)
            }
            first => FaceQuery::Or(vec![first, other]),
        }
    }

    /// Whether the query picks this face.
    pub fn matches(&self, face: &FaceGeometry, tolerance: &Tolerance) -> bool {
        let along = |d: [f64; 3]| {
            face.normal
                .is_some_and(|n| angle_deg(n, d) <= tolerance.angle_deg)
        };
        match self {
            FaceQuery::All => true,
            FaceQuery::OnPlane { point, normal } => {
                // A point on the face's plane: the plane's own origin, or for
                // a planar NURBS face the middle of its bounds.
                let on_plane = match face.surface {
                    Surface::Plane { origin, .. } => origin,
                    _ => face.bounds_centre(),
                };
                along(*normal)
                    && dot(
                        sub(on_plane, *point),
                        scale(*normal, 1.0 / dot(*normal, *normal).sqrt()),
                    )
                    .abs()
                        <= tolerance.length
            }
            FaceQuery::Facing { direction } => along(*direction),
            FaceQuery::Kind(kind) => SurfaceKind::of(&face.surface) == *kind,
            FaceQuery::Labelled(pattern) => face.label.as_ref().is_some_and(|l| pattern.matches(l)),
            FaceQuery::And(all) => all.iter().all(|q| q.matches(face, tolerance)),
            FaceQuery::Or(any) => any.iter().any(|q| q.matches(face, tolerance)),
            FaceQuery::Not(query) => !query.matches(face, tolerance),
        }
    }

    /// Why the query can't be evaluated, if it can't.
    pub(crate) fn problem(&self) -> Option<&'static str> {
        let finite = |v: &[f64]| v.iter().all(|c| c.is_finite());
        let direction = |d: [f64; 3]| finite(&d) && dot(d, d) > 0.0;
        match self {
            FaceQuery::All | FaceQuery::Kind(_) | FaceQuery::Labelled(_) => None,
            FaceQuery::OnPlane { point, normal } => (!finite(point) || !direction(*normal))
                .then_some("plane needs a finite point and a non-zero normal"),
            FaceQuery::Facing { direction: d } => {
                (!direction(*d)).then_some("direction must be finite and non-zero")
            }
            FaceQuery::And(queries) | FaceQuery::Or(queries) => {
                queries.iter().find_map(FaceQuery::problem)
            }
            FaceQuery::Not(query) => query.problem(),
        }
    }
}

impl Not for FaceQuery {
    type Output = FaceQuery;

    fn not(self) -> FaceQuery {
        match self {
            FaceQuery::Not(query) => *query,
            query => FaceQuery::Not(Box::new(query)),
        }
    }
}

fn angle_deg(a: [f64; 3], b: [f64; 3]) -> f64 {
    let cos = dot(a, b) / (dot(a, a).sqrt() * dot(b, b).sqrt());
    cos.clamp(-1.0, 1.0).acos().to_degrees()
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalised(v: [f64; 3]) -> Option<[f64; 3]> {
    let length = dot(v, v).sqrt();
    (length.is_finite() && length > 1e-12).then(|| scale(v, 1.0 / length))
}
