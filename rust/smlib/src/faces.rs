//! Faces as places to put things: their geometry, frames to place solids in,
//! prisms swept from them, and queries that pick them.

use std::ops::Not;

use smlib_sys as sys;

use crate::brep::Brep;
use crate::edges::Tolerance;
use crate::label::{Label, LabelMatch};
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

/// A closed prism swept from a planar face by
/// [`Brep::sweep_face`](crate::Brep::sweep_face).
#[derive(Debug)]
pub struct Sweep {
    /// The prism, a new solid apart from the brep it was swept from.
    pub prism: Brep,
    /// What each of the prism's faces came from, in the prism's face order.
    pub faces: Vec<SweptFace>,
}

/// What one face of a [`Sweep`] came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SweptFace {
    /// The copy of the swept face where the sweep starts, on the face itself.
    Start,
    /// The copy of the swept face the sweep distance away.
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
