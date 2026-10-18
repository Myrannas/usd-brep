//! Faces as places to put things: their geometry, frames to place solids in,
//! and prisms swept from them.

use smlib_sys as sys;

use crate::brep::Brep;
use crate::label::Label;
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
