//! What [`Brep::move_faces`](crate::Brep::move_faces) takes and gives back.

use smlib_sys as sys;

use crate::brep::Brep;

/// A rigid motion for [`Brep::move_faces`](crate::Brep::move_faces): an
/// optional [`Turn`], then a translation. The default moves nothing.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RigidMotion {
    /// Added to every point after the turn.
    pub translation: [f64; 3],
    /// A turn about an axis, before the translation.
    pub turn: Option<Turn>,
}

/// A turn about the axis through `point` along `axis`, by `degrees`,
/// anticlockwise looking down the axis from its tip (the right-hand rule),
/// clockwise when negative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Turn {
    /// A point on the axis.
    pub point: [f64; 3],
    /// The axis's direction; any length but zero.
    pub axis: [f64; 3],
    /// The angle, in degrees.
    pub degrees: f64,
}

impl RigidMotion {
    /// A translation alone.
    pub fn translation(translation: [f64; 3]) -> RigidMotion {
        RigidMotion {
            translation,
            turn: None,
        }
    }

    /// A turn alone, about the axis through `point` along `axis`.
    pub fn turn(point: [f64; 3], axis: [f64; 3], degrees: f64) -> RigidMotion {
        RigidMotion {
            translation: [0.0; 3],
            turn: Some(Turn {
                point,
                axis,
                degrees,
            }),
        }
    }

    /// This motion followed by a translation.
    pub fn then_translate(self, translation: [f64; 3]) -> RigidMotion {
        let [x, y, z] = self.translation;
        RigidMotion {
            translation: [x + translation[0], y + translation[1], z + translation[2]],
            ..self
        }
    }

    /// Where the motion takes `point`.
    pub fn apply(&self, point: [f64; 3]) -> [f64; 3] {
        let turned = match self.turn {
            Some(turn) if turn.degrees != 0.0 => turn.apply(point),
            _ => point,
        };
        [
            turned[0] + self.translation[0],
            turned[1] + self.translation[1],
            turned[2] + self.translation[2],
        ]
    }
}

impl Turn {
    /// Where the turn takes `point` (Rodrigues' rotation formula).
    pub fn apply(&self, point: [f64; 3]) -> [f64; 3] {
        let length = self.axis.iter().map(|c| c * c).sum::<f64>().sqrt();
        let k = self.axis.map(|c| c / length);
        let p = [
            point[0] - self.point[0],
            point[1] - self.point[1],
            point[2] - self.point[2],
        ];
        let (sin, cos) = self.degrees.to_radians().sin_cos();
        let along = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
        let cross = [
            k[1] * p[2] - k[2] * p[1],
            k[2] * p[0] - k[0] * p[2],
            k[0] * p[1] - k[1] * p[0],
        ];
        std::array::from_fn(|i| {
            self.point[i] + p[i] * cos + cross[i] * sin + k[i] * along * (1.0 - cos)
        })
    }
}

/// A solid whose faces [`Brep::move_faces`](crate::Brep::move_faces)
/// moved.
#[derive(Debug)]
pub struct Moved {
    /// The solid, with the faces and edges it had, in their order, each
    /// with its label.
    pub brep: Brep,
    /// What the move did to each face, in face order.
    pub faces: Vec<MovedFace>,
}

/// What [`Brep::move_faces`](crate::Brep::move_faces) did to one face.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MovedFace {
    /// Nothing: it's where it was, as it was.
    Kept,
    /// One of the faces asked for, moved rigidly.
    Moved,
    /// A blend on a moved face (a fillet, a corner between fillets), slid
    /// to stay tangent to the faces round it.
    Carried,
    /// A face sharing an edge or a corner with a moving face, extended or
    /// trimmed along its own surface to meet it.
    Followed,
}

impl MovedFace {
    pub(crate) fn from_raw(raw: u8) -> MovedFace {
        match sys::sm_move_role::from(raw) {
            sys::SM_MOVE_MOVED => MovedFace::Moved,
            sys::SM_MOVE_CARRIED => MovedFace::Carried,
            sys::SM_MOVE_FOLLOWED => MovedFace::Followed,
            _ => MovedFace::Kept,
        }
    }
}
