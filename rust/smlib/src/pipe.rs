//! What [`Brep::pipe`](crate::Brep::pipe) and
//! [`Brep::sweep_profile`](crate::Brep::sweep_profile) give back.

use smlib_sys as sys;

use crate::brep::Brep;

/// A pipe swept along a path by [`Brep::pipe`](crate::Brep::pipe).
#[derive(Debug)]
pub struct Pipe {
    /// The pipe, a new solid.
    pub solid: Brep,
    /// What each of the pipe's faces is, in its face order.
    pub faces: Vec<PipeFace>,
}

/// What one face of a [`Pipe`] is. A wall's `segment` is the index in the
/// path of the curve it was swept along, or `None` when it couldn't be told.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PipeFace {
    /// The outer wall, swept round one curve of the path.
    Outer {
        /// The path curve, by index.
        segment: Option<usize>,
    },
    /// The inner wall of a hollow pipe, swept round one curve of the path.
    Inner {
        /// The path curve, by index.
        segment: Option<usize>,
    },
    /// The cap at the path's start.
    Start,
    /// The cap at the path's end.
    End,
}

impl PipeFace {
    pub(crate) fn from_raw(raw: &sys::sm_pipe_face) -> PipeFace {
        let segment = usize::try_from(raw.segment).ok();
        match raw.role as sys::sm_pipe_role {
            sys::SM_PIPE_INNER => PipeFace::Inner { segment },
            sys::SM_PIPE_START => PipeFace::Start,
            sys::SM_PIPE_END => PipeFace::End,
            _ => PipeFace::Outer { segment },
        }
    }
}

/// A profile swept along a path by
/// [`Brep::sweep_profile`](crate::Brep::sweep_profile).
#[derive(Debug)]
pub struct ProfileSweep {
    /// The swept solid, a new one.
    pub solid: Brep,
    /// What each of its faces is, in its face order.
    pub faces: Vec<ProfileFace>,
}

/// How a profile swept along an open path ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PathEnd {
    /// A flat cap: the profile itself, square to the path.
    Square,
    /// The profile's left half turned half a turn about the frame's z
    /// through the path's end: round beyond it.
    Round,
}

/// What one face of a [`ProfileSweep`] is. `path` is the index in the path
/// of the curve a side was swept along, and `profile` the index in the
/// profile of the segment a side, or a round end's face, was swept from;
/// either is `None` when it couldn't be told, and a square end's cap has no
/// profile segment.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ProfileFace {
    /// A side, swept from one profile segment along one path curve.
    Side {
        /// The path curve, by index.
        path: Option<usize>,
        /// The profile segment, by index.
        profile: Option<usize>,
    },
    /// The path's start: its cap, or a face of its round end.
    Start {
        /// The profile segment a round end's face was turned from.
        profile: Option<usize>,
    },
    /// The path's end, likewise.
    End {
        /// The profile segment a round end's face was turned from.
        profile: Option<usize>,
    },
}

impl ProfileFace {
    pub(crate) fn from_raw(raw: &sys::sm_profile_face) -> ProfileFace {
        let path = usize::try_from(raw.path).ok();
        let profile = usize::try_from(raw.profile).ok();
        match raw.role as sys::sm_profile_role {
            sys::SM_PROFILE_START => ProfileFace::Start { profile },
            sys::SM_PROFILE_END => ProfileFace::End { profile },
            _ => ProfileFace::Side { path, profile },
        }
    }
}
