//! What [`Brep::pipe`](crate::Brep::pipe) gives back.

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
