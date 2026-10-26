//! What [`Brep::shell`](crate::Brep::shell) gives back.

use smlib_sys as sys;

use crate::brep::Brep;

/// A solid hollowed by [`Brep::shell`](crate::Brep::shell).
#[derive(Debug)]
pub struct Shell {
    /// The hollowed solid, apart from the one it was made from.
    pub brep: Brep,
    /// What each of its faces is, in its face order.
    pub faces: Vec<ShellFace>,
}

/// What one face of a [`Shell`] is. `face` is the face of the solid it came
/// from, as its index in that solid's [`faces`](crate::Brep::faces).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ShellFace {
    /// A face that stays where it was, outside, with its label.
    Outer {
        /// The face it is (or a piece of).
        face: usize,
    },
    /// What's left of an opened face round its opening, with its label.
    Rim {
        /// The face opened.
        face: usize,
    },
    /// A face's offset: the inside of the wall it stands for. No label.
    Inner {
        /// The face offset.
        face: usize,
    },
    /// A face that came from no face of the solid. No label.
    New,
}

impl ShellFace {
    pub(crate) fn from_raw(raw: &sys::sm_shell_face) -> ShellFace {
        let Ok(face) = usize::try_from(raw.source_face) else {
            return ShellFace::New;
        };
        match raw.role as sys::sm_shell_role {
            sys::SM_SHELL_OUTER => ShellFace::Outer { face },
            sys::SM_SHELL_RIM => ShellFace::Rim { face },
            sys::SM_SHELL_INNER => ShellFace::Inner { face },
            _ => ShellFace::New,
        }
    }
}
