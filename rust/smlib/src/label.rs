//! Labels: caller-chosen names on faces and edges that the kernel carries
//! through every operation.

use std::ffi::c_char;
use std::fmt;

use smlib_sys as sys;

use crate::error::Error;

/// The longest role a [`Label`] can carry, in bytes.
pub const MAX_ROLE_LEN: usize = sys::SM_LABEL_ROLE_CAPACITY as usize - 1;

/// A caller's name for a face or edge: the feature that made it, the role it
/// plays there, and an index telling apart entities one feature makes with
/// the same role.
///
/// Set labels with [`Brep::set_face_label`](crate::Brep::set_face_label) and
/// [`Brep::set_edge_label`](crate::Brep::set_edge_label); read them back from
/// [`Brep::faces`](crate::Brep::faces), [`Brep::edges`](crate::Brep::edges)
/// and [`Mesh`](crate::Mesh). The kernel carries them as attributes
/// (`SmAttribute.h`):
///
/// - When an operation splits a face or edge, both pieces keep its label.
/// - Copies ([`Brep::try_clone`](crate::Brep::try_clone)) keep every label.
/// - When two labelled faces or edges merge into one, the result keeps the
///   older label: the one with the lower `feature`.
/// - Faces and edges an operation creates from nothing, such as a fillet's
///   rounds and their rails, have no label.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Label {
    /// The feature that made the face or edge. Lower ids are older.
    pub feature: u32,
    /// What the face or edge is in that feature, such as `"wall"`: not empty,
    /// at most [`MAX_ROLE_LEN`] bytes and without NUL characters.
    pub role: String,
    /// Tells apart entities one feature makes with the same role.
    pub index: u32,
}

impl Label {
    /// A label with index 0.
    pub fn new(feature: u32, role: impl Into<String>) -> Label {
        Label {
            feature,
            role: role.into(),
            index: 0,
        }
    }

    /// The same label with another index.
    pub fn with_index(self, index: u32) -> Label {
        Label { index, ..self }
    }

    pub(crate) fn to_raw(&self, operation: &'static str) -> Result<sys::sm_label, Error> {
        let bytes = self.role.as_bytes();
        if bytes.is_empty() || bytes.len() > MAX_ROLE_LEN || bytes.contains(&0) {
            return Err(Error::InvalidArgument {
                operation,
                reason: "a label's role must be 1 to 63 bytes without NUL",
            });
        }
        let mut raw = NO_LABEL;
        raw.present = 1;
        raw.feature = self.feature;
        raw.index = self.index;
        for (slot, &byte) in raw.role.iter_mut().zip(bytes) {
            *slot = byte as c_char;
        }
        Ok(raw)
    }

    pub(crate) fn from_raw(raw: &sys::sm_label) -> Option<Label> {
        if raw.present == 0 {
            return None;
        }
        let bytes: Vec<u8> = raw
            .role
            .iter()
            .map(|&c| c as u8)
            .take_while(|&b| b != 0)
            .collect();
        Some(Label {
            feature: raw.feature,
            role: String::from_utf8_lossy(&bytes).into_owned(),
            index: raw.index,
        })
    }
}

impl fmt::Display for Label {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}.{}", self.feature, self.role, self.index)
    }
}

/// The raw form of "no label".
pub(crate) const NO_LABEL: sys::sm_label = sys::sm_label {
    present: 0,
    feature: 0,
    index: 0,
    role: [0; sys::SM_LABEL_ROLE_CAPACITY as usize],
};

/// A test on a [`Label`] for [`EdgeQuery::Between`](crate::EdgeQuery::Between):
/// the feature must match, and so must the role and index when given.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct LabelMatch {
    /// The feature the label must name.
    pub feature: u32,
    /// The role it must have, or any.
    pub role: Option<String>,
    /// The index it must have, or any.
    pub index: Option<u32>,
}

impl LabelMatch {
    /// Labels from this feature with this role, whatever their index.
    pub fn role(feature: u32, role: impl Into<String>) -> LabelMatch {
        LabelMatch {
            feature,
            role: Some(role.into()),
            index: None,
        }
    }

    /// Whether the label passes.
    pub fn matches(&self, label: &Label) -> bool {
        label.feature == self.feature
            && self.role.as_ref().is_none_or(|role| *role == label.role)
            && self.index.is_none_or(|index| index == label.index)
    }
}
