use std::fmt;

use smlib_sys as sys;

/// A status code returned by the kernel or the C shim beneath this crate.
///
/// Kernel codes are SMLib's `SM_*` values (`source/SMLib/inc/SmMessages.h`);
/// codes from 2001 up come from the shim.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Status(pub i64);

impl Status {
    /// SMLib's success code.
    pub const SUCCESS: Status = Status::from_raw(sys::SM_STATUS_SUCCESS);

    // sm_status is a C long: i64 here, but i32 on Windows.
    #[allow(clippy::unnecessary_cast)]
    pub(crate) const fn from_raw(status: sys::sm_status) -> Status {
        Status(status as i64)
    }

    /// The code's symbolic name, if it is a known one.
    pub fn name(self) -> Option<&'static str> {
        Some(match self.0 {
            1000 => "SM_SUCCESS",
            1001 => "SM_ERR",
            1002 => "SM_ERR_OUT_OF_MEMORY",
            1003 => "SM_ERR_UNKNOWN",
            1004 => "SM_ERR_FATAL",
            1005 => "SM_ERR_ASSERT_FAILURE",
            1006 => "SM_ERR_NULL_POINTER",
            1007 => "SM_ERR_INVALID_INPUT",
            1008 => "SM_ERR_NON_NULL_OUTPUT_POINTER",
            1010 => "SM_ERR_ASSERTVALID_FAILURE",
            1020 => "SM_ERR_OUTSIDE_OF_DOMAIN",
            1021 => "SM_ERR_NOT_WITHIN_TOLERANCE",
            1022 => "SM_ERR_NOT_CONVERGING",
            1023 => "SM_ERR_WARNING",
            1024 => "SM_ERR_MESSAGE",
            1025 => "SM_ERR_AXIS_INSIDE_FACE",
            1026 => "SM_ERR_LICENSE_EXPIRED",
            1027 => "SM_ERR_BAD_INTERSECTIONS",
            1028 => "SM_ERR_BAD_COINCIDENT_VERTICES",
            1030 => "SM_ERR_BAD_SURFACE_POINT",
            1032 => "SM_ERR_BAD_TANGENT_DROP",
            1034 => "SM_ERR_DEGENERATE_SURFACE",
            1036 => "SM_ERR_BAD_FIND_DEGEN_PARAM",
            1038 => "SM_ERR_METHOD_FAILURE_QUITING",
            1040 => "SM_ERR_NOTYET_HEAL_FACE",
            1041 => "SM_ERR_TESS_FAIL_FACE",
            2001 => "SHIM_ERR_NULL_ARGUMENT",
            2002 => "SHIM_ERR_EXCEPTION",
            2003 => "SHIM_ERR_OUT_OF_MEMORY",
            2004 => "SHIM_ERR_UNMAPPED_FACE",
            2005 => "SHIM_ERR_UNSUPPORTED_POLYGON",
            2006 => "SHIM_ERR_BAD_EDGE",
            2007 => "SHIM_ERR_BAD_FILLET",
            2008 => "SHIM_ERR_BROKE_SOLID",
            2020 => "SHIM_ERR_BAD_FACE",
            2021 => "SHIM_ERR_BAD_LABEL",
            2022 => "SHIM_ERR_NOT_PLANAR",
            2023 => "SHIM_ERR_BAD_FRAME",
            2024 => "SHIM_ERR_BAD_SWEEP",
            2025 => "SHIM_ERR_UNSUPPORTED_FACE",
            2026 => "SHIM_ERR_CHANGED_TOPOLOGY",
            2027 => "SHIM_ERR_OBSTRUCTED",
            2028 => "SHIM_ERR_UNSUPPORTED_STEP",
            2030 => "SHIM_ERR_SPENT_TESSELLATION",
            2031 => "SHIM_ERR_BAD_AXIS",
            2032 => "SHIM_ERR_AXIS_CROSSES_FACE",
            2033 => "SHIM_ERR_BAD_DRAFT",
            2034 => "SHIM_ERR_DRAFT_CLOSES",
            2035 => "SHIM_ERR_BAD_PATH",
            2036 => "SHIM_ERR_ARC_CORNER",
            _ => return None,
        })
    }
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.name() {
            Some(name) => write!(f, "{name} ({})", self.0),
            None => write!(f, "status {}", self.0),
        }
    }
}

/// Why a kernel operation failed.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// A kernel call returned a failure status.
    Kernel {
        /// The kernel function that failed, such as `SmApiBooleanDifference`.
        operation: &'static str,
        /// The status it returned.
        status: Status,
    },
    /// An argument was rejected before it reached the kernel.
    InvalidArgument {
        /// The operation the argument was for.
        operation: &'static str,
        /// What was wrong with it.
        reason: &'static str,
    },
    /// [`Kernel::init`](crate::Kernel::init) was called on a thread other than
    /// the one that owns the kernel.
    WrongThread,
}

impl Error {
    pub(crate) fn check(operation: &'static str, status: sys::sm_status) -> Result<(), Error> {
        if status == sys::SM_STATUS_SUCCESS {
            Ok(())
        } else {
            Err(Error::Kernel {
                operation,
                status: Status::from_raw(status),
            })
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Kernel { operation, status } => write!(f, "{operation} failed: {status}"),
            Error::InvalidArgument { operation, reason } => write!(f, "{operation}: {reason}"),
            Error::WrongThread => write!(f, "the SMLib kernel belongs to another thread"),
        }
    }
}

impl std::error::Error for Error {}
