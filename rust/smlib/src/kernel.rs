use std::marker::PhantomData;
use std::sync::OnceLock;
use std::thread::{self, ThreadId};

use smlib_sys as sys;

use crate::error::{Error, Status};

/// The thread that owns the kernel, and the status its context creation
/// returned. Set once per process.
static OWNER: OnceLock<(ThreadId, sys::sm_status)> = OnceLock::new();

/// A handle to the SMLib kernel, proof that the caller is on the kernel's
/// thread.
///
/// Every [`Brep`](crate::Brep) is created through a `Kernel`. The handle is
/// neither `Send` nor `Sync`, so it and everything made with it stay on the
/// thread that called [`Kernel::init`] first.
#[derive(Debug)]
pub struct Kernel {
    _not_send: PhantomData<*mut ()>,
}

impl Kernel {
    /// Returns a handle to the kernel, creating SMLib's process-global context
    /// (`SmApiCreateContext`) on the first call.
    ///
    /// The first thread to call this owns the kernel for the life of the
    /// process. Later calls on that thread return another handle; calls on any
    /// other thread fail with [`Error::WrongThread`]. Call it on the kernel's
    /// thread before starting threads that might race it.
    pub fn init() -> Result<Kernel, Error> {
        let me = thread::current().id();
        // SAFETY: OnceLock runs this exactly once per process, so the context
        // is created once, on the thread that becomes its owner.
        let (owner, status) = *OWNER.get_or_init(|| (me, unsafe { sys::sm_context_create() }));
        if owner != me {
            return Err(Error::WrongThread);
        }
        if status != sys::SM_STATUS_SUCCESS {
            return Err(Error::Kernel {
                operation: "SmApiCreateContext",
                status: Status::from_raw(status),
            });
        }
        Ok(Kernel {
            _not_send: PhantomData,
        })
    }
}
