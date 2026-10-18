use std::marker::PhantomData;
use std::ptr::{self, NonNull};

use smlib_sys as sys;

use crate::brep::Brep;
use crate::error::Error;
use crate::mesh::Mesh;

/// A tessellation of a [`Brep`] in progress, from
/// [`Brep::begin_tessellation`]: every edge split as
/// [`Brep::tessellate`] splits it, with the faces still to triangulate.
///
/// A face's triangles depend only on its own geometry and its boundary as
/// split, which its [`FaceKey`] captures, so they come out the same whichever
/// other faces [`mesh`](Tessellation::mesh) triangulates. A caller can keep
/// earlier triangles by key and ask only for the faces it hasn't seen; the
/// pieces stay watertight because every shared edge is split the same on
/// both sides.
pub struct Tessellation<'a> {
    raw: NonNull<sys::sm_tessellation>,
    _brep: PhantomData<&'a Brep>,
    _not_send: PhantomData<*mut ()>,
}

impl Tessellation<'_> {
    /// # Safety
    ///
    /// `raw` must be a live tessellation from `sm_brep_tessellation_begin`
    /// of a brep that outlives the returned value.
    pub(crate) unsafe fn from_raw(raw: NonNull<sys::sm_tessellation>) -> Self {
        Tessellation {
            raw,
            _brep: PhantomData,
            _not_send: PhantomData,
        }
    }

    /// The number of faces: the brep's.
    pub fn face_count(&self) -> usize {
        // SAFETY: raw is live.
        unsafe { sys::sm_tessellation_face_count(self.raw.as_ptr()) }
    }

    /// Each face's key, in face order: a 128-bit digest of everything that
    /// fixes the face's triangles, from the tessellation parameters, the
    /// face's surface, trim and tolerances, and its boundary as split. Faces
    /// with equal keys are triangulated identically, vertex for vertex.
    /// `None` for a face SMLib can't describe that way, which should always
    /// be triangulated afresh.
    pub fn face_keys(&mut self) -> Result<Vec<Option<FaceKey>>, Error> {
        (0..self.face_count())
            .map(|face| {
                let mut key = [0u64; 2];
                let mut keyed = 0;
                // SAFETY: raw is live and on the kernel thread; face is in
                // range and key holds the two words the shim writes.
                let status = unsafe {
                    sys::sm_tessellation_face_key(
                        self.raw.as_ptr(),
                        face,
                        key.as_mut_ptr(),
                        &mut keyed,
                    )
                };
                Error::check("SmApiTessellate", status)?;
                Ok((keyed != 0).then_some(FaceKey(key)))
            })
            .collect()
    }

    /// Triangulates the faces `wanted` marks, one flag per face (`None` for
    /// every face), into a [`Mesh`] whose triangles and vertices cover just
    /// those faces. Faces follow each other in face order, each with its own
    /// vertices, so a face's triangles form one run of [`Mesh::face_ids`]
    /// and its vertices one run of the vertex arrays. [`Mesh::faces`] and
    /// [`Mesh::edges`] cover the whole brep. With every face wanted, each
    /// face's triangles are [`Brep::tessellate`]'s for it.
    pub fn mesh(self, wanted: Option<&[bool]>) -> Result<Mesh, Error> {
        let flags: Option<Vec<u8>> = wanted.map(|w| w.iter().map(|&f| u8::from(f)).collect());
        let (flags_ptr, flags_len) = flags
            .as_ref()
            .map_or((ptr::null(), 0), |f| (f.as_ptr(), f.len()));
        let mut out = ptr::null_mut();
        // SAFETY: raw is live and on the kernel thread; the flags outlive the
        // call and cover flags_len bytes.
        let status =
            unsafe { sys::sm_tessellation_mesh(self.raw.as_ptr(), flags_ptr, flags_len, &mut out) };
        Error::check("SmApiTessellate", status)?;
        let raw = NonNull::new(out).ok_or(Error::Kernel {
            operation: "SmApiTessellate",
            status: crate::Status::from_raw(sys::SM_SHIM_ERR_NULL_ARGUMENT),
        })?;
        // SAFETY: raw is a mesh the shim just allocated; we free it once,
        // after copying its arrays out.
        unsafe {
            let mesh = Mesh::from_raw(raw.as_ptr());
            sys::sm_mesh_free(raw.as_ptr());
            Ok(mesh)
        }
    }
}

impl Drop for Tessellation<'_> {
    fn drop(&mut self) {
        // SAFETY: we own the tessellation, and it can't leave the kernel
        // thread.
        unsafe { sys::sm_tessellation_free(self.raw.as_ptr()) }
    }
}

/// What a face's triangles depend on, from [`Tessellation::face_keys`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FaceKey(pub(crate) [u64; 2]);
