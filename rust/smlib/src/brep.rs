use std::fmt;
use std::marker::PhantomData;
use std::mem::ManuallyDrop;
use std::ptr::{self, NonNull};

use smlib_sys as sys;

use crate::error::Error;
use crate::kernel::Kernel;
use crate::mesh::{Mesh, TessellationParams};

/// A B-rep solid owned by the kernel. It is freed on drop.
///
/// A `Brep` is neither `Send` nor `Sync`: it stays on the kernel's thread.
pub struct Brep {
    raw: NonNull<sys::sm_brep>,
    _not_send: PhantomData<*mut ()>,
}

type BooleanFn = unsafe extern "C" fn(
    *mut sys::sm_brep,
    *mut sys::sm_brep,
    *mut *mut sys::sm_brep,
) -> sys::sm_status;

impl Brep {
    /// Wraps a brep the shim just handed over.
    fn from_raw(operation: &'static str, raw: *mut sys::sm_brep) -> Result<Brep, Error> {
        match NonNull::new(raw) {
            Some(raw) => Ok(Brep {
                raw,
                _not_send: PhantomData,
            }),
            None => Err(Error::Kernel {
                operation,
                status: crate::Status::from_raw(sys::SM_SHIM_ERR_NULL_ARGUMENT),
            }),
        }
    }

    /// An axis-aligned box from its minimum corner and its size along x, y
    /// and z (`SmApiCreateBox`).
    pub fn cuboid(_kernel: &Kernel, min_corner: [f64; 3], size: [f64; 3]) -> Result<Brep, Error> {
        const OP: &str = "SmApiCreateBox";
        check_point(OP, min_corner)?;
        for extent in size {
            check_positive(OP, extent, "size must be finite and positive")?;
        }
        let mut out = ptr::null_mut();
        // SAFETY: plain values in, a fresh brep out; we're on the kernel thread.
        let status = unsafe {
            sys::sm_brep_create_box(vec3(min_corner), size[0], size[1], size[2], &mut out)
        };
        Error::check(OP, status)?;
        Brep::from_raw(OP, out)
    }

    /// A solid cylinder standing on `base_center`, along +z
    /// (`SmApiCreateCylinder`).
    pub fn cylinder(
        _kernel: &Kernel,
        base_center: [f64; 3],
        radius: f64,
        height: f64,
    ) -> Result<Brep, Error> {
        const OP: &str = "SmApiCreateCylinder";
        check_point(OP, base_center)?;
        check_positive(OP, radius, "radius must be finite and positive")?;
        check_positive(OP, height, "height must be finite and positive")?;
        let mut out = ptr::null_mut();
        // SAFETY: as in `cuboid`.
        let status =
            unsafe { sys::sm_brep_create_cylinder(vec3(base_center), radius, height, &mut out) };
        Error::check(OP, status)?;
        Brep::from_raw(OP, out)
    }

    /// A solid sphere (`SmApiCreateSphere`).
    pub fn sphere(_kernel: &Kernel, center: [f64; 3], radius: f64) -> Result<Brep, Error> {
        const OP: &str = "SmApiCreateSphere";
        check_point(OP, center)?;
        check_positive(OP, radius, "radius must be finite and positive")?;
        let mut out = ptr::null_mut();
        // SAFETY: as in `cuboid`.
        let status = unsafe { sys::sm_brep_create_sphere(vec3(center), radius, &mut out) };
        Error::check(OP, status)?;
        Brep::from_raw(OP, out)
    }

    /// A deep copy (`SmApiBrepCopy`).
    pub fn try_clone(&self) -> Result<Brep, Error> {
        const OP: &str = "SmApiBrepCopy";
        let mut out = ptr::null_mut();
        // SAFETY: self is a live brep on the kernel thread.
        let status = unsafe { sys::sm_brep_copy(self.raw.as_ptr(), &mut out) };
        Error::check(OP, status)?;
        Brep::from_raw(OP, out)
    }

    /// `self ∪ other` (`SmApiBooleanUnion`). Both inputs are consumed; on
    /// failure the kernel may have modified them, so both are dropped.
    pub fn union(self, other: Brep) -> Result<Brep, Error> {
        self.boolean(other, sys::sm_brep_union, "SmApiBooleanUnion")
    }

    /// `self − other` (`SmApiBooleanDifference`). Consumes both, as
    /// [`union`](Brep::union).
    pub fn difference(self, other: Brep) -> Result<Brep, Error> {
        self.boolean(other, sys::sm_brep_difference, "SmApiBooleanDifference")
    }

    /// `self ∩ other` (`SmApiBooleanIntersection`). Consumes both, as
    /// [`union`](Brep::union).
    pub fn intersection(self, other: Brep) -> Result<Brep, Error> {
        self.boolean(other, sys::sm_brep_intersection, "SmApiBooleanIntersection")
    }

    fn boolean(self, other: Brep, op: BooleanFn, name: &'static str) -> Result<Brep, Error> {
        let a = ManuallyDrop::new(self);
        let b = ManuallyDrop::new(other);
        let mut out = ptr::null_mut();
        // SAFETY: a and b are distinct live breps (each owned by value) on the
        // kernel thread. The kernel takes ownership of both for the call.
        let status = unsafe { op(a.raw.as_ptr(), b.raw.as_ptr(), &mut out) };
        if status == sys::SM_STATUS_SUCCESS {
            // a was modified in place and returned as out; b was deleted.
            debug_assert_eq!(
                out,
                a.raw.as_ptr(),
                "{name} should return its first operand"
            );
            return Brep::from_raw(name, out);
        }
        // Both are still allocated but may be half-modified: discard them.
        drop(ManuallyDrop::into_inner(a));
        drop(ManuallyDrop::into_inner(b));
        Err(Error::Kernel {
            operation: name,
            status: crate::Status::from_raw(status),
        })
    }

    /// Whether this is a closed, manifold solid (`SmApiBrepIsManifoldSolid`).
    pub fn is_manifold(&self) -> Result<bool, Error> {
        let mut manifold = 0;
        // SAFETY: self is a live brep on the kernel thread.
        let status = unsafe { sys::sm_brep_is_manifold_solid(self.raw.as_ptr(), &mut manifold) };
        Error::check("SmApiBrepIsManifoldSolid", status)?;
        Ok(manifold != 0)
    }

    /// The enclosed volume (`SmApiBrepComputeVolume`). The kernel clamps
    /// `relative_accuracy` to [1e-8, 1e-1]; tighter is slower but avoids the
    /// known error on faces trimmed by circles (KNOWN_ISSUES.md).
    pub fn volume(&self, relative_accuracy: f64) -> Result<f64, Error> {
        const OP: &str = "SmApiBrepComputeVolume";
        check_positive(
            OP,
            relative_accuracy,
            "relative accuracy must be finite and positive",
        )?;
        let mut volume = 0.0;
        // SAFETY: self is a live brep on the kernel thread.
        let status =
            unsafe { sys::sm_brep_volume(self.raw.as_ptr(), relative_accuracy, &mut volume) };
        Error::check(OP, status)?;
        Ok(volume)
    }

    /// The number of faces in the brep.
    pub fn face_count(&self) -> usize {
        self.counts().0
    }

    /// The number of edges in the brep.
    pub fn edge_count(&self) -> usize {
        self.counts().1
    }

    fn counts(&self) -> (usize, usize) {
        let (mut faces, mut edges) = (0, 0);
        // SAFETY: self is a live brep on the kernel thread; this only reads
        // list sizes and can't fail for a non-null brep.
        unsafe { sys::sm_brep_counts(self.raw.as_ptr(), &mut faces, &mut edges) };
        (faces, edges)
    }

    /// Tessellates the solid (`SmApiTessellate`) and samples its edges
    /// (`SmApiTessellateBoundaries`) into a [`Mesh`] that keeps its links to
    /// the brep's faces and edges.
    pub fn tessellate(&self, params: &TessellationParams) -> Result<Mesh, Error> {
        let raw_params = params.to_raw()?;
        let mut out = ptr::null_mut();
        // SAFETY: self is a live brep on the kernel thread; params outlives
        // the call.
        let status = unsafe { sys::sm_brep_tessellate(self.raw.as_ptr(), &raw_params, &mut out) };
        Error::check("SmApiTessellate", status)?;
        let Some(raw) = NonNull::new(out) else {
            return Err(Error::Kernel {
                operation: "SmApiTessellate",
                status: crate::Status::from_raw(sys::SM_SHIM_ERR_NULL_ARGUMENT),
            });
        };
        // SAFETY: raw is a mesh the shim just allocated; we free it once,
        // after copying its arrays out.
        unsafe {
            let mesh = Mesh::from_raw(raw.as_ptr());
            sys::sm_mesh_free(raw.as_ptr());
            Ok(mesh)
        }
    }
}

impl Drop for Brep {
    fn drop(&mut self) {
        // SAFETY: we own the brep, and Drop runs on the kernel thread because
        // Brep isn't Send.
        unsafe { sys::sm_brep_free(self.raw.as_ptr()) }
    }
}

impl fmt::Debug for Brep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Brep")
            .field("faces", &self.face_count())
            .field("edges", &self.edge_count())
            .finish()
    }
}

fn vec3([x, y, z]: [f64; 3]) -> sys::sm_vec3 {
    sys::sm_vec3 { x, y, z }
}

fn check_point(operation: &'static str, point: [f64; 3]) -> Result<(), Error> {
    if point.iter().all(|c| c.is_finite()) {
        Ok(())
    } else {
        Err(Error::InvalidArgument {
            operation,
            reason: "coordinates must be finite",
        })
    }
}

fn check_positive(operation: &'static str, value: f64, reason: &'static str) -> Result<(), Error> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(Error::InvalidArgument { operation, reason })
    }
}
