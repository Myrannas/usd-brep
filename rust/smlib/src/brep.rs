use std::fmt;
use std::marker::PhantomData;
use std::mem::ManuallyDrop;
use std::ptr::{self, NonNull};

use smlib_sys as sys;

use crate::edges::{EDGE_SAMPLES, EdgeGeometry, EdgeQuery, Tolerance};
use crate::error::Error;
use crate::faces::{
    BoundaryCurve, FaceCentroid, FaceGeometry, FaceQuery, Frame, PlaneCurve, Sweep,
};
use crate::fillet::{EdgeBlend, FilletOptions};
use crate::imprint::{Imprint, ImprintEdge};
use crate::kernel::Kernel;
use crate::label::{Label, NO_LABEL};
use crate::mesh::{Mesh, TessellationParams, slice};
use crate::motion::{Moved, MovedFace, RigidMotion};
use crate::pipe::{Pipe, PipeFace};
use crate::shell::{Shell, ShellFace};
use crate::tessellation::{FaceKey, SurfaceCaches, Tessellation};

/// Relative accuracy of [`Brep::face_centroid`]: placing a feature at a
/// centroid is then exact to a tenth of a micron on a metre-sized part.
const CENTROID_ACCURACY: f64 = 1e-7;

/// Relative accuracy of [`Brep::edge_length`].
const EDGE_LENGTH_ACCURACY: f64 = 1e-7;

/// How far apart, in model units, one curve's end and the next one's start
/// may be in a loop for [`Brep::planar_faces`] or a path for [`Brep::pipe`].
const LOOP_GAP: f64 = 1e-6;

/// The sharpest turn, in degrees, a [`Brep::pipe`]'s path may take where
/// two curves meet: past it a mitre's corner would reach too far.
const MAX_PIPE_TURN_DEG: f64 = 179.0;

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
    pub(crate) fn from_raw(operation: &'static str, raw: *mut sys::sm_brep) -> Result<Brep, Error> {
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

    /// A solid frustum standing on `base_center`, along +z, with the given
    /// base and top radii; one of them may be 0 for a cone
    /// (`SmApiCreateCone`).
    pub fn cone(
        _kernel: &Kernel,
        base_center: [f64; 3],
        base_radius: f64,
        top_radius: f64,
        height: f64,
    ) -> Result<Brep, Error> {
        const OP: &str = "SmApiCreateCone";
        check_point(OP, base_center)?;
        let radius_ok = |r: f64| r.is_finite() && r >= 0.0;
        if !radius_ok(base_radius) || !radius_ok(top_radius) || base_radius + top_radius <= 0.0 {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "radii must be finite, not negative, and not both zero",
            });
        }
        check_positive(OP, height, "height must be finite and positive")?;
        let mut out = ptr::null_mut();
        // SAFETY: as in `cuboid`.
        let status = unsafe {
            sys::sm_brep_create_cone(vec3(base_center), base_radius, top_radius, height, &mut out)
        };
        Error::check(OP, status)?;
        Brep::from_raw(OP, out)
    }

    /// A solid ring around the +z axis through `center` (`SmApiCreateTorus`).
    pub fn torus(
        _kernel: &Kernel,
        center: [f64; 3],
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Brep, Error> {
        const OP: &str = "SmApiCreateTorus";
        check_point(OP, center)?;
        check_positive(OP, minor_radius, "minor radius must be finite and positive")?;
        if !(major_radius.is_finite() && major_radius > minor_radius) {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "major radius must be finite and larger than the minor radius",
            });
        }
        let mut out = ptr::null_mut();
        // SAFETY: as in `cuboid`.
        let status = unsafe {
            sys::sm_brep_create_torus(vec3(center), major_radius, minor_radius, &mut out)
        };
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

    /// A prism over the closed polygon through `base` (in order, in one
    /// plane), `height` along the polygon's normal, with its walls leaning in
    /// by `draft_deg` (0 to under 90) from that direction and sharp corners
    /// (`SmApiCreateTaperExtrude`). Its tapered walls are planar NURBS
    /// surfaces. The kernel has no outward draft: angles over 90 mirror back.
    pub fn tapered_prism(
        _kernel: &Kernel,
        base: &[[f64; 3]],
        height: f64,
        draft_deg: f64,
    ) -> Result<Brep, Error> {
        const OP: &str = "SmApiCreateTaperExtrude";
        if base.len() < 3 {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "the base needs at least three points",
            });
        }
        for &point in base {
            check_point(OP, point)?;
        }
        check_positive(OP, height, "height must be finite and positive")?;
        if !(0.0..90.0).contains(&draft_deg) {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "draft angle must be from 0 to under 90 degrees",
            });
        }
        let points: Vec<sys::sm_vec3> = base.iter().map(|&p| vec3(p)).collect();
        let mut out = ptr::null_mut();
        // SAFETY: points holds points.len() values; a fresh brep comes out.
        let status = unsafe {
            sys::sm_brep_create_tapered_prism(
                points.as_ptr(),
                points.len(),
                height,
                draft_deg,
                &mut out,
            )
        };
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

    /// The geometry of every edge, in the brep's edge order, for seeing what
    /// an [`EdgeQuery`] will pick.
    pub fn edges(&self) -> Result<Vec<EdgeGeometry>, Error> {
        const OP: &str = "sm_brep_edge_geometry";
        let mut raw = vec![EMPTY_EDGE; self.edge_count()];
        let mut count = 0;
        // SAFETY: self is a live brep on the kernel thread; raw has room for
        // raw.len() edges.
        let status = unsafe {
            sys::sm_brep_edge_geometry(self.raw.as_ptr(), raw.as_mut_ptr(), raw.len(), &mut count)
        };
        Error::check(OP, status)?;
        raw.truncate(count);
        Ok(raw
            .iter()
            .enumerate()
            .map(|(index, raw)| EdgeGeometry::from_raw(index, raw))
            .collect())
    }

    /// The arc length of the edge at `edge` in [`edges`](Brep::edges)
    /// (`SmApiEdgeComputeLength`), integrated to a relative accuracy of 1e-7.
    /// This costs far more than [`edges`](Brep::edges), which leaves it out.
    pub fn edge_length(&self, edge: usize) -> Result<f64, Error> {
        let mut length = 0.0;
        // SAFETY: self is a live brep on the kernel thread.
        let status = unsafe {
            sys::sm_brep_edge_length(self.raw.as_ptr(), edge, EDGE_LENGTH_ACCURACY, &mut length)
        };
        Error::check("SmApiEdgeComputeLength", status)?;
        Ok(length)
    }

    /// The geometry of every face, in the brep's face order: surface, label,
    /// bounds, and the outward normal of planar faces.
    pub fn faces(&self) -> Result<Vec<FaceGeometry>, Error> {
        const OP: &str = "sm_brep_face_geometry";
        let mut raw = vec![EMPTY_FACE; self.face_count()];
        let mut count = 0;
        // SAFETY: self is a live brep on the kernel thread; raw has room for
        // raw.len() faces.
        let status = unsafe {
            sys::sm_brep_face_geometry(self.raw.as_ptr(), raw.as_mut_ptr(), raw.len(), &mut count)
        };
        Error::check(OP, status)?;
        raw.truncate(count);
        Ok(raw
            .iter()
            .enumerate()
            .map(|(index, raw)| FaceGeometry::from_raw(index, raw))
            .collect())
    }

    /// The area and area centroid of the face at `face` in
    /// [`faces`](Brep::faces), integrated to a relative accuracy of 1e-7.
    /// This costs far more than [`faces`](Brep::faces): up to milliseconds
    /// for a face trimmed by many curves.
    pub fn face_centroid(&self, face: usize) -> Result<FaceCentroid, Error> {
        let (mut area, mut centroid) = (0.0, [0.0; 3]);
        // SAFETY: self is a live brep on the kernel thread; centroid holds
        // three doubles.
        let status = unsafe {
            sys::sm_brep_face_centroid(
                self.raw.as_ptr(),
                face,
                CENTROID_ACCURACY,
                &mut area,
                centroid.as_mut_ptr(),
            )
        };
        Error::check("SmFace::ComputePreciseProperties", status)?;
        Ok(FaceCentroid { area, centroid })
    }

    /// An axis-aligned box enclosing the brep, as its minimum and maximum
    /// corners. It may be loose around curved faces.
    pub fn bounding_box(&self) -> Result<[[f64; 3]; 2], Error> {
        let (mut min, mut max) = ([0.0; 3], [0.0; 3]);
        // SAFETY: self is a live brep on the kernel thread; min and max hold
        // three doubles each.
        let status = unsafe {
            sys::sm_brep_bounding_box(self.raw.as_ptr(), min.as_mut_ptr(), max.as_mut_ptr())
        };
        Error::check("SmApiBrepBoundingBox", status)?;
        Ok([min, max])
    }

    /// Labels the face at `face` in [`faces`](Brep::faces), replacing any
    /// label it has. Labels don't change the brep's shape, so face and edge
    /// indices stay valid.
    pub fn set_face_label(&mut self, face: usize, label: &Label) -> Result<(), Error> {
        const OP: &str = "sm_brep_set_face_label";
        let raw = label.to_raw(OP)?;
        // SAFETY: self is a live brep on the kernel thread; raw outlives the
        // call.
        let status = unsafe { sys::sm_brep_set_face_label(self.raw.as_ptr(), face, &raw) };
        Error::check(OP, status)
    }

    /// Removes the label from the face at `face`, if it has one.
    pub fn remove_face_label(&mut self, face: usize) -> Result<(), Error> {
        // SAFETY: as in `set_face_label`.
        let status = unsafe { sys::sm_brep_set_face_label(self.raw.as_ptr(), face, &NO_LABEL) };
        Error::check("sm_brep_set_face_label", status)
    }

    /// Labels the edge at `edge` in [`edges`](Brep::edges), as
    /// [`set_face_label`](Brep::set_face_label) does for faces.
    pub fn set_edge_label(&mut self, edge: usize, label: &Label) -> Result<(), Error> {
        const OP: &str = "sm_brep_set_edge_label";
        let raw = label.to_raw(OP)?;
        // SAFETY: as in `set_face_label`.
        let status = unsafe { sys::sm_brep_set_edge_label(self.raw.as_ptr(), edge, &raw) };
        Error::check(OP, status)
    }

    /// Removes the label from the edge at `edge`, if it has one.
    pub fn remove_edge_label(&mut self, edge: usize) -> Result<(), Error> {
        // SAFETY: as in `set_face_label`.
        let status = unsafe { sys::sm_brep_set_edge_label(self.raw.as_ptr(), edge, &NO_LABEL) };
        Error::check("sm_brep_set_edge_label", status)
    }

    /// Moves the brep rigidly into `frame` (`SmApiTransform`): the point
    /// (x, y, z) goes to `frame.to_world([x, y, z])`. Build a solid at the
    /// origin, then place it, to put it anywhere at any angle.
    pub fn place(&mut self, frame: &Frame) -> Result<(), Error> {
        // SAFETY: self is a live brep on the kernel thread.
        let status = unsafe {
            sys::sm_brep_place(
                self.raw.as_ptr(),
                vec3(frame.origin),
                vec3(frame.x_axis),
                vec3(frame.y_axis),
            )
        };
        Error::check("SmApiTransform", status)
    }

    /// Moves the brep rigidly (`SmApiTransform`): the point p goes to
    /// `rotation · p + translation`, with `rotation` given row by row. It
    /// must be a rotation, orthonormal with determinant +1, to within 1e-9;
    /// otherwise this fails with `SHIM_ERR_BAD_FRAME` and leaves the brep
    /// alone.
    ///
    /// Every surface and curve keeps its type and parameterisation, so a
    /// plane stays a plane and a face keeps its (u, v) domain, and every face
    /// and edge keeps its index and [`Label`].
    pub fn transform(
        &mut self,
        rotation: [[f64; 3]; 3],
        translation: [f64; 3],
    ) -> Result<(), Error> {
        const OP: &str = "SmApiTransform";
        check_point(OP, translation)?;
        let rows = rotation.as_flattened();
        // SAFETY: self is a live brep on the kernel thread; rows holds the
        // nine doubles the shim reads.
        let status =
            unsafe { sys::sm_brep_transform(self.raw.as_ptr(), rows.as_ptr(), vec3(translation)) };
        Error::check(OP, status)
    }

    /// Moves the brep by `offset`: [`transform`](Brep::transform) with no
    /// rotation.
    pub fn translate(&mut self, offset: [f64; 3]) -> Result<(), Error> {
        self.transform(IDENTITY, offset)
    }

    /// Turns the brep by `angle_deg` degrees about the line through `point`
    /// along `axis`, counter-clockwise seen from the axis's tip:
    /// [`transform`](Brep::transform) with that rotation.
    pub fn rotate_about(
        &mut self,
        point: [f64; 3],
        axis: [f64; 3],
        angle_deg: f64,
    ) -> Result<(), Error> {
        const OP: &str = "SmApiTransform";
        check_point(OP, point)?;
        let rotation = rotation(axis, angle_deg).ok_or(Error::InvalidArgument {
            operation: OP,
            reason: "the axis must be finite and not zero, and the angle finite",
        })?;
        let turned = rotate(&rotation, point);
        self.transform(rotation, [0, 1, 2].map(|i| point[i] - turned[i]))
    }

    /// Scales the brep by `factor` about `point`, which stays put: the point
    /// p goes to `point + factor · (p - point)`. The factor must be finite
    /// and positive.
    ///
    /// A uniform scale keeps every surface and curve's type: planes,
    /// cylinders, cones, spheres and tori stay what they are. The volume
    /// scales by `factor³`. Every face and edge keeps its index and
    /// [`Label`]. The brep's tolerance, and every face's, edge's and
    /// vertex's, scale with it, so every gap stays within the tolerance it
    /// was within; other breps and the kernel's own tolerance are untouched.
    pub fn scale_about(&mut self, point: [f64; 3], factor: f64) -> Result<(), Error> {
        self.scale_axes_about(point, [factor; 3])
    }

    /// Scales the brep about `point` by a factor along each of x, y and z:
    /// the point p goes to `point + factors * (p - point)`, axis by axis.
    /// The factors must be finite and positive.
    ///
    /// Unequal factors turn every surface and curve into a NURBS first
    /// (`SmBrep::TurnToNURBS`), as SMLib does for any scale that would bend
    /// an analytic surface out of its kind: a cylinder stretched across its
    /// axis is no longer a cylinder. Equal factors are
    /// [`scale_about`](Brep::scale_about). Tolerances scale by the largest
    /// factor.
    pub fn scale_axes_about(&mut self, point: [f64; 3], factors: [f64; 3]) -> Result<(), Error> {
        const OP: &str = "SmBrep::Transform";
        check_point(OP, point)?;
        for factor in factors {
            check_positive(OP, factor, "scale factors must be finite and positive")?;
        }
        // SAFETY: self is a live brep on the kernel thread.
        let status = unsafe { sys::sm_brep_scale(self.raw.as_ptr(), vec3(point), vec3(factors)) };
        Error::check(OP, status)
    }

    /// The brep's tolerance in model units (`SmBrep::GetTolerance`): how far
    /// apart the kernel lets a vertex lie from its edges and faces, and an
    /// edge from its faces.
    pub fn tolerance(&self) -> Result<f64, Error> {
        let mut tolerance = 0.0;
        // SAFETY: self is a live brep on the kernel thread.
        let status = unsafe { sys::sm_brep_tolerance(self.raw.as_ptr(), &mut tolerance) };
        Error::check("SmBrep::GetTolerance", status)?;
        Ok(tolerance)
    }

    /// Sweeps the planar face at `face` in [`faces`](Brep::faces) along its
    /// outward normal by `distance` into a closed prism
    /// (`SmApiCreateCurveSweepFromFaces`). A negative distance sweeps into
    /// the solid. The brep is unchanged; combine the prism with it by
    /// [`union`](Brep::union) or [`difference`](Brep::difference).
    ///
    /// A side swept from a straight edge is a plane
    /// ([`Surface::Plane`](crate::Surface::Plane)) and one swept from a
    /// circular arc or circle a cylinder
    /// ([`Surface::Cylinder`](crate::Surface::Cylinder)); a side swept from
    /// any other curve is a NURBS surface.
    ///
    /// The prism's faces carry no labels. [`Sweep::faces`] says which is the
    /// start, which the end, and which edge of `face` each side was swept
    /// from, for labelling them.
    pub fn sweep_face(&self, face: usize, distance: f64) -> Result<Sweep, Error> {
        const OP: &str = "SmApiCreateCurveSweepFromFaces";
        if !distance.is_finite() || distance == 0.0 {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "sweep distance must be finite and non-zero",
            });
        }
        let mut out = ptr::null_mut();
        let mut faces = ptr::null_mut();
        let mut count = 0;
        // SAFETY: self is a live brep on the kernel thread; the outputs are
        // written only on success.
        let status = unsafe {
            sys::sm_brep_sweep_face(
                self.raw.as_ptr(),
                face,
                distance,
                &mut out,
                &mut faces,
                &mut count,
            )
        };
        Error::check(OP, status)?;
        // SAFETY: on success faces holds count entries.
        unsafe { Sweep::from_raw(OP, out, faces, count) }
    }

    /// [`sweep_face`](Brep::sweep_face) with draft: the sides lean
    /// `draft_deg` (under 90 either way) from the sweep's direction. A
    /// positive draft leans them into the face, so the end is the face shrunk
    /// all round by `distance` times the draft's tangent, its outer edge in
    /// and its holes out; a negative draft leans them out and grows it.
    /// Corners stay sharp. A draft of 0 is [`sweep_face`](Brep::sweep_face).
    ///
    /// Each loop of the face is swept by `SmPrimitiveCreation::CreateTaperExtrude`
    /// and the holes' solids cut from the outer loop's. A side swept from a
    /// straight edge is a plane and one from an arc or circle a cone
    /// ([`Surface::Cone`](crate::Surface::Cone)).
    ///
    /// [`Sweep::faces`] says what each face came from, as for
    /// [`sweep_face`](Brep::sweep_face). Where a line meets an arc at an
    /// angle, an outward draft can add a side standing on the corner, whose
    /// `source_edge` is `None`.
    ///
    /// Refuses a draft whose sides meet before the end of the sweep, where
    /// the face closes up or a hole leaning out swallows the end. The kernel
    /// fails an outward draft where a line meets an arc at an angle (a D
    /// shape).
    pub fn sweep_face_tapered(
        &self,
        face: usize,
        distance: f64,
        draft_deg: f64,
    ) -> Result<Sweep, Error> {
        const OP: &str = "SmPrimitiveCreation::CreateTaperExtrude";
        let invalid = |reason| {
            Err(Error::InvalidArgument {
                operation: OP,
                reason,
            })
        };
        if !distance.is_finite() || distance == 0.0 {
            return invalid("sweep distance must be finite and non-zero");
        }
        if !(draft_deg.is_finite() && draft_deg.abs() < 90.0) {
            return invalid("the draft angle must be finite and under 90 degrees either way");
        }
        let mut out = ptr::null_mut();
        let mut faces = ptr::null_mut();
        let mut count = 0;
        // SAFETY: self is a live brep on the kernel thread; the outputs are
        // written only on success.
        let status = unsafe {
            sys::sm_brep_sweep_face_tapered(
                self.raw.as_ptr(),
                face,
                distance,
                draft_deg,
                &mut out,
                &mut faces,
                &mut count,
            )
        };
        match status {
            sys::SM_SHIM_ERR_NOT_PLANAR => return invalid("only a planar face can be swept"),
            sys::SM_SHIM_ERR_DRAFT_CLOSES => {
                return invalid("the draft closes the face up before the end of the sweep");
            }
            _ => Error::check(OP, status)?,
        }
        // SAFETY: on success faces holds count entries.
        unsafe { Sweep::from_raw(OP, out, faces, count) }
    }

    /// Revolves the planar face at `face` in [`faces`](Brep::faces) about the
    /// axis through `axis_point` along `axis_direction` by `angle_deg`, into
    /// a closed solid (`SmTopologySweep` with `SmRotationalSweepGeometry`).
    /// It turns anticlockwise looking down the axis from its tip when the
    /// angle is positive, clockwise when negative, and at most a whole turn
    /// either way. The brep is unchanged; combine the solid with it by
    /// [`union`](Brep::union) or [`difference`](Brep::difference).
    ///
    /// The axis must lie in the face's plane. It may run along the face's
    /// boundary or touch it at a point, but not through its inside: such a
    /// solid would pass through itself. An edge lying on the axis sweeps
    /// nothing, so it has no side.
    ///
    /// A side turned from a straight edge is a plane (an edge square to the
    /// axis), a cylinder (along it) or a cone (at a slant), and one turned
    /// from an arc a sphere when the arc is centred on the axis and a torus
    /// otherwise.
    ///
    /// The solid's faces carry no labels. [`Sweep::faces`] says which is the
    /// start (the face itself), which the end (the face turned by the
    /// angle), and which edge of `face` each side was turned from. A whole
    /// turn has no start or end.
    pub fn revolve_face(
        &self,
        face: usize,
        axis_point: [f64; 3],
        axis_direction: [f64; 3],
        angle_deg: f64,
    ) -> Result<Sweep, Error> {
        const OP: &str = "SmTopologySweep (revolve)";
        let invalid = |reason| {
            Err(Error::InvalidArgument {
                operation: OP,
                reason,
            })
        };
        check_point(OP, axis_point)?;
        if !axis_direction.iter().all(|c| c.is_finite()) || axis_direction == [0.0; 3] {
            return invalid("the axis direction must be finite and not zero");
        }
        if !(angle_deg.is_finite() && angle_deg != 0.0 && angle_deg.abs() <= 360.0) {
            return invalid("the angle must be finite, not zero, and at most a whole turn");
        }
        let mut out = ptr::null_mut();
        let mut faces = ptr::null_mut();
        let mut count = 0;
        // SAFETY: self is a live brep on the kernel thread; the outputs are
        // written only on success.
        let status = unsafe {
            sys::sm_brep_revolve_face(
                self.raw.as_ptr(),
                face,
                vec3(axis_point),
                vec3(axis_direction),
                angle_deg,
                &mut out,
                &mut faces,
                &mut count,
            )
        };
        match status {
            sys::SM_SHIM_ERR_BAD_AXIS => return invalid("the axis must lie in the face's plane"),
            sys::SM_SHIM_ERR_AXIS_CROSSES_FACE => {
                return invalid("the axis runs through the face, which can't revolve across it");
            }
            sys::SM_SHIM_ERR_NOT_PLANAR => return invalid("only a planar face can be revolved"),
            _ => Error::check(OP, status)?,
        }
        // SAFETY: on success faces holds count entries.
        unsafe { Sweep::from_raw(OP, out, faces, count) }
    }

    /// Planar faces bounded by closed loops of lines, arcs and circles in
    /// `frame`'s xy plane (`SmApiCreatePlanarFaces`), for sweeping into
    /// prisms with [`sweep_face`](Brep::sweep_face).
    ///
    /// Each loop's curves run head to tail, the last back to the first; a
    /// circle is a loop on its own. Loops mustn't cross or touch each other.
    /// An outer loop bounds a face, a loop inside it a hole, and a loop inside
    /// a hole another face, however deep they nest. The result is a sheet, not
    /// a solid: one face per outer loop, with its
    /// [`normal`](FaceGeometry::normal) along the frame's z axis whichever way
    /// its loops run, so a positive sweep runs along z. Each curve's edge
    /// carries the curve's label, so each side of a prism swept from a face
    /// can be traced back to the curve it came from: through
    /// [`SweptFace::Side`](crate::SweptFace::Side)'s `source_edge` to the edge
    /// in this sheet's [`edges`](Brep::edges).
    ///
    /// The kernel's own error when the loops don't bound faces. It doesn't
    /// catch every bad arrangement: two squares that cross still make a face,
    /// so callers must hand over loops that don't.
    pub fn planar_faces(
        _kernel: &Kernel,
        frame: &Frame,
        loops: &[Vec<BoundaryCurve>],
    ) -> Result<Brep, Error> {
        const OP: &str = "SmApiCreatePlanarFaces";
        let invalid = |reason| {
            Err(Error::InvalidArgument {
                operation: OP,
                reason,
            })
        };
        if loops.is_empty() || loops.iter().any(Vec::is_empty) {
            return invalid("every loop needs at least one curve");
        }
        let mut raw = Vec::new();
        for curves in loops {
            if let Some(reason) = curves.iter().find_map(|c| c.curve.problem()) {
                return invalid(reason);
            }
            let circle = |c: &BoundaryCurve| matches!(c.curve, PlaneCurve::Circle { .. });
            if curves.len() > 1 && curves.iter().any(circle) {
                return invalid("a circle is a loop on its own");
            }
            let ends: Vec<([f64; 2], [f64; 2])> = curves
                .iter()
                .map(|c| (c.curve.start(), c.curve.end()))
                .collect();
            let next = ends.iter().cycle().skip(1);
            if ends
                .iter()
                .zip(next)
                .any(|((_, end), (start, _))| gap(*end, *start) > LOOP_GAP)
            {
                return invalid("each loop's curves must run head to tail and close");
            }
            for curve in curves {
                raw.push(curve.curve.to_raw(&curve.label, OP)?);
            }
        }
        let mut out = ptr::null_mut();
        // SAFETY: raw holds raw.len() curves; a fresh brep comes out.
        let status = unsafe {
            sys::sm_brep_create_planar_faces(
                raw.as_ptr(),
                raw.len(),
                vec3(frame.origin),
                vec3(frame.x_axis),
                vec3(frame.y_axis),
                &mut out,
            )
        };
        Error::check(OP, status)?;
        Brep::from_raw(OP, out)
    }

    /// A pipe: a circle of `radius` swept along `path`, lines and arcs in
    /// `frame`'s xy plane, hollow with `inner_radius` when that's given. Each
    /// curve of the path starts where the last ended and runs its own way
    /// (an arc's sweep may be negative); a circle is a closed path on its
    /// own. An open path is capped at both ends, a closed one isn't.
    ///
    /// Each curve gets a piece: a disc square to it at its start, swept
    /// along a line as [`sweep_face`](Brep::sweep_face) does or turned about
    /// an arc's centre as
    /// [`revolve_face`](Brep::revolve_face) does, so the walls along lines
    /// are cylinders and along arcs tori. Where curves meet tangentially the
    /// pieces meet face to face; where they meet at an angle the pipe is
    /// mitred, both pieces cut at the plane halfway between them. The pieces
    /// are joined (`SmApiBooleanUnion`), and a hollow pipe is the solid one
    /// less a solid pipe of `inner_radius`. [`Pipe::faces`] says which face
    /// is the outer wall, which the inner, which the start and end caps, and
    /// which path curve each wall was swept along: one outer wall (and one
    /// inner) per curve.
    ///
    /// Lines may meet at an angle, up to 179 degrees, but an arc must meet
    /// the curves either side of it tangentially: a plane can't mitre a
    /// cylinder against a torus. Refuses those paths, a bend tighter than the
    /// pipe (an arc whose radius isn't more than `radius`), and a hollow pipe
    /// along a closed path, whose bore would be sealed inside it (SMLib's
    /// tessellation winds some of such a cavity's walls the wrong way).
    pub fn pipe(
        _kernel: &Kernel,
        frame: &Frame,
        path: &[PlaneCurve],
        radius: f64,
        inner_radius: Option<f64>,
    ) -> Result<Pipe, Error> {
        const OP: &str = "sm_brep_create_pipe";
        let invalid = |reason| {
            Err(Error::InvalidArgument {
                operation: OP,
                reason,
            })
        };
        if path.is_empty() {
            return invalid("the path needs at least one curve");
        }
        if let Some(reason) = path.iter().find_map(PlaneCurve::problem) {
            return invalid(reason);
        }
        let circle = |c: &PlaneCurve| matches!(c, PlaneCurve::Circle { .. });
        if path.len() > 1 && path.iter().any(circle) {
            return invalid("a circle is a path on its own");
        }
        if path
            .windows(2)
            .any(|pair| gap(pair[0].end(), pair[1].start()) > LOOP_GAP)
        {
            return invalid("each curve of the path must start where the last ended");
        }
        // Where each curve meets the next, and the last the first when closed.
        let closed = path.len() > 1 && gap(path[path.len() - 1].end(), path[0].start()) <= LOOP_GAP;
        let ring = closed
            || path.iter().any(|c| match *c {
                PlaneCurve::Arc { sweep_deg, .. } => sweep_deg.abs() == 360.0,
                _ => circle(c),
            });
        if ring && inner_radius.is_some() {
            return invalid(
                "a pipe along a closed path can't be hollow: its bore would be sealed inside",
            );
        }
        let joints = path.windows(2).map(|pair| (pair[0], pair[1]));
        let closing = closed.then(|| (path[path.len() - 1], path[0]));
        for (before, after) in joints.chain(closing) {
            let ([x1, y1], [x2, y2]) = (before.end_tangent(), after.start_tangent());
            let turn = (x1 * y2 - y1 * x2).abs().atan2(x1 * x2 + y1 * y2);
            if turn > MAX_PIPE_TURN_DEG.to_radians() {
                return invalid("the path turns back on itself");
            }
            let arc = |c: PlaneCurve| !matches!(c, PlaneCurve::Line { .. });
            if turn > 1e-9 && (arc(before) || arc(after)) {
                return invalid(
                    "an arc in the path must meet the curves either side of it tangentially",
                );
            }
        }
        check_positive(OP, radius, "the radius must be finite and positive")?;
        let inner = inner_radius.unwrap_or(0.0);
        if !(inner.is_finite() && inner >= 0.0 && inner < radius) {
            return invalid(
                "the inner radius must be finite, not negative, and less than the radius",
            );
        }
        let tight = |c: &PlaneCurve| match *c {
            PlaneCurve::Arc { radius: bend, .. } | PlaneCurve::Circle { radius: bend, .. } => {
                bend <= radius
            }
            PlaneCurve::Line { .. } => false,
        };
        if path.iter().any(tight) {
            return invalid("a bend's radius must be more than the pipe's");
        }
        let raw = path
            .iter()
            .map(|curve| curve.to_raw(&None, OP))
            .collect::<Result<Vec<_>, _>>()?;
        let mut out = ptr::null_mut();
        let mut faces = ptr::null_mut();
        let mut count = 0;
        // SAFETY: raw holds raw.len() curves; the outputs are written only on
        // success.
        let status = unsafe {
            sys::sm_brep_create_pipe(
                raw.as_ptr(),
                raw.len(),
                vec3(frame.origin),
                vec3(frame.x_axis),
                vec3(frame.y_axis),
                radius,
                inner,
                &mut out,
                &mut faces,
                &mut count,
            )
        };
        match status {
            sys::SM_SHIM_ERR_BAD_PATH => return invalid("the path turns back on itself"),
            sys::SM_SHIM_ERR_ARC_CORNER => {
                return invalid(
                    "an arc in the path must meet the curves either side of it tangentially",
                );
            }
            _ => Error::check(OP, status)?,
        }
        // SAFETY: on success faces holds count entries, which we copy and then
        // free once.
        let roles = unsafe {
            let roles = slice(faces, count).iter().map(PipeFace::from_raw).collect();
            sys::sm_pipe_faces_free(faces);
            roles
        };
        Ok(Pipe {
            solid: Brep::from_raw(OP, out)?,
            faces: roles,
        })
    }

    /// The solid hollowed to a wall `thickness` thick, inward from its faces,
    /// with the faces at `open` in [`faces`](Brep::faces) left open
    /// (`SmApiShellBrepFull`: SMLib insets the solid, those faces left where
    /// they are, and cuts the inset from it). Corners inside meet sharp, as
    /// their faces extended meet. The brep is unchanged.
    ///
    /// [`Shell::faces`] says what each face is, by the face it came from: a
    /// face that stays outside, what's left of an opened face round its
    /// opening (the rim), or a face's offset, the inside of its wall. Faces
    /// that stay keep their labels, and so do the edges only they use;
    /// inner faces have none, for the caller to label after the faces they
    /// came from.
    ///
    /// At least one face must be open: a closed hollow is a sealed cavity,
    /// which a game mesh never shows, and SMLib's tessellation can wind its
    /// curved walls the wrong way.
    ///
    /// The kernel's error when it can't offset every face by `thickness`: a
    /// round, a chamfer or a boss can only be offset in by less than its own
    /// size. `SHIM_ERR_BAD_SHELL` when the result isn't a closed solid whose
    /// faces' integrals agree with their edges, as when the opened face is
    /// rounded into its neighbours.
    pub fn shell(&self, thickness: f64, open: &[usize]) -> Result<Shell, Error> {
        const OP: &str = "SmApiShellBrepFull";
        let invalid = |reason| {
            Err(Error::InvalidArgument {
                operation: OP,
                reason,
            })
        };
        check_positive(OP, thickness, "the thickness must be finite and positive")?;
        if !self.is_manifold()? {
            return invalid("only a closed, manifold solid can be hollowed");
        }
        let faces = self.face_count();
        let mut seen = vec![false; faces];
        for &face in open {
            if face >= faces || std::mem::replace(&mut seen[face], true) {
                return invalid("a face to open is out of range or repeated");
            }
        }
        if open.len() == faces {
            return invalid("opening every face leaves nothing");
        }
        if open.is_empty() {
            return invalid("open at least one face: a closed hollow is a sealed cavity");
        }
        let indices: Vec<u32> = open.iter().map(|&f| f as u32).collect();
        let mut out = ptr::null_mut();
        let mut roles = ptr::null_mut();
        let mut count = 0;
        // SAFETY: self is a live brep on the kernel thread; indices holds
        // indices.len() elements; the outputs are written only on success.
        let status = unsafe {
            sys::sm_brep_shell(
                self.raw.as_ptr(),
                thickness,
                indices.as_ptr(),
                indices.len(),
                &mut out,
                &mut roles,
                &mut count,
            )
        };
        Error::check(OP, status)?;
        // SAFETY: on success roles holds count entries, which we copy and then
        // free once.
        let faces = unsafe {
            let faces = slice(roles, count)
                .iter()
                .map(ShellFace::from_raw)
                .collect();
            sys::sm_shell_faces_free(roles);
            faces
        };
        Ok(Shell {
            brep: Brep::from_raw(OP, out)?,
            faces,
        })
    }

    /// Merges faces that lie on one surface into one face, and edges that lie
    /// on one curve into one edge (`SmBrep::RemoveTopologicalEdgesAndVertices`),
    /// without changing the shape: after a union that leaves a wall in two
    /// coplanar pieces, it's one face again. Merged faces keep the older label
    /// (see [`Label`]). Booleans already do this along their intersection
    /// curves.
    pub fn merge_faces(&mut self) -> Result<(), Error> {
        // SAFETY: self is a live brep on the kernel thread.
        let status = unsafe { sys::sm_brep_merge_faces(self.raw.as_ptr()) };
        Error::check("SmBrep::RemoveTopologicalEdgesAndVertices", status)
    }

    /// [`merge_faces`](Brep::merge_faces), except that the edges at `keep` in
    /// [`edges`](Brep::edges) stay, whether or not the faces either side lie
    /// on one surface (`SmBrep::RemoveTopologicalEdgesAndVertices` with its
    /// keep list). Keeping every edge but some merges faces only across
    /// those, as when taking an [`imprint`](Brep::imprint) back off some
    /// faces. Vertices left in the middle of one curve go, wherever they are.
    pub fn merge_faces_keeping(&mut self, keep: &[usize]) -> Result<(), Error> {
        const OP: &str = "SmBrep::RemoveTopologicalEdgesAndVertices";
        let indices = keep
            .iter()
            .map(|&edge| u32::try_from(edge))
            .collect::<Result<Vec<u32>, _>>()
            .map_err(|_| Error::InvalidArgument {
                operation: OP,
                reason: "an edge index is out of range",
            })?;
        // SAFETY: self is a live brep on the kernel thread; indices holds
        // indices.len() elements.
        let status = unsafe {
            sys::sm_brep_merge_faces_keeping(self.raw.as_ptr(), indices.as_ptr(), indices.len())
        };
        Error::check(OP, status)
    }

    /// Imprints `tool` on this brep: splits its faces and edges along the
    /// curves where `tool`'s boundary meets them, which become edges and
    /// vertices, and changes nothing else (`SmMerge::ManifoldBoolean` with
    /// `SM_BO_IMPRINT`). Unlike the booleans, it doesn't merge faces on one
    /// surface afterwards, so the pieces of a split face stay apart. The
    /// shape and volume are unchanged.
    ///
    /// The pieces of a split face or edge keep its label (see [`Label`]);
    /// the new edges have none. [`Imprint::edges`] lists every edge lying
    /// along the tool, with the labels of the tool faces it lies on.
    ///
    /// Consumes both, as [`union`](Brep::union).
    pub fn imprint(self, tool: Brep) -> Result<Imprint, Error> {
        self.imprint_with(tool, false)
    }

    /// [`imprint`](Brep::imprint), then sorts the faces as an intersection
    /// with `tool` would, without deleting any (`SM_BO_INTERSECTION` with
    /// `bImprintAndClassifyFaces`): [`Imprint::inside`] says which faces lie
    /// inside the tool.
    pub fn imprint_classified(self, tool: Brep) -> Result<Imprint, Error> {
        self.imprint_with(tool, true)
    }

    fn imprint_with(self, tool: Brep, classify: bool) -> Result<Imprint, Error> {
        const OP: &str = "SmMerge::ManifoldBoolean (imprint)";
        let a = ManuallyDrop::new(self);
        let b = ManuallyDrop::new(tool);
        let mut out = ptr::null_mut();
        let mut edges = ptr::null_mut();
        let mut edge_count = 0;
        let mut inside = ptr::null_mut();
        let mut face_count = 0;
        // SAFETY: a and b are distinct live breps (each owned by value) on the
        // kernel thread. The shim takes ownership of both for the call; the
        // outputs are written only on success.
        let status = unsafe {
            sys::sm_brep_imprint(
                a.raw.as_ptr(),
                b.raw.as_ptr(),
                i32::from(classify),
                &mut out,
                &mut edges,
                &mut edge_count,
                &mut inside,
                &mut face_count,
            )
        };
        if status != sys::SM_STATUS_SUCCESS {
            // Both are still allocated but may be half-modified: discard them.
            drop(ManuallyDrop::into_inner(a));
            drop(ManuallyDrop::into_inner(b));
            return Err(Error::Kernel {
                operation: OP,
                status: crate::Status::from_raw(status),
            });
        }
        debug_assert_eq!(out, a.raw.as_ptr(), "{OP} should return its first operand");
        // SAFETY: on success edges holds edge_count entries and inside
        // face_count, which we copy and then free once.
        let (edges, inside) = unsafe {
            let copied = (
                slice(edges, edge_count)
                    .iter()
                    .map(ImprintEdge::from_raw)
                    .collect(),
                slice(inside, face_count).iter().map(|&b| b != 0).collect(),
            );
            sys::sm_imprint_free(edges, inside);
            copied
        };
        Ok(Imprint {
            brep: Brep::from_raw(OP, out)?,
            edges,
            inside,
        })
    }

    /// The indices in [`faces`](Brep::faces) of the faces `query` picks.
    pub fn find_faces(
        &self,
        query: &FaceQuery,
        tolerance: &Tolerance,
    ) -> Result<Vec<usize>, Error> {
        if let Some(reason) = query.problem() {
            return Err(Error::InvalidArgument {
                operation: "FaceQuery",
                reason,
            });
        }
        Ok(self
            .faces()?
            .iter()
            .filter(|face| query.matches(face, tolerance))
            .map(|face| face.index)
            .collect())
    }

    /// Moves the faces at `faces` (indices in [`faces`](Brep::faces)) the
    /// signed `distance` along their outward normals, as Parasolid's face
    /// offset does: each other neighbour is extended or trimmed along its own
    /// surface to meet the moved faces, so a chamfer grows or shrinks and a
    /// drafted wall keeps its draft. A picked plane is translated, and a
    /// picked cone slides along its axis (its offset is that translation);
    /// other curved faces (cylinders, spheres, tori, NURBS) take their offset
    /// surface, so a hole or a boss changes radius. Blends on a picked plane
    /// or cone (fillets along its edges, the corners between them) are
    /// translated by the shortest move that keeps them tangent to all their
    /// neighbours: on a raised top they rise and the walls stretch, against a
    /// drafted wall they slide along it, between two picked faces they follow
    /// both.
    ///
    /// A picked plane pushed in can run into features standing on the faces
    /// beside it: a boss or a hole on a lid near the edge above a pushed
    /// wall, or under a corner round that moves with the wall, whose loops
    /// would cross the lid's new boundary. A picked plane can also swallow a
    /// face of a feature on it whole: a pocket pushed past its floor, a
    /// chamfer past its width, a boss pulled out past its top. A picked
    /// curved face moved either way crosses features beside it the same
    /// way: a boss near a hole that widens, a hole or boss beside a boss
    /// that grows, a pin on a post that thins. Each such feature is set
    /// aside, the faces are offset without it, and it comes back cut by the
    /// face's new position, as Parasolid pushes a face through: a boss
    /// joined back as far as the face leaves it, or not at all; a hole as
    /// the part with its holes less what the offset swept away, or with what
    /// it grew. [`offset_obstacles`](Brep::offset_obstacles) lists them.
    /// Where a cut boss's new face doesn't merge into the moved face it is
    /// unlabelled; everything else keeps its [`Label`], though a cut face
    /// can come back in two pieces along its surface's seam. Only the faces
    /// the push-through's own cuts split merge back: an edge that already
    /// joined two faces on one surface, as an [`imprint`](Brep::imprint)
    /// leaves them, stays.
    /// `SHIM_ERR_OBSTRUCTED` when one can't be set aside: it doesn't stand
    /// entirely on one side of a planar face, or something else stands in
    /// its way; or when it would come back on its own or sealed inside.
    ///
    /// A face pushed in past a chamfer beside it swallows the chamfer and
    /// meets the plane beyond directly, as if the edge had never been
    /// chamfered: a lid past the chamfer along its edge, a hole widened or a
    /// post thinned past the chamfer on its rim. That needs the face to run
    /// along the plane's normal, and nothing but the chamfer and what meets it
    /// near where the face passes.
    ///
    /// A picked face beside an unpicked face on its own surface, as the two
    /// pieces of a face an [`imprint`](Brep::imprint) splits are, can't
    /// stretch that neighbour to meet it, so it steps away from it: the
    /// picked faces are swept along the move into prisms, joined when pulled
    /// out and cut when pushed in, which stands a wall along each edge to
    /// such a face while every other neighbour stretches or is trimmed along
    /// its surface as before. A panel imprinted on a lid rises or sinks. The
    /// moved faces keep their [`Label`]s, and each wall takes the label of
    /// the edge it stands on, if that has one, so an imprinted edge's label
    /// tells which wall stands on it. That needs
    /// every picked face planar, all facing one way, and every other
    /// neighbour standing along the move (a plane containing it, or a
    /// cylinder along it); otherwise `SHIM_ERR_UNSUPPORTED_STEP`. A step never
    /// pushes through anything: `SHIM_ERR_OBSTRUCTED` when the prisms don't
    /// add or take exactly their own volume, as when a face pushed in runs
    /// into a pocket under it.
    ///
    /// With nothing in the way the topology doesn't change, so every face and
    /// edge keeps its [`Label`]: the moved faces and the neighbours that
    /// stretch alike. A move that would otherwise shrink a face to nothing
    /// (the walls under a lid's rounds, a whole box) fails with
    /// `SHIM_ERR_BROKE_SOLID` or `SHIM_ERR_CHANGED_TOPOLOGY`, never a broken
    /// solid: the kernel's result is checked for manifoldness, the same faces
    /// and edges, every planar face's loops apart (none crossing, every inner
    /// one inside the outer), gaps within the brep's tolerance, and a volume that can be
    /// trusted: each face's integral agrees with its edges, or else the
    /// volume agrees with the tessellation's. `SHIM_ERR_UNSUPPORTED_FACE`: a
    /// picked curved face other than a cone has blends on it (their new
    /// shape isn't a translation), no translation keeps a blend tangent to
    /// all its neighbours, or a surface has no offset. Some moves the kernel
    /// still gets wrong, and refuses (`SHIM_ERR_BROKE_SOLID`): fillets
    /// sliding along drafted walls, and some moves of faces cut by spheres
    /// (see `PATCHES.md`).
    ///
    /// Built on `SmBrep::LocalOperation`. `faces` must be non-empty and
    /// distinct (`SHIM_ERR_BAD_FACE` otherwise); a zero distance leaves the
    /// brep as it is. Consumed and dropped on failure, as
    /// [`fillet`](Brep::fillet).
    pub fn offset_faces(self, faces: &[usize], distance: f64) -> Result<Brep, Error> {
        const OP: &str = "SmBrep::LocalOperation";
        let indices = offset_indices(OP, faces, distance)?;
        // SAFETY: self is a live brep on the kernel thread; indices holds
        // indices.len() values.
        let status = unsafe {
            sys::sm_brep_offset_faces(self.raw.as_ptr(), indices.as_ptr(), indices.len(), distance)
        };
        Error::check(OP, status)?;
        Ok(self)
    }

    /// [`offset_faces`](Brep::offset_faces)'s local operation alone: the
    /// kernel's own result, neither checked nor pushed through anything in
    /// the way, so it can be a broken solid that fails later, in
    /// [`tessellate`](Brep::tessellate) for one. For diagnosing the kernel,
    /// and for testing what callers do with a broken solid.
    pub fn offset_faces_unchecked(self, faces: &[usize], distance: f64) -> Result<Brep, Error> {
        const OP: &str = "SmBrep::LocalOperation";
        let indices = offset_indices(OP, faces, distance)?;
        // SAFETY: as in offset_faces.
        let status = unsafe {
            sys::sm_brep_offset_faces_unchecked(
                self.raw.as_ptr(),
                indices.as_ptr(),
                indices.len(),
                distance,
            )
        };
        Error::check(OP, status)?;
        Ok(self)
    }

    /// What offsetting `faces` by `distance` runs into: the features standing
    /// on faces beside the moving ones, along inner loops, that reach into
    /// the space a picked plane pushed in, the blends moving with it, or a
    /// picked curved face moved either way sweep, and those on a picked
    /// plane it swallows a face of (see
    /// [`offset_faces`](Brep::offset_faces)). Each is its faces' indices, in
    /// face order, and the features come in the order of their first faces;
    /// none when nothing is in the way.
    pub fn offset_obstacles(
        &self,
        faces: &[usize],
        distance: f64,
    ) -> Result<Vec<Vec<usize>>, Error> {
        const OP: &str = "offset obstacles";
        let indices = offset_indices(OP, faces, distance)?;
        let query = |out_faces: *mut u32, out_obstacles: *mut u32, capacity: usize| {
            let mut count = 0;
            // SAFETY: self is a live brep on the kernel thread; indices holds
            // indices.len() values; the outputs hold `capacity` values each.
            let status = unsafe {
                sys::sm_brep_offset_obstacles(
                    self.raw.as_ptr(),
                    indices.as_ptr(),
                    indices.len(),
                    distance,
                    out_faces,
                    out_obstacles,
                    capacity,
                    &mut count,
                )
            };
            Error::check(OP, status).map(|()| count)
        };
        let count = query(ptr::null_mut(), ptr::null_mut(), 0)?;
        let mut face_list = vec![0u32; count];
        let mut obstacle_list = vec![0u32; count];
        query(face_list.as_mut_ptr(), obstacle_list.as_mut_ptr(), count)?;
        let mut obstacles: Vec<Vec<usize>> = Vec::new();
        for (face, obstacle) in face_list.into_iter().zip(obstacle_list) {
            let obstacle = obstacle as usize;
            if obstacles.len() <= obstacle {
                obstacles.resize_with(obstacle + 1, Vec::new);
            }
            obstacles[obstacle].push(face as usize);
        }
        Ok(obstacles)
    }

    /// Moves the faces at `faces` (indices in [`faces`](Brep::faces)) by a
    /// rigid `motion`: turned about an axis, if it turns, then translated.
    /// Every other neighbour is extended or trimmed along its own surface to
    /// meet them, as in [`offset_faces`](Brep::offset_faces), on the same
    /// local operation (`SmBrep::LocalOperation`), so a moved plane can lean
    /// at any angle and a hole or a boss can move sideways.
    ///
    /// A planar face moved within its own plane stays where it is (its
    /// boundary is its neighbours'), so only the move along its normal shows;
    /// move a feature's faces together to slide the feature. A translation
    /// carries the blends on the moved faces (fillets, the corners between
    /// them): each slides the shortest way that keeps it tangent to every
    /// neighbour, as an offset's do, and one the move leaves in place
    /// stretches. A turn can't carry a blend: a round on a turned face is
    /// refused.
    ///
    /// The faces and edges stay the same, in the same order, each with its
    /// [`Label`]. [`Moved::faces`](crate::Moved::faces) says what happened
    /// to each face: moved, carried, followed (extended or trimmed) or kept.
    ///
    /// Refused, in plain words ([`Error::InvalidArgument`]), when the faces
    /// round the moved ones can't reach them: one would shrink to nothing or
    /// cross another, the result wouldn't be a closed solid, a round can't
    /// follow, a face beside them on their own surface (as an
    /// [`imprint`](Brep::imprint) leaves) would have to stretch, or the
    /// kernel can't intersect a neighbour with them. The result is checked
    /// as an offset's is: manifold, no gaps, every planar face's loops
    /// apart, and a volume that can be trusted. Nothing in the way is cut
    /// through. Consumed and dropped on failure, as
    /// [`offset_faces`](Brep::offset_faces). A motion that moves nothing
    /// gives the brep back as it was.
    pub fn move_faces(self, faces: &[usize], motion: &RigidMotion) -> Result<Moved, Error> {
        const OP: &str = "SmBrep::LocalOperation (move)";
        let invalid = |reason| {
            Err(Error::InvalidArgument {
                operation: OP,
                reason,
            })
        };
        if faces.is_empty() {
            return invalid("no faces to move");
        }
        if !motion.translation.iter().all(|c| c.is_finite()) {
            return invalid("the translation must be finite");
        }
        let (point, axis, degrees) = match motion.turn {
            Some(turn) => {
                check_point(OP, turn.point)?;
                if !turn.degrees.is_finite() {
                    return invalid("the angle must be finite");
                }
                if turn.degrees != 0.0
                    && (!turn.axis.iter().all(|c| c.is_finite()) || turn.axis == [0.0; 3])
                {
                    return invalid("the axis direction must be finite and not zero");
                }
                (turn.point, turn.axis, turn.degrees)
            }
            None => ([0.0; 3], [0.0, 0.0, 1.0], 0.0),
        };
        let indices: Vec<u32> = faces
            .iter()
            .map(|&face| u32::try_from(face).unwrap_or(u32::MAX))
            .collect();
        let mut roles = vec![0u8; self.face_count()];
        // SAFETY: self is a live brep on the kernel thread; indices holds
        // indices.len() values and roles roles.len().
        let status = unsafe {
            sys::sm_brep_move_faces(
                self.raw.as_ptr(),
                indices.as_ptr(),
                indices.len(),
                vec3(motion.translation),
                vec3(point),
                vec3(axis),
                degrees,
                roles.as_mut_ptr(),
                roles.len(),
            )
        };
        match status {
            sys::SM_STATUS_SUCCESS => {}
            sys::SM_SHIM_ERR_BAD_FACE => {
                return invalid("a face to move is out of range or repeated");
            }
            sys::SM_SHIM_ERR_BAD_AXIS => {
                return invalid("the axis direction must be finite and not zero");
            }
            sys::SM_SHIM_ERR_CHANGED_TOPOLOGY => {
                return invalid(
                    "the faces round them can't reach them: one would shrink to nothing or cross another",
                );
            }
            sys::SM_SHIM_ERR_BROKE_SOLID => {
                return invalid("the faces round them can't reach them and still close the solid");
            }
            sys::SM_SHIM_ERR_UNSUPPORTED_FACE if degrees != 0.0 => {
                return invalid(
                    "the faces round them can't reach them: a round on a turned face can't turn with it",
                );
            }
            sys::SM_SHIM_ERR_UNSUPPORTED_FACE => {
                return invalid(
                    "the faces round them can't reach them: no slide keeps a round on them tangent",
                );
            }
            sys::SM_SHIM_ERR_UNSUPPORTED_STEP => {
                return invalid(
                    "the faces round them can't reach them: a face beside them on their own surface can't stretch",
                );
            }
            // The kernel's own failures: the local operation found no way to
            // meet a neighbour with the moved faces.
            #[allow(clippy::useless_conversion)] // sm_status is a C long, i32 on Windows
            status if (1001..2000).contains(&i64::from(status)) => {
                return invalid("the faces round them can't reach them");
            }
            status => Error::check(OP, status)?,
        }
        Ok(Moved {
            brep: self,
            faces: roles.into_iter().map(MovedFace::from_raw).collect(),
        })
    }

    /// How many edges `query` picks.
    pub fn count_edges(&self, query: &EdgeQuery, tolerance: &Tolerance) -> Result<usize, Error> {
        check_query("EdgeQuery", query)?;
        Ok(self
            .edges()?
            .iter()
            .filter(|edge| query.matches(edge, tolerance))
            .count())
    }

    /// Fillets edges picked by geometric queries
    /// (`SmApiFilletEdgesPerEdge`), each group with its own radius and
    /// cross-section.
    ///
    /// Every query must pick at least one edge. An edge picked by more than
    /// one blend takes the first. A manifold solid that the kernel leaves
    /// non-manifold is an error (`SHIM_ERR_BROKE_SOLID`), even when the kernel
    /// reports success.
    ///
    /// The brep is consumed: on failure the kernel may have half-modified it,
    /// so it is dropped. There is no undo; the kernel's `SmApiRemoveFillet`
    /// does nothing (KNOWN_ISSUES.md), so rebuild from the unfilleted solid
    /// instead.
    pub fn fillet(self, blends: &[EdgeBlend], options: &FilletOptions) -> Result<Brep, Error> {
        const OP: &str = "SmApiFilletEdgesPerEdge";
        if blends.is_empty() {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "no edges to fillet",
            });
        }
        for blend in blends {
            check_query(OP, &blend.edges)?;
            check_positive(
                OP,
                blend.radius,
                "fillet radius must be finite and positive",
            )?;
        }
        check_positive(
            OP,
            options.thumbweight,
            "thumbweight must be finite and positive",
        )?;
        let picked = self.pick(OP, blends.iter().map(|b| &b.edges), &options.tolerance)?;
        let indices: Vec<u32> = picked.iter().map(|&(edge, _)| edge).collect();
        let radii: Vec<f64> = picked.iter().map(|&(_, b)| blends[b].radius).collect();
        let sections: Vec<i32> = picked
            .iter()
            .map(|&(_, b)| blends[b].section.to_raw())
            .collect();
        // SAFETY: self is a live brep on the kernel thread; the three arrays
        // each hold indices.len() elements.
        let status = unsafe {
            sys::sm_brep_fillet_edges(
                self.raw.as_ptr(),
                indices.as_ptr(),
                radii.as_ptr(),
                sections.as_ptr(),
                indices.len(),
                options.continuity as u32,
                options.thumbweight,
            )
        };
        Error::check(OP, status)?;
        Ok(self)
    }

    /// Chamfers the edges `query` picks with a flat face `width` across, as
    /// `SmApiChamferFillet` does for every edge: a 90° edge is set back by
    /// `width / √2` on each face, and corners where chamfered edges meet are
    /// beveled.
    ///
    /// Unlike `SmApiChamferFillet`, this leaves the kernel's global merge off;
    /// with it on, chamfers on circular edges can come back as open shells.
    /// A chamfer that can't fit, such as one that would cut past a whole
    /// face, is an error. For a chamfer given by its setback on each face,
    /// use [`fillet`](Brep::fillet) with
    /// [`EdgeBlend::linear`](crate::EdgeBlend::linear).
    ///
    /// Consumed and dropped on failure, as [`fillet`](Brep::fillet).
    pub fn chamfer(
        self,
        query: &EdgeQuery,
        width: f64,
        tolerance: &Tolerance,
    ) -> Result<Brep, Error> {
        const OP: &str = "SmApiChamferFillet";
        check_query(OP, query)?;
        check_positive(OP, width, "chamfer width must be finite and positive")?;
        let picked = self.pick(OP, [query], tolerance)?;
        let indices: Vec<u32> = picked.iter().map(|&(edge, _)| edge).collect();
        // SAFETY: self is a live brep on the kernel thread; indices holds
        // indices.len() elements.
        let status = unsafe {
            sys::sm_brep_chamfer_edges(self.raw.as_ptr(), indices.as_ptr(), indices.len(), width)
        };
        Error::check(OP, status)?;
        Ok(self)
    }

    /// Resolves queries to `(edge index, query index)` pairs, each edge going
    /// to the first query that picks it. The indices are only good until the
    /// brep changes, so they never leave this crate.
    fn pick<'q>(
        &self,
        operation: &'static str,
        queries: impl IntoIterator<Item = &'q EdgeQuery>,
        tolerance: &Tolerance,
    ) -> Result<Vec<(u32, usize)>, Error> {
        let edges = self.edges()?;
        let mut owner: Vec<Option<usize>> = vec![None; edges.len()];
        for (q, query) in queries.into_iter().enumerate() {
            let mut any = false;
            for (edge, slot) in edges.iter().zip(&mut owner) {
                if query.matches(edge, tolerance) {
                    any = true;
                    slot.get_or_insert(q);
                }
            }
            if !any {
                return Err(Error::InvalidArgument {
                    operation,
                    reason: "an edge query picks no edges",
                });
            }
        }
        Ok(owner
            .iter()
            .enumerate()
            .filter_map(|(edge, q)| q.map(|q| (edge as u32, q)))
            .collect())
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

    /// Starts tessellating the solid one face at a time: splits every edge
    /// as [`tessellate`](Brep::tessellate) would, on a copy, and leaves the
    /// faces to [`Tessellation::mesh`].
    pub fn begin_tessellation(
        &self,
        params: &TessellationParams,
    ) -> Result<Tessellation<'_>, Error> {
        let raw_params = params.to_raw()?;
        let mut out = ptr::null_mut();
        // SAFETY: self is a live brep on the kernel thread; params outlives
        // the call. The tessellation borrows self, so the brep outlives it
        // and, being borrowed, can't change before its mesh is made.
        let status =
            unsafe { sys::sm_brep_tessellation_begin(self.raw.as_ptr(), &raw_params, &mut out) };
        Error::check("SmApiTessellate", status)?;
        let raw = NonNull::new(out).ok_or(Error::Kernel {
            operation: "SmApiTessellate",
            status: crate::Status::from_raw(sys::SM_SHIM_ERR_NULL_ARGUMENT),
        })?;
        // SAFETY: raw is a tessellation the shim just made of this brep.
        Ok(unsafe { Tessellation::from_raw(raw) })
    }

    /// [`tessellate`](Brep::tessellate) one face at a time, triangulating
    /// only the faces `wanted` marks: see [`Tessellation::mesh`].
    pub fn tessellate_faces(
        &self,
        params: &TessellationParams,
        wanted: Option<&[bool]>,
    ) -> Result<Mesh, Error> {
        self.begin_tessellation(params)?.mesh(wanted)
    }

    /// Each face's key for [`tessellate_some`](Brep::tessellate_some), in
    /// face order, taken before tessellating: a 128-bit digest of everything
    /// SMLib's tessellation reads for the face when no face's boundary needs
    /// refining. That is the parameters, the brep's tolerance, the face's
    /// surface, domain, tolerance and loops (each edge's curve, interval and
    /// trim curves, and their vertices), and for each of its edges, the
    /// surface, domain, tolerance and order of the other faces on it. `None`
    /// for a face SMLib can't describe that way.
    ///
    /// These keys never equal [`Tessellation::face_keys`]'.
    pub fn tessellation_keys(
        &self,
        params: &TessellationParams,
    ) -> Result<Vec<Option<FaceKey>>, Error> {
        let raw_params = params.to_raw()?;
        let count = self.face_count();
        let mut keys = vec![[0u64; 2]; count];
        let mut keyed = vec![0i32; count];
        // SAFETY: self is a live brep on the kernel thread; keys holds two
        // words and keyed one flag for each of its count faces.
        let status = unsafe {
            sys::sm_brep_tessellation_keys(
                self.raw.as_ptr(),
                &raw_params,
                keys.as_mut_ptr().cast(),
                keyed.as_mut_ptr(),
                count,
            )
        };
        Error::check("SmApiTessellate", status)?;
        Ok(keys
            .into_iter()
            .zip(keyed)
            .map(|(key, keyed)| (keyed != 0).then_some(FaceKey(key)))
            .collect())
    }

    /// Triangulates the faces `wanted` marks, one flag per face, into a mesh
    /// laid out as [`Tessellation::mesh`] lays it out, running SMLib's edge
    /// splitting only for those faces and the faces sharing their edges.
    ///
    /// Each face's triangles are [`tessellate`](Brep::tessellate)'s for it
    /// provided every face not wanted has the
    /// [`tessellation_keys`](Brep::tessellation_keys) key of a face an
    /// earlier `tessellate_some` triangulated: then no face's boundary needs
    /// refining, which this checks only for the faces wanted, and a face's
    /// key fixes its triangles. `None` when this can't vouch for its faces:
    /// a wanted face's boundary needs refining, or something SMLib's whole
    /// tessellation would recover from went wrong. Then
    /// [`begin_tessellation`](Brep::begin_tessellation) tessellates the
    /// whole solid.
    ///
    /// With `caches`, faces take surface caches kept from earlier calls for
    /// identical surfaces (wanted faces take copies), and the faces sharing
    /// wanted faces' edges leave theirs there, as do wanted faces whose
    /// surfaces recur: see [`SurfaceCaches`]. The mesh comes out the same.
    pub fn tessellate_some(
        &self,
        params: &TessellationParams,
        wanted: &[bool],
        caches: Option<&mut SurfaceCaches>,
    ) -> Result<Option<Mesh>, Error> {
        let raw_params = params.to_raw()?;
        let flags: Vec<u8> = wanted.iter().map(|&w| u8::from(w)).collect();
        let mut out = ptr::null_mut();
        let caches = caches.map_or(ptr::null_mut(), SurfaceCaches::as_ptr);
        // SAFETY: self is a live brep on the kernel thread; params and the
        // flags outlive the call, and caches is null or a live store on this
        // thread.
        let status = unsafe {
            sys::sm_brep_tessellate_some(
                self.raw.as_ptr(),
                &raw_params,
                flags.as_ptr(),
                flags.len(),
                caches,
                &mut out,
            )
        };
        Error::check("SmApiTessellate", status)?;
        let Some(raw) = NonNull::new(out) else {
            return Ok(None);
        };
        // SAFETY: raw is a mesh the shim just allocated; we free it once,
        // after copying its arrays out.
        unsafe {
            let mesh = Mesh::from_raw(raw.as_ptr());
            sys::sm_mesh_free(raw.as_ptr());
            Ok(Some(mesh))
        }
    }
}

/// `faces` as the shim's indices, after the checks every offset makes.
fn offset_indices(
    operation: &'static str,
    faces: &[usize],
    distance: f64,
) -> Result<Vec<u32>, Error> {
    if !distance.is_finite() {
        return Err(Error::InvalidArgument {
            operation,
            reason: "offset distance must be finite",
        });
    }
    if faces.is_empty() {
        return Err(Error::InvalidArgument {
            operation,
            reason: "no faces to offset",
        });
    }
    Ok(faces
        .iter()
        .map(|&face| u32::try_from(face).unwrap_or(u32::MAX))
        .collect())
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

fn gap(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1])
}

fn vec3([x, y, z]: [f64; 3]) -> sys::sm_vec3 {
    sys::sm_vec3 { x, y, z }
}

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// The rotation by `angle_deg` degrees about `axis`, counter-clockwise seen
/// from its tip (Rodrigues' formula), row by row. `None` for a zero or
/// non-finite axis or angle.
fn rotation(axis: [f64; 3], angle_deg: f64) -> Option<[[f64; 3]; 3]> {
    let length = axis.iter().map(|c| c * c).sum::<f64>().sqrt();
    if !length.is_finite() || length == 0.0 || !angle_deg.is_finite() {
        return None;
    }
    let [x, y, z] = axis.map(|c| c / length);
    let (sin, cos) = angle_deg.to_radians().sin_cos();
    let t = 1.0 - cos;
    Some([
        [cos + x * x * t, x * y * t - z * sin, x * z * t + y * sin],
        [y * x * t + z * sin, cos + y * y * t, y * z * t - x * sin],
        [z * x * t - y * sin, z * y * t + x * sin, cos + z * z * t],
    ])
}

fn rotate(rotation: &[[f64; 3]; 3], point: [f64; 3]) -> [f64; 3] {
    rotation.map(|row| row[0] * point[0] + row[1] * point[1] + row[2] * point[2])
}

const EMPTY_EDGE: sys::sm_edge_geometry = sys::sm_edge_geometry {
    points: [[0.0; 3]; EDGE_SAMPLES],
    face_a: -1,
    face_b: -1,
    use_count: 0,
    normal_a: [0.0; 3],
    normal_b: [0.0; 3],
    label: NO_LABEL,
    label_a: NO_LABEL,
    label_b: NO_LABEL,
};

const EMPTY_FACE: sys::sm_face_geometry = sys::sm_face_geometry {
    info: sys::sm_face_info {
        kind: 0,
        kernel_type: 0,
        origin: [0.0; 3],
        axis: [0.0; 3],
        ref_dir: [0.0; 3],
        radius: 0.0,
        minor_radius: 0.0,
        half_angle_rad: 0.0,
        degree_u: 0,
        degree_v: 0,
        poles_u: 0,
        poles_v: 0,
        rational: 0,
        uv_min: [0.0; 2],
        uv_max: [0.0; 2],
        label: NO_LABEL,
    },
    bounds_min: [0.0; 3],
    bounds_max: [0.0; 3],
    planar: 0,
    normal: [0.0; 3],
};

fn check_query(operation: &'static str, query: &EdgeQuery) -> Result<(), Error> {
    match query.problem() {
        Some(reason) => Err(Error::InvalidArgument { operation, reason }),
        None => Ok(()),
    }
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
