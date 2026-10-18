use std::fmt;
use std::marker::PhantomData;
use std::mem::ManuallyDrop;
use std::ptr::{self, NonNull};

use smlib_sys as sys;

use crate::edges::{EDGE_SAMPLES, EdgeGeometry, EdgeQuery, Tolerance};
use crate::error::Error;
use crate::faces::{FaceCentroid, FaceGeometry, FaceQuery, Frame, Sweep, SweptFace};
use crate::fillet::{EdgeBlend, FilletOptions};
use crate::kernel::Kernel;
use crate::label::{Label, NO_LABEL};
use crate::mesh::{Mesh, TessellationParams};
use crate::tessellation::{FaceKey, SurfaceCaches, Tessellation};

/// Relative accuracy of [`Brep::face_centroid`]: placing a feature at a
/// centroid is then exact to a tenth of a micron on a metre-sized part.
const CENTROID_ACCURACY: f64 = 1e-7;

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

    /// Sweeps the planar face at `face` in [`faces`](Brep::faces) along its
    /// outward normal by `distance` into a closed prism
    /// (`SmApiCreateCurveSweepFromFaces`). A negative distance sweeps into
    /// the solid. The brep is unchanged; combine the prism with it by
    /// [`union`](Brep::union) or [`difference`](Brep::difference).
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
        // SAFETY: on success faces holds count entries, which we copy and then
        // free once.
        let roles = unsafe {
            let roles = if faces.is_null() {
                Vec::new()
            } else {
                std::slice::from_raw_parts(faces, count)
                    .iter()
                    .map(SweptFace::from_raw)
                    .collect()
            };
            sys::sm_swept_faces_free(faces);
            roles
        };
        Ok(Sweep {
            prism: Brep::from_raw(OP, out)?,
            faces: roles,
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
    /// The topology doesn't change, so every face and edge keeps its
    /// [`Label`]: the moved faces and the neighbours that stretch alike. A
    /// move that would shrink a face to nothing fails with
    /// `SHIM_ERR_BROKE_SOLID` or `SHIM_ERR_CHANGED_TOPOLOGY`, never a broken
    /// solid: the kernel's result is checked for manifoldness, the same faces
    /// and edges, gaps within the brep's tolerance, and a volume its own
    /// tessellation agrees with. `SHIM_ERR_UNSUPPORTED_FACE`: a picked curved
    /// face other than a cone has blends on it (their new shape isn't a
    /// translation), no translation keeps a blend tangent to all its
    /// neighbours, or a surface has no offset. Some moves the kernel still
    /// gets wrong, and refuses (`SHIM_ERR_BROKE_SOLID`): fillets sliding along
    /// drafted walls, and some moves of faces cut by spheres (see
    /// `PATCHES.md`).
    ///
    /// Built on `SmBrep::LocalOperation`. `faces` must be non-empty and
    /// distinct (`SHIM_ERR_BAD_FACE` otherwise); a zero distance leaves the
    /// brep as it is. Consumed and dropped on failure, as
    /// [`fillet`](Brep::fillet).
    pub fn offset_faces(self, faces: &[usize], distance: f64) -> Result<Brep, Error> {
        const OP: &str = "SmBrep::LocalOperation";
        if !distance.is_finite() {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "offset distance must be finite",
            });
        }
        if faces.is_empty() {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "no faces to offset",
            });
        }
        let indices: Vec<u32> = faces
            .iter()
            .map(|&face| u32::try_from(face).unwrap_or(u32::MAX))
            .collect();
        // SAFETY: self is a live brep on the kernel thread; indices holds
        // indices.len() values.
        let status = unsafe {
            sys::sm_brep_offset_faces(self.raw.as_ptr(), indices.as_ptr(), indices.len(), distance)
        };
        Error::check(OP, status)?;
        Ok(self)
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
    /// identical surfaces, and the faces sharing wanted faces' edges leave
    /// theirs there: see [`SurfaceCaches`]. The mesh comes out the same.
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

const EMPTY_EDGE: sys::sm_edge_geometry = sys::sm_edge_geometry {
    points: [[0.0; 3]; EDGE_SAMPLES],
    length: 0.0,
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
