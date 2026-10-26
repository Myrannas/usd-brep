use smlib_sys as sys;

use crate::error::Error;
use crate::label::Label;

/// Tessellation quality controls (`SmTessellationParams`). Each is disabled
/// by 0.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TessellationParams {
    /// Maximum distance between the mesh and the exact surface, in model
    /// units.
    pub chord_height: f64,
    /// Maximum angle between tangents along a tessellated curve, in degrees.
    pub curve_angle_deg: f64,
    /// Maximum turning angle across a tessellated surface element, in degrees.
    pub surface_angle_deg: f64,
    /// Target edge length in model units. Not a hard cap; negative values
    /// select 0.025 × |value| × the bounding-box diagonal.
    pub max_edge_length: f64,
    /// Maximum triangle aspect ratio.
    pub max_aspect_ratio: f64,
    /// Angle tolerance for the edge polylines in [`Edge::polyline`], in
    /// degrees.
    pub edge_angle_deg: f64,
}

impl Default for TessellationParams {
    /// A 0.05 chord height (fine for millimetre models) with 15° curve,
    /// surface and edge angles. 15° is the documented workaround for
    /// inside-out triangles next to fillets (KNOWN_ISSUES.md).
    fn default() -> Self {
        TessellationParams {
            chord_height: 0.05,
            curve_angle_deg: 15.0,
            surface_angle_deg: 15.0,
            max_edge_length: 0.0,
            max_aspect_ratio: 0.0,
            edge_angle_deg: 15.0,
        }
    }
}

impl TessellationParams {
    pub(crate) fn to_raw(self) -> Result<sys::sm_tessellation_params, Error> {
        let values = [
            self.chord_height,
            self.curve_angle_deg,
            self.surface_angle_deg,
            self.max_edge_length,
            self.max_aspect_ratio,
            self.edge_angle_deg,
        ];
        if values.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidArgument {
                operation: "SmApiTessellate",
                reason: "parameters must be finite",
            });
        }
        if self.edge_angle_deg <= 0.0 {
            return Err(Error::InvalidArgument {
                operation: "SmApiTessellateBoundaries",
                reason: "edge angle must be positive",
            });
        }
        Ok(sys::sm_tessellation_params {
            chord_height_tol: self.chord_height,
            curve_angle_tol_deg: self.curve_angle_deg,
            surface_angle_tol_deg: self.surface_angle_deg,
            max_edge_length: self.max_edge_length,
            max_aspect_ratio: self.max_aspect_ratio,
            edge_angle_tol_deg: self.edge_angle_deg,
        })
    }
}

/// A triangle mesh that keeps its links back to the B-rep it came from.
///
/// Vertices are split per B-rep face (and on either side of a seam), so hard
/// edges stay hard and every vertex belongs to exactly one face. Face and edge
/// ids index [`faces`](Mesh::faces) and [`edges`](Mesh::edges), which follow
/// the brep's own face and edge lists; they hold for this tessellation only.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    /// Vertex positions.
    pub positions: Vec<[f32; 3]>,
    /// Unit surface normals from the kernel, one per vertex, pointing out of
    /// the solid.
    pub normals: Vec<[f32; 3]>,
    /// Each vertex's (u, v) on its face's surface.
    pub surface_uvs: Vec<[f32; 2]>,
    /// Triangle list, three indices per triangle, counter-clockwise seen from
    /// outside.
    pub indices: Vec<u32>,
    /// The source face of each triangle, as an index into
    /// [`faces`](Mesh::faces).
    pub face_ids: Vec<u32>,
    /// The brep's faces.
    pub faces: Vec<Face>,
    /// The brep's edges.
    pub edges: Vec<Edge>,
}

impl Mesh {
    /// The number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.face_ids.len()
    }

    /// Triangles as vertex index triples.
    pub fn triangles(&self) -> impl Iterator<Item = [u32; 3]> + '_ {
        self.indices.as_chunks::<3>().0.iter().copied()
    }

    /// Moves the mesh rigidly, as [`Brep::transform`](crate::Brep::transform)
    /// moves the brep it came from: the point p goes to `rotation · p +
    /// translation`, with `rotation` given row by row, and the normals, the
    /// edges' polylines and the faces' [`Surface`]s turn with it. Each vertex
    /// keeps its surface (u, v), as each surface keeps its parameters, and
    /// each face its domain and [`Label`].
    ///
    /// A moved brep's own tessellation is these triangles to within
    /// rounding, but its [`FaceKey`](crate::FaceKey)s differ, since they hash
    /// the geometry where it lies: a caller that keeps triangles by key can
    /// keep a body's from before it moved and move them. `rotation` must be
    /// a rotation to within 1e-9, as for
    /// [`Brep::transform`](crate::Brep::transform); otherwise this fails and
    /// leaves the mesh alone.
    pub fn transform(
        &mut self,
        rotation: [[f64; 3]; 3],
        translation: [f64; 3],
    ) -> Result<(), Error> {
        if !is_rotation(&rotation) || !translation.iter().all(|c| c.is_finite()) {
            return Err(Error::InvalidArgument {
                operation: "Mesh::transform",
                reason: "the rotation must be orthonormal with determinant +1, and the translation finite",
            });
        }
        self.map(&rotation, translation, false);
        Ok(())
    }

    /// Reflects the mesh as [`Brep::mirror`](crate::Brep::mirror) reflects
    /// the brep it came from, across the plane through `point` with
    /// `normal`: the point p goes to `p - 2 ((p - point) · n) n` for the unit
    /// normal n. The normals reflect too, so they still point out of the
    /// solid, and every triangle turns round, so it is still
    /// counter-clockwise seen from outside. Each vertex keeps its surface
    /// (u, v), as each surface keeps its parameters; a surface's angles now
    /// run the other way round its axis, which [`Surface`] doesn't record.
    /// As with [`transform`](Mesh::transform), the mirrored brep's own
    /// tessellation is these triangles to within rounding. Fails, leaving
    /// the mesh alone, for a point or normal that isn't finite or a zero
    /// normal.
    pub fn mirror(&mut self, point: [f64; 3], normal: [f64; 3]) -> Result<(), Error> {
        let length = normal.iter().map(|c| c * c).sum::<f64>().sqrt();
        if !point.iter().all(|c| c.is_finite()) || !length.is_finite() || length == 0.0 {
            return Err(Error::InvalidArgument {
                operation: "Mesh::mirror",
                reason: "the point and normal must be finite, and the normal not zero",
            });
        }
        let n = normal.map(|c| c / length);
        let reflection: [[f64; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| f64::from(i == j) - 2.0 * n[i] * n[j]));
        let along = 2.0 * (point[0] * n[0] + point[1] * n[1] + point[2] * n[2]);
        self.map(&reflection, n.map(|c| along * c), true);
        Ok(())
    }

    /// Takes every point p to `linear · p + translation` and every direction
    /// d to `linear · d`, turning the triangles round when `linear` reflects.
    fn map(&mut self, linear: &[[f64; 3]; 3], translation: [f64; 3], reflects: bool) {
        let vector = |v: [f64; 3]| linear.map(|row| row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
        let point = |p: [f64; 3]| {
            let q = vector(p);
            [q[0] + translation[0], q[1] + translation[1], q[2] + translation[2]]
        };
        let narrow = |p: [f64; 3]| p.map(|c| c as f32);
        for p in self.positions.iter_mut().chain(self.edges.iter_mut().flat_map(|e| &mut e.polyline)) {
            *p = narrow(point(p.map(f64::from)));
        }
        for n in &mut self.normals {
            *n = narrow(vector(n.map(f64::from)));
        }
        for face in &mut self.faces {
            face.surface = face.surface.mapped(point, vector);
        }
        if reflects {
            for triangle in self.indices.as_chunks_mut::<3>().0 {
                triangle.swap(1, 2);
            }
        }
    }

    /// Copies a shim mesh into Rust-owned arrays.
    ///
    /// # Safety
    ///
    /// `raw` must be a live mesh from `sm_brep_tessellate`.
    pub(crate) unsafe fn from_raw(raw: *const sys::sm_mesh) -> Mesh {
        let mut view = std::mem::MaybeUninit::<sys::sm_mesh_view>::zeroed();
        // SAFETY: raw is live; the view borrows its arrays until we free it.
        let view = unsafe {
            sys::sm_mesh_view_get(raw, view.as_mut_ptr());
            view.assume_init()
        };
        // SAFETY: each pointer covers the count the view gives for it.
        unsafe {
            let positions = slice(view.positions, view.vertex_count * 3);
            let normals = slice(view.normals, view.vertex_count * 3);
            let uvs = slice(view.surface_uvs, view.vertex_count * 2);
            let edge_points = slice(view.edge_points, view.edge_point_count * 3);
            Mesh {
                positions: positions.as_chunks::<3>().0.to_vec(),
                normals: normals.as_chunks::<3>().0.to_vec(),
                surface_uvs: uvs.as_chunks::<2>().0.to_vec(),
                indices: slice(view.indices, view.triangle_count * 3).to_vec(),
                face_ids: slice(view.face_ids, view.triangle_count).to_vec(),
                faces: slice(view.faces, view.face_count)
                    .iter()
                    .map(Face::from_raw)
                    .collect(),
                edges: slice(view.edges, view.edge_count)
                    .iter()
                    .map(|e| Edge::from_raw(e, edge_points))
                    .collect(),
            }
        }
    }
}

/// A B-rep face.
#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    /// The surface the face lies on.
    pub surface: Surface,
    /// The face's domain in its surface's (u, v): minimum corner.
    pub uv_min: [f64; 2],
    /// The face's domain in its surface's (u, v): maximum corner.
    pub uv_max: [f64; 2],
    /// The face's [`Label`], if it has one.
    pub label: Option<Label>,
}

impl Face {
    pub(crate) fn from_raw(raw: &sys::sm_face_info) -> Face {
        let (origin, axis, ref_dir, radius) = (raw.origin, raw.axis, raw.ref_dir, raw.radius);
        let surface = match raw.kind as sys::sm_surface_kind {
            sys::SM_SURFACE_PLANE => Surface::Plane {
                origin,
                normal: axis,
                u_axis: ref_dir,
            },
            sys::SM_SURFACE_CYLINDER => Surface::Cylinder {
                origin,
                axis,
                ref_dir,
                radius,
            },
            sys::SM_SURFACE_CONE => Surface::Cone {
                origin,
                axis,
                ref_dir,
                radius,
                half_angle: raw.half_angle_rad,
            },
            sys::SM_SURFACE_SPHERE => Surface::Sphere {
                center: origin,
                axis,
                ref_dir,
                radius,
            },
            sys::SM_SURFACE_TORUS => Surface::Torus {
                center: origin,
                axis,
                ref_dir,
                major_radius: radius,
                minor_radius: raw.minor_radius,
            },
            sys::SM_SURFACE_NURBS => Surface::Nurbs {
                degree: [raw.degree_u, raw.degree_v],
                control_points: [raw.poles_u, raw.poles_v],
                rational: raw.rational != 0,
            },
            _ => Surface::Other {
                kernel_type: raw.kernel_type,
            },
        };
        Face {
            surface,
            uv_min: raw.uv_min,
            uv_max: raw.uv_max,
            label: Label::from_raw(&raw.label),
        }
    }
}

/// The surface a face lies on, with its defining parameters.
///
/// `axis` and `ref_dir` are unit vectors: the axis of revolution and the
/// direction of angle zero about it.
#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    /// A plane through `origin`. `normal` is the surface normal, which may
    /// face either way; the mesh normals point out of the solid.
    Plane {
        /// A point on the plane: the origin of its (u, v).
        origin: [f64; 3],
        /// The plane's normal.
        normal: [f64; 3],
        /// The direction of increasing u.
        u_axis: [f64; 3],
    },
    /// A circular cylinder.
    Cylinder {
        /// A point on the axis.
        origin: [f64; 3],
        /// The axis direction.
        axis: [f64; 3],
        /// The direction of angle zero.
        ref_dir: [f64; 3],
        /// The radius.
        radius: f64,
    },
    /// A circular cone.
    Cone {
        /// A point on the axis.
        origin: [f64; 3],
        /// The axis direction.
        axis: [f64; 3],
        /// The direction of angle zero.
        ref_dir: [f64; 3],
        /// The radius in the plane through `origin`.
        radius: f64,
        /// The semi-angle in radians; positive when the cone widens along
        /// `axis`.
        half_angle: f64,
    },
    /// A sphere.
    Sphere {
        /// The centre.
        center: [f64; 3],
        /// The polar axis.
        axis: [f64; 3],
        /// The direction of angle zero.
        ref_dir: [f64; 3],
        /// The radius.
        radius: f64,
    },
    /// A torus.
    Torus {
        /// The centre.
        center: [f64; 3],
        /// The axis of revolution.
        axis: [f64; 3],
        /// The direction of angle zero.
        ref_dir: [f64; 3],
        /// The distance from the axis to the tube's centre.
        major_radius: f64,
        /// The tube's radius.
        minor_radius: f64,
    },
    /// A NURBS surface.
    Nurbs {
        /// Degree in u and v.
        degree: [u32; 2],
        /// Control-point count in u and v.
        control_points: [u32; 2],
        /// Whether it has weights.
        rational: bool,
    },
    /// Another surface class, such as an offset or swept surface.
    Other {
        /// SMLib's surface type id (`SmSurfTypes.h`).
        kernel_type: i32,
    },
}

impl Surface {
    /// The surface with its points taken by `point` and its directions by
    /// `vector`, which together are a rigid motion or a reflection.
    fn mapped(&self, point: impl Fn([f64; 3]) -> [f64; 3], vector: impl Fn([f64; 3]) -> [f64; 3]) -> Surface {
        match *self {
            Surface::Plane {
                origin,
                normal,
                u_axis,
            } => Surface::Plane {
                origin: point(origin),
                normal: vector(normal),
                u_axis: vector(u_axis),
            },
            Surface::Cylinder {
                origin,
                axis,
                ref_dir,
                radius,
            } => Surface::Cylinder {
                origin: point(origin),
                axis: vector(axis),
                ref_dir: vector(ref_dir),
                radius,
            },
            Surface::Cone {
                origin,
                axis,
                ref_dir,
                radius,
                half_angle,
            } => Surface::Cone {
                origin: point(origin),
                axis: vector(axis),
                ref_dir: vector(ref_dir),
                radius,
                half_angle,
            },
            Surface::Sphere {
                center,
                axis,
                ref_dir,
                radius,
            } => Surface::Sphere {
                center: point(center),
                axis: vector(axis),
                ref_dir: vector(ref_dir),
                radius,
            },
            Surface::Torus {
                center,
                axis,
                ref_dir,
                major_radius,
                minor_radius,
            } => Surface::Torus {
                center: point(center),
                axis: vector(axis),
                ref_dir: vector(ref_dir),
                major_radius,
                minor_radius,
            },
            Surface::Nurbs { .. } | Surface::Other { .. } => self.clone(),
        }
    }
}

/// Whether `r`, row by row, is a rotation to within 1e-9: unit columns,
/// perpendicular, the third the cross product of the first two, as the shim
/// checks for [`Brep::transform`](crate::Brep::transform).
fn is_rotation(r: &[[f64; 3]; 3]) -> bool {
    let column = |j: usize| [r[0][j], r[1][j], r[2][j]];
    let [x, y, z] = [column(0), column(1), column(2)];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = [
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    ];
    let unit = |v: [f64; 3]| (dot(v, v).sqrt() - 1.0).abs() <= 1e-9;
    unit(x) && unit(y) && dot(x, y).abs() <= 1e-9 && (0..3).all(|i| (cross[i] - z[i]).abs() <= 1e-9)
}

/// A B-rep edge: where two faces meet, and so a candidate UV seam.
#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
    /// Sample points along the edge (`SmApiTessellateBoundaries`).
    pub polyline: Vec<[f32; 3]>,
    /// The faces the edge joins, as indices into [`Mesh::faces`]. A seam
    /// joins a face to itself. `None` when the edge has fewer than two uses.
    pub faces: [Option<u32>; 2],
    /// How many faces use the edge: 2 on a manifold solid, more when
    /// non-manifold (only the first two are in `faces`).
    pub use_count: u32,
    /// The edge's own [`Label`], if it has one.
    pub label: Option<Label>,
}

impl Edge {
    fn from_raw(raw: &sys::sm_edge_info, points: &[f32]) -> Edge {
        let start = raw.first_point as usize * 3;
        let end = start + raw.point_count as usize * 3;
        let face = |id: i32| u32::try_from(id).ok();
        Edge {
            polyline: points
                .get(start..end)
                .unwrap_or_default()
                .as_chunks::<3>()
                .0
                .to_vec(),
            faces: [face(raw.face_a), face(raw.face_b)],
            use_count: raw.use_count,
            label: Label::from_raw(&raw.label),
        }
    }
}

/// A slice over `len` elements at `ptr`, which may be null when `len` is 0.
///
/// # Safety
///
/// When `len` is non-zero, `ptr` must point to `len` initialised elements that
/// outlive the returned slice.
pub(crate) unsafe fn slice<'a, T>(ptr: *const T, len: usize) -> &'a [T] {
    if len == 0 || ptr.is_null() {
        &[]
    } else {
        // SAFETY: per this function's contract.
        unsafe { std::slice::from_raw_parts(ptr, len) }
    }
}
