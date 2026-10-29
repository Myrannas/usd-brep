//! Free-form curves and surfaces: NURBS curves from poles or through points,
//! surfaces from poles, lofted through sections, through a curve network
//! (Gordon), on four boundaries (Coons) or filling an N-sided hole, and
//! faces of them stitched into solids. The kernel's `SmBSplineCurve` and
//! `SmBSplineSurface` beneath (`NL_SrfAdv.h`).

use std::marker::PhantomData;
use std::ptr::{self, NonNull};

use smlib_sys as sys;

use crate::brep::Brep;
use crate::error::Error;
use crate::kernel::Kernel;
use crate::label::Label;

fn vec3([x, y, z]: [f64; 3]) -> sys::sm_vec3 {
    sys::sm_vec3 { x, y, z }
}

fn array(v: sys::sm_vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// A curve the kernel holds, a NURBS for everything made here. Freed on drop;
/// it stays on the kernel's thread.
pub struct Curve {
    raw: NonNull<sys::sm_curve>,
    _not_send: PhantomData<*mut ()>,
}

/// A surface the kernel holds, a NURBS for everything made here. Freed on
/// drop; it stays on the kernel's thread.
pub struct FreeSurface {
    raw: NonNull<sys::sm_surface>,
    _not_send: PhantomData<*mut ()>,
}

/// A NURBS's degree, pole count and whether it has weights.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CurveShape {
    /// Degree.
    pub degree: u32,
    /// Control points.
    pub poles: usize,
    /// Whether it has weights.
    pub rational: bool,
}

/// A NURBS surface's degrees and pole counts along u and v.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SurfaceShape {
    /// Degrees along u and v.
    pub degree: [u32; 2],
    /// Control points along u and v.
    pub poles: [usize; 2],
    /// Whether it has weights.
    pub rational: bool,
}

/// A NURBS surface's whole definition: what `FreeSurface::nurbs` takes and
/// `FreeSurface::net` gives back.
#[derive(Clone, Debug, PartialEq)]
pub struct Net {
    /// Degrees along u and v.
    pub degree: [u32; 2],
    /// `poles[i][j]`, i along u and j along v, in Euclidean coordinates.
    pub poles: Vec<Vec<[f64; 3]>>,
    /// Each pole's weight, as `poles`; `None` for a non-rational surface.
    pub weights: Option<Vec<Vec<f64>>>,
    /// The whole knot vectors along u and v, `poles + degree + 1` knots each.
    pub knots: [Vec<f64>; 2],
}

/// The point of a surface nearest another (`FreeSurface::closest_point`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Closest {
    /// Where it is.
    pub point: [f64; 3],
    /// Its (u, v) on the surface.
    pub uv: [f64; 2],
    /// How far it is from the point asked about.
    pub distance: f64,
}

/// What stitching faces into a solid did (`SmStitch::StitchIntoSolid`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stitched {
    /// Whether the faces closed into a manifold solid.
    pub solid: bool,
    /// Edges glued.
    pub edges: u32,
    /// The largest gap between vertices it closed.
    pub vertex_gap: f64,
    /// The largest gap between edges it closed.
    pub edge_gap: f64,
}

impl Curve {
    fn wrap(operation: &'static str, raw: *mut sys::sm_curve) -> Result<Curve, Error> {
        NonNull::new(raw)
            .map(|raw| Curve {
                raw,
                _not_send: PhantomData,
            })
            .ok_or(Error::Kernel {
                operation,
                status: crate::Status::from_raw(sys::SM_SHIM_ERR_NULL_ARGUMENT),
            })
    }

    /// A NURBS curve of `degree` from its poles and whole knot vector
    /// (`poles.len() + degree + 1` knots), with weights for a rational one.
    pub fn nurbs(
        _kernel: &Kernel,
        degree: u32,
        poles: &[[f64; 3]],
        weights: Option<&[f64]>,
        knots: &[f64],
    ) -> Result<Curve, Error> {
        const OP: &str = "sm_curve_create";
        if weights.is_some_and(|w| w.len() != poles.len()) {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "a weight for every pole",
            });
        }
        let raw: Vec<_> = poles.iter().copied().map(vec3).collect();
        let mut out = ptr::null_mut();
        // SAFETY: the slices outlive the call; we're on the kernel thread.
        let status = unsafe {
            sys::sm_curve_create(
                degree as i32,
                raw.as_ptr(),
                weights.map_or(ptr::null(), |w| w.as_ptr()),
                raw.len(),
                knots.as_ptr(),
                knots.len(),
                &mut out,
            )
        };
        Error::check(OP, status)?;
        Curve::wrap(OP, out)
    }

    /// A non-rational curve of `degree` through the points, chord-length
    /// parameterised; `closed` runs it on smoothly through the first again.
    pub fn through(
        _kernel: &Kernel,
        points: &[[f64; 3]],
        degree: u32,
        closed: bool,
    ) -> Result<Curve, Error> {
        const OP: &str = "SmBSplineCurve::InterpolatePoints";
        let raw: Vec<_> = points.iter().copied().map(vec3).collect();
        let mut out = ptr::null_mut();
        // SAFETY: as in `nurbs`.
        let status = unsafe {
            sys::sm_curve_interpolate(
                raw.as_ptr(),
                raw.len(),
                degree as i32,
                closed as i32,
                &mut out,
            )
        };
        Error::check(OP, status)?;
        Curve::wrap(OP, out)
    }

    /// The brep's edge's curve as a NURBS, trimmed to the edge.
    pub fn of_edge(brep: &Brep, edge: usize) -> Result<Curve, Error> {
        const OP: &str = "sm_curve_of_edge";
        let mut out = ptr::null_mut();
        // SAFETY: the brep is live on the kernel thread.
        let status = unsafe { sys::sm_curve_of_edge(brep.as_raw(), edge, &mut out) };
        Error::check(OP, status)?;
        Curve::wrap(OP, out)
    }

    /// A curve with this one's knots and degree and every pole at `vector`:
    /// a constant cross-boundary derivative to go with it.
    pub fn constant_like(&self, vector: [f64; 3]) -> Result<Curve, Error> {
        const OP: &str = "sm_curve_constant_like";
        let mut out = ptr::null_mut();
        // SAFETY: self is live on the kernel thread.
        let status =
            unsafe { sys::sm_curve_constant_like(self.raw.as_ptr(), vec3(vector), &mut out) };
        Error::check(OP, status)?;
        Curve::wrap(OP, out)
    }

    /// A copy.
    pub fn copy(&self) -> Result<Curve, Error> {
        const OP: &str = "SmCurve::Copy";
        let mut out = ptr::null_mut();
        // SAFETY: as in `constant_like`.
        let status = unsafe { sys::sm_curve_copy(self.raw.as_ptr(), &mut out) };
        Error::check(OP, status)?;
        Curve::wrap(OP, out)
    }

    /// Turns the curve round, keeping its interval.
    pub fn reverse(&mut self) -> Result<(), Error> {
        // SAFETY: as in `constant_like`.
        let status = unsafe { sys::sm_curve_reverse(self.raw.as_ptr()) };
        Error::check("SmCurve::ReverseParameterization", status)
    }

    /// The curve's parameter interval.
    pub fn interval(&self) -> Result<[f64; 2], Error> {
        let (mut min, mut max) = (0.0, 0.0);
        // SAFETY: as in `constant_like`.
        let status = unsafe { sys::sm_curve_interval(self.raw.as_ptr(), &mut min, &mut max) };
        Error::check("sm_curve_interval", status)?;
        Ok([min, max])
    }

    /// The point and first derivative at `t`.
    pub fn evaluate(&self, t: f64) -> Result<([f64; 3], [f64; 3]), Error> {
        let mut point = vec3([0.0; 3]);
        let mut derivative = vec3([0.0; 3]);
        // SAFETY: as in `constant_like`.
        let status =
            unsafe { sys::sm_curve_evaluate(self.raw.as_ptr(), t, &mut point, &mut derivative) };
        Error::check("SmCurve::Evaluate", status)?;
        Ok((array(point), array(derivative)))
    }

    /// The point a share `s` (0 to 1) of the way along its interval.
    pub fn at(&self, s: f64) -> Result<[f64; 3], Error> {
        let [a, b] = self.interval()?;
        Ok(self.evaluate(a + (b - a) * s)?.0)
    }

    /// Degree, poles and weights, for a NURBS.
    pub fn shape(&self) -> Result<CurveShape, Error> {
        let (mut degree, mut poles, mut rational) = (0, 0, 0);
        // SAFETY: as in `constant_like`.
        let status = unsafe {
            sys::sm_curve_shape(self.raw.as_ptr(), &mut degree, &mut poles, &mut rational)
        };
        Error::check("sm_curve_shape", status)?;
        Ok(CurveShape {
            degree: degree as u32,
            poles,
            rational: rational != 0,
        })
    }

    fn as_raw(&self) -> *const sys::sm_curve {
        self.raw.as_ptr()
    }
}

impl Drop for Curve {
    fn drop(&mut self) {
        // SAFETY: we own the curve and drop it once, on the kernel thread.
        unsafe { sys::sm_curve_free(self.raw.as_ptr()) }
    }
}

fn raws(curves: &[&Curve]) -> Vec<*const sys::sm_curve> {
    curves.iter().map(|c| c.as_raw()).collect()
}

impl FreeSurface {
    pub(crate) fn wrap(operation: &'static str, raw: *mut sys::sm_surface) -> Result<FreeSurface, Error> {
        NonNull::new(raw)
            .map(|raw| FreeSurface {
                raw,
                _not_send: PhantomData,
            })
            .ok_or(Error::Kernel {
                operation,
                status: crate::Status::from_raw(sys::SM_SHIM_ERR_NULL_ARGUMENT),
            })
    }

    /// A NURBS surface from its control net: `net.poles[i][j]`, i along u
    /// and j along v, of the net's degrees, with whole knot vectors along u
    /// and v (`SmBSplineSurface::CreateCanonical`).
    pub fn nurbs(_kernel: &Kernel, net: &Net) -> Result<FreeSurface, Error> {
        const OP: &str = "SmBSplineSurface::CreateCanonical";
        let Net {
            degree,
            poles,
            weights,
            knots,
        } = net;
        let count_u = poles.len();
        let count_v = poles.first().map_or(0, Vec::len);
        if poles.iter().any(|row| row.len() != count_v) {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "every row of poles the same length",
            });
        }
        let raw: Vec<_> = poles.iter().flatten().copied().map(vec3).collect();
        let flat_weights: Option<Vec<f64>> = weights
            .as_ref()
            .map(|w| w.iter().flatten().copied().collect());
        if flat_weights.as_ref().is_some_and(|w| w.len() != raw.len()) {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "a weight for every pole",
            });
        }
        let mut out = ptr::null_mut();
        // SAFETY: the buffers outlive the call; we're on the kernel thread.
        let status = unsafe {
            sys::sm_surface_create(
                degree[0] as i32,
                degree[1] as i32,
                raw.as_ptr(),
                count_u,
                count_v,
                flat_weights.as_ref().map_or(ptr::null(), |w| w.as_ptr()),
                knots[0].as_ptr(),
                knots[0].len(),
                knots[1].as_ptr(),
                knots[1].len(),
                &mut out,
            )
        };
        Error::check(OP, status)?;
        FreeSurface::wrap(OP, out)
    }

    /// A loft through the sections in order, each a u isoparameter curve of
    /// it, of up to `degree` across them
    /// (`SmBSplineSurface::CreateSkinnedSurface`).
    pub fn loft(_kernel: &Kernel, sections: &[&Curve], degree: u32) -> Result<FreeSurface, Error> {
        const OP: &str = "SmBSplineSurface::CreateSkinnedSurface";
        let raw = raws(sections);
        let mut out = ptr::null_mut();
        // SAFETY: the curves outlive the call; we're on the kernel thread.
        let status =
            unsafe { sys::sm_surface_loft(raw.as_ptr(), raw.len(), degree as i32, &mut out) };
        Error::check(OP, status)?;
        FreeSurface::wrap(OP, out)
    }

    /// A Gordon surface through a network of u curves crossing v curves
    /// (`SmBSplineSurface::CreateGordonSurface`).
    pub fn gordon(_kernel: &Kernel, u: &[&Curve], v: &[&Curve]) -> Result<FreeSurface, Error> {
        const OP: &str = "SmBSplineSurface::CreateGordonSurface";
        let (u, v) = (raws(u), raws(v));
        let mut out = ptr::null_mut();
        // SAFETY: as in `loft`.
        let status =
            unsafe { sys::sm_surface_gordon(u.as_ptr(), u.len(), v.as_ptr(), v.len(), &mut out) };
        Error::check(OP, status)?;
        FreeSurface::wrap(OP, out)
    }

    /// A Coons patch on four boundary curves round it, each running on from
    /// the last; with cross-boundary derivatives into it, a bicubic one
    /// meeting them, else bilinear (`SmBSplineSurface::CreateCoonsPatch`).
    pub fn coons(
        _kernel: &Kernel,
        boundary: [&Curve; 4],
        derivatives: Option<[&Curve; 4]>,
        tolerance: f64,
    ) -> Result<FreeSurface, Error> {
        const OP: &str = "SmBSplineSurface::CreateCoonsPatch";
        let boundary = raws(&boundary);
        let derivatives = derivatives.map(|d| raws(&d));
        let mut out = ptr::null_mut();
        // SAFETY: as in `loft`.
        let status = unsafe {
            sys::sm_surface_coons(
                boundary.as_ptr(),
                derivatives.as_ref().map_or(ptr::null(), |d| d.as_ptr()),
                tolerance,
                &mut out,
            )
        };
        Error::check(OP, status)?;
        FreeSurface::wrap(OP, out)
    }

    /// Fills an N-sided hole bounded by non-rational curves round it, each
    /// running on from the last, with N four-sided patches meeting in its
    /// middle, meeting any cross-boundary derivatives given to within
    /// `angle_tol_deg` (`SmBSplineSurface::CreateNSidedPatch`,
    /// `N_FillNSidedHole`).
    ///
    /// The fill is smooth, not exact: where three rounds of radius r meet
    /// at a box's corner, the patches meet the rounds within 0.11° but
    /// stray up to 0.043 r (4% of the radius, 0.43 mm at r 10) from the
    /// sphere the exact corner lies on, bulging out, and cross-boundary
    /// derivatives of any length change neither
    /// (`tests/freeform.rs`, `n_sided_fill_quality`). A part with such
    /// corners is a little bigger than one rounded by `Brep::fillet`, whose
    /// corners are spheres.
    pub fn fill(
        _kernel: &Kernel,
        boundary: &[&Curve],
        derivatives: Option<&[Option<&Curve>]>,
        angle_tol_deg: f64,
    ) -> Result<Vec<FreeSurface>, Error> {
        const OP: &str = "SmBSplineSurface::CreateNSidedPatch";
        if derivatives.is_some_and(|d| d.len() != boundary.len()) {
            return Err(Error::InvalidArgument {
                operation: OP,
                reason: "a derivative, or None, for every boundary curve",
            });
        }
        let bounds = raws(boundary);
        let derivs: Option<Vec<*const sys::sm_curve>> = derivatives.map(|d| {
            d.iter()
                .map(|c| c.map_or(ptr::null(), |c| c.as_raw()))
                .collect()
        });
        let mut out = vec![ptr::null_mut(); boundary.len()];
        let mut count = 0;
        // SAFETY: as in `loft`; out has room for one surface per side.
        let status = unsafe {
            sys::sm_surface_fill(
                bounds.as_ptr(),
                derivs.as_ref().map_or(ptr::null(), |d| d.as_ptr()),
                bounds.len(),
                angle_tol_deg,
                out.as_mut_ptr(),
                &mut count,
            )
        };
        Error::check(OP, status)?;
        out.into_iter()
            .take(count)
            .map(|raw| FreeSurface::wrap(OP, raw))
            .collect()
    }

    /// A copy of the brep's face's surface.
    pub fn of_face(brep: &Brep, face: usize) -> Result<FreeSurface, Error> {
        const OP: &str = "sm_surface_of_face";
        let mut out = ptr::null_mut();
        // SAFETY: the brep is live on the kernel thread.
        let status = unsafe { sys::sm_surface_of_face(brep.as_raw(), face, &mut out) };
        Error::check(OP, status)?;
        FreeSurface::wrap(OP, out)
    }

    /// The shim's handle, for the crate's other modules.
    pub(crate) fn as_raw(&self) -> *const sys::sm_surface {
        self.raw.as_ptr()
    }

    /// The surface's (u, v) domain, its lower and upper corners.
    pub fn domain(&self) -> Result<[[f64; 2]; 2], Error> {
        let (mut min, mut max) = ([0.0; 2], [0.0; 2]);
        // SAFETY: self is live on the kernel thread.
        let status = unsafe {
            sys::sm_surface_domain(self.raw.as_ptr(), min.as_mut_ptr(), max.as_mut_ptr())
        };
        Error::check("sm_surface_domain", status)?;
        Ok([min, max])
    }

    /// The point and unit normal at (u, v).
    pub fn evaluate(&self, u: f64, v: f64) -> Result<([f64; 3], [f64; 3]), Error> {
        let mut point = vec3([0.0; 3]);
        let mut normal = vec3([0.0; 3]);
        // SAFETY: as in `domain`.
        let status =
            unsafe { sys::sm_surface_evaluate(self.raw.as_ptr(), u, v, &mut point, &mut normal) };
        Error::check("SmApiEvaluateSurfacePoint", status)?;
        Ok((array(point), array(normal)))
    }

    /// The point and normal at shares `s` and `t` (0 to 1) of its domain.
    pub fn at(&self, s: f64, t: f64) -> Result<([f64; 3], [f64; 3]), Error> {
        let [min, max] = self.domain()?;
        self.evaluate(
            min[0] + (max[0] - min[0]) * s,
            min[1] + (max[1] - min[1]) * t,
        )
    }

    /// The point of the surface nearest `point`, over its whole domain
    /// (`SmApiSurfaceClosestPoint`): where, its (u, v) and how far. Of
    /// several equally near, one.
    pub fn closest_point(&self, point: [f64; 3]) -> Result<Closest, Error> {
        let mut closest = vec3([0.0; 3]);
        let mut uv = [0.0; 2];
        let mut distance = 0.0;
        // SAFETY: as in `domain`.
        let status = unsafe {
            sys::sm_surface_closest_point(
                self.raw.as_ptr(),
                vec3(point),
                &mut closest,
                uv.as_mut_ptr(),
                &mut distance,
            )
        };
        Error::check("SmApiSurfaceClosestPoint", status)?;
        Ok(Closest {
            point: array(closest),
            uv,
            distance,
        })
    }

    /// The surface's control net, as `nurbs` takes it: a NURBS's own, an
    /// analytic surface's as a NURBS.
    pub fn net(&self) -> Result<Net, Error> {
        const OP: &str = "sm_surface_net";
        let shape = self.shape()?;
        let [count_u, count_v] = shape.poles;
        let [degree_u, degree_v] = shape.degree.map(|d| d as usize);
        let mut poles = vec![vec3([0.0; 3]); count_u * count_v];
        let mut weights = vec![1.0; count_u * count_v];
        let mut knots_u = vec![0.0; count_u + degree_u + 1];
        let mut knots_v = vec![0.0; count_v + degree_v + 1];
        // SAFETY: as in `domain`; the buffers are sized as `shape` says.
        let status = unsafe {
            sys::sm_surface_net(
                self.raw.as_ptr(),
                poles.as_mut_ptr(),
                weights.as_mut_ptr(),
                knots_u.as_mut_ptr(),
                knots_v.as_mut_ptr(),
            )
        };
        Error::check(OP, status)?;
        let row = |i: usize| i * count_v..(i + 1) * count_v;
        Ok(Net {
            degree: shape.degree,
            poles: (0..count_u)
                .map(|i| poles[row(i)].iter().copied().map(array).collect())
                .collect(),
            weights: shape
                .rational
                .then(|| (0..count_u).map(|i| weights[row(i)].to_vec()).collect()),
            knots: [knots_u, knots_v],
        })
    }

    /// Degrees and pole counts of the surface as a NURBS: a NURBS's own, an
    /// analytic surface's as a NURBS.
    pub fn shape(&self) -> Result<SurfaceShape, Error> {
        let (mut du, mut dv, mut pu, mut pv, mut rational) = (0, 0, 0, 0, 0);
        // SAFETY: as in `domain`.
        let status = unsafe {
            sys::sm_surface_shape(
                self.raw.as_ptr(),
                &mut du,
                &mut dv,
                &mut pu,
                &mut pv,
                &mut rational,
            )
        };
        Error::check("sm_surface_shape", status)?;
        Ok(SurfaceShape {
            degree: [du as u32, dv as u32],
            poles: [pu, pv],
            rational: rational != 0,
        })
    }
}

impl Drop for FreeSurface {
    fn drop(&mut self) {
        // SAFETY: we own the surface and drop it once, on the kernel thread.
        unsafe { sys::sm_surface_free(self.raw.as_ptr()) }
    }
}

/// Faces gathered to be stitched into a solid: each surface on its whole
/// domain, and planar faces bounded by loops of curves.
pub struct Sheets {
    brep: Brep,
}

impl Sheets {
    /// Nothing yet.
    pub fn new(_kernel: &Kernel) -> Result<Sheets, Error> {
        const OP: &str = "sm_brep_create_empty";
        let mut out = ptr::null_mut();
        // SAFETY: a fresh brep out; we're on the kernel thread.
        let status = unsafe { sys::sm_brep_create_empty(&mut out) };
        Error::check(OP, status)?;
        Ok(Sheets {
            brep: Brep::from_raw(OP, out)?,
        })
    }

    /// Adds the surface on its whole domain, split where it's only C0,
    /// labelled; how many faces it made.
    pub fn add(&mut self, surface: &FreeSurface, label: Option<&Label>) -> Result<usize, Error> {
        const OP: &str = "SmBrep::CreateFacesFromSurface";
        let raw = label.map(|l| l.to_raw(OP)).transpose()?;
        let mut faces = 0;
        // SAFETY: both are live on the kernel thread; raw outlives the call.
        let status = unsafe {
            sys::sm_brep_add_surface(
                self.brep.as_raw(),
                surface.raw.as_ptr(),
                raw.as_ref().map_or(ptr::null(), |r| r as *const _),
                &mut faces,
            )
        };
        Error::check(OP, status)?;
        Ok(faces)
    }

    /// Adds a planar face bounded by the closed loop of curves, labelled.
    pub fn add_planar(&mut self, loop_: &[&Curve], label: Option<&Label>) -> Result<(), Error> {
        const OP: &str = "SmBrep::CreatePlanarFacesWith3DCurves";
        let raw = label.map(|l| l.to_raw(OP)).transpose()?;
        let curves = raws(loop_);
        // SAFETY: as in `add`.
        let status = unsafe {
            sys::sm_brep_add_planar_face(
                self.brep.as_raw(),
                curves.as_ptr(),
                curves.len(),
                raw.as_ref().map_or(ptr::null(), |r| r as *const _),
            )
        };
        Error::check(OP, status)
    }

    /// The faces so far, unstitched.
    pub fn brep(&self) -> &Brep {
        &self.brep
    }

    /// Stitches the faces along the edges that lie together and orients
    /// them outward (`SmBrep::StitchAndOrient`), `tolerance` 0 for the
    /// kernel's; the brep and whether it closed into a manifold solid.
    pub fn stitch(self, tolerance: f64) -> Result<(Brep, bool), Error> {
        const OP: &str = "SmBrep::StitchAndOrient";
        let mut solid = 0;
        // SAFETY: as in `add`.
        let status = unsafe { sys::sm_brep_stitch(self.brep.as_raw(), tolerance, &mut solid) };
        Error::check(OP, status)?;
        Ok((self.brep, solid != 0))
    }

    /// Stitches with `SmStitch::StitchIntoSolid`, saying what it did.
    pub fn stitch_into_solid(self) -> Result<(Brep, Stitched), Error> {
        const OP: &str = "SmStitch::StitchIntoSolid";
        let (mut solid, mut edges, mut vertex_gap, mut edge_gap) = (0, 0, 0.0, 0.0);
        // SAFETY: as in `add`.
        let status = unsafe {
            sys::sm_brep_stitch_into_solid(
                self.brep.as_raw(),
                &mut solid,
                &mut edges,
                &mut vertex_gap,
                &mut edge_gap,
            )
        };
        Error::check(OP, status)?;
        Ok((
            self.brep,
            Stitched {
                solid: solid != 0,
                edges,
                vertex_gap,
                edge_gap,
            },
        ))
    }
}

impl Curve {
    /// A whole ellipse about `centre` with semi-axes `a` and `b` (vectors,
    /// square to each other), starting at `centre + a` and turning towards
    /// `b`: rational quadratic, clamped, 9 poles, exact.
    pub fn ellipse(
        kernel: &Kernel,
        centre: [f64; 3],
        a: [f64; 3],
        b: [f64; 3],
    ) -> Result<Curve, Error> {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let at = |i: f64, j: f64| [0, 1, 2].map(|k| centre[k] + i * a[k] + j * b[k]);
        let poles = [
            at(1.0, 0.0),
            at(1.0, 1.0),
            at(0.0, 1.0),
            at(-1.0, 1.0),
            at(-1.0, 0.0),
            at(-1.0, -1.0),
            at(0.0, -1.0),
            at(1.0, -1.0),
            at(1.0, 0.0),
        ];
        let weights = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
        let knots = [
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ];
        Curve::nurbs(kernel, 2, &poles, Some(&weights), &knots)
    }

    /// A straight line from `a` to `b`, a degree-1 NURBS.
    pub fn line(kernel: &Kernel, a: [f64; 3], b: [f64; 3]) -> Result<Curve, Error> {
        Curve::nurbs(kernel, 1, &[a, b], None, &[0.0, 0.0, 1.0, 1.0])
    }
}

/// Test solids with free-form faces, for kernel tests and spikes.
pub mod samples {
    use std::f64::consts::{FRAC_PI_2, PI};

    use super::{Curve, FreeSurface, Net, Sheets};
    use crate::brep::Brep;
    use crate::error::Error;
    use crate::kernel::Kernel;

    /// A tank lofted along x from `start` through closed elliptical sections,
    /// each `radii[i]` (across y and z) at `start.x + length * i / (n - 1)`,
    /// capped flat at both ends: a loft, closed round, and two planar faces
    /// bounded by its end sections.
    pub fn lofted_tank(
        kernel: &Kernel,
        start: [f64; 3],
        length: f64,
        radii: &[[f64; 2]],
    ) -> Result<(Brep, bool), Error> {
        let n = radii.len();
        let mut sections = Vec::with_capacity(n);
        for (i, [ry, rz]) in radii.iter().enumerate() {
            let x = start[0] + length * i as f64 / (n - 1) as f64;
            sections.push(Curve::ellipse(
                kernel,
                [x, start[1], start[2]],
                [0.0, *ry, 0.0],
                [0.0, 0.0, *rz],
            )?);
        }
        let refs: Vec<&Curve> = sections.iter().collect();
        let skin = FreeSurface::loft(kernel, &refs, 3)?;
        let mut sheets = Sheets::new(kernel)?;
        sheets.add(&skin, None)?;
        sheets.add_planar(&[&sections[0]], None)?;
        sheets.add_planar(&[&sections[n - 1]], None)?;
        sheets.stitch(0.0)
    }

    /// A dome over the rectangle from `min` (x, y) by `size`, its rim at
    /// `height` and its middle poles lifted by `lift`: 5 by 5 poles, cubic,
    /// the rim's rows straight.
    pub fn dome(
        kernel: &Kernel,
        min: [f64; 2],
        size: [f64; 2],
        height: f64,
        lift: f64,
    ) -> Result<FreeSurface, Error> {
        const F: [f64; 5] = [0.0, 0.8, 1.2, 0.8, 0.0];
        let poles: Vec<Vec<[f64; 3]>> = (0..5)
            .map(|i| {
                (0..5)
                    .map(|j| {
                        [
                            min[0] + size[0] * i as f64 / 4.0,
                            min[1] + size[1] * j as f64 / 4.0,
                            height + lift * F[i] * F[j],
                        ]
                    })
                    .collect()
            })
            .collect();
        let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        let net = Net {
            degree: [3, 3],
            poles,
            weights: None,
            knots: [knots.clone(), knots],
        };
        FreeSurface::nurbs(kernel, &net)
    }

    /// A box from `min` by `size` whose top is a NURBS dome rising `lift`
    /// (its poles; the surface rises about 0.7 of it) over its walls' tops:
    /// four planar walls, a planar floor and the dome, stitched.
    pub fn domed_box(
        kernel: &Kernel,
        min: [f64; 3],
        size: [f64; 3],
        lift: f64,
    ) -> Result<(Brep, bool), Error> {
        let [x0, y0, z0] = min;
        let [x1, y1, z1] = [x0 + size[0], y0 + size[1], z0 + size[2]];
        let dome = dome(kernel, [x0, y0], [size[0], size[1]], z1, lift)?;
        let mut sheets = Sheets::new(kernel)?;
        sheets.add(&dome, None)?;
        let corners = |z: f64| [[x0, y0, z], [x1, y0, z], [x1, y1, z], [x0, y1, z]];
        let (low, high) = (corners(z0), corners(z1));
        // The floor, then a wall under each side of the rim.
        let floor: Vec<Curve> = (0..4)
            .map(|k| Curve::line(kernel, low[k], low[(k + 1) % 4]))
            .collect::<Result<_, _>>()?;
        sheets.add_planar(&floor.iter().collect::<Vec<_>>(), None)?;
        for k in 0..4 {
            let m = (k + 1) % 4;
            let wall = [
                Curve::line(kernel, low[k], low[m])?,
                Curve::line(kernel, low[m], high[m])?,
                Curve::line(kernel, high[m], high[k])?,
                Curve::line(kernel, high[k], low[k])?,
            ];
            sheets.add_planar(&wall.iter().collect::<Vec<_>>(), None)?;
        }
        sheets.stitch(0.0)
    }

    /// A quarter arc about `centre` from `centre + r a` to `centre + r b`
    /// (a and b unit and square), through 9 points, non-rational.
    fn arc(
        kernel: &Kernel,
        centre: [f64; 3],
        r: f64,
        a: [f64; 3],
        b: [f64; 3],
    ) -> Result<Curve, Error> {
        let points: Vec<[f64; 3]> = (0..9)
            .map(|k| {
                let t = FRAC_PI_2 * k as f64 / 8.0;
                let (c, s) = (t.cos(), t.sin());
                [0, 1, 2].map(|i| centre[i] + r * (a[i] * c + b[i] * s))
            })
            .collect();
        Curve::through(kernel, &points, 3, false)
    }

    /// What a rounded box is made of, for checking: how far the corner
    /// patches stray from the exact sphere.
    pub struct RoundedBox {
        /// The solid and whether it closed.
        pub brep: Brep,
        /// Whether stitching closed it.
        pub solid: bool,
        /// The largest distance from the corner patches' sample points to
        /// the sphere each corner would be on.
        pub corner_deviation: f64,
        /// The exact volume of the rounded box.
        pub exact_volume: f64,
    }

    /// A box from `min` by `size` with every edge rounded by `r`: six planar
    /// faces, twelve rounds lofted between non-rational arcs, and each
    /// corner, where three rounds meet, filled by an N-sided patch
    /// (`FreeSurface::fill`, three four-sided patches) meeting the rounds
    /// tangentially, all stitched.
    pub fn rounded_box(
        kernel: &Kernel,
        min: [f64; 3],
        size: [f64; 3],
        r: f64,
    ) -> Result<RoundedBox, Error> {
        let lo = min;
        let hi = [0, 1, 2].map(|i| min[i] + size[i]);
        let mut sheets = Sheets::new(kernel)?;
        let mut deviation: f64 = 0.0;
        // Each corner by its signs along x, y and z.
        for corner in 0..8 {
            let sign = [0, 1, 2].map(|i| if corner >> i & 1 == 1 { 1.0 } else { -1.0 });
            let centre = [0, 1, 2].map(|i| if sign[i] > 0.0 { hi[i] - r } else { lo[i] + r });
            let axis = |i: usize| {
                let mut v = [0.0; 3];
                v[i] = sign[i];
                v
            };
            // The arcs where the rounds along each axis end at this corner:
            // the round along axis i ends in the plane through the centre
            // square to i, its arc between the other two axes. Run them
            // round the corner: (x→y), (y→z), (z→x) seen from outside,
            // turned round when the corner's handedness flips.
            let flip = sign[0] * sign[1] * sign[2] < 0.0;
            let order: [(usize, usize, usize); 3] = [(2, 0, 1), (0, 1, 2), (1, 2, 0)];
            let mut boundary = Vec::new();
            let mut derivs = Vec::new();
            for &(along, a, b) in &order {
                let (a, b) = if flip { (b, a) } else { (a, b) };
                let curve = arc(kernel, centre, r, axis(a), axis(b))?;
                // Into the patch: away from the round, towards the corner.
                derivs.push(curve.constant_like(axis(along).map(|v| v * r))?);
                boundary.push(curve);
            }
            if flip {
                boundary.reverse();
                derivs.reverse();
            }
            let bounds: Vec<&Curve> = boundary.iter().collect();
            let ders: Vec<Option<&Curve>> = derivs.iter().map(Some).collect();
            let patches = FreeSurface::fill(kernel, &bounds, Some(&ders), 1.0)?;
            for patch in &patches {
                for s in 0..=4 {
                    for t in 0..=4 {
                        let (p, _) = patch.at(s as f64 / 4.0, t as f64 / 4.0)?;
                        let d = (0..3)
                            .map(|i| (p[i] - centre[i]).powi(2))
                            .sum::<f64>()
                            .sqrt();
                        deviation = deviation.max((d - r).abs());
                    }
                }
                sheets.add(patch, None)?;
            }
        }
        // A round along each edge: axis i, at the corners' signs along the
        // other two.
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            for corner in 0..4 {
                let sj = if corner & 1 == 1 { 1.0 } else { -1.0 };
                let sk = if corner & 2 == 2 { 1.0 } else { -1.0 };
                let mut centre = [0.0; 3];
                centre[j] = if sj > 0.0 { hi[j] - r } else { lo[j] + r };
                centre[k] = if sk > 0.0 { hi[k] - r } else { lo[k] + r };
                let mut a = [0.0; 3];
                a[j] = sj;
                let mut b = [0.0; 3];
                b[k] = sk;
                let ends: Vec<Curve> = [lo[i] + r, hi[i] - r]
                    .iter()
                    .map(|&at| {
                        let mut c = centre;
                        c[i] = at;
                        arc(kernel, c, r, a, b)
                    })
                    .collect::<Result<_, _>>()?;
                let round = FreeSurface::loft(kernel, &[&ends[0], &ends[1]], 1)?;
                sheets.add(&round, None)?;
            }
        }
        // A flat face on each side, inset by r.
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            for side in [lo[i], hi[i]] {
                let corners: Vec<[f64; 3]> = [(0, 0), (1, 0), (1, 1), (0, 1)]
                    .iter()
                    .map(|&(u, v)| {
                        let mut p = [0.0; 3];
                        p[i] = side;
                        p[j] = if u == 0 { lo[j] + r } else { hi[j] - r };
                        p[k] = if v == 0 { lo[k] + r } else { hi[k] - r };
                        p
                    })
                    .collect();
                let edges: Vec<Curve> = (0..4)
                    .map(|e| Curve::line(kernel, corners[e], corners[(e + 1) % 4]))
                    .collect::<Result<_, _>>()?;
                sheets.add_planar(&edges.iter().collect::<Vec<_>>(), None)?;
            }
        }
        let (brep, solid) = sheets.stitch(0.0)?;
        let [a, b, c] = size;
        let exact_volume =
            a * b * c - 4.0 * (1.0 - PI / 4.0) * r * r * (a + b + c - 6.0 * r) - 8.0 * r.powi(3)
                + 4.0 / 3.0 * PI * r.powi(3);
        Ok(RoundedBox {
            brep,
            solid,
            corner_deviation: deviation,
            exact_volume,
        })
    }
}
