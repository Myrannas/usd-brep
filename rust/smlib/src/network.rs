//! A surface through a network of curves: sections one way, each a curve
//! the surface passes through at a constant v, and guides the other, each a
//! curve it passes through at a constant u, every guide crossing every
//! section. A loft pulled along guides is one; so is a network of curves
//! each way.
//!
//! A Gordon surface needs each family to cross the other at parameters they
//! share: section i meets guide j at u = σ_j on every section, and guide j
//! meets section i at v = τ_i on every guide. Curves drawn by hand never do,
//! and the kernel's own Gordon surface reparameterises them piecewise
//! linearly, which leaves a crease along every crossing curve. This is the
//! compatibility step instead ([`FreeSurface::network`]): each crossing is
//! found where the two curves come nearest, a gap wider than asked refused
//! by how far it is, σ_j and τ_i taken as the crossings' means, and each
//! curve reparameterised smoothly, by a monotone cubic through its own
//! crossings to the shared ones (`nurbs::Monotone`), and refitted as a
//! non-rational cubic within `nurbs::FIT_TOLERANCE`, its corners kept as
//! corners. The kernel then makes the surface through the refitted curves at
//! those parameters (`sm_surface_gordon_at`), the middle of each crossing's
//! two points its own, so a gap within the one asked is closed there.
//!
//! Closed sections, a duct's, each start where the first guide crosses
//! them, and that guide closes the network again at the far end, so the
//! surface comes round on itself. Open sections whose ends no guide runs
//! along have one added there: the curve through the sections' ends at the
//! τ the guides share, as a skin through them would have its edge.

use std::fmt;
use std::ptr;

use smlib_sys as sys;

use crate::error::Error;
use crate::freeform::{Curve, FreeSurface};
use crate::kernel::Kernel;
use crate::nurbs::{CurveNet, FIT_TOLERANCE, Monotone};

/// The surface a network makes, and the widest gap between a guide and a
/// section it closed.
pub struct Network {
    /// The surface: u round the sections, v along the guides.
    pub surface: FreeSurface,
    /// How far apart the farthest crossing's two curves were.
    pub widest: f64,
}

/// Why a network makes no surface.
#[derive(Debug)]
pub enum NetworkError {
    /// Guide `guide` comes no nearer section `section` than `by`, more than
    /// the gap allowed.
    Misses {
        /// The section, by its place among them.
        section: usize,
        /// The guide, by its place among them.
        guide: usize,
        /// How near it comes.
        by: f64,
    },
    /// Guide `guide` doesn't cross the sections in their order, one after
    /// another.
    OutOfOrder {
        /// The guide.
        guide: usize,
    },
    /// The guides cross section `section` in another order round it than
    /// they cross the first, or two of them at one place.
    Crossed {
        /// The section.
        section: usize,
    },
    /// The kernel refused it, or what it was given.
    Kernel(Error),
}

impl fmt::Display for NetworkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NetworkError::Misses { section, guide, by } => write!(f, "guide {guide} misses section {section} by {by}"),
            NetworkError::OutOfOrder { guide } => write!(f, "guide {guide} doesn't cross the sections in order"),
            NetworkError::Crossed { section } => write!(f, "the guides cross section {section} out of order"),
            NetworkError::Kernel(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for NetworkError {}

impl From<Error> for NetworkError {
    fn from(error: Error) -> Self {
        NetworkError::Kernel(error)
    }
}

/// Where two curves come nearest: the parameter on each, and how far apart.
#[derive(Clone, Copy, Debug)]
struct Crossing {
    on_section: f64,
    on_guide: f64,
    gap: f64,
}

/// Samples taken along each curve to start finding where two come nearest.
const SAMPLES: usize = 128;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    dot(sub(a, b), sub(a, b)).sqrt()
}

/// Where `section` and `guide`, both on 0 to 1, come nearest: the nearest
/// pair of samples, then Gauss and Newton's steps on the distance squared,
/// round a closed section's seam where it has one.
fn crossing(section: &CurveNet, guide: &CurveNet, closed: bool) -> Crossing {
    let a: Vec<[f64; 3]> = (0..=SAMPLES).map(|k| section.point(k as f64 / SAMPLES as f64)).collect();
    let b: Vec<[f64; 3]> = (0..=SAMPLES).map(|k| guide.point(k as f64 / SAMPLES as f64)).collect();
    let (mut s, mut g, mut best) = (0.0, 0.0, f64::INFINITY);
    for (i, &p) in a.iter().enumerate() {
        for (j, &q) in b.iter().enumerate() {
            let d = distance(p, q);
            if d < best {
                (s, g, best) = (i as f64 / SAMPLES as f64, j as f64 / SAMPLES as f64, d);
            }
        }
    }
    for _ in 0..60 {
        let (p, dp) = section.derivative(s);
        let (q, dq) = guide.derivative(g);
        let r = sub(p, q);
        // J = [dp, -dq]; solve JᵀJ x = -Jᵀ r.
        let (a11, a12, a22) = (dot(dp, dp), -dot(dp, dq), dot(dq, dq));
        let (b1, b2) = (-dot(dp, r), dot(dq, r));
        let det = a11 * a22 - a12 * a12;
        if det.abs() < 1e-30 {
            break;
        }
        let ds = (b1 * a22 - a12 * b2) / det;
        let dg = (a11 * b2 - a12 * b1) / det;
        let mut next_s = s + ds;
        if closed {
            next_s = next_s.rem_euclid(1.0);
        } else {
            next_s = next_s.clamp(0.0, 1.0);
        }
        let next_g = (g + dg).clamp(0.0, 1.0);
        let moved = (next_s - s).abs() + (next_g - g).abs();
        (s, g) = (next_s, next_g);
        if moved < 1e-13 {
            break;
        }
    }
    Crossing { on_section: s, on_guide: g, gap: distance(section.point(s), guide.point(g)) }
}

/// The curve through `points` at `at` (rising from 0 to 1), of degree up
/// to 3: the edge a skin through the sections' ends has, near enough. A
/// polyline through them where they can't take one.
fn through(points: &[[f64; 3]], at: &[f64]) -> CurveNet {
    let degree = (points.len() - 1).min(3) as u32;
    CurveNet::through_at(points, at, degree).unwrap_or_else(|| {
        let mut knots = vec![0.0];
        knots.extend_from_slice(at);
        knots.push(1.0);
        CurveNet { degree: 1, poles: points.to_vec(), weights: None, knots }
    })
}

impl FreeSurface {
    /// The surface through `sections` pulled along `guides` (each crossing
    /// every section, within `gap`), of `degree` across the sections, at
    /// most one less than their count: see the module's notes. The sections
    /// all closed or all open, each running the same way round or along,
    /// in order; a guide may run either way. Refused, saying which, when a
    /// guide misses a section by more than `gap`, crosses them out of order,
    /// or the guides cross a section in another order than the first.
    pub fn network(kernel: &Kernel, sections: &[CurveNet], guides: &[CurveNet], degree: u32, gap: f64) -> Result<Network, NetworkError> {
        const OP: &str = "FreeSurface::network";
        let invalid = |reason: &'static str| NetworkError::Kernel(Error::InvalidArgument { operation: OP, reason });
        let n = sections.len();
        if n < 2 || guides.is_empty() {
            return Err(invalid("two sections or more, and a guide"));
        }
        if degree < 1 || degree as usize >= n {
            return Err(invalid("a degree from 1 to one less than the sections"));
        }
        if sections.iter().chain(guides).any(|curve| !curve.is_valid()) {
            return Err(invalid("valid curves"));
        }
        let unit = |curve: &CurveNet| curve.reparameterised([0.0, 1.0]);
        let mut sections: Vec<CurveNet> = sections.iter().map(unit).collect();
        let mut guides: Vec<CurveNet> = guides.iter().map(unit).collect();
        let closed_at = |curve: &CurveNet| distance(curve.point(0.0), curve.point(1.0)) <= gap.max(1e-9);
        let closed = closed_at(&sections[0]);
        if sections.iter().any(|section| closed_at(section) != closed) {
            return Err(invalid("sections all closed or all open"));
        }

        // Where each guide crosses each section.
        let mut crossings: Vec<Vec<Crossing>> = sections.iter().map(|section| guides.iter().map(|guide| crossing(section, guide, closed)).collect()).collect();
        let mut widest: f64 = 0.0;
        for (i, row) in crossings.iter().enumerate() {
            for (j, crossing) in row.iter().enumerate() {
                if crossing.gap > gap {
                    return Err(NetworkError::Misses { section: i, guide: j, by: crossing.gap });
                }
                widest = widest.max(crossing.gap);
            }
        }

        // Each guide from the first section to the last, and on to no more.
        for (j, guide) in guides.iter_mut().enumerate() {
            let along: Vec<f64> = crossings.iter().map(|row| row[j].on_guide).collect();
            if along.windows(2).all(|pair| pair[1] < pair[0]) {
                *guide = guide.reversed();
                for row in &mut crossings {
                    row[j].on_guide = 1.0 - row[j].on_guide;
                }
            } else if !along.windows(2).all(|pair| pair[1] > pair[0]) {
                return Err(NetworkError::OutOfOrder { guide: j });
            }
            let (first, last) = (crossings[0][j].on_guide, crossings[n - 1][j].on_guide);
            let mut trimmed = guide.clone();
            if last < 1.0 - 1e-12 {
                trimmed = trimmed.split(last).0;
            }
            if first > 1e-12 {
                trimmed = trimmed.split(first).1;
            }
            *guide = trimmed.reparameterised([0.0, 1.0]);
            for row in &mut crossings {
                row[j].on_guide = ((row[j].on_guide - first) / (last - first)).clamp(0.0, 1.0);
            }
        }

        // The guides in the order they cross the first section; a closed
        // section starting where the first of them crosses it.
        let mut order: Vec<usize> = (0..guides.len()).collect();
        order.sort_by(|&a, &b| crossings[0][a].on_section.total_cmp(&crossings[0][b].on_section));
        guides = order.iter().map(|&j| guides[j].clone()).collect();
        crossings = crossings.iter().map(|row| order.iter().map(|&j| row[j]).collect()).collect();
        let m = guides.len();
        if closed {
            for (i, section) in sections.iter_mut().enumerate() {
                let start = crossings[i][0].on_section;
                if start > 1e-12 && start < 1.0 - 1e-12 {
                    let (before, after) = section.split(start);
                    *section = CurveNet::joined_at(&[after, before], &[0.0, 1.0 - start, 1.0]).ok_or_else(|| invalid("a section to start where a guide crosses it"))?;
                }
                for crossing in &mut crossings[i] {
                    crossing.on_section = (crossing.on_section - start).rem_euclid(1.0);
                }
                crossings[i][0].on_section = 0.0;
            }
        }
        for (i, row) in crossings.iter().enumerate() {
            if row.windows(2).any(|pair| pair[1].on_section <= pair[0].on_section) {
                return Err(NetworkError::Crossed { section: i });
            }
        }

        // The shared parameters: τ along the guides, σ round the sections.
        let tau: Vec<f64> = (0..n)
            .map(|i| match i {
                0 => 0.0,
                _ if i == n - 1 => 1.0,
                _ => crossings[i].iter().map(|c| c.on_guide).sum::<f64>() / m as f64,
            })
            .collect();
        if tau.windows(2).any(|pair| pair[1] <= pair[0]) {
            return Err(NetworkError::OutOfOrder { guide: 0 });
        }
        // Each section's crossings round it: (σ, s) pairs, then the guides,
        // with the first again at the far end of a closed one, and an edge
        // added where an open one's ends have no guide.
        let mut columns: Vec<(CurveNet, Vec<f64>)> = guides.iter().enumerate().map(|(j, guide)| (guide.clone(), crossings.iter().map(|row| row[j].on_guide).collect())).collect();
        let mut on_sections: Vec<Vec<f64>> = crossings.iter().map(|row| row.iter().map(|c| c.on_section).collect()).collect();
        if closed {
            columns.push(columns[0].clone());
            for row in &mut on_sections {
                row.push(1.0);
            }
        } else {
            let at_end = |on_sections: &[Vec<f64>], end: f64, j: usize| {
                (0..n).all(|i| distance(sections[i].point(end), sections[i].point(on_sections[i][j])) <= gap.max(1e-9))
            };
            if at_end(&on_sections, 0.0, 0) {
                for row in &mut on_sections {
                    row[0] = 0.0;
                }
            } else {
                let edge = through(&sections.iter().map(|section| section.point(0.0)).collect::<Vec<_>>(), &tau);
                columns.insert(0, (edge, tau.clone()));
                for row in &mut on_sections {
                    row.insert(0, 0.0);
                }
            }
            let last = on_sections[0].len() - 1;
            if at_end(&on_sections, 1.0, last) {
                for row in &mut on_sections {
                    row[last] = 1.0;
                }
            } else {
                let edge = through(&sections.iter().map(|section| section.point(1.0)).collect::<Vec<_>>(), &tau);
                columns.push((edge, tau.clone()));
                for row in &mut on_sections {
                    row.push(1.0);
                }
            }
        }
        let count = columns.len();
        let sigma: Vec<f64> = (0..count)
            .map(|j| match j {
                0 => 0.0,
                _ if j == count - 1 => 1.0,
                _ => on_sections.iter().map(|row| row[j]).sum::<f64>() / n as f64,
            })
            .collect();
        if sigma.windows(2).any(|pair| pair[1] <= pair[0]) {
            return Err(NetworkError::Crossed { section: 0 });
        }

        // Each curve reparameterised to cross at the shared parameters, and
        // refitted.
        let refit = |curve: &CurveNet, own: &[f64], shared: &[f64]| -> Result<CurveNet, NetworkError> {
            let points: Vec<[f64; 2]> = shared.iter().zip(own).map(|(&x, &y)| [x, y]).collect();
            let map = Monotone::through(&points).ok_or_else(|| invalid("crossings that rise along each curve"))?;
            // Corners stay corners; elsewhere it may bend sharply where the
            // curve's own spans meet and where the map's do.
            let corners = curve.corners(1e-6);
            let joints: Vec<f64> = curve.breaks().into_iter().filter(|t| !corners.contains(t)).map(|t| map.inverse(t)).chain(shared[1..shared.len() - 1].iter().copied()).collect();
            let corners: Vec<f64> = corners.into_iter().map(|c| map.inverse(c)).collect();
            Ok(CurveNet::fit(|t| curve.point(map.at(t)), &corners, &joints, FIT_TOLERANCE))
        };
        let mut u_curves = Vec::with_capacity(n);
        for (section, own) in sections.iter().zip(&on_sections) {
            u_curves.push(Curve::from_net(kernel, &refit(section, own, &sigma)?)?);
        }
        let mut v_curves = Vec::with_capacity(count);
        for (guide, own) in &columns {
            v_curves.push(Curve::from_net(kernel, &refit(guide, own, &tau)?)?);
        }
        let degree_u = (count - 1).min(3) as i32;
        let u: Vec<*const sys::sm_curve> = u_curves.iter().map(Curve::as_raw).collect();
        let v: Vec<*const sys::sm_curve> = v_curves.iter().map(Curve::as_raw).collect();
        let mut out = ptr::null_mut();
        // SAFETY: the curves and parameters outlive the call; we're on the
        // kernel thread.
        let status = unsafe {
            sys::sm_surface_gordon_at(u.as_ptr(), tau.as_ptr(), u.len(), v.as_ptr(), sigma.as_ptr(), v.len(), degree_u, degree as i32, &mut out)
        };
        Error::check("SmBSplineSurface::CreateGordonSurfaceAt", status)?;
        Ok(Network { surface: FreeSurface::wrap("SmBSplineSurface::CreateGordonSurfaceAt", out)?, widest })
    }
}
