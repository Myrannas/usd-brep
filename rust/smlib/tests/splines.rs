//! Splines in sketch profiles and paths, as the free-form parts spec has
//! them: planar faces bounded by splines, a pipe and a profile swept along a
//! spline path, and faces turned from a spline read back as NURBS. Each is
//! checked against the spline's own evaluation here, de Boor's algorithm and
//! Gauss quadrature, not the kernel's: areas and volumes against integrals of
//! the spline, walls against their distance from it, nets against the kernel's
//! evaluation and its closest point. The spec's probe cases are here too.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    BoundaryCurve, Brep, Curve, Error, Frame, FreeSurface, Kernel, Label, Mesh, Net, PathEnd,
    PipeFace, PlaneCurve, ProfileFace, Sheets, Surface, SweptFace, TessellationParams,
};

type Job = Box<dyn FnOnce(&Kernel) + Send>;

fn on_kernel<R: Send + 'static>(f: impl FnOnce(&Kernel) -> R + Send + 'static) -> R {
    static JOBS: OnceLock<mpsc::Sender<Job>> = OnceLock::new();
    let jobs = JOBS.get_or_init(|| {
        let (tx, rx) = mpsc::channel::<Job>();
        std::thread::Builder::new()
            .name("smlib-kernel".into())
            .spawn(move || {
                let kernel = Kernel::init().expect("kernel init");
                for job in rx {
                    let _ = panic::catch_unwind(AssertUnwindSafe(|| job(&kernel)));
                }
            })
            .unwrap();
        tx
    });
    let (reply, result) = mpsc::sync_channel(1);
    jobs.send(Box::new(move |k| {
        let _ = reply.send(f(k));
    }))
    .unwrap();
    result.recv().expect("kernel job panicked")
}

// ---------------------------------------------------------------------------
// The splines, evaluated here.

/// A clamped non-rational B-spline in a plane.
#[derive(Clone, Debug)]
struct Spline {
    degree: usize,
    poles: Vec<[f64; 2]>,
    knots: Vec<f64>,
}

impl Spline {
    /// Of `degree` on `poles`, with uniform clamped knots from 0 to 1.
    fn clamped(degree: usize, poles: &[[f64; 2]]) -> Spline {
        let spans = poles.len() - degree;
        let mut knots = vec![0.0; degree + 1];
        knots.extend((1..spans).map(|i| i as f64 / spans as f64));
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        Spline {
            degree,
            poles: poles.to_vec(),
            knots,
        }
    }

    fn curve(&self) -> PlaneCurve {
        PlaneCurve::Spline {
            degree: self.degree as u32,
            poles: self.poles.clone(),
            knots: self.knots.clone(),
        }
    }

    fn domain(&self) -> (f64, f64) {
        (self.knots[0], self.knots[self.knots.len() - 1])
    }

    /// The point at `t`, by de Boor's algorithm.
    fn point(&self, t: f64) -> [f64; 2] {
        let p = self.degree;
        let n = self.poles.len();
        let t = t.clamp(self.knots[p], self.knots[n]);
        let k = (p..n)
            .rev()
            .find(|&i| self.knots[i] <= t && self.knots[i] < self.knots[i + 1])
            .unwrap_or(p);
        let mut d: Vec<[f64; 2]> = (0..=p).map(|j| self.poles[j + k - p]).collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = j + k - p;
                let span = self.knots[i + p + 1 - r] - self.knots[i];
                let alpha = if span > 0.0 {
                    (t - self.knots[i]) / span
                } else {
                    0.0
                };
                d[j] = [
                    (1.0 - alpha) * d[j - 1][0] + alpha * d[j][0],
                    (1.0 - alpha) * d[j - 1][1] + alpha * d[j][1],
                ];
            }
        }
        d[p]
    }

    /// Its derivative, a spline of one degree less.
    fn hodograph(&self) -> Spline {
        let p = self.degree;
        let poles = (0..self.poles.len() - 1)
            .map(|i| {
                let span = self.knots[i + p + 1] - self.knots[i + 1];
                let f = p as f64 / span;
                [
                    f * (self.poles[i + 1][0] - self.poles[i][0]),
                    f * (self.poles[i + 1][1] - self.poles[i][1]),
                ]
            })
            .collect();
        Spline {
            degree: p - 1,
            poles,
            knots: self.knots[1..self.knots.len() - 1].to_vec(),
        }
    }

    /// The knot spans that aren't empty.
    fn spans(&self) -> Vec<(f64, f64)> {
        self.knots
            .windows(2)
            .filter(|w| w[1] > w[0])
            .map(|w| (w[0], w[1]))
            .collect()
    }

    /// ∫ f(t) dt over the whole spline, by 16-point Gauss–Legendre on each
    /// knot span cut in `pieces`.
    fn integrate(&self, pieces: usize, f: impl Fn(f64) -> f64) -> f64 {
        let mut total = 0.0;
        for (a, b) in self.spans() {
            for k in 0..pieces {
                let (lo, hi) = (
                    a + (b - a) * k as f64 / pieces as f64,
                    a + (b - a) * (k + 1) as f64 / pieces as f64,
                );
                let (mid, half) = ((lo + hi) / 2.0, (hi - lo) / 2.0);
                for (x, w) in GAUSS_16 {
                    for sign in [-1.0, 1.0] {
                        total += half * w * f(mid + sign * half * x);
                    }
                }
            }
        }
        total
    }

    /// The radius of its tightest bend, sampled finely.
    fn tightest(&self) -> f64 {
        let d1 = self.hodograph();
        let d2 = d1.hodograph();
        let (lo, hi) = self.domain();
        (0..=20000)
            .map(|k| {
                let t = lo + (hi - lo) * k as f64 / 20000.0;
                let ([x1, y1], [x2, y2]) = (d1.point(t), d2.point(t));
                x1.hypot(y1).powi(3) / (x1 * y2 - y1 * x2).abs()
            })
            .fold(f64::INFINITY, f64::min)
    }

    fn length(&self) -> f64 {
        let d = self.hodograph();
        self.integrate(32, |t| {
            let [x, y] = d.point(t);
            x.hypot(y)
        })
    }

    /// Its parameter nearest `q` and how far `q` is from it.
    fn closest(&self, q: [f64; 2]) -> (f64, f64) {
        let d1 = self.hodograph();
        let d2 = d1.hodograph();
        let (lo, hi) = self.domain();
        let distance = |t: f64| {
            let [x, y] = self.point(t);
            (x - q[0]).hypot(y - q[1])
        };
        let samples = 400 * self.spans().len();
        let mut best = lo;
        for k in 0..=samples {
            let t = lo + (hi - lo) * k as f64 / samples as f64;
            if distance(t) < distance(best) {
                best = t;
            }
        }
        // Newton on (C - q)·C' = 0.
        let mut t = best;
        for _ in 0..50 {
            let c = self.point(t);
            let [dx, dy] = d1.point(t);
            let [ex, ey] = d2.point(t);
            let (rx, ry) = (c[0] - q[0], c[1] - q[1]);
            let f = rx * dx + ry * dy;
            let df = dx * dx + dy * dy + rx * ex + ry * ey;
            if df.abs() < 1e-300 {
                break;
            }
            let next = (t - f / df).clamp(lo, hi);
            if (next - t).abs() < 1e-15 {
                break;
            }
            t = next;
        }
        if distance(t) > distance(best) {
            t = best;
        }
        (t, distance(t))
    }

    /// ∮ terms of Green's theorem along the spline: ∫ x dy, ∫ x²/2 dy and
    /// ∫ y²/2 dx. Exact, since each integrand is a polynomial in t of degree
    /// under 32.
    fn green(&self) -> [f64; 3] {
        let d = self.hodograph();
        let terms = |f: &dyn Fn([f64; 2], [f64; 2]) -> f64| {
            self.integrate(1, |t| f(self.point(t), d.point(t)))
        };
        [
            terms(&|[x, _], [_, dy]| x * dy),
            terms(&|[x, _], [_, dy]| x * x / 2.0 * dy),
            terms(&|[_, y], [dx, _]| y * y / 2.0 * dx),
        ]
    }
}

/// 16-point Gauss–Legendre nodes (positive half) and weights.
const GAUSS_16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_5),
    (0.617_876_244_402_643_7, 0.149_595_988_816_576_7),
    (0.755_404_408_355_003, 0.124_628_971_255_533_9),
    (0.865_631_202_387_831_7, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_89),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_09),
];

/// The probe's D: a 60 mm line and a cubic spline arching back to it, meeting
/// it at an angle at both ends.
fn d_arch() -> Spline {
    Spline::clamped(
        3,
        &[
            [60.0, 0.0],
            [66.0, 18.0],
            [50.0, 42.0],
            [24.0, 30.0],
            [8.0, 40.0],
            [-4.0, 14.0],
            [0.0, 0.0],
        ],
    )
}

/// The probe's other D, whose spline leaves the line and meets it again
/// tangentially.
fn tangent_arch() -> Spline {
    Spline::clamped(
        3,
        &[
            [60.0, 0.0],
            [75.0, 0.0],
            [75.0, 40.0],
            [-15.0, 40.0],
            [-15.0, 0.0],
            [0.0, 0.0],
        ],
    )
}

/// The probe's closed blob, clamped at its seam and tangent across it.
fn blob() -> Spline {
    Spline::clamped(
        3,
        &[
            [50.0, 0.0],
            [50.0, 30.0],
            [10.0, 45.0],
            [-40.0, 30.0],
            [-45.0, -10.0],
            [-10.0, -40.0],
            [50.0, -30.0],
            [50.0, 0.0],
        ],
    )
}

/// The base line of both Ds.
fn base() -> PlaneCurve {
    line([0.0, 0.0], [60.0, 0.0])
}

fn line(start: [f64; 2], end: [f64; 2]) -> PlaneCurve {
    PlaneCurve::Line { start, end }
}

/// The area and centroid of a region whose boundary is the base line (which
/// adds nothing to Green's integrals, lying on y = 0) and `arch`, closing it
/// anticlockwise.
fn d_area(arch: &Spline) -> (f64, [f64; 2]) {
    let [area, mx, my] = arch.green();
    (area, [mx / area, -my / area])
}

fn flat() -> Frame {
    Frame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap()
}

/// A plane tilted 30° about x, away from the origin.
fn tilted() -> Frame {
    let (sin, cos) = 30f64.to_radians().sin_cos();
    Frame::new([10.0, 20.0, 30.0], [1.0, 0.0, 0.0], [0.0, cos, sin]).unwrap()
}

/// Where a point in the frame's own (x, y, z) lies, and back.
fn world(frame: &Frame, [x, y]: [f64; 2], z: f64) -> [f64; 3] {
    frame.to_world([x, y, z])
}

fn local(frame: &Frame, p: [f64; 3]) -> [f64; 3] {
    let d = sub(p, frame.origin);
    [dot(d, frame.x_axis), dot(d, frame.y_axis), dot(d, frame.z_axis())]
}

// ---------------------------------------------------------------------------
// NURBS surfaces, evaluated here.

/// The non-zero basis functions of `degree` at `t` and the first one's
/// index (Piegl and Tiller's A2.2).
fn basis(knots: &[f64], degree: usize, t: f64) -> (usize, Vec<f64>) {
    let n = knots.len() - degree - 1;
    let t = t.clamp(knots[degree], knots[n]);
    let span = (degree..n)
        .rev()
        .find(|&i| knots[i] <= t && knots[i] < knots[i + 1])
        .unwrap_or(degree);
    let mut values = vec![0.0; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    values[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = values[r] / (right[r + 1] + left[j - r]);
            values[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        values[j] = saved;
    }
    (span - degree, values)
}

/// The point of a (rational) NURBS surface at (u, v).
fn net_point(net: &Net, u: f64, v: f64) -> [f64; 3] {
    let [pu, pv] = net.degree.map(|d| d as usize);
    let (iu, bu) = basis(&net.knots[0], pu, u);
    let (iv, bv) = basis(&net.knots[1], pv, v);
    let mut sum = [0.0; 4];
    for (a, nu) in bu.iter().enumerate() {
        for (b, nv) in bv.iter().enumerate() {
            let (i, j) = (iu + a, iv + b);
            let w = net.weights.as_ref().map_or(1.0, |w| w[i][j]);
            let p = net.poles[i][j];
            let f = nu * nv * w;
            sum[0] += f * p[0];
            sum[1] += f * p[1];
            sum[2] += f * p[2];
            sum[3] += f;
        }
    }
    [sum[0] / sum[3], sum[1] / sum[3], sum[2] / sum[3]]
}

/// That a face's surface is NURBS and its net, read back, is the surface:
/// it evaluates to the kernel's own points over the face's (u, v) box within
/// 1e-9, every point it gives is its own closest point within 1e-9, and the
/// tessellation's vertices on the face are its points at their (u, v) to f32
/// rounding. Returns points on the face, by the net, for the caller's checks.
fn assert_net(solid: &Brep, mesh: &Mesh, face: usize) -> Vec<[f64; 3]> {
    let faces = solid.faces().unwrap();
    assert!(matches!(faces[face].surface, Surface::Nurbs { .. }), "{:?}", faces[face].surface);
    let surface = FreeSurface::of_face(solid, face).unwrap();
    let net = surface.net().unwrap();
    let [u0, v0] = faces[face].uv_min;
    let [u1, v1] = faces[face].uv_max;
    let mut points = Vec::new();
    for a in 0..=12 {
        for b in 0..=12 {
            let (u, v) = (u0 + (u1 - u0) * a as f64 / 12.0, v0 + (v1 - v0) * b as f64 / 12.0);
            let ours = net_point(&net, u, v);
            let (theirs, _) = surface.evaluate(u, v).unwrap();
            assert_near("net against the kernel", ours, theirs, 1e-9);
            let closest = surface.closest_point(ours).unwrap();
            assert!(closest.distance <= 1e-9, "{closest:?}");
            points.push(ours);
        }
    }
    let mut vertices = 0;
    for (i, &id) in vertex_faces(mesh).iter().enumerate() {
        if id != Some(face) {
            continue;
        }
        vertices += 1;
        let [u, v] = mesh.surface_uvs[i].map(f64::from);
        let at = net_point(&net, u, v);
        let vertex = mesh.positions[i].map(f64::from);
        // The vertex and its (u, v) are both rounded to f32.
        assert_near("vertex against the net", at, vertex, 2e-3);
    }
    assert!(vertices > 0);
    points
}

/// The face each vertex of the mesh lies on, from its triangles.
fn vertex_faces(mesh: &Mesh) -> Vec<Option<usize>> {
    let mut faces = vec![None; mesh.positions.len()];
    for (t, tri) in mesh.triangles().enumerate() {
        for i in tri {
            faces[i as usize] = Some(mesh.face_ids[t] as usize);
        }
    }
    faces
}

// ---------------------------------------------------------------------------
// Planar faces.

/// What the tests look at in a sheet of planar faces and its prism.
fn d_sheet(k: &Kernel, frame: &Frame, arch: &Spline) -> Result<Brep, Error> {
    Brep::planar_faces(
        k,
        frame,
        &[vec![
            BoundaryCurve::labelled(base(), Label::new(1, "base")),
            BoundaryCurve::labelled(arch.curve(), Label::new(1, "arch")),
        ]],
    )
}

#[test]
fn a_spline_bounds_a_planar_face() {
    for frame in [flat(), tilted()] {
        let arch = d_arch();
        let (area, centroid) = d_area(&arch);
        on_kernel(move |k| {
            let sheet = d_sheet(k, &frame, &arch).unwrap();
            assert_eq!(sheet.face_count(), 1);
            let faces = sheet.faces().unwrap();
            assert_near("normal", faces[0].normal.unwrap(), frame.z_axis(), 1e-12);
            let measured = sheet.face_centroid(0).unwrap();
            assert!((measured.area - area).abs() < 1e-7 * area, "{} vs {area}", measured.area);
            assert_near("centroid", measured.centroid, world(&frame, centroid, 0.0), 1e-6);
            // Each curve's edge carries its label, and the arch's edge is the
            // spline: every point of it on the spline here.
            let edges = sheet.edges().unwrap();
            assert_eq!(edges.len(), 2);
            let arched = edges.iter().position(|e| e.label.as_ref().is_some_and(|l| l.role == "arch")).unwrap();
            assert!(edges.iter().any(|e| e.label.as_ref().is_some_and(|l| l.role == "base")));
            let curve = Curve::of_edge(&sheet, arched).unwrap();
            let [t0, t1] = curve.interval().unwrap();
            for i in 0..=50 {
                let (p, _) = curve.evaluate(t0 + (t1 - t0) * i as f64 / 50.0).unwrap();
                let [x, y, z] = local(&frame, p);
                assert!(z.abs() < 1e-9);
                assert!(arch.closest([x, y]).1 < 1e-9, "{:?}", arch.closest([x, y]));
            }
        });
    }
}

#[test]
fn a_spline_face_extrudes_and_its_wall_is_nurbs() {
    let arch = d_arch();
    let (area, _) = d_area(&arch);
    on_kernel(move |k| {
        let frame = tilted();
        let sheet = d_sheet(k, &frame, &arch).unwrap();
        let prism = sheet.sweep_face(0, 30.0).unwrap();
        let solid = &prism.prism;
        assert!(solid.is_manifold().unwrap());
        let volume = solid.volume(1e-9).unwrap();
        assert!((volume - area * 30.0).abs() < 1e-8 * volume, "{volume} vs {}", area * 30.0);
        assert_eq!(solid.face_count(), 4);
        // The arch's wall is one NURBS face standing on the spline.
        let edges = sheet.edges().unwrap();
        let arched = edges.iter().position(|e| e.label.as_ref().is_some_and(|l| l.role == "arch")).unwrap();
        let wall = prism
            .faces
            .iter()
            .position(|f| *f == SweptFace::Side { source_edge: Some(arched) })
            .unwrap();
        let mesh = solid.tessellate(&TessellationParams::default()).unwrap();
        for p in assert_net(solid, &mesh, wall) {
            let [x, y, z] = local(&frame, p);
            assert!((-1e-9..=30.0 + 1e-9).contains(&z));
            assert!(arch.closest([x, y]).1 < 1e-9);
        }
    });
}

#[test]
fn a_closed_spline_is_a_loop_on_its_own() {
    let blob = blob();
    let [area, ..] = blob.green();
    on_kernel(move |k| {
        // A plate with the blob as a hole, and the blob as a face.
        let rectangle = vec![
            line([-60.0, -50.0], [60.0, -50.0]),
            line([60.0, -50.0], [60.0, 50.0]),
            line([60.0, 50.0], [-60.0, 50.0]),
            line([-60.0, 50.0], [-60.0, -50.0]),
        ];
        let loops = [rectangle.clone(), vec![blob.curve()]]
            .map(|curves| curves.into_iter().map(|curve| BoundaryCurve { curve, label: None }).collect::<Vec<_>>());
        let plate = Brep::planar_faces(k, &flat(), &loops).unwrap();
        assert_eq!(plate.face_count(), 1);
        let plate_area = plate.face_centroid(0).unwrap().area;
        assert!((plate_area - (120.0 * 100.0 - area)).abs() < 1e-7 * plate_area, "{plate_area}");
        let face = Brep::planar_faces(k, &flat(), &loops[1..]).unwrap();
        let prism = face.sweep_face(0, 30.0).unwrap().prism;
        let volume = prism.volume(1e-9).unwrap();
        assert!((volume - 30.0 * area).abs() < 1e-8 * volume);
        // A closed spline can't be one curve of a longer loop.
        let crossed = Brep::planar_faces(
            k,
            &flat(),
            &[vec![
                BoundaryCurve { curve: blob.curve(), label: None },
                BoundaryCurve { curve: line([50.0, 0.0], [60.0, 0.0]), label: None },
            ]],
        );
        assert!(matches!(crossed, Err(Error::InvalidArgument { .. })), "{crossed:?}");
    });
}

#[test]
fn a_bad_spline_is_refused() {
    let arch = d_arch();
    on_kernel(move |k| {
        let mut wrong = arch.clone();
        for (what, change) in [
            ("knots that fall", Box::new(|s: &mut Spline| s.knots.reverse()) as Box<dyn Fn(&mut Spline)>),
            ("a knot short", Box::new(|s: &mut Spline| s.knots.truncate(s.knots.len() - 1))),
            ("unclamped knots", Box::new(|s: &mut Spline| s.knots[0] = -0.5)),
            ("a degree past its poles", Box::new(|s: &mut Spline| s.degree = 7)),
            ("a pole not finite", Box::new(|s: &mut Spline| s.poles[3][0] = f64::NAN)),
        ] {
            wrong.clone_from(&arch);
            change(&mut wrong);
            let refused = d_sheet(k, &flat(), &wrong);
            assert!(matches!(refused, Err(Error::InvalidArgument { .. })), "{what}: {refused:?}");
        }
    });
}

// ---------------------------------------------------------------------------
// The probe's cases: what the spec's probe found the kernel does with a
// spline profile.

#[test]
fn probe_drafts_lean_in_and_out_where_they_can() {
    // Drafted volumes are integrated to 1e-6 alone: the kernel took 96 s to
    // hold the tangent D's, drafted 15° out, to 1e-8.
    let (arch, tangent) = (d_arch(), tangent_arch());
    on_kernel(move |k| {
        let sheet = d_sheet(k, &flat(), &arch).unwrap();
        let straight = sheet.sweep_face(0, 30.0).unwrap().prism.volume(1e-6).unwrap();
        // In, it builds, and comes out smaller than the straight prism.
        for draft in [5.0, 15.0] {
            let drafted = sheet.sweep_face_tapered(0, 30.0, draft).unwrap().prism;
            assert!(drafted.is_manifold().unwrap());
            assert!(drafted.volume(1e-6).unwrap() < straight);
        }
        // Out, the corners where the spline meets its line at an angle stop
        // it, as an arc's do.
        for draft in [-5.0, -15.0] {
            assert!(sheet.sweep_face_tapered(0, 30.0, draft).is_err());
        }
        // Where they meet tangentially, and on a closed spline, it builds
        // either way.
        let tangent_sheet = d_sheet(k, &flat(), &tangent).unwrap();
        let loop_ = vec![BoundaryCurve { curve: blob().curve(), label: None }];
        let blob_sheet = Brep::planar_faces(k, &flat(), &[loop_]).unwrap();
        for sheet in [&tangent_sheet, &blob_sheet] {
            let straight = sheet.sweep_face(0, 30.0).unwrap().prism.volume(1e-6).unwrap();
            for draft in [5.0, -5.0, -15.0] {
                let drafted = sheet.sweep_face_tapered(0, 30.0, draft).unwrap().prism;
                assert!(drafted.is_manifold().unwrap());
                let volume = drafted.volume(1e-6).unwrap();
                assert!((volume < straight) == (draft > 0.0), "{draft}: {volume} vs {straight}");
            }
        }
    });
}

#[test]
fn probe_cuts_joins_and_imprints_exactly() {
    let arch = d_arch();
    let (area, _) = d_area(&arch);
    let blob = blob();
    on_kernel(move |k| {
        let sheet = d_sheet(k, &flat(), &arch).unwrap();
        let mut prism = sheet.sweep_face(0, 30.0).unwrap().prism;
        prism.translate([0.0, 0.0, -10.0]).unwrap();
        let block = Brep::cuboid(k, [-20.0, -20.0, 0.0], [100.0, 80.0, 10.0]).unwrap();
        let cut = block.try_clone().unwrap().difference(prism.try_clone().unwrap()).unwrap();
        assert!(cut.is_manifold().unwrap());
        let expected = 100.0 * 80.0 * 10.0 - area * 10.0;
        assert!((cut.volume(1e-9).unwrap() - expected).abs() < 1e-7 * expected);
        let joined = block.union(prism).unwrap();
        assert!(joined.is_manifold().unwrap());
        let expected = 100.0 * 80.0 * 10.0 + area * 20.0;
        assert!((joined.volume(1e-9).unwrap() - expected).abs() < 1e-7 * expected);
        // A closed spline imprinted onto a block's top splits it in two.
        let loop_ = vec![BoundaryCurve { curve: blob.curve(), label: None }];
        let mut tool = Brep::planar_faces(k, &flat(), &[loop_]).unwrap().sweep_face(0, 30.0).unwrap().prism;
        tool.translate([0.0, 0.0, -10.0]).unwrap();
        let block = Brep::cuboid(k, [-60.0, -50.0, 0.0], [120.0, 100.0, 10.0]).unwrap();
        let faces = block.face_count();
        let imprinted = block.imprint(tool).unwrap().brep;
        assert!(imprinted.is_manifold().unwrap());
        assert_eq!(imprinted.face_count(), faces + 2);
        assert!((imprinted.volume(1e-9).unwrap() - 120.0 * 100.0 * 10.0).abs() < 1e-6);
    });
}

#[test]
fn probe_sweeps_along_the_splines_edge() {
    // The probe swept along the sheet's spline edge (`sweep_profile_along`),
    // before a path could hold a spline.
    let arch = d_arch();
    let length = arch.length();
    on_kernel(move |k| {
        let sheet = d_sheet(k, &flat(), &arch).unwrap();
        let edges = sheet.edges().unwrap();
        let arched = edges.iter().position(|e| e.label.as_ref().is_some_and(|l| l.role == "arch")).unwrap();
        assert!((sheet.edge_length(arched).unwrap() - length).abs() < 1e-6 * length);
        for (profile, area) in [(vec![PlaneCurve::Circle { centre: [0.0, 0.0], radius: 3.0 }], PI * 9.0), (square(3.0), 36.0)] {
            let swept = sheet.sweep_profile_along(&[arched], &profile, 32, 1e-3).unwrap();
            assert!(swept.solid.is_manifold().unwrap());
            let volume = swept.solid.volume(1e-8).unwrap();
            assert!((volume - area * length).abs() < 1e-4 * volume, "{volume} vs {}", area * length);
        }
    });
}

// ---------------------------------------------------------------------------
// Revolving: a face turned from a spline reads back as NURBS.

#[test]
fn a_turned_spline_reads_as_nurbs() {
    let arch = d_arch();
    let (area, centroid) = d_area(&arch);
    for (frame, axis_point, axis, angle, radius) in [
        // About its line, a whole turn and three quarters.
        (flat(), [0.0, 0.0], [1.0, 0.0], 360.0, centroid[1]),
        (tilted(), [0.0, 0.0], [1.0, 0.0], 270.0, centroid[1]),
        // A whole turn about an axis off it, x = -10.
        (tilted(), [-10.0, 0.0], [0.0, 1.0], 360.0, centroid[0] + 10.0),
    ] {
        let arch = arch.clone();
        on_kernel(move |k| {
            let sheet = d_sheet(k, &frame, &arch).unwrap();
            let edges = sheet.edges().unwrap();
            let arched = edges.iter().position(|e| e.label.as_ref().is_some_and(|l| l.role == "arch")).unwrap();
            let base = world(&frame, axis_point, 0.0);
            let direction = sub(world(&frame, [axis_point[0] + axis[0], axis_point[1] + axis[1]], 0.0), base);
            let turned = sheet.revolve_face(0, base, direction, angle).unwrap();
            let solid = &turned.prism;
            assert!(solid.is_manifold().unwrap());
            // The kernel can't hold a turned NURBS side's volume to 1e-8
            // (`SM_ERR`), as `transform_bodies.rs` finds of the drafted D.
            let volume = solid.volume(1e-7).unwrap();
            let pappus = area * radius * angle.to_radians();
            assert!((volume - pappus).abs() < 1e-6 * pappus, "{volume} vs {pappus}");
            let wall = turned
                .faces
                .iter()
                .position(|f| *f == SweptFace::Side { source_edge: Some(arched) })
                .unwrap();
            let mesh = solid.tessellate(&TessellationParams::default()).unwrap();
            // Every point of it turned from the spline: its height along the
            // axis and its distance from it a point of the spline, in the
            // spline's own plane measured from the axis.
            for p in assert_net(solid, &mesh, wall) {
                let d = sub(p, base);
                let along = dot(d, direction);
                let out = norm(sub(d, scale(direction, along)));
                // The spline's points in the same terms: along the axis from
                // its point, and out from it.
                let (_, gap) = closest_turned(&arch, axis_point, axis, along, out);
                assert!(gap < 1e-8, "{gap}");
            }
        });
    }
}

/// How far (along, out) is from the spline measured from an axis in its
/// plane through `point` along unit `axis`: the nearest of its points, taken
/// to the axis's terms, by sampling and bisecting on the spline.
fn closest_turned(arch: &Spline, point: [f64; 2], axis: [f64; 2], along: f64, out: f64) -> (f64, f64) {
    let terms = |t: f64| {
        let [x, y] = arch.point(t);
        let (dx, dy) = (x - point[0], y - point[1]);
        let a = dx * axis[0] + dy * axis[1];
        let o = (dx * axis[1] - dy * axis[0]).abs();
        (a, o)
    };
    let gap = |t: f64| {
        let (a, o) = terms(t);
        (a - along).hypot(o - out)
    };
    let (lo, hi) = arch.domain();
    let samples = 4000;
    let mut best = lo;
    for k in 0..=samples {
        let t = lo + (hi - lo) * k as f64 / samples as f64;
        if gap(t) < gap(best) {
            best = t;
        }
    }
    // Golden-section about the best sample.
    let step = (hi - lo) / samples as f64;
    let (mut a, mut b) = ((best - step).max(lo), (best + step).min(hi));
    let ratio = (5f64.sqrt() - 1.0) / 2.0;
    for _ in 0..200 {
        let c = b - ratio * (b - a);
        let d = a + ratio * (b - a);
        if gap(c) < gap(d) {
            b = d;
        } else {
            a = c;
        }
    }
    let t = (a + b) / 2.0;
    (t, gap(t))
}

// ---------------------------------------------------------------------------
// Pipes and profiles along splines.

/// An S-shaped rail: a cubic leaving along +x from the origin and arriving
/// along +x at (90, 30).
fn rail() -> Spline {
    Spline::clamped(
        3,
        &[[0.0, 0.0], [30.0, 0.0], [45.0, 15.0], [60.0, 30.0], [90.0, 30.0]],
    )
}

/// A square profile of half-side `h`, round its middle.
fn square(h: f64) -> Vec<PlaneCurve> {
    vec![
        line([-h, -h], [h, -h]),
        line([h, -h], [h, h]),
        line([h, h], [-h, h]),
        line([-h, h], [-h, -h]),
    ]
}

/// What the tests look at in a pipe.
struct Piped {
    roles: Vec<PipeFace>,
    surfaces: Vec<Surface>,
    manifold: bool,
    volume: f64,
    mesh: Mesh,
}

fn pipe_along(frame: Frame, path: Vec<PlaneCurve>, radius: f64, inner: Option<f64>) -> Result<Piped, Error> {
    on_kernel(move |k| {
        let pipe = Brep::pipe(k, &frame, &path, radius, inner)?;
        let solid = &pipe.solid;
        let surfaces: Vec<Surface> = solid.faces()?.iter().map(|f| f.surface.clone()).collect();
        let mesh = solid.tessellate(&TessellationParams::default())?;
        for (face, surface) in surfaces.iter().enumerate() {
            if matches!(surface, Surface::Nurbs { .. }) {
                assert_net(solid, &mesh, face);
            }
        }
        Ok(Piped {
            roles: pipe.faces.clone(),
            surfaces,
            manifold: solid.is_manifold()?,
            volume: solid.volume(1e-8)?,
            mesh,
        })
    })
}

/// That a pipe is closed, has one outer wall (and inner, when hollow) per
/// path curve and nothing else but its caps, holds `volume` within `relative`
/// and winds its mesh outward round it.
fn assert_piped(piped: &Piped, curves: usize, hollow: bool, closed: bool, volume: f64, relative: f64) {
    assert!(piped.manifold);
    assert!((piped.volume - volume).abs() <= relative * volume, "{} vs {volume}", piped.volume);
    for segment in 0..curves {
        let count = |role: PipeFace| piped.roles.iter().filter(|r| **r == role).count();
        assert_eq!(count(PipeFace::Outer { segment: Some(segment) }), 1, "{:?}", piped.roles);
        if hollow {
            assert_eq!(count(PipeFace::Inner { segment: Some(segment) }), 1, "{:?}", piped.roles);
        }
    }
    let walls = curves * if hollow { 2 } else { 1 };
    let caps = if closed { 0 } else { 2 };
    assert_eq!(piped.roles.len(), walls + caps, "{:?}", piped.roles);
    assert_mesh(&piped.mesh, piped.volume);
}

/// That every vertex of the wall swept along `spline` in `frame` lies
/// `radius` from it, within `tolerance`: the sweep is approximated to 1e-5,
/// and the mesh's vertices are f32, 1e-5 apart at 100 mm.
fn assert_round(piped: &Piped, wall: usize, frame: &Frame, spline: &Spline, radius: f64, tolerance: f64) {
    assert!(matches!(piped.surfaces[wall], Surface::Nurbs { .. }), "{:?}", piped.surfaces[wall]);
    let mut worst: f64 = 0.0;
    for (i, &face) in vertex_faces(&piped.mesh).iter().enumerate() {
        if face != Some(wall) {
            continue;
        }
        let [x, y, z] = local(frame, piped.mesh.positions[i].map(f64::from));
        let (_, d) = spline.closest([x, y]);
        worst = worst.max((d.hypot(z) - radius).abs());
    }
    assert!(worst <= tolerance, "{worst}");
}

#[test]
fn a_pipe_follows_a_spline() {
    let rail = rail();
    let length = rail.length();
    for frame in [flat(), tilted()] {
        let piped = pipe_along(frame, vec![rail.curve()], 3.0, None).unwrap();
        assert_piped(&piped, 1, false, false, PI * 9.0 * length, 1e-6);
        let wall = piped.roles.iter().position(|r| *r == PipeFace::Outer { segment: Some(0) }).unwrap();
        assert_round(&piped, wall, &frame, &rail, 3.0, 2e-5);
        // Hollow, its bore runs out through both caps.
        let hollow = pipe_along(frame, vec![rail.curve()], 3.0, Some(1.5)).unwrap();
        assert_piped(&hollow, 1, true, false, PI * (9.0 - 2.25) * length, 1e-6);
        let bore = hollow.roles.iter().position(|r| *r == PipeFace::Inner { segment: Some(0) }).unwrap();
        assert_round(&hollow, bore, &frame, &rail, 1.5, 2e-5);
    }
}

#[test]
fn a_spline_meets_lines_and_arcs_face_to_face() {
    // A line in along +x, the rail, a quarter turn anticlockwise, and a
    // line on up.
    let rail = rail();
    let path = vec![
        line([-40.0, 0.0], [0.0, 0.0]),
        rail.curve(),
        PlaneCurve::Arc { centre: [90.0, 50.0], radius: 20.0, start_deg: -90.0, sweep_deg: 90.0 },
        line([110.0, 50.0], [110.0, 90.0]),
    ];
    let length = 40.0 + rail.length() + PI * 10.0 + 40.0;
    for frame in [flat(), tilted()] {
        let piped = pipe_along(frame, path.clone(), 3.0, None).unwrap();
        assert_piped(&piped, 4, false, false, PI * 9.0 * length, 1e-6);
        let wall = piped.roles.iter().position(|r| *r == PipeFace::Outer { segment: Some(1) }).unwrap();
        assert_round(&piped, wall, &frame, &rail, 3.0, 2e-5);
        // The lines' walls are still cylinders, the arc's a torus.
        for (segment, analytic) in [(0, "cylinder"), (2, "torus"), (3, "cylinder")] {
            let wall = piped.roles.iter().position(|r| *r == PipeFace::Outer { segment: Some(segment) }).unwrap();
            let kind = match piped.surfaces[wall] {
                Surface::Cylinder { .. } => "cylinder",
                Surface::Torus { .. } => "torus",
                _ => "other",
            };
            assert_eq!(kind, analytic);
        }
    }
}

#[test]
fn a_pipe_rings_a_closed_spline() {
    let blob = blob();
    let length = blob.length();
    for frame in [flat(), tilted()] {
        let piped = pipe_along(frame, vec![blob.curve()], 3.0, None).unwrap();
        assert_piped(&piped, 1, false, true, PI * 9.0 * length, 1e-6);
        assert_round(&piped, 0, &frame, &blob, 3.0, 2e-5);
        // A closed path's bore would be sealed inside.
        assert!(pipe_along(frame, vec![blob.curve()], 3.0, Some(1.5)).is_err());
    }
}

#[test]
fn a_spline_path_is_refused_where_a_pipe_cant_follow() {
    let rail = rail();
    let refused = |path: Vec<PlaneCurve>, radius: f64| match pipe_along(flat(), path, radius, None) {
        Err(Error::InvalidArgument { reason, .. }) => reason,
        other => panic!("{:?}", other.map(|p| p.roles)),
    };
    // Meeting a line at an angle.
    let reason = refused(vec![line([0.0, -40.0], [0.0, 0.0]), rail.curve()], 3.0);
    assert!(reason.contains("spline"), "{reason}");
    // A pipe wider than the rail's tightest bend.
    let reason = refused(vec![rail.curve()], 1.01 * rail.tightest());
    assert!(reason.contains("bend"), "{reason}");
    // A corner inside it: a degree 1 spline that turns.
    let kinked = PlaneCurve::Spline {
        degree: 1,
        poles: vec![[0.0, 0.0], [40.0, 0.0], [40.0, 40.0]],
        knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
    };
    let reason = refused(vec![kinked], 3.0);
    assert!(reason.contains("corner"), "{reason}");
}

/// What the tests look at in a profile swept along a path.
struct Swept {
    roles: Vec<ProfileFace>,
    manifold: bool,
    volume: f64,
    mesh: Mesh,
}

fn sweep_along(frame: Frame, path: Vec<PlaneCurve>, profile: Vec<PlaneCurve>, ends: [PathEnd; 2]) -> Result<Swept, Error> {
    on_kernel(move |k| {
        let swept = Brep::sweep_profile(k, &frame, &path, &profile, ends)?;
        let mesh = swept.solid.tessellate(&TessellationParams::default())?;
        for (face, geometry) in swept.solid.faces()?.iter().enumerate() {
            if matches!(geometry.surface, Surface::Nurbs { .. }) {
                assert_net(&swept.solid, &mesh, face);
            }
        }
        Ok(Swept {
            roles: swept.faces.clone(),
            manifold: swept.solid.is_manifold()?,
            volume: swept.solid.volume(1e-8)?,
            mesh,
        })
    })
}

#[test]
fn a_profile_follows_a_spline() {
    let rail = rail();
    let path = vec![line([-40.0, 0.0], [0.0, 0.0]), rail.curve(), line([90.0, 30.0], [130.0, 30.0])];
    let length = 80.0 + rail.length();
    for frame in [flat(), tilted()] {
        // The probe's 6 mm square bar.
        let swept = sweep_along(frame, path.clone(), square(3.0), [PathEnd::Square; 2]).unwrap();
        assert!(swept.manifold);
        assert!((swept.volume - 36.0 * length).abs() < 1e-6 * swept.volume, "{} vs {}", swept.volume, 36.0 * length);
        assert_mesh(&swept.mesh, swept.volume);
        // The walls along the spline are traced to their profile segments.
        // The level sides, the bar's top and bottom, are one plane each all
        // along, merged by the unions into the first line's, as a level side
        // along lines and arcs is.
        let count = |side: ProfileFace| swept.roles.iter().filter(|r| **r == side).count();
        for segment in [1, 3] {
            assert_eq!(count(ProfileFace::Side { path: Some(1), profile: Some(segment) }), 1, "{:?}", swept.roles);
        }
        for segment in [0, 2] {
            assert_eq!(count(ProfileFace::Side { path: Some(0), profile: Some(segment) }), 1, "{:?}", swept.roles);
        }
        assert_eq!(swept.roles.len(), 3 * 2 + 2 + 2, "{:?}", swept.roles);
        assert!(swept.roles.contains(&ProfileFace::Start { profile: None }));
        assert!(swept.roles.contains(&ProfileFace::End { profile: None }));
        // Round ends where the path is the spline alone: each end the
        // profile's left half turned half a turn, half a cylinder 6 tall.
        let rounded = sweep_along(frame, vec![rail.curve()], square(3.0), [PathEnd::Round; 2]).unwrap();
        assert!(rounded.manifold);
        let expected = 36.0 * rail.length() + PI * 9.0 * 6.0;
        assert!((rounded.volume - expected).abs() < 1e-6 * expected, "{} vs {expected}", rounded.volume);
    }
}

#[test]
fn a_profile_is_refused_a_spline_or_a_bend_it_cant_turn() {
    let rail = rail();
    // A spline in a profile.
    let profile = vec![line([-3.0, -3.0], [3.0, -3.0]), PlaneCurve::Spline {
        degree: 2,
        poles: vec![[3.0, -3.0], [3.0, 3.0], [-3.0, -3.0]],
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    }];
    assert!(matches!(
        sweep_along(flat(), vec![rail.curve()], profile, [PathEnd::Square; 2]),
        Err(Error::InvalidArgument { .. })
    ));
    // A profile reaching past the rail's tightest bend.
    assert!(matches!(
        sweep_along(flat(), vec![rail.curve()], square(1.01 * rail.tightest()), [PathEnd::Square; 2]),
        Err(Error::InvalidArgument { .. })
    ));
}

// ---------------------------------------------------------------------------
// Lofts, the probe's last case: unlike sections stitch into a solid.

/// A closed section through `corners`, as a degree 1 NURBS from and back to
/// the first, at chord-length knots.
fn polygon(k: &Kernel, corners: &[[f64; 3]]) -> Curve {
    let mut poles = corners.to_vec();
    poles.push(corners[0]);
    let mut lengths = vec![0.0];
    for w in poles.windows(2) {
        lengths.push(lengths.last().unwrap() + norm(sub(w[1], w[0])));
    }
    let total = *lengths.last().unwrap();
    let mut knots = vec![0.0];
    knots.extend(lengths.iter().map(|l| l / total));
    knots.push(1.0);
    Curve::nurbs(k, 1, &poles, None, &knots).unwrap()
}

#[test]
fn probe_lofts_unlike_sections_into_a_solid() {
    on_kernel(|k| {
        // A 60 × 40 rectangle at z 0 to a circle r 20 at z 80, straight and
        // through an ellipse 26 × 18 at z 40.
        let rect = polygon(k, &[[30.0, 0.0, 0.0], [30.0, 20.0, 0.0], [-30.0, 20.0, 0.0], [-30.0, -20.0, 0.0], [30.0, -20.0, 0.0]]);
        let mid = Curve::ellipse(k, [0.0, 0.0, 40.0], [26.0, 0.0, 0.0], [0.0, 18.0, 0.0]).unwrap();
        let top = Curve::ellipse(k, [0.0, 0.0, 80.0], [20.0, 0.0, 0.0], [0.0, 20.0, 0.0]).unwrap();
        for sections in [vec![&rect, &top], vec![&rect, &mid, &top]] {
            let skin = FreeSurface::loft(k, &sections, 3).unwrap();
            let mut sheets = Sheets::new(k).unwrap();
            // The skin comes in pieces, split at the rectangle's corners and
            // the circle's quarters.
            let pieces = sheets.add(&skin, None).unwrap();
            assert_eq!(pieces, 8);
            sheets.add_planar(&[sections[0]], None).unwrap();
            sheets.add_planar(&[sections[sections.len() - 1]], None).unwrap();
            let (solid, closed) = sheets.stitch(0.0).unwrap();
            assert!(closed && solid.is_manifold().unwrap());
            let volume = solid.volume(1e-6).unwrap();
            // Between the prism on the rectangle and the cylinder on the
            // circle.
            assert!(volume > PI * 400.0 * 80.0 && volume < 2400.0 * 80.0, "{volume}");
        }
    });
}

// ---------------------------------------------------------------------------

fn assert_mesh(mesh: &Mesh, volume: f64) {
    assert!(mesh.triangle_count() > 0);
    let widen = |v: [f32; 3]| v.map(f64::from);
    let mut enclosed = 0.0;
    for tri in mesh.triangles() {
        let [a, b, c] = tri.map(|i| widen(mesh.positions[i as usize]));
        let normal = tri.iter().fold([0.0; 3], |n, &i| add(n, widen(mesh.normals[i as usize])));
        assert!(dot(cross(sub(b, a), sub(c, a)), normal) > 0.0);
        enclosed += dot(a, cross(b, c)) / 6.0;
    }
    assert!((enclosed - volume).abs() < 0.01 * volume, "{enclosed} vs {volume}");
}

fn assert_near(what: &str, actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    let gap = norm(sub(actual, expected));
    assert!(gap <= tolerance, "{what}: {actual:?} isn't within {tolerance} of {expected:?} ({gap})");
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
