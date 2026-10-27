//! Profiles swept along paths of lines and arcs: a panel line's groove. Every
//! side must trace back to the path curve and profile segment it was swept
//! from, the solid must be closed, and a profile symmetric across the path
//! must take its area times the path's length, round its bends and mitres
//! too (Pappus, with its centroid on the path).
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, Error, Frame, Kernel, Mesh, PathEnd, PlaneCurve, ProfileFace, Surface, TessellationParams,
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

/// The grooves' width, depth and the clearance their walls run on above the
/// surface.
const W: f64 = 4.0;
const D: f64 = 2.0;
const C: f64 = 1.0;

fn flat() -> Frame {
    Frame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap()
}

/// A plane tilted 30° about x, away from the origin.
fn tilted() -> Frame {
    let (sin, cos) = 30f64.to_radians().sin_cos();
    Frame::new([10.0, 20.0, 30.0], [1.0, 0.0, 0.0], [0.0, cos, sin]).unwrap()
}

fn line(start: [f64; 2], end: [f64; 2]) -> PlaneCurve {
    PlaneCurve::Line { start, end }
}

/// A closed polyline through `points`.
fn polygon(points: &[[f64; 2]]) -> Vec<PlaneCurve> {
    (0..points.len())
        .map(|i| line(points[i], points[(i + 1) % points.len()]))
        .collect()
}

/// A square groove: floor, right wall, top, left wall.
fn square() -> (Vec<PlaneCurve>, f64) {
    let h = W / 2.0;
    (polygon(&[[-h, -D], [h, -D], [h, C], [-h, C]]), W * (D + C))
}

/// A V: its two slopes, then its walls and top above the surface.
fn vee() -> (Vec<PlaneCurve>, f64) {
    let h = W / 2.0;
    (
        polygon(&[[-h, 0.0], [0.0, -D], [h, 0.0], [h, C], [-h, C]]),
        W * D / 2.0 + W * C,
    )
}

/// A round floor, a half circle as deep as half the width, then walls and
/// top.
fn round() -> (Vec<PlaneCurve>, f64) {
    let h = W / 2.0;
    let mut curves = vec![PlaneCurve::Arc {
        centre: [0.0, 0.0],
        radius: h,
        start_deg: 180.0,
        sweep_deg: 180.0,
    }];
    curves.extend([line([h, 0.0], [h, C]), line([h, C], [-h, C]), line([-h, C], [-h, 0.0])]);
    (curves, PI * h * h / 2.0 + W * C)
}

/// A stepped groove: a shoulder half as deep and a quarter as wide each
/// side, and the narrow floor below.
fn stepped() -> (Vec<PlaneCurve>, f64) {
    let (h, q) = (W / 2.0, W / 4.0);
    let points = [
        [-h, -D / 2.0],
        [-q, -D / 2.0],
        [-q, -D],
        [q, -D],
        [q, -D / 2.0],
        [h, -D / 2.0],
        [h, C],
        [-h, C],
    ];
    (polygon(&points), W * (D / 2.0 + C) + 2.0 * q * D / 2.0)
}

/// A rounded rectangle, anticlockwise from the middle of its bottom: lines
/// and quarter turns, 60 by 40 with corners of 6.
fn rounded_rectangle() -> (Vec<PlaneCurve>, f64) {
    let (a, b, r) = (30.0, 20.0, 6.0);
    let arc = |centre: [f64; 2], start_deg: f64| PlaneCurve::Arc {
        centre,
        radius: r,
        start_deg,
        sweep_deg: 90.0,
    };
    let path = vec![
        line([0.0, -b], [a - r, -b]),
        arc([a - r, -b + r], -90.0),
        line([a, -b + r], [a, b - r]),
        arc([a - r, b - r], 0.0),
        line([a - r, b], [-a + r, b]),
        arc([-a + r, b - r], 90.0),
        line([-a, b - r], [-a, -b + r]),
        arc([-a + r, -b + r], 180.0),
        line([-a + r, -b], [0.0, -b]),
    ];
    (path, 4.0 * (a + b - 2.0 * r) + 2.0 * PI * r)
}

/// What the tests look at in a swept profile.
#[derive(Debug)]
struct Swept {
    roles: Vec<ProfileFace>,
    surfaces: Vec<Surface>,
    labels_absent: bool,
    manifold: bool,
    volume: f64,
    mesh: Mesh,
}

impl Swept {
    /// The faces with `role`, in face order.
    fn faces(&self, role: ProfileFace) -> Vec<usize> {
        (0..self.roles.len())
            .filter(|&i| self.roles[i] == role)
            .collect()
    }

    /// The one face with `role`.
    fn face(&self, role: ProfileFace) -> usize {
        let faces = self.faces(role);
        assert_eq!(faces.len(), 1, "{role:?}: {:?}", self.roles);
        faces[0]
    }
}

fn sweep(
    frame: Frame,
    path: Vec<PlaneCurve>,
    profile: Vec<PlaneCurve>,
    ends: [PathEnd; 2],
) -> Result<Swept, Error> {
    on_kernel(move |k| {
        let swept = Brep::sweep_profile(k, &frame, &path, &profile, ends)?;
        let solid = &swept.solid;
        let faces = solid.faces().unwrap();
        Ok(Swept {
            roles: swept.faces.clone(),
            surfaces: faces.iter().map(|f| f.surface.clone()).collect(),
            labels_absent: faces.iter().all(|f| f.label.is_none()),
            manifold: solid.is_manifold().unwrap(),
            volume: solid.volume(1e-9).unwrap(),
            mesh: solid.tessellate(&TessellationParams::default()).unwrap(),
        })
    })
}

/// What every swept profile must be: closed, unlabelled, and a mesh that
/// winds outwards round its volume, with a side for each path curve and
/// profile segment. The union merges sides on one surface, which leaves a
/// level segment's side (a floor, a step or the top, all on one plane) one
/// face for the whole path, and a line's sides in the sides of the line
/// before it when it carries straight on: those come from the first curve.
fn assert_swept(swept: &Swept, path: &[PlaneCurve], profile: &[PlaneCurve], volume: f64) {
    assert!(swept.manifold);
    assert!(close(swept.volume, volume), "{} vs {volume}", swept.volume);
    assert!(swept.labels_absent);
    let sides: BTreeSet<(usize, usize)> = swept
        .roles
        .iter()
        .filter_map(|role| match *role {
            ProfileFace::Side { path, profile } => Some((path.unwrap(), profile.unwrap())),
            _ => None,
        })
        .collect();
    let level = |s: usize| matches!(profile[s], PlaneCurve::Line { start, end } if start[1] == end[1]);
    // A line carrying straight on from or into an earlier one.
    let closed = {
        let ([x1, y1], [x2, y2]) = (path[path.len() - 1].end(), path[0].start());
        (x1 - x2).hypot(y1 - y2) < 1e-6
    };
    let collinear = |a: PlaneCurve, b: PlaneCurve| match (a, b) {
        (PlaneCurve::Line { start: p, end: q }, PlaneCurve::Line { start: r, end: t }) => [r, t]
            .iter()
            .all(|u| ((q[0] - p[0]) * (u[1] - p[1]) - (q[1] - p[1]) * (u[0] - p[0])).abs() < 1e-9),
        _ => false,
    };
    let straight_on = |c: usize| {
        let n = path.len();
        let neighbours = [(c + n - 1) % n, (c + 1) % n];
        neighbours.iter().any(|&m| {
            let adjacent = closed || m.abs_diff(c) == 1;
            m < c && adjacent && collinear(path[m], path[c])
        })
    };
    for s in 0..profile.len() {
        let curves: Vec<usize> = sides.iter().filter(|side| side.1 == s).map(|side| side.0).collect();
        let expected: Vec<usize> = if level(s) {
            vec![0]
        } else {
            (0..path.len()).filter(|&c| !straight_on(c)).collect()
        };
        assert_eq!(curves, expected, "segment {s}: {:?}", swept.roles);
    }
    assert_mesh(&swept.mesh, swept.volume);
}

#[test]
fn a_straight_groove_is_its_profile_times_its_length() {
    for frame in [flat(), tilted()] {
        for (profile, area) in [square(), vee(), round(), stepped()] {
            let segments = profile.len();
            let path = vec![line([0.0, 0.0], [50.0, 0.0])];
            let swept = sweep(frame, path.clone(), profile.clone(), [PathEnd::Square; 2]).unwrap();
            assert_swept(&swept, &path, &profile, area * 50.0);
            // One side per segment, and the two caps.
            assert_eq!(swept.roles.len(), segments + 2, "{:?}", swept.roles);
            swept.face(ProfileFace::Start { profile: None });
            swept.face(ProfileFace::End { profile: None });
        }
    }
}

#[test]
fn a_square_grooves_floor_lies_its_depth_down() {
    let frame = tilted();
    let swept = sweep(frame, vec![line([0.0, 0.0], [50.0, 0.0])], square().0, [PathEnd::Square; 2]).unwrap();
    let floor = swept.face(ProfileFace::Side {
        path: Some(0),
        profile: Some(0),
    });
    let Surface::Plane { origin, normal, .. } = swept.surfaces[floor] else {
        panic!("{:?} isn't a plane", swept.surfaces[floor]);
    };
    let z = frame.z_axis();
    assert!(1.0 - dot(normal, z).abs() < 1e-12, "{normal:?}");
    assert!((dot(sub(origin, frame.origin), z) + D).abs() < 1e-9);
}

#[test]
fn bends_are_turned_and_corners_mitred() {
    for (profile, area) in [square(), vee(), round(), stepped()] {
        // A rounded rectangle: its bends' sides are planes, cylinders and
        // cones, each analytic.
        let (path, length) = rounded_rectangle();
        let swept = sweep(flat(), path.clone(), profile.clone(), [PathEnd::Square; 2]).unwrap();
        assert_swept(&swept, &path, &profile, area * length);
        assert!(!swept.roles.iter().any(|r| !matches!(r, ProfileFace::Side { .. })));
        for (i, role) in swept.roles.iter().enumerate() {
            let analytic = matches!(
                swept.surfaces[i],
                Surface::Plane { .. } | Surface::Cylinder { .. } | Surface::Cone { .. } | Surface::Torus { .. }
            );
            // SMLib's torus recogniser misses a half circle turned round a
            // bend, the round floor's, which stays an exact surface of
            // revolution.
            let half_torus = matches!(role, ProfileFace::Side { profile: Some(0), .. })
                && matches!(profile[0], PlaneCurve::Arc { .. })
                && matches!(swept.surfaces[i], Surface::Other { kernel_type: 15020 });
            assert!(analytic || half_torus, "{role:?}: {:?}", swept.surfaces[i]);
        }
        // A frame of four lines, mitred at each corner, and an open zig-zag.
        let frame = polygon(&[[0.0, 0.0], [60.0, 0.0], [60.0, 40.0], [0.0, 40.0]]);
        let swept = sweep(tilted(), frame.clone(), profile.clone(), [PathEnd::Square; 2]).unwrap();
        assert_swept(&swept, &frame, &profile, area * 200.0);
        let zigzag = vec![line([0.0, 0.0], [30.0, 0.0]), line([30.0, 0.0], [50.0, 25.0]), line([50.0, 25.0], [80.0, 5.0])];
        let length = 30.0 + 20f64.hypot(25.0) + 30f64.hypot(20.0);
        let swept = sweep(flat(), zigzag.clone(), profile.clone(), [PathEnd::Square; 2]).unwrap();
        assert_swept(&swept, &zigzag, &profile, area * length);
    }
}

#[test]
fn round_ends_turn_half_the_profile() {
    // Each end is the left half turned half a turn: the two together the
    // half turned a whole one, by Pappus.
    let (profile, area) = square();
    let half = (W / 2.0) * (D + C);
    let ends = 2.0 * PI * half * (W / 4.0);
    let path = vec![line([0.0, 0.0], [50.0, 0.0])];
    let swept = sweep(tilted(), path.clone(), profile.clone(), [PathEnd::Round; 2]).unwrap();
    assert_swept(&swept, &path, &profile, area * 50.0 + ends);
    for end in [ProfileFace::Start { profile: Some(1) }, ProfileFace::End { profile: Some(1) }] {
        // The wall turns into a half cylinder as wide as the groove.
        let wall = swept.face(end);
        let Surface::Cylinder { radius, .. } = swept.surfaces[wall] else {
            panic!("{:?} isn't a cylinder", swept.surfaces[wall]);
        };
        assert!((radius - W / 2.0).abs() < 1e-9);
    }
    // A V's round end is a cone; one end round and the other square.
    let (profile, area) = vee();
    let cone = PI * (W / 2.0).powi(2) * D / 3.0 / 2.0 + PI * (W / 2.0).powi(2) * C / 2.0;
    let path = vec![line([0.0, 0.0], [30.0, 0.0])];
    let swept = sweep(flat(), path.clone(), profile.clone(), [PathEnd::Square, PathEnd::Round]).unwrap();
    assert_swept(&swept, &path, &profile, area * 30.0 + cone);
    swept.face(ProfileFace::Start { profile: None });
    let slope = swept.face(ProfileFace::End { profile: Some(1) });
    assert!(matches!(swept.surfaces[slope], Surface::Cone { .. }), "{:?}", swept.surfaces[slope]);
}

#[test]
fn a_circle_is_a_pipe() {
    let disc = vec![PlaneCurve::Circle {
        centre: [0.0, 0.0],
        radius: 2.0,
    }];
    let (path, length) = rounded_rectangle();
    let swept = sweep(flat(), path.clone(), disc.clone(), [PathEnd::Square; 2]).unwrap();
    assert_swept(&swept, &path, &disc, PI * 4.0 * length);
}

#[test]
fn refuses_what_it_cant_sweep() {
    let reason = |path: Vec<PlaneCurve>, profile: Vec<PlaneCurve>, ends: [PathEnd; 2]| match sweep(flat(), path, profile, ends) {
        Err(Error::InvalidArgument { reason, .. }) => reason,
        other => panic!("{other:?}"),
    };
    let straight = || vec![line([0.0, 0.0], [50.0, 0.0])];
    // A bend of 2 can't take a groove reaching 2 towards its centre.
    let tight = vec![
        line([0.0, 0.0], [10.0, 0.0]),
        PlaneCurve::Arc {
            centre: [10.0, 2.0],
            radius: 2.0,
            start_deg: -90.0,
            sweep_deg: 90.0,
        },
    ];
    assert_eq!(
        reason(tight, square().0, [PathEnd::Square; 2]),
        "a bend's radius must be more than the profile reaches towards its centre"
    );
    // Wholly to the left of the path, it has no half to turn round.
    let aside = polygon(&[[1.0, -D], [3.0, -D], [3.0, C], [1.0, C]]);
    assert_eq!(
        reason(straight(), aside.clone(), [PathEnd::Round, PathEnd::Square]),
        "a round end needs a profile reaching both sides of the path"
    );
    assert!(sweep(flat(), straight(), aside, [PathEnd::Square; 2]).is_ok());
    let mut open = square().0;
    open.pop();
    assert_eq!(
        reason(straight(), open, [PathEnd::Square; 2]),
        "the profile's curves must run head to tail and close"
    );
    let corner = vec![
        line([0.0, 0.0], [10.0, 0.0]),
        PlaneCurve::Arc {
            centre: [20.0, 0.0],
            radius: 10.0,
            start_deg: 180.0,
            sweep_deg: 90.0,
        },
    ];
    assert_eq!(
        reason(corner, square().0, [PathEnd::Square; 2]),
        "an arc in the path must meet the curves either side of it tangentially"
    );
}

/// That the mesh's triangles wind with the kernel's normals and enclose the
/// solid's volume to 1%.
fn assert_mesh(mesh: &Mesh, volume: f64) {
    assert!(mesh.triangle_count() > 0);
    let widen = |v: [f32; 3]| v.map(f64::from);
    let mut enclosed = 0.0;
    for tri in mesh.triangles() {
        let [a, b, c] = tri.map(|i| widen(mesh.positions[i as usize]));
        let normal = tri
            .iter()
            .fold([0.0; 3], |n, &i| add(n, widen(mesh.normals[i as usize])));
        assert!(dot(cross(sub(b, a), sub(c, a)), normal) > 0.0);
        enclosed += dot(a, cross(b, c)) / 6.0;
    }
    assert!(
        (enclosed - volume).abs() < 0.01 * volume,
        "{enclosed} vs {volume}"
    );
}

/// Within the accuracy volumes are integrated to.
fn close(actual: f64, expected: f64) -> bool {
    (actual - expected).abs() <= 1e-6 * expected.abs()
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
