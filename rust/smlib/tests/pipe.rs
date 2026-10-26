//! Pipes swept along paths of lines and arcs: what a sketch curve's Pipe
//! needs. Every wall must trace back to the path curve it was swept along,
//! the solid must be closed, and its walls must be cylinders and tori.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{Brep, Error, Frame, Kernel, Mesh, PipeFace, PlaneCurve, Surface, TessellationParams};

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

/// The pipe's radius, and its inner radius when hollow.
const R: f64 = 5.0;
const INNER: f64 = 3.0;

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

/// The area of the pipe's cross-section.
fn section(inner: Option<f64>) -> f64 {
    PI * (R * R - inner.map_or(0.0, |r| r * r))
}

/// What the tests look at in a pipe.
#[derive(Debug)]
struct Piped {
    roles: Vec<PipeFace>,
    surfaces: Vec<Surface>,
    labels_absent: bool,
    centroids: Vec<[f64; 3]>,
    normals: Vec<Option<[f64; 3]>>,
    manifold: bool,
    volume: f64,
    mesh: Mesh,
}

impl Piped {
    /// The faces with each role, in face order.
    fn faces(&self, role: PipeFace) -> Vec<usize> {
        (0..self.roles.len())
            .filter(|&i| self.roles[i] == role)
            .collect()
    }

    /// The one face with `role`.
    fn face(&self, role: PipeFace) -> usize {
        let faces = self.faces(role);
        assert_eq!(faces.len(), 1, "{role:?}: {:?}", self.roles);
        faces[0]
    }
}

fn pipe(frame: Frame, path: Vec<PlaneCurve>, inner: Option<f64>) -> Piped {
    on_kernel(move |k| {
        let pipe = Brep::pipe(k, &frame, &path, R, inner).unwrap();
        let solid = &pipe.solid;
        let faces = solid.faces().unwrap();
        Piped {
            roles: pipe.faces.clone(),
            surfaces: faces.iter().map(|f| f.surface.clone()).collect(),
            labels_absent: faces.iter().all(|f| f.label.is_none()),
            centroids: (0..faces.len())
                .map(|i| solid.face_centroid(i).unwrap().centroid)
                .collect(),
            normals: faces.iter().map(|f| f.normal).collect(),
            manifold: solid.is_manifold().unwrap(),
            volume: solid.volume(1e-9).unwrap(),
            mesh: solid.tessellate(&TessellationParams::default()).unwrap(),
        }
    })
}

/// What every pipe must be: closed, unlabelled, one outer wall (and inner,
/// when hollow) per curve of the path and nothing else but its caps, and a
/// mesh that winds outwards round its volume.
fn assert_piped(piped: &Piped, curves: usize, inner: Option<f64>, closed: bool, volume: f64) {
    assert!(piped.manifold);
    assert!(close(piped.volume, volume), "{} vs {volume}", piped.volume);
    assert!(piped.labels_absent);
    for segment in 0..curves {
        piped.face(PipeFace::Outer {
            segment: Some(segment),
        });
        if inner.is_some() {
            piped.face(PipeFace::Inner {
                segment: Some(segment),
            });
        }
    }
    let walls = curves * if inner.is_some() { 2 } else { 1 };
    let caps = if closed { 0 } else { 2 };
    assert_eq!(piped.roles.len(), walls + caps, "{:?}", piped.roles);
    assert_mesh(&piped.mesh, piped.volume);
}

/// That the cap with `role` lies at `point`, facing along `outward`.
fn assert_cap(piped: &Piped, role: PipeFace, point: [f64; 3], outward: [f64; 3]) {
    let cap = piped.face(role);
    assert!(matches!(piped.surfaces[cap], Surface::Plane { .. }));
    // Integrated to 1e-7 of the face's size.
    assert_near("cap", piped.centroids[cap], point, 1e-5);
    assert_near("cap normal", piped.normals[cap].unwrap(), outward, 1e-9);
}

/// That the wall is a cylinder of `radius` round the line through `point`
/// along `along`.
fn assert_cylinder(surface: &Surface, point: [f64; 3], along: [f64; 3], radius: f64) {
    let Surface::Cylinder {
        origin,
        axis,
        radius: r,
        ..
    } = *surface
    else {
        panic!("{surface:?} isn't a cylinder");
    };
    assert!(1.0 - dot(axis, along).abs() < 1e-12, "{axis:?}");
    let d = sub(origin, point);
    let off = sub(d, scale(along, dot(d, along)));
    assert!(dot(off, off).sqrt() < 1e-9, "{origin:?}");
    assert!((r - radius).abs() < 1e-9, "radius {r}");
}

/// That the wall is a torus round `centre` and `axis`, of radii `major` and
/// `minor`.
fn assert_torus(surface: &Surface, centre: [f64; 3], axis: [f64; 3], major: f64, minor: f64) {
    let Surface::Torus {
        center,
        axis: along,
        major_radius,
        minor_radius,
        ..
    } = *surface
    else {
        panic!("{surface:?} isn't a torus");
    };
    assert!(1.0 - dot(along, axis).abs() < 1e-12, "{along:?}");
    assert_near("centre", center, centre, 1e-9);
    assert!((major_radius - major).abs() < 1e-9 && (minor_radius - minor).abs() < 1e-9);
}

#[test]
fn a_straight_pipe_is_a_cylinder() {
    for frame in [flat(), tilted()] {
        for inner in [None, Some(INNER)] {
            let piped = pipe(frame, vec![line([0.0, 0.0], [50.0, 0.0])], inner);
            assert_piped(&piped, 1, inner, false, section(inner) * 50.0);
            let x = frame.x_axis;
            assert_cap(&piped, PipeFace::Start, frame.origin, scale(x, -1.0));
            assert_cap(&piped, PipeFace::End, frame.to_world([50.0, 0.0, 0.0]), x);
            let outer = piped.face(PipeFace::Outer { segment: Some(0) });
            assert_cylinder(&piped.surfaces[outer], frame.origin, x, R);
            if inner.is_some() {
                let wall = piped.face(PipeFace::Inner { segment: Some(0) });
                assert_cylinder(&piped.surfaces[wall], frame.origin, x, INNER);
            }
        }
    }
}

#[test]
fn a_bend_is_a_torus() {
    // Along x, a quarter turn either way round a 10 radius, then 30 on.
    for turn in [1.0, -1.0] {
        let path = vec![
            line([0.0, 0.0], [40.0, 0.0]),
            PlaneCurve::Arc {
                centre: [40.0, 10.0 * turn],
                radius: 10.0,
                start_deg: -90.0 * turn,
                sweep_deg: 90.0 * turn,
            },
            line([50.0, 10.0 * turn], [50.0, 40.0 * turn]),
        ];
        for frame in [flat(), tilted()] {
            for inner in [None, Some(INNER)] {
                let piped = pipe(frame, path.clone(), inner);
                let length = 40.0 + 10.0 * PI / 2.0 + 30.0;
                assert_piped(&piped, 3, inner, false, section(inner) * length);
                let centre = frame.to_world([40.0, 10.0 * turn, 0.0]);
                let z = frame.z_axis();
                let bend = piped.face(PipeFace::Outer { segment: Some(1) });
                assert_torus(&piped.surfaces[bend], centre, z, 10.0, R);
                if inner.is_some() {
                    let wall = piped.face(PipeFace::Inner { segment: Some(1) });
                    assert_torus(&piped.surfaces[wall], centre, z, 10.0, INNER);
                }
                let y = scale(frame.y_axis, turn);
                let last = piped.face(PipeFace::Outer { segment: Some(2) });
                assert_cylinder(
                    &piped.surfaces[last],
                    frame.to_world([50.0, 0.0, 0.0]),
                    y,
                    R,
                );
                assert_cap(
                    &piped,
                    PipeFace::End,
                    frame.to_world([50.0, 40.0 * turn, 0.0]),
                    y,
                );
            }
        }
    }
}

#[test]
fn corners_are_mitred() {
    // Two lines turning by up to 165°: each wall whole up to the plane
    // halfway between them, and the pipe as long as its path.
    let frame = flat();
    for turn in [30.0f64, 90.0, 150.0, 165.0] {
        let (sin, cos) = turn.to_radians().sin_cos();
        let path = vec![
            line([0.0, 0.0], [50.0, 0.0]),
            line([50.0, 0.0], [50.0 + 40.0 * cos, 40.0 * sin]),
        ];
        for inner in [None, Some(INNER)] {
            let piped = pipe(frame, path.clone(), inner);
            assert_piped(&piped, 2, inner, false, section(inner) * 90.0);
            let wall = piped.face(PipeFace::Outer { segment: Some(1) });
            assert_cylinder(&piped.surfaces[wall], [50.0, 0.0, 0.0], [cos, sin, 0.0], R);
        }
    }
}

#[test]
fn closed_paths_have_no_caps() {
    let frame = tilted();
    // A ring.
    let ring = vec![PlaneCurve::Circle {
        centre: [0.0, 0.0],
        radius: 30.0,
    }];
    let piped = pipe(frame, ring, None);
    assert_piped(&piped, 1, None, true, section(None) * 2.0 * PI * 30.0);
    let wall = piped.face(PipeFace::Outer { segment: Some(0) });
    assert_torus(&piped.surfaces[wall], frame.origin, frame.z_axis(), 30.0, R);

    // A slot: tangent all round.
    let slot = vec![
        line([0.0, 0.0], [40.0, 0.0]),
        PlaneCurve::Arc {
            centre: [40.0, 15.0],
            radius: 15.0,
            start_deg: -90.0,
            sweep_deg: 180.0,
        },
        line([40.0, 30.0], [0.0, 30.0]),
        PlaneCurve::Arc {
            centre: [0.0, 15.0],
            radius: 15.0,
            start_deg: 90.0,
            sweep_deg: 180.0,
        },
    ];
    let length = 80.0 + 2.0 * PI * 15.0;
    assert_piped(
        &pipe(frame, slot, None),
        4,
        None,
        true,
        section(None) * length,
    );

    // A square: mitred all round.
    let corners = [[0.0, 0.0], [40.0, 0.0], [40.0, 40.0], [0.0, 40.0]];
    let square = (0..4)
        .map(|i| line(corners[i], corners[(i + 1) % 4]))
        .collect();
    assert_piped(
        &pipe(frame, square, None),
        4,
        None,
        true,
        section(None) * 160.0,
    );
}

#[test]
fn a_pipe_joins_a_block() {
    // A handle: a bent pipe whose ends sink into a block's top.
    let (manifold, volume, pipe_volume) = on_kernel(|k| {
        let frame = Frame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
        let path = [
            line([0.0, -5.0], [0.0, 20.0]),
            PlaneCurve::Arc {
                centre: [15.0, 20.0],
                radius: 15.0,
                start_deg: 180.0,
                sweep_deg: -180.0,
            },
            line([30.0, 20.0], [30.0, -5.0]),
        ];
        let handle = Brep::pipe(k, &frame, &path, 4.0, None).unwrap().solid;
        let pipe_volume = handle.volume(1e-9).unwrap();
        let block = Brep::cuboid(k, [-20.0, -20.0, -30.0], [70.0, 40.0, 30.0]).unwrap();
        let part = block.union(handle).unwrap();
        (
            part.is_manifold().unwrap(),
            part.volume(1e-9).unwrap(),
            pipe_volume,
        )
    });
    assert!(manifold);
    // The block plus the pipe less its two 5 long ends sunk in the block.
    let sunk = 2.0 * PI * 16.0 * 5.0;
    let expected = 70.0 * 40.0 * 30.0 + pipe_volume - sunk;
    assert!(close(volume, expected), "{volume} vs {expected}");
}

#[test]
fn refuses_paths_it_cant_follow() {
    let reasons = on_kernel(|k| {
        let frame = flat();
        let bend = |radius| PlaneCurve::Arc {
            centre: [0.0, radius],
            radius,
            start_deg: -90.0,
            sweep_deg: 90.0,
        };
        let circle = PlaneCurve::Circle {
            centre: [0.0, 0.0],
            radius: 20.0,
        };
        let results = [
            Brep::pipe(k, &frame, &[], R, None),
            Brep::pipe(
                k,
                &frame,
                &[
                    line([0.0, 0.0], [10.0, 0.0]),
                    line([11.0, 0.0], [20.0, 0.0]),
                ],
                R,
                None,
            ),
            Brep::pipe(
                k,
                &frame,
                &[circle, line([20.0, 0.0], [30.0, 0.0])],
                R,
                None,
            ),
            Brep::pipe(k, &frame, &[line([0.0, 0.0], [10.0, 0.0])], 0.0, None),
            Brep::pipe(k, &frame, &[line([0.0, 0.0], [10.0, 0.0])], R, Some(R)),
            Brep::pipe(k, &frame, &[bend(R)], R, None),
            // Out and straight back.
            Brep::pipe(
                k,
                &frame,
                &[line([0.0, 0.0], [10.0, 0.0]), line([10.0, 0.0], [0.0, 0.0])],
                R,
                None,
            ),
            // A hollow ring.
            Brep::pipe(k, &frame, &[circle], R, Some(INNER)),
            // Along x, then a bend heading down y.
            Brep::pipe(
                k,
                &frame,
                &[
                    line([0.0, 0.0], [40.0, 0.0]),
                    PlaneCurve::Arc {
                        centre: [70.0, 0.0],
                        radius: 30.0,
                        start_deg: 180.0,
                        sweep_deg: 60.0,
                    },
                ],
                R,
                None,
            ),
        ];
        results
            .into_iter()
            .map(|r| match r {
                Err(Error::InvalidArgument { reason, .. }) => reason,
                other => panic!("{other:?}"),
            })
            .collect::<Vec<_>>()
    });
    assert_eq!(
        reasons,
        [
            "the path needs at least one curve",
            "each curve of the path must start where the last ended",
            "a circle is a path on its own",
            "the radius must be finite and positive",
            "the inner radius must be finite, not negative, and less than the radius",
            "a bend's radius must be more than the pipe's",
            "the path turns back on itself",
            "a pipe along a closed path can't be hollow: its bore would be sealed inside",
            "an arc in the path must meet the curves either side of it tangentially",
        ]
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

fn assert_near(what: &str, actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    let gap = dot(sub(actual, expected), sub(actual, expected)).sqrt();
    assert!(
        gap <= tolerance,
        "{what}: {actual:?} isn't within {tolerance} of {expected:?}"
    );
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

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
