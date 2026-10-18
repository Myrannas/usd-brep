//! Planar faces built from lines, arcs and circles in a plane, and the prisms
//! swept from them: what a sketch's extrude needs. Every side of a prism must
//! trace back to the curve it was swept from, through the curve's label.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{BoundaryCurve, Brep, Error, Frame, Kernel, Label, PlaneCurve, Surface, SweptFace};

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

const SKETCH: u32 = 7;

/// A plane tilted 30° about x, away from the origin.
fn tilted() -> Frame {
    let (sin, cos) = 30f64.to_radians().sin_cos();
    Frame::new([10.0, 20.0, 30.0], [1.0, 0.0, 0.0], [0.0, cos, sin]).unwrap()
}

fn side(index: u32) -> Label {
    Label::new(SKETCH, "side").with_index(index)
}

/// The rectangle from `min` to `max` as four lines labelled `first`
/// onwards, counter-clockwise or clockwise.
fn rectangle(min: [f64; 2], max: [f64; 2], first: u32, clockwise: bool) -> Vec<BoundaryCurve> {
    let mut corners = [min, [max[0], min[1]], max, [min[0], max[1]]];
    if clockwise {
        corners.reverse();
    }
    (0..4)
        .map(|i| {
            let line = PlaneCurve::Line {
                start: corners[i],
                end: corners[(i + 1) % 4],
            };
            BoundaryCurve::labelled(line, side(first + i as u32))
        })
        .collect()
}

fn circle(centre: [f64; 2], radius: f64, index: u32) -> Vec<BoundaryCurve> {
    vec![BoundaryCurve::labelled(
        PlaneCurve::Circle { centre, radius },
        side(index),
    )]
}

/// What the tests look at in a sheet and a prism swept from its face.
#[derive(Debug)]
struct Swept {
    faces: usize,
    normal: [f64; 3],
    area: f64,
    edge_labels: Vec<Option<Label>>,
    edges_in_plane: bool,
    roles: Vec<SweptFace>,
    /// Each side's surface, in the prism's face order.
    sides: Vec<Surface>,
    start: [f64; 3],
    end: [f64; 3],
    volume: f64,
}

/// Builds the loops on `frame` and sweeps the sheet's one face by
/// `distance`.
fn sweep(loops: Vec<Vec<BoundaryCurve>>, frame: Frame, distance: f64) -> Swept {
    on_kernel(move |k| {
        let sheet = Brep::planar_faces(k, &frame, &loops).unwrap();
        let faces = sheet.faces().unwrap();
        let edges = sheet.edges().unwrap();
        let z = frame.z_axis();
        let edges_in_plane = edges.iter().all(|e| {
            e.samples
                .iter()
                .all(|&p| dot(sub(p, frame.origin), z).abs() < 1e-9)
        });
        let area = sheet.face_centroid(0).unwrap().area;
        let swept = sheet.sweep_face(0, distance).unwrap();
        let prism_faces = swept.prism.faces().unwrap();
        let centroid = |role| {
            let index = swept.faces.iter().position(|r| *r == role).unwrap();
            swept.prism.face_centroid(index).unwrap().centroid
        };
        let sides = swept
            .faces
            .iter()
            .zip(&prism_faces)
            .filter(|(role, _)| matches!(role, SweptFace::Side { .. }))
            .map(|(_, face)| {
                // A plane has a normal; nothing else does.
                assert_eq!(
                    face.normal.is_some(),
                    matches!(face.surface, Surface::Plane { .. })
                );
                face.surface.clone()
            })
            .collect();
        assert!(swept.prism.is_manifold().unwrap());
        Swept {
            faces: faces.len(),
            normal: faces[0].normal.unwrap(),
            area,
            edge_labels: edges.iter().map(|e| e.label.clone()).collect(),
            edges_in_plane,
            sides,
            start: centroid(SweptFace::Start),
            end: centroid(SweptFace::End),
            volume: swept.prism.volume(1e-9).unwrap(),
            roles: swept.faces,
        }
    })
}

/// Each side's source edge's label, which must be there.
fn side_labels(swept: &Swept) -> Vec<Label> {
    swept
        .roles
        .iter()
        .filter_map(|role| match role {
            SweptFace::Side { source_edge } => {
                let edge = source_edge.expect("every side has a source edge");
                Some(
                    swept.edge_labels[edge]
                        .clone()
                        .expect("every edge is labelled"),
                )
            }
            _ => None,
        })
        .collect()
}

#[test]
fn a_rectangle_with_a_hole_faces_its_frame_and_traces_its_sides() {
    let frame = tilted();
    for clockwise in [false, true] {
        for distance in [12.0, -12.0] {
            let loops = vec![
                rectangle([-20.0, -10.0], [20.0, 10.0], 0, clockwise),
                circle([5.0, 0.0], 5.0, 4),
            ];
            let swept = sweep(loops, frame, distance);
            assert_eq!(swept.faces, 1);
            assert_near("normal", swept.normal, frame.z_axis(), 1e-9);
            let area = 800.0 - 25.0 * PI;
            assert!(close(swept.area, area), "{}", swept.area);
            assert!(swept.edges_in_plane);
            assert!(swept.edge_labels.iter().all(Option::is_some));
            assert_eq!(swept.edge_labels.len(), 5);

            // Two caps and a side for each curve, each traced to its own.
            assert_eq!(swept.roles.len(), 7);
            let traced: BTreeSet<u32> = side_labels(&swept).iter().map(|l| l.index).collect();
            assert_eq!(traced, (0..5).collect());
            // The lines sweep into planes and the circle into a cylinder
            // round the hole, along the frame's z.
            let planes = swept
                .sides
                .iter()
                .filter(|s| matches!(s, Surface::Plane { .. }));
            assert_eq!(planes.count(), 4, "{:?}", swept.sides);
            let centre = frame.to_world([5.0, 0.0, 0.0]);
            assert_cylinders(&swept.sides, &[(centre, 5.0)], frame.z_axis());

            // The end is `distance` along the frame's z from the start.
            let along = scale(frame.z_axis(), distance);
            assert_near("end", swept.end, add(swept.start, along), 1e-6);
            assert!(close(swept.volume, area * 12.0), "{}", swept.volume);
        }
    }
}

#[test]
fn arcs_run_either_way() {
    // A slot 30 between centres and 10 wide, round its ends by arcs given
    // counter-clockwise, then the same loop the other way round.
    let line = |start, end| PlaneCurve::Line { start, end };
    let arc = |centre, start_deg, sweep_deg| PlaneCurve::Arc {
        centre,
        radius: 5.0,
        start_deg,
        sweep_deg,
    };
    let forward = [
        line([0.0, -5.0], [30.0, -5.0]),
        arc([30.0, 0.0], -90.0, 180.0),
        line([30.0, 5.0], [0.0, 5.0]),
        arc([0.0, 0.0], 90.0, 180.0),
    ];
    let backward = [
        arc([0.0, 0.0], 270.0, -180.0),
        line([0.0, 5.0], [30.0, 5.0]),
        arc([30.0, 0.0], 90.0, -180.0),
        line([30.0, -5.0], [0.0, -5.0]),
    ];
    let area = 300.0 + 25.0 * PI;
    for curves in [forward, backward] {
        let curves = curves
            .iter()
            .enumerate()
            .map(|(i, &c)| BoundaryCurve::labelled(c, side(i as u32)))
            .collect();
        let swept = sweep(
            vec![curves],
            Frame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap(),
            4.0,
        );
        assert_near("normal", swept.normal, [0.0, 0.0, 1.0], 1e-9);
        assert!(close(swept.area, area), "{}", swept.area);
        let traced: BTreeSet<u32> = side_labels(&swept).iter().map(|l| l.index).collect();
        assert_eq!(traced, (0..4).collect());
        // Two planes, and the ends round their centres.
        let planes = swept
            .sides
            .iter()
            .filter(|s| matches!(s, Surface::Plane { .. }));
        assert_eq!(planes.count(), 2, "{:?}", swept.sides);
        let ends = [([0.0; 3], 5.0), ([30.0, 0.0, 0.0], 5.0)];
        assert_cylinders(&swept.sides, &ends, [0.0, 0.0, 1.0]);
        assert!(close(swept.volume, area * 4.0), "{}", swept.volume);
    }
}

#[test]
fn loops_nest_into_faces_holes_and_islands() {
    let (faces, areas, ends) = on_kernel(|k| {
        let loops = vec![
            circle([0.0, 0.0], 10.0, 8),
            rectangle([-50.0, -50.0], [50.0, 50.0], 0, false),
            rectangle([-30.0, -30.0], [30.0, 30.0], 4, true),
        ];
        let frame = Frame::new([0.0, 0.0, 5.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
        let sheet = Brep::planar_faces(k, &frame, &loops).unwrap();
        let faces = sheet.faces().unwrap();
        let mut areas: Vec<f64> = (0..faces.len())
            .map(|f| sheet.face_centroid(f).unwrap().area)
            .collect();
        areas.sort_by(f64::total_cmp);
        // Each face sweeps along the frame's z: from x = 0 to 3.
        let ends: Vec<[f64; 3]> = (0..faces.len())
            .map(|f| {
                let swept = sheet.sweep_face(f, 3.0).unwrap();
                let end = swept.faces.iter().position(|r| *r == SweptFace::End);
                swept.prism.face_centroid(end.unwrap()).unwrap().centroid
            })
            .collect();
        (faces, areas, ends)
    });
    assert_eq!(faces.len(), 2, "the frame and the island");
    for face in &faces {
        assert_near("normal", face.normal.unwrap(), [1.0, 0.0, 0.0], 1e-9);
    }
    assert!(close(areas[0], 100.0 * PI), "{areas:?}");
    assert!(close(areas[1], 6400.0), "{areas:?}");
    for end in ends {
        assert!((end[0] - 3.0).abs() < 1e-9, "{end:?}");
    }
}

#[test]
fn bad_loops_are_refused() {
    let frame = Frame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    let results = on_kernel(move |k| {
        let mut open = rectangle([0.0, 0.0], [10.0, 10.0], 0, false);
        open.pop();
        let mut with_circle = rectangle([0.0, 0.0], [10.0, 10.0], 0, false);
        with_circle.extend(circle([5.0, 5.0], 1.0, 4));
        let flat = vec![BoundaryCurve::labelled(
            PlaneCurve::Circle {
                centre: [0.0, 0.0],
                radius: 0.0,
            },
            side(0),
        )];
        let bad_label = vec![BoundaryCurve::labelled(
            PlaneCurve::Circle {
                centre: [0.0, 0.0],
                radius: 1.0,
            },
            Label::new(SKETCH, ""),
        )];
        [
            vec![open],
            vec![with_circle],
            vec![flat],
            vec![bad_label],
            Vec::new(),
        ]
        .map(|loops| Brep::planar_faces(k, &frame, &loops).map(|sheet| sheet.face_count()))
    });
    for result in results {
        assert!(
            matches!(result, Err(Error::InvalidArgument { .. })),
            "{result:?}"
        );
    }
}

/// That the sides other than planes are cylinders along `axis`, one round
/// each (centre, radius), in any order.
fn assert_cylinders(sides: &[Surface], expected: &[([f64; 3], f64)], axis: [f64; 3]) {
    let mut left = expected.to_vec();
    for side in sides {
        if matches!(side, Surface::Plane { .. }) {
            continue;
        }
        let Surface::Cylinder {
            origin,
            axis: along,
            radius,
            ..
        } = *side
        else {
            panic!("{side:?} isn't a plane or a cylinder");
        };
        assert!(1.0 - dot(along, axis).abs() < 1e-12, "{along:?}");
        // The centre's distance from the axis.
        let offset = |centre| {
            let d = sub(centre, origin);
            let d = sub(d, scale(along, dot(d, along)));
            dot(d, d).sqrt()
        };
        let found = left
            .iter()
            .position(|&(centre, r)| offset(centre) < 1e-9 && (radius - r).abs() < 1e-9)
            .unwrap_or_else(|| panic!("{side:?} isn't one of {left:?}"));
        left.remove(found);
    }
    assert!(left.is_empty(), "no cylinders for {left:?}");
}

/// Within the accuracy areas and volumes are integrated to.
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
