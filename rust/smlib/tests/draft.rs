//! Planar faces swept with draft: what a sketch's Extrude with a draft angle
//! needs. Every side must trace back to the edge it was swept from, the
//! solid must be closed, and its sides must be planes and cones.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    BoundaryCurve, Brep, Error, Frame, Kernel, Label, Mesh, PlaneCurve, Surface, SweptFace,
    TessellationParams,
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

const SKETCH: u32 = 7;

fn flat() -> Frame {
    Frame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap()
}

/// A plane tilted 30° about x, away from the origin.
fn tilted() -> Frame {
    let (sin, cos) = 30f64.to_radians().sin_cos();
    Frame::new([10.0, 20.0, 30.0], [1.0, 0.0, 0.0], [0.0, cos, sin]).unwrap()
}

fn edge(index: u32) -> Label {
    Label::new(SKETCH, "edge").with_index(index)
}

fn line(start: [f64; 2], end: [f64; 2], index: u32) -> BoundaryCurve {
    BoundaryCurve::labelled(PlaneCurve::Line { start, end }, edge(index))
}

fn arc(centre: [f64; 2], radius: f64, start_deg: f64, sweep_deg: f64, index: u32) -> BoundaryCurve {
    let arc = PlaneCurve::Arc {
        centre,
        radius,
        start_deg,
        sweep_deg,
    };
    BoundaryCurve::labelled(arc, edge(index))
}

fn circle(centre: [f64; 2], radius: f64, index: u32) -> Vec<BoundaryCurve> {
    vec![BoundaryCurve::labelled(
        PlaneCurve::Circle { centre, radius },
        edge(index),
    )]
}

/// The square of side `2 * half` round the origin, its lines labelled 0 to 3.
fn square(half: f64) -> Vec<BoundaryCurve> {
    let corners = [[-half, -half], [half, -half], [half, half], [-half, half]];
    (0..4)
        .map(|i| line(corners[i], corners[(i + 1) % 4], i as u32))
        .collect()
}

/// A 30 × 30 square round the origin with its corners rounded to 5: lines
/// 0, 2, 4, 6 and quarter arcs 1, 3, 5, 7.
fn rounded() -> Vec<BoundaryCurve> {
    vec![
        line([-10.0, -15.0], [10.0, -15.0], 0),
        arc([10.0, -10.0], 5.0, -90.0, 90.0, 1),
        line([15.0, -10.0], [15.0, 10.0], 2),
        arc([10.0, 10.0], 5.0, 0.0, 90.0, 3),
        line([10.0, 15.0], [-10.0, 15.0], 4),
        arc([-10.0, 10.0], 5.0, 90.0, 90.0, 5),
        line([-15.0, 10.0], [-15.0, -10.0], 6),
        arc([-10.0, -10.0], 5.0, 180.0, 90.0, 7),
    ]
}

/// The rounded square's area, shrunk all round by `s` (grown when negative).
fn rounded_area(s: f64) -> f64 {
    (30.0 - 2.0 * s).powi(2) - (4.0 - PI) * (5.0 - s).powi(2)
}

/// What the tests look at in a sheet's face swept with draft.
#[derive(Debug)]
struct Drafted {
    /// The sheet's edges' labels, by edge index.
    edge_labels: Vec<Option<Label>>,
    roles: Vec<SweptFace>,
    surfaces: Vec<Surface>,
    labels: Vec<Option<Label>>,
    /// Each face's outward normal, for planes.
    normals: Vec<Option<[f64; 3]>>,
    centroids: Vec<[f64; 3]>,
    manifold: bool,
    volume: f64,
    mesh: Mesh,
}

impl Drafted {
    /// The label index of the sheet edge each side came from, by face.
    fn side_labels(&self) -> Vec<(usize, u32)> {
        self.roles
            .iter()
            .enumerate()
            .filter_map(|(face, role)| match role {
                SweptFace::Side { source_edge } => {
                    let edge = source_edge.expect("every side has a source edge");
                    let label = self.edge_labels[edge]
                        .clone()
                        .expect("every edge is labelled");
                    Some((face, label.index))
                }
                _ => None,
            })
            .collect()
    }

    /// The face swept from the edge labelled `index`.
    fn side_from(&self, index: u32) -> usize {
        self.side_labels()
            .into_iter()
            .find(|&(_, label)| label == index)
            .unwrap_or_else(|| panic!("no side from edge {index}"))
            .0
    }

    fn face_with(&self, role: SweptFace) -> usize {
        self.roles.iter().position(|r| *r == role).unwrap()
    }
}

/// Builds the loops on `frame` and sweeps the sheet's one face `distance`
/// with `draft_deg` of draft.
fn draft(frame: Frame, loops: Vec<Vec<BoundaryCurve>>, distance: f64, draft_deg: f64) -> Drafted {
    on_kernel(move |k| {
        let sheet = Brep::planar_faces(k, &frame, &loops).unwrap();
        let edge_labels = sheet
            .edges()
            .unwrap()
            .into_iter()
            .map(|e| e.label)
            .collect();
        let swept = sheet.sweep_face_tapered(0, distance, draft_deg).unwrap();
        let solid = &swept.prism;
        let faces = solid.faces().unwrap();
        Drafted {
            edge_labels,
            surfaces: faces.iter().map(|f| f.surface.clone()).collect(),
            labels: faces.iter().map(|f| f.label.clone()).collect(),
            normals: faces.iter().map(|f| f.normal).collect(),
            centroids: (0..faces.len())
                .map(|i| solid.face_centroid(i).unwrap().centroid)
                .collect(),
            manifold: solid.is_manifold().unwrap(),
            volume: solid.volume(1e-9).unwrap(),
            mesh: solid.tessellate(&TessellationParams::default()).unwrap(),
            roles: swept.faces,
        }
    })
}

/// What every drafted sweep must be: closed, unlabelled, a start on the face
/// and an end `distance` along the frame's z, each a plane, a side traced to
/// each of the sheet's edges `expected` (by label index), and a mesh that
/// winds outwards round the solid's volume.
fn assert_drafted(drafted: &Drafted, frame: Frame, distance: f64, expected: &[u32], volume: f64) {
    assert!(drafted.manifold);
    assert!(
        close(drafted.volume, volume),
        "{} vs {volume}",
        drafted.volume
    );
    assert!(drafted.labels.iter().all(Option::is_none));
    let traced: Vec<u32> = drafted.side_labels().iter().map(|&(_, l)| l).collect();
    let unique: BTreeSet<u32> = traced.iter().copied().collect();
    assert_eq!(
        traced.len(),
        unique.len(),
        "a side traced twice: {traced:?}"
    );
    assert_eq!(unique, expected.iter().copied().collect());
    assert_eq!(drafted.roles.len(), expected.len() + 2);

    let z = frame.z_axis();
    let (start, end) = (
        drafted.face_with(SweptFace::Start),
        drafted.face_with(SweptFace::End),
    );
    let along = |face: usize| dot(sub(drafted.centroids[face], frame.origin), z);
    assert!(along(start).abs() < 1e-9);
    assert!((along(end) - distance).abs() < 1e-9);
    let away = scale(z, distance.signum());
    assert_near(
        "start",
        drafted.normals[start].unwrap(),
        scale(away, -1.0),
        1e-9,
    );
    assert_near("end", drafted.normals[end].unwrap(), away, 1e-9);
    assert_mesh(&drafted.mesh, drafted.volume);
}

/// That the face is a plane whose outward normal leans `draft_deg` from
/// `across` (its straight side's) towards the sweep's direction `away`.
fn assert_leaning_plane(
    drafted: &Drafted,
    face: usize,
    across: [f64; 3],
    away: [f64; 3],
    draft_deg: f64,
) {
    assert!(
        matches!(drafted.surfaces[face], Surface::Plane { .. }),
        "{:?}",
        drafted.surfaces[face]
    );
    let (sin, cos) = draft_deg.to_radians().sin_cos();
    let expected = add(scale(across, cos), scale(away, sin));
    assert_near("normal", drafted.normals[face].unwrap(), expected, 1e-9);
}

/// That the face is a cone along `axis` through `centre` whose radius is
/// `radius` on the face's plane and leans `draft_deg` (either sign).
fn assert_cone(surface: &Surface, centre: [f64; 3], axis: [f64; 3], radius: f64, draft_deg: f64) {
    let Surface::Cone {
        origin,
        axis: along,
        radius: r,
        half_angle,
        ..
    } = *surface
    else {
        panic!("{surface:?} isn't a cone");
    };
    assert!(1.0 - dot(along, axis).abs() < 1e-12, "{along:?}");
    let d = sub(origin, centre);
    let height = dot(d, along);
    let off_axis = sub(d, scale(along, height));
    assert!(dot(off_axis, off_axis).sqrt() < 1e-9, "{origin:?}");
    // The radius where the axis meets the face's plane.
    let at_face = r - height * half_angle.tan();
    assert!((at_face - radius).abs() < 1e-9, "radius {at_face}");
    assert!(
        (half_angle.abs() - draft_deg.abs().to_radians()).abs() < 1e-9,
        "half-angle {half_angle}"
    );
}

#[test]
fn a_square_leans_in_and_out() {
    let frame = flat();
    for (distance, draft_deg) in [(10.0, 10.0), (10.0, -10.0), (-10.0, 10.0), (25.0, 20.0)] {
        let drafted = draft(frame, vec![square(10.0)], distance, draft_deg);
        // A frustum: the far square is 20 less twice the lean on a side.
        let far = 20.0 - 2.0 * distance.abs() * draft_deg.to_radians().tan();
        let (a, b) = (400.0, far * far);
        let volume = distance.abs() / 3.0 * (a + b + (a * b).sqrt());
        assert_drafted(&drafted, frame, distance, &[0, 1, 2, 3], volume);
        let away = scale(frame.z_axis(), distance.signum());
        let across = [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ];
        for (index, across) in across.into_iter().enumerate() {
            let face = drafted.side_from(index as u32);
            assert_leaning_plane(&drafted, face, across, away, draft_deg);
        }
    }
}

#[test]
fn circles_and_rounded_corners_lean_into_cones() {
    for frame in [flat(), tilted()] {
        let z = frame.z_axis();
        for draft_deg in [10.0f64, -10.0] {
            let t = draft_deg.to_radians().tan();
            let drafted = draft(frame, vec![circle([0.0, 0.0], 10.0, 0)], 10.0, draft_deg);
            let far = 10.0 - 10.0 * t;
            let volume = PI * 10.0 / 3.0 * (100.0 + far * far + 10.0 * far);
            assert_drafted(&drafted, frame, 10.0, &[0], volume);
            assert_cone(
                &drafted.surfaces[drafted.side_from(0)],
                frame.origin,
                z,
                10.0,
                draft_deg,
            );

            // The area shrinks by the lean all round: integrate it up.
            let drafted = draft(frame, vec![rounded()], 10.0, draft_deg);
            let volume = integrate(|h| rounded_area(h * t), 10.0);
            assert_drafted(&drafted, frame, 10.0, &(0..8).collect::<Vec<_>>(), volume);
            for (index, centre) in [(1, [10.0, -10.0]), (3, [10.0, 10.0]), (5, [-10.0, 10.0])] {
                let surface = &drafted.surfaces[drafted.side_from(index)];
                let centre = frame.to_world([centre[0], centre[1], 0.0]);
                assert_cone(surface, centre, z, 5.0, draft_deg);
            }
            let away = frame.z_axis();
            let across = scale(frame.x_axis, -1.0);
            assert_leaning_plane(&drafted, drafted.side_from(6), across, away, draft_deg);
        }
    }
}

#[test]
fn holes_lean_the_other_way() {
    let frame = flat();
    for draft_deg in [10.0f64, -10.0] {
        let t = draft_deg.to_radians().tan();
        let loops = vec![square(10.0), circle([0.0, 0.0], 3.0, 4)];
        let drafted = draft(frame, loops, 10.0, draft_deg);
        let far = 20.0 - 20.0 * t;
        let block = 10.0 / 3.0 * (400.0 + far * far + 20.0 * far);
        // The hole grows as the face shrinks.
        let wide = 3.0 + 10.0 * t;
        let hole = PI * 10.0 / 3.0 * (9.0 + wide * wide + 3.0 * wide);
        assert_drafted(&drafted, frame, 10.0, &[0, 1, 2, 3, 4], block - hole);
        assert_cone(
            &drafted.surfaces[drafted.side_from(4)],
            frame.origin,
            frame.z_axis(),
            3.0,
            draft_deg,
        );
    }
}

#[test]
fn a_solids_face_with_a_hole_in_it_leans() {
    // The top of a block with a round hole through it, drafted out of it.
    let (roles, volume, sides_from, manifold) = on_kernel(|k| {
        let block = Brep::cuboid(k, [-10.0, -10.0, -10.0], [20.0, 20.0, 10.0]).unwrap();
        let pin = Brep::cylinder(k, [0.0, 0.0, -20.0], 3.0, 40.0).unwrap();
        let part = block.difference(pin).unwrap();
        let faces = part.faces().unwrap();
        let top = faces
            .iter()
            .position(|f| f.normal == Some([0.0, 0.0, 1.0]))
            .unwrap();
        let edges = part.edges().unwrap();
        let swept = part.sweep_face_tapered(top, 10.0, 10.0).unwrap();
        let sides_from: Vec<[f64; 3]> = swept
            .faces
            .iter()
            .filter_map(|role| match role {
                SweptFace::Side { source_edge } => Some(edges[source_edge.unwrap()].samples[2]),
                _ => None,
            })
            .collect();
        let manifold = swept.prism.is_manifold().unwrap();
        let volume = swept.prism.volume(1e-9).unwrap();
        // And joined back on, as an extrude does.
        let before = part.volume(1e-9).unwrap();
        let joined = part.union(swept.prism).unwrap();
        assert!(joined.is_manifold().unwrap());
        assert!(close(joined.volume(1e-9).unwrap(), before + volume));
        (swept.faces, volume, sides_from, manifold)
    });
    assert!(manifold);
    assert_eq!(roles.len(), 7);
    let t = 10f64.to_radians().tan();
    let far = 20.0 - 20.0 * t;
    let wide = 3.0 + 10.0 * t;
    let expected = 10.0 / 3.0 * (400.0 + far * far + 20.0 * far)
        - PI * 10.0 / 3.0 * (9.0 + wide * wide + 3.0 * wide);
    assert!(close(volume, expected), "{volume} vs {expected}");
    // Four sides from the top's edges and one from the hole's rim.
    let mut midpoints: Vec<[i64; 3]> = sides_from
        .iter()
        .map(|p| p.map(|c| c.round() as i64))
        .collect();
    midpoints.sort();
    assert_eq!(midpoints.len(), 5);
    assert!(midpoints.iter().all(|p| p[2] == 0));
}

#[test]
fn an_arc_meeting_a_line_at_an_angle_comes_out_ruled() {
    // A D: a line and a half circle meeting square to it. The kernel offsets
    // the arc and trims it at the corners, and rules its wall between the
    // arc and its trimmed offset: a NURBS wall near the drafted cone, whose
    // corners are straight rulings, so the solid comes out a little small.
    let d = vec![
        line([0.0, -10.0], [0.0, 10.0], 0),
        arc([0.0, 0.0], 10.0, 90.0, -180.0, 1),
    ];
    let frame = flat();
    let drafted = draft(frame, vec![d.clone()], 10.0, 10.0);
    let t = 10f64.to_radians().tan();
    let exact = integrate(
        |h| {
            let (s, r) = (h * t, 10.0 - h * t);
            r * r * (s / r).acos() - s * (r * r - s * s).sqrt()
        },
        10.0,
    );
    assert!(drafted.manifold);
    assert!((drafted.volume - exact) / exact > -0.004 && drafted.volume < exact);
    assert_eq!(drafted.side_labels().len(), 2);
    assert!(matches!(
        drafted.surfaces[drafted.side_from(0)],
        Surface::Plane { .. }
    ));
    assert!(matches!(
        drafted.surfaces[drafted.side_from(1)],
        Surface::Nurbs { .. }
    ));
    // Leaning out, the kernel's corners leave no solid.
    let failed = on_kernel(move |k| {
        let sheet = Brep::planar_faces(k, &frame, &[d]).unwrap();
        sheet.sweep_face_tapered(0, 10.0, -10.0).err()
    });
    assert!(matches!(failed, Some(Error::Kernel { .. })), "{failed:?}");
}

#[test]
fn refuses_drafts_that_close_up_or_lean_too_far() {
    let reasons = on_kernel(|k| {
        let frame = flat();
        let sheet = Brep::planar_faces(k, &frame, &[square(10.0)]).unwrap();
        let holed =
            Brep::planar_faces(k, &frame, &[square(10.0), circle([0.0, 0.0], 3.0, 4)]).unwrap();
        let ball = Brep::sphere(k, [0.0; 3], 5.0).unwrap();
        let results = [
            // The square closes to a point 17.3 up.
            sheet.sweep_face_tapered(0, 30.0, 30.0),
            // The hole leans out through the sides.
            holed.sweep_face_tapered(0, 10.0, 30.0),
            sheet.sweep_face_tapered(0, 10.0, 90.0),
            sheet.sweep_face_tapered(0, 10.0, f64::NAN),
            sheet.sweep_face_tapered(0, 0.0, 10.0),
            ball.sweep_face_tapered(0, 10.0, 10.0),
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
            "the draft closes the face up before the end of the sweep",
            "the draft closes the face up before the end of the sweep",
            "the draft angle must be finite and under 90 degrees either way",
            "the draft angle must be finite and under 90 degrees either way",
            "sweep distance must be finite and non-zero",
            "only a planar face can be swept",
        ]
    );
}

/// ∫ f over [0, length], by Simpson's rule: exact for the quadratics a
/// shrinking face's area makes.
fn integrate(f: impl Fn(f64) -> f64, length: f64) -> f64 {
    let steps = 2000;
    let h = length / steps as f64;
    (0..steps)
        .map(|i| {
            let a = i as f64 * h;
            h / 6.0 * (f(a) + 4.0 * f(a + h / 2.0) + f(a + h))
        })
        .sum()
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
