//! Planar faces revolved about an axis in their plane: what a sketch's
//! Revolve needs. Every side must trace back to the edge it was turned from,
//! the solid must be closed, and its sides must keep analytic surfaces.
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

/// A plane tilted 30° about x, away from the origin. Profiles are drawn in
/// its (x, y) and turned about its y axis.
fn tilted() -> Frame {
    let (sin, cos) = 30f64.to_radians().sin_cos();
    Frame::new([10.0, 20.0, 30.0], [1.0, 0.0, 0.0], [0.0, cos, sin]).unwrap()
}

fn edge(index: u32) -> Label {
    Label::new(SKETCH, "edge").with_index(index)
}

/// The closed polygon through `points`, its lines labelled from 0.
fn polygon(points: &[[f64; 2]]) -> Vec<BoundaryCurve> {
    (0..points.len())
        .map(|i| {
            let line = PlaneCurve::Line {
                start: points[i],
                end: points[(i + 1) % points.len()],
            };
            BoundaryCurve::labelled(line, edge(i as u32))
        })
        .collect()
}

fn rectangle(min: [f64; 2], max: [f64; 2]) -> Vec<BoundaryCurve> {
    polygon(&[min, [max[0], min[1]], max, [min[0], max[1]]])
}

/// What the tests look at in a sheet's face revolved about the frame's y
/// axis.
#[derive(Debug)]
struct Revolved {
    /// The sheet's edges' labels, by edge index.
    edge_labels: Vec<Option<Label>>,
    roles: Vec<SweptFace>,
    /// Each face's surface and label, in the solid's face order.
    surfaces: Vec<Surface>,
    labels: Vec<Option<Label>>,
    centroids: Vec<[f64; 3]>,
    manifold: bool,
    volume: f64,
    mesh: Mesh,
}

impl Revolved {
    /// The label of the sheet edge each side came from, by side.
    fn side_labels(&self) -> Vec<(usize, Label)> {
        self.roles
            .iter()
            .enumerate()
            .filter_map(|(face, role)| match role {
                SweptFace::Side { source_edge } => {
                    let edge = source_edge.expect("every side has a source edge");
                    let label = self.edge_labels[edge]
                        .clone()
                        .expect("every edge is labelled");
                    Some((face, label))
                }
                _ => None,
            })
            .collect()
    }

    /// The surface of the side turned from the edge labelled `index`.
    fn side_from(&self, index: u32) -> &Surface {
        let (face, _) = self
            .side_labels()
            .into_iter()
            .find(|(_, label)| label.index == index)
            .unwrap_or_else(|| panic!("no side from edge {index}"));
        &self.surfaces[face]
    }

    fn face_with(&self, role: SweptFace) -> Option<usize> {
        self.roles.iter().position(|r| *r == role)
    }
}

/// Builds the loops on the tilted frame and revolves the sheet's one face by
/// `angle_deg` about the frame's y axis.
fn revolve(loops: Vec<Vec<BoundaryCurve>>, angle_deg: f64) -> Revolved {
    let frame = tilted();
    on_kernel(move |k| {
        let sheet = Brep::planar_faces(k, &frame, &loops).unwrap();
        let edge_labels = sheet
            .edges()
            .unwrap()
            .into_iter()
            .map(|e| e.label)
            .collect();
        let swept = sheet
            .revolve_face(0, frame.origin, frame.y_axis, angle_deg)
            .unwrap();
        let solid = &swept.prism;
        let faces = solid.faces().unwrap();
        Revolved {
            edge_labels,
            surfaces: faces.iter().map(|f| f.surface.clone()).collect(),
            labels: faces.iter().map(|f| f.label.clone()).collect(),
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

/// What every revolve must be: closed, unlabelled, its sides traced to the
/// sheet's edges `expected` (by label index), with caps only for a partial
/// turn, the end being the start turned by the angle, and a mesh that winds
/// outwards round the solid's volume.
fn assert_revolved(revolved: &Revolved, angle_deg: f64, expected: &[u32], volume: f64) {
    assert!(revolved.manifold);
    assert!(
        close(revolved.volume, volume),
        "{} vs {volume}",
        revolved.volume
    );
    assert!(revolved.labels.iter().all(Option::is_none));
    let traced: Vec<u32> = revolved
        .side_labels()
        .iter()
        .map(|(_, l)| l.index)
        .collect();
    let unique: BTreeSet<u32> = traced.iter().copied().collect();
    assert_eq!(
        traced.len(),
        unique.len(),
        "a side traced twice: {traced:?}"
    );
    assert_eq!(unique, expected.iter().copied().collect());

    let whole = angle_deg.abs() == 360.0;
    let (start, end) = (
        revolved.face_with(SweptFace::Start),
        revolved.face_with(SweptFace::End),
    );
    assert_eq!(
        revolved.roles.len(),
        expected.len() + if whole { 0 } else { 2 }
    );
    if whole {
        assert_eq!((start, end), (None, None));
    } else {
        let (start, end) = (start.unwrap(), end.unwrap());
        for cap in [start, end] {
            assert!(matches!(revolved.surfaces[cap], Surface::Plane { .. }));
        }
        let frame = tilted();
        let turned = rotate(
            revolved.centroids[start],
            frame.origin,
            frame.y_axis,
            angle_deg,
        );
        assert_near("end", revolved.centroids[end], turned, 1e-6);
    }
    assert_mesh(&revolved.mesh, revolved.volume);
}

#[test]
fn a_rectangle_beside_the_axis_turns_into_a_tube() {
    for angle in [360.0, 90.0, -90.0, 180.0, 270.0, -360.0] {
        let revolved = revolve(vec![rectangle([5.0, 0.0], [10.0, 20.0])], angle);
        let volume = PI * (100.0 - 25.0) * 20.0 * angle.abs() / 360.0;
        assert_revolved(&revolved, angle, &[0, 1, 2, 3], volume);
        // The lines across the axis turn into planes square to it, the
        // lines along it into cylinders round it.
        let frame = tilted();
        for index in [0, 2] {
            let Surface::Plane { normal, .. } = *revolved.side_from(index) else {
                panic!("{:?}", revolved.side_from(index));
            };
            assert!(1.0 - dot(normal, frame.y_axis).abs() < 1e-12);
        }
        for (index, radius) in [(1, 10.0), (3, 5.0)] {
            assert_round(revolved.side_from(index), "cylinder", radius, 0.0);
        }
    }
}

#[test]
fn an_edge_on_the_axis_turns_into_nothing() {
    for angle in [360.0, 45.0] {
        let revolved = revolve(vec![rectangle([0.0, 0.0], [10.0, 20.0])], angle);
        // The edge on the axis, 3, has no side: the ends are discs.
        let volume = PI * 100.0 * 20.0 * angle / 360.0;
        assert_revolved(&revolved, angle, &[0, 1, 2], volume);
        assert_round(revolved.side_from(1), "cylinder", 10.0, 0.0);
    }
}

#[test]
fn slanted_edges_turn_into_cones() {
    // A wall leaning in towards the top, and a triangle whose slant meets
    // the axis: a cone up to its apex.
    let wall = polygon(&[[5.0, 0.0], [10.0, 0.0], [8.0, 20.0], [5.0, 20.0]]);
    let triangle = polygon(&[[0.0, 0.0], [10.0, 0.0], [0.0, 20.0]]);
    for angle in [360.0, 120.0] {
        let wall_volume = f64::to_radians(angle) * frustum_moment(5.0, 10.0, 8.0, 20.0);
        let revolved = revolve(vec![wall.clone()], angle);
        assert_revolved(&revolved, angle, &[0, 1, 2, 3], wall_volume);
        assert_round(
            revolved.side_from(1),
            "cone",
            f64::NAN,
            (2.0f64 / 20.0).atan(),
        );
        assert_round(revolved.side_from(3), "cylinder", 5.0, 0.0);

        let revolved = revolve(vec![triangle.clone()], angle);
        let cone_volume = PI * 100.0 * 20.0 / 3.0 * angle / 360.0;
        assert_revolved(&revolved, angle, &[0, 1], cone_volume);
        assert_round(
            revolved.side_from(1),
            "cone",
            f64::NAN,
            (10.0f64 / 20.0).atan(),
        );
    }
}

#[test]
fn arcs_turn_into_spheres_and_tori() {
    // A half disc on the axis turns into a ball.
    let half = vec![
        BoundaryCurve::labelled(
            PlaneCurve::Line {
                start: [0.0, -10.0],
                end: [0.0, 10.0],
            },
            edge(0),
        ),
        BoundaryCurve::labelled(
            PlaneCurve::Arc {
                centre: [0.0, 0.0],
                radius: 10.0,
                start_deg: 90.0,
                sweep_deg: -180.0,
            },
            edge(1),
        ),
    ];
    // A circle off the axis turns into a ring, and a round hole in a
    // rectangle into a ring-shaped hole in a tube.
    let ring = vec![vec![BoundaryCurve::labelled(
        PlaneCurve::Circle {
            centre: [20.0, 0.0],
            radius: 5.0,
        },
        edge(0),
    )]];
    let holed = vec![
        rectangle([5.0, 0.0], [10.0, 20.0]),
        vec![BoundaryCurve::labelled(
            PlaneCurve::Circle {
                centre: [7.5, 10.0],
                radius: 1.0,
            },
            edge(4),
        )],
    ];
    for angle in [360.0, 90.0] {
        let share = angle / 360.0;
        let revolved = revolve(vec![half.clone()], angle);
        assert_revolved(&revolved, angle, &[1], 4.0 / 3.0 * PI * 1000.0 * share);
        let Surface::Sphere { radius, center, .. } = *revolved.side_from(1) else {
            panic!("{:?}", revolved.side_from(1));
        };
        assert!((radius - 10.0).abs() < 1e-9);
        assert_near("centre", center, tilted().origin, 1e-9);

        let revolved = revolve(ring.clone(), angle);
        let ring_volume = 2.0 * PI * PI * 20.0 * 25.0 * share;
        assert_revolved(&revolved, angle, &[0], ring_volume);
        assert_torus(revolved.side_from(0), 20.0, 5.0);

        let revolved = revolve(holed.clone(), angle);
        let tube = PI * (100.0 - 25.0) * 20.0 * share;
        let hole = 2.0 * PI * PI * 7.5 * share;
        assert_revolved(&revolved, angle, &[0, 1, 2, 3, 4], tube - hole);
        assert_torus(revolved.side_from(4), 7.5, 1.0);
    }
}

#[test]
fn a_face_touching_the_axis_at_a_point_turns_into_one_solid() {
    // Two cones meeting at their apex on the axis, and a cylinder round
    // them.
    let triangle = polygon(&[[0.0, 0.0], [10.0, 5.0], [10.0, -5.0]]);
    for angle in [360.0, 90.0] {
        let revolved = revolve(vec![triangle.clone()], angle);
        let volume = 50.0 * (20.0 / 3.0) * angle.to_radians();
        assert_revolved(&revolved, angle, &[0, 1, 2], volume);
        assert_round(revolved.side_from(1), "cylinder", 10.0, 0.0);
        for index in [0, 2] {
            assert_round(revolved.side_from(index), "cone", f64::NAN, (2.0f64).atan());
        }
    }
}

#[test]
fn a_solids_face_turns_about_its_own_edge() {
    let (roles, volume, sides_from) = on_kernel(|k| {
        let block = Brep::cuboid(k, [0.0, 0.0, 0.0], [10.0, 20.0, 30.0]).unwrap();
        let faces = block.faces().unwrap();
        let face = faces
            .iter()
            .position(|f| f.normal == Some([1.0, 0.0, 0.0]))
            .unwrap();
        let edges = block.edges().unwrap();
        // About the face's edge along z at y = 0, clockwise looking down z:
        // out of the block.
        let swept = block
            .revolve_face(face, [10.0, 0.0, 0.0], [0.0, 0.0, 1.0], -90.0)
            .unwrap();
        assert!(swept.prism.is_manifold().unwrap());
        // Joined back on, it adds its own volume.
        let joined = block
            .try_clone()
            .unwrap()
            .union(swept.prism.try_clone().unwrap())
            .unwrap();
        assert!(joined.is_manifold().unwrap());
        let added = joined.volume(1e-9).unwrap() - 6000.0;
        assert!(close(added, swept.prism.volume(1e-9).unwrap()), "{added}");
        let sides_from: Vec<[f64; 3]> = swept
            .faces
            .iter()
            .filter_map(|role| match role {
                SweptFace::Side { source_edge } => Some(edges[source_edge.unwrap()].samples[2]),
                _ => None,
            })
            .collect();
        (swept.faces, swept.prism.volume(1e-9).unwrap(), sides_from)
    });
    assert_eq!(roles.len(), 5);
    // Pappus: a 20 × 30 face whose centroid is 10 from the axis.
    assert!(close(volume, 600.0 * 10.0 * PI / 2.0), "{volume}");
    // The sides come from the three edges off the axis: the top, the
    // bottom and the far side.
    let mut midpoints: Vec<[i64; 3]> = sides_from
        .iter()
        .map(|p| p.map(|c| c.round() as i64))
        .collect();
    midpoints.sort();
    assert_eq!(midpoints, vec![[10, 10, 0], [10, 10, 30], [10, 20, 15]]);
}

#[test]
fn refuses_axes_off_the_plane_or_through_the_face() {
    let reasons = on_kernel(|k| {
        let frame = tilted();
        let sheet = Brep::planar_faces(k, &frame, &[rectangle([-5.0, 0.0], [10.0, 20.0])]).unwrap();
        let ball = Brep::sphere(k, [0.0; 3], 5.0).unwrap();
        let z = frame.z_axis();
        let off = add(frame.origin, z);
        let results = [
            // Through the face's inside.
            sheet.revolve_face(0, frame.origin, frame.y_axis, 90.0),
            // Out of its plane: leaning out of it, or beside it.
            sheet.revolve_face(0, frame.origin, add(frame.y_axis, z), 90.0),
            sheet.revolve_face(0, off, frame.y_axis, 90.0),
            // No direction, no angle, more than a turn, not finite.
            sheet.revolve_face(0, frame.origin, [0.0; 3], 90.0),
            sheet.revolve_face(0, frame.origin, frame.y_axis, 0.0),
            sheet.revolve_face(0, frame.origin, frame.y_axis, 360.5),
            sheet.revolve_face(0, [f64::NAN, 0.0, 0.0], frame.y_axis, 90.0),
            // Not a planar face.
            ball.revolve_face(0, [20.0, 0.0, 0.0], [0.0, 0.0, 1.0], 90.0),
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
            "the axis runs through the face, which can't revolve across it",
            "the axis must lie in the face's plane",
            "the axis must lie in the face's plane",
            "the axis direction must be finite and not zero",
            "the angle must be finite, not zero, and at most a whole turn",
            "the angle must be finite, not zero, and at most a whole turn",
            "coordinates must be finite",
            "only a planar face can be revolved",
        ]
    );
}

/// The first moment about the axis (x = 0) of the quadrilateral (r0, 0),
/// (r1, 0), (r2, h), (r0, h), which has a side along the axis' direction at
/// r0: times the angle turned, the volume it sweeps (Pappus).
fn frustum_moment(r0: f64, r1: f64, r2: f64, h: f64) -> f64 {
    // ∫∫ x dA over 0 ≤ y ≤ h, r0 ≤ x ≤ r(y), with r running r1 to r2.
    h / 2.0 * ((r1 * r1 + r1 * r2 + r2 * r2) / 3.0 - r0 * r0)
}

/// That `surface` is a cylinder or cone round the frame's y axis through its
/// origin: of `radius` (unchecked when NaN) and half-angle `half_angle`
/// (radians, either sign).
fn assert_round(surface: &Surface, kind: &str, radius: f64, half_angle: f64) {
    let frame = tilted();
    let (origin, axis, r, angle) = match (kind, surface) {
        (
            "cylinder",
            Surface::Cylinder {
                origin,
                axis,
                radius,
                ..
            },
        ) => (*origin, *axis, *radius, 0.0),
        (
            "cone",
            Surface::Cone {
                origin,
                axis,
                radius,
                half_angle,
                ..
            },
        ) => (*origin, *axis, *radius, *half_angle),
        _ => panic!("{surface:?} isn't a {kind}"),
    };
    assert!(1.0 - dot(axis, frame.y_axis).abs() < 1e-12, "{axis:?}");
    let d = sub(origin, frame.origin);
    let off_axis = sub(d, scale(frame.y_axis, dot(d, frame.y_axis)));
    assert!(dot(off_axis, off_axis).sqrt() < 1e-9, "{origin:?}");
    if !radius.is_nan() {
        assert!((r - radius).abs() < 1e-9, "radius {r}");
    }
    assert!(
        (angle.abs() - half_angle).abs() < 1e-9,
        "half-angle {angle}"
    );
}

fn assert_torus(surface: &Surface, major: f64, minor: f64) {
    let frame = tilted();
    let Surface::Torus {
        center,
        axis,
        major_radius,
        minor_radius,
        ..
    } = *surface
    else {
        panic!("{surface:?} isn't a torus");
    };
    assert!(1.0 - dot(axis, frame.y_axis).abs() < 1e-12, "{axis:?}");
    let d = sub(center, frame.origin);
    let off_axis = sub(d, scale(frame.y_axis, dot(d, frame.y_axis)));
    assert!(dot(off_axis, off_axis).sqrt() < 1e-9, "{center:?}");
    assert!((major_radius - major).abs() < 1e-9 && (minor_radius - minor).abs() < 1e-9);
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

/// `point` turned by `angle_deg` about the axis through `origin` along the
/// unit `axis`, anticlockwise looking down it from its tip.
fn rotate(point: [f64; 3], origin: [f64; 3], axis: [f64; 3], angle_deg: f64) -> [f64; 3] {
    let (sin, cos) = angle_deg.to_radians().sin_cos();
    let p = sub(point, origin);
    let along = scale(axis, dot(p, axis));
    let across = sub(p, along);
    let turned = add(scale(across, cos), scale(cross(axis, across), sin));
    add(origin, add(along, turned))
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
