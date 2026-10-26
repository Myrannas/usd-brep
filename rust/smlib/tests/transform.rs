//! Rigid transforms and scales on the real kernel: what a solid keeps when it
//! moves, turns or changes size (its labels, its surfaces' kinds, its volume)
//! and what a non-uniform scale turns its surfaces into.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, EdgeBlend, EdgeQuery, Error, FaceGeometry, FilletOptions, Kernel, Label, Surface,
    TessellationParams, Tolerance,
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

/// A 400 × 300 × 250 box centred on the origin with its vertical edges
/// rounded (cylinders), a 25 mm hole down through its top whose rim is
/// chamfered (a NURBS), a countersunk bore right through (a cone), and a ball
/// cut from one top corner (a sphere). Every
/// face and edge is labelled by its index, so labels can be followed through
/// a transform.
fn part(k: &Kernel) -> Brep {
    let block = Brep::cuboid(k, [-200.0, -150.0, -125.0], [400.0, 300.0, 250.0]).unwrap();
    let vertical = EdgeQuery::ParallelTo {
        direction: [0.0, 0.0, 1.0],
    };
    let block = block
        .fillet(
            &[EdgeBlend::round(vertical, 20.0)],
            &FilletOptions::default(),
        )
        .unwrap();
    let hole = Brep::cylinder(k, [30.0, 20.0, 40.0], 25.0, 200.0).unwrap();
    let ball = Brep::sphere(k, [-200.0, -150.0, 125.0], 60.0).unwrap();
    let bore = Brep::cylinder(k, [-80.0, 40.0, -200.0], 15.0, 400.0).unwrap();
    let countersink = Brep::cone(k, [-80.0, 40.0, 100.0], 10.0, 45.0, 35.0).unwrap();
    let part = [ball, bore, countersink]
        .into_iter()
        .fold(block.difference(hole).unwrap(), |part, tool| {
            part.difference(tool).unwrap()
        });
    let rim = EdgeQuery::OnCircle {
        center: [30.0, 20.0, 125.0],
        axis: [0.0, 0.0, 1.0],
        radius: 25.0,
    };
    let mut part = part.chamfer(&rim, 4.0, &Tolerance::default()).unwrap();
    for face in 0..part.face_count() {
        part.set_face_label(face, &Label::new(1, "face").with_index(face as u32))
            .unwrap();
    }
    for edge in 0..part.edge_count() {
        part.set_edge_label(edge, &Label::new(1, "edge").with_index(edge as u32))
            .unwrap();
    }
    part
}

/// A torus, standing on the origin.
fn ring(k: &Kernel) -> Brep {
    Brep::torus(k, [0.0, 0.0, 0.0], 40.0, 10.0).unwrap()
}

fn kind(surface: &Surface) -> &'static str {
    match surface {
        Surface::Plane { .. } => "plane",
        Surface::Cylinder { .. } => "cylinder",
        Surface::Cone { .. } => "cone",
        Surface::Sphere { .. } => "sphere",
        Surface::Torus { .. } => "torus",
        Surface::Nurbs { .. } => "nurbs",
        Surface::Other { .. } => "other",
    }
}

fn kinds(faces: &[FaceGeometry]) -> Vec<&'static str> {
    faces.iter().map(|f| kind(&f.surface)).collect()
}

fn edge_labels(solid: &Brep) -> Vec<Option<Label>> {
    solid
        .edges()
        .unwrap()
        .into_iter()
        .map(|e| e.label)
        .collect()
}

fn face_labels(faces: &[FaceGeometry]) -> Vec<Option<Label>> {
    faces.iter().map(|f| f.label.clone()).collect()
}

/// The rotation by `deg` degrees about z, row by row.
fn about_z(deg: f64) -> [[f64; 3]; 3] {
    let (s, c) = deg.to_radians().sin_cos();
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

fn apply(rotation: &[[f64; 3]; 3], translation: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let turned = rotation.map(|row| row[0] * p[0] + row[1] * p[1] + row[2] * p[2]);
    [0, 1, 2].map(|i| turned[i] + translation[i])
}

fn assert_near(what: &str, actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    let gap = (0..3)
        .map(|i| (actual[i] - expected[i]).abs())
        .fold(0.0, f64::max);
    assert!(
        gap <= tolerance,
        "{what}: {actual:?} is {gap:e} from {expected:?}"
    );
}

/// Checks a brep's bounding box, which SMLib pads by twice its tolerance.
fn assert_box(what: &str, solid: &Brep, [min, max]: [[f64; 3]; 2]) {
    let [low, high] = solid.bounding_box().unwrap();
    let pad = 2.0 * solid.tolerance().unwrap();
    assert_near(what, low, min.map(|c| c - pad), 1e-9);
    assert_near(what, high, max.map(|c| c + pad), 1e-9);
}

fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

#[test]
fn a_rigid_transform_keeps_labels_surfaces_and_volume() {
    on_kernel(|k| {
        let before = part(k);
        let faces = before.faces().unwrap();
        let volume = before.volume(1e-9).unwrap();
        let rotation = about_z(30.0);
        let translation = [100.0, 0.0, 50.0];

        let mut moved = part(k);
        moved.transform(rotation, translation).unwrap();
        let after = moved.faces().unwrap();
        assert!(moved.is_manifold().unwrap());
        assert_eq!(
            face_labels(&after),
            face_labels(&faces),
            "every face keeps its label"
        );
        assert_eq!(
            edge_labels(&moved),
            edge_labels(&before),
            "every edge keeps its label"
        );
        assert_eq!(kinds(&after), kinds(&faces), "every surface keeps its kind");
        assert!(
            kinds(&faces).contains(&"cone") && kinds(&faces).contains(&"sphere"),
            "{:?}",
            kinds(&faces)
        );
        // The integrator adapts its steps to where the faces are, so the
        // volume agrees to its accuracy, not bit for bit.
        let error = relative(moved.volume(1e-9).unwrap(), volume);
        println!("volume moved: {error:e} relative");
        assert!(error < 1e-8, "{error:e}");
        for (was, is) in faces.iter().zip(&after) {
            // A rigid motion moves a surface's placement, not its (u, v).
            assert_eq!(
                (was.uv_min, was.uv_max),
                (is.uv_min, is.uv_max),
                "face {}",
                was.index
            );
            match (&was.surface, &is.surface) {
                (
                    Surface::Cylinder {
                        origin,
                        axis,
                        radius,
                        ..
                    },
                    Surface::Cylinder {
                        origin: o,
                        axis: a,
                        radius: r,
                        ..
                    },
                ) => {
                    assert_near(
                        "cylinder origin",
                        *o,
                        apply(&rotation, translation, *origin),
                        1e-9,
                    );
                    assert_near(
                        "cylinder axis",
                        *a,
                        apply(&rotation, [0.0; 3], *axis),
                        1e-12,
                    );
                    assert_eq!(radius, r);
                }
                (Surface::Plane { .. }, Surface::Plane { .. }) => {
                    let [n, m] = [was.normal.unwrap(), is.normal.unwrap()];
                    assert_near("normal", m, apply(&rotation, [0.0; 3], n), 1e-12);
                }
                _ => {}
            }
        }
        // A rigid motion changes no face's triangulation, only where it is.
        let params = TessellationParams::default();
        let (was, is) = (
            before.tessellate(&params).unwrap(),
            moved.tessellate(&params).unwrap(),
        );
        assert_eq!(was.indices, is.indices);
        assert_eq!(was.face_ids, is.face_ids);
        let far = was
            .positions
            .iter()
            .zip(&is.positions)
            .map(|(p, q)| {
                let expected = apply(&rotation, translation, p.map(f64::from));
                (0..3)
                    .map(|i| (f64::from(q[i]) - expected[i]).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max);
        assert!(
            far < 1e-3,
            "vertices land {far:e} from where the motion takes them"
        );

        // Later operations find the moved part's edges by label: round the
        // bore's bottom rim, found on the part before it moved.
        let back = EdgeQuery::OnCircle {
            center: [-80.0, 40.0, -125.0],
            axis: [0.0, 0.0, 1.0],
            radius: 15.0,
        };
        let tolerance = Tolerance::default();
        let picked: Vec<Label> = before
            .edges()
            .unwrap()
            .into_iter()
            .filter(|e| back.matches(e, &tolerance))
            .filter_map(|e| e.label)
            .collect();
        assert_eq!(picked.len(), 1);
        let edges = moved.edges().unwrap();
        let found = edges
            .iter()
            .position(|e| e.label.as_ref() == Some(&picked[0]))
            .unwrap();
        let turned = EdgeQuery::OnCircle {
            center: apply(&rotation, translation, [-80.0, 40.0, -125.0]),
            axis: [0.0, 0.0, 1.0],
            radius: 15.0,
        };
        assert!(
            turned.matches(&edges[found], &tolerance),
            "it turned with the part"
        );
        let rounded = moved
            .fillet(
                &[EdgeBlend::round(EdgeQuery::Indices(vec![found]), 3.0)],
                &FilletOptions::default(),
            )
            .unwrap();
        assert!(rounded.is_manifold().unwrap());
        assert_eq!(rounded.face_count(), before.face_count() + 1);
    });
}

#[test]
fn translate_and_rotate_about_a_line() {
    on_kernel(|k| {
        let mut block = Brep::cuboid(k, [0.0, 0.0, 0.0], [40.0, 20.0, 10.0]).unwrap();
        block.translate([5.0, -3.0, 2.0]).unwrap();
        assert_box("moved", &block, [[5.0, -3.0, 2.0], [45.0, 17.0, 12.0]]);

        // A quarter turn about the vertical line through the block's corner
        // at (5, -3) swings the block from +x of it to +y of it.
        block
            .rotate_about([5.0, -3.0, 0.0], [0.0, 0.0, 2.0], 90.0)
            .unwrap();
        assert_box("turned", &block, [[-15.0, -3.0, 2.0], [5.0, 37.0, 12.0]]);
        block
            .rotate_about([5.0, -3.0, 0.0], [0.0, 0.0, 1.0], -90.0)
            .unwrap();
        assert_box(
            "turned back",
            &block,
            [[5.0, -3.0, 2.0], [45.0, 17.0, 12.0]],
        );
        assert!(relative(block.volume(1e-9).unwrap(), 8000.0) < 1e-9);
    });
}

#[test]
fn transforms_that_arent_rigid_are_refused() {
    on_kernel(|k| {
        let mut block = Brep::cuboid(k, [0.0, 0.0, 0.0], [40.0, 20.0, 10.0]).unwrap();
        let bad_frame = |result: Result<(), Error>| match result {
            Err(Error::Kernel { status, .. }) => status.name() == Some("SHIM_ERR_BAD_FRAME"),
            _ => false,
        };
        let stretched = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let sheared = [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mirrored = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]];
        for rotation in [stretched, sheared, mirrored, [[f64::NAN; 3]; 3]] {
            assert!(
                bad_frame(block.transform(rotation, [0.0; 3])),
                "{rotation:?}"
            );
        }
        assert!(matches!(
            block.transform(about_z(10.0), [f64::INFINITY, 0.0, 0.0]),
            Err(Error::InvalidArgument { .. })
        ));
        assert!(matches!(
            block.rotate_about([0.0; 3], [0.0; 3], 10.0),
            Err(Error::InvalidArgument { .. })
        ));
        assert_box("left alone", &block, [[0.0; 3], [40.0, 20.0, 10.0]]);
    });
}

#[test]
fn a_uniform_scale_keeps_surfaces_analytic() {
    on_kernel(|k| {
        let centre = [12.0, -7.0, 30.0];
        for (name, make) in [("part", part as fn(&Kernel) -> Brep), ("ring", ring)] {
            let before = make(k);
            let faces = before.faces().unwrap();
            let volume = before.volume(1e-9).unwrap();
            let [min, max] = before.bounding_box().unwrap();
            let mut scaled = make(k);
            scaled.scale_about(centre, 2.0).unwrap();
            let after = scaled.faces().unwrap();
            assert!(scaled.is_manifold().unwrap(), "{name}");
            assert_eq!(
                kinds(&after),
                kinds(&faces),
                "{name}: every surface keeps its kind"
            );
            assert_eq!(face_labels(&after), face_labels(&faces), "{name}");
            assert_eq!(edge_labels(&scaled), edge_labels(&before), "{name}");
            // Doubling is exact in binary, and the tolerances double too, so
            // the volume integrates to exactly eight times the old one.
            let error = relative(scaled.volume(1e-9).unwrap(), 8.0 * volume);
            println!("{name} ×2: volume {error:e} from ×8");
            assert!(error < 1e-8, "{name}: {error:e}");
            assert_eq!(
                scaled.tolerance().unwrap(),
                2.0 * before.tolerance().unwrap(),
                "{name}"
            );
            let [low, high] = scaled.bounding_box().unwrap();
            let grown = |p: [f64; 3]| [0, 1, 2].map(|i| centre[i] + 2.0 * (p[i] - centre[i]));
            assert_near(name, low, grown(min), 1e-9);
            assert_near(name, high, grown(max), 1e-9);
            for (was, is) in faces.iter().zip(&after) {
                match (&was.surface, &is.surface) {
                    (Surface::Cylinder { radius, .. }, Surface::Cylinder { radius: r, .. })
                    | (Surface::Sphere { radius, .. }, Surface::Sphere { radius: r, .. }) => {
                        assert_eq!(*r, 2.0 * radius, "{name}")
                    }
                    (
                        Surface::Cone { half_angle, .. },
                        Surface::Cone {
                            half_angle: angle, ..
                        },
                    ) => assert_eq!(angle, half_angle, "{name}"),
                    (
                        Surface::Torus {
                            major_radius,
                            minor_radius,
                            ..
                        },
                        Surface::Torus {
                            major_radius: major,
                            minor_radius: minor,
                            ..
                        },
                    ) => assert_eq!([*major, *minor], [2.0 * major_radius, 2.0 * minor_radius]),
                    _ => {}
                }
            }
        }

        // Shrinking works the same way, and the scaled part takes later
        // operations: a boolean and a fillet.
        let mut shrunk = part(k);
        let volume = shrunk.volume(1e-9).unwrap();
        shrunk.scale_about([0.0; 3], 0.5).unwrap();
        assert!(relative(shrunk.volume(1e-9).unwrap(), volume / 8.0) < 1e-8);
        let cut = shrunk
            .difference(Brep::cylinder(k, [60.0, -60.0, -100.0], 5.0, 200.0).unwrap())
            .unwrap();
        let rounded = cut
            .fillet(
                &[EdgeBlend::round(
                    EdgeQuery::OnCircle {
                        center: [60.0, -60.0, 62.5],
                        axis: [0.0, 0.0, 1.0],
                        radius: 5.0,
                    },
                    1.0,
                )],
                &FilletOptions::default(),
            )
            .unwrap();
        assert!(rounded.is_manifold().unwrap());
    });
}

#[test]
fn a_non_uniform_scale_stretches_the_part() {
    on_kernel(|k| {
        let before = part(k);
        let volume = before.volume(1e-9).unwrap();
        let mut stretched = part(k);
        stretched
            .scale_axes_keeping([0.0; 3], [2.0, 1.0, 1.0], false)
            .unwrap();
        let faces = stretched.faces().unwrap();
        assert!(stretched.is_manifold().unwrap());
        assert_eq!(face_labels(&faces), face_labels(&before.faces().unwrap()));
        assert_eq!(edge_labels(&stretched), edge_labels(&before));
        assert!(relative(stretched.volume(1e-9).unwrap(), 2.0 * volume) < 1e-6);
        // A NURBS's bounding box is its control points', so measure the mesh.
        let mesh = stretched
            .tessellate(&TessellationParams::default())
            .unwrap();
        let (mut min, mut max) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in &mesh.positions {
            for i in 0..3 {
                min[i] = min[i].min(f64::from(p[i]));
                max[i] = max[i].max(f64::from(p[i]));
            }
        }
        assert_near("min", min, [-400.0, -150.0, -125.0], 1e-3);
        assert_near("max", max, [400.0, 150.0, 125.0], 1e-3);
        assert!(faces.iter().filter(|f| f.normal.is_some()).count() >= 6);
    });
}

#[test]
fn scales_that_arent_positive_are_refused() {
    on_kernel(|k| {
        let mut block = Brep::cuboid(k, [0.0, 0.0, 0.0], [40.0, 20.0, 10.0]).unwrap();
        for factor in [0.0, -2.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                block.scale_about([0.0; 3], factor),
                Err(Error::InvalidArgument { .. })
            ));
        }
        assert!(matches!(
            block.scale_axes_keeping([0.0; 3], [1.0, -1.0, 1.0], false),
            Err(Error::InvalidArgument { .. })
        ));
        assert!(matches!(
            block.scale_about([f64::NAN, 0.0, 0.0], 2.0),
            Err(Error::InvalidArgument { .. })
        ));
        assert_box("left alone", &block, [[0.0; 3], [40.0, 20.0, 10.0]]);
    });
}

/// The tangent of a cone face's half-angle, or `None` for another face.
fn cone_slope(face: &FaceGeometry) -> Option<f64> {
    match face.surface {
        Surface::Cone { half_angle, .. } => Some(half_angle.tan()),
        _ => None,
    }
}

#[test]
fn a_stretch_along_the_cylinders_keeps_them_analytic() {
    on_kernel(|k| {
        let before = part(k);
        let faces = before.faces().unwrap();
        let volume = before.volume(1e-9).unwrap();
        let mut stretched = part(k);
        let report = stretched
            .scale_axes_keeping([0.0; 3], [1.0, 1.0, 2.0], false)
            .unwrap();
        let after = stretched.faces().unwrap();
        assert!(stretched.is_manifold().unwrap());
        // Every plane, cylinder and the cone stays; the ball doesn't, and the
        // chamfer was a NURBS already.
        assert_eq!(
            (report.faces_nurbs, report.faces_were_nurbs),
            (1, 1),
            "{report:?}"
        );
        assert_eq!(report.faces_kept as usize, faces.len() - 2);
        let expected: Vec<&str> = kinds(&faces)
            .into_iter()
            .map(|kind| if kind == "sphere" { "nurbs" } else { kind })
            .collect();
        assert_eq!(kinds(&after), expected);
        assert_eq!(face_labels(&after), face_labels(&faces));
        assert_eq!(edge_labels(&stretched), edge_labels(&before));
        assert!(relative(stretched.volume(1e-9).unwrap(), 2.0 * volume) < 1e-8);
        for (was, is) in faces.iter().zip(&after) {
            if let (Surface::Cylinder { radius, .. }, Surface::Cylinder { radius: r, .. }) =
                (&was.surface, &is.surface)
            {
                assert!((r - radius).abs() < 1e-9, "radius {radius} became {r}");
            }
            // Twice as tall, as wide: half as steep.
            if let (Some(slope), Some(now)) = (cone_slope(was), cone_slope(is)) {
                assert!(
                    (now.abs() - slope.abs() / 2.0).abs() < 1e-9,
                    "{slope} became {now}"
                );
            }
        }
        // The bore's bottom rim, now 250 below the origin, still rounds.
        let rim = EdgeQuery::OnCircle {
            center: [-80.0, 40.0, -250.0],
            axis: [0.0, 0.0, 1.0],
            radius: 15.0,
        };
        let rounded = stretched
            .try_clone()
            .unwrap()
            .fillet(&[EdgeBlend::round(rim, 2.0)], &FilletOptions::default())
            .unwrap();
        assert!(rounded.is_manifold().unwrap());
    });
}

#[test]
fn a_stretch_across_the_cylinders_keeps_the_planes() {
    on_kernel(|k| {
        let before = part(k);
        let faces = before.faces().unwrap();
        let volume = before.volume(1e-9).unwrap();
        let planes = kinds(&faces)
            .iter()
            .filter(|&&kind| kind == "plane")
            .count();
        let cylinders = kinds(&faces)
            .iter()
            .filter(|&&kind| kind == "cylinder")
            .count();
        for extrusions in [false, true] {
            let mut stretched = part(k);
            let report = stretched
                .scale_axes_keeping([0.0; 3], [2.0, 1.0, 1.0], extrusions)
                .unwrap();
            let after = stretched.faces().unwrap();
            let now = kinds(&after);
            assert!(stretched.is_manifold().unwrap());
            assert_eq!(now.iter().filter(|&&kind| kind == "plane").count(), planes);
            // Stretched cylinders are exact ellipses swept along a line, or
            // NURBS; the cone and the ball are NURBS.
            let swept = now.iter().filter(|&&kind| kind == "other").count();
            assert_eq!(swept, if extrusions { cylinders } else { 0 }, "{report:?}");
            assert_eq!(report.faces_extruded as usize, swept);
            assert_eq!(report.faces_kept as usize, planes);
            assert!(report.curves_stretched > 0, "{report:?}");
            assert_eq!(face_labels(&after), face_labels(&faces));
            assert_eq!(edge_labels(&stretched), edge_labels(&before));
            let error = relative(stretched.volume(1e-9).unwrap(), 2.0 * volume);
            assert!(error < 1e-8, "{error:e}");

            // Later operations take it: a hole down through a stretched
            // corner round, and a fillet on the bore's rim, now an ellipse.
            let hole = Brep::cylinder(k, [370.0, 135.0, -200.0], 10.0, 400.0).unwrap();
            let cut = stretched.try_clone().unwrap().difference(hole).unwrap();
            assert!(cut.is_manifold().unwrap());
            let rim = EdgeQuery::MidpointOnPlane {
                point: [0.0, 0.0, -125.0],
                normal: [0.0, 0.0, 1.0],
            }
            .and(EdgeQuery::Not(Box::new(EdgeQuery::ParallelTo {
                direction: [1.0, 0.0, 0.0],
            })))
            .and(EdgeQuery::Not(Box::new(EdgeQuery::ParallelTo {
                direction: [0.0, 1.0, 0.0],
            })));
            let ellipse = stretched
                .edges()
                .unwrap()
                .into_iter()
                .filter(|e| rim.matches(e, &Tolerance::default()))
                .find(|e| {
                    let x = e.samples.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
                    (x - (-190.0)).abs() < 1e-6
                })
                .expect("the bore's bottom rim")
                .index;
            let rounded = stretched
                .try_clone()
                .unwrap()
                .fillet(
                    &[EdgeBlend::round(EdgeQuery::Indices(vec![ellipse]), 2.0)],
                    &FilletOptions::default(),
                )
                .unwrap();
            assert!(rounded.is_manifold().unwrap());
        }
    });
}

#[test]
fn a_scale_keeps_its_trims_so_tessellating_is_repeatable() {
    on_kernel(|k| {
        let params = TessellationParams::default();
        let mesh_after = |factors: [f64; 3], volume_first: bool| {
            let mut solid = part(k);
            solid.scale_axes_keeping([0.0; 3], factors, true).unwrap();
            if volume_first {
                solid.volume(1e-9).unwrap();
            }
            solid.tessellate(&params).unwrap()
        };
        // SmBrep::Transform drops the UV trim curves under a non-uniform
        // scale, and the first call that needs them rebuilds them its own
        // way, so what the mesh is depended on what was asked first.
        for factors in [[2.0, 1.0, 1.0], [1.5; 3], [1.0, 0.5, 3.0]] {
            assert_eq!(mesh_after(factors, true), mesh_after(factors, false));
        }
    });
}

#[test]
fn keeping_analytics_under_equal_factors_keeps_every_surface() {
    on_kernel(|k| {
        for make in [part as fn(&Kernel) -> Brep, ring] {
            let before = make(k);
            let mut scaled = make(k);
            let report = scaled
                .scale_axes_keeping([12.0, -7.0, 30.0], [2.0; 3], false)
                .unwrap();
            assert_eq!(report.faces_nurbs, 0, "{report:?}");
            assert_eq!(
                kinds(&scaled.faces().unwrap()),
                kinds(&before.faces().unwrap())
            );
            let error = relative(
                scaled.volume(1e-9).unwrap(),
                8.0 * before.volume(1e-9).unwrap(),
            );
            assert!(error < 1e-8, "{error:e}");
        }
    });
}

#[test]
fn ellipses_swept_along_a_line_read_back_from_a_file() {
    on_kernel(|k| {
        let mut stretched = part(k);
        stretched
            .scale_axes_keeping([0.0; 3], [2.0, 1.0, 1.0], true)
            .unwrap();
        let path = std::env::temp_dir().join(format!("smlib-transform-{}.smb", std::process::id()));
        stretched.write_file(&path).unwrap();
        let back = Brep::read_file(k, &path);
        let _ = std::fs::remove_file(&path);
        let back = back.unwrap();
        assert!(back.is_manifold().unwrap());
        assert_eq!(
            kinds(&back.faces().unwrap()),
            kinds(&stretched.faces().unwrap())
        );
        let error = relative(back.volume(1e-9).unwrap(), stretched.volume(1e-9).unwrap());
        assert!(error < 1e-9, "{error:e}");
    });
}
