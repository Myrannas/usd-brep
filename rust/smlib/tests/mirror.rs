//! Reflections on the real kernel: what a solid keeps when it is mirrored
//! across a plane (its labels, its surfaces' kinds, its volume, its outward
//! normals) and how the mirrored copy joins the original.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, EdgeBlend, EdgeQuery, FaceGeometry, FilletOptions, Kernel, Label, Mesh, Surface,
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

/// The transform tests' part: a 400 × 300 × 250 box centred on the origin
/// with its vertical edges rounded (cylinders), a 25 mm hole down through its
/// top whose rim is chamfered (a NURBS), a countersunk bore right through (a
/// cone), and a ball cut from one top corner (a sphere). Every face and edge
/// is labelled by its index.
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
    label(&mut part);
    part
}

fn label(solid: &mut Brep) {
    for face in 0..solid.face_count() {
        solid
            .set_face_label(face, &Label::new(1, "face").with_index(face as u32))
            .unwrap();
    }
    for edge in 0..solid.edge_count() {
        solid
            .set_edge_label(edge, &Label::new(1, "edge").with_index(edge as u32))
            .unwrap();
    }
}

/// A torus, standing on the origin.
fn ring(k: &Kernel) -> Brep {
    let mut ring = Brep::torus(k, [0.0, 0.0, 0.0], 40.0, 10.0).unwrap();
    label(&mut ring);
    ring
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

fn face_labels(faces: &[FaceGeometry]) -> Vec<Option<Label>> {
    faces.iter().map(|f| f.label.clone()).collect()
}

fn edge_labels(solid: &Brep) -> Vec<Option<Label>> {
    solid
        .edges()
        .unwrap()
        .into_iter()
        .map(|e| e.label)
        .collect()
}

fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

/// The reflection of p across the plane through `point` with unit `normal`.
fn reflect(point: [f64; 3], normal: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let d: f64 = (0..3).map(|i| (p[i] - point[i]) * normal[i]).sum();
    [0, 1, 2].map(|i| p[i] - 2.0 * d * normal[i])
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let length = v.iter().map(|c| c * c).sum::<f64>().sqrt();
    v.map(|c| c / length)
}

/// The volume a mesh encloses, from its triangles' winding: positive when
/// they run counter-clockwise seen from outside.
fn enclosed(mesh: &Mesh) -> f64 {
    mesh.triangles()
        .map(|[a, b, c]| {
            let [p, q, r] = [a, b, c].map(|i| mesh.positions[i as usize].map(f64::from));
            (p[0] * (q[1] * r[2] - q[2] * r[1]) - p[1] * (q[0] * r[2] - q[2] * r[0])
                + p[2] * (q[0] * r[1] - q[1] * r[0]))
                / 6.0
        })
        .sum()
}

/// Compares a mirrored body's mesh with its source's, moved by `map`: each
/// mirrored vertex is matched to the nearest moved source vertex on the same
/// face. Returns the farthest match, and how many mirrored triangles have no
/// source triangle with the same corners, and how many of those that match
/// wind the same way as their source (rather than the other way round, as
/// a reflection turns them).
fn compare_meshes(
    source: &Mesh,
    mirrored: &Mesh,
    map: impl Fn([f64; 3]) -> [f64; 3],
) -> (f64, usize, usize) {
    use std::collections::HashMap;
    let mut by_face: HashMap<u32, Vec<u32>> = HashMap::new();
    let mut face_of = vec![u32::MAX; source.positions.len()];
    for (t, &face) in source.triangles().zip(&source.face_ids) {
        for i in t {
            if face_of[i as usize] == u32::MAX {
                face_of[i as usize] = face;
                by_face.entry(face).or_default().push(i);
            }
        }
    }
    // A seam's two copies of a vertex count as one.
    let mut canonical: Vec<u32> = (0..source.positions.len() as u32).collect();
    let mut seen: HashMap<(u32, [u32; 3]), u32> = HashMap::new();
    for (i, p) in source.positions.iter().enumerate() {
        if face_of[i] != u32::MAX {
            canonical[i] = *seen
                .entry((face_of[i], p.map(f32::to_bits)))
                .or_insert(i as u32);
        }
    }
    let moved: Vec<[f64; 3]> = source
        .positions
        .iter()
        .map(|p| map(p.map(f64::from)))
        .collect();
    let mut matched = vec![u32::MAX; mirrored.positions.len()];
    let mut worst = 0.0f64;
    for (t, &face) in mirrored.triangles().zip(&mirrored.face_ids) {
        for i in t {
            if matched[i as usize] != u32::MAX {
                continue;
            }
            let p = mirrored.positions[i as usize].map(f64::from);
            let (best, gap) = by_face[&face]
                .iter()
                .map(|&j| {
                    let q = moved[j as usize];
                    (j, (0..3).map(|k| (p[k] - q[k]).abs()).fold(0.0, f64::max))
                })
                .min_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap();
            matched[i as usize] = canonical[best as usize];
            worst = worst.max(gap);
        }
    }
    let sorted = |t: [u32; 3]| {
        let mut s = t;
        s.sort();
        s
    };
    let mut wanted: HashMap<[u32; 3], Vec<[u32; 3]>> = HashMap::new();
    for t in source.triangles() {
        let t = t.map(|i| canonical[i as usize]);
        wanted.entry(sorted(t)).or_default().push(t);
    }
    let (mut missing, mut same_winding) = (0, 0);
    for t in mirrored.triangles() {
        let t = t.map(|i| matched[i as usize]);
        match wanted.get_mut(&sorted(t)).and_then(|v| v.pop()) {
            Some(was) => {
                let cyclic = [was, [was[1], was[2], was[0]], [was[2], was[0], was[1]]];
                if cyclic.contains(&t) {
                    same_winding += 1;
                }
            }
            None => missing += 1,
        }
    }
    (worst, missing, same_winding)
}

/// Triangles whose winding disagrees with the kernel's normals at their
/// corners.
fn inside_out(mesh: &Mesh) -> usize {
    mesh.triangles()
        .filter(|&[a, b, c]| {
            let [p, q, r] = [a, b, c].map(|i| mesh.positions[i as usize].map(f64::from));
            let e = [0, 1, 2].map(|i| q[i] - p[i]);
            let f = [0, 1, 2].map(|i| r[i] - p[i]);
            let n = [
                e[1] * f[2] - e[2] * f[1],
                e[2] * f[0] - e[0] * f[2],
                e[0] * f[1] - e[1] * f[0],
            ];
            let normal = [a, b, c]
                .map(|i| mesh.normals[i as usize].map(f64::from))
                .into_iter()
                .fold([0.0; 3], |s, m| [s[0] + m[0], s[1] + m[1], s[2] + m[2]]);
            n[0] * normal[0] + n[1] * normal[1] + n[2] * normal[2] < 0.0
        })
        .count()
}

const PLANES: [([f64; 3], [f64; 3]); 3] = [
    // Beside the part, a separate copy.
    ([300.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
    // Through it, at an angle, with a normal of no particular length.
    ([10.0, 20.0, 30.0], [1.0, 2.0, 3.0]),
    // Under it.
    ([0.0, 0.0, -125.0], [0.0, 0.0, -1.0]),
];

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
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

/// Asserts that `axis` is the reflection of `was`, either way along it.
fn assert_axis(what: &str, axis: [f64; 3], was: [f64; 3], normal: [f64; 3]) {
    let expected = reflect([0.0; 3], normal, was);
    let along = dot(axis, expected);
    assert!(
        (along.abs() - 1.0).abs() < 1e-12,
        "{what}: {axis:?} isn't along {expected:?}"
    );
}

#[test]
fn a_mirror_keeps_every_surface_analytic_and_every_label() {
    on_kernel(|k| {
        let makers = [
            ("part", part as fn(&Kernel) -> Brep),
            ("ring", ring),
            ("bracket", bracket),
        ];
        for (name, make) in makers {
            let before = make(k);
            let faces = before.faces().unwrap();
            let volume = before.volume(1e-9).unwrap();
            for (point, normal) in PLANES {
                let n = unit(normal);
                let mut mirrored = make(k);
                let start = std::time::Instant::now();
                let report = mirrored.mirror(point, normal).unwrap();
                println!(
                    "{name} across {normal:?}: {report:?} in {:?}",
                    start.elapsed()
                );
                let after = mirrored.faces().unwrap();
                assert!(mirrored.is_manifold().unwrap(), "{name}");
                // Nothing analytic became a NURBS: the planes, cylinders,
                // cone, sphere and torus, and every line and circle.
                assert_eq!(report.faces_nurbs, 0, "{name}: {report:?}");
                assert_eq!(report.curves_nurbs, 0, "{name}: {report:?}");
                assert_eq!(kinds(&after), kinds(&faces), "{name}");
                assert_eq!(face_labels(&after), face_labels(&faces), "{name}");
                assert_eq!(edge_labels(&mirrored), edge_labels(&before), "{name}");
                assert_eq!(mirrored.tolerance().unwrap(), before.tolerance().unwrap());
                let error = relative(mirrored.volume(1e-9).unwrap(), volume);
                println!("  volume {error:e} relative");
                assert!(error < 1e-6, "{name}: {error:e}");
                for (was, is) in faces.iter().zip(&after) {
                    // Each surface keeps its parameters.
                    assert_eq!(
                        (was.uv_min, was.uv_max),
                        (is.uv_min, is.uv_max),
                        "{name}: face {}",
                        was.index
                    );
                    let at = |p: [f64; 3]| reflect(point, n, p);
                    match (&was.surface, &is.surface) {
                        (Surface::Plane { .. }, Surface::Plane { .. }) => {
                            // The outward normal is the old one reflected.
                            let [m, o] = [was.normal.unwrap(), is.normal.unwrap()];
                            assert_near("normal", o, reflect([0.0; 3], n, m), 1e-12);
                        }
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
                            assert_axis("cylinder axis", *a, *axis, n);
                            // The new origin lies on the reflected axis.
                            let d = [0, 1, 2].map(|i| o[i] - at(*origin)[i]);
                            let off = [0, 1, 2].map(|i| d[i] - dot(d, *a) * a[i]);
                            assert_near("cylinder axis line", off, [0.0; 3], 1e-9);
                            assert!((r - radius).abs() < 1e-12, "{radius} became {r}");
                        }
                        (
                            Surface::Cone {
                                axis, half_angle, ..
                            },
                            Surface::Cone {
                                axis: a,
                                half_angle: h,
                                ..
                            },
                        ) => {
                            assert_axis("cone axis", *a, *axis, n);
                            assert!((h.abs() - half_angle.abs()).abs() < 1e-12);
                        }
                        (
                            Surface::Sphere { center, radius, .. },
                            Surface::Sphere {
                                center: c,
                                radius: r,
                                ..
                            },
                        ) => {
                            assert_near("sphere centre", *c, at(*center), 1e-9);
                            assert!((r - radius).abs() < 1e-12);
                        }
                        (
                            Surface::Torus {
                                center,
                                axis,
                                major_radius,
                                minor_radius,
                                ..
                            },
                            Surface::Torus {
                                center: c,
                                axis: a,
                                major_radius: major,
                                minor_radius: minor,
                                ..
                            },
                        ) => {
                            assert_near("torus centre", *c, at(*center), 1e-9);
                            assert_axis("torus axis", *a, *axis, n);
                            assert!((major - major_radius).abs() < 1e-12);
                            assert!((minor - minor_radius).abs() < 1e-12);
                        }
                        (Surface::Nurbs { .. }, Surface::Nurbs { .. }) => {}
                        // The bracket's rounded rims: surfaces of revolution.
                        (Surface::Other { kernel_type: a }, Surface::Other { kernel_type: b }) => {
                            assert_eq!(a, b, "{name}")
                        }
                        (was, is) => panic!("{name}: {was:?} became {is:?}"),
                    }
                }
            }
        }
    });
}

#[test]
fn a_mirrored_solid_faces_out() {
    on_kernel(|k| {
        let params = TessellationParams::default();
        for (name, make) in [("part", part as fn(&Kernel) -> Brep), ("ring", ring)] {
            let before = make(k);
            let volume = before.volume(1e-9).unwrap();
            let mesh = before.tessellate(&params).unwrap();
            for (point, normal) in PLANES {
                let n = unit(normal);
                let mut mirrored = make(k);
                mirrored.mirror(point, normal).unwrap();
                let m = mirrored.tessellate(&params).unwrap();
                // Every triangle runs counter-clockwise seen from outside,
                // and the kernel's normals point out.
                assert_eq!(inside_out(&m), 0, "{name}");
                let error = relative(enclosed(&m), volume);
                assert!(error < 1e-2, "{name}: the mesh encloses {error:e} off");
                assert!(relative(enclosed(&m), enclosed(&mesh)) < 1e-6, "{name}");
                // The same triangles as the original's, reflected: every
                // vertex where the reflection takes one of the original's,
                // and every triangle wound the other way round. A few quads
                // split along the other diagonal.
                let (worst, missing, same) = compare_meshes(&mesh, &m, |p| reflect(point, n, p));
                println!(
                    "{name} across {normal:?}: vertices within {worst:e}, {missing} of {} triangles split otherwise",
                    m.triangle_count()
                );
                assert_eq!(m.triangle_count(), mesh.triangle_count(), "{name}");
                assert!(worst < 1e-4, "{name}: {worst:e}");
                assert_eq!(same, 0, "{name}");
                assert!(missing * 20 < m.triangle_count(), "{name}: {missing}");
            }
        }
    });
}

#[test]
fn mirroring_twice_gives_the_solid_back() {
    on_kernel(|k| {
        let params = TessellationParams::default();
        for (name, make) in [("part", part as fn(&Kernel) -> Brep), ("ring", ring)] {
            let before = make(k);
            let mesh = before.tessellate(&params).unwrap();
            for (point, normal) in PLANES {
                let mut back = make(k);
                back.mirror(point, normal).unwrap();
                back.mirror(point, normal).unwrap();
                assert_eq!(
                    kinds(&back.faces().unwrap()),
                    kinds(&before.faces().unwrap())
                );
                let m = back.tessellate(&params).unwrap();
                // The same vertices, where they were, and the same triangles
                // wound the same way, but for one planar face of the part
                // with many boundary points (see a_mirrored_solid_faces_out).
                let (worst, missing, same) = compare_meshes(&mesh, &m, |p| p);
                println!("{name} across {normal:?} and back: {missing} triangles split otherwise");
                assert!(worst < 1e-4, "{name}: {worst:e}");
                assert_eq!(same + missing, m.triangle_count(), "{name}");
                assert!(missing * 20 < m.triangle_count(), "{name}: {missing}");
            }
        }
    });
}

/// A bracket symmetric about x = 0: a 200 × 120 × 60 block centred on the
/// origin with its vertical edges rounded, a hole through either side, a bore
/// and a ball pocket on the plane of symmetry, and rounded rims on the side
/// holes (SMLib makes those surfaces of revolution).
fn bracket(k: &Kernel) -> Brep {
    let block = Brep::cuboid(k, [-100.0, -60.0, -30.0], [200.0, 120.0, 60.0]).unwrap();
    let vertical = EdgeQuery::ParallelTo {
        direction: [0.0, 0.0, 1.0],
    };
    let block = block
        .fillet(
            &[EdgeBlend::round(vertical, 10.0)],
            &FilletOptions::default(),
        )
        .unwrap();
    let tools = [
        Brep::cylinder(k, [-60.0, 0.0, -40.0], 12.0, 80.0).unwrap(),
        Brep::cylinder(k, [60.0, 0.0, -40.0], 12.0, 80.0).unwrap(),
        Brep::cylinder(k, [0.0, 30.0, -40.0], 8.0, 80.0).unwrap(),
        Brep::sphere(k, [0.0, -25.0, 30.0], 18.0).unwrap(),
    ];
    let cut = tools
        .into_iter()
        .fold(block, |part, tool| part.difference(tool).unwrap());
    let rims = [-60.0, 60.0].map(|x| EdgeQuery::OnCircle {
        center: [x, 0.0, 30.0],
        axis: [0.0, 0.0, 1.0],
        radius: 12.0,
    });
    let [a, b] = rims;
    let mut cut = cut
        .fillet(&[EdgeBlend::round(a.or(b), 3.0)], &FilletOptions::default())
        .unwrap();
    label(&mut cut);
    cut
}

fn sorted_kinds(solid: &Brep) -> Vec<&'static str> {
    let mut kinds = kinds(&solid.faces().unwrap());
    kinds.sort();
    kinds
}

#[test]
fn half_mirrored_and_joined_is_the_whole() {
    on_kernel(|k| {
        let whole = bracket(k);
        let volume = whole.volume(1e-9).unwrap();
        assert!(sorted_kinds(&whole).contains(&"other"));
        // The half on +x, mirrored across x = 0 and joined back.
        let keep = Brep::cuboid(k, [0.0, -100.0, -100.0], [200.0, 200.0, 200.0]).unwrap();
        let half = bracket(k).intersection(keep).unwrap();
        let mut other = half.try_clone().unwrap();
        let start = std::time::Instant::now();
        other.mirror([0.0; 3], [1.0, 0.0, 0.0]).unwrap();
        let mirrored = start.elapsed();
        let start = std::time::Instant::now();
        let mut joined = half.try_clone().unwrap().union(other).unwrap();
        let unioned = start.elapsed();
        let pieces = joined.face_count();
        joined.merge_faces().unwrap();
        println!(
            "half {} faces; joined {pieces} faces, {} merged, whole {}; mirror {mirrored:?}, union {unioned:?}",
            joined.face_count(),
            joined.face_count(),
            whole.face_count()
        );
        assert!(joined.is_manifold().unwrap());
        let error = relative(joined.volume(1e-9).unwrap(), volume);
        println!("volume {error:e} from the whole");
        assert!(error < 1e-8, "{error:e}");
        // The faces cut at x = 0 met their mirrored copies and closed up: one
        // bore, one pocket, one top, as in the part made whole.
        // The faces cut at x = 0 met their mirrored copies and closed up:
        // the top, bottom, front and back and the bore are whole again. The
        // ball pocket isn't: SMLib doesn't merge a sphere's pieces across its
        // seam, which the cut put in the +x half, so the pocket comes back as
        // two pieces from each half. Joining the two halves cut from the
        // bracket, no mirror at all, leaves it in three.
        let count = |solid: &Brep, wanted: &str| {
            sorted_kinds(solid)
                .into_iter()
                .filter(|&kind| kind == wanted)
                .count()
        };
        for kind in ["plane", "cylinder", "other"] {
            assert_eq!(count(&joined, kind), count(&whole, kind), "{kind}");
        }
        assert_eq!((count(&whole, "sphere"), count(&joined, "sphere")), (1, 4));
        let left = bracket(k)
            .intersection(Brep::cuboid(k, [-200.0, -100.0, -100.0], [200.0, 200.0, 200.0]).unwrap())
            .unwrap();
        let mut halves = half.try_clone().unwrap().union(left).unwrap();
        halves.merge_faces().unwrap();
        assert_eq!(count(&halves, "sphere"), 3);
    });
}

#[test]
fn later_operations_take_a_mirrored_solid() {
    on_kernel(|k| {
        let mut mirrored = part(k);
        mirrored.mirror([0.0; 3], [1.0, 0.0, 0.0]).unwrap();
        let volume = mirrored.volume(1e-9).unwrap();
        // The bore, now at x = 80, takes a fillet on its bottom rim, found
        // by the label it had before.
        let rim = EdgeQuery::OnCircle {
            center: [-80.0, 40.0, -125.0],
            axis: [0.0, 0.0, 1.0],
            radius: 15.0,
        };
        let tolerance = Tolerance::default();
        let label = part(k)
            .edges()
            .unwrap()
            .into_iter()
            .find(|e| rim.matches(e, &tolerance))
            .and_then(|e| e.label)
            .unwrap();
        let edges = mirrored.edges().unwrap();
        let found = edges
            .iter()
            .position(|e| e.label.as_ref() == Some(&label))
            .unwrap();
        let moved = EdgeQuery::OnCircle {
            center: [80.0, 40.0, -125.0],
            axis: [0.0, 0.0, 1.0],
            radius: 15.0,
        };
        assert!(moved.matches(&edges[found], &tolerance));
        let rounded = mirrored
            .try_clone()
            .unwrap()
            .fillet(
                &[EdgeBlend::round(EdgeQuery::Indices(vec![found]), 3.0)],
                &FilletOptions::default(),
            )
            .unwrap();
        assert!(rounded.is_manifold().unwrap());
        assert!(rounded.volume(1e-9).unwrap() < volume);
        // A chamfer round the mirrored top's edges, a hole through a
        // mirrored round, and a face pushed out.
        let top = EdgeQuery::MidpointOnPlane {
            point: [0.0, 0.0, 125.0],
            normal: [0.0, 0.0, 1.0],
        };
        let chamfered = mirrored
            .try_clone()
            .unwrap()
            .chamfer(
                &top.and(EdgeQuery::Sharp {
                    min_angle_deg: 60.0,
                }),
                2.0,
                &tolerance,
            )
            .unwrap();
        assert!(chamfered.is_manifold().unwrap());
        let hole = Brep::cylinder(k, [190.0, 140.0, -200.0], 6.0, 400.0).unwrap();
        let cut = mirrored.try_clone().unwrap().difference(hole).unwrap();
        assert!(cut.is_manifold().unwrap());
        assert!(cut.volume(1e-9).unwrap() < volume);
        let back = mirrored
            .find_faces(
                &smlib::FaceQuery::Facing {
                    direction: [0.0, 1.0, 0.0],
                },
                &tolerance,
            )
            .unwrap();
        let pushed = mirrored
            .try_clone()
            .unwrap()
            .offset_faces(&back, 5.0)
            .unwrap();
        // As much as the same push on the part as it was.
        let was = part(k);
        let back = was
            .find_faces(
                &smlib::FaceQuery::Facing {
                    direction: [0.0, 1.0, 0.0],
                },
                &tolerance,
            )
            .unwrap();
        let expected = was.offset_faces(&back, 5.0).unwrap().volume(1e-9).unwrap();
        let error = relative(pushed.volume(1e-9).unwrap(), expected);
        assert!(error < 1e-7, "{error:e}");
        // And a copy mirrored across a plane beside the part joins it as a
        // second lump.
        let mut beside = part(k);
        beside.mirror([300.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
        let both = part(k).union(beside).unwrap();
        let error = relative(both.volume(1e-9).unwrap(), 2.0 * volume);
        assert!(error < 1e-8, "{error:e}");
    });
}

#[test]
fn mirrors_across_bad_planes_are_refused() {
    on_kernel(|k| {
        let mut block = Brep::cuboid(k, [0.0, 0.0, 0.0], [40.0, 20.0, 10.0]).unwrap();
        for (point, normal) in [
            ([0.0; 3], [0.0; 3]),
            ([0.0; 3], [f64::NAN, 0.0, 1.0]),
            ([f64::INFINITY, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ] {
            assert!(matches!(
                block.mirror(point, normal),
                Err(smlib::Error::InvalidArgument { .. })
            ));
            assert!(matches!(
                block.mirror_stock(point, normal),
                Err(smlib::Error::InvalidArgument { .. })
            ));
        }
        let [low, high] = block.bounding_box().unwrap();
        assert!(low[0] > -1e-3 && high[0] < 40.0 + 1e-3);
    });
}

#[test]
fn the_stock_mirror_rebuilds_trims_and_turns_curves_into_nurbs() {
    on_kernel(|k| {
        let params = TessellationParams::default();
        for (name, make) in [("part", part as fn(&Kernel) -> Brep), ("ring", ring)] {
            let before = make(k);
            let faces = before.faces().unwrap();
            let volume = before.volume(1e-9).unwrap();
            let mut stock = make(k);
            let start = std::time::Instant::now();
            stock
                .mirror_stock([300.0, 0.0, 0.0], [1.0, 0.0, 0.0])
                .unwrap();
            let took = start.elapsed();
            let after = stock.faces().unwrap();
            assert!(stock.is_manifold().unwrap());
            assert_eq!(face_labels(&after), face_labels(&faces));
            assert!(relative(stock.volume(1e-9).unwrap(), volume) < 1e-8);
            assert_eq!(inside_out(&stock.tessellate(&params).unwrap()), 0);
            let moved = faces
                .iter()
                .zip(&after)
                .filter(|(a, b)| (a.uv_min, a.uv_max) != (b.uv_min, b.uv_max))
                .count();
            println!(
                "{name}: stock {took:?}, kinds {:?}, {moved} domains moved",
                kinds(&after)
            );
            // Mirroring it back with mirror() says what its curves became.
            let report = stock.mirror([300.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
            println!("{name}: {report:?}");
            assert_eq!(report.curves_were_nurbs as usize, before.edge_count());
            if name == "ring" {
                // SMLib finds no torus in the reflected NURBS.
                assert_eq!(kinds(&after), ["other"]);
            } else {
                assert!(moved > 0);
            }
        }
    });
}
