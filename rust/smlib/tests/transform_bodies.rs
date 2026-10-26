//! Every kind of body this crate makes, moved, scaled, stretched, mirrored
//! and combined: what each keeps (its surfaces' kinds, its faces' domains, its
//! labels, its volume, its outward faces), what its triangles do, and how
//! long each takes. Run with `--nocapture` for the table of which surfaces
//! stay analytic under which operation, and the timings.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::collections::BTreeMap;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;
use std::time::{Duration, Instant};

use smlib::{
    BoundaryCurve, Brep, EdgeBlend, EdgeQuery, Error, FaceQuery, FilletOptions, Frame, Kernel,
    Label, Mesh, PlaneCurve, Relationship, RigidMotion, ScaleReport, Surface, TessellationParams,
    Tolerance,
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

fn line(start: [f64; 2], end: [f64; 2]) -> BoundaryCurve {
    BoundaryCurve {
        curve: PlaneCurve::Line { start, end },
        label: None,
    }
}

fn arc(centre: [f64; 2], radius: f64, start_deg: f64, sweep_deg: f64) -> BoundaryCurve {
    BoundaryCurve {
        curve: PlaneCurve::Arc {
            centre,
            radius,
            start_deg,
            sweep_deg,
        },
        label: None,
    }
}

/// Every kind of body: from the primitives, rounded, revolved, piped,
/// drafted, hollowed, with faces moved and offset, and stretched into an
/// ellipse swept along a line. Between them they have every kind of surface
/// the kernel makes. Every face and edge is labelled by its index.
fn bodies(k: &Kernel) -> Vec<(&'static str, Brep)> {
    let flat = Frame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    let tolerance = Tolerance::default();
    let block = || Brep::cuboid(k, [-50.0, -30.0, -20.0], [100.0, 60.0, 40.0]).unwrap();
    let up = |solid: &Brep| {
        solid
            .find_faces(
                &FaceQuery::Facing {
                    direction: [0.0, 0.0, 1.0],
                },
                &tolerance,
            )
            .unwrap()
    };
    let mut bodies = Vec::new();

    // Planes, cylinders along every axis and spheres at the corners.
    let rounded = block()
        .fillet(
            &[EdgeBlend::round(EdgeQuery::All, 5.0)],
            &FilletOptions::default(),
        )
        .unwrap();
    bodies.push(("filleted box", rounded));

    // A boss whose rim is rounded: a fillet on a circular edge.
    let boss = Brep::cylinder(k, [0.0, 0.0, 0.0], 20.0, 30.0)
        .unwrap()
        .fillet(
            &[EdgeBlend::round(
                EdgeQuery::OnCircle {
                    center: [0.0, 0.0, 30.0],
                    axis: [0.0, 0.0, 1.0],
                    radius: 20.0,
                },
                6.0,
            )],
            &FilletOptions::default(),
        )
        .unwrap();
    bodies.push(("rounded boss", boss));

    // Revolved whole: a flat ring, a torus off the axis and a cone.
    let profile = vec![
        line([0.0, 0.0], [20.0, 0.0]),
        arc([20.0, 10.0], 10.0, -90.0, 180.0),
        line([20.0, 20.0], [10.0, 30.0]),
        line([10.0, 30.0], [0.0, 30.0]),
        line([0.0, 30.0], [0.0, 0.0]),
    ];
    let sheet = Brep::planar_faces(k, &flat, &[profile]).unwrap();
    let knob = sheet
        .revolve_face(0, [0.0; 3], [0.0, 1.0, 0.0], 360.0)
        .unwrap();
    bodies.push(("revolved knob", knob.prism));

    // Three quarters of a ball: a sphere and two planes.
    let half_disc = vec![arc([0.0, 0.0], 10.0, -90.0, 180.0), line([0.0, 10.0], [0.0, -10.0])];
    let sheet = Brep::planar_faces(k, &flat, &[half_disc]).unwrap();
    let ball = sheet
        .revolve_face(0, [0.0; 3], [0.0, 1.0, 0.0], 270.0)
        .unwrap();
    bodies.push(("revolved ball", ball.prism));

    // A hollow pipe round a bend: cylinders, tori and its two end rings.
    let path = [
        PlaneCurve::Line {
            start: [0.0, 0.0],
            end: [40.0, 0.0],
        },
        PlaneCurve::Arc {
            centre: [40.0, 20.0],
            radius: 20.0,
            start_deg: -90.0,
            sweep_deg: 90.0,
        },
        PlaneCurve::Line {
            start: [60.0, 20.0],
            end: [60.0, 60.0],
        },
    ];
    let pipe = Brep::pipe(k, &flat, &path, 5.0, Some(3.0)).unwrap();
    bodies.push(("hollow pipe", pipe.solid));

    // Drafted: a disc into a cone, and a D into a NURBS side.
    let disc = vec![arc([0.0, 0.0], 12.0, 0.0, 360.0)];
    let sheet = Brep::planar_faces(k, &flat, &[disc]).unwrap();
    let frustum = sheet.sweep_face_tapered(0, 15.0, 10.0).unwrap();
    bodies.push(("drafted disc", frustum.prism));
    let d = vec![line([-10.0, 0.0], [10.0, 0.0]), arc([0.0, 0.0], 10.0, 0.0, 180.0)];
    let sheet = Brep::planar_faces(k, &flat, &[d]).unwrap();
    let drafted = sheet.sweep_face_tapered(0, 10.0, 10.0).unwrap();
    bodies.push(("drafted D", drafted.prism));

    // Hollowed: a box with rounded verticals open at the top, and one with
    // a chamfered hole (a NURBS chamfer) open at the bottom.
    let rounded = block()
        .fillet(
            &[EdgeBlend::round(
                EdgeQuery::ParallelTo {
                    direction: [0.0, 0.0, 1.0],
                },
                8.0,
            )],
            &FilletOptions::default(),
        )
        .unwrap();
    let open = up(&rounded);
    bodies.push(("shelled box", rounded.shell(2.0, &open).unwrap().brep));
    let chamfered = || {
        let hole = Brep::cylinder(k, [0.0, 0.0, -30.0], 10.0, 60.0).unwrap();
        let rim = EdgeQuery::OnCircle {
            center: [0.0, 0.0, 20.0],
            axis: [0.0, 0.0, 1.0],
            radius: 10.0,
        };
        block()
            .difference(hole)
            .unwrap()
            .chamfer(&rim, 3.0, &tolerance)
            .unwrap()
    };
    let solid = chamfered();
    let bottom = solid
        .find_faces(
            &FaceQuery::Facing {
                direction: [0.0, 0.0, -1.0],
            },
            &tolerance,
        )
        .unwrap();
    bodies.push(("shelled chamfer", solid.shell(2.0, &bottom).unwrap().brep));

    // The chamfer's NURBS offset: an offset surface.
    let solid = chamfered();
    let nurbs: Vec<usize> = solid
        .faces()
        .unwrap()
        .iter()
        .filter(|f| matches!(f.surface, Surface::Nurbs { .. }))
        .map(|f| f.index)
        .collect();
    bodies.push(("offset chamfer", solid.offset_faces(&nurbs, 1.0).unwrap()));

    // A box with its right face turned in.
    let solid = block();
    let right = solid
        .find_faces(
            &FaceQuery::Facing {
                direction: [1.0, 0.0, 0.0],
            },
            &tolerance,
        )
        .unwrap();
    let motion = RigidMotion::turn([50.0, 0.0, -20.0], [0.0, 1.0, 0.0], 15.0);
    bodies.push(("moved face", solid.move_faces(&right, &motion).unwrap().brep));

    // A cylinder stretched across its axis: an ellipse swept along a line.
    let mut oval = Brep::cylinder(k, [0.0, 0.0, 0.0], 10.0, 20.0).unwrap();
    oval.scale_axes_keeping([0.0; 3], [2.0, 1.0, 1.0], true)
        .unwrap();
    bodies.push(("stretched cylinder", oval));

    for (_, solid) in &mut bodies {
        label(solid);
    }
    bodies
}

fn kind(surface: &Surface) -> &'static str {
    match surface {
        Surface::Plane { .. } => "plane",
        Surface::Cylinder { .. } => "cylinder",
        Surface::Cone { .. } => "cone",
        Surface::Sphere { .. } => "sphere",
        Surface::Torus { .. } => "torus",
        Surface::Nurbs { .. } => "nurbs",
        // SmSurfTypes.h.
        Surface::Other { kernel_type: 15020 } => "revolution",
        Surface::Other { kernel_type: 15030 } => "extrusion",
        Surface::Other { kernel_type: 15400 } => "offset",
        Surface::Other { .. } => "other",
    }
}

/// Everything about a solid that a transform should keep or scale.
#[derive(Debug)]
struct Shape {
    kinds: Vec<&'static str>,
    axes: Vec<Option<[f64; 3]>>,
    domains: Vec<([f64; 2], [f64; 2])>,
    face_labels: Vec<Option<Label>>,
    edge_labels: Vec<Option<Label>>,
    volume: f64,
    manifold: bool,
    inside_out: usize,
    enclosed: f64,
}

fn shape(solid: &Brep) -> Shape {
    let faces = solid.faces().unwrap();
    let mesh = solid.tessellate(&TessellationParams::default()).unwrap();
    Shape {
        kinds: faces.iter().map(|f| kind(&f.surface)).collect(),
        axes: faces
            .iter()
            .map(|f| match f.surface {
                Surface::Cylinder { axis, .. } | Surface::Cone { axis, .. } => Some(axis),
                _ => None,
            })
            .collect(),
        domains: faces.iter().map(|f| (f.uv_min, f.uv_max)).collect(),
        face_labels: faces.iter().map(|f| f.label.clone()).collect(),
        edge_labels: solid.edges().unwrap().into_iter().map(|e| e.label).collect(),
        volume: solid.volume(1e-9).unwrap(),
        manifold: solid.is_manifold().unwrap(),
        inside_out: inside_out(&mesh),
        enclosed: enclosed(&mesh),
    }
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

fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

const SCALE: f64 = 1.7;

/// The operations, each with how much it scales the volume.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Operation {
    /// Turned 37° about a skew line and moved.
    Rigid,
    /// Scaled by 1.7 about a point off the body.
    Uniform,
    /// Stretched twice as long along z, then along x.
    StretchZ,
    StretchX,
    /// Mirrored across a skew plane through the body.
    Mirror,
}

const OPERATIONS: [Operation; 5] = [
    Operation::Rigid,
    Operation::Uniform,
    Operation::StretchZ,
    Operation::StretchX,
    Operation::Mirror,
];

impl Operation {
    fn name(self) -> &'static str {
        match self {
            Operation::Rigid => "rigid",
            Operation::Uniform => "uniform",
            Operation::StretchZ => "stretch z",
            Operation::StretchX => "stretch x",
            Operation::Mirror => "mirror",
        }
    }

    fn apply(self, solid: &mut Brep) -> Result<Option<ScaleReport>, Error> {
        match self {
            Operation::Rigid => {
                solid.rotate_about([10.0, -20.0, 5.0], [1.0, 2.0, 3.0], 37.0)?;
                solid.translate([120.0, -35.0, 60.0])?;
                Ok(None)
            }
            Operation::Uniform => solid.scale_about([7.0, 11.0, -13.0], SCALE).map(|()| None),
            Operation::StretchZ => solid
                .scale_axes_keeping([3.0, -4.0, 5.0], [1.0, 1.0, 2.0], false)
                .map(Some),
            Operation::StretchX => solid
                .scale_axes_keeping([3.0, -4.0, 5.0], [2.0, 1.0, 1.0], false)
                .map(Some),
            Operation::Mirror => solid.mirror([4.0, 5.0, 6.0], [1.0, -2.0, 3.0]).map(Some),
        }
    }

    fn volume_factor(self) -> f64 {
        match self {
            Operation::Rigid | Operation::Mirror => 1.0,
            Operation::Uniform => SCALE.powi(3),
            Operation::StretchZ | Operation::StretchX => 2.0,
        }
    }

    /// The axis a stretch lengthens.
    fn stretched(self) -> Option<[f64; 3]> {
        match self {
            Operation::StretchZ => Some([0.0, 0.0, 1.0]),
            Operation::StretchX => Some([1.0, 0.0, 0.0]),
            _ => None,
        }
    }
}

#[test]
fn every_body_keeps_its_labels_volume_and_outward_faces() {
    on_kernel(|k| {
        // (kind before, operation) -> kind after -> faces.
        let mut table: BTreeMap<(&str, &str), BTreeMap<&str, usize>> = BTreeMap::new();
        let mut timings: BTreeMap<&str, Vec<(&str, Duration)>> = BTreeMap::new();
        let mut refused = Vec::new();
        for (name, body) in bodies(k) {
            let before = shape(&body);
            assert!(before.manifold, "{name}");
            println!("{name}: {:?}", {
                let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
                for kind in &before.kinds {
                    *counts.entry(kind).or_default() += 1;
                }
                counts
            });
            for operation in OPERATIONS {
                let what = format!("{name}, {}", operation.name());
                let mut solid = body.try_clone().unwrap();
                let start = Instant::now();
                let applied = operation.apply(&mut solid);
                timings
                    .entry(operation.name())
                    .or_default()
                    .push((name, start.elapsed()));
                if let Err(error) = applied {
                    // Only an offset surface under unequal factors, which has
                    // no exact image, and the brep is left alone.
                    assert!(operation.stretched().is_some(), "{what}: {error}");
                    assert!(before.kinds.contains(&"offset"), "{what}: {error}");
                    assert_eq!(solid.faces().unwrap().len(), before.kinds.len());
                    assert!(relative(solid.volume(1e-9).unwrap(), before.volume) < 1e-9);
                    refused.push(what);
                    continue;
                }
                let after = shape(&solid);
                assert!(after.manifold, "{what}");
                assert_eq!(after.face_labels, before.face_labels, "{what}");
                assert_eq!(after.edge_labels, before.edge_labels, "{what}");
                assert_eq!(after.domains, before.domains, "{what}");
                let error = relative(after.volume, operation.volume_factor() * before.volume);
                assert!(error < 1e-6, "{what}: volume off by {error:e}");
                // Faces out: the mesh's winding encloses the volume, and no
                // more triangles disagree with the kernel's normals than did.
                assert!(
                    relative(after.enclosed, after.volume) < 1e-2,
                    "{what}: the mesh encloses {} of {}",
                    after.enclosed,
                    after.volume
                );
                assert!(after.inside_out <= before.inside_out, "{what}");
                for (face, (&was, &now)) in before.kinds.iter().zip(&after.kinds).enumerate() {
                    *table
                        .entry((was, operation.name()))
                        .or_default()
                        .entry(now)
                        .or_default() += 1;
                    match operation.stretched() {
                        None => assert_eq!(now, was, "{what}: face {face}"),
                        Some(stretched) => {
                            let along = before.axes[face].is_some_and(|axis| {
                                (axis.iter().zip(stretched).map(|(a, b)| a * b).sum::<f64>().abs() - 1.0)
                                    .abs()
                                    < 1e-9
                            });
                            if was == "plane" || along {
                                assert_eq!(now, was, "{what}: face {face}");
                            }
                        }
                    }
                }
            }
        }

        println!("\nWhat each kind of surface becomes:");
        let kinds: Vec<&str> = table.keys().map(|(kind, _)| *kind).collect::<std::collections::BTreeSet<_>>().into_iter().collect();
        print!("{:<12}", "");
        for operation in OPERATIONS {
            print!("{:<26}", operation.name());
        }
        println!();
        for kind in kinds {
            print!("{kind:<12}");
            for operation in OPERATIONS {
                let cell = table.get(&(kind, operation.name())).map_or(String::from("-"), |after| {
                    after
                        .iter()
                        .map(|(now, count)| {
                            if *now == kind {
                                format!("kept {count}")
                            } else {
                                format!("{now} {count}")
                            }
                        })
                        .collect::<Vec<_>>()
                        .join(", ")
                });
                print!("{cell:<26}");
            }
            println!();
        }
        println!("\nRefused, the brep unchanged: {refused:?}");
        println!("\nTimings:");
        for (operation, times) in timings {
            let total: Duration = times.iter().map(|(_, t)| *t).sum();
            let slowest = times.iter().max_by_key(|(_, t)| *t).unwrap();
            println!(
                "  {operation}: {total:?} for {} bodies, slowest {} {:?}",
                times.len(),
                slowest.0,
                slowest.1
            );
        }
    });
}

/// Compares two meshes of one solid: each vertex of `b` is matched to the
/// nearest vertex of `a` on the same face. Returns the farthest match, how
/// many of `b`'s triangles have no triangle of `a` with the same corners (a
/// quad split along its other diagonal), and how many of those that match
/// wind the same way.
fn compare(a: &Mesh, b: &Mesh) -> (f64, usize, usize) {
    use std::collections::HashMap;
    let mut by_face: HashMap<u32, Vec<u32>> = HashMap::new();
    let mut face_of = vec![u32::MAX; a.positions.len()];
    for (t, &face) in a.triangles().zip(&a.face_ids) {
        for i in t {
            if face_of[i as usize] == u32::MAX {
                face_of[i as usize] = face;
                by_face.entry(face).or_default().push(i);
            }
        }
    }
    // A seam's two copies of a vertex count as one.
    let mut canonical: Vec<u32> = (0..a.positions.len() as u32).collect();
    let mut seen: HashMap<(u32, [u32; 3]), u32> = HashMap::new();
    for (i, p) in a.positions.iter().enumerate() {
        if face_of[i] != u32::MAX {
            canonical[i] = *seen
                .entry((face_of[i], p.map(f32::to_bits)))
                .or_insert(i as u32);
        }
    }
    let mut matched = vec![u32::MAX; b.positions.len()];
    let mut worst = 0.0f64;
    for (t, &face) in b.triangles().zip(&b.face_ids) {
        for i in t {
            if matched[i as usize] != u32::MAX {
                continue;
            }
            let p = b.positions[i as usize].map(f64::from);
            let (best, gap) = by_face[&face]
                .iter()
                .map(|&j| {
                    let q = a.positions[j as usize].map(f64::from);
                    (j, (0..3).map(|k| (p[k] - q[k]).abs()).fold(0.0, f64::max))
                })
                .min_by(|x, y| x.1.total_cmp(&y.1))
                .unwrap();
            matched[i as usize] = canonical[best as usize];
            worst = worst.max(gap);
        }
    }
    let sorted = |mut t: [u32; 3]| {
        t.sort();
        t
    };
    let mut wanted: HashMap<[u32; 3], Vec<[u32; 3]>> = HashMap::new();
    for t in a.triangles() {
        let t = t.map(|i| canonical[i as usize]);
        wanted.entry(sorted(t)).or_default().push(t);
    }
    let (mut missing, mut same_winding) = (0, 0);
    for t in b.triangles() {
        let t = t.map(|i| matched[i as usize]);
        match wanted.get_mut(&sorted(t)).and_then(|v| v.pop()) {
            Some(was) => {
                if [was, [was[1], was[2], was[0]], [was[2], was[0], was[1]]].contains(&t) {
                    same_winding += 1;
                }
            }
            None => missing += 1,
        }
    }
    (worst, missing, same_winding)
}

#[test]
fn meshes_move_and_mirror_with_their_solids() {
    on_kernel(|k| {
        let params = TessellationParams::default();
        let rotation = {
            // 37° about (1, 2, 3), row by row.
            let length = 14f64.sqrt();
            let [x, y, z] = [1.0 / length, 2.0 / length, 3.0 / length];
            let (sin, cos) = 37f64.to_radians().sin_cos();
            let t = 1.0 - cos;
            [
                [cos + x * x * t, x * y * t - z * sin, x * z * t + y * sin],
                [y * x * t + z * sin, cos + y * y * t, y * z * t - x * sin],
                [z * x * t - y * sin, z * y * t + x * sin, cos + z * z * t],
            ]
        };
        let translation = [120.0, -35.0, 60.0];
        let (point, normal) = ([4.0, 5.0, 6.0], [1.0, -2.0, 3.0]);
        for (name, body) in bodies(k) {
            let mesh = body.tessellate(&params).unwrap();
            let keys = body.begin_tessellation(&params).unwrap().face_keys().unwrap();

            let mut moved = body.try_clone().unwrap();
            moved.transform(rotation, translation).unwrap();
            let mut carried = mesh.clone();
            carried.transform(rotation, translation).unwrap();
            let fresh = moved.tessellate(&params).unwrap();
            let moved_keys = moved.begin_tessellation(&params).unwrap().face_keys().unwrap();

            let mut mirrored = body.try_clone().unwrap();
            mirrored.mirror(point, normal).unwrap();
            let mut reflected = mesh.clone();
            reflected.mirror(point, normal).unwrap();
            let fresh_mirror = mirrored.tessellate(&params).unwrap();

            let kept = keys
                .iter()
                .zip(&moved_keys)
                .filter(|(a, b)| a.is_some() && a == b)
                .count();
            println!(
                "{name}: moved {:?}, mirrored {:?} of {} triangles; {kept} of {} face keys kept by the move",
                compare(&carried, &fresh),
                compare(&reflected, &fresh_mirror),
                mesh.triangle_count(),
                keys.len()
            );
            // The moved and mirrored meshes have the moved and mirrored
            // solids' own vertices, to within rounding, and their triangles
            // wound the same way, but for some quads split along the other
            // diagonal, where rounding breaks a tie the other way (more of
            // them for a mirror, whose (u, v) runs the other way round). So
            // they are meshes of the moved and mirrored solids, watertight
            // and facing out, but not bit for bit the ones they'd make.
            for (what, carried, fresh) in [
                ("moved", &carried, &fresh),
                ("mirrored", &reflected, &fresh_mirror),
            ] {
                let (gap, differing, same) = compare(carried, fresh);
                assert!(gap < 1e-4, "{name}: the {what} mesh is {gap:e} off");
                assert_eq!(same + differing, fresh.triangle_count(), "{name}: {what}");
                assert!(differing * 5 < mesh.triangle_count(), "{name}: {differing}");
                assert_eq!(carried.faces.len(), fresh.faces.len());
                assert_eq!(carried.edges.len(), fresh.edges.len());
                assert_eq!(inside_out(carried), inside_out(&mesh), "{name}: {what}");
                assert!(relative(enclosed(carried), enclosed(&mesh)) < 1e-6, "{name}: {what}");
            }
            // The moved faces' keys differ, since they hash the geometry
            // where it lies.
            assert_eq!(kept, 0, "{name}");
        }
    });
}

#[test]
fn copies_apart_combine_and_copies_that_touch_are_refused() {
    on_kernel(|k| {
        for (name, body) in bodies(k) {
            let faces = body.face_count();
            let volume = body.volume(1e-9).unwrap();
            // The kernel's bounding box is its faces' surfaces', which can
            // stand off a rounded body; the mesh's vertices lie on it, and
            // reach no further than it does.
            let mesh = body.tessellate(&TessellationParams::default()).unwrap();
            let along = || mesh.positions.iter().map(|p| f64::from(p[0]));
            let (low, high) = (along().fold(f64::INFINITY, f64::min), along().fold(f64::NEG_INFINITY, f64::max));
            let width = high - low;
            let copy = |offset: f64| {
                let mut copy = body.try_clone().unwrap();
                copy.translate([offset, 0.0, 0.0]).unwrap();
                copy
            };

            // Apart, and mirrored apart.
            let mut pair = body.try_clone().unwrap();
            let start = Instant::now();
            pair.combine(&copy(width + 10.0)).unwrap();
            let took = start.elapsed();
            let mut image = body.try_clone().unwrap();
            image.mirror([high + width + 15.0, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
            pair.combine(&image).unwrap();
            assert!(pair.is_manifold().unwrap(), "{name}");
            assert_eq!(pair.face_count(), 3 * faces, "{name}");
            let error = relative(pair.volume(1e-9).unwrap(), 3.0 * volume);
            assert!(error < 1e-6, "{name}: {error:e}");
            let labels: Vec<_> = pair.faces().unwrap().into_iter().map(|f| f.label).collect();
            assert_eq!(labels[..faces], labels[faces..2 * faces], "{name}");
            assert_eq!(labels[..faces], labels[2 * faces..], "{name}");

            // Touching its mirror image across a plane on its far side,
            // overlapping a copy moved 1 mm, and on itself: refused, the
            // brep unchanged.
            let mut touching = body.try_clone().unwrap();
            touching.mirror([high, 0.0, 0.0], [1.0, 0.0, 0.0]).unwrap();
            for (what, other) in [
                ("touching", touching),
                ("overlapping", copy(1.0)),
                ("on itself", copy(0.0)),
            ] {
                let mut solid = body.try_clone().unwrap();
                let start = Instant::now();
                let combined = solid.combine(&other);
                let refused_in = start.elapsed();
                assert!(
                    matches!(&combined, Err(Error::Kernel { status, .. }) if status.name() == Some("SHIM_ERR_NOT_APART")),
                    "{name}, {what}: {combined:?}"
                );
                assert_eq!(solid.face_count(), faces, "{name}, {what}");
                println!("{name}: {what} refused in {refused_in:?}");
            }
            println!("{name}: combined apart in {took:?}");
        }
    });
}

#[test]
fn the_relationship_query_is_exact_for_planes_and_spheres_and_cautious_for_cylinders() {
    on_kernel(|k| {
        let block = |x: f64| Brep::cuboid(k, [x, 0.0, 0.0], [10.0, 10.0, 10.0]).unwrap();
        let ball = |x: f64| Brep::sphere(k, [x, 0.0, 0.0], 5.0).unwrap();
        // The crate's bolt: a circle of radius 8 raised 6.
        let bolt = |x: f64, y: f64| Brep::cylinder(k, [x, y, 140.0], 8.0, 6.0).unwrap();
        for gap in [0.01, 0.1, 1.0, 10.0, 300.0] {
            for (what, found) in [
                ("blocks", block(0.0).relationship(&block(10.0 + gap)).unwrap()),
                ("balls", ball(0.0).relationship(&ball(10.0 + gap)).unwrap()),
                (
                    "bolts one above the other",
                    bolt(0.0, 0.0)
                        .relationship(&Brep::cylinder(k, [0.0, 0.0, 146.0 + gap], 8.0, 6.0).unwrap())
                        .unwrap(),
                ),
            ] {
                let Relationship::Apart { gap: found } = found else {
                    panic!("{what} {gap} apart: {found:?}")
                };
                assert!((found - gap).abs() < 1e-9, "{what} {gap} apart: {found}");
            }
        }
        // Side by side, a cylinder's distance falls short: never more than
        // the truth, so what reads as apart is.
        for gap in [0.01, 1.0, 5.0, 10.0, 300.0] {
            let start = Instant::now();
            let found = bolt(0.0, 0.0).relationship(&bolt(16.0 + gap, 0.0)).unwrap();
            let took = start.elapsed();
            match found {
                Relationship::Apart { gap: found } => assert!(found <= gap, "{gap}: {found}"),
                Relationship::Touching => {}
                other => panic!("{gap}: {other:?}"),
            }
            println!("bolts {gap} apart side by side: {found:?} in {took:?}");
        }
        let start = Instant::now();
        let found = bolt(0.0, 0.0).relationship(&bolt(16.0, 0.0)).unwrap();
        println!("bolts touching: {found:?} in {:?}", start.elapsed());
        assert!(!matches!(found, Relationship::Apart { .. }));
        let start = Instant::now();
        let found = bolt(0.0, 0.0).relationship(&bolt(10.0, 0.0)).unwrap();
        println!("bolts overlapping: {found:?} in {:?}", start.elapsed());
        assert_eq!(found, Relationship::Overlapping);
        let start = Instant::now();
        let found = bolt(165.0, 115.0).relationship(&bolt(-165.0, 115.0)).unwrap();
        println!("the crate's bolts: {found:?} in {:?}", start.elapsed());
        assert!(matches!(found, Relationship::Apart { .. }));
        // One inside the other.
        let inner = Brep::cuboid(k, [2.0, 2.0, 2.0], [6.0, 6.0, 6.0]).unwrap();
        assert_eq!(block(0.0).relationship(&inner).unwrap(), Relationship::Contains { gap: 2.0 });
        assert_eq!(inner.relationship(&block(0.0)).unwrap(), Relationship::Inside { gap: 2.0 });
    });
}

#[test]
fn a_combined_brep_is_asked_about_all_of_itself() {
    on_kernel(|k| {
        let bolt = |x: f64| Brep::cylinder(k, [x, 0.0, 140.0], 8.0, 6.0).unwrap();
        let mut tool = bolt(0.0);
        let far = bolt(100.0);
        // Built now, the query's spatial trees cover the first bolt alone.
        assert!(matches!(tool.relationship(&far).unwrap(), Relationship::Apart { .. }));
        tool.combine(&bolt(50.0)).unwrap();
        let Relationship::Apart { gap } = tool.relationship(&far).unwrap() else {
            panic!()
        };
        assert!(gap < 50.0 - 16.0 + 1e-9, "{gap}");
        assert!(tool.combine(&bolt(60.0)).is_err());
    });
}

#[test]
fn a_solid_inside_another_or_its_cavity_isnt_combined() {
    on_kernel(|k| {
        let cube = |at: f64, size: f64| Brep::cuboid(k, [at; 3], [size; 3]).unwrap();
        let not_apart = |result: Result<(), Error>| {
            matches!(result, Err(Error::Kernel { status, .. }) if status.name() == Some("SHIM_ERR_NOT_APART"))
        };
        // In its material.
        let mut solid = cube(0.0, 100.0);
        assert_eq!(solid.relationship(&cube(10.0, 5.0)).unwrap(), Relationship::Contains { gap: 10.0 });
        assert!(not_apart(solid.combine(&cube(10.0, 5.0))));
        // In a sealed cavity: apart to the kernel's query, but not to
        // combine, since the copy belongs in the cavity's region.
        let mut hollow = cube(0.0, 100.0).difference(cube(20.0, 60.0)).unwrap();
        assert!(matches!(hollow.relationship(&cube(40.0, 20.0)).unwrap(), Relationship::Apart { .. }));
        assert!(not_apart(hollow.combine(&cube(40.0, 20.0))));
        assert_eq!(hollow.face_count(), 12);
        // Round it, apart.
        hollow.combine(&cube(150.0, 20.0)).unwrap();
        assert_eq!(hollow.face_count(), 18);
        assert!(solid.is_manifold().unwrap() && hollow.is_manifold().unwrap());
    });
}

#[test]
fn patterns_combine_copy_by_copy() {
    on_kernel(|k| {
        let bolt = |[x, y]: [f64; 2]| Brep::cylinder(k, [x, y, 140.0], 8.0, 6.0).unwrap();
        let round = |n: usize, radius: f64| -> Vec<[f64; 2]> {
            (0..n)
                .map(|i| {
                    let (sin, cos) = (i as f64 * 360.0 / n as f64).to_radians().sin_cos();
                    [radius * cos, radius * sin]
                })
                .collect()
        };
        let patterns = [
            ("the crate's 4 bolts", vec![[165.0, 115.0], [-165.0, 115.0], [165.0, -115.0], [-165.0, -115.0]]),
            ("12 round a circle", round(12, 100.0)),
            ("100 round a circle", round(100, 600.0)),
            ("a 10 by 10 grid", (0..100).map(|i| [30.0 * (i % 10) as f64, 30.0 * (i / 10) as f64]).collect()),
        ];
        for (name, places) in patterns {
            let copies: Vec<Brep> = places.iter().map(|&p| bolt(p)).collect();
            let start = Instant::now();
            let mut tool = copies[0].try_clone().unwrap();
            let mut slowest = Duration::ZERO;
            for copy in &copies[1..] {
                let one = Instant::now();
                tool.combine(copy).unwrap();
                slowest = slowest.max(one.elapsed());
            }
            let combined = start.elapsed();
            assert_eq!(tool.face_count(), 3 * places.len());
            let plate = || Brep::cuboid(k, [-700.0, -700.0, 130.0], [1400.0, 1400.0, 10.0]).unwrap();
            let start = Instant::now();
            let once = plate().union(tool).unwrap();
            let union = start.elapsed();
            let start = Instant::now();
            let each = copies.into_iter().fold(plate(), |plate, copy| plate.union(copy).unwrap());
            let in_turn = start.elapsed();
            assert_eq!(once.face_count(), each.face_count());
            assert!(relative(once.volume(1e-9).unwrap(), each.volume(1e-9).unwrap()) < 1e-9);
            println!(
                "{name}: combined in {combined:?} (slowest copy {slowest:?}), then one union in {union:?}; \
                 a union per copy {in_turn:?}"
            );
        }
    });
}
