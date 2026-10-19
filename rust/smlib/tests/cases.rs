//! The demo's cases on a 100 × 60 × 40 mm box centred on the origin, Z up,
//! in kernel millimetres: four booleans, then fillets and chamfers on edges
//! picked by geometric queries, then offset faces checked against the same
//! part built directly at the new size.
//!
//! The test harness runs each test on its own thread, but the kernel belongs
//! to one thread, so every test sends its kernel work to a shared kernel
//! thread and checks the plain-data results on its own.

use std::collections::HashSet;
use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, Continuity, EdgeBlend, EdgeQuery, Error, FaceQuery, FilletOptions, Frame, Kernel, Label,
    Mesh, Surface, SurfaceKind, TessellationParams, Tolerance,
};

const L: f64 = 100.0;
const W: f64 = 60.0;
const H: f64 = 40.0;
const BOX_VOLUME: f64 = L * W * H;
const RELATIVE_ACCURACY: f64 = 1e-6;

type Job = Box<dyn FnOnce(&Kernel) + Send>;

/// Runs `f` on the shared kernel thread and returns its result.
fn on_kernel<R: Send + 'static>(f: impl FnOnce(&Kernel) -> R + Send + 'static) -> R {
    static JOBS: OnceLock<mpsc::Sender<Job>> = OnceLock::new();
    let jobs = JOBS.get_or_init(|| {
        let (tx, rx) = mpsc::channel::<Job>();
        std::thread::Builder::new()
            .name("smlib-kernel".into())
            .spawn(move || {
                let kernel = Kernel::init().expect("kernel init");
                for job in rx {
                    // A failed job drops its reply sender; the test sees that.
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

/// What a case reports back from the kernel thread.
struct Built {
    manifold: bool,
    volume: f64,
    brep_faces: usize,
    brep_edges: usize,
    mesh: Mesh,
}

fn build(make: impl FnOnce(&Kernel) -> Result<Brep, Error> + Send + 'static) -> Built {
    on_kernel(move |k| {
        let brep = make(k).expect("building the case");
        Built {
            manifold: brep.is_manifold().unwrap(),
            volume: brep.volume(RELATIVE_ACCURACY).unwrap(),
            brep_faces: brep.face_count(),
            brep_edges: brep.edge_count(),
            mesh: brep.tessellate(&TessellationParams::default()).unwrap(),
        }
    })
}

fn block(k: &Kernel) -> Result<Brep, Error> {
    Brep::cuboid(k, [-L / 2.0, -W / 2.0, -H / 2.0], [L, W, H])
}

fn assert_volume(name: &str, actual: f64, expected: f64) {
    let error = (actual - expected).abs() / expected;
    println!(
        "{name}: volume {actual:.3} vs {expected:.3} ({:.5}%)",
        error * 100.0
    );
    assert!(
        error < 1e-3,
        "{name}: volume {actual} is {:.4}% off {expected}",
        error * 100.0
    );
}

fn sub(a: [f32; 3], b: [f32; 3]) -> [f64; 3] {
    [
        (a[0] - b[0]) as f64,
        (a[1] - b[1]) as f64,
        (a[2] - b[2]) as f64,
    ]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn widen(v: [f32; 3]) -> [f64; 3] {
    [v[0] as f64, v[1] as f64, v[2] as f64]
}

/// The brief's face-link criterion plus mesh sanity: every triangle maps to a
/// real face, every face has triangles, every vertex has a (u, v) and a unit
/// normal, triangles wind outward and agree with the kernel's normals, and
/// the mesh encloses the solid's volume.
fn assert_mesh(name: &str, built: &Built) {
    let mesh = &built.mesh;
    let vertices = mesh.positions.len();
    assert!(mesh.triangle_count() > 0, "{name}: empty mesh");
    assert_eq!(mesh.normals.len(), vertices, "{name}: normals");
    assert_eq!(mesh.surface_uvs.len(), vertices, "{name}: uvs");
    assert_eq!(
        mesh.indices.len(),
        mesh.triangle_count() * 3,
        "{name}: indices"
    );
    assert_eq!(mesh.faces.len(), built.brep_faces, "{name}: faces");
    assert_eq!(mesh.edges.len(), built.brep_edges, "{name}: edges");

    assert!(
        mesh.face_ids
            .iter()
            .all(|&f| (f as usize) < mesh.faces.len()),
        "{name}: triangle with an invalid face id"
    );
    let used: HashSet<u32> = mesh.face_ids.iter().copied().collect();
    assert_eq!(
        used.len(),
        mesh.faces.len(),
        "{name}: a face has no triangles"
    );
    assert!(
        mesh.indices.iter().all(|&i| (i as usize) < vertices),
        "{name}: index out of range"
    );
    assert!(
        mesh.surface_uvs.iter().flatten().all(|c| c.is_finite()),
        "{name}: a vertex has no (u, v)"
    );
    assert!(
        mesh.positions.iter().flatten().all(|c| c.is_finite()),
        "{name}: non-finite position"
    );
    for n in &mesh.normals {
        let len = dot(widen(*n), widen(*n)).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-3,
            "{name}: normal {n:?} isn't unit length"
        );
    }

    // Vertices are split per face: no vertex is used by two faces.
    let mut vertex_face = vec![u32::MAX; vertices];
    for (t, tri) in mesh.triangles().enumerate() {
        for i in tri {
            let owner = &mut vertex_face[i as usize];
            assert!(
                *owner == u32::MAX || *owner == mesh.face_ids[t],
                "{name}: vertex {i} shared by two faces"
            );
            *owner = mesh.face_ids[t];
        }
    }

    // Winding agrees with the kernel's normals, and the signed volume of the
    // closed mesh matches the solid's.
    let mut signed_volume = 0.0;
    let mut disagreeing = 0;
    for tri in mesh.triangles() {
        let [a, b, c] = tri.map(|i| mesh.positions[i as usize]);
        let geometric = cross(sub(b, a), sub(c, a));
        let kernel: [f64; 3] = tri.iter().fold([0.0; 3], |acc, &i| {
            let n = widen(mesh.normals[i as usize]);
            [acc[0] + n[0], acc[1] + n[1], acc[2] + n[2]]
        });
        if dot(geometric, kernel) <= 0.0 {
            disagreeing += 1;
        }
        signed_volume += dot(widen(a), cross(widen(b), widen(c))) / 6.0;
    }
    assert_eq!(
        disagreeing, 0,
        "{name}: {disagreeing} triangles wind against the kernel's normals"
    );
    let error = (signed_volume - built.volume).abs() / built.volume;
    println!(
        "{name}: {} triangles, {vertices} vertices, {} faces, {} edges; mesh volume {signed_volume:.1} ({:.3}% off)",
        mesh.triangle_count(),
        mesh.faces.len(),
        mesh.edges.len(),
        error * 100.0
    );
    assert!(
        error < 5e-3,
        "{name}: mesh encloses {signed_volume}, solid is {}",
        built.volume
    );

    for (i, edge) in mesh.edges.iter().enumerate() {
        assert!(edge.polyline.len() >= 2, "{name}: edge {i} has no polyline");
        assert_eq!(edge.use_count, 2, "{name}: edge {i} isn't manifold");
        assert!(
            edge.faces
                .iter()
                .all(|f| f.is_some_and(|f| (f as usize) < mesh.faces.len())),
            "{name}: edge {i} faces"
        );
    }
}

fn surface_kinds(mesh: &Mesh) -> Vec<&'static str> {
    let mut kinds: Vec<_> = mesh
        .faces
        .iter()
        .map(|f| match f.surface {
            Surface::Plane { .. } => "plane",
            Surface::Cylinder { .. } => "cylinder",
            Surface::Cone { .. } => "cone",
            Surface::Sphere { .. } => "sphere",
            Surface::Torus { .. } => "torus",
            Surface::Nurbs { .. } => "nurbs",
            Surface::Other { .. } => "other",
        })
        .collect();
    kinds.sort();
    kinds
}

#[test]
fn box_case() {
    let built = build(block);
    assert!(built.manifold);
    assert_volume("box", built.volume, BOX_VOLUME);
    assert_mesh("box", &built);

    let mesh = &built.mesh;
    let ids: HashSet<u32> = mesh.face_ids.iter().copied().collect();
    assert_eq!(ids.len(), 6, "the box mesh should carry exactly 6 face ids");
    assert_eq!(surface_kinds(mesh), ["plane"; 6]);
    assert_eq!(mesh.edges.len(), 12);
    // Hard edges: every vertex of a face carries that face's axis-aligned
    // normal.
    for (t, tri) in mesh.triangles().enumerate() {
        let first = mesh.normals[tri[0] as usize];
        for i in tri {
            assert_eq!(
                mesh.normals[i as usize], first,
                "face {} isn't flat-shaded",
                mesh.face_ids[t]
            );
        }
        assert_eq!(
            first.iter().filter(|c| c.abs() > 0.999).count(),
            1,
            "box normal {first:?} isn't axis-aligned"
        );
    }
}

#[test]
fn through_hole() {
    let r = 10.0;
    let built =
        build(move |k| block(k)?.difference(Brep::cylinder(k, [0.0, 0.0, -H / 2.0], r, H)?));
    assert!(built.manifold);
    assert_volume("hole", built.volume, BOX_VOLUME - PI * r * r * H);
    assert_mesh("hole", &built);
    assert!(surface_kinds(&built.mesh).contains(&"cylinder"));
}

#[test]
fn boss() {
    let (r, h) = (12.0, 20.0);
    let built = build(move |k| block(k)?.union(Brep::cylinder(k, [0.0, 0.0, H / 2.0], r, h)?));
    assert!(built.manifold);
    assert_volume("boss", built.volume, BOX_VOLUME + PI * r * r * h);
    assert_mesh("boss", &built);
    assert!(surface_kinds(&built.mesh).contains(&"cylinder"));
}

#[test]
fn rounded_block() {
    let built = build(|k| block(k)?.intersection(Brep::sphere(k, [0.0; 3], 55.0)?));
    assert!(built.manifold);
    println!("rounded: volume {:.3}", built.volume);
    assert!(
        built.volume > 0.0 && built.volume < BOX_VOLUME,
        "rounded block volume {}",
        built.volume
    );
    assert_mesh("rounded", &built);
    assert!(surface_kinds(&built.mesh).contains(&"sphere"));
}

#[test]
fn kernel_belongs_to_one_thread() {
    // Make sure the kernel thread has claimed the kernel, then try from here.
    on_kernel(|_| ());
    assert_eq!(Kernel::init().unwrap_err(), Error::WrongThread);
}

#[test]
fn rejects_bad_dimensions() {
    let err = on_kernel(|k| {
        Brep::cuboid(k, [0.0; 3], [1.0, f64::NAN, 1.0])
            .map(|_| ())
            .unwrap_err()
    });
    assert!(matches!(err, Error::InvalidArgument { .. }), "{err:?}");
}

const UP: [f64; 3] = [0.0, 0.0, 1.0];
const X: [f64; 3] = [1.0, 0.0, 0.0];
const Y: [f64; 3] = [0.0, 1.0, 0.0];
const Z: [f64; 3] = UP;

fn neg([x, y, z]: [f64; 3]) -> [f64; 3] {
    [-x, -y, -z]
}

/// The box's four vertical edges: between a ±X face and a ±Y face.
fn vertical_edges() -> EdgeQuery {
    let mut any = Vec::new();
    for x in [X, neg(X)] {
        for y in [Y, neg(Y)] {
            any.push(EdgeQuery::BetweenFaces { a: x, b: y });
        }
    }
    EdgeQuery::Or(any)
}

/// The rim of a hole of radius `r` along Z at (x, y), where it meets the top
/// face.
fn top_rim(x: f64, y: f64, r: f64) -> EdgeQuery {
    EdgeQuery::OnCircle {
        center: [x, y, H / 2.0],
        axis: UP,
        radius: r,
    }
}

fn hole(k: &Kernel, r: f64) -> Result<Brep, Error> {
    block(k)?.difference(Brep::cylinder(k, [0.0, 0.0, -H / 2.0], r, H)?)
}

/// A box with every edge and corner rounded by `r`: the inner box, slabs on
/// its faces, quarter cylinders on its edges and sphere octants at its
/// corners.
fn rounded_box_volume(r: f64) -> f64 {
    let (a, b, c) = (L - 2.0 * r, W - 2.0 * r, H - 2.0 * r);
    a * b * c
        + 2.0 * r * (a * b + b * c + c * a)
        + PI * r * r * (a + b + c)
        + 4.0 / 3.0 * PI * r.powi(3)
}

/// The material a 45° chamfer with legs `d` removes from the convex rim of a
/// hole of radius `r`: a triangle of area d²/2 swept around the axis at the
/// radius of its centroid, r + d/3 (Pappus).
fn rim_chamfer_volume(r: f64, d: f64) -> f64 {
    PI * d * d * (r + d / 3.0)
}

fn options() -> FilletOptions {
    FilletOptions::default()
}

#[test]
fn edge_queries_pick_by_geometry() {
    let tol = Tolerance::default();
    let counts = on_kernel(move |k| {
        let block = block(k).unwrap();
        let count = |q: EdgeQuery| block.count_edges(&q, &tol).unwrap();
        let top = EdgeQuery::MidpointOnPlane {
            point: [0.0, 0.0, H / 2.0],
            normal: UP,
        };
        let bottom = EdgeQuery::InPlane {
            point: [0.0, 0.0, -H / 2.0],
            normal: neg(UP),
        };
        let box_counts = [
            count(EdgeQuery::All),
            count(EdgeQuery::Sharp {
                min_angle_deg: 30.0,
            }),
            count(top.clone()),
            count(bottom.clone()),
            count(EdgeQuery::BetweenFaces { a: Z, b: X }),
            count(EdgeQuery::BetweenFaces {
                a: neg(Y),
                b: neg(Z),
            }),
            count(vertical_edges()),
            count(!top.clone().or(bottom.clone())),
            count(top.and(EdgeQuery::BetweenFaces { a: Z, b: neg(X) })),
        ];

        let part = hole(k, 10.0).unwrap();
        let edges = part.edges().unwrap();
        let lengths: Vec<f64> = (0..edges.len())
            .map(|edge| part.edge_length(edge).unwrap())
            .collect();
        let seams = edges
            .iter()
            .filter(|e| e.bend_deg().is_some_and(|b| b < 1.0))
            .count();
        let count = |q: EdgeQuery| part.count_edges(&q, &tol).unwrap();
        let hole_counts = [
            edges.len(),
            seams,
            count(EdgeQuery::Sharp {
                min_angle_deg: 30.0,
            }),
            count(top_rim(0.0, 0.0, 10.0)),
            count(top_rim(0.0, 0.0, 10.5)),
            count(EdgeQuery::OnCircle {
                center: [0.0, 0.0, -H / 2.0],
                axis: neg(UP),
                radius: 10.0,
            }),
        ];
        (box_counts, hole_counts, edges, lengths)
    });
    let (box_counts, [hole_edges, seams, sharp, rim, wrong_rim, bottom_rim], edges, lengths) =
        counts;
    println!(
        "box picks {box_counts:?}; hole: {hole_edges} edges, {seams} seams, {sharp} sharp, rim {rim}"
    );
    assert_eq!(box_counts, [12, 12, 4, 4, 1, 1, 4, 4, 1]);
    assert!(seams >= 1, "the hole's cylinder should have a seam");
    assert_eq!(sharp, hole_edges - seams, "seams aren't sharp");
    assert!(rim >= 1, "the top rim");
    assert_eq!(wrong_rim, 0, "a circle 0.5 mm too wide picks nothing");
    assert_eq!(bottom_rim, rim, "top and bottom rims split alike");
    let rim_length: f64 = edges
        .iter()
        .filter(|e| top_rim(0.0, 0.0, 10.0).matches(e, &tol))
        .map(|e| lengths[e.index])
        .sum();
    assert!(
        (rim_length - 20.0 * std::f64::consts::PI).abs() < 1e-6,
        "the rim is {rim_length} long"
    );
    for (edge, length) in edges.iter().zip(&lengths) {
        assert_eq!(edge.use_count, 2);
        assert!(*length > 0.0);
        for n in edge.normals.iter().flatten() {
            assert!(
                (dot(*n, *n).sqrt() - 1.0).abs() < 1e-9,
                "normal {n:?} isn't unit length"
            );
        }
    }
}

#[test]
fn filleted_box() {
    let r = 8.0;
    let built = build(move |k| {
        let block = block(k)?;
        assert_eq!(
            block.count_edges(
                &EdgeQuery::Sharp {
                    min_angle_deg: 30.0
                },
                &Tolerance::default()
            )?,
            12
        );
        block.fillet(
            &[EdgeBlend::round(
                EdgeQuery::Sharp {
                    min_angle_deg: 30.0,
                },
                r,
            )],
            &options(),
        )
    });
    assert!(built.manifold);
    assert_volume("filleted box", built.volume, rounded_box_volume(r));
    assert_mesh("filleted box", &built);
    println!("filleted box surfaces: {:?}", surface_kinds(&built.mesh));
    assert_eq!(
        built.brep_faces,
        6 + 12 + 8,
        "6 faces, 12 edge fillets, 8 corners"
    );
}

#[test]
fn filleted_box_across_radii() {
    for r in [0.5, 2.0, 5.0, 12.0, 19.0] {
        let built =
            build(move |k| block(k)?.fillet(&[EdgeBlend::round(EdgeQuery::All, r)], &options()));
        let name = format!("filleted box r={r}");
        assert!(built.manifold, "{name}");
        assert_volume(&name, built.volume, rounded_box_volume(r));
        assert_mesh(&name, &built);
    }
}

#[test]
fn chamfered_hole_rim() {
    // A linear section: a 45° chamfer with legs equal to the radius.
    let r = 10.0;
    for d in [0.5, 1.0, 2.0, 3.0, 4.0, 6.0] {
        let built = build(move |k| {
            hole(k, r)?.fillet(&[EdgeBlend::linear(top_rim(0.0, 0.0, r), d)], &options())
        });
        let name = format!("chamfered rim d={d}");
        assert!(built.manifold, "{name}");
        assert_volume(
            &name,
            built.volume,
            BOX_VOLUME - PI * r * r * H - rim_chamfer_volume(r, d),
        );
        assert_mesh(&name, &built);
        assert_eq!(
            built.brep_faces, 8,
            "{name}: the hole's 7 faces and the chamfer"
        );
    }
}

#[test]
fn chamfered_hole_off_centre() {
    let (x, y, r, d) = (25.0, -12.0, 6.0, 2.0);
    let built = build(move |k| {
        let part = block(k)?.difference(Brep::cylinder(k, [x, y, -H / 2.0], r, H)?)?;
        part.fillet(&[EdgeBlend::linear(top_rim(x, y, r), d)], &options())
    });
    assert!(built.manifold);
    let expected = BOX_VOLUME - PI * r * r * H - rim_chamfer_volume(r, d);
    assert_volume("off-centre chamfered rim", built.volume, expected);
    assert_mesh("off-centre chamfered rim", &built);
}

#[test]
fn rounded_hole_rim() {
    // A round removes the spandrel between the corner and the arc, area
    // (1 − π/4)·d², whose centroid sits (10 − 3π) / (12 − 3π)·d from each
    // face.
    let (r, d) = (10.0, 3.0);
    let built = build(move |k| {
        hole(k, r)?.fillet(&[EdgeBlend::round(top_rim(0.0, 0.0, r), d)], &options())
    });
    let spandrel = (1.0 - PI / 4.0) * d * d;
    let centroid = r + (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI) * d;
    assert!(built.manifold);
    assert_volume(
        "rounded rim",
        built.volume,
        BOX_VOLUME - PI * r * r * H - 2.0 * PI * centroid * spandrel,
    );
    assert_mesh("rounded rim", &built);
}

#[test]
fn constant_width_chamfer_on_a_circle() {
    // With SmApiChamferFillet's global merge on, legs of about 2.5 to 4 on
    // the 10 mm rim (and 3 on the 20 mm one) came back as open shells. The
    // shim leaves it off; every size should now be exact.
    for (r, legs) in [
        (10.0, &[0.5, 1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0][..]),
        (20.0, &[3.0, 6.0][..]),
    ] {
        for &d in legs {
            let built = build(move |k| {
                hole(k, r)?.chamfer(
                    &top_rim(0.0, 0.0, r),
                    d * 2f64.sqrt(),
                    &Tolerance::default(),
                )
            });
            let name = format!("constant-width rim r={r} d={d}");
            assert!(built.manifold, "{name}");
            assert_volume(
                &name,
                built.volume,
                BOX_VOLUME - PI * r * r * H - rim_chamfer_volume(r, d),
            );
            assert_mesh(&name, &built);
            assert_eq!(built.brep_faces, 8, "{name}");
        }
    }
}

#[test]
fn chamfer_too_big_to_fit_fails_cleanly() {
    // A chamfer with legs of 3 on the top rim of a 2 mm boss would cut past
    // the boss's whole side. The kernel refuses it; the kernel keeps working.
    let error = on_kernel(|k| {
        let part = block(k)
            .unwrap()
            .union(Brep::cylinder(k, [0.0, 0.0, H / 2.0], 12.0, 2.0).unwrap())
            .unwrap();
        let rim = EdgeQuery::OnCircle {
            center: [0.0, 0.0, H / 2.0 + 2.0],
            axis: UP,
            radius: 12.0,
        };
        part.chamfer(&rim, 3.0 * 2f64.sqrt(), &Tolerance::default())
            .map(|_| ())
            .unwrap_err()
    });
    println!("chamfer past a 2 mm boss: {error}");
    assert!(matches!(error, Error::Kernel { .. }), "{error:?}");
    let built = build(block);
    assert_volume("box after a refused chamfer", built.volume, BOX_VOLUME);
}

#[test]
fn per_edge_radii_and_sections() {
    // Vertical edges only, so no corner joins three fillets and the volume is
    // exact: two rounds of 5 and two flat cuts set back 4.
    let (round, flat) = (5.0, 4.0);
    let plus_x = EdgeQuery::Or(vec![
        EdgeQuery::BetweenFaces { a: X, b: Y },
        EdgeQuery::BetweenFaces { a: X, b: neg(Y) },
    ]);
    let built = build(move |k| {
        block(k)?.fillet(
            &[
                EdgeBlend::round(plus_x, round),
                EdgeBlend::linear(vertical_edges(), flat),
            ],
            &options(),
        )
    });
    let expected =
        BOX_VOLUME - 2.0 * (1.0 - PI / 4.0) * round * round * H - 2.0 * flat * flat / 2.0 * H;
    assert!(built.manifold);
    assert_volume("mixed vertical fillets", built.volume, expected);
    assert_mesh("mixed vertical fillets", &built);
    assert_eq!(built.brep_faces, 10);
}

#[test]
fn blend_sections() {
    let r = 6.0;
    let flat = BOX_VOLUME - 4.0 * r * r / 2.0 * H;
    let round = BOX_VOLUME - 4.0 * (1.0 - PI / 4.0) * r * r * H;
    for continuity in [Continuity::G1, Continuity::G2, Continuity::G3] {
        let built = build(move |k| {
            let options = FilletOptions {
                continuity,
                ..options()
            };
            block(k)?.fillet(&[EdgeBlend::blend(vertical_edges(), r)], &options)
        });
        let name = format!("{continuity:?} blend");
        assert!(built.manifold, "{name}");
        println!(
            "{name}: volume {:.1}; flat {flat:.1}, round {round:.1}",
            built.volume
        );
        // A blend between the same lines bulges past the flat cut.
        assert!(
            built.volume > flat && built.volume < BOX_VOLUME,
            "{name}: {}",
            built.volume
        );
        assert_mesh(&name, &built);
    }
}

#[test]
fn chamfered_box_edges() {
    let w = 4.0;
    // Vertical edges: each loses a right triangle with legs w/√2.
    let built = build(move |k| block(k)?.chamfer(&vertical_edges(), w, &Tolerance::default()));
    assert!(built.manifold);
    assert_volume(
        "chamfered vertical edges",
        built.volume,
        BOX_VOLUME - 4.0 * w * w / 4.0 * H,
    );
    assert_mesh("chamfered vertical edges", &built);

    // Every edge, with beveled corners. Each edge's prism is w²/4 in section;
    // summed over full edge lengths they overcount the corners.
    let built = build(move |k| block(k)?.chamfer(&EdgeQuery::All, w, &Tolerance::default()));
    assert!(built.manifold);
    let prisms = w * w / 4.0 * 4.0 * (L + W + H);
    println!(
        "chamfered box: volume {:.1}, faces {}",
        built.volume, built.brep_faces
    );
    assert!(
        built.volume > BOX_VOLUME - prisms && built.volume < BOX_VOLUME - 4.0 * w * w / 4.0 * H
    );
    assert_mesh("chamfered box", &built);
}

#[test]
fn rejects_bad_fillets() {
    let results = on_kernel(|k| {
        let tol = Tolerance::default();
        let nowhere = EdgeQuery::MidpointOnPlane {
            point: [0.0, 0.0, 100.0],
            normal: UP,
        };
        let no_normal = EdgeQuery::BetweenFaces { a: [0.0; 3], b: X };
        [
            block(k)
                .unwrap()
                .fillet(&[EdgeBlend::round(nowhere.clone(), 2.0)], &options()),
            block(k)
                .unwrap()
                .fillet(&[EdgeBlend::round(EdgeQuery::All, f64::NAN)], &options()),
            block(k)
                .unwrap()
                .fillet(&[EdgeBlend::round(no_normal, 2.0)], &options()),
            block(k).unwrap().fillet(&[], &options()),
            block(k).unwrap().chamfer(&nowhere, 2.0, &tol),
            block(k).unwrap().chamfer(&EdgeQuery::All, -1.0, &tol),
        ]
        .map(|r| r.map(|_| ()).unwrap_err())
    });
    for error in results {
        assert!(matches!(error, Error::InvalidArgument { .. }), "{error:?}");
    }
}

#[test]
fn oversized_fillet_fails_cleanly() {
    // Radius 25 can't fit the 40 mm height: the kernel either refuses or
    // returns something that isn't a valid solid. Either way, no crash, and
    // the kernel keeps working.
    let outcome = on_kernel(|k| {
        match block(k)
            .unwrap()
            .fillet(&[EdgeBlend::round(EdgeQuery::All, 25.0)], &options())
        {
            Err(error) => format!("refused: {error}"),
            Ok(brep) => format!("returned; manifold {:?}", brep.is_manifold()),
        }
    });
    println!("radius 25 on a 40 mm box: {outcome}");
    let built = build(block);
    assert_volume("box after a failed fillet", built.volume, BOX_VOLUME);
}

// Offset faces. Each case moves the top face of a part and compares the
// result with the same part built directly at the new height (the reference)
// and with its analytic volume.

/// Volumes for the offset oracles: tighter than the kernel's default, since
/// the criterion is 1e-6.
const PRECISE: f64 = 1e-8;

fn box_of_height(k: &Kernel, h: f64) -> Result<Brep, Error> {
    Brep::cuboid(k, [-L / 2.0, -W / 2.0, -H / 2.0], [L, W, h])
}

/// The top face of a part `h` tall standing on the box's base.
fn top_at(h: f64) -> FaceQuery {
    FaceQuery::OnPlane {
        point: [0.0, 0.0, -H / 2.0 + h],
        normal: UP,
    }
}

/// The four edges around the top of a part `h` tall.
fn top_edges_at(h: f64) -> EdgeQuery {
    EdgeQuery::MidpointOnPlane {
        point: [0.0, 0.0, -H / 2.0 + h],
        normal: UP,
    }
}

/// A box `h` tall with 45° chamfers of leg `c` on its four top edges.
fn chamfered_box(k: &Kernel, h: f64, c: f64) -> Result<Brep, Error> {
    box_of_height(k, h)?.chamfer(&top_edges_at(h), c * 2f64.sqrt(), &Tolerance::default())
}

/// The chamfered box's volume: each top edge loses a prism of section c²/2,
/// and adjacent prisms overlap by c³/3 at each mitred corner.
fn chamfered_box_volume(h: f64, c: f64) -> f64 {
    L * W * h - c * c * (L + W) + 4.0 * c.powi(3) / 3.0
}

/// A box `h` tall with rounds of radius `r` on its four top edges.
fn top_rounded_box(k: &Kernel, h: f64, r: f64) -> Result<Brep, Error> {
    box_of_height(k, h)?.fillet(&[EdgeBlend::round(top_edges_at(h), r)], &options())
}

/// The filleted box's volume: each top edge loses a spandrel of (1 − π/4)r²,
/// and adjacent spandrels overlap by (5/3 − π/2)r³ at each mitred corner.
fn top_rounded_box_volume(h: f64, r: f64) -> f64 {
    L * W * h - (1.0 - PI / 4.0) * r * r * 2.0 * (L + W) + 4.0 * r.powi(3) * (5.0 / 3.0 - PI / 2.0)
}

/// A frustum on the box's base, `h` tall, its walls leaning in by
/// `draft_deg`.
fn drafted_block(k: &Kernel, h: f64, draft_deg: f64) -> Result<Brep, Error> {
    let (x, y, z) = (L / 2.0, W / 2.0, -H / 2.0);
    let base = [[-x, -y, z], [x, -y, z], [x, y, z], [-x, y, z]];
    Brep::tapered_prism(k, &base, h, draft_deg)
}

fn drafted_block_volume(h: f64, draft_deg: f64) -> f64 {
    let t = draft_deg.to_radians().tan();
    L * W * h - (L + W) * t * h * h + 4.0 / 3.0 * t * t * h.powi(3)
}

/// A box `h` tall with a hole of radius `r` through it at the origin.
fn holed_box(k: &Kernel, h: f64, r: f64) -> Result<Brep, Error> {
    box_of_height(k, h)?.difference(Brep::cylinder(k, [0.0, 0.0, -H / 2.0], r, h)?)
}

struct Offset {
    built: Built,
    precise_volume: f64,
}

fn build_offset(make: impl FnOnce(&Kernel) -> Result<Brep, Error> + Send + 'static) -> Offset {
    on_kernel(move |k| {
        let brep = make(k).expect("building the case");
        Offset {
            precise_volume: brep.volume(PRECISE).unwrap(),
            built: Built {
                manifold: brep.is_manifold().unwrap(),
                volume: brep.volume(RELATIVE_ACCURACY).unwrap(),
                brep_faces: brep.face_count(),
                brep_edges: brep.edge_count(),
                mesh: brep.tessellate(&TessellationParams::default()).unwrap(),
            },
        }
    })
}

/// Offsets the top face of `part` (`H` tall) by `d` and holds it to the
/// brief's criteria against the part built `H + d` tall: see
/// [`assert_offset_case`].
fn assert_offset(name: &str, d: f64, part: fn(&Kernel, f64) -> Result<Brep, Error>, analytic: f64) {
    assert_offset_case(
        &format!("{name} d={d:+}"),
        move |k| part(k, H)?.offset_where(&top_at(H), d),
        move |k| part(k, H + d),
        Some(analytic),
    );
}

/// Holds an offset to the brief's criteria against `reference`, the same part
/// built directly at its new size, and its analytic volume when there is one:
/// manifold, the same face and edge counts, volume within 1e-6, and the mesh
/// checks.
fn assert_offset_case(
    name: &str,
    moved: impl FnOnce(&Kernel) -> Result<Brep, Error> + Send + 'static,
    reference: impl FnOnce(&Kernel) -> Result<Brep, Error> + Send + 'static,
    analytic: Option<f64>,
) -> Built {
    let moved = build_offset(moved);
    let reference = build_offset(reference);
    let vs_reference = (moved.precise_volume - reference.precise_volume) / reference.precise_volume;
    let vs_analytic = analytic.map(|a| (moved.precise_volume - a) / a);
    println!(
        "| {name} | {}/{} faces, {}/{} edges | {:.4} | {:.4} | {} | {vs_reference:+.1e} | {} |",
        moved.built.brep_faces,
        reference.built.brep_faces,
        moved.built.brep_edges,
        reference.built.brep_edges,
        moved.precise_volume,
        reference.precise_volume,
        analytic.map_or("-".into(), |a| format!("{a:.4}")),
        vs_analytic.map_or("-".into(), |e| format!("{e:+.1e}")),
    );
    assert!(moved.built.manifold, "{name}: not manifold");
    assert!(reference.built.manifold, "{name}: reference not manifold");
    assert_eq!(
        moved.built.brep_faces, reference.built.brep_faces,
        "{name}: faces"
    );
    assert_eq!(
        moved.built.brep_edges, reference.built.brep_edges,
        "{name}: edges"
    );
    assert!(
        vs_reference.abs() < 1e-6,
        "{name}: {vs_reference:e} off the reference"
    );
    if let Some(error) = vs_analytic {
        assert!(error.abs() < 1e-6, "{name}: {error:e} off analytic");
    }
    assert_mesh(name, &moved.built);
    moved.built
}

fn status_name(error: &Error) -> Option<&'static str> {
    match error {
        Error::Kernel { status, .. } => status.name(),
        _ => None,
    }
}

/// Offsets the faces a query picks: [`Brep::find_faces`] then
/// [`Brep::offset_faces`].
trait OffsetWhere: Sized {
    fn offset_where(self, query: &FaceQuery, distance: f64) -> Result<Brep, Error>;
}

impl OffsetWhere for Brep {
    fn offset_where(self, query: &FaceQuery, distance: f64) -> Result<Brep, Error> {
        let faces = self.find_faces(query, &Tolerance::default())?;
        self.offset_faces(&faces, distance)
    }
}

#[test]
fn face_queries_pick_by_geometry() {
    let tol = Tolerance::default();
    let counts = on_kernel(move |k| {
        let block = block(k).unwrap();
        let count = |b: &Brep, q: FaceQuery| b.find_faces(&q, &tol).unwrap().len();
        let rounded = top_rounded_box(k, H, 8.0).unwrap();
        let drafted = drafted_block(k, H, 5.0).unwrap();
        [
            count(&block, FaceQuery::All),
            count(&block, top_at(H)),
            count(&block, top_at(H + 1.0)),
            count(
                &block,
                FaceQuery::OnPlane {
                    point: [0.0, 0.0, H / 2.0],
                    normal: neg(UP),
                },
            ),
            count(&block, FaceQuery::Facing { direction: X }),
            count(&block, FaceQuery::Kind(SurfaceKind::Plane)),
            count(
                &block,
                !FaceQuery::Facing { direction: UP }.or(FaceQuery::Facing { direction: neg(UP) }),
            ),
            count(&rounded, FaceQuery::Kind(SurfaceKind::Cylinder)),
            count(&rounded, top_at(H)),
            // A drafted block's walls are planar NURBS surfaces: they have
            // normals, so Facing finds them.
            count(
                &drafted,
                FaceQuery::Facing {
                    direction: [5f64.to_radians().cos(), 0.0, 5f64.to_radians().sin()],
                },
            ),
        ]
    });
    assert_eq!(counts, [6, 1, 0, 0, 1, 6, 4, 4, 1, 1]);
    let error = on_kernel(move |k| {
        block(k)
            .unwrap()
            .find_faces(
                &FaceQuery::Facing {
                    direction: [0.0; 3],
                },
                &tol,
            )
            .unwrap_err()
    });
    assert!(matches!(error, Error::InvalidArgument { .. }), "{error:?}");
}

#[test]
fn offset_face_plain_walls() {
    for d in [10.0, -10.0] {
        assert_offset("plain box", d, box_of_height, L * W * (H + d));
    }
}

#[test]
fn offset_face_grows_and_shrinks_chamfers() {
    // The chamfer planes extend along their own geometry: legs of 2 become
    // 2 + d. The reference is chamfered with those legs.
    for d in [10.0, -1.0] {
        assert_offset(
            "chamfered box",
            d,
            |k, h| chamfered_box(k, h, 2.0 + (h - H)),
            chamfered_box_volume(H + d, 2.0 + d),
        );
    }
}

#[test]
fn offset_face_carries_fillets() {
    // The rounds ride on the face; the walls they lean on stretch to follow.
    for d in [10.0, -10.0] {
        assert_offset(
            "filleted box",
            d,
            |k, h| top_rounded_box(k, h, 8.0),
            top_rounded_box_volume(H + d, 8.0),
        );
    }
}

#[test]
fn offset_face_carries_fillets_at_every_distance() {
    // The offset executive's iteration order follows pointer order, so a
    // fragile result shows up as a scattering of failures across distances.
    let results = on_kernel(|k| {
        (-31..=31)
            .filter(|&i| i != 0)
            .map(|i| {
                let d = i as f64;
                let moved = top_rounded_box(k, H, 8.0)
                    .and_then(|part| part.offset_where(&top_at(H), d))
                    .and_then(|b| Ok((b.is_manifold()?, b.face_count(), b.volume(PRECISE)?)));
                (d, moved)
            })
            .collect::<Vec<_>>()
    });
    for (d, moved) in results {
        let (manifold, faces, volume) = moved.unwrap_or_else(|e| panic!("d={d}: {e}"));
        let error =
            (volume - top_rounded_box_volume(H + d, 8.0)) / top_rounded_box_volume(H + d, 8.0);
        assert!(
            manifold && faces == 10 && error.abs() < 1e-6,
            "d={d}: manifold {manifold}, {faces} faces, {error:e}"
        );
    }
}

#[test]
fn offset_face_extends_drafted_walls() {
    for d in [10.0, -10.0] {
        assert_offset(
            "drafted block 5°",
            d,
            |k, h| drafted_block(k, h, 5.0),
            drafted_block_volume(H + d, 5.0),
        );
    }
}

#[test]
fn offset_face_around_a_hole() {
    // The hole's wall runs along the move, so it simply extends; a round on
    // its rim rides along like the box's fillets.
    for d in [10.0, -10.0] {
        let r = 10.0;
        assert_offset(
            "holed box",
            d,
            |k, h| holed_box(k, h, 10.0),
            (L * W - PI * r * r) * (H + d),
        );
        assert_offset(
            "holed box, rounded rim",
            d,
            |k, h| {
                let rim = EdgeQuery::OnCircle {
                    center: [0.0, 0.0, -H / 2.0 + h],
                    axis: UP,
                    radius: 10.0,
                };
                holed_box(k, h, 10.0)?.fillet(&[EdgeBlend::round(rim, 3.0)], &options())
            },
            {
                // As in rounded_hole_rim.
                let spandrel = (1.0 - PI / 4.0) * 9.0;
                let centroid = r + (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI) * 3.0;
                (L * W - PI * r * r) * (H + d) - 2.0 * PI * centroid * spandrel
            },
        );
    }
}

#[test]
fn offset_face_fails_cleanly_when_a_face_would_vanish() {
    // Each move would shrink a face to nothing that the moved face can't
    // simply pass: the walls under the fillets from d = -32, the box's walls
    // at d = -40, the drafted block's top past its apex (about d = +303).
    // The kernel hands back an open or inside-out solid; the shim refuses
    // it. (A chamfer along the moved face's edge it can pass: see
    // a_pushed_face_meets_the_face_beyond_one_it_swallows.)
    let outcomes = on_kernel(move |k| {
        let attempts: Vec<(&str, f64, Result<Brep, Error>)> = vec![
            ("filleted box", -32.0, top_rounded_box(k, H, 8.0)),
            ("filleted box", -35.0, top_rounded_box(k, H, 8.0)),
            ("plain box", -40.0, box_of_height(k, H)),
            ("plain box", -45.0, box_of_height(k, H)),
            ("drafted block", 340.0, drafted_block(k, H, 5.0)),
        ];
        attempts
            .into_iter()
            .map(|(name, d, part)| {
                let outcome = part
                    .unwrap()
                    .offset_where(&top_at(H), d)
                    .map(|b| b.face_count());
                (name, d, outcome)
            })
            .collect::<Vec<_>>()
    });
    for (name, d, outcome) in outcomes {
        println!("{name} d={d:+}: {outcome:?}");
        let error = outcome.expect_err("a vanishing face must not come back as a solid");
        assert!(
            matches!(
                status_name(&error),
                Some("SHIM_ERR_BROKE_SOLID" | "SHIM_ERR_CHANGED_TOPOLOGY")
            ),
            "{name} d={d}: {error}"
        );
    }
    // The kernel still works, and a move just short of vanishing is exact.
    assert_offset(
        "chamfered box",
        -1.99,
        |k, h| chamfered_box(k, h, 2.0 + (h - H)),
        chamfered_box_volume(H - 1.99, 0.01),
    );
    assert_offset(
        "drafted block 5°",
        300.0,
        |k, h| drafted_block(k, h, 5.0),
        drafted_block_volume(H + 300.0, 5.0),
    );
}

#[test]
fn offset_face_rejects_what_it_cannot_move() {
    let errors = on_kernel(move |k| {
        let rim = EdgeQuery::OnCircle {
            center: [0.0, 0.0, H / 2.0],
            axis: UP,
            radius: 10.0,
        };
        let rounded_rim = holed_box(k, H, 10.0)
            .unwrap()
            .fillet(&[EdgeBlend::round(rim, 3.0)], &options())
            .unwrap();
        let top = block(k)
            .unwrap()
            .find_faces(&top_at(H), &Tolerance::default())
            .unwrap()[0];
        [
            block(k).unwrap().offset_faces(&[], 5.0),
            block(k).unwrap().offset_where(&top_at(H + 1.0), 5.0),
            block(k).unwrap().offset_faces(&[top], f64::NAN),
            block(k).unwrap().offset_faces(&[6], 5.0),
            block(k).unwrap().offset_faces(&[top, top], 5.0),
            // A round on a picked curved face: its new shape isn't a
            // translation of the old one.
            rounded_rim.offset_where(&FaceQuery::Kind(SurfaceKind::Cylinder), 1.0),
        ]
        .map(|r| r.map(|_| ()).unwrap_err())
    });
    let [none, none_picked, nan, out_of_range, repeated, rounded_rim] = errors;
    for error in [none, none_picked, nan] {
        assert!(matches!(error, Error::InvalidArgument { .. }), "{error:?}");
    }
    for error in [out_of_range, repeated] {
        assert_eq!(status_name(&error), Some("SHIM_ERR_BAD_FACE"), "{error}");
    }
    assert_eq!(
        status_name(&rounded_rim),
        Some("SHIM_ERR_UNSUPPORTED_FACE"),
        "{rounded_rim}"
    );

    // Zero is a no-op.
    let built = build(|k| block(k)?.offset_where(&top_at(H), 0.0));
    assert_volume("offset by zero", built.volume, BOX_VOLUME);
}

/// A box from its minimum and maximum corners.
fn span(k: &Kernel, min: [f64; 3], max: [f64; 3]) -> Result<Brep, Error> {
    Brep::cuboid(k, min, [max[0] - min[0], max[1] - min[1], max[2] - min[2]])
}

/// The face of a box whose outward normal is `normal`.
fn side(normal: [f64; 3]) -> FaceQuery {
    FaceQuery::Facing { direction: normal }.and(FaceQuery::Kind(SurfaceKind::Plane))
}

/// The box's walls: its planar faces facing sideways.
fn walls() -> FaceQuery {
    side(X).or(side(neg(X))).or(side(Y)).or(side(neg(Y)))
}

/// A box spanning `min` to `max` with rounds of radius `r` on its top edges.
fn rounded_span(k: &Kernel, min: [f64; 3], max: [f64; 3], r: f64) -> Result<Brep, Error> {
    let top = EdgeQuery::MidpointOnPlane {
        point: [0.0, 0.0, max[2]],
        normal: UP,
    };
    span(k, min, max)?.fillet(&[EdgeBlend::round(top, r)], &options())
}

const LO: [f64; 3] = [-L / 2.0, -W / 2.0, -H / 2.0];

#[test]
fn offset_faces_moves_several_planes() {
    assert_offset_case(
        "plain box, top and +X d=+10",
        move |k| block(k)?.offset_where(&top_at(H).or(side(X)), 10.0),
        |k| span(k, LO, [60.0, 30.0, 30.0]),
        Some(110.0 * 60.0 * 50.0),
    );
    assert_offset_case(
        "plain box, every face d=-5",
        move |k| block(k)?.offset_where(&FaceQuery::All, -5.0),
        |k| span(k, [-45.0, -25.0, -15.0], [45.0, 25.0, 15.0]),
        Some(90.0 * 50.0 * 30.0),
    );
    assert_offset_case(
        "chamfered box, top and bottom d=+10",
        move |k| chamfered_box(k, H, 2.0)?.offset_where(&top_at(H).or(side(neg(UP))), 10.0),
        |k| {
            let top = EdgeQuery::MidpointOnPlane {
                point: [0.0, 0.0, 30.0],
                normal: UP,
            };
            span(k, [-50.0, -30.0, -30.0], [50.0, 30.0, 30.0])?.chamfer(
                &top,
                12.0 * 2f64.sqrt(),
                &Tolerance::default(),
            )
        },
        Some(chamfered_box_volume(H + 20.0, 12.0)),
    );
    assert_offset_case(
        "filleted box, top and bottom d=+10",
        move |k| top_rounded_box(k, H, 8.0)?.offset_where(&top_at(H).or(side(neg(UP))), 10.0),
        |k| rounded_span(k, [-50.0, -30.0, -30.0], [50.0, 30.0, 30.0], 8.0),
        Some(top_rounded_box_volume(H + 20.0, 8.0)),
    );
}

/// The volume of a box `l` × `w` × `h` with rounds of radius `r` on its top
/// edges, as top_rounded_box_volume.
fn rounded_span_volume(l: f64, w: f64, h: f64, r: f64) -> f64 {
    l * w * h - (1.0 - PI / 4.0) * r * r * 2.0 * (l + w) + 4.0 * r.powi(3) * (5.0 / 3.0 - PI / 2.0)
}

#[test]
fn offset_faces_pushes_filleted_walls() {
    // A wall under a fillet: the fillet slides out with it and the top and the
    // other fillets stretch. Two picked faces: the fillet between them moves
    // with both.
    for d in [10.0, -5.0] {
        assert_offset_case(
            &format!("filleted box, +X wall d={d:+}"),
            move |k| top_rounded_box(k, H, 8.0)?.offset_where(&side(X), d),
            move |k| rounded_span(k, LO, [50.0 + d, 30.0, 20.0], 8.0),
            Some(rounded_span_volume(L + d, W, H, 8.0)),
        );
        assert_offset_case(
            &format!("filleted box, top and +X wall d={d:+}"),
            move |k| top_rounded_box(k, H, 8.0)?.offset_where(&top_at(H).or(side(X)), d),
            move |k| rounded_span(k, LO, [50.0 + d, 30.0, 20.0 + d], 8.0),
            Some(rounded_span_volume(L + d, W, H + d, 8.0)),
        );
    }
    assert_offset_case(
        "filleted box, top and every wall d=+5",
        move |k| top_rounded_box(k, H, 8.0)?.offset_where(&top_at(H).or(walls()), 5.0),
        |k| rounded_span(k, [-55.0, -35.0, -20.0], [55.0, 35.0, 25.0], 8.0),
        Some(rounded_span_volume(L + 10.0, W + 10.0, H + 5.0, 8.0)),
    );
    // One offset after another, as a feature tree applies them: the fillets
    // the first offset stretched are still rounds the next one carries.
    let fillets = on_kernel(|k| {
        let part = top_rounded_box(k, H, 8.0)?.offset_where(&side(X), 10.0)?;
        part.find_faces(
            &FaceQuery::Kind(SurfaceKind::Cylinder),
            &Tolerance::default(),
        )
    })
    .unwrap();
    assert_eq!(fillets.len(), 4, "the stretched fillets stay cylinders");
    assert_offset_case(
        "filleted box, +X wall d=+10, then top d=+5, then -Y wall d=+4",
        move |k| {
            top_rounded_box(k, H, 8.0)?
                .offset_where(&side(X), 10.0)?
                .offset_where(&top_at(H), 5.0)?
                .offset_where(&side(neg(Y)), 4.0)
        },
        |k| rounded_span(k, [-50.0, -34.0, -20.0], [60.0, 30.0, 25.0], 8.0),
        Some(rounded_span_volume(L + 10.0, W + 4.0, H + 5.0, 8.0)),
    );
    for d in [5.0, -5.0] {
        let all = |k: &Kernel, h: f64| {
            box_of_height(k, h)?.fillet(&[EdgeBlend::round(EdgeQuery::All, 8.0)], &options())
        };
        assert_offset_case(
            &format!("fully rounded box, top d={d:+}"),
            move |k| all(k, H)?.offset_where(&top_at(H), d),
            move |k| all(k, H + d),
            None,
        );
    }
}

#[test]
fn offset_face_moves_a_drafted_wall() {
    // A drafted wall itself: its base edge moves out by d / cos(5°).
    let wall = FaceQuery::Facing {
        direction: [5f64.to_radians().cos(), 0.0, 5f64.to_radians().sin()],
    };
    assert_offset_case(
        "drafted block, +X wall d=+5",
        move |k| drafted_block(k, H, 5.0)?.offset_where(&wall, 5.0),
        |k| {
            let (x, y, z) = (L / 2.0 + 5.0 / 5f64.to_radians().cos(), W / 2.0, -H / 2.0);
            Brep::tapered_prism(
                k,
                &[[-L / 2.0, -y, z], [x, -y, z], [x, y, z], [-L / 2.0, y, z]],
                H,
                5.0,
            )
        },
        None,
    );
}

#[test]
fn offset_faces_on_curved_faces() {
    // A curved face takes its offset surface: holes and bosses change radius,
    // a cone slides along its axis, a sphere or torus changes radius.
    let cylinder = FaceQuery::Kind(SurfaceKind::Cylinder);
    for (d, r) in [(2.0, 8.0), (-2.0, 12.0)] {
        let c = cylinder.clone();
        assert_offset_case(
            &format!("hole wall d={d:+}"),
            move |k| holed_box(k, H, 10.0)?.offset_where(&c, d),
            move |k| holed_box(k, H, r),
            Some((L * W - PI * r * r) * H),
        );
    }
    let c = cylinder.clone();
    assert_offset_case(
        "hole wall and top d=+2",
        move |k| holed_box(k, H, 10.0)?.offset_where(&c.or(top_at(H)), 2.0),
        |k| holed_box(k, H + 2.0, 8.0),
        Some((L * W - PI * 64.0) * (H + 2.0)),
    );
    let boss =
        |k: &Kernel, r: f64, h: f64| block(k)?.union(Brep::cylinder(k, [0.0, 0.0, H / 2.0], r, h)?);
    for (d, r) in [(2.0, 14.0), (-2.0, 10.0)] {
        let c = cylinder.clone();
        assert_offset_case(
            &format!("boss side d={d:+}"),
            move |k| boss(k, 12.0, 20.0)?.offset_where(&c, d),
            move |k| boss(k, r, 20.0),
            Some(BOX_VOLUME + PI * r * r * 20.0),
        );
    }
    assert_offset_case(
        "boss top d=+5",
        move |k| boss(k, 12.0, 20.0)?.offset_where(&top_at(H + 20.0), 5.0),
        move |k| boss(k, 12.0, 25.0),
        Some(BOX_VOLUME + PI * 144.0 * 25.0),
    );
    // A frustum, radius 20 to 10 over 30: its side moves out d / cos α, and
    // stays a cone (the offset is the cone slid along its axis).
    for d in [2.0, -2.0] {
        let grow = d / (10.0f64 / 30.0).atan().cos();
        let (r1, r2) = (20.0 + grow, 10.0 + grow);
        let built = assert_offset_case(
            &format!("frustum side d={d:+}"),
            move |k| {
                Brep::cone(k, [0.0; 3], 20.0, 10.0, 30.0)?
                    .offset_where(&FaceQuery::Kind(SurfaceKind::Cone), d)
            },
            move |k| Brep::cone(k, [0.0; 3], r1, r2, 30.0),
            Some(PI * 30.0 / 3.0 * (r1 * r1 + r1 * r2 + r2 * r2)),
        );
        assert!(
            surface_kinds(&built.mesh).contains(&"cone"),
            "frustum d={d}: {:?}",
            surface_kinds(&built.mesh)
        );
    }
    // A 45° countersink widened by 1: the cone slides down by √2.
    let sink = |k: &Kernel, apex: f64| {
        block(k)?.difference(Brep::cone(
            k,
            [0.0, 0.0, apex],
            0.0,
            H / 2.0 + 1.0 - apex,
            H / 2.0 + 1.0 - apex,
        )?)
    };
    assert_offset_case(
        "countersink d=-1",
        move |k| sink(k, 10.0)?.offset_where(&FaceQuery::Kind(SurfaceKind::Cone), -1.0),
        move |k| sink(k, 10.0 - 2f64.sqrt()),
        None,
    );
    let ball = |k: &Kernel, r: f64, h: f64| {
        box_of_height(k, h)?.intersection(Brep::sphere(k, [0.0; 3], r)?)
    };
    let sphere = FaceQuery::Kind(SurfaceKind::Sphere);
    assert_offset_case(
        "box cut by a sphere, sphere faces d=+2",
        move |k| ball(k, 55.0, H)?.offset_where(&sphere, 2.0),
        move |k| ball(k, 57.0, H),
        None,
    );
    // A 10 mm tube: thinner ones mesh coarser than assert_mesh allows at the
    // default 15°, offset or not.
    for (d, r) in [(1.0, 11.0), (-1.0, 9.0)] {
        assert_offset_case(
            &format!("torus d={d:+}"),
            move |k| Brep::torus(k, [0.0; 3], 30.0, 10.0)?.offset_where(&FaceQuery::All, d),
            move |k| Brep::torus(k, [0.0; 3], 30.0, r),
            Some(2.0 * PI * PI * 30.0 * r * r),
        );
    }
}

#[test]
fn offsets_come_out_the_same_every_time() {
    // The walls' edges with the vertical fillets stretch when the top rises.
    // Left to the executive, they varied in the last digits with the heap's
    // addresses; the same part must rebuild to the same bits.
    let runs = on_kernel(|k| {
        let part = box_of_height(k, H)?.fillet(
            &[EdgeBlend::round(
                EdgeQuery::ParallelTo { direction: UP },
                8.0,
            )],
            &options(),
        )?;
        let lowest = |b: &Brep| -> Result<f64, Error> {
            Ok(b.edges()?
                .iter()
                .flat_map(|e| e.samples)
                .map(|p| p[2])
                .fold(f64::INFINITY, f64::min))
        };
        let before = lowest(&part)?;
        let mut kept = Vec::new();
        let mut runs = Vec::new();
        for _ in 0..8 {
            let moved = part.try_clone()?.offset_where(&top_at(H), 15.0)?;
            let edges: Vec<_> = moved.edges()?.iter().map(|e| e.samples).collect();
            runs.push((edges, lowest(&moved)? == before));
            kept.push(moved);
        }
        Ok::<_, Error>(runs)
    })
    .unwrap();
    for run in &runs[1..] {
        assert!(run.0 == runs[0].0, "an offset came out different");
    }
    assert!(
        runs[0].1,
        "the walls' bottom ends stay exactly where they were"
    );
}

#[test]
fn offsets_mate_edges_in_the_same_order_every_time() {
    // The executive mated edges in the order of a map keyed by their
    // addresses, which also picked the side each was built from: these
    // parts came out different in the last digits from run to run.
    type Make = fn(&Kernel) -> Result<Brep, Error>;
    let parts: [(&str, Make); 2] = [
        ("box cut by a sphere, sphere faces d=+2", |k| {
            box_of_height(k, H)?
                .intersection(Brep::sphere(k, [0.0; 3], 55.0)?)?
                .offset_where(&FaceQuery::Kind(SurfaceKind::Sphere), 2.0)
        }),
        ("filleted box, top and every wall d=+5", |k| {
            top_rounded_box(k, H, 8.0)?.offset_where(&top_at(H).or(walls()), 5.0)
        }),
    ];
    for (name, make) in parts {
        let runs = on_kernel(move |k| {
            let mut kept = Vec::new();
            let mut runs = Vec::new();
            for _ in 0..8 {
                let moved = make(k)?;
                let edges: Vec<_> = moved.edges()?.iter().map(|e| e.samples).collect();
                runs.push((edges, moved.volume(PRECISE)?.to_bits()));
                kept.push(moved);
            }
            Ok::<_, Error>(runs)
        })
        .unwrap();
        for run in &runs[1..] {
            assert!(*run == runs[0], "{name}: an offset came out different");
        }
    }
}

#[test]
fn stretched_seams_are_exact() {
    // A hole's seam runs from the box's base up to the floor of a pocket cut
    // around the hole, so it stretches when the floor moves. Mated, it came
    // from one side's lifted trim curve and moved its base end by about 1e-14;
    // it must come out as the part built at its new depth has it.
    let pocketed = |k: &Kernel, depth: f64| {
        holed_box(k, H, 10.0)?.difference(Brep::cuboid(
            k,
            [-20.0, -20.0, H / 2.0 - depth],
            [40.0, 40.0, depth + 1.0],
        )?)
    };
    let seam = |b: &Brep| -> Result<_, Error> {
        let seams: Vec<_> = b
            .edges()?
            .into_iter()
            .filter(|e| e.faces[0].is_some() && e.faces[0] == e.faces[1])
            .map(|e| e.samples)
            .collect();
        assert_eq!(seams.len(), 1, "the hole has one seam");
        Ok(seams[0])
    };
    for d in [-5.0, 5.0] {
        let (before, moved, reference) = on_kernel(move |k| {
            let part = pocketed(k, 15.0)?;
            let before = seam(&part)?;
            let floor = FaceQuery::OnPlane {
                point: [0.0, 0.0, H / 2.0 - 15.0],
                normal: UP,
            };
            let moved = part.offset_where(&floor, d)?;
            assert!(moved.is_manifold()?, "floor d={d:+}: not manifold");
            Ok::<_, Error>((before, seam(&moved)?, seam(&pocketed(k, 15.0 - d)?)?))
        })
        .unwrap();
        assert_eq!(moved, reference, "floor d={d:+}: the seam is off");
        assert_eq!(
            moved[0], before[0],
            "floor d={d:+}: the seam's base end moved"
        );
    }
}

#[test]
fn edges_read_after_an_offset_and_a_fillet() {
    // A boss raised by a sweep, the far end trimmed, then the boss's rim
    // rounded: the fillet leaves the seam's range an ulp past its line's
    // domain, which SmApiEvaluateCurve refuses to evaluate.
    let (edges, manifold) = on_kernel(|k| {
        let mut part = Brep::cuboid(k, [-100.0, -60.0, -40.0], [200.0, 120.0, 80.0])?
            .union(Brep::cylinder(k, [-50.0, 0.0, 40.0], 15.0, 20.0)?)?;
        let boss_top = FaceQuery::OnPlane {
            point: [0.0, 0.0, 60.0],
            normal: UP,
        };
        let boss_top = part.find_faces(&boss_top, &Tolerance::default())?;
        let sweep = part.sweep_face(boss_top[0], 10.0)?;
        part = part.union(sweep.prism)?;
        part.merge_faces()?;
        let end = part.find_faces(&side(neg(X)), &Tolerance::default())?;
        let rim = EdgeQuery::OnCircle {
            center: [-50.0, 0.0, 40.0],
            axis: UP,
            radius: 15.0,
        };
        let part = part
            .offset_faces(&end, -10.0)?
            .fillet(&[EdgeBlend::round(rim, 3.0)], &options())?;
        Ok::<_, Error>((part.edges()?.len(), part.is_manifold()?))
    })
    .unwrap();
    assert_eq!(edges, 17);
    assert!(manifold);
}

#[test]
fn offset_faces_refuses_rather_than_return_a_wrong_solid() {
    // Moves the kernel gets wrong today. Each must either come back exact or
    // fail with a clean error, never as a broken or inconsistent solid:
    //
    // - fillets sliding along drafted walls: the executive's result
    //   tessellates right but its own volume is 20-30% off;
    // - a box cut by a sphere, its sphere faces moved in, or its top moved
    //   up: the executive leaves the solid open.
    type Make = Box<dyn FnOnce(&Kernel) -> Result<Brep, Error> + Send>;
    let rounded_drafted = |k: &Kernel, h: f64| {
        drafted_block(k, h, 5.0)?.fillet(&[EdgeBlend::round(top_edges_at(h), 3.0)], &options())
    };
    let ball = |k: &Kernel, r: f64, h: f64| {
        box_of_height(k, h)?.intersection(Brep::sphere(k, [0.0; 3], r)?)
    };
    let cases: Vec<(String, Make, Make)> = vec![
        (
            "drafted block with rounded top, d=+10".into(),
            Box::new(move |k| rounded_drafted(k, H)?.offset_where(&top_at(H), 10.0)),
            Box::new(move |k| rounded_drafted(k, H + 10.0)),
        ),
        (
            "drafted block with rounded top, d=-10".into(),
            Box::new(move |k| rounded_drafted(k, H)?.offset_where(&top_at(H), -10.0)),
            Box::new(move |k| rounded_drafted(k, H - 10.0)),
        ),
        (
            "box cut by a sphere, sphere faces d=-2".into(),
            Box::new(move |k| {
                ball(k, 55.0, H)?.offset_where(&FaceQuery::Kind(SurfaceKind::Sphere), -2.0)
            }),
            Box::new(move |k| ball(k, 53.0, H)),
        ),
        (
            "box cut by a sphere, top d=+5".into(),
            Box::new(move |k| ball(k, 55.0, H)?.offset_where(&top_at(H), 5.0)),
            Box::new(move |k| ball(k, 55.0, H + 5.0)),
        ),
    ];
    for (name, moved, reference) in cases {
        let outcome = on_kernel(move |k| {
            let reference = reference(k).unwrap().volume(PRECISE).unwrap();
            moved(k).map(|b| {
                (
                    b.is_manifold().unwrap(),
                    b.volume(PRECISE).unwrap(),
                    reference,
                )
            })
        });
        match outcome {
            Ok((manifold, volume, reference)) => {
                let error = (volume - reference) / reference;
                println!("{name}: exact, {error:+.1e}");
                assert!(
                    manifold && error.abs() < 1e-6,
                    "{name}: manifold {manifold}, {error:e} off"
                );
            }
            Err(error) => {
                println!("{name}: refused, {error}");
                assert!(
                    matches!(
                        status_name(&error),
                        Some("SHIM_ERR_BROKE_SOLID" | "SHIM_ERR_CHANGED_TOPOLOGY")
                    ),
                    "{name}: {error}"
                );
            }
        }
    }
}

/// A crate-sized block, its vertical edges rounded, with 16 mm bosses 6 high
/// (`r` 8) or holes 20 deep on its top at each (x, y), as in game-cad's
/// bugs.gcad: two stand near the front edge.
fn crate_with(k: &Kernel, bosses: &[[f64; 3]], holes: &[[f64; 3]]) -> Result<Brep, Error> {
    let mut part = Brep::cuboid(k, [-200.0, -150.0, -125.0], [400.0, 300.0, 265.0])?.fillet(
        &[EdgeBlend::round(
            EdgeQuery::ParallelTo { direction: UP },
            20.0,
        )],
        &options(),
    )?;
    for &[x, y, r] in bosses {
        part = part.union(Brep::cylinder(k, [x, y, 140.0], r, 6.0)?)?;
    }
    for &[x, y, r] in holes {
        part = part.difference(Brep::cylinder(k, [x, y, 120.0], r, 20.0)?)?;
    }
    Ok(part)
}

const FRONT_PAIR: [[f64; 3]; 2] = [[-165.0, -115.0, 8.0], [165.0, -115.0, 8.0]];

fn crate_front() -> FaceQuery {
    FaceQuery::OnPlane {
        point: [0.0, -150.0, 0.0],
        normal: neg(Y),
    }
}

/// The area of a disc of radius 8 at y = -115 left behind the front face
/// moved to `y`: what remains of a front boss or hole.
fn left_behind(y: f64) -> f64 {
    let r: f64 = 8.0;
    let segment = |d: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let from_centre = y + 115.0;
    if from_centre <= -r {
        PI * r * r
    } else if from_centre >= r {
        0.0
    } else if from_centre >= 0.0 {
        segment(from_centre)
    } else {
        PI * r * r - segment(-from_centre)
    }
}

#[test]
fn a_pushed_face_cuts_through_what_stands_in_its_way() {
    // The front face pushed in under two bosses: short of them, through
    // them, and past them. The bosses come back trimmed flush with the moved
    // face, keeping their labels, and the rest is as if they weren't there;
    // so with holes, cut back as notches.
    for (holes, height) in [(false, 6.0), (true, -20.0)] {
        for d in [-26.0, -30.0, -46.0] {
            let (with, without, obstacles, faces, labels) = on_kernel(move |k| {
                let (bosses, cuts): (&[[f64; 3]], &[[f64; 3]]) = if holes {
                    (&[], &FRONT_PAIR)
                } else {
                    (&FRONT_PAIR, &[])
                };
                let mut part = crate_with(k, bosses, cuts)?;
                for face in part.faces()? {
                    if matches!(face.surface, Surface::Cylinder { radius, .. } if (radius - 8.0).abs() < 1e-6)
                    {
                        part.set_face_label(face.index, &Label::new(9, "side"))?;
                    }
                }
                let front = part.find_faces(&crate_front(), &Tolerance::default())?;
                part.set_face_label(front[0], &Label::new(1, "-Y"))?;
                let obstacles = part.offset_obstacles(&front, d)?.len();
                let moved = part.offset_faces(&front, d)?;
                let plain = crate_with(k, &[], &[])?.offset_where(&crate_front(), d)?;
                let labels: Vec<String> = moved
                    .faces()?
                    .iter()
                    .filter_map(|f| f.label.as_ref().map(ToString::to_string))
                    .collect();
                Ok::<_, Error>((
                    moved.volume(PRECISE)?,
                    plain.volume(PRECISE)?,
                    obstacles,
                    moved.face_count(),
                    labels,
                ))
            })
            .unwrap();
            let left = 2.0 * height * left_behind(-150.0 - d);
            let expected = without + left;
            println!(
                "holes {holes} d={d:+}: {faces} faces, {with:.4} against {expected:.4}, {obstacles} in the way"
            );
            assert!(
                ((with - expected) / expected).abs() < PRECISE,
                "holes {holes} d={d}: {with} against {expected}"
            );
            assert_eq!(obstacles, if d == -26.0 { 0 } else { 2 }, "d={d}");
            assert!(
                labels.iter().any(|l| l == "1/-Y.0"),
                "the moved face keeps its label"
            );
            let pieces = labels.iter().filter(|l| *l == "9/side.0").count();
            match d {
                -46.0 => assert_eq!(pieces, 0, "gone past the plane"),
                _ => assert!(
                    pieces >= 2,
                    "the trimmed {} keep their labels: {labels:?}",
                    if holes { "holes" } else { "bosses" }
                ),
            }
        }
    }
}

fn crate_lid() -> FaceQuery {
    FaceQuery::OnPlane {
        point: [0.0, 0.0, 140.0],
        normal: UP,
    }
}

/// The part of a disc of radius `r` beyond a chord `d` from its centre, the
/// part a cut there leaves on the far side.
fn disc_beyond(r: f64, d: f64) -> f64 {
    if d >= r {
        return 0.0;
    }
    PI * r * r - (r * r * (d / r).acos() - d * (r * r - d * d).sqrt())
}

#[test]
fn a_moved_face_swallows_or_cuts_what_it_passes() {
    // What stands on the moved face itself: the front pushed past the floor
    // of a pocket 30 deep in it, which then goes, and the lid pulled up past
    // the 6 mm bosses on it, which it swallows. Short of that, the faces just
    // move: the pocket is shallower and the bosses shorter. Either way the
    // part is the plain crate moved the same, less what's left of the pocket
    // or plus what's left of the bosses.
    const BOSSES: [[f64; 3]; 4] = [
        [-165.0, -115.0, 8.0],
        [165.0, -115.0, 8.0],
        [-165.0, 115.0, 8.0],
        [165.0, 115.0, 8.0],
    ];
    let pocket = |k: &Kernel| {
        crate_with(k, &[], &[])?.difference(Brep::cuboid(
            k,
            [-60.0, -160.0, 40.0],
            [120.0, 40.0, 40.0],
        )?)
    };
    type Make = Box<dyn Fn(&Kernel) -> Result<Brep, Error> + Send>;
    let cases: Vec<(&str, Make, FaceQuery, f64, usize, f64)> = vec![
        (
            "front short of the pocket's floor",
            Box::new(pocket),
            crate_front(),
            -20.0,
            0,
            -120.0 * 40.0 * 10.0,
        ),
        (
            "front past the pocket's floor",
            Box::new(pocket),
            crate_front(),
            -36.0,
            1,
            0.0,
        ),
        (
            "lid short of the bosses' tops",
            Box::new(|k| crate_with(k, &BOSSES, &[])),
            crate_lid(),
            3.0,
            0,
            4.0 * PI * 64.0 * 3.0,
        ),
        (
            "lid past the bosses' tops",
            Box::new(|k| crate_with(k, &BOSSES, &[])),
            crate_lid(),
            10.0,
            4,
            0.0,
        ),
    ];
    for (name, make, face, d, in_the_way, left) in cases {
        let (volume, plain, obstacles) = on_kernel(move |k| {
            let part = make(k)?;
            let picked = part.find_faces(&face, &Tolerance::default())?;
            let obstacles = part.offset_obstacles(&picked, d)?.len();
            let moved = part.offset_faces(&picked, d)?;
            assert!(moved.is_manifold()?);
            let plain = crate_with(k, &[], &[])?.offset_where(&face, d)?;
            Ok::<_, Error>((moved.volume(PRECISE)?, plain.volume(PRECISE)?, obstacles))
        })
        .unwrap();
        let expected = plain + left;
        println!("{name}: {volume:.4} against {expected:.4}, {obstacles} in the way");
        assert!(
            ((volume - expected) / expected).abs() < PRECISE,
            "{name}: {volume} against {expected}"
        );
        assert_eq!(obstacles, in_the_way, "{name}");
    }

    // A vent 25 in radius through the lid, its rim chamfered 3, as in
    // game-cad's sample crate: the front pushed in to 22.2 short of its axis
    // cuts both its wall and its chamfer, leaving what's behind the plane.
    // Cutting the vent's space back out of the moved part instead used to
    // leave a solid 45 mm out of true; copying the result back with
    // CopyFaces, one inside out.
    let (volume, plain) = on_kernel(|k| {
        let rim = EdgeQuery::OnCircle {
            center: [0.0, 0.0, 140.0],
            axis: UP,
            radius: 25.0,
        };
        let part = crate_with(k, &[], &[])?
            .difference(Brep::cylinder(k, [0.0, 0.0, 60.0], 25.0, 90.0)?)?
            .chamfer(&rim, 3.0 * 2.0_f64.sqrt(), &Tolerance::default())?;
        let moved = part.offset_where(&crate_front(), -127.8)?;
        assert!(moved.is_manifold()?);
        moved.tessellate(&TessellationParams::default())?;
        let plain = crate_with(k, &[], &[])?.offset_where(&crate_front(), -127.8)?;
        Ok::<_, Error>((moved.volume(PRECISE)?, plain.volume(PRECISE)?))
    })
    .unwrap();
    let steps = 300;
    let chamfer: f64 = (0..steps)
        .map(|i| {
            disc_beyond(25.0 + (f64::from(i) + 0.5) * 3.0 / f64::from(steps), 22.2) * 3.0
                / f64::from(steps)
        })
        .sum();
    let expected = plain - 77.0 * disc_beyond(25.0, 22.2) - chamfer;
    println!("chamfered vent: {volume:.4} against {expected:.4}");
    assert!(
        ((volume - expected) / expected).abs() < 1e-6,
        "chamfered vent: {volume} against {expected}"
    );
}

#[test]
fn rounds_pushed_over_pockets_cut_them_beside_a_chamfered_vent() {
    // The front pushed in takes its corner rounds over the ends of pockets 30
    // deep, 120 long and 40 high in the sides, which then open into the
    // rounds, while the chamfered vent through the lid stays as it was: at
    // -115 the front's new plane comes within 7 of the chamfer. The part is
    // the crate with only the vent, pushed the same, less what's left of the
    // pockets.
    const DEPTH: f64 = 30.0;
    for d in [-80.0, -115.0] {
        let (volume, plain) = on_kernel(move |k| {
            let rim = EdgeQuery::OnCircle {
                center: [0.0, 0.0, 140.0],
                axis: UP,
                radius: 25.0,
            };
            let vented = crate_with(k, &[], &[])?
                .difference(Brep::cylinder(k, [0.0, 0.0, 60.0], 25.0, 90.0)?)?
                .chamfer(&rim, 3.0 * 2.0_f64.sqrt(), &Tolerance::default())?;
            let part = vented
                .try_clone()?
                .difference(Brep::cuboid(
                    k,
                    [200.0 - DEPTH, -60.0, 0.0],
                    [DEPTH + 10.0, 120.0, 40.0],
                )?)?
                .difference(Brep::cuboid(
                    k,
                    [-210.0, -60.0, 0.0],
                    [DEPTH + 10.0, 120.0, 40.0],
                )?)?;
            let moved = part.offset_where(&crate_front(), d)?;
            assert!(moved.is_manifold()?);
            moved.tessellate(&TessellationParams::default())?;
            let plain = vented.offset_where(&crate_front(), d)?;
            Ok::<_, Error>((moved.volume(PRECISE)?, plain.volume(PRECISE)?))
        })
        .unwrap();
        // The rounds, radius 20, end at the moved front and centre 20 behind
        // it: a pocket is as deep as a round leaves it over the first stretch
        // past the front, and DEPTH from the rounds' centres on.
        let front = -150.0 - d;
        let centre = front + 20.0;
        let under = centre - front.max(-60.0);
        let round = under / 2.0 * (400.0 - under * under).sqrt() + 200.0 * (under / 20.0).asin();
        let pocket = 40.0 * (DEPTH * (60.0 - centre) + (DEPTH - 20.0) * under + round);
        let expected = plain - 2.0 * pocket;
        println!("d={d}: {volume:.4} against {expected:.4}");
        assert!(
            ((volume - expected) / expected).abs() < 1e-6,
            "d={d}: {volume} against {expected}"
        );
    }
}

/// The area two discs of radii `a` and `b`, `d` apart, share.
fn lens(a: f64, b: f64, d: f64) -> f64 {
    if d >= a + b {
        return 0.0;
    }
    if d <= (a - b).abs() {
        return PI * a.min(b).powi(2);
    }
    let alpha = ((d * d + a * a - b * b) / (2.0 * d * a)).acos();
    let beta = ((d * d + b * b - a * a) / (2.0 * d * b)).acos();
    let kite = ((-d + a + b) * (d + a - b) * (d - a + b) * (d + a + b)).sqrt() / 2.0;
    a * a * alpha + b * b * beta - kite
}

/// The cylindrical face of radius `r` whose axis runs through (x, y).
fn cylinder_at(part: &Brep, x: f64, y: f64, r: f64) -> Result<usize, Error> {
    Ok(part
        .faces()?
        .into_iter()
        .find(|f| {
            matches!(f.surface, Surface::Cylinder { origin, radius, .. }
                if (radius - r).abs() < 1e-6 && (origin[0] - x).hypot(origin[1] - y) < 1e-6)
        })
        .unwrap_or_else(|| panic!("no cylinder of radius {r} at ({x}, {y})"))
        .index)
}

#[test]
fn a_moved_curved_face_cuts_what_it_passes_along_its_new_surface() {
    // Parts on the crate's lid, at z = 140, each with the curved face moved,
    // how many features are in its way, and its volume as it should come out.
    type Make = Box<dyn Fn(&Kernel) -> Result<Brep, Error> + Send>;
    type Expected = Box<dyn Fn(&Kernel) -> Result<f64, Error> + Send>;
    // The crate with bosses (x, y, radius, height) on its lid and holes (x,
    // y, radius) 20 deep in it.
    fn lid(k: &Kernel, bosses: &[[f64; 4]], holes: &[[f64; 3]]) -> Result<Brep, Error> {
        let mut part = crate_with(k, &[], &[])?;
        for &[x, y, r, h] in bosses {
            part = part.union(Brep::cylinder(k, [x, y, 140.0], r, h)?)?;
        }
        for &[x, y, r] in holes {
            part = part.difference(Brep::cylinder(k, [x, y, 120.0], r, 30.0)?)?;
        }
        Ok(part)
    }
    fn volume(k: &Kernel, bosses: &[[f64; 4]], holes: &[[f64; 3]]) -> Result<f64, Error> {
        lid(k, bosses, holes)?.volume(PRECISE)
    }
    struct Case {
        name: &'static str,
        part: Make,
        // The cylinder moved: its axis's (x, y) and radius.
        face: [f64; 3],
        d: f64,
        in_the_way: usize,
        expected: Expected,
    }
    let cases = vec![
        // A hole of radius 20 widened to 38 and to 50 past a boss of radius 8
        // beside it, 40 out: the boss is cut along the new wall, then gone.
        // Short of the boss, nothing's in the way.
        Case {
            name: "hole widened short of the boss",
            part: Box::new(move |k| lid(k, &[[40.0, 0.0, 8.0, 6.0]], &[[0.0, 0.0, 20.0]])),
            face: [0.0, 0.0, 20.0],
            d: -5.0,
            in_the_way: 0,
            expected: Box::new(move |k| volume(k, &[[40.0, 0.0, 8.0, 6.0]], &[[0.0, 0.0, 25.0]])),
        },
        Case {
            name: "hole widened through the boss",
            part: Box::new(move |k| lid(k, &[[40.0, 0.0, 8.0, 6.0]], &[[0.0, 0.0, 20.0]])),
            face: [0.0, 0.0, 20.0],
            d: -18.0,
            in_the_way: 1,
            expected: Box::new(move |k| {
                Ok(
                    volume(k, &[], &[[0.0, 0.0, 38.0]])?
                        + 6.0 * (PI * 64.0 - lens(38.0, 8.0, 40.0)),
                )
            }),
        },
        Case {
            name: "hole widened past the boss",
            part: Box::new(move |k| lid(k, &[[40.0, 0.0, 8.0, 6.0]], &[[0.0, 0.0, 20.0]])),
            face: [0.0, 0.0, 20.0],
            d: -30.0,
            in_the_way: 1,
            expected: Box::new(move |k| volume(k, &[], &[[0.0, 0.0, 50.0]])),
        },
        // A post 20 in radius and 30 tall thinned to 12 past a pin on its
        // top, 15 out: the pin is cut along the post's new side.
        Case {
            name: "post thinned through the pin on it",
            part: Box::new(move |k| {
                lid(k, &[[0.0, 0.0, 20.0, 30.0]], &[])?.union(Brep::cylinder(
                    k,
                    [15.0, 0.0, 170.0],
                    3.0,
                    5.0,
                )?)
            }),
            face: [0.0, 0.0, 20.0],
            d: -8.0,
            in_the_way: 1,
            expected: Box::new(move |k| {
                Ok(volume(k, &[[0.0, 0.0, 12.0, 30.0]], &[])? + 5.0 * lens(12.0, 3.0, 15.0))
            }),
        },
        // A boss 10 in radius widened to 22 over a hole 5 in radius 25 out,
        // whose mouth it then half covers, and over a taller boss there,
        // which it joins.
        Case {
            name: "boss widened over a hole",
            part: Box::new(move |k| lid(k, &[[0.0, 0.0, 10.0, 6.0]], &[[25.0, 0.0, 5.0]])),
            face: [0.0, 0.0, 10.0],
            d: 12.0,
            in_the_way: 1,
            expected: Box::new(move |k| {
                Ok(volume(k, &[[0.0, 0.0, 10.0, 6.0]], &[[25.0, 0.0, 5.0]])?
                    + 6.0 * PI * (22.0 * 22.0 - 100.0))
            }),
        },
        Case {
            name: "boss widened over a boss",
            part: Box::new(move |k| lid(k, &[[0.0, 0.0, 10.0, 6.0], [25.0, 0.0, 5.0, 10.0]], &[])),
            face: [0.0, 0.0, 10.0],
            d: 12.0,
            in_the_way: 1,
            expected: Box::new(move |k| {
                let both = volume(k, &[[0.0, 0.0, 10.0, 6.0], [25.0, 0.0, 5.0, 10.0]], &[])?;
                Ok(both + 6.0 * (PI * (22.0 * 22.0 - 100.0) - lens(22.0, 5.0, 25.0)))
            }),
        },
    ];
    for case in cases {
        let Case {
            name,
            part,
            face,
            d,
            in_the_way,
            expected,
        } = case;
        let (volume, expected, obstacles) = on_kernel(move |k| {
            let part = part(k)?;
            let [x, y, r] = face;
            let picked = [cylinder_at(&part, x, y, r)?];
            let obstacles = part.offset_obstacles(&picked, d)?.len();
            let moved = part.offset_faces(&picked, d)?;
            assert!(moved.is_manifold()?, "{name}");
            moved.tessellate(&TessellationParams::default())?;
            Ok::<_, Error>((moved.volume(PRECISE)?, expected(k)?, obstacles))
        })
        .unwrap();
        println!("{name}: {volume:.4} against {expected:.4}, {obstacles} in the way");
        assert!(
            ((volume - expected) / expected).abs() < 1e-7,
            "{name}: {volume} against {expected}"
        );
        assert_eq!(obstacles, in_the_way, "{name}");
    }

    // Widened to 32, the boss covers the hole's mouth whole, which would
    // leave it a void inside the part: refused.
    let outcome = on_kernel(move |k| {
        let part = lid(k, &[[0.0, 0.0, 10.0, 6.0]], &[[25.0, 0.0, 5.0]])?;
        let picked = [cylinder_at(&part, 0.0, 0.0, 10.0)?];
        Ok::<_, Error>(part.offset_faces(&picked, 22.0).map(|b| b.face_count()))
    })
    .unwrap();
    let error = outcome.expect_err("the hole would be sealed in");
    assert_eq!(status_name(&error), Some("SHIM_ERR_OBSTRUCTED"), "{error}");
}

#[test]
fn a_pushed_face_refuses_what_it_cannot_set_aside() {
    // A small boss stands 0.5 mm behind a front one, inside the room the
    // front one is set aside with: the front one can't be taken out alone.
    let (outcome, obstacles) = on_kernel(|k| {
        let part = crate_with(k, &[[-165.0, -115.0, 8.0], [-165.0, -104.5, 2.0]], &[])?;
        let front = part.find_faces(&crate_front(), &Tolerance::default())?;
        let obstacles = part.offset_obstacles(&front, -30.0)?;
        let faces = part.faces()?;
        let radii: Vec<Vec<f64>> = obstacles
            .iter()
            .map(|o| {
                o.iter()
                    .filter_map(|&i| match faces[i].surface {
                        Surface::Cylinder { radius, .. } => Some(radius),
                        _ => None,
                    })
                    .collect()
            })
            .collect();
        Ok::<_, Error>((
            part.offset_faces(&front, -30.0).map(|b| b.face_count()),
            radii,
        ))
    })
    .unwrap();
    let error = outcome.expect_err("the small boss is in the way");
    assert_eq!(status_name(&error), Some("SHIM_ERR_OBSTRUCTED"), "{error}");
    assert_eq!(
        obstacles.len(),
        1,
        "only the front boss is reached: {obstacles:?}"
    );
    assert!(
        obstacles[0].iter().all(|r| (r - 8.0).abs() < 1e-6),
        "{obstacles:?}"
    );

    // A boss on the side wall near the front: the round between the front
    // and the wall moves in with the front, so the wall pulls back from under
    // the boss's remains, which would come back on their own. Past the new
    // plane altogether it simply goes.
    let outcomes = on_kernel(|k| {
        let mut outcomes = Vec::new();
        for d in [-30.0, -46.0] {
            let mut boss = Brep::cylinder(k, [0.0; 3], 8.0, 6.0)?;
            boss.place(&Frame {
                origin: [-200.0, -115.0, 0.0],
                x_axis: UP,
                y_axis: Y,
            })?;
            let part = crate_with(k, &[], &[])?.union(boss)?;
            let front = part.find_faces(&crate_front(), &Tolerance::default())?;
            outcomes.push(part.offset_faces(&front, d).map(|b| b.face_count()));
        }
        Ok::<_, Error>(outcomes)
    })
    .unwrap();
    let error = outcomes[0]
        .as_ref()
        .expect_err("the boss's remains would float");
    assert_eq!(status_name(error), Some("SHIM_ERR_OBSTRUCTED"), "{error}");
    assert_eq!(outcomes[1].as_ref().ok(), Some(&10), "gone past the plane");
}

#[test]
fn an_unchecked_offset_returns_what_the_kernel_made() {
    // Without the push through and the checks, the local operation leaves the
    // bosses' rims outside the lid: a solid that passes for closed but can't
    // be tessellated. The checked offset never returns it.
    let (unchecked, checked) = on_kernel(|k| {
        let part = crate_with(k, &FRONT_PAIR, &[])?;
        let front = part.find_faces(&crate_front(), &Tolerance::default())?;
        let raw = part.try_clone()?.offset_faces_unchecked(&front, -46.0)?;
        let unchecked = raw.tessellate(&TessellationParams::default()).map(|_| ());
        let checked = part
            .offset_faces(&front, -46.0)?
            .tessellate(&TessellationParams::default())
            .map(|_| ());
        Ok::<_, Error>((unchecked, checked))
    })
    .unwrap();
    assert!(unchecked.is_err(), "the kernel's own result tessellates");
    checked.unwrap();
}

#[test]
fn a_pushed_face_meets_the_face_beyond_one_it_swallows() {
    // A face pushed in past a chamfer along its edge swallows it, and meets
    // the face beyond directly, as if the edge had never been chamfered: the
    // chamfered box's top pushed past its chamfers (legs 2), to exactly
    // their width and beyond, is a plain box, and so is its +X wall pushed
    // past the chamfer above it, less that chamfer, the others running into
    // the wall's new plane.
    type Make = Box<dyn Fn(&Kernel) -> Result<Brep, Error> + Send>;
    let c: f64 = 2.0;
    let mut cases: Vec<(String, Make, f64, usize)> = Vec::new();
    for d in [-2.0, -2.5, -10.0] {
        cases.push((
            format!("chamfered box, top d={d:+}"),
            Box::new(move |k| {
                let mut part = chamfered_box(k, H, c)?;
                let top = part.find_faces(&top_at(H), &Tolerance::default())?;
                part.set_face_label(top[0], &Label::new(1, "top"))?;
                part.offset_faces(&top, d)
            }),
            L * W * (H + d),
            6,
        ));
    }
    let short = L - 5.0;
    cases.push((
        "chamfered box, +X wall d=-5".into(),
        Box::new(move |k| chamfered_box(k, H, c)?.offset_where(&side(X), -5.0)),
        short * W * H - c * c / 2.0 * (2.0 * short + W) + 2.0 * c.powi(3) / 3.0,
        9,
    ));
    for (name, make, expected, faces) in cases {
        let (volume, count, manifold, labels) = on_kernel(move |k| {
            let moved = make(k)?;
            moved.tessellate(&TessellationParams::default())?;
            let labels: Vec<String> = moved
                .faces()?
                .iter()
                .filter_map(|f| f.label.as_ref().map(ToString::to_string))
                .collect();
            Ok::<_, Error>((
                moved.volume(PRECISE)?,
                moved.face_count(),
                moved.is_manifold()?,
                labels,
            ))
        })
        .unwrap_or_else(|e| panic!("{name}: {e}"));
        println!("{name}: {volume:.6} against {expected:.6}, {count} faces");
        assert!(manifold, "{name}");
        assert!(
            ((volume - expected) / expected).abs() < 1e-9,
            "{name}: {volume} against {expected}"
        );
        assert_eq!(count, faces, "{name}");
        if name.contains("top") {
            assert_eq!(
                labels,
                ["1/top.0"],
                "{name}: the moved face keeps its label"
            );
        }
    }

    // Curved faces the same: game-cad's sample crate's vent, 25 in radius
    // and 80 deep with its rim chamfered 3, widened to the chamfer's width
    // and past it, is a plain hole up to the lid; and a post 20 in radius and
    // 30 tall, its top rim chamfered 3, thinned so is a plain post.
    let rim = |z: f64, r: f64| EdgeQuery::OnCircle {
        center: [0.0, 0.0, z],
        axis: UP,
        radius: r,
    };
    type Expected = Box<dyn Fn(&Kernel) -> Result<f64, Error> + Send>;
    let mut cases: Vec<(String, Make, f64, Expected)> = Vec::new();
    for d in [-3.0, -3.5, -10.0] {
        cases.push((
            format!("vent widened d={d:+}"),
            Box::new(move |k| {
                crate_with(k, &[], &[])?
                    .difference(Brep::cylinder(k, [0.0, 0.0, 60.0], 25.0, 90.0)?)?
                    .chamfer(
                        &rim(140.0, 25.0),
                        3.0 * 2.0_f64.sqrt(),
                        &Tolerance::default(),
                    )
            }),
            d,
            Box::new(move |k| {
                crate_with(k, &[], &[])?
                    .difference(Brep::cylinder(k, [0.0, 0.0, 60.0], 25.0 - d, 90.0)?)?
                    .volume(PRECISE)
            }),
        ));
    }
    for d in [-3.0, -5.0] {
        cases.push((
            format!("post thinned d={d:+}"),
            Box::new(move |k| {
                crate_with(k, &[], &[])?
                    .union(Brep::cylinder(k, [0.0, 0.0, 140.0], 20.0, 30.0)?)?
                    .chamfer(
                        &rim(170.0, 20.0),
                        3.0 * 2.0_f64.sqrt(),
                        &Tolerance::default(),
                    )
            }),
            d,
            Box::new(move |k| {
                crate_with(k, &[], &[])?
                    .union(Brep::cylinder(k, [0.0, 0.0, 140.0], 20.0 + d, 30.0)?)?
                    .volume(PRECISE)
            }),
        ));
    }
    for (name, make, d, expected) in cases {
        let r = if name.starts_with("vent") { 25.0 } else { 20.0 };
        let (volume, expected, faces, plain, manifold) = on_kernel(move |k| {
            let part = make(k)?;
            let picked = [cylinder_at(&part, 0.0, 0.0, r)?];
            let moved = part.offset_faces(&picked, d)?;
            moved.tessellate(&TessellationParams::default())?;
            let plain = crate_with(k, &[], &[])?.face_count();
            Ok::<_, Error>((
                moved.volume(PRECISE)?,
                expected(k)?,
                moved.face_count(),
                plain,
                moved.is_manifold()?,
            ))
        })
        .unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(manifold, "{name}");
        println!("{name}: {volume:.6} against {expected:.6}, {faces} faces");
        assert!(
            ((volume - expected) / expected).abs() < 1e-8,
            "{name}: {volume} against {expected}"
        );
        assert_eq!(
            faces,
            plain + 2,
            "{name}: the crate's faces, a wall and a floor or top"
        );
    }
}
