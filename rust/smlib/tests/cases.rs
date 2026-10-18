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
    Brep, Continuity, EdgeBlend, EdgeQuery, Error, FaceQuery, FilletOptions, Kernel, Mesh, Surface,
    SurfaceKind, TessellationParams, Tolerance,
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
    // Each move would shrink a face to nothing: the chamfers (legs 2) from
    // d = -2, the walls under the fillets from d = -32, the box's walls at
    // d = -40, the drafted block's top past its apex (about d = +303). The
    // kernel hands back an open or inside-out solid; the shim refuses it.
    let outcomes = on_kernel(move |k| {
        let attempts: Vec<(&str, f64, Result<Brep, Error>)> = vec![
            ("chamfered box", -2.5, chamfered_box(k, H, 2.0)),
            ("chamfered box", -10.0, chamfered_box(k, H, 2.0)),
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
