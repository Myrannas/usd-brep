//! The demo's cases on a 100 × 60 × 40 mm box centred on the origin, Z up,
//! in kernel millimetres: four booleans, then fillets and chamfers on edges
//! picked by geometric queries.
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
    Brep, Continuity, EdgeBlend, EdgeQuery, Error, FilletOptions, Kernel, Mesh, Surface,
    TessellationParams, Tolerance,
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
        (box_counts, hole_counts, edges)
    });
    let (box_counts, [hole_edges, seams, sharp, rim, wrong_rim, bottom_rim], edges) = counts;
    println!(
        "box picks {box_counts:?}; hole: {hole_edges} edges, {seams} seams, {sharp} sharp, rim {rim}"
    );
    assert_eq!(box_counts, [12, 12, 4, 4, 1, 1, 4, 4, 1]);
    assert!(seams >= 1, "the hole's cylinder should have a seam");
    assert_eq!(sharp, hole_edges - seams, "seams aren't sharp");
    assert!(rim >= 1, "the top rim");
    assert_eq!(wrong_rim, 0, "a circle 0.5 mm too wide picks nothing");
    assert_eq!(bottom_rim, rim, "top and bottom rims split alike");
    for edge in &edges {
        assert_eq!(edge.use_count, 2);
        assert!(edge.length > 0.0);
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
    // SmApiChamferFillet's solver handles small rims, but from about d = 2.5
    // to 4 on this one it reports success and leaves an open shell. The shim
    // turns that into an error rather than hand back a broken solid.
    let r = 10.0;
    for d in [1.0, 2.0] {
        let built = build(move |k| {
            hole(k, r)?.chamfer(
                &top_rim(0.0, 0.0, r),
                d * 2f64.sqrt(),
                &Tolerance::default(),
            )
        });
        let name = format!("constant-width rim d={d}");
        assert!(built.manifold, "{name}");
        assert_volume(
            &name,
            built.volume,
            BOX_VOLUME - PI * r * r * H - rim_chamfer_volume(r, d),
        );
        assert_mesh(&name, &built);
    }
    let error = on_kernel(move |k| {
        hole(k, r)
            .unwrap()
            .chamfer(
                &top_rim(0.0, 0.0, r),
                3.0 * 2f64.sqrt(),
                &Tolerance::default(),
            )
            .map(|_| ())
            .unwrap_err()
    });
    println!("constant-width rim d=3: {error}");
    assert!(
        matches!(error, Error::Kernel { status, .. } if status.name() == Some("SHIM_ERR_BROKE_SOLID")),
        "{error:?}"
    );
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
