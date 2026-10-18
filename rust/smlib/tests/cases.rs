//! The demo's four cases on a 100 × 60 × 40 mm box centred on the origin,
//! Z up, in kernel millimetres.
//!
//! The test harness runs each test on its own thread, but the kernel belongs
//! to one thread, so every test sends its kernel work to a shared kernel
//! thread and checks the plain-data results on its own.

use std::collections::HashSet;
use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{Brep, Error, Kernel, Mesh, Surface, TessellationParams};

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
