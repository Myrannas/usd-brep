//! Solids hollowed to a wall thickness, opened at chosen faces: what a
//! part's Hollow needs. Faces that stay must keep their labels, every inner
//! face must say which face it came from, and the result must be closed.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, EdgeBlend, EdgeQuery, Error, FaceGeometry, FilletOptions, Kernel, Label, Mesh, ShellFace,
    Status, Surface, TessellationParams,
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

const BLOCK: u32 = 1;

/// A 100 × 60 × 40 block from the origin, each face labelled by its index.
fn block(k: &Kernel) -> Brep {
    let mut block = Brep::cuboid(k, [0.0; 3], [100.0, 60.0, 40.0]).unwrap();
    for face in 0..block.face_count() {
        let label = Label::new(BLOCK, "face").with_index(face as u32);
        block.set_face_label(face, &label).unwrap();
    }
    block
}

/// The block with its vertical edges rounded 8 and a 6 radius hole down
/// through it.
fn rounded(k: &Kernel) -> Brep {
    let vertical = EdgeQuery::ParallelTo {
        direction: [0.0, 0.0, 1.0],
    };
    let rounded = block(k)
        .fillet(
            &[EdgeBlend::round(vertical, 8.0)],
            &FilletOptions::default(),
        )
        .unwrap();
    let pin = Brep::cylinder(k, [50.0, 30.0, -10.0], 6.0, 60.0).unwrap();
    rounded.difference(pin).unwrap()
}

/// The faces facing along `normal`.
fn facing(faces: &[FaceGeometry], normal: [f64; 3]) -> Vec<usize> {
    faces
        .iter()
        .filter(|f| f.normal == Some(normal))
        .map(|f| f.index)
        .collect()
}

/// What the tests look at in a hollowed solid, with its source's faces.
#[derive(Debug)]
struct Hollowed {
    sources: Vec<FaceGeometry>,
    roles: Vec<ShellFace>,
    faces: Vec<FaceGeometry>,
    manifold: bool,
    volume: f64,
    mesh: Mesh,
}

fn hollow(make: fn(&Kernel) -> Brep, thickness: f64, open: Vec<[f64; 3]>) -> Hollowed {
    on_kernel(move |k| {
        let solid = make(k);
        let sources = solid.faces().unwrap();
        let open: Vec<usize> = open.iter().flat_map(|&n| facing(&sources, n)).collect();
        let shell = solid.shell(thickness, &open).unwrap();
        Hollowed {
            sources,
            roles: shell.faces,
            faces: shell.brep.faces().unwrap(),
            manifold: shell.brep.is_manifold().unwrap(),
            volume: shell.brep.volume(1e-9).unwrap(),
            mesh: shell
                .brep
                .tessellate(&TessellationParams::default())
                .unwrap(),
        }
    })
}

/// What every hollowed solid must be: closed, of `volume`, its mesh winding
/// outwards round it; faces that stay keep their source's label, inner
/// faces have none, and each lies `thickness` in from its source when that's
/// a plane, on the same axis when a cylinder.
fn assert_hollowed(hollowed: &Hollowed, thickness: f64, volume: f64) {
    assert!(hollowed.manifold);
    assert!(
        close(hollowed.volume, volume),
        "{} vs {volume}",
        hollowed.volume
    );
    for (face, role) in hollowed.faces.iter().zip(&hollowed.roles) {
        match *role {
            ShellFace::Outer { face: source } | ShellFace::Rim { face: source } => {
                assert_eq!(face.label, hollowed.sources[source].label);
            }
            ShellFace::Inner { face: source } => {
                assert_eq!(face.label, None);
                assert_offset(face, &hollowed.sources[source], thickness);
            }
            ShellFace::New => panic!("a new face: {face:?}"),
        }
    }
    assert_mesh(&hollowed.mesh, hollowed.volume);
}

/// That `inner` lies `thickness` inside `source`: a parallel plane, or a
/// cylinder on the same axis.
fn assert_offset(inner: &FaceGeometry, source: &FaceGeometry, thickness: f64) {
    match (&inner.surface, &source.surface) {
        (Surface::Plane { origin, .. }, Surface::Plane { .. }) => {
            let normal = source.normal.unwrap();
            let [a, b] = [inner.normal.unwrap(), normal];
            // The inner face faces back into the wall it lines.
            assert!((dot(a, b) + 1.0).abs() < 1e-9, "{a:?} against {b:?}");
            let on_source = source.bounds_centre();
            let depth = -dot(sub(*origin, on_source), normal);
            assert!((depth - thickness).abs() < 1e-6, "{depth}");
        }
        (
            Surface::Cylinder { origin, radius, .. },
            Surface::Cylinder {
                origin: o,
                radius: r,
                axis,
                ..
            },
        ) => {
            let d = sub(*origin, *o);
            let off = sub(d, scale(*axis, dot(d, *axis)));
            assert!(dot(off, off).sqrt() < 1e-6, "{origin:?}");
            assert!((radius - r).abs() - thickness < 1e-6, "{radius} from {r}");
        }
        (inner, source) => panic!("{inner:?} from {source:?}"),
    }
}

#[test]
fn a_box_opened_at_its_top() {
    let hollowed = hollow(block, 5.0, vec![[0.0, 0.0, 1.0]]);
    // The inset: 90 × 50, from 5 up to the top.
    assert_hollowed(&hollowed, 5.0, 240_000.0 - 90.0 * 50.0 * 35.0);
    let count = |f: fn(&ShellFace) -> bool| hollowed.roles.iter().filter(|r| f(r)).count();
    assert_eq!(count(|r| matches!(r, ShellFace::Outer { .. })), 5);
    assert_eq!(count(|r| matches!(r, ShellFace::Rim { .. })), 1);
    assert_eq!(count(|r| matches!(r, ShellFace::Inner { .. })), 5);
    // Every face but the top has an inner face, and the top's rim is a
    // 5 wide ring round the opening.
    let top = facing(&hollowed.sources, [0.0, 0.0, 1.0])[0];
    let inner: Vec<usize> = hollowed
        .roles
        .iter()
        .filter_map(|r| match r {
            ShellFace::Inner { face } => Some(*face),
            _ => None,
        })
        .collect();
    assert_eq!(inner.len(), 5);
    assert!(!inner.contains(&top));
}

#[test]
fn a_box_opened_at_both_ends_is_a_tube() {
    let open = vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
    let hollowed = hollow(block, 5.0, open);
    assert_hollowed(&hollowed, 5.0, 240_000.0 - 90.0 * 50.0 * 40.0);
    let rims = hollowed
        .roles
        .iter()
        .filter(|r| matches!(r, ShellFace::Rim { .. }))
        .count();
    assert_eq!(rims, 2);
}

#[test]
fn rounds_and_holes_offset_round_their_axes() {
    // The rounds shrink to 8 - t round the same axes, the hole grows to 6 + t.
    let area = |w: f64, d: f64, round: f64, hole: f64| {
        w * d - (4.0 - PI) * round * round - PI * hole * hole
    };
    let outer = area(100.0, 60.0, 8.0, 6.0) * 40.0;
    for t in [2.0, 5.0, 7.9] {
        let hollowed = hollow(rounded, t, vec![[0.0, 0.0, 1.0]]);
        let inner = area(100.0 - 2.0 * t, 60.0 - 2.0 * t, 8.0 - t, 6.0 + t) * (40.0 - t);
        assert_hollowed(&hollowed, t, outer - inner);
        let cylinders = hollowed
            .faces
            .iter()
            .zip(&hollowed.roles)
            .filter(|(f, r)| {
                matches!(r, ShellFace::Inner { .. })
                    && matches!(f.surface, Surface::Cylinder { .. })
            })
            .count();
        // Four rounds and the hole, each in one piece.
        assert_eq!(cylinders, 5);
    }
}

#[test]
fn a_wall_as_thick_as_a_round_fails() {
    // A round offset in by its radius or more has nothing left.
    let failures = on_kernel(|k| {
        let solid = rounded(k);
        let top = facing(&solid.faces().unwrap(), [0.0, 0.0, 1.0]);
        [8.0, 10.0].map(|t| solid.shell(t, &top).err())
    });
    for failure in failures {
        assert!(matches!(failure, Some(Error::Kernel { .. })), "{failure:?}");
    }
}

#[test]
fn opening_a_face_rounded_into_its_neighbours_fails() {
    // The top's rim rounded 4: the offset comes back with faces whose own
    // integrals disagree with their edges, which the shim refuses.
    let failed = on_kernel(|k| {
        let rim = EdgeQuery::MidpointOnPlane {
            point: [0.0, 0.0, 40.0],
            normal: [0.0, 0.0, 1.0],
        };
        let solid = block(k)
            .fillet(&[EdgeBlend::round(rim, 4.0)], &FilletOptions::default())
            .unwrap();
        let top = facing(&solid.faces().unwrap(), [0.0, 0.0, 1.0]);
        solid.shell(2.0, &top).err()
    });
    assert_eq!(
        failed,
        Some(Error::Kernel {
            operation: "SmApiShellBrepFull",
            status: Status(2037),
        })
    );
}

#[test]
fn the_same_shell_comes_out_the_same() {
    // SMLib's offset mated edges in the order of their addresses (patch
    // 0afdfdb puts them in the brep's): the same shell gave different volumes
    // from one process to the next.
    let results = on_kernel(|k| {
        let solid = rounded(k);
        let top = facing(&solid.faces().unwrap(), [0.0, 0.0, 1.0]);
        (0..3)
            .map(|_| {
                let shell = solid.shell(5.0, &top).unwrap();
                (shell.brep.volume(1e-9).unwrap(), shell.brep.face_count())
            })
            .collect::<Vec<_>>()
    });
    assert!(results.windows(2).all(|w| w[0] == w[1]), "{results:?}");
}

#[test]
fn refuses_what_it_cant_hollow() {
    let reasons = on_kernel(|k| {
        let solid = block(k);
        let sheet = Brep::planar_faces(
            k,
            &smlib::Frame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap(),
            &[vec![smlib::BoundaryCurve {
                curve: smlib::PlaneCurve::Circle {
                    centre: [0.0, 0.0],
                    radius: 5.0,
                },
                label: None,
            }]],
        )
        .unwrap();
        let results = [
            solid.shell(0.0, &[0]),
            solid.shell(f64::NAN, &[0]),
            solid.shell(5.0, &[6]),
            solid.shell(5.0, &[1, 1]),
            solid.shell(5.0, &[0, 1, 2, 3, 4, 5]),
            solid.shell(5.0, &[]),
            sheet.shell(1.0, &[0]),
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
            "the thickness must be finite and positive",
            "the thickness must be finite and positive",
            "a face to open is out of range or repeated",
            "a face to open is out of range or repeated",
            "opening every face leaves nothing",
            "open at least one face: a closed hollow is a sealed cavity",
            "only a closed, manifold solid can be hollowed",
        ]
    );
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
