//! `Brep::move_faces`: faces moved and turned rigidly, their neighbours
//! extended or trimmed to follow, on a 100 × 60 × 40 mm box centred on the
//! origin, Z up. Each result is held to a manifold solid with the faces and
//! edges it had, their labels, a volume exact to 1e-7 against the shape's
//! own formula, and a mesh that encloses it and winds outward.
//!
//! The kernel belongs to one thread, so every test sends its kernel work to
//! a shared kernel thread and checks the plain-data results on its own.

use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, EdgeBlend, EdgeQuery, Error, FaceQuery, FilletOptions, Kernel, Label, Mesh, MovedFace,
    RigidMotion, Surface, SurfaceKind, TessellationParams, Tolerance,
};

const L: f64 = 100.0;
const W: f64 = 60.0;
const H: f64 = 40.0;
const PRECISE: f64 = 1e-9;

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

/// The box, each face labelled by its outward normal, as a Block's are.
fn block(k: &Kernel) -> Result<Brep, Error> {
    labelled(Brep::cuboid(k, [-L / 2.0, -W / 2.0, -H / 2.0], [L, W, H])?)
}

/// Labels every planar face by its outward normal ("+X", ...), and every
/// other face "curved" by its index.
fn labelled(mut brep: Brep) -> Result<Brep, Error> {
    for face in brep.faces()? {
        let role = match face.normal {
            Some(n) if n[0] > 0.99 => "+X",
            Some(n) if n[0] < -0.99 => "-X",
            Some(n) if n[1] > 0.99 => "+Y",
            Some(n) if n[1] < -0.99 => "-Y",
            Some(n) if n[2] > 0.99 => "+Z",
            Some(n) if n[2] < -0.99 => "-Z",
            Some(_) => "slant",
            None => "curved",
        };
        brep.set_face_label(
            face.index,
            &Label::new(1, role).with_index(face.index as u32),
        )?;
    }
    Ok(brep)
}

fn facing(direction: [f64; 3]) -> FaceQuery {
    FaceQuery::Facing { direction }
}

fn find(brep: &Brep, query: &FaceQuery) -> Vec<usize> {
    brep.find_faces(query, &Tolerance::default()).unwrap()
}

/// What a move gives back, as plain data.
#[derive(Debug)]
struct Result_ {
    manifold: bool,
    volume: f64,
    faces: usize,
    edges: usize,
    labels: Vec<Option<Label>>,
    kinds: Vec<&'static str>,
    roles: Vec<MovedFace>,
    surfaces: Vec<Surface>,
    mesh: Mesh,
}

struct Case {
    before: Result_,
    after: Result_,
}

fn measure(brep: &Brep, roles: Vec<MovedFace>) -> Result_ {
    let faces = brep.faces().unwrap();
    Result_ {
        manifold: brep.is_manifold().unwrap(),
        volume: brep.volume(PRECISE).unwrap(),
        faces: brep.face_count(),
        edges: brep.edge_count(),
        labels: faces.iter().map(|f| f.label.clone()).collect(),
        kinds: faces.iter().map(|f| kind(&f.surface)).collect(),
        roles,
        surfaces: faces.iter().map(|f| f.surface.clone()).collect(),
        mesh: brep.tessellate(&TessellationParams::default()).unwrap(),
    }
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

/// Builds a part, moves the faces `pick` finds by `motion`, and reports both.
fn moved(
    part: fn(&Kernel) -> Result<Brep, Error>,
    pick: impl FnOnce(&Brep) -> Vec<usize> + Send + 'static,
    motion: RigidMotion,
) -> Case {
    on_kernel(move |k| {
        let brep = part(k).unwrap();
        let faces = pick(&brep);
        let before = measure(&brep, vec![MovedFace::Kept; brep.face_count()]);
        let moved = brep.move_faces(&faces, &motion).expect("the move");
        let after = measure(&moved.brep, moved.faces);
        Case { before, after }
    })
}

/// The move's error.
fn refused(
    part: fn(&Kernel) -> Result<Brep, Error>,
    pick: impl FnOnce(&Brep) -> Vec<usize> + Send + 'static,
    motion: RigidMotion,
) -> Error {
    on_kernel(move |k| {
        let brep = part(k).unwrap();
        let faces = pick(&brep);
        brep.move_faces(&faces, &motion).expect_err("a refusal")
    })
}

fn reason(error: &Error) -> &'static str {
    match error {
        Error::InvalidArgument { reason, .. } => reason,
        other => panic!("not a refusal in words: {other:?}"),
    }
}

fn assert_close(name: &str, actual: f64, expected: f64) {
    let error = (actual - expected).abs() / expected.abs();
    println!("{name}: volume {actual:.6} against {expected:.6} ({error:.1e})");
    assert!(
        error < 1e-7,
        "{name}: volume {actual} is {error:e} off {expected}"
    );
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

/// The same faces and edges, labels kept face by face, a manifold solid,
/// and a mesh that encloses its volume, winds outward and is closed.
fn assert_sound(name: &str, case: &Case) {
    let (before, after) = (&case.before, &case.after);
    assert!(after.manifold, "{name}: not manifold");
    assert_eq!(after.faces, before.faces, "{name}: faces");
    assert_eq!(after.edges, before.edges, "{name}: edges");
    assert_eq!(after.labels, before.labels, "{name}: labels by face");
    let mesh = &after.mesh;
    assert_eq!(mesh.faces.len(), after.faces, "{name}: mesh faces");
    let mut enclosed = 0.0;
    let mut against = 0;
    for tri in mesh.triangles() {
        let [a, b, c] = tri.map(|i| mesh.positions[i as usize].map(f64::from));
        let normal = cross(
            [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
            [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
        );
        let kernel = tri.iter().fold([0.0; 3], |acc, &i| {
            let n = mesh.normals[i as usize].map(f64::from);
            [acc[0] + n[0], acc[1] + n[1], acc[2] + n[2]]
        });
        if dot(normal, kernel) <= 0.0 {
            against += 1;
        }
        enclosed += dot(a, cross(b, c)) / 6.0;
    }
    assert_eq!(against, 0, "{name}: {against} triangles wind inward");
    let error = (enclosed - after.volume).abs() / after.volume;
    assert!(
        error < 5e-3,
        "{name}: the mesh encloses {enclosed}, the solid is {}",
        after.volume
    );
    for (i, edge) in mesh.edges.iter().enumerate() {
        assert_eq!(edge.use_count, 2, "{name}: edge {i} isn't manifold");
    }
}

fn top(brep: &Brep) -> Vec<usize> {
    find(brep, &facing([0.0, 0.0, 1.0]))
}

fn right(brep: &Brep) -> Vec<usize> {
    find(brep, &facing([1.0, 0.0, 0.0]))
}

fn role_of(case: &Case, label: &str) -> MovedFace {
    let index = case
        .after
        .labels
        .iter()
        .position(|l| l.as_ref().is_some_and(|l| l.role == label))
        .unwrap_or_else(|| panic!("no face labelled {label}"));
    case.after.roles[index]
}

#[test]
fn a_face_moved_along_its_normal() {
    for z in [12.0, -15.0] {
        let case = moved(block, top, RigidMotion::translation([0.0, 0.0, z]));
        let name = format!("top moved {z:+}");
        assert_sound(&name, &case);
        assert_close(&name, case.after.volume, L * W * (H + z));
        assert_eq!(role_of(&case, "+Z"), MovedFace::Moved);
        for wall in ["+X", "-X", "+Y", "-Y"] {
            assert_eq!(role_of(&case, wall), MovedFace::Followed, "{name}: {wall}");
        }
        assert_eq!(role_of(&case, "-Z"), MovedFace::Kept, "{name}: bottom");
        assert!(case.after.kinds.iter().all(|k| *k == "plane"));
    }
}

#[test]
fn a_face_moved_across_its_normal() {
    // A plane moved within itself stays where it is: only the part of the
    // move along its normal shows.
    let case = moved(block, top, RigidMotion::translation([7.0, -4.0, 10.0]));
    assert_sound("top moved across", &case);
    assert_close("top moved across", case.after.volume, L * W * (H + 10.0));
    // The right face moved straight across its normal changes nothing.
    let case = moved(block, right, RigidMotion::translation([0.0, 9.0, -3.0]));
    assert_sound("right moved within its plane", &case);
    assert_close("right moved within its plane", case.after.volume, L * W * H);
    // A top with a boss standing on it, moved up and sideways, is moved up:
    // the boss 3 shorter. Its plane's own sideways shift once carried the
    // trims round the boss's foot off it.
    let boss = std::f64::consts::PI * 100.0;
    let lid = |brep: &Brep| {
        find(
            brep,
            &FaceQuery::OnPlane {
                point: [0.0, 0.0, H / 2.0],
                normal: [0.0, 0.0, 1.0],
            },
        )
    };
    for translation in [[0.0, 0.0, 3.0], [7.0, -4.0, 3.0]] {
        let case = moved(bossed, lid, RigidMotion::translation(translation));
        let name = format!("bossed top moved {translation:?}");
        assert_sound(&name, &case);
        assert_close(&name, case.after.volume, L * W * (H + 3.0) + boss * 12.0);
    }
}

/// The box with a boss 10 mm in radius and 15 tall standing on its top at
/// (20, 5), labelled.
fn bossed(k: &Kernel) -> Result<Brep, Error> {
    let boss = Brep::cylinder(k, [20.0, 5.0, H / 2.0], 10.0, 15.0)?;
    labelled(Brep::cuboid(k, [-L / 2.0, -W / 2.0, -H / 2.0], [L, W, H])?.union(boss)?)
}

/// The boss's side and top.
fn boss_faces(brep: &Brep) -> Vec<usize> {
    let mut faces = find(brep, &FaceQuery::Kind(SurfaceKind::Cylinder));
    faces.extend(find(
        brep,
        &FaceQuery::OnPlane {
            point: [0.0, 0.0, H / 2.0 + 15.0],
            normal: [0.0, 0.0, 1.0],
        },
    ));
    faces
}

#[test]
fn a_feature_slid_across_the_face_it_stands_on() {
    let case = moved(
        bossed,
        boss_faces,
        RigidMotion::translation([-15.0, 8.0, 0.0]),
    );
    assert_sound("boss slid", &case);
    let boss = std::f64::consts::PI * 100.0 * 15.0;
    assert_close("boss slid", case.after.volume, L * W * H + boss);
    let axis = case
        .after
        .surfaces
        .iter()
        .find_map(|s| match s {
            Surface::Cylinder { origin, .. } => Some(*origin),
            _ => None,
        })
        .unwrap();
    assert!(
        (axis[0] - 5.0).abs() < 1e-9 && (axis[1] - 13.0).abs() < 1e-9,
        "the boss's axis is at {axis:?}"
    );
    // Lifted along with it, the boss grows taller.
    let case = moved(
        bossed,
        boss_faces,
        RigidMotion::translation([-15.0, 8.0, 6.0]),
    );
    assert_sound("boss lifted", &case);
    let lifted = std::f64::consts::PI * 100.0 * 21.0;
    assert_close("boss lifted", case.after.volume, L * W * H + lifted);
}

/// The right face turned `degrees` about its bottom edge, the line
/// x = L/2, z = -H/2 along +Y: positive leans it out.
fn tilt_right(degrees: f64) -> RigidMotion {
    RigidMotion::turn([L / 2.0, 0.0, -H / 2.0], [0.0, 1.0, 0.0], degrees)
}

#[test]
fn a_face_tilted_about_an_edge() {
    for degrees in [15.0, -15.0, 30.0] {
        let name = format!("right face tilted {degrees}°");
        let case = moved(block, right, tilt_right(degrees));
        assert_sound(&name, &case);
        // A trapezoid in XZ: the top edge moves out by H tan(angle).
        let lean = H * degrees.to_radians().tan();
        assert_close(&name, case.after.volume, W * H * (L + lean / 2.0));
        assert_eq!(role_of(&case, "+X"), MovedFace::Moved);
        assert_eq!(role_of(&case, "-X"), MovedFace::Kept);
        let normal = case
            .after
            .surfaces
            .iter()
            .zip(&case.after.labels)
            .find_map(|(s, l)| match (s, l) {
                (Surface::Plane { normal, .. }, Some(l)) if l.role == "+X" => Some(*normal),
                _ => None,
            })
            .unwrap();
        let tilt = normal[2].abs().asin().to_degrees();
        assert!(
            (tilt - degrees.abs()).abs() < 1e-9,
            "{name}: the face leans {tilt}°"
        );
        assert!(case.after.kinds.iter().all(|k| *k == "plane"));
    }
}

#[test]
fn every_wall_turned_in_makes_a_frustum() {
    // Each wall turned 10° in about its bottom edge, one after another, is a
    // draft: angled facets on every side.
    let case = on_kernel(|k| {
        let mut brep = block(k).unwrap();
        let before = measure(&brep, vec![MovedFace::Kept; brep.face_count()]);
        let walls = [
            ([1.0, 0.0, 0.0], [L / 2.0, 0.0, -H / 2.0], [0.0, -1.0, 0.0]),
            ([-1.0, 0.0, 0.0], [-L / 2.0, 0.0, -H / 2.0], [0.0, 1.0, 0.0]),
            ([0.0, 1.0, 0.0], [0.0, W / 2.0, -H / 2.0], [1.0, 0.0, 0.0]),
            (
                [0.0, -1.0, 0.0],
                [0.0, -W / 2.0, -H / 2.0],
                [-1.0, 0.0, 0.0],
            ),
        ];
        let mut roles = Vec::new();
        for (normal, point, axis) in walls {
            let face = find(&brep, &facing(normal));
            let moved = brep
                .move_faces(&face, &RigidMotion::turn(point, axis, 10.0))
                .unwrap();
            brep = moved.brep;
            roles = moved.faces;
        }
        let after = measure(&brep, roles);
        Case { before, after }
    });
    assert_sound("frustum", &case);
    let t = 10f64.to_radians().tan();
    let frustum = L * W * H - (L + W) * t * H * H + 4.0 / 3.0 * t * t * H.powi(3);
    assert_close("frustum", case.after.volume, frustum);
}

/// The box with its four vertical edges rounded, 8 mm.
fn rounded(k: &Kernel) -> Result<Brep, Error> {
    let vertical = EdgeQuery::ParallelTo {
        direction: [0.0, 0.0, 1.0],
    };
    labelled(
        Brep::cuboid(k, [-L / 2.0, -W / 2.0, -H / 2.0], [L, W, H])?.fillet(
            &[EdgeBlend::round(vertical, 8.0)],
            &FilletOptions::default(),
        )?,
    )
}

fn rounded_volume(l: f64) -> f64 {
    (l * W - (4.0 - std::f64::consts::PI) * 64.0) * H
}

#[test]
fn a_face_beside_a_round() {
    // The rounds between the right face and the front and back slide out
    // with it, staying tangent to both; the front and back stretch.
    let case = moved(rounded, right, RigidMotion::translation([10.0, 0.0, 0.0]));
    assert_sound("rounded right face moved out", &case);
    assert_close(
        "rounded right face moved out",
        case.after.volume,
        rounded_volume(L + 10.0),
    );
    let carried = case
        .after
        .roles
        .iter()
        .filter(|r| **r == MovedFace::Carried)
        .count();
    assert_eq!(carried, 2, "the two rounds on the right face slide");
    assert!(
        case.after
            .kinds
            .iter()
            .filter(|k| **k == "cylinder")
            .count()
            == 4
    );
    // The top, moved up, stretches the rounds standing under its edges.
    let case = moved(rounded, top, RigidMotion::translation([0.0, 0.0, 5.0]));
    assert_sound("rounded top raised", &case);
    assert_close(
        "rounded top raised",
        case.after.volume,
        rounded_volume(L) * (H + 5.0) / H,
    );
    // A turn can't carry the rounds.
    let error = refused(rounded, right, tilt_right(10.0));
    assert!(
        reason(&error).contains("a round on a turned face"),
        "{error}"
    );
}

#[test]
fn labels_are_kept() {
    let case = moved(
        block,
        right,
        tilt_right(20.0).then_translate([5.0, 0.0, 0.0]),
    );
    assert_sound("tilted and moved", &case);
    let roles: Vec<&str> = case
        .after
        .labels
        .iter()
        .map(|l| l.as_ref().unwrap().role.as_str())
        .collect();
    let before: Vec<&str> = case
        .before
        .labels
        .iter()
        .map(|l| l.as_ref().unwrap().role.as_str())
        .collect();
    assert_eq!(roles, before);
    let lean = H * 20f64.to_radians().tan();
    assert_close(
        "tilted and moved",
        case.after.volume,
        W * H * (L + 5.0 + lean / 2.0),
    );
}

#[test]
fn a_move_the_faces_round_it_cant_reach_is_refused() {
    // Down past the bottom: the walls would shrink to nothing.
    let error = refused(block, top, RigidMotion::translation([0.0, 0.0, -H - 5.0]));
    assert!(
        reason(&error).starts_with("the faces round them can't reach them"),
        "{error}"
    );
    // Turned so far in that its top edge crosses the left face.
    let error = refused(block, right, tilt_right(-75.0));
    assert!(
        reason(&error).starts_with("the faces round them can't reach them"),
        "{error}"
    );
}

#[test]
fn bad_moves_are_refused_before_the_kernel() {
    let errors = on_kernel(|k| {
        let reasons = [
            block(k)
                .unwrap()
                .move_faces(&[], &RigidMotion::translation([1.0, 0.0, 0.0])),
            block(k)
                .unwrap()
                .move_faces(&[0], &RigidMotion::translation([f64::NAN, 0.0, 0.0])),
            block(k)
                .unwrap()
                .move_faces(&[0], &RigidMotion::turn([0.0; 3], [0.0; 3], 10.0)),
            block(k).unwrap().move_faces(
                &[0],
                &RigidMotion::turn([0.0; 3], [0.0, 0.0, 1.0], f64::INFINITY),
            ),
            block(k)
                .unwrap()
                .move_faces(&[0, 0], &RigidMotion::translation([1.0, 0.0, 0.0])),
            block(k)
                .unwrap()
                .move_faces(&[6], &RigidMotion::translation([1.0, 0.0, 0.0])),
        ];
        reasons.map(|r| r.map(|_| ()).unwrap_err())
    });
    let words: Vec<&str> = errors.iter().map(reason).collect();
    assert_eq!(
        words,
        [
            "no faces to move",
            "the translation must be finite",
            "the axis direction must be finite and not zero",
            "the angle must be finite",
            "a face to move is out of range or repeated",
            "a face to move is out of range or repeated",
        ]
    );
}

#[test]
fn a_move_of_nothing_keeps_the_brep() {
    let case = moved(block, top, RigidMotion::default());
    assert_sound("no move", &case);
    assert_close("no move", case.after.volume, L * W * H);
    assert_eq!(
        case.after
            .roles
            .iter()
            .filter(|r| **r == MovedFace::Moved)
            .count(),
        1
    );
    assert!(
        case.after
            .roles
            .iter()
            .all(|r| matches!(r, MovedFace::Moved | MovedFace::Kept))
    );
}

#[test]
fn a_turned_cylinder_is_exact_or_refused() {
    // The kernel-solids notes found the local operation leaving an open
    // solid once a turned side was a cylinder: the move must either come
    // out whole or be refused in words, never broken.
    let result = on_kernel(|k| {
        let brep = bossed(k).unwrap();
        let faces = boss_faces(&brep);
        let turn = RigidMotion::turn([20.0, 5.0, H / 2.0], [1.0, 0.0, 0.0], 8.0);
        match brep.move_faces(&faces, &turn) {
            Ok(moved) => Ok((
                moved.brep.is_manifold().unwrap(),
                moved.brep.volume(PRECISE).unwrap(),
            )),
            Err(error) => Err(error),
        }
    });
    match result {
        Ok((manifold, volume)) => {
            println!("turned boss: manifold {manifold}, volume {volume}");
            assert!(manifold);
        }
        Err(error) => {
            println!("turned boss refused: {error}");
            assert!(
                reason(&error).starts_with("the faces round them can't reach them"),
                "{error}"
            );
        }
    }
}
