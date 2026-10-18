//! Face-by-face tessellation (`Brep::begin_tessellation`): each face comes
//! out as a whole-solid tessellation gives it, whichever other faces are
//! triangulated, and a face's key fixes its triangles, so triangles kept by
//! key from one solid stand in exactly for a face of the next.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::collections::HashMap;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, EdgeBlend, EdgeQuery, Error, FaceKey, FilletOptions, Kernel, Mesh, SurfaceCaches,
    SurfaceKind, TessellationParams,
};

const L: f64 = 100.0;
const W: f64 = 60.0;
const H: f64 = 40.0;

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

fn block(k: &Kernel) -> Result<Brep, Error> {
    Brep::cuboid(k, [-L / 2.0, -W / 2.0, -H / 2.0], [L, W, H])
}

/// A box with rounded vertical edges and a hole at (x, 0) whose top rim is
/// rounded too: planes, cylinders, tori and their seams.
fn holed(k: &Kernel, x: f64, r: f64) -> Result<Brep, Error> {
    let options = FilletOptions::default();
    let mut vertical = Vec::new();
    for sx in [1.0, -1.0] {
        for sy in [1.0, -1.0] {
            vertical.push(EdgeQuery::BetweenFaces {
                a: [sx, 0.0, 0.0],
                b: [0.0, sy, 0.0],
            });
        }
    }
    let rounded = block(k)?.fillet(&[EdgeBlend::round(EdgeQuery::Or(vertical), 6.0)], &options)?;
    let drilled = rounded.difference(Brep::cylinder(k, [x, 0.0, -H / 2.0], r, H)?)?;
    let rim = EdgeQuery::OnCircle {
        center: [x, 0.0, H / 2.0],
        axis: [0.0, 0.0, 1.0],
        radius: r,
    };
    drilled.fillet(&[EdgeBlend::round(rim, 2.0)], &options)
}

/// The box cut down by a sphere: spherical patches meeting planes.
fn rounded_block(k: &Kernel, radius: f64) -> Result<Brep, Error> {
    block(k)?.intersection(Brep::sphere(k, [0.0; 3], radius)?)
}

/// A boss on the box with every sharp edge rounded.
fn bossed(k: &Kernel, boss_radius: f64) -> Result<Brep, Error> {
    let part = block(k)?.union(Brep::cylinder(k, [0.0, 0.0, H / 2.0], boss_radius, 15.0)?)?;
    part.fillet(
        &[EdgeBlend::round(
            EdgeQuery::Sharp {
                min_angle_deg: 30.0,
            },
            1.5,
        )],
        &FilletOptions::default(),
    )
}

/// One face's share of a mesh: its vertices in the order its triangles first
/// use them, and its triangles indexing those.
#[derive(Clone, Debug, PartialEq)]
struct Piece {
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    uvs: Vec<[f32; 2]>,
    triangles: Vec<[u32; 3]>,
}

/// Every face's piece, numbering each face's vertices afresh, so meshes
/// that order faces differently compare.
fn pieces(mesh: &Mesh) -> Vec<Piece> {
    let mut pieces = vec![
        Piece {
            positions: Vec::new(),
            normals: Vec::new(),
            uvs: Vec::new(),
            triangles: Vec::new()
        };
        mesh.faces.len()
    ];
    let mut local: Vec<HashMap<u32, u32>> = vec![HashMap::new(); mesh.faces.len()];
    for (triangle, &face) in mesh.triangles().zip(&mesh.face_ids) {
        let piece = &mut pieces[face as usize];
        let map = &mut local[face as usize];
        let corners = triangle.map(|v| {
            *map.entry(v).or_insert_with(|| {
                let i = v as usize;
                piece.positions.push(mesh.positions[i]);
                piece.normals.push(mesh.normals[i]);
                piece.uvs.push(mesh.surface_uvs[i]);
                (piece.positions.len() - 1) as u32
            })
        });
        piece.triangles.push(corners);
    }
    pieces
}

/// Checks a face-by-face mesh's layout: faces in order, each with one run of
/// triangles and one run of vertices of its own, every vertex used.
fn assert_face_runs(name: &str, mesh: &Mesh) {
    let mut next_vertex = 0u32;
    let mut last_face = None;
    let mut first = 0usize;
    let face_ids = &mesh.face_ids;
    while first < face_ids.len() {
        let face = face_ids[first];
        assert!(
            last_face.is_none_or(|last| face > last),
            "{name}: face {face} out of order"
        );
        last_face = Some(face);
        let end = first + face_ids[first..].iter().take_while(|&&f| f == face).count();
        let used: Vec<u32> = mesh.indices[first * 3..end * 3].to_vec();
        let low = *used.iter().min().unwrap();
        let high = *used.iter().max().unwrap();
        assert_eq!(
            low, next_vertex,
            "{name}: face {face}'s vertices don't follow the last face's"
        );
        let mut seen = vec![false; (high - low + 1) as usize];
        for v in used {
            seen[(v - low) as usize] = true;
        }
        assert!(
            seen.iter().all(|&s| s),
            "{name}: face {face} has unused vertices"
        );
        next_vertex = high + 1;
        first = end;
    }
    assert_eq!(
        next_vertex as usize,
        mesh.positions.len(),
        "{name}: vertices past the last face"
    );
}

fn final_params() -> TessellationParams {
    TessellationParams::default()
}

fn preview_params() -> TessellationParams {
    TessellationParams {
        chord_height: 0.5,
        ..TessellationParams::default()
    }
}

type Make = fn(&Kernel) -> Result<Brep, Error>;

fn parts() -> Vec<(&'static str, Make)> {
    vec![
        ("holed", |k| holed(k, 10.0, 9.0)),
        ("rounded block", |k| rounded_block(k, 55.0)),
        ("bossed", |k| bossed(k, 12.0)),
    ]
}

#[test]
fn face_by_face_matches_the_whole_tessellation() {
    for (name, make) in parts() {
        for params in [final_params(), preview_params()] {
            let (whole, by_face) = on_kernel(move |k| {
                let part = make(k).unwrap();
                (
                    part.tessellate(&params).unwrap(),
                    part.tessellate_faces(&params, None).unwrap(),
                )
            });
            assert_face_runs(name, &by_face);
            assert_eq!(by_face.faces, whole.faces, "{name}: faces");
            assert_eq!(by_face.edges, whole.edges, "{name}: edges");
            assert_eq!(
                by_face.triangle_count(),
                whole.triangle_count(),
                "{name}: triangles"
            );
            assert_eq!(
                by_face.positions.len(),
                whole.positions.len(),
                "{name}: vertices"
            );
            for (face, (a, b)) in pieces(&by_face).iter().zip(pieces(&whole)).enumerate() {
                assert_eq!(
                    *a, b,
                    "{name}: face {face} differs from the whole tessellation's"
                );
            }
        }
    }
}

#[test]
fn partial_meshes_hold_the_same_faces() {
    for (name, make) in parts() {
        let (full, partial, wanted) = on_kernel(move |k| {
            let part = make(k).unwrap();
            let params = final_params();
            let full = part.tessellate_faces(&params, None).unwrap();
            let wanted: Vec<bool> = (0..full.faces.len()).map(|f| f % 3 != 1).collect();
            let partial = part.tessellate_faces(&params, Some(&wanted)).unwrap();
            (full, partial, wanted)
        });
        assert_face_runs(name, &partial);
        assert_eq!(partial.faces, full.faces, "{name}: faces");
        assert_eq!(partial.edges, full.edges, "{name}: edges");
        for (face, (a, b)) in pieces(&partial).iter().zip(pieces(&full)).enumerate() {
            if wanted[face] {
                assert_eq!(*a, b, "{name}: face {face}");
            } else {
                assert!(
                    a.triangles.is_empty() && a.positions.is_empty(),
                    "{name}: face {face} wasn't wanted"
                );
            }
        }
    }
}

#[test]
fn keys_repeat() {
    for (name, make) in parts() {
        let (first, second, third) = on_kernel(move |k| {
            let part = make(k).unwrap();
            let again = make(k).unwrap();
            let params = final_params();
            let keys = |brep: &Brep| {
                brep.begin_tessellation(&params)
                    .unwrap()
                    .face_keys()
                    .unwrap()
            };
            (keys(&part), keys(&part), keys(&again))
        });
        assert!(
            first.iter().all(Option::is_some),
            "{name}: every face has a key"
        );
        assert_eq!(first, second, "{name}: keys of one solid");
        assert_eq!(first, third, "{name}: keys of the same solid built twice");
    }
}

#[test]
fn keys_follow_the_parameters() {
    let (fine, coarse) = on_kernel(|k| {
        let part = holed(k, 10.0, 9.0).unwrap();
        let keys = |params: &TessellationParams| {
            part.begin_tessellation(params)
                .unwrap()
                .face_keys()
                .unwrap()
        };
        (keys(&final_params()), keys(&preview_params()))
    });
    for (face, (a, b)) in fine.iter().zip(&coarse).enumerate() {
        assert_ne!(
            a, b,
            "face {face}: a coarser chord height changes every key"
        );
    }
}

/// A face's key from one solid, with the triangles it was given there.
type Kept = HashMap<FaceKey, Piece>;

/// Tessellates `part` reusing the pieces `kept` has keys for, as a caller
/// keeping triangles would, and checks the result against a whole
/// tessellation. Returns how many faces were reused, and adds the new ones.
fn reuse(name: &str, part: &Brep, params: &TessellationParams, kept: &mut Kept) -> usize {
    let mut tessellation = part.begin_tessellation(params).unwrap();
    let keys = tessellation.face_keys().unwrap();
    let wanted: Vec<bool> = keys
        .iter()
        .map(|key| key.as_ref().is_none_or(|key| !kept.contains_key(key)))
        .collect();
    let fresh = pieces(&tessellation.mesh(Some(&wanted)).unwrap());
    let whole = pieces(&part.tessellate(params).unwrap());
    let mut reused = 0;
    for (face, key) in keys.into_iter().enumerate() {
        let piece = if wanted[face] {
            fresh[face].clone()
        } else {
            reused += 1;
            kept[key.as_ref().unwrap()].clone()
        };
        assert_eq!(
            piece, whole[face],
            "{name}: face {face} (reused: {})",
            !wanted[face]
        );
        if let Some(key) = key {
            kept.entry(key).or_insert(piece);
        }
    }
    reused
}

#[test]
fn kept_faces_stand_in_exactly_across_edits() {
    // A hole moving and growing, as a drag would move it, at both qualities.
    let steps = [
        (10.0, 9.0),
        (12.0, 9.0),
        (14.0, 9.0),
        (14.0, 10.0),
        (10.0, 9.0),
        (-20.0, 9.0),
    ];
    let reused = on_kernel(move |k| {
        let mut kept = Kept::new();
        let mut reused = Vec::new();
        for (x, r) in steps {
            let part = holed(k, x, r).unwrap();
            for params in [final_params(), preview_params()] {
                reused.push(reuse(
                    &format!("hole at {x}, radius {r}"),
                    &part,
                    &params,
                    &mut kept,
                ));
            }
        }
        reused
    });
    // The first step finds nothing kept; every later one keeps the faces
    // away from the hole (the sides and the bottom's rounded edges).
    assert_eq!(&reused[..2], &[0, 0]);
    assert!(
        reused[2..].iter().all(|&n| n >= 8),
        "faces reused per step: {reused:?}"
    );
    // Going back to the first position reuses every face.
    let faces = on_kernel(|k| holed(k, 10.0, 9.0).unwrap().face_count());
    assert_eq!(
        reused[8], faces,
        "faces reused when the hole comes back: {reused:?}"
    );
}

#[test]
fn kept_faces_stand_in_across_other_parts() {
    let reused = on_kernel(|k| {
        let mut kept = Kept::new();
        let mut reused = Vec::new();
        for (name, part) in [
            ("boss 12", bossed(k, 12.0)),
            ("boss 14", bossed(k, 14.0)),
            ("sphere 55", rounded_block(k, 55.0)),
            ("sphere 57", rounded_block(k, 57.0)),
        ] {
            reused.push(reuse(name, &part.unwrap(), &final_params(), &mut kept));
        }
        reused
    });
    // A wider boss keeps the box's sides and bottom.
    assert!(reused[1] >= 5, "faces reused: {reused:?}");
}

#[test]
fn rejects_a_flag_count_other_than_the_face_count() {
    let (faces, too_few, too_many) = on_kernel(|k| {
        let part = block(k).unwrap();
        let params = final_params();
        let faces = part.begin_tessellation(&params).unwrap().face_count();
        let too_few = part.tessellate_faces(&params, Some(&[true; 5])).map(|_| ());
        let too_many = part.tessellate_faces(&params, Some(&[true; 7])).map(|_| ());
        (faces, too_few, too_many)
    });
    assert_eq!(faces, 6);
    assert!(too_few.is_err() && too_many.is_err());
}

/// Tessellates `part` with [`Brep::tessellate_some`], reusing the pieces
/// `kept` has [`Brep::tessellation_keys`] keys for, and checks every face
/// against a whole tessellation. Returns how many faces were reused, and adds
/// the new ones.
fn reuse_some(name: &str, part: &Brep, params: &TessellationParams, kept: &mut Kept) -> usize {
    reuse_some_with(name, part, params, kept, None)
}

/// [`reuse_some`], keeping surface caches in `caches`.
fn reuse_some_with(
    name: &str,
    part: &Brep,
    params: &TessellationParams,
    kept: &mut Kept,
    caches: Option<&mut SurfaceCaches>,
) -> usize {
    let keys = part.tessellation_keys(params).unwrap();
    let wanted: Vec<bool> = keys
        .iter()
        .map(|key| key.as_ref().is_none_or(|key| !kept.contains_key(key)))
        .collect();
    let mesh = part
        .tessellate_some(params, &wanted, caches)
        .unwrap()
        .unwrap_or_else(|| panic!("{name}: no face of these parts refines"));
    assert_face_runs(name, &mesh);
    let fresh = pieces(&mesh);
    let whole_mesh = part.tessellate(params).unwrap();
    assert_eq!(mesh.faces, whole_mesh.faces, "{name}: faces");
    assert_eq!(mesh.edges, whole_mesh.edges, "{name}: edges");
    let whole = pieces(&whole_mesh);
    let mut reused = 0;
    for (face, key) in keys.into_iter().enumerate() {
        let piece = if wanted[face] {
            fresh[face].clone()
        } else {
            reused += 1;
            assert!(
                fresh[face].triangles.is_empty(),
                "{name}: face {face} wasn't wanted"
            );
            kept[key.as_ref().unwrap()].clone()
        };
        assert_eq!(
            piece, whole[face],
            "{name}: face {face} (reused: {})",
            !wanted[face]
        );
        if let Some(key) = key {
            kept.entry(key).or_insert(piece);
        }
    }
    reused
}

#[test]
fn some_faces_match_the_whole_tessellation() {
    for (name, make) in parts() {
        let reused = on_kernel(move |k| {
            let part = make(k).unwrap();
            let mut kept = Kept::new();
            let mut reused = Vec::new();
            for params in [final_params(), preview_params()] {
                // Every face, then every face again from what was kept.
                reused.push(reuse_some(name, &part, &params, &mut kept));
                reused.push(reuse_some(name, &part, &params, &mut kept));
            }
            (reused, part.face_count())
        });
        let (reused, faces) = reused;
        assert_eq!(reused, [0, faces, 0, faces], "{name}: faces reused");
    }
}

#[test]
fn some_faces_alone_match_the_whole_tessellation() {
    for (name, make) in parts() {
        on_kernel(move |k| {
            let part = make(k).unwrap();
            let params = final_params();
            let keys = part.tessellation_keys(&params).unwrap();
            let mut kept = Kept::new();
            reuse_some(name, &part, &params, &mut kept);
            // Drop every third face, so the rest are triangulated alone.
            for (face, key) in keys.iter().enumerate() {
                if face % 3 == 1 {
                    kept.remove(key.as_ref().unwrap());
                }
            }
            let dropped = (keys.len() + 1) / 3;
            let reused = reuse_some(name, &part, &params, &mut kept);
            assert_eq!(reused, keys.len() - dropped, "{name}: faces reused");
        });
    }
}

#[test]
fn kept_faces_stand_in_exactly_when_tessellated_alone() {
    // As kept_faces_stand_in_exactly_across_edits, but splitting edges only
    // for the faces that changed.
    let steps = [
        (10.0, 9.0),
        (12.0, 9.0),
        (14.0, 9.0),
        (14.0, 10.0),
        (10.0, 9.0),
        (-20.0, 9.0),
    ];
    let reused = on_kernel(move |k| {
        let mut kept = Kept::new();
        let mut reused = Vec::new();
        for (x, r) in steps {
            let part = holed(k, x, r).unwrap();
            for params in [final_params(), preview_params()] {
                reused.push(reuse_some(
                    &format!("hole at {x}, radius {r}"),
                    &part,
                    &params,
                    &mut kept,
                ));
            }
        }
        reused
    });
    assert_eq!(&reused[..2], &[0, 0]);
    assert!(
        reused[2..].iter().all(|&n| n >= 8),
        "faces reused per step: {reused:?}"
    );
    let faces = on_kernel(|k| holed(k, 10.0, 9.0).unwrap().face_count());
    assert_eq!(
        reused[8], faces,
        "faces reused when the hole comes back: {reused:?}"
    );
}

#[test]
fn kept_faces_stand_in_across_other_parts_when_tessellated_alone() {
    let reused = on_kernel(|k| {
        let mut kept = Kept::new();
        let mut reused = Vec::new();
        for (name, part) in [
            ("boss 12", bossed(k, 12.0)),
            ("boss 14", bossed(k, 14.0)),
            ("sphere 55", rounded_block(k, 55.0)),
            ("sphere 57", rounded_block(k, 57.0)),
            ("boss 12 again", bossed(k, 12.0)),
        ] {
            let part = part.unwrap();
            for params in [final_params(), preview_params()] {
                reused.push(reuse_some(name, &part, &params, &mut kept));
            }
        }
        reused
    });
    assert!(reused[2] >= 5, "faces reused: {reused:?}");
    let faces = on_kernel(|k| bossed(k, 12.0).unwrap().face_count());
    assert_eq!(&reused[8..], &[faces, faces], "faces reused: {reused:?}");
}

#[test]
fn tessellation_keys_repeat_and_follow_the_parameters() {
    for (name, make) in parts() {
        let (first, second, third, coarse, whole) = on_kernel(move |k| {
            let part = make(k).unwrap();
            let again = make(k).unwrap();
            let keys =
                |brep: &Brep, params: &TessellationParams| brep.tessellation_keys(params).unwrap();
            (
                keys(&part, &final_params()),
                keys(&part, &final_params()),
                keys(&again, &final_params()),
                keys(&part, &preview_params()),
                part.begin_tessellation(&final_params())
                    .unwrap()
                    .face_keys()
                    .unwrap(),
            )
        });
        assert!(
            first.iter().all(Option::is_some),
            "{name}: every face has a key"
        );
        assert_eq!(first, second, "{name}: keys of one solid");
        assert_eq!(first, third, "{name}: keys of the same solid built twice");
        for (face, key) in first.iter().enumerate() {
            assert!(
                !coarse.contains(key) && !whole.contains(key),
                "{name}: face {face}'s key is another's"
            );
        }
    }
}

#[test]
fn tessellate_some_rejects_a_flag_count_other_than_the_face_count() {
    let (too_few, too_many, none) = on_kernel(|k| {
        let part = block(k).unwrap();
        let params = final_params();
        (
            part.tessellate_some(&params, &[true; 5], None).map(|_| ()),
            part.tessellation_keys(&params).map(|keys| keys.len()),
            part.tessellate_some(&params, &[false; 6], None).unwrap(),
        )
    });
    assert!(too_few.is_err());
    assert_eq!(too_many.unwrap(), 6);
    let none = none.expect("nothing to refine");
    assert_eq!(none.faces.len(), 6);
    assert_eq!(none.triangle_count(), 0);
}

/// A disk with a hole half a millimetre from its rim, away from the rim's
/// seam. At a coarse chord height the rim's polygon cuts across the hole's,
/// so SMLib refines the top and bottom faces' boundaries.
fn crowded(k: &Kernel) -> Result<Brep, Error> {
    let disk = Brep::cylinder(k, [0.0, 0.0, 0.0], 30.0, 10.0)?;
    let (r, angle) = (20.5, 22.5f64.to_radians());
    disk.difference(Brep::cylinder(
        k,
        [r * angle.cos(), r * angle.sin(), -1.0],
        9.0,
        12.0,
    )?)
}

fn coarse_params() -> TessellationParams {
    TessellationParams {
        chord_height: 5.0,
        curve_angle_deg: 45.0,
        surface_angle_deg: 45.0,
        ..TessellationParams::default()
    }
}

#[test]
fn tessellate_some_declines_boundaries_that_need_refining() {
    let (alone, together, whole) = on_kernel(|k| {
        let part = crowded(k).unwrap();
        let params = coarse_params();
        let faces = part.face_count();
        let alone: Vec<bool> = (0..faces)
            .map(|face| {
                let wanted: Vec<bool> = (0..faces).map(|f| f == face).collect();
                part.tessellate_some(&params, &wanted, None)
                    .unwrap()
                    .is_some()
            })
            .collect();
        let together = part
            .tessellate_some(&params, &vec![true; faces], None)
            .unwrap();
        (alone, together.is_some(), part.tessellate(&params).unwrap())
    });
    // The planar ends refine; the cylinders don't.
    assert_eq!(alone, [true, false, false, true]);
    assert!(!together, "the ends refine with every face wanted too");
    assert!(whole.triangle_count() > 0);
}

#[test]
fn kept_surface_caches_leave_faces_as_they_were() {
    // The hole's steps and the other parts again, with surface caches kept
    // from step to step.
    let (reused, caches_kept) = on_kernel(|k| {
        let mut kept = Kept::new();
        let mut caches = SurfaceCaches::new(k);
        let mut reused = Vec::new();
        let mut caches_kept = Vec::new();
        let mut parts: Vec<(String, Brep)> = [(10.0, 9.0), (12.0, 9.0), (14.0, 10.0), (10.0, 9.0)]
            .into_iter()
            .map(|(x, r)| (format!("hole at {x}, radius {r}"), holed(k, x, r).unwrap()))
            .collect();
        parts.push(("boss 12".into(), bossed(k, 12.0).unwrap()));
        parts.push(("boss 14".into(), bossed(k, 14.0).unwrap()));
        for (name, part) in &parts {
            for params in [final_params(), preview_params()] {
                reused.push(reuse_some_with(
                    name,
                    part,
                    &params,
                    &mut kept,
                    Some(&mut caches),
                ));
                caches_kept.push(caches.len());
            }
        }
        (reused, caches_kept)
    });
    assert_eq!(&reused[..2], &[0, 0]);
    assert!(
        reused[2..4].iter().all(|&n| n >= 8),
        "faces reused per step: {reused:?}"
    );
    // Faces next to the moving hole and the growing boss leave their
    // caches for the next step; going back to the first hole reuses every
    // face, which leaves the store as it was.
    assert_eq!(
        caches_kept[..2],
        [0, 0],
        "nothing is kept beside every face"
    );
    assert!(
        [2, 3, 4, 5, 6, 7, 10, 11]
            .iter()
            .all(|&i| caches_kept[i] > 0),
        "caches kept per step: {caches_kept:?}"
    );
}

/// A prism with drafted walls and every edge rounded, notched by a cylinder
/// at (x, 20) through its front wall and top front round: NURBS faces (the
/// walls are planar NURBS, the rounds curved ones) that keep their surfaces
/// while the notch moves.
fn notched_prism(k: &Kernel, x: f64) -> Result<Brep, Error> {
    let base = [
        [-30.0, -20.0, 0.0],
        [30.0, -20.0, 0.0],
        [30.0, 20.0, 0.0],
        [-30.0, 20.0, 0.0],
    ];
    let rounded = Brep::tapered_prism(k, &base, 25.0, 10.0)?.fillet(
        &[EdgeBlend::round(
            EdgeQuery::Sharp {
                min_angle_deg: 30.0,
            },
            4.0,
        )],
        &FilletOptions::default(),
    )?;
    rounded.difference(Brep::cylinder(k, [x, 20.0, -1.0], 8.0, 30.0)?)
}

#[test]
fn kept_surface_caches_are_copied_for_faces_tessellated_again() {
    // Every face is tessellated at every step, so the store keeps only
    // copies of the caches of faces it tessellates: a face whose surface
    // stays takes a copy of the cache kept for it, B-spline subdivision and
    // all, and leaves the original for the next step.
    let (caches_kept, nurbs) = on_kernel(|k| {
        let mut caches = SurfaceCaches::new(k);
        let mut caches_kept = Vec::new();
        let mut nurbs = Vec::new();
        for x in [0.0, 3.0, 6.0, 9.0] {
            let part = notched_prism(k, x).unwrap();
            let faces = part.faces().unwrap();
            nurbs.push(
                faces
                    .iter()
                    .filter(|face| SurfaceKind::of(&face.surface) == SurfaceKind::Nurbs)
                    .count(),
            );
            let name = format!("notch at {x}");
            let reused = reuse_some_with(
                &name,
                &part,
                &final_params(),
                &mut Kept::new(),
                Some(&mut caches),
            );
            assert_eq!(reused, 0, "{name}: nothing kept to reuse");
            caches_kept.push(caches.len());
        }
        (caches_kept, nurbs)
    });
    assert!(nurbs.iter().all(|&n| n >= 20), "NURBS faces: {nurbs:?}");
    // The first step builds every cache and keeps none; the second keeps
    // copies of the caches it builds for surfaces the first built too; from
    // the third, faces take copies and the store keeps what they copied.
    assert_eq!(caches_kept[0], 0, "caches kept per step: {caches_kept:?}");
    assert!(
        caches_kept[1] >= 20
            && caches_kept[2] >= caches_kept[1]
            && caches_kept[3] >= caches_kept[1],
        "caches kept per step: {caches_kept:?}"
    );
}
