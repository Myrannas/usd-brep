//! Imprints on the real kernel: a tool's boundary split into a solid's faces
//! as edges, and taken back off with a keep-list merge.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{Brep, Imprint, Kernel, Label, TessellationParams};

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
const TOOL: u32 = 2;

/// A 100 × 60 × 40 block centred on the origin, each face labelled by its
/// index, so faces split later are told apart from faces made later.
fn block(k: &Kernel) -> Brep {
    let mut solid = Brep::cuboid(k, [-50.0, -30.0, -20.0], [100.0, 60.0, 40.0]).unwrap();
    for face in solid.faces().unwrap() {
        solid
            .set_face_label(
                face.index,
                &Label::new(BLOCK, "side").with_index(face.index as u32),
            )
            .unwrap();
    }
    solid
}

/// A box from `min` to `max`, each face labelled `wall` by its index.
fn tool(k: &Kernel, min: [f64; 3], max: [f64; 3]) -> Brep {
    let size = [0, 1, 2].map(|i| max[i] - min[i]);
    let mut solid = Brep::cuboid(k, min, size).unwrap();
    for face in solid.faces().unwrap() {
        solid
            .set_face_label(
                face.index,
                &Label::new(TOOL, "wall").with_index(face.index as u32),
            )
            .unwrap();
    }
    solid
}

/// A 20 × 10 square tool standing up through the block's top face only:
/// its bottom lies inside the block, so only its walls meet the block.
fn through_top(k: &Kernel) -> Brep {
    tool(k, [-10.0, -5.0, 0.0], [10.0, 5.0, 60.0])
}

fn top_face_pieces(brep: &Brep) -> Vec<usize> {
    brep.faces()
        .unwrap()
        .into_iter()
        .filter(|f| f.normal.is_some_and(|n| n[2] > 0.999) && f.bounds[1][2] > 19.999)
        .map(|f| f.index)
        .collect()
}

#[test]
fn a_tool_splits_the_faces_it_crosses_and_changes_nothing_else() {
    on_kernel(|k| {
        let before = block(k);
        let volume = before.volume(1e-9).unwrap();
        let Imprint {
            brep,
            edges,
            inside,
        } = before.imprint(through_top(k)).unwrap();

        assert!(brep.is_manifold().unwrap());
        assert!((brep.volume(1e-9).unwrap() - volume).abs() < 1e-6 * volume);
        assert!(inside.is_empty(), "a plain imprint doesn't classify");
        // The top is split into the square and the rest round it; the square's
        // four sides are new edges.
        assert_eq!(brep.face_count(), 7);
        assert_eq!(brep.edge_count(), 16);
        let top = top_face_pieces(&brep);
        assert_eq!(top.len(), 2);
        let faces = brep.faces().unwrap();
        for &piece in &top {
            assert_eq!(
                faces[piece].label, faces[top[0]].label,
                "both pieces keep the top's label"
            );
        }

        // Each new edge lies on one of the tool's walls, and says which.
        assert_eq!(edges.len(), 4);
        let all = brep.edges().unwrap();
        let mut walls = Vec::new();
        for edge in &edges {
            let geometry = &all[edge.edge];
            assert!(edge.new);
            assert_eq!(geometry.label, None, "new edges have no label");
            assert!(
                geometry
                    .faces
                    .iter()
                    .all(|f| f.is_some_and(|f| top.contains(&(f as usize))))
            );
            assert_eq!(edge.tool_faces.len(), 1);
            assert_eq!(edge.tool_faces[0].feature, TOOL);
            walls.push(edge.tool_faces[0].index);
        }
        walls.sort();
        walls.dedup();
        assert_eq!(walls.len(), 4, "one edge per wall");
    });
}

#[test]
fn classifying_says_which_faces_lie_inside_the_tool() {
    on_kernel(|k| {
        let Imprint {
            brep,
            edges,
            inside,
        } = block(k).imprint_classified(through_top(k)).unwrap();
        assert!(brep.is_manifold().unwrap());
        assert_eq!(
            (brep.face_count(), brep.edge_count(), edges.len()),
            (7, 16, 4)
        );
        assert_eq!(inside.len(), 7);
        let faces = brep.faces().unwrap();
        let inner: Vec<usize> = (0..7).filter(|&f| inside[f]).collect();
        assert_eq!(inner.len(), 1, "only the square on the top is inside");
        let [min, max] = faces[inner[0]].bounds;
        for (axis, (low, high)) in [(-10.0, 10.0), (-5.0, 5.0), (20.0, 20.0)]
            .into_iter()
            .enumerate()
        {
            assert!((min[axis] - low).abs() < 1e-4 && (max[axis] - high).abs() < 1e-4);
        }
    });
}

#[test]
fn a_tool_through_the_part_imprints_both_sides_and_its_edges_split_the_side_faces() {
    on_kernel(|k| {
        // A slab across the whole block, from y = -5 up: its wall at y = -5
        // runs round the block, its far side outside it.
        let slab = tool(k, [-60.0, -5.0, -30.0], [60.0, 40.0, 30.0]);
        let Imprint { brep, edges, .. } = block(k).imprint_classified(slab).unwrap();
        assert!(brep.is_manifold().unwrap());
        // Top, bottom, +X and -X are each split in two.
        assert_eq!(brep.face_count(), 10);
        assert_eq!(edges.len(), 4);
        for edge in &edges {
            assert_eq!(edge.tool_faces.len(), 1);
        }
        // The old edges the wall crosses are split at new vertices: 12 + 4
        // split + 4 new.
        assert_eq!(brep.edge_count(), 20);
        assert!(edges.iter().all(|edge| edge.new));
    });
}

#[test]
fn a_tool_face_along_a_face_finds_its_edges_and_splits_nothing() {
    on_kernel(|k| {
        // Outside the block, its wall at y = 30 on the block's +Y face, which
        // it covers; and the block inside one whose wall lies on that face.
        for slab in [
            tool(k, [-60.0, 30.0, -30.0], [60.0, 60.0, 30.0]),
            tool(k, [-60.0, -40.0, -30.0], [60.0, 30.0, 30.0]),
        ] {
            let Imprint {
                brep,
                edges,
                inside,
            } = block(k).imprint_classified(slab).unwrap();
            assert!(brep.is_manifold().unwrap());
            assert_eq!((brep.face_count(), brep.edge_count()), (6, 12));
            // The +Y face's four edges lie along the slab, and were there.
            assert_eq!(edges.len(), 4);
            assert!(
                edges
                    .iter()
                    .all(|edge| !edge.new && edge.tool_faces.len() == 1)
            );
            let all = brep.edges().unwrap();
            for edge in &edges {
                assert!(
                    all[edge.edge]
                        .samples
                        .iter()
                        .all(|p| (p[1] - 30.0).abs() < 1e-9)
                );
            }
            assert!(inside.iter().all(|&i| i) || inside.iter().all(|&i| !i));
        }
    });
}

/// The faces come back as they were, but for their tessellation keys: the
/// merge extends the surface of one piece over the others, so the merged
/// face and every face sharing an edge with it key differently.
#[test]
fn taking_the_edges_back_off_restores_the_faces() {
    on_kernel(|k| {
        let params = TessellationParams::default();
        let original = block(k);
        let keys = original.tessellation_keys(&params).unwrap();
        let mesh = original.tessellate(&params).unwrap();
        let Imprint {
            mut brep, edges, ..
        } = original
            .try_clone()
            .unwrap()
            .imprint(through_top(k))
            .unwrap();
        let imprinted: Vec<usize> = edges.iter().map(|e| e.edge).collect();
        let keep: Vec<usize> = (0..brep.edge_count())
            .filter(|e| !imprinted.contains(e))
            .collect();
        brep.merge_faces_keeping(&keep).unwrap();

        assert!(brep.is_manifold().unwrap());
        assert_eq!((brep.face_count(), brep.edge_count()), (6, 12));
        let faces = |b: &Brep| {
            let faces = b.faces().unwrap().into_iter();
            faces
                .map(|f| (f.label, f.surface, f.bounds))
                .collect::<Vec<_>>()
        };
        assert_eq!(faces(&brep), faces(&original));
        assert_eq!(brep.tessellate(&params).unwrap(), mesh);
        let now = brep.tessellation_keys(&params).unwrap();
        let top = top_face_pieces(&brep)[0];
        let bottom = (0..6).find(|&f| {
            f != top
                && !brep.edges().unwrap().iter().any(|e| {
                    e.faces.contains(&Some(f as u32)) && e.faces.contains(&Some(top as u32))
                })
        });
        assert_eq!(
            bottom.map(|f| now[f] == keys[f]),
            Some(true),
            "a face apart from the merge keys alike"
        );
        assert_ne!(now[top], keys[top]);
    });
}

#[test]
fn kept_edges_stay_through_a_merge() {
    on_kernel(|k| {
        let Imprint {
            mut brep, edges, ..
        } = block(k).imprint(through_top(k)).unwrap();
        // Keep one imprinted edge: the faces either side are still one face
        // round it, so the other three go but that one stays.
        let kept = edges[0].edge;
        let point = brep.edges().unwrap()[kept].samples[0];
        brep.merge_faces_keeping(&[kept]).unwrap();
        assert!(brep.is_manifold().unwrap());
        assert_eq!(top_face_pieces(&brep).len(), 1);
        let near = |p: [f64; 3]| (0..3).all(|i| (p[i] - point[i]).abs() < 1e-9);
        assert!(
            brep.edges()
                .unwrap()
                .iter()
                .any(|e| e.samples.iter().copied().any(near))
        );

        // An unrestricted merge takes every imprint off.
        let Imprint { mut brep, .. } = block(k).imprint(through_top(k)).unwrap();
        brep.merge_faces().unwrap();
        assert_eq!((brep.face_count(), brep.edge_count()), (6, 12));
    });
}

#[test]
fn bad_keep_lists_are_refused() {
    on_kernel(|k| {
        let mut solid = block(k);
        assert!(solid.merge_faces_keeping(&[12]).is_err());
        assert!(solid.merge_faces_keeping(&[1, 1]).is_err());
        assert!(solid.merge_faces_keeping(&[usize::MAX]).is_err());
    });
}
