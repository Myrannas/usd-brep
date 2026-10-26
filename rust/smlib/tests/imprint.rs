//! Imprints on the real kernel: a tool's boundary split into a solid's faces
//! as edges, taken back off with a keep-list merge, and the pieces it leaves
//! moved with Offset Face, which steps away from the rest of their face.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{Brep, EdgeQuery, Error, Imprint, Kernel, Label, TessellationParams, Tolerance};

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

fn status_name(error: &Error) -> Option<&'static str> {
    match error {
        Error::Kernel { status, .. } => status.name(),
        _ => None,
    }
}

/// The block with the 20 × 10 square imprinted on its top, and the square.
fn panel(k: &Kernel) -> (Brep, usize) {
    let Imprint { brep, .. } = block(k).imprint(through_top(k)).unwrap();
    let square = square_of(&brep);
    (brep, square)
}

fn square_of(brep: &Brep) -> usize {
    let faces = brep.faces().unwrap();
    top_face_pieces(brep)
        .into_iter()
        .find(|&f| faces[f].bounds[1][0] < 10.1)
        .unwrap()
}

/// The block with a slab imprinted from y = -5 up and taken back off every
/// face but the top, and the top's piece from y = -5 to its +Y edge.
fn edge_piece(solid: Brep, k: &Kernel) -> (Brep, usize) {
    let slab = tool(k, [-60.0, -5.0, -30.0], [60.0, 40.0, 30.0]);
    let Imprint {
        brep: mut solid,
        edges,
        ..
    } = solid.imprint(slab).unwrap();
    let all = solid.edges().unwrap();
    let on_top =
        |f: Option<u32>| f.is_some_and(|f| top_face_pieces(&solid).contains(&(f as usize)));
    let imprinted: Vec<usize> = edges.iter().map(|e| e.edge).collect();
    let keep: Vec<usize> = (0..all.len())
        .filter(|e| !imprinted.contains(e) || all[*e].faces.iter().any(|&f| on_top(f)))
        .collect();
    solid.merge_faces_keeping(&keep).unwrap();
    let faces = solid.faces().unwrap();
    let piece = top_face_pieces(&solid)
        .into_iter()
        .find(|&f| faces[f].bounds[0][1] > -5.1)
        .unwrap();
    (solid, piece)
}

#[test]
fn a_panel_steps_out_and_in_from_the_rest_of_its_face() {
    on_kernel(|k| {
        for distance in [3.0, -3.0] {
            let (solid, square) = panel(k);
            let volume = solid.volume(1e-9).unwrap();
            let label = solid.faces().unwrap()[square].label.clone();
            let moved = solid.offset_faces(&[square], distance).unwrap();

            assert!(moved.is_manifold().unwrap());
            let change = moved.volume(1e-9).unwrap() - volume;
            assert!(
                (change - 200.0 * distance).abs() < 1e-3,
                "{distance}: {change}"
            );
            // The square moved, four walls round it: 7 - 1 + 1 + 4 faces, and
            // 16 edges + 4 along its new place + 4 up the corners.
            assert_eq!((moved.face_count(), moved.edge_count()), (11, 24));
            let faces = moved.faces().unwrap();
            let height = 20.0 + distance;
            let at = |f: &&smlib::FaceGeometry| {
                f.bounds[1][0] < 10.1
                    && (f.bounds[0][2] - height).abs() < 1e-4
                    && f.normal.is_some()
            };
            let square: Vec<_> = faces
                .iter()
                .filter(|f| f.normal.is_some_and(|n| n[2] > 0.999))
                .filter(at)
                .collect();
            assert_eq!(square.len(), 1);
            assert_eq!(square[0].label, label, "the panel keeps its label");
            let walls: Vec<_> = faces.iter().filter(|f| f.label.is_none()).collect();
            assert_eq!(walls.len(), 4);
            for wall in walls {
                assert!(wall.normal.unwrap()[2].abs() < 1e-9);
                assert!((wall.bounds[1][2] - wall.bounds[0][2] - distance.abs()).abs() < 1e-4);
            }
        }
    });
}

#[test]
fn a_piece_reaching_the_edge_steps_beside_the_imprint_and_stretches_the_sides() {
    on_kernel(|k| {
        for distance in [3.0, -3.0] {
            let (solid, piece) = edge_piece(block(k), k);
            assert_eq!(solid.face_count(), 7);
            let labels = |b: &Brep| {
                let mut labels: Vec<Label> = b
                    .faces()
                    .unwrap()
                    .into_iter()
                    .filter_map(|f| f.label)
                    .collect();
                labels.sort();
                labels
            };
            let before = labels(&solid);
            let volume = solid.volume(1e-9).unwrap();
            let moved = solid.offset_faces(&[piece], distance).unwrap();

            assert!(moved.is_manifold().unwrap());
            let change = moved.volume(1e-9).unwrap() - volume;
            assert!(
                (change - 3500.0 * distance).abs() < 1e-3,
                "{distance}: {change}"
            );
            // One wall, along the imprint; +X, -X and +Y stretch or are
            // trimmed, and stay one face each.
            assert_eq!(moved.face_count(), 8);
            assert_eq!(labels(&moved), before);
            let faces = moved.faces().unwrap();
            let walls: Vec<_> = faces.iter().filter(|f| f.label.is_none()).collect();
            assert_eq!(walls.len(), 1);
            assert!(walls[0].normal.unwrap()[1].abs() > 0.999);
        }
    });
}

#[test]
fn both_pieces_moved_together_need_no_step() {
    on_kernel(|k| {
        let (solid, _) = panel(k);
        let top = top_face_pieces(&solid);
        let volume = solid.volume(1e-9).unwrap();
        let moved = solid.offset_faces(&top, 3.0).unwrap();
        assert!(moved.is_manifold().unwrap());
        assert_eq!(moved.face_count(), 7);
        assert!((moved.volume(1e-9).unwrap() - volume - 18000.0).abs() < 1e-3);
    });
}

#[test]
fn steps_that_cant_be_built_are_refused() {
    on_kernel(|k| {
        // A pocket under the panel, up to 2 below the top: pushed in 3, the
        // panel would cut into it; pushed in 1, it clears it.
        let (solid, _) = panel(k);
        let pocket = Brep::cuboid(k, [-5.0, -3.0, -30.0], [10.0, 6.0, 48.0]).unwrap();
        let holed = solid.difference(pocket).unwrap();
        let square = square_of(&holed);
        let error = holed
            .try_clone()
            .unwrap()
            .offset_faces(&[square], -3.0)
            .unwrap_err();
        assert_eq!(status_name(&error), Some("SHIM_ERR_OBSTRUCTED"), "{error}");
        assert!(
            holed
                .offset_faces(&[square], -1.0)
                .unwrap()
                .is_manifold()
                .unwrap()
        );

        // A chamfer along the top's +Y edge, beside the piece: it doesn't
        // stand along the move.
        let mut chamfered = block(k);
        let faces = chamfered.faces().unwrap();
        let facing = |axis: usize| {
            faces
                .iter()
                .find(|f| f.normal.is_some_and(|n| n[axis] > 0.999))
                .unwrap()
                .index as u32
        };
        let (top, back) = (facing(2), facing(1));
        let edge = chamfered
            .edges()
            .unwrap()
            .into_iter()
            .find(|e| e.faces.contains(&Some(top)) && e.faces.contains(&Some(back)))
            .unwrap()
            .index;
        chamfered = chamfered
            .chamfer(&EdgeQuery::Indices(vec![edge]), 2.0, &Tolerance::default())
            .unwrap();
        let (solid, piece) = edge_piece(chamfered, k);
        let error = solid.offset_faces(&[piece], 3.0).unwrap_err();
        assert_eq!(
            status_name(&error),
            Some("SHIM_ERR_UNSUPPORTED_STEP"),
            "{error}"
        );

        // A rod's side split round its middle: the piece is curved.
        let rod = Brep::cylinder(k, [0.0, 0.0, -20.0], 10.0, 40.0).unwrap();
        let Imprint { brep: rod, .. } = rod
            .imprint(tool(k, [-20.0, -20.0, 0.0], [20.0, 20.0, 30.0]))
            .unwrap();
        let faces = rod.faces().unwrap();
        let upper = faces
            .iter()
            .find(|f| f.normal.is_none() && f.bounds[0][2] > -1e-3)
            .unwrap()
            .index;
        let error = rod.offset_faces(&[upper], 1.0).unwrap_err();
        assert_eq!(
            status_name(&error),
            Some("SHIM_ERR_UNSUPPORTED_STEP"),
            "{error}"
        );
    });
}

#[test]
fn step_walls_take_the_labels_of_the_edges_they_stand_on() {
    on_kernel(|k| {
        let Imprint {
            brep: mut solid,
            edges,
            ..
        } = block(k).imprint(through_top(k)).unwrap();
        for (n, edge) in edges.iter().enumerate() {
            solid
                .set_edge_label(edge.edge, &Label::new(3, "outline").with_index(n as u32))
                .unwrap();
        }
        let square = square_of(&solid);
        let all = solid.edges().unwrap();
        // Each new edge's label, by where it runs.
        let along: Vec<([f64; 3], Label)> = edges
            .iter()
            .map(|e| (all[e.edge].midpoint(), all[e.edge].label.clone().unwrap()))
            .collect();
        for distance in [3.0, -3.0] {
            let moved = solid.try_clone().unwrap().offset_faces(&[square], distance).unwrap();
            assert!(moved.is_manifold().unwrap());
            let walls: Vec<_> = moved
                .faces()
                .unwrap()
                .into_iter()
                .filter(|f| f.label.as_ref().is_some_and(|l| l.role == "outline"))
                .collect();
            assert_eq!(walls.len(), 4, "{distance}");
            // Each wall stands on the edge whose label it took.
            for wall in walls {
                let (at, label) = along
                    .iter()
                    .find(|(at, _)| (0..2).all(|i| at[i] >= wall.bounds[0][i] - 1e-4 && at[i] <= wall.bounds[1][i] + 1e-4))
                    .unwrap();
                assert_eq!(wall.label.as_ref(), Some(label), "{distance}: the wall at {at:?}");
                assert!(wall.normal.unwrap()[2].abs() < 1e-9);
            }
        }
        // The edges' own labels stay on them.
        let moved = solid.offset_faces(&[square], 3.0).unwrap();
        let labelled = moved.edges().unwrap().into_iter().filter(|e| e.label.as_ref().is_some_and(|l| l.role == "outline"));
        assert_eq!(labelled.count(), 4);
    });
}

#[test]
fn a_push_through_keeps_an_imprint_on_the_faces_it_merges() {
    on_kernel(|k| {
        // A crate's lid with two bosses near its front edge and a panel
        // imprinted between them, then the front pushed in under the bosses.
        let mut part = Brep::cuboid(k, [-200.0, -150.0, -125.0], [400.0, 300.0, 265.0]).unwrap();
        for x in [-165.0, 165.0] {
            part = part.union(Brep::cylinder(k, [x, -115.0, 140.0], 8.0, 6.0).unwrap()).unwrap();
        }
        let panel = tool(k, [-50.0, -40.0, 100.0], [50.0, 40.0, 200.0]);
        let Imprint { brep: mut part, edges, .. } = part.imprint(panel).unwrap();
        for edge in &edges {
            part.set_edge_label(edge.edge, &Label::new(3, "outline")).unwrap();
        }
        let lid = |b: &Brep| {
            b.faces()
                .unwrap()
                .into_iter()
                .filter(|f| f.normal.is_some_and(|n| n[2] > 0.999) && (f.bounds[0][2] - 140.0).abs() < 1e-4)
                .count()
        };
        assert_eq!(lid(&part), 2, "the lid and the panel");
        let front = part
            .find_faces(
                &smlib::FaceQuery::OnPlane {
                    point: [0.0, -150.0, 0.0],
                    normal: [0.0, -1.0, 0.0],
                },
                &Tolerance::default(),
            )
            .unwrap();
        assert_eq!(part.offset_obstacles(&front, -30.0).unwrap().len(), 2);
        let volume = part.volume(1e-9).unwrap();
        let moved = part.offset_faces(&front, -30.0).unwrap();
        assert!(moved.is_manifold().unwrap());
        assert_eq!(lid(&moved), 2, "the panel is still apart from the lid");
        let outline = moved.edges().unwrap();
        let outline = outline.iter().filter(|e| e.label.as_ref().is_some_and(|l| l.role == "outline"));
        assert_eq!(outline.count(), 4);
        // What the push took: a slab of the box, less its rounds, and what it
        // cut off the bosses.
        assert!(moved.volume(1e-9).unwrap() < volume - 30.0 * 400.0 * 265.0 + 1.0);
    });
}

#[test]
fn merging_an_imprint_off_one_round_leaves_a_panel_on_another_as_it_was() {
    on_kernel(|k| {
        // A crate-sized block with its vertical edges rounded, and a 10 ×
        // 60 panel imprinted through it along X: on the round at each end of
        // its +Y side.
        let mut part = Brep::cuboid(k, [-200.0, -150.0, -125.0], [400.0, 300.0, 265.0])
            .unwrap()
            .fillet(
                &[smlib::EdgeBlend::round(
                    EdgeQuery::ParallelTo {
                        direction: [0.0, 0.0, 1.0],
                    },
                    20.0,
                )],
                &smlib::FilletOptions::default(),
            )
            .unwrap();
        for face in part.faces().unwrap() {
            part.set_face_label(face.index, &Label::new(BLOCK, "side").with_index(face.index as u32))
                .unwrap();
        }
        let tool = Brep::cuboid(k, [-500.0, 135.0, 20.0], [1000.0, 10.0, 60.0]).unwrap();
        let Imprint {
            brep: mut part,
            edges,
            ..
        } = part.imprint(tool).unwrap();
        for edge in edges.iter().filter(|e| e.new) {
            part.set_edge_label(edge.edge, &Label::new(TOOL, "outline")).unwrap();
        }
        let round_at = |face: &smlib::FaceGeometry, x: f64| {
            matches!(face.surface, smlib::Surface::Cylinder { origin, .. }
                if (origin[0] - x).abs() < 1e-6 && (origin[1] - 130.0).abs() < 1e-6)
        };
        // The panel on the +X round: the piece of it that doesn't reach
        // the lid.
        let panel = |b: &Brep| {
            b.faces()
                .unwrap()
                .into_iter()
                .find(|f| round_at(f, 180.0) && f.bounds[1][2] < 100.0)
                .unwrap()
        };
        let before = panel(&part);
        let area = part.face_centroid(before.index).unwrap().area;
        // The arc from 5 to 15 off the round's axis along Y, 60 high.
        let arc = 20.0 * ((15.0f64 / 20.0).asin() - (5.0f64 / 20.0).asin());
        assert!((area - arc * 60.0).abs() < 1e-6 * area, "{area}");

        // Merged back off the -X round, keeping every other edge.
        let faces = part.faces().unwrap();
        let on_far = |f: Option<u32>| f.is_some_and(|f| round_at(&faces[f as usize], -180.0));
        let all = part.edges().unwrap();
        let keep: Vec<usize> = all
            .iter()
            .filter(|e| e.label.is_none() || !e.faces.iter().all(|&f| on_far(f)))
            .map(|e| e.index)
            .collect();
        assert_eq!(all.len() - keep.len(), 4);
        part.merge_faces_keeping(&keep).unwrap();
        let far: Vec<_> = part.faces().unwrap().into_iter().filter(|f| round_at(f, -180.0)).collect();
        assert_eq!(far.len(), 1, "the -X round is one face again");
        // The kept panel is as it was: its box, and its area.
        let after = panel(&part);
        assert_eq!(after.bounds, before.bounds);
        assert_eq!(after.uv_min, before.uv_min);
        assert_eq!(after.uv_max, before.uv_max);
        assert!((part.face_centroid(after.index).unwrap().area - area).abs() < 1e-9 * area);
        // And a later join over it works.
        let volume = part.volume(1e-9).unwrap();
        let joined = part
            .union(Brep::cuboid(k, [190.0, 140.0, 40.0], [30.0, 30.0, 20.0]).unwrap())
            .unwrap();
        assert!(joined.is_manifold().unwrap());
        assert!(joined.volume(1e-9).unwrap() > volume);
    });
}
