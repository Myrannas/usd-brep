//! Faces taken out of a solid and rebuilt from the rim, sewn back in one
//! call (`Brep::rebuild_faces`): what the Free-form parts spec's Rebuild
//! needs, from its continuity spike. A rebuilt side asked to be Tangent
//! meets the face beyond it with no crease, the stitch says how far apart
//! what it left open is, and the measures that check it read true.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, EdgeBlend, EdgeQuery, ExtendShape, FilletOptions, FreeSurface, Kernel, Label, Meets, RebuildSurface, Side,
    Surface, SurfaceSide, Tool,
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

/// Hollow's limit: faces meeting further off smooth than this, and under
/// a degree, stop SMLib's shell (game-cad's `hollow.rs`).
const CREASE_DEG: f64 = 0.0005;
/// The spec's stitch tolerance.
const STITCH: f64 = 0.01;
const REBUILT: u32 = 9000;

fn label() -> Option<Label> {
    Some(Label::new(REBUILT, "patch"))
}

/// A 100 × 80 × 60 block from the origin with its top front edge rounded
/// r 10, and the round's index.
fn rounded_edge(k: &Kernel) -> (Brep, usize) {
    let block = Brep::cuboid(k, [0.0; 3], [100.0, 80.0, 60.0]).unwrap();
    let edge = EdgeQuery::BetweenFaces { a: [0.0, 0.0, 1.0], b: [0.0, -1.0, 0.0] };
    let solid = block.fillet(&[EdgeBlend::round(edge, 10.0)], &FilletOptions::default()).unwrap();
    let round = solid.faces().unwrap().into_iter().find(|f| f.normal.is_none()).unwrap().index;
    (solid, round)
}

/// The round's rim: its edges to the top and the front, Tangent, and to
/// the ends, Position, in order round the loop.
fn rim(solid: &Brep, round: usize) -> Vec<Side> {
    let faces = solid.faces().unwrap();
    let mut rim: Vec<(usize, [f64; 3])> = solid
        .edges()
        .unwrap()
        .into_iter()
        .filter(|e| e.faces.contains(&Some(round as u32)))
        .map(|e| {
            let kept = e.faces.iter().flatten().map(|&f| f as usize).find(|&f| f != round).unwrap();
            (e.index, faces[kept].normal.unwrap())
        })
        .collect();
    // Top, an end, front, the other end: round the loop.
    let order = |n: [f64; 3]| match n {
        n if n[2] > 0.99 => 0,
        n if n[0] < -0.99 => 1,
        n if n[1] < -0.99 => 2,
        _ => 3,
    };
    rim.sort_by_key(|&(_, n)| order(n));
    rim.into_iter()
        .map(|(edge, n)| Side { edge, meets: if n[0].abs() > 0.99 { Meets::Position } else { Meets::Tangent } })
        .collect()
}

/// The worst crease on each rebuilt face's edges to a kept face.
fn rim_creases(brep: &Brep) -> Vec<f64> {
    let faces = brep.faces().unwrap();
    let rebuilt = |f: Option<u32>| f.is_some_and(|f| faces[f as usize].label.as_ref().is_some_and(|l| l.feature == REBUILT));
    let creases = brep.edge_creases().unwrap();
    brep.edges()
        .unwrap()
        .into_iter()
        .filter(|e| e.faces.iter().filter(|&&f| rebuilt(f)).count() == 1 && e.faces.iter().flatten().count() == 2)
        .filter_map(|e| creases[e.index])
        .collect()
}

#[test]
fn a_coons_patch_closes_a_rounded_edge_tangent_to_both_faces() {
    on_kernel(|k| {
        let (solid, round) = rounded_edge(k);
        let sides = rim(&solid, round);
        let surface = RebuildSurface { tool: Tool::Coons, sides: sides.clone(), tangent_tol_deg: 1e-4, fit_tolerance: 0.0, label: label() };
        let rebuilt = solid.rebuild_faces(&[round], &[surface], STITCH).unwrap();
        assert!(rebuilt.report.solid);
        assert_eq!(rebuilt.report.open_edges, 0);
        assert_eq!(rebuilt.report.faces_made, 1);
        assert_eq!(rebuilt.brep.face_count(), solid.face_count());
        // The two Tangent sides meet with no crease, the ends square.
        let creases = rim_creases(&rebuilt.brep);
        assert_eq!(creases.len(), 4);
        assert_eq!(creases.iter().filter(|&&c| c < 1e-6).count(), 2, "{creases:?}");
        assert_eq!(creases.iter().filter(|&&c| (c - 90.0).abs() < 1e-6).count(), 2, "{creases:?}");
        // It's the round it replaced, near enough.
        let before = solid.volume(1e-6).unwrap();
        let after = rebuilt.brep.volume(1e-6).unwrap();
        assert!((after - before).abs() < 1e-6 * before, "{before} → {after}");
        // And it shells across it.
        let bottom = rebuilt.brep.faces().unwrap().into_iter().find(|f| f.normal.is_some_and(|n| n[2] < -0.99)).unwrap();
        assert!(rebuilt.brep.shell(2.0, &[bottom.index]).is_ok());
    });
}

#[test]
fn a_fill_closes_it_in_a_patch_a_side() {
    on_kernel(|k| {
        let (solid, round) = rounded_edge(k);
        let surface = RebuildSurface { tool: Tool::Fill, sides: rim(&solid, round), tangent_tol_deg: 0.01, fit_tolerance: 0.0, label: label() };
        let rebuilt = solid.rebuild_faces(&[round], &[surface], STITCH).unwrap();
        assert!(rebuilt.report.solid);
        // Four patches meeting in the middle, the rim's edges split where
        // they meet it.
        assert_eq!(rebuilt.report.faces_made, 4);
        let creases = rim_creases(&rebuilt.brep);
        assert!(creases.iter().filter(|&&c| c < 45.0).all(|&c| c < CREASE_DEG), "{creases:?}");
    });
}

#[test]
fn an_edge_blend_between_two_sides_leaves_its_ends_open_and_says_how_far() {
    on_kernel(|k| {
        let (solid, round) = rounded_edge(k);
        let tangent: Vec<Side> = rim(&solid, round).into_iter().filter(|s| s.meets == Meets::Tangent).collect();
        let surface = RebuildSurface { tool: Tool::EdgeBlend { order: 1 }, sides: tangent, tangent_tol_deg: 0.0, fit_tolerance: 0.0, label: label() };
        let Err(failure) = solid.rebuild_faces(&[round], std::slice::from_ref(&surface), STITCH) else {
            panic!("an edge blend's ends closed");
        };
        // Its ends' curves aren't the round's ends' arcs, so they stay 0.2
        // mm off the ends' faces.
        assert_eq!(failure.report.open_edges, 4);
        let gap = failure.report.open_edge_gap.unwrap();
        assert!((0.1..0.3).contains(&gap), "{gap}");
        // StitchIntoSolid glues it shut, widening its tolerance, and says
        // only how wide what it glued was.
        let glued = solid.rebuild_faces(&[round], &[surface], 0.0).unwrap();
        assert!(glued.report.solid);
        assert!(glued.report.edge_gap.max(glued.report.vertex_gap) > 0.1, "{:?}", glued.report);
    });
}

#[test]
fn a_corner_blend_meets_a_rounded_box_corner_within_hollows_limit() {
    on_kernel(|k| {
        // A box's vertical edges rounded r 20, then its top's edges r 10:
        // at each corner a blend between the top, the vertical round and
        // two top rounds.
        let block = Brep::cuboid(k, [0.0; 3], [100.0, 80.0, 60.0]).unwrap();
        let vertical = EdgeQuery::ParallelTo { direction: [0.0, 0.0, 1.0] };
        let block = block.fillet(&[EdgeBlend::round(vertical, 20.0)], &FilletOptions::default()).unwrap();
        let top = block.faces().unwrap().into_iter().find(|f| f.normal.is_some_and(|n| n[2] > 0.99)).unwrap().index;
        let edges: Vec<usize> =
            block.edges().unwrap().into_iter().filter(|e| e.faces.contains(&Some(top as u32))).map(|e| e.index).collect();
        let solid = block.fillet(&[EdgeBlend::round(EdgeQuery::Indices(edges), 10.0)], &FilletOptions::default()).unwrap();
        let corner = solid
            .faces()
            .unwrap()
            .into_iter()
            // SMLib's corner blend: neither the rounds' cylinders nor a
            // plane, read as NURBS or a surface of its own.
            .filter(|f| !matches!(f.surface, Surface::Plane { .. } | Surface::Cylinder { .. }))
            .find(|f| f.bounds[0][0] > 79.0 && f.bounds[0][1] > 59.0 && f.bounds[0][2] > 49.0)
            .unwrap()
            .index;
        let sides: Vec<Side> = solid
            .edges()
            .unwrap()
            .into_iter()
            .filter(|e| e.faces.contains(&Some(corner as u32)))
            .map(|e| Side { edge: e.index, meets: Meets::Tangent })
            .collect();
        assert_eq!(sides.len(), 4);
        // Fitted to the brep's tolerance, the rounds' arcs leave it a hair
        // past Hollow's limit; fitted to a millionth of a millimetre, well
        // inside it.
        for (fit, within) in [(0.0, 1e-3), (1e-6, CREASE_DEG)] {
            let surface = RebuildSurface { tool: Tool::Corner, sides: sides.clone(), tangent_tol_deg: 1e-4, fit_tolerance: fit, label: label() };
            let rebuilt = solid.rebuild_faces(&[corner], &[surface], STITCH).unwrap();
            assert!(rebuilt.report.solid);
            let creases = rim_creases(&rebuilt.brep);
            assert_eq!(creases.len(), 4);
            assert!(creases.iter().all(|&c| c < within), "fitted to {fit}: {creases:?}");
        }
    });
}

#[test]
fn edge_curvatures_read_a_round_against_the_plane_beside_it() {
    on_kernel(|k| {
        let (solid, round) = rounded_edge(k);
        let curvatures = solid.edge_curvatures().unwrap();
        let faces = solid.faces().unwrap();
        for edge in solid.edges().unwrap().into_iter().filter(|e| e.faces.contains(&Some(round as u32))) {
            let kept = edge.faces.iter().flatten().map(|&f| f as usize).find(|&f| f != round).unwrap();
            let curvature = curvatures[edge.index].unwrap();
            let (on_round, on_kept) = if edge.faces[0] == Some(round as u32) {
                (curvature.a, curvature.b)
            } else {
                (curvature.b, curvature.a)
            };
            if faces[kept].normal.unwrap()[0].abs() > 0.99 {
                // Across an end the round runs straight.
                assert!(on_round.iter().all(|c| c.abs() < 1e-9), "{on_round:?}");
            } else {
                // Across the top and front it's convex, 1/10, and the
                // planes flat.
                assert!(on_round.iter().all(|c| (c - 0.1).abs() < 1e-9), "{on_round:?}");
                assert!(on_kept.iter().all(|c| c.abs() < 1e-9), "{on_kept:?}");
                assert!((curvature.worst_difference() - 0.1).abs() < 1e-9);
            }
        }
    });
}

#[test]
fn a_surface_carried_on_meets_itself_where_it_was_cut_off() {
    on_kernel(|k| {
        let block = Brep::cuboid(k, [0.0; 3], [10.0, 10.0, 10.0]).unwrap();
        let plane = FreeSurface::of_face(&block, 0).unwrap();
        let [min, max] = plane.domain().unwrap();
        for shape in [ExtendShape::Linear, ExtendShape::Smooth, ExtendShape::Natural, ExtendShape::Ruled] {
            let extended = plane.extended(SurfaceSide::UMax, 20.0, shape).unwrap();
            let [emin, emax] = extended.domain().unwrap();
            assert!(emax[0] > max[0] && (emin[0] - min[0]).abs() < 1e-9, "{shape:?}: {emin:?} {emax:?}");
            for k in 0..=4 {
                let v = min[1] + (max[1] - min[1]) * k as f64 / 4.0;
                let (p, n) = plane.evaluate(max[0], v).unwrap();
                let near = extended.closest_point(p).unwrap();
                assert!(near.distance < 1e-9, "{shape:?}: {}", near.distance);
                let (_, m) = extended.evaluate(near.uv[0], near.uv[1]).unwrap();
                assert!((n[0] * m[0] + n[1] * m[1] + n[2] * m[2] - 1.0).abs() < 1e-12, "{shape:?}");
            }
        }
    });
}
