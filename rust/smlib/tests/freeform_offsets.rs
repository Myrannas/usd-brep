//! Offsets and shells of free-form faces: a dome's top and a lofted tank's
//! skin offset out and in, each lying its distance off the face it came
//! from, and the tank opened at either cap. Before the free-form offset
//! patches every offset of the dome failed, and the tank's outward ones,
//! its seam 5 mm in and, after 11 s, 10 mm in (see `PATCHES.md`). Also how
//! far off smooth edges are (`Brep::edge_creases`), which says which solids
//! the kernel's shell can't take.
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{Brep, EdgeBlend, EdgeQuery, FilletOptions, FreeSurface, Kernel, Surface, samples};

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

const DISTANCES: [f64; 8] = [-10.0, -5.0, -2.0, -0.5, 0.5, 2.0, 5.0, 10.0];

fn dome(k: &Kernel) -> Brep {
    samples::domed_box(k, [-60.0, -40.0, -30.0], [120.0, 80.0, 42.0], 36.0).unwrap().0
}

/// The robustness survey's tank, its last ring `last` across.
fn tank(k: &Kernel, last: f64) -> Brep {
    let radii = [[27.0, 24.0], [41.4, 36.8], [45.0, 40.0], [40.5, 36.0], [last, 22.0]];
    samples::lofted_tank(k, [-80.0, 0.0, 0.0], 160.0, &radii).unwrap().0
}

/// The indices of the faces that aren't planar.
fn curved(brep: &Brep) -> Vec<usize> {
    brep.faces().unwrap().iter().filter(|f| f.normal.is_none()).map(|f| f.index).collect()
}

fn add(a: [f64; 3], b: [f64; 3], by: f64) -> [f64; 3] {
    [a[0] + b[0] * by, a[1] + b[1] * by, a[2] + b[2] * by]
}

/// How far, at most, the surface of `moved`'s face `face` lies from where
/// `from`'s surface moved `distance` along `outward` would be, over the
/// middle of `from`'s domain.
fn off_by(from: &FreeSurface, moved: &Brep, face: usize, distance: f64, outward: impl Fn([f64; 3], [f64; 3]) -> bool) -> f64 {
    let moved = FreeSurface::of_face(moved, face).unwrap();
    let mut worst: f64 = 0.0;
    for i in 0..=4 {
        for j in 0..=4 {
            let (point, normal) = from.at(0.2 + 0.15 * i as f64, 0.2 + 0.15 * j as f64).unwrap();
            let normal = if outward(point, normal) { normal } else { normal.map(|n| -n) };
            let target = add(point, normal, distance);
            worst = worst.max(moved.closest_point(target).unwrap().distance);
        }
    }
    worst
}

#[test]
fn a_domes_top_offsets_both_ways() {
    on_kernel(|k| {
        let part = dome(k);
        let volume = part.volume(1e-6).unwrap();
        let top = curved(&part);
        assert_eq!(top.len(), 1);
        let surface = FreeSurface::of_face(&part, top[0]).unwrap();
        for distance in DISTANCES {
            let moved = part.try_clone().unwrap().offset_faces(&top, distance).unwrap_or_else(|e| panic!("{distance}: {e}"));
            assert!(moved.is_manifold().unwrap(), "{distance}");
            let grown = moved.volume(1e-6).unwrap() - volume;
            assert!(grown * distance > 0.0, "{distance}: the volume changed by {grown}");
            let new_top = curved(&moved);
            assert_eq!(new_top.len(), 1, "{distance}");
            assert!(matches!(moved.faces().unwrap()[new_top[0]].surface, Surface::Nurbs { .. }), "{distance}");
            let off = off_by(&surface, &moved, new_top[0], distance, |_, n| n[2] > 0.0);
            assert!(off < 1e-4, "{distance}: the new top is {off} mm off the offset");
        }
    });
}

#[test]
fn a_tanks_skin_offsets_both_ways_across_its_seam() {
    on_kernel(|k| {
        let part = tank(k, 24.8);
        let volume = part.volume(1e-6).unwrap();
        let skin = curved(&part);
        assert_eq!(skin.len(), 1);
        let surface = FreeSurface::of_face(&part, skin[0]).unwrap();
        // Out from the tank's axis, along x.
        let outward = |p: [f64; 3], n: [f64; 3]| p[1] * n[1] + p[2] * n[2] > 0.0;
        for distance in DISTANCES {
            let moved = part.try_clone().unwrap().offset_faces(&skin, distance).unwrap_or_else(|e| panic!("{distance}: {e}"));
            assert!(moved.is_manifold().unwrap(), "{distance}");
            assert_eq!(moved.face_count(), 3, "{distance}");
            let grown = moved.volume(1e-6).unwrap() - volume;
            assert!(grown * distance > 0.0, "{distance}: the volume changed by {grown}");
            let new_skin = curved(&moved);
            assert_eq!(new_skin.len(), 1, "{distance}");
            let off = off_by(&surface, &moved, new_skin[0], distance, outward);
            assert!(off < 1e-4, "{distance}: the new skin is {off} mm off the offset");
        }
    });
}

#[test]
fn a_tank_opens_at_either_cap() {
    on_kernel(|k| {
        // The backend's tank is 24.75 across at its last ring, where the far
        // cap's ring failed to trim at 3 mm.
        for last in [24.8, 24.75] {
            let part = tank(k, last);
            let volume = part.volume(1e-6).unwrap();
            let caps: Vec<usize> = part.faces().unwrap().iter().filter(|f| f.normal.is_some()).map(|f| f.index).collect();
            assert_eq!(caps.len(), 2);
            for cap in caps {
                for thickness in [1.0, 2.0, 3.0, 5.0, 8.0] {
                    let shell = part.shell(thickness, &[cap]).unwrap_or_else(|e| panic!("{last} cap {cap} {thickness}: {e}"));
                    assert!(shell.brep.is_manifold().unwrap(), "{last} cap {cap} {thickness}");
                    let left = shell.brep.volume(1e-6).unwrap();
                    assert!(left > 0.0 && left < volume, "{last} cap {cap} {thickness}: {left} of {volume}");
                }
            }
        }
    });
}

/// The worst of the edges' creases for which `which` holds, by index.
fn worst(brep: &Brep, which: impl Fn(usize, f64) -> bool) -> f64 {
    brep.edge_creases().unwrap().into_iter().enumerate().filter_map(|(i, c)| c.filter(|&c| which(i, c))).fold(0.0, f64::max)
}

#[test]
fn edge_creases_find_where_faces_part_from_smooth() {
    on_kernel(|k| {
        // A box: every edge square.
        let block = Brep::cuboid(k, [0.0; 3], [100.0, 60.0, 40.0]).unwrap();
        let creases = block.edge_creases().unwrap();
        assert_eq!(creases.len(), 12);
        assert!(creases.iter().all(|c| c.is_some_and(|c| (c - 90.0).abs() < 1e-6)), "{creases:?}");

        // A Block with every edge rounded: smooth everywhere, its rounds
        // and corner spheres exact.
        let rounded = block.fillet(&[EdgeBlend::round(EdgeQuery::Sharp { min_angle_deg: 30.0 }, 8.0)], &FilletOptions::default()).unwrap();
        let creases = rounded.edge_creases().unwrap();
        assert!(creases.iter().flatten().all(|&c| c < 5e-4), "{creases:?}");

        // N-sided fills meet the rounds within 0.11° (`FreeSurface::fill`),
        // but not exactly.
        let filled = samples::rounded_box(k, [-50.0, -30.0, -20.0], [100.0, 60.0, 40.0], 8.0).unwrap().brep;
        let off = worst(&filled, |_, c| c < 1.0);
        assert!(off > 0.05 && off < 0.11, "{off}");

        // A tank's cap rims rounded r 2 by SMLib's fillet on its NURBS skin:
        // off smooth near the rims' ends, though not at their middles.
        let part = tank(k, 24.8);
        let skin = curved(&part)[0];
        let rims: Vec<usize> = part
            .edges()
            .unwrap()
            .iter()
            .filter(|e| e.faces.contains(&Some(skin as u32)) && e.faces.iter().flatten().any(|&f| f as usize != skin))
            .map(|e| e.index)
            .collect();
        assert_eq!(rims.len(), 2);
        let rounded = part.fillet(&[EdgeBlend::round(EdgeQuery::Indices(rims), 2.0)], &FilletOptions::default()).unwrap();
        let edges = rounded.edges().unwrap();
        let off = worst(&rounded, |i, c| c < 1.0 && edges[i].bend_deg().is_some_and(|middle| middle < 5e-4));
        assert!(off > 5e-4, "{off}");
    });
}
