//! Free-form surfaces made into solids: a loft closed round and capped, a box
//! whose top is a NURBS dome, and a rounded box whose corners are N-sided
//! fills. Each must close into a manifold solid of the right volume and
//! tessellate. Also each surface's control net and closest point, split
//! free-form faces merging back, and, run
//! by hand, how booleans, fillets, offsets and shells fare on them
//! (`freeform_robustness`).
//!
//! As in `cases.rs`, every test sends its kernel work to one shared kernel
//! thread.

use std::f64::consts::PI;
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::{
    Brep, Curve, EdgeBlend, EdgeQuery, FilletOptions, FreeSurface, Kernel, Sheets, Surface,
    TessellationParams, Tolerance, samples,
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

fn kinds(brep: &Brep) -> String {
    let mut counts = std::collections::BTreeMap::new();
    for face in brep.faces().unwrap() {
        let kind = match face.surface {
            Surface::Plane { .. } => "plane",
            Surface::Cylinder { .. } => "cylinder",
            Surface::Cone { .. } => "cone",
            Surface::Sphere { .. } => "sphere",
            Surface::Torus { .. } => "torus",
            Surface::Nurbs { .. } if face.normal.is_some() => "flat nurbs",
            Surface::Nurbs { .. } => "nurbs",
            _ => "other",
        };
        *counts.entry(kind).or_insert(0) += 1;
    }
    format!("{counts:?}")
}

fn report(name: &str, brep: &Brep) -> f64 {
    let volume = brep.volume(1e-6).unwrap();
    let mesh = brep.tessellate(&TessellationParams::default()).unwrap();
    println!(
        "{name}: {} faces {} edges, {}, volume {volume:.3}, {} triangles",
        brep.face_count(),
        brep.edge_count(),
        kinds(brep),
        mesh.triangle_count()
    );
    volume
}

#[test]
fn curves_from_poles_and_through_points() {
    on_kernel(|k| {
        // A quadratic rational quarter circle, exact.
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let arc = Curve::nurbs(
            k,
            2,
            &[[10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
            Some(&[1.0, w, 1.0]),
            &[0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        )
        .unwrap();
        let shape = arc.shape().unwrap();
        assert_eq!((shape.degree, shape.poles, shape.rational), (2, 3, true));
        for s in [0.0, 0.3, 0.5, 0.9, 1.0] {
            let p = arc.at(s).unwrap();
            assert!(
                ((p[0] * p[0] + p[1] * p[1]).sqrt() - 10.0).abs() < 1e-9,
                "{p:?}"
            );
        }
        let through = Curve::through(
            k,
            &[[0.0; 3], [1.0, 2.0, 0.0], [3.0, 1.0, 1.0], [4.0, 4.0, 2.0]],
            3,
            false,
        )
        .unwrap();
        assert!(through.at(0.0).unwrap().iter().all(|v| v.abs() < 1e-9));
        let end = through.at(1.0).unwrap();
        assert!((end[0] - 4.0).abs() < 1e-9 && (end[2] - 2.0).abs() < 1e-9);
        let mut turned = through.copy().unwrap();
        turned.reverse().unwrap();
        assert!((turned.at(0.0).unwrap()[0] - 4.0).abs() < 1e-9);
    });
}

#[test]
fn surfaces_from_curves() {
    on_kernel(|k| {
        // Coons: four straight sides of a twisted quad give a bilinear patch.
        let c = [
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, 10.0, 5.0],
            [0.0, 10.0, 0.0],
        ];
        let sides: Vec<Curve> = (0..4)
            .map(|i| Curve::line(k, c[i], c[(i + 1) % 4]).unwrap())
            .collect();
        let coons = FreeSurface::coons(k, [&sides[0], &sides[1], &sides[2], &sides[3]], None, 1e-4)
            .unwrap();
        println!("coons, bilinear: {:?}", coons.shape().unwrap());
        // Its corners are the quad's, and its middle the corners' mean.
        let (middle, _) = coons.at(0.5, 0.5).unwrap();
        assert!((middle[2] - 1.25).abs() < 1e-6, "{middle:?}");
        // Gordon: a 3 by 3 network on a saddle z = (x² - y²) / 20.
        let saddle = |x: f64, y: f64| [x, y, (x * x - y * y) / 20.0];
        let along = |fixed_y: f64| {
            let p: Vec<_> = (0..=8)
                .map(|i| saddle(-10.0 + 2.5 * i as f64, fixed_y))
                .collect();
            Curve::through(k, &p, 3, false).unwrap()
        };
        let across = |fixed_x: f64| {
            let p: Vec<_> = (0..=8)
                .map(|j| saddle(fixed_x, -10.0 + 2.5 * j as f64))
                .collect();
            Curve::through(k, &p, 3, false).unwrap()
        };
        let u: Vec<Curve> = [-10.0, 0.0, 10.0].iter().map(|&y| along(y)).collect();
        let v: Vec<Curve> = [-10.0, 0.0, 10.0].iter().map(|&x| across(x)).collect();
        let gordon = FreeSurface::gordon(
            k,
            &u.iter().collect::<Vec<_>>(),
            &v.iter().collect::<Vec<_>>(),
        );
        let gordon = gordon.unwrap();
        let mut worst: f64 = 0.0;
        for i in 0..=10 {
            for j in 0..=10 {
                let (p, _) = gordon.at(i as f64 / 10.0, j as f64 / 10.0).unwrap();
                worst = worst.max((p[2] - (p[0] * p[0] - p[1] * p[1]) / 20.0).abs());
            }
        }
        println!(
            "gordon: {:?}, off the saddle by up to {worst:.2e}",
            gordon.shape().unwrap()
        );
        // Chord-length curves aren't parameter-compatible, so between them
        // it strays a little.
        assert!(worst < 0.3);
        // A loft through three sections.
        let section = |x: f64, r: f64| {
            let p: Vec<_> = (0..12)
                .map(|i| {
                    let a = 2.0 * PI * i as f64 / 12.0;
                    [x, r * a.cos(), r * a.sin()]
                })
                .collect();
            Curve::through(k, &p, 3, true).unwrap()
        };
        let s = [section(0.0, 5.0), section(10.0, 8.0), section(20.0, 4.0)];
        let loft = FreeSurface::loft(k, &[&s[0], &s[1], &s[2]], 3).unwrap();
        let shape = loft.shape().unwrap();
        println!("loft: {shape:?}");
        assert!(!shape.rational && shape.degree[1] <= 3);
    });
}

#[test]
fn lofted_tank_closes() {
    on_kernel(|k| {
        let radii = [
            [30.0, 24.0],
            [45.0, 36.0],
            [50.0, 40.0],
            [44.0, 35.0],
            [28.0, 22.0],
        ];
        let (tank, solid) = samples::lofted_tank(k, [0.0, 0.0, 50.0], 160.0, &radii).unwrap();
        let volume = report("tank", &tank);
        // Elliptical sections, area pi ry rz, integrated roughly by Simpson.
        let areas: Vec<f64> = radii.iter().map(|[a, b]| PI * a * b).collect();
        let h = 40.0;
        let simpson =
            h / 3.0 * (areas[0] + 4.0 * areas[1] + 2.0 * areas[2] + 4.0 * areas[3] + areas[4]);
        println!("tank: solid {solid}, Simpson's estimate {simpson:.0}");
        assert!(solid && tank.is_manifold().unwrap());
        assert!((volume - simpson).abs() / simpson < 0.02);
    });
}

#[test]
fn domed_box_closes() {
    on_kernel(|k| {
        let (part, solid) =
            samples::domed_box(k, [-60.0, -40.0, 0.0], [120.0, 80.0, 50.0], 30.0).unwrap();
        let volume = report("domed box", &part);
        println!(
            "domed box: solid {solid}, box alone {}",
            120.0 * 80.0 * 50.0
        );
        assert!(solid && part.is_manifold().unwrap());
        assert!(volume > 120.0 * 80.0 * 50.0);
        // The walls and floor are flat, whatever surface they're on; the dome
        // isn't.
        let faces = part.faces().unwrap();
        assert_eq!(faces.iter().filter(|f| f.normal.is_some()).count(), 5);
        let dome = faces.iter().find(|f| f.normal.is_none()).unwrap();
        assert!(matches!(dome.surface, Surface::Nurbs { .. }));
    });
}

#[test]
fn rounded_box_corners_fill() {
    on_kernel(|k| {
        let rounded =
            samples::rounded_box(k, [-50.0, -30.0, 0.0], [100.0, 60.0, 40.0], 10.0).unwrap();
        let volume = report("rounded box", &rounded.brep);
        println!(
            "rounded box: solid {}, exact volume {:.3}, off by {:.2e}, corners off the sphere by up to {:.3e}",
            rounded.solid,
            rounded.exact_volume,
            (volume - rounded.exact_volume) / rounded.exact_volume,
            rounded.corner_deviation
        );
        assert!(rounded.solid && rounded.brep.is_manifold().unwrap());
        assert!((volume - rounded.exact_volume).abs() / rounded.exact_volume < 1e-3);
        // `FreeSurface::fill`'s docs: up to 4% of the radius off the sphere.
        assert!(
            rounded.corner_deviation < 0.045 * 10.0,
            "{}",
            rounded.corner_deviation
        );
        assert_eq!(
            rounded
                .brep
                .faces()
                .unwrap()
                .iter()
                .filter(|f| f.normal.is_some())
                .count(),
            6
        );
    });
}

/// A slab imprinted through a rounded box's corner, splitting its fill
/// patches and its rounds, merges back to the 42 faces it had, whether the
/// new edges alone are merged away (`merge_faces_keeping`) or every edge
/// that can be (`merge_faces`). Each piece of a split NURBS face lies on its
/// own copy of the surface, and some pieces of a round come back as surfaces
/// of extrusion, which the kernel's merge left apart before the fork's
/// `bFreeForm`.
#[test]
fn split_freeform_faces_merge_back() {
    on_kernel(|k| {
        let rounded = samples::rounded_box(k, [-50.0, -30.0, -20.0], [100.0, 60.0, 40.0], 8.0)
            .unwrap()
            .brep;
        let whole = kinds(&rounded);
        let slab = || Brep::cuboid(k, [44.0, -100.0, -100.0], [1.0, 200.0, 200.0]).unwrap();
        for keeping in [true, false] {
            let imprint = rounded.try_clone().unwrap().imprint_classified(slab()).unwrap();
            let mut brep = imprint.brep;
            let split = brep.face_count();
            assert!(split > 42, "the slab splits faces: {split}");
            if keeping {
                let new: Vec<usize> =
                    imprint.edges.iter().filter(|e| e.new).map(|e| e.edge).collect();
                let keep: Vec<usize> =
                    (0..brep.edge_count()).filter(|e| !new.contains(e)).collect();
                brep.merge_faces_keeping(&keep).unwrap();
            } else {
                brep.merge_faces().unwrap();
            }
            println!(
                "keeping {keeping}: 42 faces, {split} imprinted, {} merged back",
                brep.face_count()
            );
            assert_eq!(brep.face_count(), 42, "keeping {keeping}");
            assert_eq!(kinds(&brep), whole, "keeping {keeping}");
            assert!(brep.is_manifold().unwrap());
        }
    });
}

/// A box whose top is a flat NURBS, a dome with no lift, the rest planar
/// faces bounded by lines: since 1.0.1 the top reports its
/// normal as a plane's does (`SmSurface::IsPlanar`), though it stays a
/// NURBS, and the planar faces bounded by curves are planes.
#[test]
fn flat_nurbs_faces_are_flat() {
    on_kernel(|k| {
        let (x, y, z) = (30.0, 20.0, 10.0);
        let low = [[0.0, 0.0, 0.0], [x, 0.0, 0.0], [x, y, 0.0], [0.0, y, 0.0]];
        let high = low.map(|[a, b, _]| [a, b, z]);
        let lines = |c: &[[f64; 3]]| -> Vec<Curve> {
            (0..c.len())
                .map(|i| Curve::line(k, c[i], c[(i + 1) % c.len()]).unwrap())
                .collect()
        };
        let mut sheets = Sheets::new(k).unwrap();
        let flat = samples::dome(k, [0.0, 0.0], [x, y], z, 0.0).unwrap();
        assert_eq!(sheets.add(&flat, None).unwrap(), 1);
        sheets
            .add_planar(&lines(&low).iter().collect::<Vec<_>>(), None)
            .unwrap();
        for i in 0..4 {
            let j = (i + 1) % 4;
            let wall = lines(&[low[i], low[j], high[j], high[i]]);
            sheets
                .add_planar(&wall.iter().collect::<Vec<_>>(), None)
                .unwrap();
        }
        let (part, solid) = sheets.stitch(0.0).unwrap();
        let volume = report("flat-topped box", &part);
        assert!(solid && part.is_manifold().unwrap());
        assert!((volume - x * y * z).abs() < 1e-6 * x * y * z, "{volume}");
        let faces = part.faces().unwrap();
        let top = faces
            .iter()
            .find(|f| matches!(f.surface, Surface::Nurbs { .. }))
            .expect("the NURBS top");
        let normal = top.normal.expect("a flat NURBS has a normal");
        assert!((normal[2].abs() - 1.0).abs() < 1e-9, "{normal:?}");
        let planes = faces
            .iter()
            .filter(|f| matches!(f.surface, Surface::Plane { .. }))
            .count();
        assert_eq!(planes, 5);
    });
}

#[test]
fn sheets_report_stitching() {
    on_kernel(|k| {
        let dome = samples::dome(k, [0.0, 0.0], [10.0, 10.0], 5.0, 3.0).unwrap();
        let mut sheets = Sheets::new(k).unwrap();
        sheets.add(&dome, None).unwrap();
        let (_, stitched) = sheets.stitch_into_solid().unwrap();
        println!("a lone dome: {stitched:?}");
        assert!(!stitched.solid);
    });
}

/// One corner of three rounds of radius 10 meeting, filled with and without
/// cross-boundary derivatives of several lengths: how far the fill strays
/// from the sphere the exact corner is on, and how far its normals turn from
/// the rounds' along their shared edges (0 for tangent).
#[test]
fn n_sided_fill_quality() {
    on_kernel(|k| {
        let r = 10.0;
        let arc = |a: [f64; 3], b: [f64; 3]| {
            let p: Vec<[f64; 3]> = (0..9)
                .map(|i| {
                    let t = std::f64::consts::FRAC_PI_2 * i as f64 / 8.0;
                    [0, 1, 2].map(|j| r * (a[j] * t.cos() + b[j] * t.sin()))
                })
                .collect();
            Curve::through(k, &p, 3, false).unwrap()
        };
        let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        // Round along z ends in x→y, along x in y→z, along y in z→x.
        let bounds = [arc(x, y), arc(y, z), arc(z, x)];
        let into = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        for scale in [
            None,
            Some(0.5),
            Some(1.0),
            Some(std::f64::consts::FRAC_PI_2),
            Some(2.0),
        ] {
            let derivs: Option<Vec<Curve>> = scale.map(|s| {
                (0..3)
                    .map(|i| bounds[i].constant_like(into[i].map(|v| v * r * s)).unwrap())
                    .collect()
            });
            let ders: Option<Vec<Option<&Curve>>> =
                derivs.as_ref().map(|d| d.iter().map(Some).collect());
            let patches = match FreeSurface::fill(
                k,
                &bounds.iter().collect::<Vec<_>>(),
                ders.as_deref(),
                1.0,
            ) {
                Ok(p) => p,
                Err(e) => {
                    println!("derivatives × {scale:?}: {e}");
                    continue;
                }
            };
            let (mut off, mut turn): (f64, f64) = (0.0, 0.0);
            for patch in &patches {
                for i in 0..=8 {
                    for j in 0..=8 {
                        let (p, n) = patch.at(i as f64 / 8.0, j as f64 / 8.0).unwrap();
                        let d = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                        off = off.max((d - r).abs());
                        // On the boundary, the rounds' normal is the sphere's.
                        let on_edge = [p[0], p[1], p[2]].iter().any(|v| v.abs() < 1e-6);
                        if on_edge {
                            let cos = (n[0] * p[0] + n[1] * p[1] + n[2] * p[2]).abs() / d;
                            turn = turn.max(cos.clamp(-1.0, 1.0).acos().to_degrees());
                        }
                    }
                }
            }
            println!(
                "derivatives × {scale:?}: {} patches, off the sphere by up to {off:.3} mm, normals turned up to {turn:.2}° along the rounds",
                patches.len()
            );
            // As `FreeSurface::fill`'s docs say: tangent to the rounds, but
            // up to 4% of the radius off the sphere, whatever the
            // derivatives.
            assert_eq!(patches.len(), 3);
            assert!(off < 0.045 * r, "{off}");
            assert!(turn < 0.2, "{turn}");
        }
    });
}

/// The points of two surfaces on a grid over their domains, which must be
/// the same.
fn same_surface(a: &FreeSurface, b: &FreeSurface) {
    assert_eq!(a.domain().unwrap(), b.domain().unwrap());
    let [min, max] = a.domain().unwrap();
    for i in 0..=8 {
        for j in 0..=8 {
            let u = min[0] + (max[0] - min[0]) * i as f64 / 8.0;
            let v = min[1] + (max[1] - min[1]) * j as f64 / 8.0;
            let (p, _) = a.evaluate(u, v).unwrap();
            let (q, _) = b.evaluate(u, v).unwrap();
            assert!(distance(p, q) < 1e-9, "({u}, {v}): {p:?} against {q:?}");
        }
    }
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
}

/// A surface's control net comes back as it was made, and a surface made
/// from the net another gives is that surface: a NURBS's own, a rational
/// loft's with its weights, and an analytic cylinder's as a NURBS.
#[test]
fn nets_round_trip() {
    on_kernel(|k| {
        let dome = samples::dome(k, [0.0, 0.0], [40.0, 20.0], 5.0, 8.0).unwrap();
        let net = dome.net().unwrap();
        assert_eq!(net.degree, [3, 3]);
        assert_eq!((net.poles.len(), net.poles[0].len()), (5, 5));
        assert!(net.weights.is_none());
        let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        assert_eq!(net.knots, [knots.clone(), knots]);
        assert_eq!(net.poles[4][0], [40.0, 0.0, 5.0]);
        assert_eq!(net.poles[2][2], [20.0, 10.0, 5.0 + 8.0 * 1.2 * 1.2]);
        same_surface(&dome, &FreeSurface::nurbs(k, &net).unwrap());

        // The tank's skin, lofted through rational ellipses.
        let radii = [[30.0, 24.0], [45.0, 36.0], [50.0, 40.0]];
        let (tank, _) = samples::lofted_tank(k, [0.0, 0.0, 0.0], 100.0, &radii).unwrap();
        let skin = tank
            .faces()
            .unwrap()
            .into_iter()
            .find(|f| f.normal.is_none())
            .unwrap();
        let skin = FreeSurface::of_face(&tank, skin.index).unwrap();
        let net = skin.net().unwrap();
        let weights = net.weights.as_ref().expect("a rational skin");
        assert_eq!(weights.len(), net.poles.len());
        assert!(weights.iter().flatten().any(|&w| (w - 1.0).abs() > 1e-3));
        for (along, knots) in net.knots.iter().enumerate() {
            let count = if along == 0 {
                net.poles.len()
            } else {
                net.poles[0].len()
            };
            assert_eq!(knots.len(), count + net.degree[along] as usize + 1);
        }
        same_surface(&skin, &FreeSurface::nurbs(k, &net).unwrap());

        // A cylinder's side, converted: every point of the surface its net
        // makes lies on the cylinder.
        let cylinder = Brep::cylinder(k, [0.0, 0.0, 0.0], 7.0, 20.0).unwrap();
        let side = cylinder
            .faces()
            .unwrap()
            .into_iter()
            .find(|f| matches!(f.surface, Surface::Cylinder { .. }))
            .unwrap();
        let side = FreeSurface::of_face(&cylinder, side.index).unwrap();
        let made = FreeSurface::nurbs(k, &side.net().unwrap()).unwrap();
        for i in 0..=8 {
            for j in 0..=8 {
                let (p, _) = made.at(i as f64 / 8.0, j as f64 / 8.0).unwrap();
                assert!(
                    ((p[0] * p[0] + p[1] * p[1]).sqrt() - 7.0).abs() < 1e-9,
                    "{p:?}"
                );
                assert!((-1e-9..=20.0 + 1e-9).contains(&p[2]), "{p:?}");
            }
        }
    });
}

/// The nearest point of a surface to points off it along its normals is
/// the foot of each, at its (u, v): on the dome, the rational tank skin and
/// a kernel cylinder.
#[test]
fn closest_points_are_the_feet_of_normals() {
    on_kernel(|k| {
        let radii = [[30.0, 24.0], [45.0, 36.0], [50.0, 40.0]];
        let (tank, _) = samples::lofted_tank(k, [0.0, 0.0, 0.0], 100.0, &radii).unwrap();
        let skin = tank
            .faces()
            .unwrap()
            .into_iter()
            .find(|f| f.normal.is_none())
            .unwrap();
        let cylinder = Brep::cylinder(k, [0.0, 0.0, 0.0], 7.0, 20.0).unwrap();
        let side = cylinder
            .faces()
            .unwrap()
            .into_iter()
            .find(|f| matches!(f.surface, Surface::Cylinder { .. }))
            .unwrap();
        let surfaces = [
            (
                "dome",
                samples::dome(k, [0.0, 0.0], [40.0, 20.0], 5.0, 8.0).unwrap(),
            ),
            ("tank", FreeSurface::of_face(&tank, skin.index).unwrap()),
            (
                "cylinder",
                FreeSurface::of_face(&cylinder, side.index).unwrap(),
            ),
        ];
        for (name, surface) in &surfaces {
            let [min, max] = surface.domain().unwrap();
            let mut worst: f64 = 0.0;
            for i in 1..8 {
                for j in 1..8 {
                    let u = min[0] + (max[0] - min[0]) * i as f64 / 8.0;
                    let v = min[1] + (max[1] - min[1]) * j as f64 / 8.0;
                    let (foot, normal) = surface.evaluate(u, v).unwrap();
                    for d in [0.0, 0.5, -0.5] {
                        let off = [0, 1, 2].map(|c| foot[c] + d * normal[c]);
                        let closest = surface.closest_point(off).unwrap();
                        let gap = distance(closest.point, foot);
                        worst = worst.max(gap);
                        assert!(
                            gap < 1e-5,
                            "{name} at ({u}, {v}) off {d}: {closest:?}, foot {foot:?}"
                        );
                        assert!(
                            (closest.distance - d.abs()).abs() < 1e-5,
                            "{name}: {closest:?}"
                        );
                        let (back, _) = surface.evaluate(closest.uv[0], closest.uv[1]).unwrap();
                        assert!(distance(back, closest.point) < 1e-6, "{name}: {closest:?}");
                    }
                }
            }
            println!("{name}: closest points within {worst:.1e} of the feet");
        }
    });
}

/// How a kernel operation went: built a manifold solid, refused, or came
/// back broken (not a manifold solid).
#[derive(Default, Debug)]
struct Tally {
    ok: u32,
    broken: u32,
    errors: std::collections::BTreeMap<String, u32>,
    ms: f64,
}

impl Tally {
    fn record(&mut self, result: Result<Brep, smlib::Error>, start: std::time::Instant) {
        self.ms += start.elapsed().as_secs_f64() * 1000.0;
        match result {
            Ok(brep) if brep.is_manifold().unwrap_or(false) => self.ok += 1,
            Ok(_) => self.broken += 1,
            Err(error) => {
                let key = match &error {
                    smlib::Error::Kernel { operation, status } => {
                        format!("{operation} {}", status.name().unwrap_or("?"))
                    }
                    other => other.to_string(),
                };
                *self.errors.entry(key).or_insert(0) += 1;
            }
        }
    }

    fn line(&self, what: &str) -> String {
        let n = self.ok + self.broken + self.errors.values().sum::<u32>();
        format!(
            "{what}: {}/{n} ok, {} broken, refused {:?}, {:.0} ms each",
            self.ok,
            self.broken,
            self.errors,
            self.ms / n.max(1) as f64
        )
    }
}

/// The indices of the faces that aren't planar.
fn curved(brep: &Brep) -> Vec<usize> {
    brep.faces()
        .unwrap()
        .iter()
        .filter(|f| f.normal.is_none())
        .map(|f| f.index)
        .collect()
}

/// The edges between face `a` and a face for which `other` holds.
fn edges_between(
    brep: &Brep,
    a: usize,
    other: impl Fn(&smlib::FaceGeometry) -> bool,
) -> Vec<usize> {
    let faces = brep.faces().unwrap();
    brep.edges()
        .unwrap()
        .iter()
        .filter(|e| {
            let f: Vec<usize> = e.faces.iter().flatten().map(|&f| f as usize).collect();
            f.len() == 2 && f.contains(&a) && f.iter().any(|&b| b != a && other(&faces[b]))
        })
        .map(|e| e.index)
        .collect()
}

/// Booleans, fillets, offsets and shells on free-form faces, varied, and how
/// often each fails and how: a survey, not a check, that takes about five
/// minutes in release. Run it by hand:
///
/// ```text
/// cargo test --release -p smlib --test freeform freeform_robustness -- --ignored --nocapture
/// ```
#[test]
#[ignore]
fn freeform_robustness() {
    on_kernel(|k| {
        use std::time::Instant;
        let dome = || {
            samples::domed_box(k, [-60.0, -40.0, -30.0], [120.0, 80.0, 42.0], 36.0)
                .unwrap()
                .0
        };
        let tank = || {
            let radii = [
                [27.0, 24.0],
                [41.4, 36.8],
                [45.0, 40.0],
                [40.5, 36.0],
                [24.8, 22.0],
            ];
            samples::lofted_tank(k, [-80.0, 0.0, 0.0], 160.0, &radii)
                .unwrap()
                .0
        };
        let rounded = || {
            samples::rounded_box(k, [-50.0, -30.0, -20.0], [100.0, 60.0, 40.0], 8.0)
                .unwrap()
                .brep
        };

        // Booleans: a pin through the dome on a grid, cut and joined.
        let (mut cut, mut join, mut common) =
            (Tally::default(), Tally::default(), Tally::default());
        for i in 0..9 {
            for j in 0..5 {
                for r in [3.0, 12.0] {
                    let (x, y) = (-52.0 + 13.0 * i as f64, -32.0 + 16.0 * j as f64);
                    let pin = || Brep::cylinder(k, [x, y, -60.0], r, 160.0).unwrap();
                    let t = Instant::now();
                    cut.record(dome().difference(pin()), t);
                    let short = Brep::cylinder(k, [x, y, 0.0], r, 60.0).unwrap();
                    let t = Instant::now();
                    join.record(dome().union(short), t);
                    let t = Instant::now();
                    common.record(dome().intersection(pin()), t);
                }
            }
        }
        println!("{}", cut.line("dome − pin, 90 places"));
        println!("{}", join.line("dome ∪ pin, 90 places"));
        println!("{}", common.line("dome ∩ pin, 90 places"));

        // Free-form against free-form: the tank cut from, joined to and
        // intersected with the dome, moved about.
        let (mut ff_cut, mut ff_join, mut ff_common) =
            (Tally::default(), Tally::default(), Tally::default());
        for dx in [-40.0, 0.0, 30.0] {
            for dz in [0.0, 20.0, 40.0, 60.0] {
                for turn in [0.0, 30.0] {
                    let moved = || {
                        let mut t = tank();
                        t.scale_about([0.0; 3], 0.5).unwrap();
                        t.rotate_about([0.0; 3], [0.0, 0.0, 1.0], turn).unwrap();
                        t.translate([dx, 0.0, dz]).unwrap();
                        t
                    };
                    let t = Instant::now();
                    ff_cut.record(dome().difference(moved()), t);
                    let t = Instant::now();
                    ff_join.record(dome().union(moved()), t);
                    let t = Instant::now();
                    ff_common.record(dome().intersection(moved()), t);
                }
            }
        }
        println!("{}", ff_cut.line("dome − tank, 24 placings"));
        println!("{}", ff_join.line("dome ∪ tank, 24 placings"));
        println!("{}", ff_common.line("dome ∩ tank, 24 placings"));
        let mut corners = Tally::default();
        for dz in [-10.0, -4.0, 0.0, 6.0] {
            for dx in [0.0, 20.0, 45.0] {
                let mut t = tank();
                t.scale_about([0.0; 3], 0.4).unwrap();
                t.translate([dx, 20.0, 20.0 + dz]).unwrap();
                let start = Instant::now();
                corners.record(rounded().difference(t), start);
            }
        }
        println!(
            "{}",
            corners.line("rounded box (N-sided corners) − tank, 12 placings")
        );

        // Fillets and chamfers: the dome's rim edges, one and all four, and
        // the edge a cut pin leaves on the dome.
        let mut rim = Tally::default();
        let mut rim_all = Tally::default();
        let mut rim_chamfer = Tally::default();
        for r in [1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0] {
            let d = dome();
            let top = curved(&d)[0];
            let edges = edges_between(&d, top, |f| f.normal.is_some());
            for &e in &edges {
                let t = Instant::now();
                rim.record(
                    dome().fillet(
                        &[EdgeBlend::round(EdgeQuery::Indices(vec![e]), r)],
                        &FilletOptions::default(),
                    ),
                    t,
                );
            }
            let t = Instant::now();
            rim_all.record(
                dome().fillet(
                    &[EdgeBlend::round(EdgeQuery::Indices(edges.clone()), r)],
                    &FilletOptions::default(),
                ),
                t,
            );
            let t = Instant::now();
            rim_chamfer.record(
                dome().chamfer(&EdgeQuery::Indices(edges.clone()), r, &Tolerance::default()),
                t,
            );
        }
        println!("{}", rim.line("fillet one dome rim edge, r 1–16 (28)"));
        println!("{}", rim_all.line("fillet all four rim edges, r 1–16 (7)"));
        println!(
            "{}",
            rim_chamfer.line("chamfer all four rim edges, 1–16 (7)")
        );
        let mut pin_rim = Tally::default();
        for (x, y) in [(0.0, 0.0), (20.0, 10.0), (-40.0, -25.0), (50.0, 30.0)] {
            for r in [0.5, 1.0, 2.0, 3.0] {
                let holed = dome()
                    .difference(Brep::cylinder(k, [x, y, -60.0], 6.0, 160.0).unwrap())
                    .unwrap();
                let top = curved(&holed)
                    .into_iter()
                    .find(|&f| matches!(holed.faces().unwrap()[f].surface, Surface::Nurbs { .. }))
                    .unwrap();
                let edges = edges_between(&holed, top, |f| {
                    matches!(f.surface, Surface::Cylinder { .. })
                });
                let t = Instant::now();
                pin_rim.record(
                    holed.fillet(
                        &[EdgeBlend::round(EdgeQuery::Indices(edges), r)],
                        &FilletOptions::default(),
                    ),
                    t,
                );
            }
        }
        println!(
            "{}",
            pin_rim.line("fillet a hole's rim in the dome, 4 places × r 0.5–3 (16)")
        );
        let mut tank_caps = Tally::default();
        for r in [1.0, 2.0, 4.0, 8.0, 12.0] {
            let t0 = tank();
            let skin = curved(&t0)[0];
            let edges = edges_between(&t0, skin, |f| f.normal.is_some());
            let t = Instant::now();
            tank_caps.record(
                tank().fillet(
                    &[EdgeBlend::round(EdgeQuery::Indices(edges), r)],
                    &FilletOptions::default(),
                ),
                t,
            );
        }
        println!(
            "{}",
            tank_caps.line("fillet the tank's cap rims, r 1–12 (5)")
        );

        // Offsets of the free-form face, out and in.
        let mut dome_offset = Tally::default();
        let mut tank_offset = Tally::default();
        for d in [-10.0, -5.0, -2.0, -0.5, 0.5, 2.0, 5.0, 10.0] {
            let dd = dome();
            let top = curved(&dd);
            let t = Instant::now();
            dome_offset.record(dd.offset_faces(&top, d), t);
            let tt = tank();
            let skin = curved(&tt);
            let t = Instant::now();
            tank_offset.record(tt.offset_faces(&skin, d), t);
        }
        println!("{}", dome_offset.line("offset the dome ±0.5–10 (8)"));
        println!("{}", tank_offset.line("offset the tank's skin ±0.5–10 (8)"));

        // Shells: hollowed inward, closed or open at a face.
        let (mut shell_dome, mut shell_tank, mut shell_rounded) =
            (Tally::default(), Tally::default(), Tally::default());
        for thickness in [1.0, 2.0, 3.0, 5.0, 8.0] {
            let d = dome();
            let top = curved(&d)[0];
            let floor = d
                .faces()
                .unwrap()
                .iter()
                .find(|f| f.normal.is_some_and(|n| n[2] < -0.9))
                .unwrap()
                .index;
            for open in [vec![], vec![top], vec![floor]] {
                let t = Instant::now();
                shell_dome.record(d.shell(thickness, &open).map(|s| s.brep), t);
            }
            let tt = tank();
            let cap = tt
                .faces()
                .unwrap()
                .iter()
                .find(|f| f.normal.is_some())
                .unwrap()
                .index;
            for open in [vec![], vec![cap]] {
                let t = Instant::now();
                shell_tank.record(tt.shell(thickness, &open).map(|s| s.brep), t);
            }
            let rr = rounded();
            let floor = rr
                .faces()
                .unwrap()
                .iter()
                .find(|f| f.normal.is_some_and(|n| n[2] < -0.9))
                .unwrap()
                .index;
            for open in [vec![], vec![floor]] {
                let t = Instant::now();
                shell_rounded.record(rr.shell(thickness, &open).map(|s| s.brep), t);
            }
        }
        println!(
            "{}",
            shell_dome.line("shell the domed box 1–8, closed / open at dome / at floor (15)")
        );
        println!(
            "{}",
            shell_tank.line("shell the tank 1–8, closed / open at a cap (10)")
        );
        println!(
            "{}",
            shell_rounded.line("shell the rounded box 1–8, closed / open at floor (10)")
        );
        // The same rounded box made by the kernel's own fillet (cylinders and
        // spheres), for comparison.
        let mut shell_analytic = Tally::default();
        for thickness in [1.0, 2.0, 3.0, 5.0, 8.0] {
            let b = Brep::cuboid(k, [-50.0, -30.0, -20.0], [100.0, 60.0, 40.0]).unwrap();
            let b = b
                .fillet(
                    &[EdgeBlend::round(
                        EdgeQuery::Sharp {
                            min_angle_deg: 30.0,
                        },
                        8.0,
                    )],
                    &FilletOptions::default(),
                )
                .unwrap();
            let floor = b
                .faces()
                .unwrap()
                .iter()
                .find(|f| f.normal.is_some_and(|n| n[2] < -0.9))
                .unwrap()
                .index;
            for open in [vec![], vec![floor]] {
                let t = Instant::now();
                shell_analytic.record(b.shell(thickness, &open).map(|s| s.brep), t);
            }
        }
        println!(
            "{}",
            shell_analytic.line("shell the analytic rounded box 1–8, closed / open at floor (10)")
        );
    });
}
