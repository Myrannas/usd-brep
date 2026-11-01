//! Networks of curves: a curve's net read back, and a loft pulled along
//! guides through the compatibility step (`FreeSurface::network`): through
//! its sections and its guides, closed round or open, a gap within the one
//! asked closed and a wider one refused by how far it is.

use std::f64::consts::{PI, TAU};
use std::panic::{self, AssertUnwindSafe};
use std::sync::OnceLock;
use std::sync::mpsc;

use smlib::nurbs::CurveNet;
use smlib::{Brep, Curve, FreeSurface, Kernel, NetworkError, Sheets};

type Job = Box<dyn FnOnce(&Kernel) + Send>;

/// Runs `f` on the one kernel thread the tests share.
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

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// A circle of `radius` about the z axis at `z`, from +x anticlockwise.
fn circle(z: f64, radius: f64) -> CurveNet {
    CurveNet::arc([0.0, 0.0, z], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], radius, 0.0, TAU)
}

/// The guide on the +x side of a duct from r 20 at z 0 through r 30 at z
/// 40 to r 20 at z 80, bulging past a skin's curve between them.
fn guide(shift: f64) -> CurveNet {
    let points = [[20.0, 0.0, 0.0], [29.0, 0.0, 20.0], [30.0, 0.0, 40.0], [29.0, 0.0, 60.0], [20.0, 0.0, 80.0]];
    let points = points.map(|[x, y, z]| [x + shift, y, z]);
    CurveNet::through_at(&points, &[0.0, 0.25, 0.5, 0.75, 1.0], 3).unwrap()
}

#[test]
fn a_curves_net_reads_back() {
    on_kernel(|k| {
        let ellipse = Curve::ellipse(k, [1.0, 2.0, 3.0], [26.0, 0.0, 0.0], [0.0, 18.0, 0.0]).unwrap();
        let net = ellipse.net().unwrap();
        assert!(net.is_valid() && net.rational());
        assert_eq!((net.degree, net.poles.len(), net.knots.len()), (2, 9, 12));
        let again = Curve::from_net(k, &net).unwrap();
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            assert!(distance(again.at(t).unwrap(), ellipse.at(t).unwrap()) < 1e-12);
            assert!(distance(net.point(t), ellipse.at(t).unwrap()) < 1e-12);
        }
        // An edge's curve: a cylinder's rim, a circle of 10.
        let pin = Brep::cylinder(k, [0.0, 0.0, 0.0], 10.0, 30.0).unwrap();
        let rims: Vec<usize> = pin.edges().unwrap().iter().filter(|edge| distance(edge.samples[0], edge.samples[edge.samples.len() - 1]) < 1e-9).map(|edge| edge.index).collect();
        assert!(!rims.is_empty());
        for rim in rims {
            let net = Curve::of_edge(&pin, rim).unwrap().net().unwrap();
            for i in 0..=40 {
                let p = net.point(net.domain()[0] + (net.domain()[1] - net.domain()[0]) * i as f64 / 40.0);
                assert!(((p[0] * p[0] + p[1] * p[1]).sqrt() - 10.0).abs() < 1e-9, "{p:?}");
            }
        }
    });
}

#[test]
fn a_loft_follows_its_guide_and_closes() {
    on_kernel(|k| {
        let sections = [circle(0.0, 20.0), circle(40.0, 30.0), circle(80.0, 20.0)];
        let network = FreeSurface::network(k, &sections, &[guide(0.0)], 2, 0.01).unwrap();
        assert!(network.widest < 1e-6, "{}", network.widest);
        let surface = &network.surface;
        // Through every section and along the guide.
        for section in &sections {
            for i in 0..36 {
                let p = section.point(i as f64 / 36.0);
                let near = surface.closest_point(p).unwrap();
                assert!(near.distance < 1e-3, "{p:?} is {} off", near.distance);
            }
        }
        let along = guide(0.0);
        for i in 0..=40 {
            let p = along.point(i as f64 / 40.0);
            assert!(surface.closest_point(p).unwrap().distance < 1e-3, "{p:?}");
        }
        // The guide pulls it out between the sections, where a loft through
        // the sections alone runs nearer the axis.
        let skin = FreeSurface::loft(k, &sections.iter().map(|s| Curve::from_net(k, s).unwrap()).collect::<Vec<_>>().iter().collect::<Vec<_>>(), 2).unwrap();
        let pulled = surface.closest_point([29.0, 0.0, 20.0]).unwrap().distance;
        let plain = skin.closest_point([29.0, 0.0, 20.0]).unwrap().distance;
        assert!(pulled < 1e-3 && plain > 0.5, "{pulled} {plain}");
        // Capped, it closes into a solid.
        let mut sheets = Sheets::new(k).unwrap();
        sheets.add(surface, None).unwrap();
        for end in [&sections[0], &sections[2]] {
            sheets.add_planar(&[&Curve::from_net(k, end).unwrap()], None).unwrap();
        }
        let (solid, closed) = sheets.stitch(0.01).unwrap();
        assert!(closed && solid.is_manifold().unwrap());
        let volume = solid.volume(1e-6).unwrap();
        assert!(volume > PI * 400.0 * 80.0 && volume < PI * 900.0 * 80.0, "{volume}");
    });
}

#[test]
fn a_guide_that_misses_says_by_how_much() {
    on_kernel(|k| {
        let sections = [circle(0.0, 20.0), circle(40.0, 30.0), circle(80.0, 20.0)];
        // Within the gap, it's closed; past it, refused.
        let near = FreeSurface::network(k, &sections, &[guide(0.008)], 2, 0.01).unwrap();
        assert!((near.widest - 0.008).abs() < 1e-6, "{}", near.widest);
        match FreeSurface::network(k, &sections, &[guide(1.4)], 2, 0.01) {
            Err(NetworkError::Misses { guide: 0, by, .. }) => assert!((by - 1.4).abs() < 1e-6, "{by}"),
            Err(other) => panic!("{other}"),
            Ok(_) => panic!("a guide 1.4 off built"),
        }
    });
}

#[test]
fn open_sections_take_a_guide_across_their_middles() {
    on_kernel(|k| {
        // Three arcs over the x axis, at y 0, 50 and 100, and a guide over
        // their tops rising in the middle.
        let arc = |y: f64, radius: f64| CurveNet::arc([0.0, y, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], radius, 0.0, PI);
        let sections = [arc(0.0, 30.0), arc(50.0, 30.0), arc(100.0, 30.0)];
        let top = CurveNet::through_at(&[[0.0, 0.0, 30.0], [0.0, 25.0, 34.0], [0.0, 50.0, 30.0], [0.0, 75.0, 34.0], [0.0, 100.0, 30.0]], &[0.0, 0.25, 0.5, 0.75, 1.0], 3).unwrap();
        let network = FreeSurface::network(k, &sections, std::slice::from_ref(&top), 2, 0.01).unwrap();
        for i in 0..=20 {
            let p = top.point(i as f64 / 20.0);
            assert!(network.surface.closest_point(p).unwrap().distance < 1e-3);
        }
        // Its edges stay on the sections' ends, along the ground.
        for (x, y) in [(30.0, 25.0), (-30.0, 75.0)] {
            let near = network.surface.closest_point([x, y, 0.0]).unwrap();
            assert!(near.distance < 0.05, "{x} {y}: {}", near.distance);
        }
        // A guide running the other way is turned round.
        let back = top.reversed();
        assert!(FreeSurface::network(k, &sections, &[back], 2, 0.01).is_ok());
    });
}

#[test]
fn a_rectangle_to_a_circle_along_a_guide_keeps_its_corners() {
    on_kernel(|k| {
        let rectangle = {
            let corners = [[30.0, -20.0, 0.0], [30.0, 20.0, 0.0], [-30.0, 20.0, 0.0], [-30.0, -20.0, 0.0], [30.0, -20.0, 0.0]];
            let lines: Vec<CurveNet> = corners.windows(2).map(|pair| CurveNet::line(pair[0], pair[1])).collect();
            CurveNet::joined(&lines).unwrap().0
        };
        assert_eq!(rectangle.corners(1e-6).len(), 3);
        let top = circle(80.0, 20.0);
        let rail = CurveNet::through_at(&[[30.0, 0.0, 0.0], [34.0, 0.0, 40.0], [20.0, 0.0, 80.0]], &[0.0, 0.5, 1.0], 2).unwrap();
        let network = FreeSurface::network(k, &[rectangle.clone(), top.clone()], &[rail], 1, 0.01).unwrap();
        let mut sheets = Sheets::new(k).unwrap();
        // Split at the rectangle's corners, at least.
        let pieces = sheets.add(&network.surface, None).unwrap();
        assert!(pieces >= 4, "{pieces}");
        for end in [&rectangle, &top] {
            sheets.add_planar(&[&Curve::from_net(k, end).unwrap()], None).unwrap();
        }
        let (solid, closed) = sheets.stitch(0.01).unwrap();
        assert!(closed && solid.is_manifold().unwrap());
    });
}

#[test]
fn closed_sections_take_two_guides_and_close_round() {
    on_kernel(|k| {
        let sections = [circle(0.0, 20.0), circle(40.0, 30.0), circle(80.0, 20.0)];
        let mirrored = {
            let one = guide(0.0);
            CurveNet { poles: one.poles.iter().map(|p| [-p[0], p[1], p[2]]).collect(), ..one }
        };
        let network = FreeSurface::network(k, &sections, &[guide(0.0), mirrored.clone()], 2, 0.01).unwrap();
        for i in 0..=20 {
            let p = mirrored.point(i as f64 / 20.0);
            assert!(network.surface.closest_point(p).unwrap().distance < 1e-3, "{p:?}");
        }
        let mut sheets = Sheets::new(k).unwrap();
        sheets.add(&network.surface, None).unwrap();
        for end in [&sections[0], &sections[2]] {
            sheets.add_planar(&[&Curve::from_net(k, end).unwrap()], None).unwrap();
        }
        let (solid, closed) = sheets.stitch(0.01).unwrap();
        assert!(closed && solid.is_manifold().unwrap());
    });
}
