//! Times `Brep::offset_face` / `offset_faces` on the test cases, against
//! rebuilding each part at its new size from scratch.
//!
//! ```sh
//! cargo run --release -p smlib --example offset_face            # timings
//! cargo run --release -p smlib --example offset_face -- repeat 500
//! ```
//!
//! `repeat N` offsets the chamfered box N times each way, once to a valid
//! height and once to one where the chamfer would vanish, for watching memory
//! with `/usr/bin/time -l`.

use std::time::{Duration, Instant};

use smlib::{
    Brep, EdgeBlend, EdgeQuery, Error, FaceQuery, FilletOptions, Kernel, SurfaceKind, Tolerance,
};

const L: f64 = 100.0;
const W: f64 = 60.0;
const H: f64 = 40.0;
const RUNS: usize = 25;

fn block(k: &Kernel, h: f64) -> Result<Brep, Error> {
    Brep::cuboid(k, [-L / 2.0, -W / 2.0, -H / 2.0], [L, W, h])
}

fn top(h: f64) -> FaceQuery {
    FaceQuery::OnPlane {
        point: [0.0, 0.0, -H / 2.0 + h],
        normal: [0.0, 0.0, 1.0],
    }
}

fn top_edges(h: f64) -> EdgeQuery {
    EdgeQuery::MidpointOnPlane {
        point: [0.0, 0.0, -H / 2.0 + h],
        normal: [0.0, 0.0, 1.0],
    }
}

fn plain(k: &Kernel, h: f64) -> Result<Brep, Error> {
    block(k, h)
}

fn chamfered(k: &Kernel, h: f64) -> Result<Brep, Error> {
    let legs = 2.0 + (h - H);
    block(k, h)?.chamfer(&top_edges(h), legs * 2f64.sqrt(), &Tolerance::default())
}

fn filleted(k: &Kernel, h: f64) -> Result<Brep, Error> {
    block(k, h)?.fillet(
        &[EdgeBlend::round(top_edges(h), 8.0)],
        &FilletOptions::default(),
    )
}

fn drafted(k: &Kernel, h: f64) -> Result<Brep, Error> {
    let (x, y, z) = (L / 2.0, W / 2.0, -H / 2.0);
    Brep::tapered_prism(k, &[[-x, -y, z], [x, -y, z], [x, y, z], [-x, y, z]], h, 5.0)
}

fn rim_rounded(k: &Kernel, h: f64) -> Result<Brep, Error> {
    let part = block(k, h)?.difference(Brep::cylinder(k, [0.0, 0.0, -H / 2.0], 10.0, h)?)?;
    let rim = EdgeQuery::OnCircle {
        center: [0.0, 0.0, -H / 2.0 + h],
        axis: [0.0, 0.0, 1.0],
        radius: 10.0,
    };
    part.fillet(&[EdgeBlend::round(rim, 3.0)], &FilletOptions::default())
}

fn holed(k: &Kernel, h: f64, r: f64) -> Result<Brep, Error> {
    block(k, h)?.difference(Brep::cylinder(k, [0.0, 0.0, -H / 2.0], r, h)?)
}

fn plus_x() -> FaceQuery {
    FaceQuery::Facing {
        direction: [1.0, 0.0, 0.0],
    }
    .and(FaceQuery::Kind(SurfaceKind::Plane))
}

fn rounded_span(k: &Kernel, length: f64, h: f64) -> Result<Brep, Error> {
    Brep::cuboid(k, [-L / 2.0, -W / 2.0, -H / 2.0], [length, W, h])?.fillet(
        &[EdgeBlend::round(top_edges(h), 8.0)],
        &FilletOptions::default(),
    )
}

fn fully_rounded(k: &Kernel, h: f64) -> Result<Brep, Error> {
    block(k, h)?.fillet(
        &[EdgeBlend::round(EdgeQuery::All, 8.0)],
        &FilletOptions::default(),
    )
}

/// The faces `query` picks, offset: Brep::find_faces then Brep::offset_faces.
fn offset_where(
    brep: Brep,
    query: &FaceQuery,
    distance: f64,
    tol: &Tolerance,
) -> Result<Brep, Error> {
    let faces = brep.find_faces(query, tol)?;
    brep.offset_faces(&faces, distance)
}

type Part = fn(&Kernel, f64) -> Result<Brep, Error>;
type Case = (
    &'static str,
    Box<dyn Fn(&Kernel) -> Result<Brep, Error>>,
    Box<dyn Fn(Brep) -> Result<Brep, Error>>,
    Box<dyn Fn(&Kernel) -> Result<Brep, Error>>,
);

struct Stats {
    min: Duration,
    median: Duration,
    max: Duration,
}

fn stats(mut samples: Vec<Duration>) -> Stats {
    samples.sort();
    Stats {
        min: samples[0],
        median: samples[samples.len() / 2],
        max: samples[samples.len() - 1],
    }
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1e3
}

fn main() -> Result<(), Error> {
    let kernel = Kernel::init()?;
    let k = &kernel;
    let tol = Tolerance::default();
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("repeat") {
        let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(500);
        let part = chamfered(k, H)?;
        let (mut ok, mut refused) = (0, 0);
        for _ in 0..n {
            for d in [5.0, -10.0] {
                match offset_where(part.try_clone()?, &top(H), d, &tol) {
                    Ok(_) => ok += 1,
                    Err(_) => refused += 1,
                }
            }
        }
        println!("{n} rounds: {ok} offsets, {refused} refused");
        return Ok(());
    }

    let cases: [(&str, Part, f64); 9] = [
        ("plain box", plain, 10.0),
        ("plain box", plain, -10.0),
        ("chamfered box", chamfered, 10.0),
        ("chamfered box", chamfered, -1.0),
        ("filleted box", filleted, 10.0),
        ("filleted box", filleted, -10.0),
        ("drafted block", drafted, 10.0),
        ("drafted block", drafted, -10.0),
        ("holed box, rounded rim", rim_rounded, 10.0),
    ];
    println!(
        "| case | d | offset_face min / median / max (ms) | rebuild at new height, median (ms) |"
    );
    println!("|---|---|---|---|");
    for (name, part, d) in cases {
        let base = part(k, H)?;
        // Warm up once, then time the offset alone on fresh copies.
        offset_where(base.try_clone()?, &top(H), d, &tol)?;
        let mut offsets = Vec::with_capacity(RUNS);
        let mut rebuilds = Vec::with_capacity(RUNS);
        for _ in 0..RUNS {
            let copy = base.try_clone()?;
            let start = Instant::now();
            let moved = offset_where(copy, &top(H), d, &tol)?;
            offsets.push(start.elapsed());
            drop(moved);
            let start = Instant::now();
            let rebuilt = part(k, H + d)?;
            rebuilds.push(start.elapsed());
            drop(rebuilt);
        }
        let o = stats(offsets);
        let r = stats(rebuilds);
        println!(
            "| {name} | {d:+} | {:.2} / {:.2} / {:.2} | {:.2} |",
            ms(o.min),
            ms(o.median),
            ms(o.max),
            ms(r.median)
        );
    }
    let tol = Tolerance::default();
    let cases: Vec<Case> = vec![
        (
            "filleted box, +X wall +10",
            Box::new(|k| rounded_span(k, L, H)),
            Box::new(move |b| offset_where(b, &plus_x(), 10.0, &tol)),
            Box::new(|k| rounded_span(k, L + 10.0, H)),
        ),
        (
            "filleted box, top and +X wall +10",
            Box::new(|k| rounded_span(k, L, H)),
            Box::new(move |b| offset_where(b, &top(H).or(plus_x()), 10.0, &tol)),
            Box::new(|k| rounded_span(k, L + 10.0, H + 10.0)),
        ),
        (
            "fully rounded box, top +5",
            Box::new(|k| fully_rounded(k, H)),
            Box::new(move |b| offset_where(b, &top(H), 5.0, &tol)),
            Box::new(|k| fully_rounded(k, H + 5.0)),
        ),
        (
            "plain box, all six faces -5",
            Box::new(|k| block(k, H)),
            Box::new(move |b| offset_where(b, &FaceQuery::All, -5.0, &tol)),
            Box::new(|k| Brep::cuboid(k, [-45.0, -25.0, -15.0], [90.0, 50.0, 30.0])),
        ),
        (
            "hole wall +2",
            Box::new(|k| holed(k, H, 10.0)),
            Box::new(move |b| offset_where(b, &FaceQuery::Kind(SurfaceKind::Cylinder), 2.0, &tol)),
            Box::new(|k| holed(k, H, 8.0)),
        ),
        (
            "frustum side +2",
            Box::new(|k| Brep::cone(k, [0.0; 3], 20.0, 10.0, 30.0)),
            Box::new(move |b| offset_where(b, &FaceQuery::Kind(SurfaceKind::Cone), 2.0, &tol)),
            Box::new(|k| {
                let grow = 2.0 / (10.0f64 / 30.0).atan().cos();
                Brep::cone(k, [0.0; 3], 20.0 + grow, 10.0 + grow, 30.0)
            }),
        ),
        (
            "torus +1",
            Box::new(|k| Brep::torus(k, [0.0; 3], 30.0, 10.0)),
            Box::new(move |b| offset_where(b, &FaceQuery::All, 1.0, &tol)),
            Box::new(|k| Brep::torus(k, [0.0; 3], 30.0, 11.0)),
        ),
    ];
    for (name, make, offset, rebuild) in &cases {
        let base = make(k)?;
        offset(base.try_clone()?)?;
        let mut offsets = Vec::with_capacity(RUNS);
        let mut rebuilds = Vec::with_capacity(RUNS);
        for _ in 0..RUNS {
            let copy = base.try_clone()?;
            let start = Instant::now();
            let moved = offset(copy)?;
            offsets.push(start.elapsed());
            drop(moved);
            let start = Instant::now();
            let rebuilt = rebuild(k)?;
            rebuilds.push(start.elapsed());
            drop(rebuilt);
        }
        let o = stats(offsets);
        let r = stats(rebuilds);
        println!(
            "| {name} | | {:.2} / {:.2} / {:.2} | {:.2} |",
            ms(o.min),
            ms(o.median),
            ms(o.max),
            ms(r.median)
        );
    }
    Ok(())
}
