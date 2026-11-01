//! NURBS curves as plain data, worked on in Rust: a curve's whole net
//! ([`CurveNet`]), as [`Curve::net`](crate::Curve::net) reads it back and
//! [`Curve::from_net`](crate::Curve::from_net) makes one of it, with what a
//! loft's sections and a network's curves need done to them before the
//! kernel sees them: evaluated, reversed, reparameterised, split at a knot,
//! raised in degree exactly, joined end to end, and fitted smoothly to a
//! function of their parameter.
//!
//! Rational curves are worked on in homogeneous coordinates, where they're
//! polynomial, so everything but a fit is exact up to rounding. The
//! algorithms are the textbook ones (Piegl and Tiller, The NURBS Book): the
//! basis functions (A2.1, A2.2), knot insertion (A5.1), and raising a degree
//! by solving for the poles in the higher-degree space, which holds the curve
//! exactly, at its Greville abscissae.

/// A NURBS curve's whole definition: its degree, its poles in Euclidean
/// coordinates, their weights (`None` for a non-rational curve) and its whole
/// knot vector, `poles + degree + 1` knots, clamped at both ends.
#[derive(Clone, Debug, PartialEq)]
pub struct CurveNet {
    /// Its degree, 1 or more.
    pub degree: u32,
    /// Its poles.
    pub poles: Vec<[f64; 3]>,
    /// Each pole's weight, or `None` for a non-rational curve.
    pub weights: Option<Vec<f64>>,
    /// Its whole knot vector.
    pub knots: Vec<f64>,
}

/// How near a fit keeps to what it fits, in model units, unless asked for
/// nearer: a micron in millimetres.
pub const FIT_TOLERANCE: f64 = 1e-3 * 1e-3;

/// The most spans a fit takes between two of its corners before giving up
/// on getting nearer.
const MOST_SPANS: usize = 512;

type H = [f64; 4];

fn hom(point: [f64; 3], weight: f64) -> H {
    [point[0] * weight, point[1] * weight, point[2] * weight, weight]
}

fn euclid(h: H) -> [f64; 3] {
    [h[0] / h[3], h[1] / h[3], h[2] / h[3]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// The knot span `t` lies in: the `i` with `knots[i] <= t < knots[i + 1]`,
/// the last span for the domain's end (A2.1).
fn span(knots: &[f64], degree: usize, count: usize, t: f64) -> usize {
    if t >= knots[count] {
        // The last non-empty span.
        let mut i = count - 1;
        while i > degree && knots[i] >= knots[count] {
            i -= 1;
        }
        return i;
    }
    if t <= knots[degree] {
        let mut i = degree;
        while i + 1 < count && knots[i + 1] <= knots[degree] {
            i += 1;
        }
        return i;
    }
    let (mut low, mut high) = (degree, count);
    let mut middle = (low + high) / 2;
    while t < knots[middle] || t >= knots[middle + 1] {
        if t < knots[middle] {
            high = middle;
        } else {
            low = middle;
        }
        middle = (low + high) / 2;
    }
    middle
}

/// The `degree + 1` basis functions not zero at `t` in span `i` (A2.2).
fn basis(knots: &[f64], degree: usize, i: usize, t: f64) -> Vec<f64> {
    let mut values = vec![0.0; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    values[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[i + 1 - j];
        right[j] = knots[i + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let denominator = right[r + 1] + left[j - r];
            let temp = if denominator == 0.0 { 0.0 } else { values[r] / denominator };
            values[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        values[j] = saved;
    }
    values
}

/// The poles in `knots` of `degree` that interpolate `values` at `sites`
/// (rising, each near its pole, as Greville abscissae and averaged knots
/// are), or `None` when the sites can't take them (Schoenberg and
/// Whitney). The collocation matrix is banded and totally positive, so it's
/// solved by banded elimination without pivoting (de Boor).
fn interpolate(knots: &[f64], degree: usize, sites: &[f64], values: Vec<H>) -> Option<Vec<H>> {
    let count = knots.len() - degree - 1;
    if sites.len() != count {
        return None;
    }
    let width = 2 * degree + 1;
    // band[r][k] is the matrix's (r, r - degree + k).
    let mut band = vec![vec![0.0; width]; count];
    for (row, &t) in sites.iter().enumerate() {
        let i = span(knots, degree, count, t);
        for (j, value) in basis(knots, degree, i, t).into_iter().enumerate() {
            let k = (i - degree + j + degree).checked_sub(row).filter(|&k| k < width)?;
            band[row][k] = value;
        }
    }
    let mut rhs = values;
    for c in 0..count {
        let pivot = band[c][degree];
        if pivot.abs() < 1e-14 {
            return None;
        }
        for r in c + 1..count.min(c + degree + 1) {
            let factor = band[r][c + degree - r] / pivot;
            if factor == 0.0 {
                continue;
            }
            for col in c..count.min(c + degree + 1) {
                band[r][col + degree - r] -= factor * band[c][col + degree - c];
            }
            let above = rhs[c];
            for (value, top) in rhs[r].iter_mut().zip(above) {
                *value -= factor * top;
            }
        }
    }
    let mut x = vec![[0.0; 4]; count];
    for r in (0..count).rev() {
        let mut value = rhs[r];
        for col in r + 1..count.min(r + degree + 1) {
            for (c, v) in value.iter_mut().enumerate() {
                *v -= band[r][col + degree - r] * x[col][c];
            }
        }
        x[r] = value.map(|v| v / band[r][degree]);
    }
    Some(x)
}

/// The Greville abscissae of a knot vector of `degree`: where each pole
/// pulls hardest.
fn greville(knots: &[f64], degree: usize) -> Vec<f64> {
    let count = knots.len() - degree - 1;
    (0..count).map(|i| knots[i + 1..=i + degree].iter().sum::<f64>() / degree as f64).collect()
}

impl CurveNet {
    /// A straight line from `a` to `b`, of degree 1 on 0 to 1.
    pub fn line(a: [f64; 3], b: [f64; 3]) -> CurveNet {
        CurveNet { degree: 1, poles: vec![a, b], weights: None, knots: vec![0.0, 0.0, 1.0, 1.0] }
    }

    /// An arc about `centre` of `radius`, from `start` radians to `start +
    /// sweep` measured from `x` towards `y` (unit vectors square to each
    /// other): rational quadratic, a span to each quarter turn or less, on 0
    /// to 1 in proportion to its angle. Exact.
    pub fn arc(centre: [f64; 3], x: [f64; 3], y: [f64; 3], radius: f64, start: f64, sweep: f64) -> CurveNet {
        let spans = ((sweep.abs() / std::f64::consts::FRAC_PI_2) - 1e-9).ceil().max(1.0) as usize;
        let step = sweep / spans as f64;
        let at = |angle: f64, reach: f64| {
            let (s, c) = angle.sin_cos();
            [0, 1, 2].map(|k| centre[k] + reach * (c * x[k] + s * y[k]))
        };
        let middle = (step / 2.0).cos();
        let mut poles = vec![at(start, radius)];
        let mut weights = vec![1.0];
        let mut knots = vec![0.0; 3];
        for k in 0..spans {
            let a = start + step * k as f64;
            poles.push(at(a + step / 2.0, radius / middle));
            weights.push(middle);
            poles.push(at(a + step, radius));
            weights.push(1.0);
            if k + 1 < spans {
                let t = (k + 1) as f64 / spans as f64;
                knots.extend([t, t]);
            }
        }
        knots.extend([1.0; 3]);
        CurveNet { degree: 2, poles, weights: Some(weights), knots }
    }

    /// Whether it has weights other than 1.
    pub fn rational(&self) -> bool {
        self.weights.as_ref().is_some_and(|weights| weights.iter().any(|&w| (w - 1.0).abs() > 1e-15))
    }

    fn weight(&self, i: usize) -> f64 {
        self.weights.as_ref().map_or(1.0, |weights| weights[i])
    }

    fn homogeneous(&self) -> Vec<H> {
        self.poles.iter().enumerate().map(|(i, &pole)| hom(pole, self.weight(i))).collect()
    }

    fn of_homogeneous(degree: u32, poles: Vec<H>, knots: Vec<f64>) -> CurveNet {
        let rational = poles.iter().any(|h| (h[3] - 1.0).abs() > 1e-15);
        CurveNet {
            degree,
            poles: poles.iter().map(|&h| euclid(h)).collect(),
            weights: rational.then(|| poles.iter().map(|h| h[3]).collect()),
            knots,
        }
    }

    /// Whether its degree, counts and knots fit together: rising, clamped,
    /// no inner knot more often than its degree, and finite.
    pub fn is_valid(&self) -> bool {
        let p = self.degree as usize;
        let n = self.poles.len();
        p >= 1
            && n > p
            && self.knots.len() == n + p + 1
            && self.weights.as_ref().is_none_or(|w| w.len() == n && w.iter().all(|&w| w > 0.0 && w.is_finite()))
            && self.knots.windows(2).all(|pair| pair[0] <= pair[1])
            && self.knots[..=p].iter().all(|&k| k == self.knots[0])
            && self.knots[n..].iter().all(|&k| k == self.knots[n + p])
            && self.knots[0] < self.knots[n + p]
            && self.poles.iter().flatten().chain(&self.knots).all(|v| v.is_finite())
    }

    /// The non-rational curve of `degree` through `points` at `at` (rising),
    /// its knots by averaging the parameters (Piegl and Tiller's 9.8);
    /// `None` for fewer points than the degree takes, or parameters that
    /// don't rise.
    pub fn through_at(points: &[[f64; 3]], at: &[f64], degree: u32) -> Option<CurveNet> {
        let p = degree as usize;
        let n = points.len();
        if p < 1 || n <= p || at.len() != n || at.windows(2).any(|pair| pair[1] <= pair[0]) {
            return None;
        }
        let mut knots = vec![at[0]; p + 1];
        for j in 1..n - p {
            knots.push(at[j..j + p].iter().sum::<f64>() / p as f64);
        }
        knots.extend(vec![at[n - 1]; p + 1]);
        let values = points.iter().map(|&point| hom(point, 1.0)).collect();
        let poles = interpolate(&knots, p, at, values)?;
        Some(CurveNet { degree, poles: poles.iter().map(|&h| euclid(h)).collect(), weights: None, knots })
    }

    /// Its parameter's range.
    pub fn domain(&self) -> [f64; 2] {
        [self.knots[0], self.knots[self.knots.len() - 1]]
    }

    fn homogeneous_at(&self, poles: &[H], t: f64) -> H {
        let p = self.degree as usize;
        let [low, high] = self.domain();
        let t = t.clamp(low, high);
        let i = span(&self.knots, p, poles.len(), t);
        let mut out = [0.0; 4];
        for (j, value) in basis(&self.knots, p, i, t).into_iter().enumerate() {
            for c in 0..4 {
                out[c] += value * poles[i - p + j][c];
            }
        }
        out
    }

    /// Its point at `t`, clamped to its domain.
    pub fn point(&self, t: f64) -> [f64; 3] {
        euclid(self.homogeneous_at(&self.homogeneous(), t))
    }

    /// Its point and first derivative at `t`.
    pub fn derivative(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        let p = self.degree as usize;
        let poles = self.homogeneous();
        let a = self.homogeneous_at(&poles, t);
        // The derivative curve, of one degree less on the inner knots.
        let n = poles.len();
        let mut d = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let gap = self.knots[i + p + 1] - self.knots[i + 1];
            let factor = if gap > 0.0 { p as f64 / gap } else { 0.0 };
            d.push([0, 1, 2, 3].map(|c| factor * (poles[i + 1][c] - poles[i][c])));
        }
        let da = if p == 1 {
            let knots = &self.knots[1..self.knots.len() - 1];
            let i = span(knots, 0, n - 1, t.clamp(self.domain()[0], self.domain()[1]));
            d[i]
        } else {
            let derived = CurveNet { degree: (p - 1) as u32, poles: vec![[0.0; 3]; n - 1], weights: None, knots: self.knots[1..self.knots.len() - 1].to_vec() };
            derived.homogeneous_at(&d, t)
        };
        let point = euclid(a);
        let w = a[3];
        let derivative = [0, 1, 2].map(|c| (da[c] - da[3] * point[c]) / w);
        (point, derivative)
    }

    /// The same curve run the other way, on the same domain.
    pub fn reversed(&self) -> CurveNet {
        let [low, high] = self.domain();
        let mut knots: Vec<f64> = self.knots.iter().rev().map(|&k| low + high - k).collect();
        knots[0] = low;
        *knots.last_mut().unwrap() = high;
        CurveNet {
            degree: self.degree,
            poles: self.poles.iter().rev().copied().collect(),
            weights: self.weights.as_ref().map(|w| w.iter().rev().copied().collect()),
            knots,
        }
    }

    /// The same curve on `[a, b]`, its knots mapped linearly.
    pub fn reparameterised(&self, [a, b]: [f64; 2]) -> CurveNet {
        let [low, high] = self.domain();
        let scale = (b - a) / (high - low);
        let mut knots: Vec<f64> = self.knots.iter().map(|&k| a + (k - low) * scale).collect();
        let p = self.degree as usize;
        let n = self.poles.len();
        knots[..=p].fill(a);
        knots[n..].fill(b);
        CurveNet { knots, ..self.clone() }
    }

    /// How many times `t` appears among its knots.
    pub fn multiplicity(&self, t: f64) -> usize {
        self.knots.iter().filter(|&&k| k == t).count()
    }

    /// The same curve with the knot `t` inserted `times` times more, but no
    /// more often than its degree in all (A5.1).
    pub fn with_knot(&self, t: f64, times: usize) -> CurveNet {
        let p = self.degree as usize;
        let [low, high] = self.domain();
        if t <= low || t >= high {
            return self.clone();
        }
        let already = self.multiplicity(t);
        let times = times.min(p.saturating_sub(already));
        let mut poles = self.homogeneous();
        let mut knots = self.knots.clone();
        for _ in 0..times {
            let n = poles.len();
            let k = span(&knots, p, n, t);
            let s = knots[..=k].iter().rev().take_while(|&&each| each == t).count();
            let mut next = Vec::with_capacity(n + 1);
            next.extend_from_slice(&poles[..=k - p]);
            for i in k - p + 1..=k - s {
                let alpha = (t - knots[i]) / (knots[i + p] - knots[i]);
                next.push([0, 1, 2, 3].map(|c| alpha * poles[i][c] + (1.0 - alpha) * poles[i - 1][c]));
            }
            next.extend_from_slice(&poles[k - s..]);
            knots.insert(k + 1, t);
            poles = next;
        }
        CurveNet::of_homogeneous(self.degree, poles, knots)
    }

    /// The curve cut in two at `t`, inside its domain: the piece before it
    /// and the piece after, each on its own part of the domain.
    pub fn split(&self, t: f64) -> (CurveNet, CurveNet) {
        let p = self.degree as usize;
        let full = self.with_knot(t, p);
        let poles = full.homogeneous();
        // The first index of t among the knots.
        let first = full.knots.iter().position(|&k| k == t).unwrap_or(p + 1);
        let mut left_knots = full.knots[..first + p].to_vec();
        left_knots.push(t);
        let left = CurveNet::of_homogeneous(self.degree, poles[..first].to_vec(), left_knots);
        let mut right_knots = vec![t];
        right_knots.extend_from_slice(&full.knots[first..]);
        let right = CurveNet::of_homogeneous(self.degree, poles[first - 1..].to_vec(), right_knots);
        (left, right)
    }

    /// The same curve raised to `degree`, each knot once more for each degree
    /// it's raised by, so its continuity is as it was. Exact: the poles are
    /// solved for at the Greville abscissae of the new knots, a space that
    /// holds the curve.
    pub fn elevated(&self, degree: u32) -> CurveNet {
        if degree <= self.degree {
            return self.clone();
        }
        let raise = (degree - self.degree) as usize;
        let mut knots = Vec::with_capacity(self.knots.len() * 2);
        for (i, &k) in self.knots.iter().enumerate() {
            knots.push(k);
            if i + 1 == self.knots.len() || self.knots[i + 1] != k {
                knots.extend(std::iter::repeat_n(k, raise));
            }
        }
        let p = degree as usize;
        let sites = greville(&knots, p);
        let poles = self.homogeneous();
        let values = sites.iter().map(|&t| self.homogeneous_at(&poles, t)).collect();
        let solved = interpolate(&knots, p, &sites, values).expect("the Greville abscissae take any degree's poles");
        let mut raised = CurveNet::of_homogeneous(degree, solved, knots);
        // The ends exactly as they were.
        let last = raised.poles.len() - 1;
        raised.poles[0] = self.poles[0];
        raised.poles[last] = self.poles[self.poles.len() - 1];
        raised
    }

    /// Joins `pieces`, each running on from the last's end, into one curve
    /// on `breaks` (one more than the pieces, rising): piece `k` on
    /// `breaks[k]` to `breaks[k + 1]`, all raised to the highest degree among
    /// them and their weights scaled to meet, each join a knot as often as
    /// the degree. Where a piece's end and the next's start differ, the join
    /// takes their middle. `None` for no pieces or breaks that don't fit.
    pub fn joined_at(pieces: &[CurveNet], breaks: &[f64]) -> Option<CurveNet> {
        if pieces.is_empty() || breaks.len() != pieces.len() + 1 || breaks.windows(2).any(|pair| pair[0] >= pair[1]) {
            return None;
        }
        let degree = pieces.iter().map(|piece| piece.degree).max()?;
        let p = degree as usize;
        let mut poles: Vec<H> = Vec::new();
        let mut knots = vec![breaks[0]; p + 1];
        for (k, piece) in pieces.iter().enumerate() {
            let piece = piece.elevated(degree).reparameterised([breaks[k], breaks[k + 1]]);
            let mut own = piece.homogeneous();
            if let Some(last) = poles.last_mut() {
                // Scaled so the first weight is the last's.
                let scale = last[3] / own[0][3];
                for h in &mut own {
                    for c in &mut *h {
                        *c *= scale;
                    }
                }
                let joint = [0, 1, 2].map(|c| (last[c] / last[3] + own[0][c] / own[0][3]) / 2.0);
                *last = hom(joint, last[3]);
                poles.extend_from_slice(&own[1..]);
                knots.extend(std::iter::repeat_n(breaks[k], p));
            } else {
                poles.extend_from_slice(&own);
            }
            knots.extend_from_slice(&piece.knots[p + 1..piece.knots.len() - p - 1]);
        }
        knots.extend(std::iter::repeat_n(breaks[pieces.len()], p + 1));
        Some(CurveNet::of_homogeneous(degree, poles, knots))
    }

    /// [`joined_at`](CurveNet::joined_at), each piece taking a share of 0 to
    /// 1 as long as its share of their length; and where each piece starts
    /// and the last ends.
    pub fn joined(pieces: &[CurveNet]) -> Option<(CurveNet, Vec<f64>)> {
        let lengths: Vec<f64> = pieces.iter().map(|piece| piece.length().max(1e-12)).collect();
        let total: f64 = lengths.iter().sum();
        let mut breaks = vec![0.0];
        for length in &lengths {
            breaks.push(breaks.last().unwrap() + length / total);
        }
        *breaks.last_mut()? = 1.0;
        let joined = CurveNet::joined_at(pieces, &breaks)?;
        Some((joined, breaks))
    }

    /// Its length, along 32 chords to each knot span.
    pub fn length(&self) -> f64 {
        let mut length = 0.0;
        let mut last = self.point(self.domain()[0]);
        for pair in self.knots.windows(2).filter(|pair| pair[0] < pair[1]) {
            for k in 1..=32 {
                let next = self.point(pair[0] + (pair[1] - pair[0]) * k as f64 / 32.0);
                length += norm(sub(next, last));
                last = next;
            }
        }
        length
    }

    /// The knots inside its domain where it may turn a corner: those as
    /// often as its degree, rising.
    pub fn breaks(&self) -> Vec<f64> {
        let p = self.degree as usize;
        let [low, high] = self.domain();
        let mut out: Vec<f64> = Vec::new();
        for &k in &self.knots {
            if k > low && k < high && self.multiplicity(k) >= p && out.last() != Some(&k) {
                out.push(k);
            }
        }
        out
    }

    /// The knots inside its domain where it turns a corner: [`breaks`]
    /// whose tangents either side are more than `angle` radians apart.
    pub fn corners(&self, angle: f64) -> Vec<f64> {
        self.breaks()
            .into_iter()
            .filter(|&t| {
                let (before, after) = self.split(t);
                let a = before.derivative(t).1;
                let b = after.derivative(t).1;
                let (la, lb) = (norm(a), norm(b));
                la == 0.0 || lb == 0.0 || (dot(a, b) / (la * lb)).clamp(-1.0, 1.0).acos() > angle
            })
            .collect()
    }

    /// The parameter of its point nearest `point`, from `start`, by Newton's
    /// method on the distance squared, clamped to its domain.
    pub fn nearest_from(&self, point: [f64; 3], start: f64) -> f64 {
        let [low, high] = self.domain();
        let mut t = start.clamp(low, high);
        for _ in 0..50 {
            let (at, d) = self.derivative(t);
            let off = sub(at, point);
            let slope = dot(d, d);
            if slope == 0.0 {
                break;
            }
            let step = dot(off, d) / slope;
            let next = (t - step).clamp(low, high);
            if (next - t).abs() < 1e-14 * (high - low) {
                t = next;
                break;
            }
            t = next;
        }
        t
    }

    /// A non-rational cubic through `curve(t)` for `t` from 0 to 1, within
    /// `tolerance` of it, turning corners only at `corners` and changing its
    /// bend sharply only there and at `joints` (each inside 0 to 1, rising),
    /// where it's C0 and C1 respectively and smooth elsewhere: uniform spans
    /// between those, doubled until the fit is near enough, fitted by
    /// interpolating at its Greville abscissae, its ends and corners on
    /// `curve`'s.
    pub fn fit(curve: impl Fn(f64) -> [f64; 3], corners: &[f64], joints: &[f64], tolerance: f64) -> CurveNet {
        let mut ends: Vec<(f64, usize)> = corners.iter().map(|&t| (t, 3)).chain(joints.iter().map(|&t| (t, 2))).filter(|&(t, _)| t > 1e-12 && t < 1.0 - 1e-12).collect();
        ends.sort_by(|a, b| a.0.total_cmp(&b.0));
        // One end at each place, a corner where either says so.
        ends.dedup_by(|later, earlier| {
            let same = (later.0 - earlier.0).abs() < 1e-12;
            if same {
                earlier.1 = earlier.1.max(later.1);
            }
            same
        });
        let mut spans = 4;
        loop {
            let piece = fit_spans(&curve, &ends, spans);
            let mut bounds = vec![0.0];
            bounds.extend(ends.iter().map(|&(t, _)| t));
            bounds.push(1.0);
            let near = bounds.windows(2).all(|pair| {
                (0..spans * 4).all(|k| {
                    let t = pair[0] + (pair[1] - pair[0]) * (k as f64 + 0.5) / (spans * 4) as f64;
                    norm(sub(piece.point(t), curve(t))) <= tolerance
                })
            });
            if near || spans >= MOST_SPANS {
                return piece;
            }
            spans *= 2;
        }
    }
}

/// A cubic through `curve` at its Greville abscissae, on knots of `spans`
/// uniform spans between each of `ends` (with how often each is a knot) and
/// 0 and 1.
fn fit_spans(curve: &impl Fn(f64) -> [f64; 3], ends: &[(f64, usize)], spans: usize) -> CurveNet {
    const P: usize = 3;
    let mut knots = vec![0.0; P + 1];
    let mut bounds = vec![(0.0, 0)];
    bounds.extend_from_slice(ends);
    bounds.push((1.0, P + 1));
    for pair in bounds.windows(2) {
        let [(a, _), (b, times)] = [pair[0], pair[1]];
        for k in 1..spans {
            knots.push(a + (b - a) * k as f64 / spans as f64);
        }
        knots.extend(std::iter::repeat_n(b, times));
    }
    let sites = greville(&knots, P);
    let values = sites.iter().map(|&t| hom(curve(t), 1.0)).collect();
    let poles = interpolate(&knots, P, &sites, values).expect("the Greville abscissae take a cubic's poles");
    CurveNet { degree: 3, poles: poles.iter().map(|&h| euclid(h)).collect(), weights: None, knots }
}

/// A function rising from `points[0]` through each of `points` to the
/// last, smooth and never falling where they rise (Fritsch and Carlson's
/// monotone cubic Hermite interpolation): how a network's curve is
/// reparameterised so it crosses the others at the parameters they share.
#[derive(Clone, Debug)]
pub struct Monotone {
    x: Vec<f64>,
    y: Vec<f64>,
    slopes: Vec<f64>,
}

impl Monotone {
    /// Through `points`, each `[x, y]`, both rising strictly; `None` for
    /// fewer than two or points that don't.
    pub fn through(points: &[[f64; 2]]) -> Option<Monotone> {
        if points.len() < 2 || points.windows(2).any(|pair| pair[1][0] <= pair[0][0] || pair[1][1] <= pair[0][1]) {
            return None;
        }
        let x: Vec<f64> = points.iter().map(|point| point[0]).collect();
        let y: Vec<f64> = points.iter().map(|point| point[1]).collect();
        let n = x.len();
        let secants: Vec<f64> = (0..n - 1).map(|k| (y[k + 1] - y[k]) / (x[k + 1] - x[k])).collect();
        let mut slopes = vec![0.0; n];
        slopes[0] = secants[0];
        slopes[n - 1] = secants[n - 2];
        for k in 1..n - 1 {
            slopes[k] = (secants[k - 1] + secants[k]) / 2.0;
        }
        for k in 0..n - 1 {
            let (a, b) = (slopes[k] / secants[k], slopes[k + 1] / secants[k]);
            let reach = a * a + b * b;
            if reach > 9.0 {
                let tau = 3.0 / reach.sqrt();
                slopes[k] = tau * a * secants[k];
                slopes[k + 1] = tau * b * secants[k];
            }
        }
        Some(Monotone { x, y, slopes })
    }

    /// Its value at `t`, clamped to its first and last x.
    pub fn at(&self, t: f64) -> f64 {
        let n = self.x.len();
        let t = t.clamp(self.x[0], self.x[n - 1]);
        let k = self.x.windows(2).position(|pair| t <= pair[1]).unwrap_or(n - 2);
        let h = self.x[k + 1] - self.x[k];
        let s = (t - self.x[k]) / h;
        let (s2, s3) = (s * s, s * s * s);
        (2.0 * s3 - 3.0 * s2 + 1.0) * self.y[k] + (s3 - 2.0 * s2 + s) * h * self.slopes[k] + (-2.0 * s3 + 3.0 * s2) * self.y[k + 1] + (s3 - s2) * h * self.slopes[k + 1]
    }

    /// The `t` it takes `value` at, by bisection.
    pub fn inverse(&self, value: f64) -> f64 {
        let n = self.x.len();
        let (mut low, mut high) = (self.x[0], self.x[n - 1]);
        for _ in 0..100 {
            let middle = (low + high) / 2.0;
            if self.at(middle) < value {
                low = middle;
            } else {
                high = middle;
            }
        }
        (low + high) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: [f64; 3], b: [f64; 3]) -> f64 {
        norm(sub(a, b))
    }

    fn quarter_circle_and_line() -> Vec<CurveNet> {
        vec![
            CurveNet::line([0.0, 0.0, 0.0], [10.0, 0.0, 0.0]),
            CurveNet::arc([10.0, 5.0, 0.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0], 5.0, 0.0, std::f64::consts::PI),
            CurveNet::line([10.0, 10.0, 0.0], [0.0, 10.0, 0.0]),
        ]
    }

    #[test]
    fn an_arc_lies_on_its_circle() {
        let arc = CurveNet::arc([1.0, 2.0, 3.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 7.0, 0.3, 4.0);
        assert!(arc.is_valid());
        for k in 0..=50 {
            let p = arc.point(k as f64 / 50.0);
            assert!((near(p, [1.0, 2.0, 3.0]) - 7.0).abs() < 1e-12, "{p:?}");
        }
        let end = arc.point(1.0);
        let expected = [1.0 + 7.0 * 4.3f64.cos(), 2.0, 3.0 + 7.0 * 4.3f64.sin()];
        assert!(near(end, expected) < 1e-12);
    }

    #[test]
    fn elevating_keeps_the_curve() {
        for piece in quarter_circle_and_line() {
            let raised = piece.elevated(4);
            assert!(raised.is_valid());
            assert_eq!(raised.degree, 4);
            for k in 0..=40 {
                let t = k as f64 / 40.0;
                assert!(near(raised.point(t), piece.point(t)) < 1e-11, "{t}");
            }
        }
    }

    #[test]
    fn joining_keeps_each_piece_and_its_corners() {
        let pieces = quarter_circle_and_line();
        let (joined, breaks) = CurveNet::joined(&pieces).unwrap();
        assert!(joined.is_valid());
        assert_eq!(breaks.len(), 4);
        for (k, piece) in pieces.iter().enumerate() {
            for s in 0..=10 {
                let s = s as f64 / 10.0;
                let t = breaks[k] + (breaks[k + 1] - breaks[k]) * s;
                assert!(near(joined.point(t), piece.point(s)) < 1e-10);
            }
        }
        // The line meets the arc square, a corner, and the arc meets the
        // line along it.
        let corners = joined.corners(1e-6);
        assert_eq!(corners.len(), 0, "the half circle meets both lines tangentially: {corners:?}");
        assert_eq!(joined.breaks().len(), 2 + 1);
    }

    #[test]
    fn splitting_and_reversing_keep_the_curve() {
        let (joined, _) = CurveNet::joined(&quarter_circle_and_line()).unwrap();
        let (a, b) = joined.split(0.37);
        assert!(a.is_valid() && b.is_valid());
        for k in 0..=20 {
            let t = 0.37 * k as f64 / 20.0;
            assert!(near(a.point(t), joined.point(t)) < 1e-10);
            let t = 0.37 + 0.63 * k as f64 / 20.0;
            assert!(near(b.point(t), joined.point(t)) < 1e-10);
        }
        let turned = joined.reversed();
        for k in 0..=20 {
            let t = k as f64 / 20.0;
            assert!(near(turned.point(t), joined.point(1.0 - t)) < 1e-10);
        }
    }

    #[test]
    fn a_fit_follows_its_curve_and_keeps_its_corners() {
        let square = |t: f64| {
            let s = (t * 4.0).min(3.999_999_999);
            let k = s.floor();
            let f = s - k;
            match k as u32 {
                0 => [f, 0.0, 0.0],
                1 => [1.0, f, 0.0],
                2 => [1.0 - f, 1.0, 0.0],
                _ => [0.0, 1.0 - f, 0.0],
            }
        };
        let fit = CurveNet::fit(square, &[0.25, 0.5, 0.75], &[], 1e-9);
        assert!(fit.is_valid() && !fit.rational());
        assert_eq!(fit.corners(1e-3), [0.25, 0.5, 0.75]);
        for k in 0..=100 {
            let t = k as f64 / 100.0;
            assert!(near(fit.point(t), square(t)) < 1e-8);
        }
        let circle = CurveNet::arc([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 20.0, 0.0, std::f64::consts::TAU);
        let fit = CurveNet::fit(|t| circle.point(t), &[], &circle.breaks(), 1e-6);
        let worst = (0..=2000).map(|k| near(fit.point(k as f64 / 2000.0), circle.point(k as f64 / 2000.0))).fold(0.0, f64::max);
        assert!(worst < 2e-6, "{worst} with {} poles", fit.poles.len());
    }

    #[test]
    fn a_monotone_map_goes_through_its_points_and_never_falls() {
        let map = Monotone::through(&[[0.0, 0.0], [0.3, 0.1], [0.5, 0.8], [1.0, 1.0]]).unwrap();
        assert!((map.at(0.3) - 0.1).abs() < 1e-15 && (map.at(0.5) - 0.8).abs() < 1e-15);
        let mut last = -1.0;
        for k in 0..=1000 {
            let value = map.at(k as f64 / 1000.0);
            assert!(value >= last);
            last = value;
        }
        assert!((map.at(map.inverse(0.42)) - 0.42).abs() < 1e-12);
    }
}
