//! A plane spline's own arithmetic, for checking a
//! [`PlaneCurve::Spline`](crate::PlaneCurve::Spline) before the kernel sees
//! it: whether it can be made, which way it runs at its ends, where it turns a
//! corner inside and how tightly it bends. A clamped, non-rational B-spline,
//! evaluated by its basis functions and their derivatives (Piegl and Tiller's
//! A2.1 and A2.3).

/// How many points each knot span is sampled at when looking for its
/// tightest bend, before the tightest is refined.
const BEND_SAMPLES: usize = 32;

/// Why a clamped non-rational B-spline of `degree` on `poles` and the whole
/// knot vector `knots` can't be made, if it can't.
pub(crate) fn problem(degree: u32, poles: &[[f64; 2]], knots: &[f64]) -> Option<&'static str> {
    let p = degree as usize;
    if !(1..=MAX_DEGREE).contains(&degree) {
        return Some("a spline's degree must be from 1 to 9");
    }
    if poles.len() < p + 1 {
        return Some("a spline needs more poles than its degree");
    }
    if knots.len() != poles.len() + p + 1 {
        return Some("a spline needs as many knots as its poles and its degree and one more");
    }
    if !poles.iter().flatten().chain(knots).all(|c| c.is_finite()) {
        return Some("a spline's poles and knots must be finite");
    }
    if knots.windows(2).any(|pair| pair[1] < pair[0]) {
        return Some("a spline's knots must rise");
    }
    let (first, last) = (knots[0], knots[knots.len() - 1]);
    let clamped = knots[..=p].iter().all(|&k| k == first) && knots[knots.len() - 1 - p..].iter().all(|&k| k == last);
    if !clamped || first == last {
        return Some("a spline's knots must be clamped: its first and last degree + 1 alike");
    }
    if multiplicities(knots).any(|(knot, times)| knot != first && knot != last && times > p) {
        return Some("a spline's knots inside it can't repeat more than its degree");
    }
    if poles.iter().all(|&pole| pole == poles[0]) {
        return Some("a spline's poles can't all be one point");
    }
    None
}

/// The highest degree a spline may have.
const MAX_DEGREE: u32 = 9;

/// Each distinct knot and how many times it comes.
fn multiplicities(knots: &[f64]) -> impl Iterator<Item = (f64, usize)> + '_ {
    let mut at = 0;
    std::iter::from_fn(move || {
        let knot = *knots.get(at)?;
        let times = knots[at..].iter().take_while(|&&k| k == knot).count();
        at += times;
        Some((knot, times))
    })
}

/// The unit direction a clamped spline runs in where it starts, or where it
/// ends: towards the first pole that differs from the first (from the last
/// pole that differs from the last), which is where a clamped B-spline's
/// lowest derivative that doesn't vanish points.
pub(crate) fn end_tangent(poles: &[[f64; 2]], at_end: bool) -> [f64; 2] {
    let ordered: Box<dyn Iterator<Item = &[f64; 2]>> = if at_end {
        Box::new(poles.iter().rev())
    } else {
        Box::new(poles.iter())
    };
    let mut ordered = ordered.peekable();
    let Some(&&end) = ordered.peek() else {
        return [0.0, 0.0];
    };
    for &[x, y] in ordered {
        let (dx, dy) = (x - end[0], y - end[1]);
        let length = dx.hypot(dy);
        if length > 1e-12 {
            let sign = if at_end { -1.0 } else { 1.0 };
            return [sign * dx / length, sign * dy / length];
        }
    }
    [0.0, 0.0]
}

/// A spline to evaluate: a [`problem`]-free clamped non-rational B-spline.
pub(crate) struct Spline<'a> {
    degree: usize,
    poles: &'a [[f64; 2]],
    knots: &'a [f64],
}

impl<'a> Spline<'a> {
    pub(crate) fn new(degree: u32, poles: &'a [[f64; 2]], knots: &'a [f64]) -> Spline<'a> {
        Spline {
            degree: degree as usize,
            poles,
            knots,
        }
    }

    /// The knot span `t` lies in: the `i` with `knots[i] <= t < knots[i + 1]`,
    /// the last span that isn't empty at the very end.
    fn span(&self, t: f64) -> usize {
        let (p, n) = (self.degree, self.poles.len());
        if t >= self.knots[n] {
            return (p..n).rev().find(|&i| self.knots[i] < self.knots[i + 1]).unwrap_or(p);
        }
        (p..n).rev().find(|&i| self.knots[i] <= t).unwrap_or(p)
    }

    /// The point at `t` and its first two derivatives along `t`.
    pub(crate) fn derivatives(&self, t: f64) -> [[f64; 2]; 3] {
        let span = self.span(t);
        let basis = basis_derivatives(span, t, self.degree, self.knots, 2);
        let mut out = [[0.0; 2]; 3];
        for (k, row) in basis.iter().enumerate() {
            for (j, &weight) in row.iter().enumerate() {
                let pole = self.poles[span - self.degree + j];
                out[k][0] += weight * pole[0];
                out[k][1] += weight * pole[1];
            }
        }
        out
    }

    /// Its knots where it may turn a corner: those inside it repeated as many
    /// times as its degree, where it's only continuous.
    fn joints(&self) -> impl Iterator<Item = f64> + '_ {
        let (first, last) = (self.knots[0], self.knots[self.knots.len() - 1]);
        multiplicities(self.knots)
            .filter(move |&(knot, times)| knot != first && knot != last && times >= self.degree)
            .map(|(knot, _)| knot)
    }

    /// Whether it turns a corner inside: runs on at one of its
    /// [`joints`](Spline::joints) in a different direction, by more than
    /// `turn_rad`, or stops dead there.
    pub(crate) fn has_corner(&self, turn_rad: f64) -> bool {
        self.joints().any(|knot| {
            let before = self.leaving(knot, false);
            let after = self.leaving(knot, true);
            match (before, after) {
                (Some([x1, y1]), Some([x2, y2])) => (x1 * y2 - y1 * x2).abs().atan2(x1 * x2 + y1 * y2) > turn_rad,
                _ => true,
            }
        })
    }

    /// Which way it runs into (`after` false) or out of a knot, from the
    /// span on that side: its first derivative there, unit.
    fn leaving(&self, knot: f64, after: bool) -> Option<[f64; 2]> {
        let span = if after {
            self.span(knot)
        } else {
            (self.degree..self.poles.len()).rev().find(|&i| self.knots[i] < knot)?
        };
        let basis = basis_derivatives(span, knot, self.degree, self.knots, 1);
        let mut d = [0.0; 2];
        for (j, &weight) in basis[1].iter().enumerate() {
            let pole = self.poles[span - self.degree + j];
            d[0] += weight * pole[0];
            d[1] += weight * pole[1];
        }
        let length = d[0].hypot(d[1]);
        (length > 1e-12).then(|| [d[0] / length, d[1] / length])
    }

    /// Its signed curvature at `t`, positive where it turns anticlockwise,
    /// to its left; infinite where it stops dead.
    fn curvature(&self, t: f64) -> f64 {
        let [_, [x1, y1], [x2, y2]] = self.derivatives(t);
        let speed = x1.hypot(y1);
        if speed < 1e-12 {
            return f64::INFINITY;
        }
        (x1 * y2 - y1 * x2) / speed.powi(3)
    }

    /// The most it bends to its left and to its right, each as a curvature,
    /// 0 when it doesn't bend that way: one over the tightest bend's radius.
    /// Each knot span is sampled and the tightest sample refined, so it can
    /// fall a hair short of the true greatest.
    pub(crate) fn sharpest_bends(&self) -> [f64; 2] {
        let mut sharpest = [0.0f64; 2];
        let p = self.degree;
        for i in p..self.poles.len() {
            let (a, b) = (self.knots[i], self.knots[i + 1]);
            if a >= b {
                continue;
            }
            let at = |k: usize| a + (b - a) * k as f64 / BEND_SAMPLES as f64;
            for side in 0..2 {
                // Left bends are positive curvature, right ones negative.
                let bend = |t: f64| {
                    let curvature = self.curvature(t);
                    if curvature.is_infinite() {
                        return f64::INFINITY;
                    }
                    if side == 0 { curvature } else { -curvature }
                };
                let best = (0..=BEND_SAMPLES).max_by(|&j, &k| bend(at(j)).total_cmp(&bend(at(k)))).unwrap_or(0);
                let mut most = bend(at(best));
                if most.is_finite() {
                    let low = at(best.saturating_sub(1));
                    let high = at((best + 1).min(BEND_SAMPLES));
                    most = most.max(golden_max(&bend, low, high));
                }
                sharpest[side] = sharpest[side].max(most);
            }
        }
        sharpest
    }
}

/// The greatest of `f` over `[low, high]` near one peak, by golden-section
/// search.
fn golden_max(f: &impl Fn(f64) -> f64, mut low: f64, mut high: f64) -> f64 {
    let ratio = (5f64.sqrt() - 1.0) / 2.0;
    let mut a = high - ratio * (high - low);
    let mut b = low + ratio * (high - low);
    let (mut fa, mut fb) = (f(a), f(b));
    for _ in 0..60 {
        if fa > fb {
            high = b;
            b = a;
            fb = fa;
            a = high - ratio * (high - low);
            fa = f(a);
        } else {
            low = a;
            a = b;
            fa = fb;
            b = low + ratio * (high - low);
            fb = f(b);
        }
    }
    fa.max(fb)
}

/// The non-zero basis functions at `t` in `span` and their derivatives up to
/// `order`: `out[k][j]` is the `k`th derivative of `N(span - degree + j)`
/// (Piegl and Tiller's A2.3).
fn basis_derivatives(span: usize, t: f64, degree: usize, knots: &[f64], order: usize) -> Vec<Vec<f64>> {
    let p = degree;
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let order = order.min(p);
    let mut out = vec![vec![0.0; p + 1]; order + 1];
    for j in 0..=p {
        out[0][j] = ndu[j][p];
    }
    let mut a = vec![vec![0.0; p + 1]; 2];
    for r in 0..=p {
        let (mut s1, mut s2) = (0, 1);
        a[0][0] = 1.0;
        for k in 1..=order {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = p - k;
            if r >= k {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk as usize];
                d = a[s2][0] * ndu[rk as usize][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk as isize { k - 1 } else { p - r };
            for j in j1..=j2 {
                let index = (rk + j as isize) as usize;
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][index];
                d += a[s2][j] * ndu[index][pk];
            }
            if r <= pk {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            out[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut factor = p as f64;
    for k in 1..=order {
        for value in out[k].iter_mut() {
            *value *= factor;
        }
        factor *= (p - k) as f64;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The parabola y = x² from x = -1 to 1 as a quadratic Bézier: it bends
    /// left, most at its vertex, where its curvature is 2.
    #[test]
    fn a_parabola_bends_as_its_curvature() {
        let poles = [[-1.0, 1.0], [0.0, -1.0], [1.0, 1.0]];
        let knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        assert_eq!(problem(2, &poles, &knots), None);
        let spline = Spline::new(2, &poles, &knots);
        let [point, d1, _] = spline.derivatives(0.5);
        assert!((point[0]).abs() < 1e-15 && point[1].abs() < 1e-15);
        assert!(d1[1].abs() < 1e-15);
        let [left, right] = spline.sharpest_bends();
        assert!((left - 2.0).abs() < 1e-9, "{left}");
        assert_eq!(right, 0.0);
        assert!(!spline.has_corner(1e-9));
        assert_eq!(end_tangent(&poles, false), [1.0 / 5f64.sqrt(), -2.0 / 5f64.sqrt()]);
    }

    #[test]
    fn a_doubled_knot_can_turn_a_corner() {
        // Two lines as a degree 1 spline, turning a right angle at (1, 0).
        let poles = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
        let knots = [0.0, 0.0, 0.5, 1.0, 1.0];
        assert_eq!(problem(1, &poles, &knots), None);
        assert!(Spline::new(1, &poles, &knots).has_corner(1e-9));
        // Straight on, it doesn't.
        let poles = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        assert!(!Spline::new(1, &poles, &knots).has_corner(1e-9));
    }

    #[test]
    fn problems() {
        let poles = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0], [3.0, 1.0]];
        let cubic = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        assert_eq!(problem(3, &poles, &cubic), None);
        assert!(problem(0, &poles, &cubic).is_some());
        assert!(problem(4, &poles, &cubic).is_some());
        assert!(problem(3, &poles, &cubic[1..]).is_some());
        assert!(problem(3, &poles, &[0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]).is_some());
        assert!(problem(3, &poles, &[1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]).is_some());
        assert!(problem(3, &[[1.0, 2.0]; 4], &cubic).is_some());
        assert!(problem(3, &[[f64::NAN, 0.0]; 4], &cubic).is_some());
        let poles = [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0], [3.0, 1.0], [4.0, 0.0]];
        assert!(problem(1, &poles, &[0.0, 0.0, 0.2, 0.2, 0.6, 1.0, 1.0]).is_some());
    }
}
