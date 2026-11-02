//! The R computations of the strain clustering: the kernel density estimate
//! and its local maxima, the negative binomial quantile, and the way R
//! prints a numeric vector.
//!
//! Written from the definitions (R's documentation and the mathematics),
//! not translated from R's GPL code; checked against R 3.5.1
//! on random inputs (`tests/r_cases/`). The density follows `density()` with
//! its defaults: Gaussian kernel, bandwidth `bw.nrd0`, 512 points, `cut =
//! 3`, linear binning on a grid extended by 4 bandwidths, convolution with
//! the kernel sampled on R's kernel grid, linear interpolation to the
//! output grid. The convolution is computed directly instead of by FFT, so
//! where the density is below R's FFT rounding noise (~1e-17 of its
//! maximum, e.g. in a wide empty coverage range) R can report spurious
//! local maxima that this version doesn't.

use std::f64::consts::PI;

/// Points of the density grid (`n = 512`).
const DENSITY_POINTS: usize = 512;
/// Window of the local-maximum search (`span = 11` in `compute_mode`).
const SPAN: usize = 11;

/// Sum with Neumaier's compensation.
fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut compensation) = (0.0_f64, 0.0_f64);
    for v in values {
        let t = sum + v;
        if sum.abs() >= v.abs() {
            compensation += (sum - t) + v;
        } else {
            compensation += (v - t) + sum;
        }
        sum = t;
    }
    sum + compensation
}

/// Sample standard deviation (`sd`, denominator n - 1).
pub fn sd(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mean = compensated_sum(values.iter().copied()) / n;
    let mean = mean + compensated_sum(values.iter().map(|v| v - mean)) / n;
    let squares = compensated_sum(values.iter().map(|v| (v - mean) * (v - mean)));
    (squares / (n - 1.0)).sqrt()
}

/// `quantile(x, p, type = 7)`: `(1 - h) x[j] + h x[j+1]` on the sorted
/// values, with `j + h = 1 + (n - 1) p` (1-based).
pub fn quantile7(sorted: &[f64], p: f64) -> f64 {
    let index = (sorted.len() - 1) as f64 * p;
    let lo = index.floor();
    let hi = index.ceil();
    let (a, b) = (sorted[lo as usize], sorted[hi as usize]);
    let h = index - lo;
    if index > lo && b != a {
        (1.0 - h) * a + h * b
    } else {
        a
    }
}

/// `bw.nrd0`: 0.9 min(sd, IQR / 1.34) n^-0.2 (falling back on sd, |x1|,
/// then 1 when the minimum is 0).
pub fn bw_nrd0(values: &[f64]) -> f64 {
    let hi = sd(values);
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let iqr = quantile7(&sorted, 0.75) - quantile7(&sorted, 0.25);
    let mut lo = hi.min(iqr / 1.34);
    if lo == 0.0 {
        lo = if hi != 0.0 {
            hi
        } else if values[0] != 0.0 {
            values[0].abs()
        } else {
            1.0
        };
    }
    0.9 * lo * (values.len() as f64).powf(-0.2)
}

/// `seq.int(from, to, length.out = n)` as R 3.5.1 computes it
/// (`from + i * by`, checked on 2,000 random grids).
fn seq(from: f64, to: f64, n: usize) -> Vec<f64> {
    let by = (to - from) / (n - 1) as f64;
    (0..n).map(|i| from + i as f64 * by).collect()
}

/// Gaussian density with standard deviation `sd`.
fn gaussian(x: f64, sd: f64) -> f64 {
    let z = x / sd;
    (-0.5 * z * z).exp() / ((2.0 * PI).sqrt() * sd)
}

/// `density(values)` with R's defaults: the output grid `x` and the
/// density `y` (512 points each).
pub fn density(values: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = DENSITY_POINTS;
    let count = values.len() as f64;
    let bw = bw_nrd0(values);
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let from = min - 3.0 * bw;
    let to = max + 3.0 * bw;
    let lo = from - 4.0 * bw;
    let up = to + 4.0 * bw;

    // linear binning on n points from lo to up, weight 1/count each
    let mut bins = vec![0.0_f64; n];
    let delta = (up - lo) / (n - 1) as f64;
    let weight = 1.0 / count;
    for &v in values {
        let position = (v - lo) / delta;
        let i = position.floor();
        let fraction = position - i;
        if i >= 0.0 && i <= (n - 2) as f64 {
            let i = i as usize;
            bins[i] += (1.0 - fraction) * weight;
            bins[i + 1] += fraction * weight;
        } else if i == -1.0 {
            bins[0] += fraction * weight;
        } else if i == (n - 1) as f64 {
            bins[n - 1] += (1.0 - fraction) * weight;
        }
    }

    // the kernel on R's kernel grid: 2n points from 0 to 2 (up - lo)
    let kernel_step = 2.0 * (up - lo) / (2 * n - 1) as f64;
    let kernel: Vec<f64> = (0..n)
        .map(|j| gaussian(j as f64 * kernel_step, bw))
        .collect();
    let smoothed: Vec<f64> = (0..n)
        .map(|m| {
            let value = compensated_sum(
                bins.iter()
                    .enumerate()
                    .filter(|(_, b)| **b != 0.0)
                    .map(|(a, b)| b * kernel[a.abs_diff(m)]),
            );
            value.max(0.0)
        })
        .collect();

    let grid = seq(lo, up, n);
    let x = seq(from, to, n);
    let y = x
        .iter()
        .map(|&v| interpolate(&grid, &smoothed, v))
        .collect();
    (x, y)
}

/// Linear interpolation of (xs, ys) at v (xs increasing, v within).
fn interpolate(xs: &[f64], ys: &[f64], v: f64) -> f64 {
    let j = xs.partition_point(|&x| x <= v);
    if j == 0 {
        return ys[0];
    }
    if j >= xs.len() {
        return ys[xs.len() - 1];
    }
    let (x0, x1) = (xs[j - 1], xs[j]);
    if v == x0 {
        return ys[j - 1];
    }
    ys[j - 1] + (ys[j] - ys[j - 1]) * ((v - x0) / (x1 - x0))
}

/// Relative level below which the density counts as 0: R's FFT gives
/// rounding noise there (~1e-17 of the maximum), the direct convolution
/// tiny exact values; neither is meaningful.
const DENSITY_FLOOR: f64 = 1e-14;

/// The local maxima of `compute_mode`'s R code: with `embed(y, 11)` and
/// `which.max` on each row (reversed windows), a point is a mode when it
/// is strictly above the 5 following points and not below the 5 preceding
/// ones. Nightly fix: every point with 5 neighbours on each side is a
/// candidate, once (the R code's logical index was 10 shorter than `x` and
/// R recycled it: modes at points 6-10, 1-based, were repeated at points
/// 508-512, and points 503-507 were never candidates).
pub fn density_modes(values: &[f64]) -> Vec<f64> {
    let (x, mut y) = density(values);
    let floor = y.iter().copied().fold(0.0, f64::max) * DENSITY_FLOOR;
    for v in &mut y {
        if *v < floor {
            *v = 0.0;
        }
    }
    let half = SPAN / 2;
    (half..y.len().saturating_sub(half))
        .filter(|&c| {
            y[c + 1..=c + half].iter().all(|&v| y[c] > v)
                && y[c - half..c].iter().all(|&v| y[c] >= v)
        })
        .map(|c| x[c])
        .collect()
}

/// The decimal exponent and significant digits (1 to `digits`) of `v`
/// written with `digits` significant digits.
fn significant(v: f64, digits: usize) -> (i32, usize) {
    if v == 0.0 {
        return (0, 1);
    }
    let text = format!("{:.*e}", digits - 1, v.abs());
    let (mantissa, exponent) = text.split_once('e').unwrap_or((&text, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits_used = mantissa.replace('.', "").trim_end_matches('0').len().max(1);
    (exponent, digits_used)
}

/// The values of a numeric vector as `print()` shows them (7 significant
/// digits): one common format for the whole vector, fixed notation with
/// the largest number of decimals any element needs, or scientific
/// notation when that is narrower.
pub fn r_printed(values: &[f64]) -> Vec<String> {
    const DIGITS: usize = 7;
    let finite: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    let mut decimals = 0usize;
    let mut left = 1usize;
    let mut max_sig = 1usize;
    let (mut min_exp, mut max_exp) = (i32::MAX, i32::MIN);
    let negative = finite.iter().any(|v| *v < 0.0);
    for &v in &finite {
        let (exp, sig) = significant(v, DIGITS);
        decimals = decimals.max((sig as i32 - 1 - exp).max(0) as usize);
        left = left.max((exp + 1).max(1) as usize);
        max_sig = max_sig.max(sig);
        min_exp = min_exp.min(exp);
        max_exp = max_exp.max(exp);
    }
    let fixed_width = usize::from(negative) + left + decimals + usize::from(decimals > 0);
    let exp_width = if max_exp.abs() >= 100 || min_exp.abs() >= 100 {
        5
    } else {
        4
    };
    let sci_width = usize::from(negative) + max_sig + usize::from(max_sig > 1) + exp_width;
    values
        .iter()
        .map(|&v| {
            if v.is_nan() {
                "NaN".to_string()
            } else if v.is_infinite() {
                if v > 0.0 { "Inf" } else { "-Inf" }.to_string()
            } else if finite.is_empty() || fixed_width <= sci_width {
                format!("{v:.decimals$}")
            } else {
                let text = format!("{:.*e}", max_sig - 1, v);
                // R writes the exponent with a sign and at least 2 digits
                match text.split_once('e') {
                    Some((m, e)) => {
                        let e: i32 = e.parse().unwrap_or(0);
                        format!("{m}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
                    }
                    None => text,
                }
            }
        })
        .collect()
}

/// `qnbinom(p, size, mu)`: the smallest y with P(Y <= y) >= p for a
/// negative binomial of mean `mu` and dispersion `size` (NaN for invalid
/// parameters).
pub fn qnbinom(p: f64, size: f64, mu: f64) -> f64 {
    if !(p.is_finite() && size.is_finite() && mu.is_finite())
        || size <= 0.0
        || mu < 0.0
        || !(0.0..=1.0).contains(&p)
    {
        return f64::NAN;
    }
    if mu == 0.0 || p == 0.0 {
        return 0.0;
    }
    if p + 1.01 * f64::EPSILON >= 1.0 {
        return f64::INFINITY;
    }
    let target = p * (1.0 - 64.0 * f64::EPSILON);
    // P(Y = k) = C(k + size - 1, k) prob^size (1 - prob)^k, prob = size / (size + mu)
    let log_prob = (size / (size + mu)).ln();
    let log_q = (mu / (size + mu)).ln();
    let mut log_term = size * log_prob;
    let (mut sum, mut compensation) = (0.0_f64, 0.0_f64);
    let mut k = 0.0_f64;
    loop {
        let term = log_term.exp();
        let t = sum + term;
        if sum.abs() >= term.abs() {
            compensation += (sum - t) + term;
        } else {
            compensation += (term - t) + sum;
        }
        sum = t;
        if sum + compensation >= target {
            return k;
        }
        log_term += ((k + size) / (k + 1.0)).ln() + log_q;
        k += 1.0;
        if k > 1e9 {
            return f64::INFINITY;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn printing_like_r() {
        let values = [3.053728111, 26.27694912, 158.4944861];
        assert_eq!(r_printed(&values), ["3.053728", "26.276949", "158.494486"]);
        assert_eq!(r_printed(&[100000.0]), ["1e+05"]);
        assert_eq!(r_printed(&[12.0, 1.5]), ["12.0", "1.5"]);
        assert_eq!(r_printed(&[123456.0]), ["123456"]);
    }

    #[test]
    fn quantiles_and_bandwidth() {
        let sorted = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(quantile7(&sorted, 0.25), 1.75);
        assert_eq!(quantile7(&sorted, 0.75), 3.25);
        // sd = 1.2909944, IQR / 1.34 = 1.119403
        let bw = bw_nrd0(&[1.0, 2.0, 3.0, 4.0]);
        assert!((bw - 0.9 * 1.5 / 1.34 * 4f64.powf(-0.2)).abs() < 1e-15);
    }

    #[test]
    fn negative_binomial_quantiles() {
        // values printed by R 3.5.1 in the E. coli reference run
        assert_eq!(qnbinom(0.02, 79.657812, 160.50947322085), 118.0);
        assert_eq!(qnbinom(0.98, 79.657812, 160.50947322085), 208.0);
        assert_eq!(qnbinom(0.02, 846.450133355214, 26.3080294244547), 16.0);
        assert_eq!(qnbinom(0.98, 846.450133355214, 26.3080294244547), 38.0);
        assert_eq!(qnbinom(0.02, 79.657812, 3.09255094615569), 0.0);
        assert_eq!(qnbinom(0.98, 79.657812, 13396.1714716397), 16666.0);
        assert_eq!(qnbinom(0.98, 1.0, 0.0), 0.0);
    }

    /// Opt-in: `RSTATS_CASES=<file from tests/r_cases/rstats_cases.R>`.
    #[test]
    fn against_r_cases() {
        let Ok(path) = std::env::var("RSTATS_CASES") else {
            return;
        };
        let text = std::fs::read_to_string(path).unwrap();
        let (mut density_cases, mut q_cases, mut noise_only) = (0, 0, 0);
        let mut failures = Vec::new();
        for line in text.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            match fields[0] {
                "modes" => {
                    density_cases += 1;
                    let values: Vec<f64> =
                        fields[1].split(',').map(|v| v.parse().unwrap()).collect();
                    let modes = density_modes(&values);
                    let mine = if modes.is_empty() {
                        "numeric(0)".to_string()
                    } else {
                        r_printed(&modes).join(" ")
                    };
                    if mine == fields[2] {
                        continue;
                    }
                    // R's modes above the noise floor (its FFT noise makes
                    // spurious modes where the density is ~1e-17 of its max)
                    let levels: Vec<f64> = fields.get(3).map_or(vec![], |l| {
                        l.split(' ').filter_map(|v| v.parse().ok()).collect()
                    });
                    let significant: Vec<&str> = fields[2]
                        .split(' ')
                        .zip(levels.iter())
                        .filter(|(_, l)| **l >= DENSITY_FLOOR)
                        .map(|(m, _)| m)
                        .collect();
                    let mine_values: Vec<f64> =
                        mine.split(' ').filter_map(|v| v.parse().ok()).collect();
                    let r_values: Vec<f64> =
                        significant.iter().filter_map(|v| v.parse().ok()).collect();
                    if mine_values == r_values {
                        noise_only += 1;
                    } else {
                        failures.push(format!(
                            "modes n={}:\n  R    {}\n  R>fl {}\n  mine {mine}",
                            values.len(),
                            fields[2],
                            significant.join(" ")
                        ));
                    }
                }
                "qnbinom" => {
                    q_cases += 1;
                    let [p, size, mu]: [f64; 3] = [1, 2, 3].map(|i| fields[i].parse().unwrap());
                    let r: f64 = fields[4].parse().unwrap();
                    let mine = qnbinom(p, size, mu);
                    if mine != r {
                        failures.push(format!("qnbinom({p}, {size}, {mu}): R {r}, mine {mine}"));
                    }
                }
                _ => {}
            }
        }
        for f in failures.iter().take(20) {
            opera_ms_core::note!("{f}");
        }
        opera_ms_core::note!(
            "{density_cases} density cases ({noise_only} differ from R only by its noise modes), {q_cases} qnbinom cases, {} failures",
            failures.len()
        );
        assert!(failures.is_empty());
    }
}
