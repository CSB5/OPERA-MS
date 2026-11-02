//! Statistics computed as Perl's `Statistics::Basic` computes them,
//! operation for operation, so that results match OPERA-MS 0.9 exactly.
//!
//! `Statistics::Basic` 1.6611 (`$UNBIAS` off), checked on the installed
//! module: `mean` = left-to-right sum / n; `variance` = sum of
//! `($x - $mean) ** 2` / n (population); `stddev` = `sqrt(variance)`;
//! `median` = middle of the numerically sorted values, or the mean of the
//! two middle ones. Empty input gives `n/a`, which is 0 as a number.
//! (Its objects print rounded to 2 decimals but compute with the full
//! value; these functions return the full value.)

/// How a `Statistics::Basic` object prints (`"$median"`): Number::Format
/// 1.75 `format_number("$v", 2)` with its defaults (thousands separator
/// `,`, decimal point `.`, no trailing zeros): the value is printed with
/// Perl's `%.15g`, rounded as `int(|v| * 100 + 0.5 + 1e-14) / 100`, and
/// its integer part grouped by 3. (Perl's output depends on the
/// locale through `localeconv()`; this is the C / English one.)
pub fn basic_string(value: f64) -> String {
    let number: f64 = crate::files::perl_number(value).parse().unwrap_or(value);
    let negative = number < 0.0;
    let product = number.abs() * 100.0;
    let rounded = (product + 0.5 + 1e-14).trunc() / 100.0;
    let text = crate::files::perl_number(rounded);
    if text.contains('e') {
        return text;
    }
    let integer = crate::files::perl_number(rounded.trunc());
    let decimal = if integer.len() < text.len() {
        &text[integer.len() + 1..]
    } else {
        ""
    };
    let mut grouped = String::new();
    for (i, c) in integer.chars().enumerate() {
        if i > 0 && (integer.len() - i).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(c);
    }
    let result = if decimal.is_empty() {
        grouped
    } else {
        format!("{grouped}.{decimal}")
    };
    if negative {
        format!("-{result}")
    } else {
        result
    }
}

/// `Statistics::Basic::mean`.
pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().fold(0.0, |sum, v| sum + v) / values.len() as f64
}

/// `Statistics::Basic::variance` (population).
pub fn variance(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let m = mean(values);
    // Perl's `$d ** 2` calls C pow(d, 2), which glibc rounds as d * d
    values.iter().fold(0.0, |sum, v| sum + (v - m) * (v - m)) / values.len() as f64
}

/// `Statistics::Basic::stddev` (population).
pub fn stddev(values: &[f64]) -> f64 {
    variance(values).sqrt()
}

/// `Statistics::Basic::median`.
pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let center = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        sorted[center]
    } else {
        (sorted[center] + sorted[center - 1]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Values printed by the installed Statistics::Basic (perl -e ...).
    #[test]
    fn same_as_statistics_basic() {
        let a = [1.0, 2.0, 3.0, 4.0, 10.123456789];
        assert_eq!(median(&a), 3.0);
        assert_eq!(stddev(&a), 3.2091642130311797);
        assert_eq!(median(&a) - 1.5 * stddev(&a), -1.8137463195467696);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        let c = [0.123456789, 0.987654321, 0.5];
        assert_eq!((median(&c), stddev(&c)), (0.5, 0.3537778494502906));
        assert_eq!((median(&[]), stddev(&[])), (0.0, 0.0));
        assert_eq!(stddev(&[5.0]), 0.0);
    }

    // `perl -e 'print median($v)'` with the installed modules.
    #[test]
    fn prints_as_statistics_basic() {
        for (value, printed) in [
            (0.0234, "0.02"),
            (1234.5678, "1,234.57"),
            (25.0, "25"),
            (0.005, "0.01"),
            (0.015, "0.02"),
            (0.025, "0.03"),
            (2.675, "2.68"),
            (1e-5, "0"),
            (123456.789, "123,456.79"),
            (-12.3456, "-12.35"),
            (0.0, "0"),
            (1.005, "1.01"),
            (0.125, "0.13"),
            (0.135, "0.14"),
            (99.995, "100"),
        ] {
            assert_eq!(basic_string(value), printed, "{value}");
        }
    }
}
