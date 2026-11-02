//! Durations given on the command line (`--opera-lg-timeout`).

use std::time::Duration;

/// Parses `300` (seconds), `300s`, `5m` or `1h`; `0` means no timeout
/// (`None`).
pub fn parse_timeout(text: &str) -> Result<Option<Duration>, String> {
    let text = text.trim();
    let (number, unit) = match text.find(|c: char| !c.is_ascii_digit()) {
        Some(pos) => text.split_at(pos),
        None => (text, "s"),
    };
    let value: u64 = number
        .parse()
        .map_err(|_| format!("invalid duration {text:?} (expected e.g. 300, 5m or 1h)"))?;
    let seconds = match unit {
        "s" => Some(value),
        "m" => value.checked_mul(60),
        "h" => value.checked_mul(3600),
        _ => return Err(format!("invalid duration unit in {text:?} (use s, m or h)")),
    }
    .ok_or_else(|| format!("duration {text:?} is too large"))?;
    Ok((seconds > 0).then(|| Duration::from_secs(seconds)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats() {
        assert_eq!(parse_timeout("300"), Ok(Some(Duration::from_secs(300))));
        assert_eq!(parse_timeout("300s"), Ok(Some(Duration::from_secs(300))));
        assert_eq!(parse_timeout("5m"), Ok(Some(Duration::from_secs(300))));
        assert_eq!(parse_timeout("1h"), Ok(Some(Duration::from_secs(3600))));
        assert_eq!(parse_timeout("0"), Ok(None));
        assert_eq!(parse_timeout("0m"), Ok(None));
    }

    #[test]
    fn errors() {
        assert!(parse_timeout("").is_err());
        assert!(parse_timeout("5x").is_err());
        assert!(parse_timeout("m").is_err());
        assert!(parse_timeout("-1").is_err());
        assert!(parse_timeout("1.5h").is_err());
    }
}
