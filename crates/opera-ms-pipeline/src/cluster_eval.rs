//! Coverage clustering parameter: checks sigma's clusters against the contig coverage, and re-runs sigma
//! with other values of R if too much of the assembly is in outlier
//! contigs or a cluster is too big.
//!
//! In the diagnostic file `r_estimate_<R>.dat`, clusters with the same
//! length (or outlier length) are listed by cluster ID. The summary line (the last one, the only one the
//! pipeline reads) doesn't depend on that order.

use std::collections::HashMap;
use std::fs;

use anyhow::{Context, Result, ensure};
use opera_ms_core::bench::StepRecorder;
use opera_ms_core::parallel::run_parallel;
use opera_ms_core::runner::{Cmd, Output, Pipeline};

use crate::edge_filter::perl_number as perl_value;
use crate::files::perl_number as perl_format;

const NB_STEP: usize = 10;
const FRACTION_THRESHOLD: f64 = 0.01;
const CLUSTER_FRACTION_THRESHOLD: f64 = 0.01;
const MIN_LENGTH_THRESHOLD: f64 = 200_000.0;
const TOP_OUTLIER_CLUSTERS: usize = 5;

fn chop(text: &str) -> &str {
    let mut chars = text.chars();
    chars.next_back();
    chars.as_str()
}

#[derive(Debug, Default, Clone)]
struct ClusterInfo {
    length: f64,
    nb_element: f64,
    nb_outlier: f64,
    outlier_length: f64,
}

/// The cluster evaluation (written to `r_estimate_<R>.dat`).
pub fn cluster_evaluation(coverage_file: &str, clusters_file: &str) -> Result<String> {
    // contig -> (coverage, length)
    let mut contigs: HashMap<String, (f64, f64)> = HashMap::new();
    let mut lines = coverage_file.split_inclusive('\n');
    let header = chop(lines.next().context("empty coverage file")?);
    let window_size = perl_value(header.split(' ').nth(3).unwrap_or(""));
    while let Some(line) = lines.next() {
        let line = line.strip_suffix('\n').unwrap_or(line);
        let fields: Vec<&str> = line.split('\t').collect();
        let contig = fields.first().copied().unwrap_or("");
        let length = perl_value(fields.get(1).copied().unwrap_or(""));
        let nb_window = perl_value(fields.get(4).copied().unwrap_or(""));
        let read_count = perl_value(chop(lines.next().unwrap_or("")));
        lines.next();
        let denominator = nb_window * window_size;
        ensure!(
            denominator != 0.0,
            "contig {contig}: no window (division by zero)"
        );
        contigs.insert(contig.to_string(), (read_count / denominator, length));
    }

    let mut clusters: HashMap<String, ClusterInfo> = HashMap::new();
    let (mut outlier_length_sum, mut length_sum) = (0.0, 0.0);
    for line in clusters_file.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let contig = fields.first().copied().unwrap_or("");
        let cluster = fields.get(1).copied().unwrap_or("");
        let cluster_mean = perl_value(fields.get(3).copied().unwrap_or(""));
        // a contig missing from the coverage file has coverage and length 0
        let (cov, length) = contigs.get(contig).copied().unwrap_or((0.0, 0.0));
        let info = clusters.entry(cluster.to_string()).or_default();
        info.nb_element += 1.0;
        if cluster_mean / 2.0 > cov || cov > cluster_mean * 2.0 {
            info.nb_outlier += 1.0;
            info.outlier_length += length;
            outlier_length_sum += length;
        }
        info.length += length;
        length_sum += length;
    }

    // Ties: by cluster ID (numeric when possible), see the module docs.
    let mut ids: Vec<&String> = clusters.keys().collect();
    ids.sort_by(|a, b| match (a.parse::<u64>(), b.parse::<u64>()) {
        (Ok(x), Ok(y)) => x.cmp(&y),
        _ => a.cmp(b),
    });
    let line = |id: &str, c: &ClusterInfo| -> Result<String> {
        ensure!(
            c.length != 0.0 && c.nb_element != 0.0,
            "cluster {id}: length 0 (division by zero)"
        );
        Ok(format!(
            "{id}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            perl_format(c.length),
            perl_format(c.outlier_length),
            perl_format(c.outlier_length / c.length),
            perl_format(c.nb_outlier),
            perl_format(c.nb_element),
            perl_format(c.nb_outlier / c.nb_element)
        ))
    };

    let mut out = String::from("\n\n########## SORTED BY LENGTH\n");
    let mut by_length = ids.clone();
    by_length.sort_by(|a, b| clusters[*b].length.total_cmp(&clusters[*a].length));
    let mut max_length = -1.0;
    let (mut max_fraction, mut max_fraction_length) = (0.0, 0.0);
    for id in by_length {
        let c = &clusters[id];
        let fraction = c.outlier_length / c.length;
        if max_length == -1.0 {
            max_length = c.length;
        }
        out.push_str(&line(id, c)?);
        if fraction > max_fraction {
            max_fraction = fraction;
            max_fraction_length = c.length;
        }
        if c.length < MIN_LENGTH_THRESHOLD {
            break;
        }
    }
    out.push_str("\n\n########## SORTED BY OUTLIER_LENGTH\n");
    let mut by_outlier = ids;
    by_outlier.sort_by(|a, b| {
        clusters[*b]
            .outlier_length
            .total_cmp(&clusters[*a].outlier_length)
    });
    for id in by_outlier.into_iter().take(TOP_OUTLIER_CLUSTERS) {
        out.push_str(&line(id, &clusters[id])?);
    }
    out.push_str(&format!(
        "\n{}\t{}\t{}\t{}\t{}\n",
        perl_format(length_sum),
        perl_format(outlier_length_sum),
        perl_format(max_length),
        perl_format(max_fraction),
        perl_format(max_fraction_length)
    ));
    Ok(out)
}

/// Evaluates the clusters of one R: writes `out_file`, adds a
/// line to the summary and says whether the clusters are good.
fn evaluate(
    coverage: &str,
    clusters_file: &str,
    out_file: &str,
    r_value: &str,
    cluster_length_threshold: &str,
    summary: &mut String,
) -> Result<bool> {
    let clusters = fs::read_to_string(clusters_file)
        .with_context(|| format!("cannot read {clusters_file}"))?;
    let text = cluster_evaluation(coverage, &clusters)
        .with_context(|| format!("Error in during cluster evaluation of {clusters_file}"))?;
    fs::write(out_file, &text).with_context(|| format!("cannot write {out_file}"))?;
    let last = text
        .trim_end_matches('\n')
        .rsplit('\n')
        .next()
        .unwrap_or("");
    let m: Vec<&str> = last.split('\t').collect();
    let field = |i: usize| m.get(i).copied().unwrap_or("");
    let total = perl_value(field(0));
    ensure!(
        total != 0.0,
        "{out_file}: total length 0 (division by zero)"
    );
    let outlier_fraction = perl_value(field(1)) / total;
    let max_cluster_length = field(2);
    let max_cluster_outlier_fraction = field(3);
    summary.push_str(&format!(
        "{r_value}\t{}\t{}\t{max_cluster_length}\t{cluster_length_threshold}\t\
         {max_cluster_outlier_fraction}\t{}\t{}\n",
        perl_format(outlier_fraction),
        perl_format(FRACTION_THRESHOLD),
        field(4),
        perl_format(CLUSTER_FRACTION_THRESHOLD)
    ));
    let threshold = perl_value(cluster_length_threshold);
    Ok(perl_value(r_value) > 1.0
        && outlier_fraction < FRACTION_THRESHOLD
        && (threshold == -1.0 || perl_value(max_cluster_length) < threshold)
        && perl_value(max_cluster_outlier_fraction) < CLUSTER_FRACTION_THRESHOLD)
}

/// The R estimate: `sigma` runs through `sigma_worker` (a command
/// to which the config file is appended).
pub fn refine_r_estimate(
    rec: &StepRecorder<'_>,
    coverage_file: &str,
    sigma_dir: &str,
    num_processor: usize,
    cluster_length_threshold: &str,
    sigma_cmd_display: &str,
    sigma_worker: &Cmd,
) -> Result<()> {
    let coverage = fs::read_to_string(coverage_file)
        .with_context(|| format!("cannot read {coverage_file}"))?;
    let mut tmp = fs::read_to_string(format!("{sigma_dir}/r_estimate_value.dat"))
        .with_context(|| format!("cannot read {sigma_dir}/r_estimate_value.dat"))?;
    tmp.pop();
    let values: Vec<&str> = tmp.split('\t').collect();
    ensure!(
        values.len() >= 2,
        "{sigma_dir}/r_estimate_value.dat: expected two values"
    );
    let (r0, r1) = (perl_value(values[0]), perl_value(values[1]));
    let step = (r1 - r0) / NB_STEP as f64;
    opera_ms_core::note!(
        " *** R refinement \t{}\t{}\t{}",
        values[0],
        values[1],
        perl_format(step)
    );
    let summary_file = format!("{sigma_dir}/r_value_evaluation_summary.dat");
    let mut summary = String::new();
    let r_first = values[0].to_string();
    let good = evaluate(
        &coverage,
        &format!("{sigma_dir}/clusters"),
        &format!("{sigma_dir}/r_estimate_{r_first}.dat"),
        &r_first,
        cluster_length_threshold,
        &mut summary,
    )?;
    if !good {
        // sigma with other R values
        let config = fs::read_to_string(format!("{sigma_dir}/sigma.config"))
            .with_context(|| format!("cannot read {sigma_dir}/sigma.config"))?;
        let cmd_file = format!("{sigma_dir}/cmd.sh");
        let mut cmd_text = String::new();
        let mut jobs = Vec::new();
        for i in 1..=NB_STEP {
            let r_value = perl_format(r0 + i as f64 * step);
            let new_dir = format!("{sigma_dir}/sigma_{r_value}");
            let new_config = format!("{new_dir}/sigma.config");
            fs::create_dir_all(&new_dir).with_context(|| format!("cannot create {new_dir}"))?;
            let mut text: String = config
                .split_inclusive('\n')
                .filter(|l| !l.contains("R_VALUE") && !l.contains("output_dir"))
                .collect();
            text.push_str(&format!("output_dir={new_dir}\nR_VALUE={r_value}\n"));
            fs::write(&new_config, text).with_context(|| format!("cannot write {new_config}"))?;
            let err = format!("{cmd_file}.{i}.err");
            cmd_text.push_str(&format!("{sigma_cmd_display} {new_config} 2> {err}\n"));
            let pipeline = Pipeline::new(sigma_worker.clone().arg(&new_config))
                .stdout(Output::Stderr)
                .stderr(Output::File(err.into()));
            jobs.push(move || {
                rec.substep(&format!("sigma R={r_value}"), |sub| {
                    let report = pipeline.run()?;
                    sub.record(&report);
                    Ok::<(), anyhow::Error>(())
                })
            });
        }
        fs::write(&cmd_file, cmd_text).with_context(|| format!("cannot write {cmd_file}"))?;
        run_parallel(jobs, num_processor.max(1))
            .map_err(|e| anyhow::anyhow!("{e}"))
            .with_context(|| {
                format!("Error in during sigma. Please see {cmd_file}.*.err for details")
            })?;
        let mut chosen = String::new();
        for i in 1..=NB_STEP {
            chosen = perl_format(r0 + i as f64 * step);
            let clusters = format!("{sigma_dir}/sigma_{chosen}/clusters");
            let out = format!("{sigma_dir}/r_estimate_{chosen}.dat");
            if evaluate(
                &coverage,
                &clusters,
                &out,
                &chosen,
                cluster_length_threshold,
                &mut summary,
            )? {
                break;
            }
        }
        // the clusters used for the next steps
        let link = format!("{sigma_dir}/clusters");
        fs::remove_file(&link).with_context(|| format!("cannot remove {link}"))?;
        std::os::unix::fs::symlink(format!("{sigma_dir}/sigma_{chosen}/clusters"), &link)
            .with_context(|| format!("cannot link {link}"))?;
    }
    fs::write(&summary_file, summary).with_context(|| format!("cannot write {summary_file}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const COVERAGE: &str = "1 500 80 340\n\
        c1\t300000\tx\tx\t10\n3400\n1\n\
        c2\t300000\tx\tx\t10\n3400\n1\n\
        c3\t1000\tx\tx\t10\n34000\n1\n\
        c4\t2000\tx\tx\t10\n3400\n1\n";

    #[test]
    fn evaluation_lines_and_summary() {
        // cov = 1 for c1, c2, c4; 10 for c3 (an outlier of cluster 1)
        let clusters = "c1\t1\t0\t1.0\nc3\t1\t0\t1.0\nc2\t2\t0\t1.0\nc4\t3\t0\t1.0\n";
        let text = cluster_evaluation(COVERAGE, clusters).unwrap();
        assert_eq!(
            text,
            "\n\n########## SORTED BY LENGTH\n\
             1\t301000\t1000\t0.00332225913621262\t1\t2\t0.5\n\
             2\t300000\t0\t0\t0\t1\t0\n\
             3\t2000\t0\t0\t0\t1\t0\n\
             \n\n########## SORTED BY OUTLIER_LENGTH\n\
             1\t301000\t1000\t0.00332225913621262\t1\t2\t0.5\n\
             2\t300000\t0\t0\t0\t1\t0\n\
             3\t2000\t0\t0\t0\t1\t0\n\
             \n603000\t1000\t301000\t0.00332225913621262\t301000\n"
        );
    }

    #[test]
    fn summary_matches_the_reference_format() {
        let clusters = "c1\t1\t0\t1.0\nc2\t1\t0\t1.0\n";
        let mut summary = String::new();
        let dir = tempfile::tempdir().unwrap();
        let cl = dir.path().join("clusters");
        std::fs::write(&cl, clusters).unwrap();
        let out = dir.path().join("r.dat");
        let good = evaluate(
            COVERAGE,
            cl.to_str().unwrap(),
            out.to_str().unwrap(),
            "79.657812",
            "6000000",
            &mut summary,
        )
        .unwrap();
        assert!(good);
        assert_eq!(summary, "79.657812\t0\t0.01\t600000\t6000000\t0\t0\t0.01\n");
    }
}
