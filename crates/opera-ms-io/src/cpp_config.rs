//! Configuration files for the C++ programs (sigma, OPERA-LG), written
//! in the formats of OPERA-MS 0.9.
//!
//! The functions here only build text; creating directories, links and
//! writing the files is left to the caller. Inputs are processed line by
//! line keeping each line's own ending, and a line's key is what comes
//! before its first `=` (as Perl's `split(/\=/, $_)`), so blank lines and
//! a last line without a newline come out unchanged.

/// Number of edge sets (`pairedEdges_i0` … `pairedEdges_i5`).
pub const EDGE_SETS: usize = 6;

/// A configuration error (a required value is missing from the input).
#[derive(Debug, thiserror::Error)]
#[error("{0}")]
pub struct ConfigError(String);

/// The key of a config line: the text before the first `=`, or the
/// whole line (with its newline) if it has none.
fn key(line: &str) -> &str {
    line.split('=').next().unwrap_or(line)
}

/// `grep -A2 <pattern> <text> | cut -d "=" -f 2`, line by line, each line
/// with its newline (as Perl's backticks in list context return them).
fn grep_after2_cut2(text: &str, pattern: &str) -> Vec<String> {
    let lines: Vec<&str> = text.split_inclusive('\n').collect();
    let mut output = Vec::new();
    let mut last_printed: Option<usize> = None;
    let mut index = 0;
    while index < lines.len() {
        if !lines[index].contains(pattern) {
            index += 1;
            continue;
        }
        // grep separates non-adjacent groups of context with "--"
        if last_printed.is_some_and(|last| index > last + 1) {
            output.push("--\n".to_string());
        }
        let mut end = (index + 2).min(lines.len() - 1);
        let mut line = index;
        while line <= end {
            if line > index && lines[line].contains(pattern) {
                end = (line + 2).min(lines.len() - 1);
            }
            output.push(cut_field2(lines[line]));
            last_printed = Some(line);
            line += 1;
        }
        index = end + 1;
    }
    output
}

/// `cut -d "=" -f 2` of one line (keeping its newline); a line without
/// `=` is printed whole.
fn cut_field2(line: &str) -> String {
    let (body, newline) = match line.strip_suffix('\n') {
        Some(body) => (body, "\n"),
        None => (line, ""),
    };
    match body.split('=').nth(1) {
        Some(field) => format!("{field}{newline}"),
        None => line.to_string(),
    }
}

/// Mean and standard deviation of edge set `i` from the long-read mapping
/// config (`$dir_mapping/config`), each with its trailing newline.
fn mean_std(mapping_config: &str, edge_set: usize) -> Result<(String, String), ConfigError> {
    let values = grep_after2_cut2(mapping_config, &format!("pairedEdges_i{edge_set}"));
    match (values.get(1), values.get(2)) {
        (Some(mean), Some(std)) => Ok((mean.clone(), std.clone())),
        _ => Err(ConfigError(format!(
            "no lib_mean / lib_std for pairedEdges_i{edge_set} in the long-read mapping config"
        ))),
    }
}

/// `generate_sigma_config_file`: the sigma config of hierarchical
/// clustering, from the coverage estimation's sigma config.
///
/// The caller also creates `<dir_sigma>/pairedEdges_i<i>/` with a link
/// `pairedEdges_i<i>` to `<dir_mapping>/pairedEdges_i<i>` and the file
/// `lib.txt` ([`sigma_lib_txt`]).
pub fn hierarchical_sigma_config(coverage_sigma_config: &str, dir_sigma: &str) -> String {
    let mut out = String::new();
    for line in coverage_sigma_config.split_inclusive('\n') {
        match key(line) {
            "mapping_files" => {}
            "output_dir" => out.push_str(&format!("output_dir={dir_sigma}\n")),
            "edges_files" => {
                let edges: Vec<String> = (0..EDGE_SETS)
                    .map(|i| format!("{dir_sigma}/pairedEdges_i{i}/pairedEdges_i{i}"))
                    .collect();
                out.push_str(&format!("edges_files={}\n", edges.join(",")));
            }
            _ => out.push_str(line),
        }
    }
    out
}

/// The `lib.txt` written next to each edge set of hierarchical clustering.
pub fn sigma_lib_txt(mapping_config: &str, edge_set: usize) -> Result<String, ConfigError> {
    let (mean, std) = mean_std(mapping_config, edge_set)?;
    Ok(format!(
        "For the mapping file: \n\
         Mean length of the library is: {mean}\
         Standard deviation of the library is: {std}\
         Orientation of paired-reads is: forward ->...->  1st read...2nd read\n\
         Number of used paired reads for estimating the mean and the standard deviation of the library is: NA\n\
         Percentage of outliers is: NA%"
    ))
}

/// `generate_opera_config_file`: the OPERA-LG config of the long-read
/// scaffolding, from the bundling config (`lib_1_bundles/lib_1.ebconfig`).
pub fn long_read_opera_config(
    bundle_config: &str,
    mapping_config: &str,
    dir_ref_clustering: &str,
    dir_opera_lr: &str,
    contig_file: &str,
) -> Result<String, ConfigError> {
    let mut out = String::new();
    for line in bundle_config.split_inclusive('\n') {
        if line == "[LIB]\n" {
            for i in 0..EDGE_SETS {
                // mean and std keep their newline: a blank line follows
                // each of them
                let (mean, std) = mean_std(mapping_config, i)?;
                out.push_str(&format!(
                    "[LIB]\nmap_type=opera\n\
                     map_file={dir_ref_clustering}/NO_REPEAT/filtered_pairedEdges_i{i}\n\
                     lib_mean={mean}\nlib_std={std}\n"
                ));
            }
            continue;
        }
        match key(line) {
            "mapping_files" | "map_file" | "filter_repeat" => {}
            "output_folder" => out.push_str(&format!("output_folder={dir_opera_lr}\n")),
            "contig_file" => out.push_str(&format!("contig_file={contig_file}\n")),
            "cluster_threshold" => out.push_str("cluster_threshold=1\n"),
            _ => out.push_str(line),
        }
    }
    out.push_str("keep_repeat=yes\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const RUN: &str = "/data/ecoli_run/intermediate_files";

    fn fixture(name: &str) -> String {
        let path = format!(
            "{}/tests/fixtures/cpp_config/{name}",
            env!("CARGO_MANIFEST_DIR")
        );
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
    }

    #[test]
    fn hierarchical_sigma_config_matches_the_reference_run() {
        let generated = hierarchical_sigma_config(
            &fixture("coverage_estimation.sigma.config"),
            &format!("{RUN}/hierarchical_clustering"),
        );
        assert_eq!(generated, fixture("hierarchical_clustering.sigma.config"));
    }

    #[test]
    fn lib_txt_matches_the_reference_run() {
        let mapping = fixture("read_mapping.config");
        for i in 0..EDGE_SETS {
            assert_eq!(
                sigma_lib_txt(&mapping, i).unwrap(),
                fixture(&format!("lib_i{i}.txt")),
                "i{i}"
            );
        }
    }

    #[test]
    fn long_read_opera_config_matches_the_reference_run() {
        let generated = long_read_opera_config(
            &fixture("lib_1.ebconfig"),
            &fixture("read_mapping.config"),
            &format!("{RUN}/reference_clustering"),
            &format!("{RUN}/opera_long_read"),
            &format!("{RUN}/reference_clustering/single_strain_species_with_excluded.fa"),
        )
        .unwrap();
        assert_eq!(generated, fixture("opera_long_read.opera.config"));
    }

    #[test]
    fn missing_library_values_are_an_error() {
        assert!(sigma_lib_txt("map_file=x/pairedEdges_i0\nlib_mean=3\n", 0).is_err());
        assert!(long_read_opera_config("[LIB]\n", "", "r", "o", "c").is_err());
    }

    #[test]
    fn grep_context_like_grep() {
        let text = "a=1\nX=2\nb=3\nc=4\nd=5\nX=6\n";
        // two non-adjacent groups: "--" between them
        assert_eq!(
            grep_after2_cut2(text, "X"),
            ["2\n", "3\n", "4\n", "--\n", "6\n"]
        );
        // a match inside the context extends it
        assert_eq!(
            grep_after2_cut2("X=1\nX=2\ny=3\nz=4\nw=5\n", "X"),
            ["1\n", "2\n", "3\n", "4\n"]
        );
        assert_eq!(cut_field2("no equals sign\n"), "no equals sign\n");
        assert_eq!(cut_field2("a=b=c"), "b");
    }
}
