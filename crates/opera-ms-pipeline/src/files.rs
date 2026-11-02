//! The small files the pipeline writes itself, in the formats of OPERA-MS
//! 0.9.
//!
//! Pure functions from input text to output text. Some Perl text handling
//! is kept on purpose, see the comments (e.g. `split(/\t/, $_)` keeps the
//! newline in the last field).

use std::collections::HashMap;

/// The settings part of `cov_estimate.config`. `values` are the run's
/// settings by key name; only the keys of the list below are written, in
/// its order, and only if set.
pub fn utils_config(values: &HashMap<&str, String>) -> String {
    const ORDER: [&str; 12] = [
        "ILLUMINA_READ_1",
        "ILLUMINA_READ_2",
        "LONG_READ",
        "OUTPUT_DIR",
        "CONTIG_EDGE_LEN",
        "CONTIG_WINDOW_LEN",
        "CONTIG_LEN_THR",
        "KMER_SIZE",
        "CONTIGS_FILE",
        "LONG_READ_MAPPER",
        "SHORT_READ_ASSEMBLER",
        "NUM_PROCESSOR",
    ];
    ORDER
        .iter()
        .filter_map(|key| values.get(key).map(|value| format!("{key} {value}\n")))
        .collect()
}

/// The lines `assembly_graph_creation` appends to the copy of the utils
/// config to make `cov_estimate.config`.
pub fn cov_estimate_config_additions(
    assembled_contig_file: Option<&str>,
    illumina_read_1: &str,
    mapping_dir: &str,
    samtools_dir: &str,
    tools_dir: &str,
) -> String {
    let mut out = String::new();
    if let Some(contigs) = assembled_contig_file {
        out.push_str(&format!("\nCONTIGS_FILE {contigs}\n"));
    }
    if !illumina_read_1.contains(',') {
        out.push_str(&format!("\nLIB {mapping_dir}/opera.bam\n"));
        out.push_str(&format!("\nMAPPING_FILES {mapping_dir}/opera.bam\n"));
    } else {
        // Legacy bug kept (multi-sample isn't supported): the BAM list is
        // built from an option that is never set, so it is empty.
        out.push_str(&format!("\nLIB {mapping_dir}/opera_1.bam\n"));
        out.push_str("\nMAPPING_FILES \n");
    }
    out.push_str(&format!("SAMTOOLS_DIR {samtools_dir}\n"));
    out.push_str(&format!("TOOLS_DIR {tools_dir}\n"));
    out
}

/// Lines with their newline, as Perl's `while(<FILE>)` reads them.
fn perl_lines(text: &str) -> impl Iterator<Item = &str> {
    text.split_inclusive('\n')
}

fn chomp(line: &str) -> &str {
    line.strip_suffix('\n').unwrap_or(line)
}

/// `multi_to_linear_fa`: one line per sequence. A header is any line
/// containing `>`; the output has no newline after the last sequence.
pub fn multi_to_linear_fa(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    let mut first = true;
    for line in perl_lines(input) {
        let line = chomp(line);
        if line.contains('>') {
            if !first {
                out.push('\n');
            }
            out.push_str(line);
            out.push('\n');
            first = false;
        } else {
            out.push_str(line);
        }
    }
    // end with a newline
    if !out.is_empty() && !out.ends_with('\n') {
        out.push('\n');
    }
    out
}

/// `split(/\t/, $_)` (trailing empty fields dropped, as Perl does).
fn split_tabs(line: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = line.split('\t').collect();
    while fields.last().is_some_and(|f| f.is_empty()) {
        fields.pop();
    }
    fields
}

/// `split(/\s+/, $_)`: a leading empty field if the line starts with
/// whitespace, trailing empty fields dropped.
fn split_whitespace_perl(line: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut rest = line;
    loop {
        match rest.find(|c: char| c.is_whitespace()) {
            Some(pos) => {
                fields.push(&rest[..pos]);
                let after = &rest[pos..];
                rest = after.trim_start_matches(|c: char| c.is_whitespace());
                if rest.is_empty() {
                    break;
                }
            }
            None => {
                fields.push(rest);
                break;
            }
        }
    }
    while fields.last().is_some_and(|f| f.is_empty()) {
        fields.pop();
    }
    fields
}

/// `generate_assembly_stats`: the stats text (printed on stdout and
/// written to `assembly.stats`, each followed by one more newline) from
/// `contig_info.txt`. `time` is Perl's `localtime` string.
pub fn assembly_stats(contig_info: &str, time: &str) -> String {
    let mut sizes: Vec<u64> = Vec::new();
    let (mut c1mb, mut c500kb, mut c100kb, mut total) = (0, 0, 0, 0u64);
    for line in perl_lines(contig_info).skip(1) {
        let fields = split_tabs(line);
        let size: u64 = fields.get(1).map_or(0, |s| perl_numeric(s));
        c1mb += u64::from(size >= 1_000_000);
        c500kb += u64::from(size >= 500_000);
        c100kb += u64::from(size >= 100_000);
        total += size;
        sizes.push(size);
    }
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    // compute_Nx(50, ...): first index where the running total reaches
    // half the assembly (Perl's $tab[-1] on an empty list is undef)
    let half = total as f64 * 0.5;
    let mut running = 0u64;
    let mut index = 0usize;
    while (running as f64) < half && index < sizes.len() {
        running += sizes[index];
        index += 1;
    }
    let n50 = if index == 0 {
        sizes.last().map(u64::to_string).unwrap_or_default()
    } else {
        sizes[index - 1].to_string()
    };
    let max = sizes.first().map(u64::to_string).unwrap_or_default();
    format!(
        "\n[{time}]\tAssembly stats\nNumber of contigs: {} \nAssembly size: {total} bp\n\
         Max contig size: {max} bp\nContig(s) longer than 1Mbp: {c1mb} \n\
         Contig(s) longer than 500kbp: {c500kb}\nContig(s) longer than 100kbp: {c100kb}\n\
         Contig N50: {n50} bp\n",
        sizes.len()
    )
}

/// Perl's numeric value of a string: leading digits (after optional
/// whitespace), 0 if none.
fn perl_numeric(text: &str) -> u64 {
    let text = text.trim_start();
    let digits: &str = &text[..text
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(text.len())];
    digits.parse().unwrap_or(0)
}

/// `write_final_assembly`: the final `contigs.fasta` from `contig_info.txt`,
/// `scaffolds.scaf.cname` and the gap-filled scaffolds (2-line FASTA).
pub fn final_contigs(contig_info: &str, cname: &str, filled: &str) -> String {
    let mut info: HashMap<String, String> = HashMap::new();
    for line in perl_lines(contig_info) {
        let fields = split_tabs(line);
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        info.insert(
            get(0).to_string(),
            format!(
                "length: {} cov_short: {} cov_long: {} cluster: {} ",
                get(1),
                get(2),
                get(3),
                get(4)
            ),
        );
    }
    let mut scaff_to_contig: HashMap<String, String> = HashMap::new();
    for line in perl_lines(cname) {
        // `chop`: remove the last character, whatever it is
        let mut chars = line.chars();
        chars.next_back();
        let fields = split_tabs(chars.as_str());
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        scaff_to_contig.insert(get(0).to_string(), get(1).to_string());
    }
    let mut out = String::with_capacity(filled.len() + filled.len() / 50);
    let mut lines = perl_lines(filled);
    while let Some(header) = lines.next() {
        let first = split_whitespace_perl(header).first().copied().unwrap_or("");
        let scaff = first.get(1..).unwrap_or("");
        let contig = scaff_to_contig.get(scaff).map(String::as_str).unwrap_or("");
        let contig_info = info.get(contig).map(String::as_str).unwrap_or("");
        out.push_str(&format!(">{contig} {contig_info}\n"));
        if let Some(seq) = lines.next() {
            out.push_str(seq);
        }
    }
    out
}

/// One output file of `get_cluster_sequence`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClusterFile {
    /// `<cluster id>.fasta`, except `NA.fasta`, renamed `unclustered.fasta`.
    pub file_name: String,
    pub content: String,
}

/// `get_cluster_sequence`: the contigs of `contig_file` split by the
/// cluster of `contig_info.txt` (column 5). Files are listed in order of
/// first appearance in `contig_info.txt`.
pub fn cluster_sequences(contig_info: &str, contig_file: &str) -> Vec<ClusterFile> {
    let mut order: Vec<String> = Vec::new();
    let mut index_of: HashMap<String, usize> = HashMap::new();
    let mut cluster_of: HashMap<String, usize> = HashMap::new();
    for line in perl_lines(contig_info).skip(1) {
        let fields = split_tabs(line);
        let contig = fields.first().copied().unwrap_or("");
        let cluster = fields.get(4).copied().unwrap_or("");
        let index = *index_of.entry(cluster.to_string()).or_insert_with(|| {
            order.push(cluster.to_string());
            order.len() - 1
        });
        cluster_of.insert(format!(">{contig}"), index);
    }
    let mut contents = vec![String::new(); order.len()];
    let mut current: Option<usize> = None;
    for line in perl_lines(contig_file) {
        if line.contains('>') {
            let id = split_whitespace_perl(line).first().copied().unwrap_or("");
            current = cluster_of.get(id).copied();
            if let Some(index) = current {
                contents[index].push_str(id);
                contents[index].push('\n');
            }
        } else if let Some(index) = current {
            contents[index].push_str(line);
        }
    }
    order
        .into_iter()
        .zip(contents)
        .map(|(cluster, content)| ClusterFile {
            file_name: if cluster == "NA" {
                "unclustered.fasta".to_string()
            } else {
                format!("{cluster}.fasta")
            },
            content,
        })
        .collect()
}

/// `convert_sigma_to_opera_lg_coverage_file` (single-genome mode).
pub fn sigma_to_opera_lg_coverage(sigma_contigs: &str) -> String {
    const READ_SIZE_FACTOR: f64 = 101.0;
    let mut out = String::new();
    let mut lines = perl_lines(sigma_contigs).skip(1);
    while let Some(line) = lines.next() {
        let fields = split_tabs(line);
        let id = fields.first().copied().unwrap_or("");
        let length = fields.get(1).copied().unwrap_or("");
        let read_count = lines.next().unwrap_or("");
        let coverage = READ_SIZE_FACTOR * perl_float(read_count) / perl_float(length);
        out.push_str(&format!("{id}\t{length}\t{}\n", perl_number(coverage)));
        lines.next();
    }
    out
}

fn perl_float(text: &str) -> f64 {
    let text = text.trim();
    let end = text
        .char_indices()
        .take_while(|(i, c)| {
            c.is_ascii_digit() || *c == '.' || (*i == 0 && (*c == '-' || *c == '+'))
        })
        .map(|(i, c)| i + c.len_utf8())
        .last()
        .unwrap_or(0);
    text[..end].parse().unwrap_or(0.0)
}

/// Perl's default number formatting (C's `%.15g`): 15 significant digits,
/// scientific notation when the (rounded) exponent is below -4 or at least
/// 15, trailing zeros removed. Limitation: Perl prints an integral value it
/// holds as an integer (IV) with all its digits (e.g. 76715762914796208
/// rather than 7.67157629147962e+16); that depends on how Perl computed the
/// value and only matters from 1e15 up, far above the lengths and counts
/// of this pipeline.
pub fn perl_number(value: f64) -> String {
    if value.is_infinite() {
        return if value > 0.0 {
            "Inf".into()
        } else {
            "-Inf".into()
        };
    }
    if value.is_nan() {
        return "NaN".into();
    }
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0".into()
        } else {
            "0".into()
        };
    }
    // rounded to 15 significant digits: d.dddddddddddddde<exp>
    let sci = format!("{value:.14e}");
    let (mantissa, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    if !(-4..15).contains(&exp) {
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        return format!(
            "{mantissa}e{}{:02}",
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        );
    }
    let decimals = usize::try_from(14 - exp).unwrap_or(0);
    let fixed = format!("{value:.decimals$}");
    if fixed.contains('.') {
        fixed
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    } else {
        fixed
    }
}

/// A line of `time.log` (`write_time`).
pub fn time_log_line(step: &str, seconds: u64) -> String {
    format!("{step}\t{seconds}\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn utils_config_order_and_missing_keys() {
        let values: HashMap<&str, String> = [
            ("NUM_PROCESSOR", "4".to_string()),
            ("ILLUMINA_READ_1", "/r1".to_string()),
            ("OUTPUT_DIR", "/out".to_string()),
            ("KRAKEN_DB", "ignored".to_string()),
        ]
        .into_iter()
        .collect();
        assert_eq!(
            utils_config(&values),
            "ILLUMINA_READ_1 /r1\nOUTPUT_DIR /out\nNUM_PROCESSOR 4\n"
        );
    }

    #[test]
    fn cov_estimate_additions_match_the_reference() {
        let m = "/run/intermediate_files/read_mapping";
        assert_eq!(
            cov_estimate_config_additions(
                None,
                "/r1.fq",
                m,
                "/oms//tools_opera_ms/",
                "/oms//tools_opera_ms/"
            ),
            format!(
                "\nLIB {m}/opera.bam\n\nMAPPING_FILES {m}/opera.bam\n\
                 SAMTOOLS_DIR /oms//tools_opera_ms/\nTOOLS_DIR /oms//tools_opera_ms/\n"
            )
        );
        assert!(
            cov_estimate_config_additions(Some("/c.fa"), "a,b", m, "s", "t")
                .starts_with("\nCONTIGS_FILE /c.fa\n\nLIB ")
        );
    }

    #[test]
    fn linear_fasta() {
        assert_eq!(
            multi_to_linear_fa(">a x\nAC\nGT\n>b\nTT\n"),
            // with the final newline
            ">a x\nACGT\n>b\nTT\n"
        );
        assert_eq!(multi_to_linear_fa(""), "");
    }

    #[test]
    fn stats() {
        let info = "CONTIG_ID\tLENGTH\n\
                    c1\t600000\tx\n\
                    c2\t300000\tx\n\
                    c3\t1200000\tx\n";
        let text = assembly_stats(info, "Sun Oct 25 08:45:47 2026");
        assert_eq!(
            text,
            "\n[Sun Oct 25 08:45:47 2026]\tAssembly stats\nNumber of contigs: 3 \n\
             Assembly size: 2100000 bp\nMax contig size: 1200000 bp\n\
             Contig(s) longer than 1Mbp: 1 \nContig(s) longer than 500kbp: 2\n\
             Contig(s) longer than 100kbp: 3\nContig N50: 1200000 bp\n"
        );
    }

    #[test]
    fn final_contigs_from_scaffolds() {
        let info = "CONTIG_ID\tLENGTH\tS\tL\tCLUSTER\tX\n\
                    opera_contig_1\t4\t1.5\t2.5\tcluster_1\tE\n";
        let cname = "opera_scaffold_1\topera_contig_1\n";
        let filled = ">opera_scaffold_1 length: 4\nACGT\n>unknown\nGG\n";
        assert_eq!(
            final_contigs(info, cname, filled),
            ">opera_contig_1 length: 4 cov_short: 1.5 cov_long: 2.5 cluster: cluster_1 \nACGT\n\
             > \nGG\n"
        );
    }

    #[test]
    fn clusters_in_order_of_contig_info() {
        // 8 columns as in contig_info.txt (with 5, the cluster would keep
        // the newline, as in Perl)
        let info = "ID\tL\tS\tL\tCLUSTER\tSP\tN\tREF\n\
                    c1\t1\t1\t1\tcluster_2\tNA\tNA\tNA\n\
                    c2\t1\t1\t1\tNA\tNA\tNA\tNA\n\
                    c3\t1\t1\t1\tcluster_2\tNA\tNA\tNA\n";
        let fasta = ">c1 info\nAA\n>c9 x\nTT\n>c2\nCC\n>c3\nGG";
        let files = cluster_sequences(info, fasta);
        assert_eq!(
            files,
            [
                ClusterFile {
                    file_name: "cluster_2.fasta".into(),
                    content: ">c1\nAA\n>c3\nGG".into()
                },
                ClusterFile {
                    file_name: "unclustered.fasta".into(),
                    content: ">c2\nCC\n".into()
                },
            ]
        );
    }

    #[test]
    fn perl_numbers() {
        assert_eq!(perl_number(25.0), "25");
        assert_eq!(perl_number(0.1 + 0.2), "0.3");
        assert_eq!(perl_number(2.0 / 3.0), "0.666666666666667");
        assert_eq!(perl_number(101.0 * 7.0 / 3.0), "235.666666666667");
        assert_eq!(perl_number(1e-7), "1e-07");
        // printf "%.15g" (C) for the same values
        assert_eq!(perl_number(9.32036804977833e-05), "9.32036804977833e-05");
        assert_eq!(perl_number(0.0001), "0.0001");
        assert_eq!(perl_number(0.00012345), "0.00012345");
        assert_eq!(perl_number(1e15), "1e+15");
        assert_eq!(perl_number(123456789012345.0), "123456789012345");
        assert_eq!(perl_number(-2.5), "-2.5");
        assert_eq!(perl_number(1e100), "1e+100");
    }

    /// Opt-in: `PERL_NUMBER_CASES=<file>` with lines `<value>\t<Perl's "$v">`
    /// (generated with `perl -ne 'chomp; my $v = $_ + 0.0; $v = $v * 1.0;
    /// printf("%.17g\t%s\n", $v, "$v")'`: the first column is the exact double
    /// Perl holds, its own string-to-number conversion being inexact at
    /// extreme exponents).
    #[test]
    fn perl_number_matches_perl_on_cases() {
        let Ok(path) = std::env::var("PERL_NUMBER_CASES") else {
            return;
        };
        let text = std::fs::read_to_string(path).unwrap();
        let mut failures = Vec::new();
        for line in text.lines() {
            let (input, expected) = line.split_once('\t').unwrap();
            let value: f64 = input.parse().unwrap();
            // Perl may hold an integral value as an integer (IV / UV, up to
            // 2^64) and print
            // all its digits: not decidable from the value (see perl_number)
            if value.fract() == 0.0 && value.abs() >= 1e15 && value.abs() < 18446744073709551616.0 {
                continue;
            }
            if perl_number(value) != expected {
                failures.push(format!(
                    "{input}: {} vs Perl {expected}",
                    perl_number(value)
                ));
            }
        }
        assert!(
            failures.is_empty(),
            "{} failures, e.g. {:?}",
            failures.len(),
            &failures[..failures.len().min(5)]
        );
    }

    #[test]
    fn perl_whitespace_split() {
        assert_eq!(split_whitespace_perl(" a b\n"), ["", "a", "b"]);
        assert_eq!(split_whitespace_perl(">c1\tx\n"), [">c1", "x"]);
        assert_eq!(split_tabs("a\tb\t\n"), ["a", "b", "\n"]);
    }
}
