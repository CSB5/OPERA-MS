# Changes

## Unreleased (2.0.1-dev)

- **Speed:** the hierarchical clustering of SIGMA (step 4) is near linear
  in the size of the cluster tree instead of quadratic: hours or days on
  very large assemblies (millions of contigs, hundreds of genomes), now
  seconds or minutes. Same outputs.
- **Fix:** cluster length and read count totals are 64-bit (overflow on
  very large clusters).
- **Reproducibility:** the choice of reference genomes for a cluster is
  now the same from run to run. When more genomes than kept (5) were tied
  at the same distance, skani chose among them at random between runs;
  ties now go to the first genomes by accession. Final contigs can differ
  from 2.0.0 on samples with such ties. Sylph's genomes are also listed in
  sorted order.
- **Fix:** the reference clustering (step 5) no longer writes one file
  per cluster all at once (tens of thousands on large assemblies, beyond
  some file quotas): the cluster files go to skani 1,000 at a time and
  are removed after each batch; the reference genomes are sketched once.
  The alignments of the linked cluster pairs are kept in memory instead
  of two files per pair in `NUCMER_OUT/` (tens of thousands, never
  removed). Same outputs.
- **Output:** the per-cluster FASTA files (`opera_ms_clusters/all/`) are
  only written with the new option `--cluster-fasta`; by default the
  clusters are given by `contig_info.txt` (column `CLUSTER`) and described
  in `cluster_info.txt` (unchanged).
- **Disk:** the large files of `intermediate_files/read_mapping/` are
  gzipped (`opera.map.sort`, `opera.map.sort.status`, `edge_read_info.dat`,
  `pairedEdges`, `pairedEdges_above_thresh`, `prefiltered_*`,
  `below_thresh_*`: up to about 100 GB on large samples), about 3 times
  smaller, compressed in parallel. A run started with 2.0.0 can be
  resumed: the plain files are still read.
- **Disk:** the strain clustering (step 6) no longer writes `wrong` per
  species (the edges outside the species, never read: about the size of
  `pairedEdges` each time, hundreds of GB on large samples). A species'
  edges are one gzipped file (`pairedEdges.gz`) instead of six directories;
  a strain's edge files are written directly in `STRAIN_<n>/`, without the
  empty `lib.txt`, and removed once OPERA-LG has run: about half the files
  of `strain_analysis/`. Same outputs.
- **Fix:** the gap filling (step 8) no longer makes a directory per
  scaffold for the whole step (tens of thousands of files on large
  samples, beyond some file quotas): the read pools go to at most 256
  files (`READ_POOLS/`), the pre-consensus sequences to one file, and a
  scaffold's directory only exists while its racon job runs. Same outputs.
- **Reproducibility:** OPERA-LG no longer depends on memory addresses. Sets
  of contigs and edges were ordered by address and iterated where the order
  reaches the scaffolds, so on ties the results could depend on the memory
  allocator, the C library or the memory layout (reversing that order
  changed a test's scaffolds). They are ordered by creation now, which
  gives the same results as 2.0.0 on the samples tested.
- **Speed:** OPERA-LG's scaffolding search (steps 6 and 7) is about 1.7
  times faster: a contig's scaffold text is no longer parsed again on each
  search step, the unassigned contigs are saved and restored without one
  allocation per contig, and a set that was filled but never read is gone.
  Same outputs. A search that hits the time limit (`--opera-lg-timeout`)
  goes further before it.

## 2.0.0

A rewrite of OPERA-MS 0.9 in Rust, as a single `opera-ms` executable (the
C++ of SIGMA and OPERA-LG is compiled in). The algorithm is the same; on
simulated gut metagenomes (CAMI III) the assembly quality is on par with
0.9 + Pilon, about 7 times faster with a quarter of the memory.

- **Installation:** `make`; the external tools come from a conda environment
  (`environment.yml`) instead of being bundled. No Perl, R or Python.
- **Tools:** minibwa replaces bwa; minimap2 replaces blasr and MUMmer;
  polishing uses POLCA's method with parallel freebayes (`--polisher
  polypolish-pypolca`: Polypolish first) instead of Pilon.
- **Reference genomes:** sylph (`--sylph-db`) and genomes downloaded from
  NCBI replace the OPERA-MS database (`install-db`, `--genome-db`, Mash);
  `--genome-cache`, `--prepare-genomes`, `--offline`, `--strict-genomes`,
  `NCBI_API_KEY`.
- **New options:** `--opera-lg-timeout` (fallback to OPERA-LG-fast),
  `--polisher`.
- **Removed:** the configuration-file mode (and its start-stage argument),
  `--long-read-mapper` (blasr), `--short-read-assembler` (SPAdes),
  OPERA-MS-UTILS, the Docker image.
- **Behaviour:** every command is checked and any failure stops the run;
  gzipped long reads are accepted; the output doesn't depend on the number
  of threads (MEGAHIT aside); a run started again in the same output
  directory resumes at the first incomplete stage; `opera-ms.log` records every command; `benchmark.tsv`
  gives time and memory per step.
- **Fixes:** many small bugs of 0.9, e.g. the last read or cluster of a file
  being ignored, unfilled gaps written as `NA`, divisions by zero on empty
  clusters.
