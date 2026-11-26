# OPERA-MS

OPERA-MS is a hybrid metagenomic assembler: it combines short reads, for
base accuracy, with long reads, for contiguity. It starts from a short-read
assembly (MEGAHIT), clusters the contigs into species and strains using
coverage, read links and reference genomes, then scaffolds, gap-fills and
polishes each genome. It can assemble near-complete genomes with as little
as 9x long-read coverage, and separates strains of the same species.

**IMPORTANT**. This version is a rewrite of OPERA-MS 0.9 in Rust; see [CHANGES.md](CHANGES.md). **This version is still experimental**. While our tests indicate similar or better results on test datasets, we can't yet confirm performance. If you want to use the previous stable version, please switch to the `Master` release on this GitHub page.

## Installation

Linux x86_64. `opera-ms` is a single executable; the tools it runs come
from a conda environment.

### With conda (Bioconda)

```
conda create -n opera-ms -c conda-forge -c bioconda opera-ms
conda activate opera-ms
opera-ms check-dependency              # lists the tools found, with their versions
```

### Installation from source

```
git clone https://github.com/CSB5/OPERA-MS.git
cd OPERA-MS
conda env create -f environment.yml    # the opera-ms environment
conda activate opera-ms
make                                   # needs Rust (https://rustup.rs) and a C++ compiler
./opera-ms check-dependency            # lists the tools found, with their versions
```

## Quick test

About a minute with 4 threads. The test data is in the repository's
`test_files/` (see its `README.md`), not in the conda package:

```
git clone https://github.com/CSB5/OPERA-MS.git    # if not done already
cd OPERA-MS/test_files
opera-ms --contig-file contigs.fasta --short-read1 R1.fastq.gz --short-read2 R2.fastq.gz \
    --long-read long_read.fastq.gz --no-ref-clustering --num-processors 4 --out-dir RESULTS
```

With a build from source, use `../opera-ms` instead of `opera-ms`.

## Reference database

Reference-based clustering finds the species of the sample with
[sylph](https://github.com/bluenote-1577/sylph), against a sylph database
of GTDB genomes: download a pre-built one (e.g. GTDB r232, `.syl2db`) from
[sylph's pre-built databases](https://sylph-docs.github.io/pre%E2%80%90built-databases/)
and give it with `--sylph-db`. Without a database, run with
`--no-ref-clustering`.

The genomes of the species found are downloaded from NCBI into
`--genome-cache`, which keeps them between runs (an NCBI API key in
`NCBI_API_KEY` raises NCBI's request limit). The cache is a folder of
`<accession>.fna.gz` files, so it can also be filled beforehand, for
example from GTDB's archive of representative genomes (same GTDB release
as the sylph database):

```
tar -xzf gtdb_genomes_reps_r232.tar.gz
mkdir -p GENOMES
find gtdb_genomes_reps_r232 -name '*_genomic.fna.gz' | while read -r f; do
    ln -s "$PWD/$f" "GENOMES/$(basename "$f" _genomic.fna.gz).fna.gz"
done
```

On machines without internet, either fill the cache that way, or run
`opera-ms --prepare-genomes` (same reads, `--sylph-db` and `--genome-cache`)
where there is internet; then run the assembly with `--offline`.

## Usage

```
opera-ms --short-read1 R1.fastq.gz --short-read2 R2.fastq.gz --long-read long_reads.fastq.gz \
    --sylph-db gtdb.syl2db --genome-cache GENOMES --num-processors 16 --out-dir RESULTS
```

Reads are FASTA or FASTQ, plain or gzip-compressed. `--contig-file` gives
a short-read assembly made beforehand (otherwise MEGAHIT runs).
`opera-ms --help` lists every option.

**Reproducibility.** The output doesn't depend on `--num-processors`,
except through MEGAHIT, whose assembly can change with the number of
threads: give the contigs with `--contig-file` for fully reproducible runs.

## Output

In the output directory:

- `contigs.polished.fasta`: the assembly (`contigs.fasta` before polishing);
- `assembly.stats`: assembly size, N50, longest contig, etc.;
- [`contig_info.txt`](https://github.com/CSB5/OPERA-MS/wiki/Contig-info-file-description): one line per contig;
- [`cluster_info.txt`](https://github.com/CSB5/OPERA-MS/wiki/Cluster-info-file-description):
  one line per strain-level cluster; the cluster of each contig is the
  `CLUSTER` column of `contig_info.txt`. With `--cluster-fasta`, also
  `opera_ms_clusters/all/`: one FASTA per cluster (`unclustered.fasta` for
  the contigs without a cluster). The clusters are conservative; contigs can
  be binned further with tools such as MetaBAT2;
- `opera-ms.log`: every step and command; `benchmark.tsv`: time and memory
  per step.
- `intermediate_files/`: the files of each step. Some are removed once
  used; `--keep-intermediate` keeps them, for inspection: all the files of
  the strains' OPERA-LG runs (`strain_analysis/*/STRAIN_<n>/`, otherwise
  only the scaffolds, contigs and logs) and the polishing alignments
  (`polished_assembly/assembly.bam`).

## Resources

On simulated human gut samples (CAMI III, 5 Gb of short reads and 4.5 Gb of
long reads), with 32 threads: about 15 minutes and 13 GB of RAM, plus the
short-read assembly when MEGAHIT runs.

## Citation

Bertrand D., Shaw J., Kalathiyappan M. et al. Hybrid metagenomic assembly
enables high-resolution analysis of resistance determinants and mobile
elements in human microbiomes. *Nature Biotechnology* 37, 937-944 (2019).
[doi:10.1038/s41587-019-0191-2](https://doi.org/10.1038/s41587-019-0191-2)

## Contact

- Jean-Sébastien Gounot <jean_sebastien@a-star.edu.sg>
- Niranjan Nagarajan <niranjan@nus.edu.sg>
