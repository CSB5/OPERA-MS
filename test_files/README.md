# Test data

A small simulated metagenome for a quick check of an installation:
three genomes of the CAMI III C3Toy sample_0 gold standard (8.1 Mb): one
single-contig genome (ASV579.0) and a pair of strains of one species
(ASV507.11 and ASV507.9, 97 % ANI).

- `R1.fastq.gz`, `R2.fastq.gz`: 267,124 read pairs, 2x150 (art_modern,
  HiSeq2500 profile, fragments 270 +- 27), depth 12x, 12x and 6x.
- `long_read.fastq.gz`: 3,186 nanopore reads (Badread, identity about 87 %),
  depth 6x, 6x and 4x.
- `contigs.fasta`: the MEGAHIT 1.2.9 assembly of the short reads.

Random seed: 7.
