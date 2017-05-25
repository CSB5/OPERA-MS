# Introduction - Under development
OPERA-MS is a long read metagenomic scaffolding pipline that takes in a set of contigs in addition with __both short reads and long reads__ to output near-complete individual microbial genomes in your environmental sample. 

It uses the following strategy: a graph partitioning tool called Sigma is used to decompose the metagenomic scaffolding problem into distinct single genome scaffolding problems that are then solved by the single genome scaffolder [OPERA-LG](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0951-y).

# Installation - Under development
### Dependencies

We require the following dependencies:

1) [samtools](https://github.com/samtools/samtools) - (version 0.1.19 or below).
2) [bwa](https://github.com/lh3/bwa).
3) [blasr](https://github.com/PacificBiosciences/blasr) - (version 5.1 and above).
4) EMBOSS water - Smith-Waterman alignment tool.
5) [vsearch](https://github.com/torognes/vsearch).
6) [graphmap](https://github.com/isovic/graphmap).
7) [pbdagcon](https://github.com/PacificBiosciences/pbdagcon).

These must either be in your PATH or the path must be specified (see the Config file).

### Building

- Download the OPERA-MS and unzip it into a specified directory or `git clone https://github.com/CSB5/OPERA-MS.git`.

- `cd` into the OPERA-MS/ folder. For example, `cd /home/usr/OPERA-MS`.

- Run the `make` command.

# Running OPERA-MS

### Inputs
OPERA-MS takes in 3 inputs

1) An assembled contigs/scaffolds file in multi-fasta format (e.g. test_dataset/contigs.fa).
2) A long read file to be used in scaffolding (e.g. test_dataset/long_read_1.fa).
3) Paired end reads to be used in scaffolding (e.g. test_dataset/lib_1_1.fa, test_dataset/lib_1_2.fa).

## Executing OPERA-MS

To run OPERA-MS, simply run `perl /path/to/OPERA-MS/OPERA-MS.pl <config file>`. The location of the test data should be specified in the config file. See below for more information on the format of the config file. An example .config file is bundled with OPERA-MS.

## Specification for the Configuration file

Configuration files must be formatted in the form :

~~~~
<OPTION1> <VALUE1>
<OPTION2> <VALUE2>
...
<OPTION2> <VALUE3>
~~~~

This is an example of a config file :

~~~~
#This is a comment. 

CONTIGS_FILE /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/SOFWARE/OPERA-MS/test_files/final.contigs_soap.fa #We can use absolute or relative paths. #This is an absolute path.
OUTPUT_DIR opera_ms_output/
LONG_READ_FILE test_files/POOL.fa
ILLUMINA_READ_1 test_files/mock1000.R1.fastq.gz
ILLUMINA_READ_2 test_files/mock1000.R2.fastq.gz

#SAMTOOLS_DIR/home/usr/software/ # comment out if samtools is in PATH.
#BLASR_DIR /home/usr/software # comment out if blasr is in PATH.
#SHORT_READ_TOOL_DIR /home/usr/software # comment out if bwa is in PATH.
#GRAPHMAP_DIR /home/usr/software  # comment out if graphmap is in PATH.
#WATER_DIR /home/usr/software  # comment out if water is in PATH.
#PBDAGCON_DIR /home/usr/software  # comment out if pbdagcon is in PATH.
#VSEARCH_DIR /home/usr/software  # comment out if vsearch is in PATH.

OPERA_VERSION OPERA-LG_v2.1.0
NUM_PROCESSOR 20
CONTIG_LEN_THR 500
CONTIG_EDGE_LEN 80
CONTIG_WINDOW_LEN 340
PDIST_TYPE NegativeBinomial 
KMER_SIZE 60
~~~~

### Options 
All relative paths are relative to the current working directory of your terminal. All paths can be chosen to be either relative or absolute.

- **CONTIGS_FILE** : `path/to/contigs.fa` - A path to the contigs file.

- **OUTPUT_DIR** : `path/to/results` - Where the final results of OPERA-MS will go.

- **LONG_READ_FILE** : `path/to/long-read.fa` - A path to the long read file.

- **ILLUMINA_READ_1** : `path/to/illum_read1.fa` - A path to the first illumina read file.

- **ILLUMINA_READ_2** : `path/to/illum_read2.fa` - A path to the second complement illumina read file.

- **(tool)_DIR** : `path/to/tool_directory` - A path to the __directory containing__ the executable file of the specific tool : e.g. blasr, bwa, vsearch. If commented out, the tool will be assumed to be in PATH.

- **OPERA_VERSION** : `OPERA-LG_v2.x.x` - The version of OPERA-LG used. User should not need to change this.

- **NUM_PROCESSOR** : The number of processors that this pipeline will use.

- **CONTIG_LEN_THR** : `default : 0` - Threshold for contig clustering.

- **CONTIG_EDGE_LEN** : `default : 0` - TODO

- **CONTIG_WINDOW_LEN** : `default : 0` - Used in SIGMA for scoring. If set to 0, contig based scoring is used instead of window scoring.

- **PDIST_TYPE** : `Poisson` or `NegativeBinomial`- The read count probability distribution. Only Poisson and NegativeBinomial are supported currently.

- **KMER_SIZE** : `default : 60` - The value of kmer used to produce the assembled contigs/scaffolds.


## Outputs

The outputs will be in OUTPUT_DIR/contigs/scaffolds-long-reads. TODO




