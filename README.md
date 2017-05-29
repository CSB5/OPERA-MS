# UNDER DEVELOPMENT
# Introduction 
OPERA-MS is a long read metagenomic scaffolding pipline that takes in a set of contigs in addition with __both short reads and long reads__ to output near-complete individual microbial genomes in your environmental sample. 

It uses the following strategy: a graph partitioning tool called Sigma is used to decompose the metagenomic scaffolding problem into distinct single genome scaffolding problems that are then solved by the single genome scaffolder [OPERA-LG](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0951-y).

The long reads are then used fill the gaps between contigs to produce a final set of scaffolds.

# Installation

To install OPERA-MS, download and unzip OPERA-MS to a specified directory or `git clone https://github.com/CSB5/OPERA-MS.git`. All of the dependencies are distributed pre-built.

The test script test-install.pl is provided. To test all dependencies, 
~~~~
cd /path/to/OPERA-MS
perl test-install.pl
~~~~
This will assemble a small mock-genome and will display any problems regarding dependencies.

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

CONTIGS_FILE /home/usr/OPERA-MS/test_files/final.contigs_soap.fa #We can use absolute or relative paths. #This is an absolute path.
OUTPUT_DIR opera_ms_output/
LONG_READ_FILE test_files/POOL.fa
ILLUMINA_READ_1 test_files/mock1000.R1.fastq.gz
ILLUMINA_READ_2 test_files/mock1000.R2.fastq.gz

NUM_PROCESSOR 20
CONTIG_LEN_THR 500
CONTIG_EDGE_LEN 80
CONTIG_WINDOW_LEN 340
KMER_SIZE 60
~~~~

### Options 
All relative paths are relative to the current working directory of your terminal. All paths can be chosen to be either relative or absolute.

- **CONTIGS_FILE** : `path/to/contigs.fa` - A path to the contigs file.

- **OUTPUT_DIR** : `path/to/results` - Where the final results of OPERA-MS will go.

- **LONG_READ_FILE** : `path/to/long-read.fa` - A path to the long read file.

- **ILLUMINA_READ_1** : `path/to/illum_read1.fa` - A path to the first illumina read file.

- **ILLUMINA_READ_2** : `path/to/illum_read2.fa` - A path to the second complement illumina read file.

- **(tool)_DIR** : `path/to/tool_directory` - A path to the __directory containing__ the executable file of the specific tool : e.g. blasr, bwa, vsearch. If commented out the tool within the utils/ directory will be used. 

- **NUM_PROCESSOR** : The number of processors that this pipeline will use.

- **CONTIG_LEN_THR** : `default : 500` - Threshold for contig clustering. Smaller contigs will not be considered for clustering.

- **CONTIG_EDGE_LEN** : `default : 80` - When calculating coverage of contigs using SIGMA, this number of bases will not be considered from each end of the contig. This is to remove to biases due to reads not matching the edges of contigs. 

- **CONTIG_WINDOW_LEN** : `(positive integer)` - Used in SIGMA for scoring contigs. We recommend using CONTIG_LEN_THR - 2 * CONTIG_EDGE_LEN as a value.

- **KMER_SIZE** : `(positive integer)` - The value of kmer used to produce the assembled contigs/scaffolds.


## Outputs

The outputs will be in OUTPUT_DIR. The scaffolds file before the secondary gap filling procedure is denoted __scaffoldSeq.fasta__. The scaffolds file after filling is denoted __scaffoledSeq.fasta.filled__. Intermediate files are are inside foe folder __intermediate_files__. 

## Dependencies

We require the following dependencies:

1) [samtools](https://github.com/samtools/samtools) - (version 0.1.19 or below).
2) [bwa](https://github.com/lh3/bwa).
3) [blasr](https://github.com/PacificBiosciences/blasr) - (version 5.1 and above).
4) [EMBOSS](http://emboss.sourceforge.net/download/) water - Smith-Waterman alignment tool.
5) [vsearch](https://github.com/torognes/vsearch).
6) [graphmap](https://github.com/isovic/graphmap).
7) [pbdagcon](https://github.com/PacificBiosciences/pbdagcon).

These are pre-built and packaged with OPERA-MS. Each binary is placed inside of the utils folder.

If the pre-built tools do not work on the user's machine, then the user can build the dependencies themselves. To use a different dependency, edit the config file and add the line `(tool)_DIR /path/to/dir`. For example, `BLASR_DIR /usr/home/blasr_dir`.




