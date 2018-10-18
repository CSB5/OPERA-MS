# Introduction 
OPERA-MS is first in class hybrid metagenomic assembler which combines the advantages of short and long-read technologies, providing highly contiguous assemblies with low base-pair error rate.

OPERA-MS uses the following strategy: (i) Assembling preliminary contigs using the short-read metagenomic assembler Megahit and overlaying long-read information to construct an assembly graph of all genomes, (ii) disentangling genomes at the species level using a reference-free clustering approach and augmenting clusters using a reference guided approach, (iii) identifying sub-species level clusters and scaffolding and gap-filling them using [OPERA-LG](https://github.com/CSB5/OPERA-MS/tree/master/OPERA-LG). 

# Installation

#### Perl Modules
- [Switch](http://search.cpan.org/~chorny/Switch-2.17/Switch.pm)

- [Statistics::Basic](http://search.cpan.org/~jettero/Statistics-Basic-1.6611/lib/Statistics/Basic.pod)

- [Statistics::R](https://metacpan.org/pod/Statistics::R)

To install OPERA-MS, download and unzip OPERA-MS to a specified directory or `git clone https://github.com/CSB5/OPERA-MS.git`.

All required tools are distributed pre-built. See below for more information on dependencies.

~~~~
cd /path/to/OPERA-MS
make
~~~~

This will build SIGMA (reference-free clustering approach), [OPERA-LG](https://github.com/CSB5/OPERA-MS/tree/master/OPERA-LG) (single genome scaffolder/gap-filler), and install the required perl modules.

### Testing Installation

A set of test files and a sample configuration file is provided to perform a test run of OPERA-MS. To test it, simply write the following commands: 
~~~~
cd /path/to/OPERA-MS
perl OPERA-MS.pl sample_config.config 2> log.err
~~~~
This will assemble a low diversity mock community in the folder __OPERA-MS/sample_output__. One should see the following statistics (name and size of the 10 longest assembled sequences) after completion:

~~~~
>strain1_opera_scaffold_1	5180338
>opera_scaffold_1	199681
>opera_scaffold_2	188904
>opera_scaffold_3	159657
>opera_scaffold_4	147456
>opera_scaffold_5	140923
>opera_scaffold_6	128548
>opera_scaffold_7	124292
>opera_scaffold_8	119336
>opera_scaffold_9	115152
~~~~

# Running OPERA-MS

### Inputs
OPERA-MS takes short and long reads fastq files as inputs.

1) A long read fastq files to be used in scaffolding and gap-filling (e.g. sample_files/sample_long_read.fastq).
2) Paired end reads fastq to be used for contig assembly and scaffolding (e.g. test_dataset/sample_R1.fastq.gz, test_dataset/sample_R2.fastq.gz).
3) Optionally, contigs assembled from short reads can be provided (e.g. sample_files/sample_contigs.fasta).

### Executing OPERA-MS

OPERA-MS requires the specfication of a configuration file that will indicate the path to the input files and the options used during the assembly steps. To run OPERA-MS, write `perl /path/to/OPERA-MS/OPERA-MS.pl <config file> 2> log.err`. 

### Configuration File Specification

The configuration file must be formatted as follow:

~~~~
#One space between OPTION and VALUE
<OPTION1> <VALUE1> 
<OPTION2> <VALUE2>
...
<OPTION2> <VALUE3>
~~~~

Here an example of a configuration file:

~~~~
#This is a comment. 
#We can use absolute or relative paths
OUTPUT_DIR opera_ms_output/ #Relative path from current working directory.
LONG_READ_FILE test_files/sample_long_read.fastq
ILLUMINA_READ_1 test_files/sample_R1.fastq.gz
ILLUMINA_READ_2 test_files/sample_R2.fastq.gz

STRAIN_CLUSTERING YES
NUM_PROCESSOR 20
CONTIG_LEN_THR 500
CONTIG_EDGE_LEN 80
CONTIG_WINDOW_LEN 340
KMER_SIZE 60
~~~~

### Options 
All paths are relative to the current working directory of your terminal. All paths can be chosen to be either relative or absolute.

- **OUTPUT_DIR**: `path/to/results` - Path to the directory where the final results of OPERA-MS will be outputted.

- **LONG_READ_FILE**: `path/to/long-read.fastq` - Path to the long read file.

- **ILLUMINA_READ_1**: `path/to/illum_read1.fastq.gz` - Path to the first illumina read file.

- **ILLUMINA_READ_2**: `path/to/illum_read2.fastq.gz` - Path to the second illumina read file.

- **NUM_PROCESSOR**: `default : 1` - The number of processors used.

- **CONTIG_LEN_THR**: `default: 500` - Threshold for contig clustering. Smaller contigs will not be considered for clustering.

- **STRAIN_CLUSTERING**: `YES/NO` - Indicate if the strain level clustering step is performed (YES) or skipped (NO) 

- **CONTIG_EDGE_LEN**: `default: 80` - When calculating coverage of contigs using SIGMA, this number of bases will not be considered from each end of the contig to avoid biases due to lower mapping efficiency at contig edges. 

- **CONTIG_WINDOW_LEN**: `340` - The window size in which the coverage estimation is performed. We recommend using CONTIG_LEN_THR - 2 * CONTIG_EDGE_LEN as the value.

- **KMER_SIZE**: `(positive integer)` - The value of kmer used to produce the assembled contigs/scaffolds.

- **CONTIGS_FILE**: `path/to/contigs.fa` - Path to the contigs file, if the short reads have already been assembled previously.

OPERA-MS can be used as complete pipeline, or can be used in a stepwise fashion:
- **SHORT_READ_ASSEMBLY**: short-read assembly using Megahit.
- **CONTIG_GRAPH**: sort-read and long-read mapping to the contig assembly to construct the assembly graph and estimate the contigs coverage.
- **CLUSTERING**: hierarchical clustering based on contig connections and coverage.
- **REF_CLUSTERING**: cluster merging based on similarity to a data-base of reference genomes (computation of `mash` distance and contig alignment to best reference using `mummer`).
- **STRAIN_CLUSTERING**: identification of species with multiple strains, strains deconvolution and independent assembly of each strain using OPERA-LG. Strain level scaffold assemblies are provided in the following files: __OUT_DIR/intermediate_files/strain_analysis/\*/\*/scaffoldSeq.fasta__.
- **ASSEMBLY**: assembly of single strain species clusters using OPERA-LG.
- **GAP_FILLING**: gap-filling of all the scaffold assemblies using long-read data and short contigs.
- **INFO**: generation of metagenome assembly statistics (see below for a detailed description).

Concatenating a **+** to a step name will allow to perfom all the steps following it. For example **REF_CLUSTERING+** will allow to run all the following steps after **REF_CLUSTERING**.

### Outputs

Outputs can be found in the specified OUTPUT_DIR and will contains all assembled sequences before and after gap-filling: __scaffoldSeq.fasta__ and __scaffoldSeq.fasta.filled__ respectively.
The file __scaffold_info.txt__ provides an overview of the assembled sequences according to the following features:
- **SEQ_ID**: sequence identifier. For single strain species sequences are name `opera_scaffold_X`. Sequences from multi-strain species are named `strainY_opera_scaffold_X` where `Y` indicate the strain ID.
- **LENGTH**: sequence length.
- **ARRIVAL_RATE**: median sequence arrival rate.
- **SPECIES**: putative species to which the assembled sequence belong to.
- **NB_STRAIN**: number of strains detected for the species.
- **REFERENCE_GENOME**: path to the closest reference genome present in the OPERA-MS data-base.
- **FRACTION_GENOME_COVERED**: fraction of the sequence that maps to the reference genome.
- **IDENTITY**: % identity between the sequence and the reference genome.

### Dependencies

We require the following dependencies:

1) [MegaHit](https://github.com/voutcn/megahit) - (version 1.0.4-beta)
2) [samtools](https://github.com/samtools/samtools) - (version 0.1.19 or below).
3) [bwa](https://github.com/lh3/bwa).
4) [blasr](https://github.com/PacificBiosciences/blasr) - (version 5.1 and above which uses '-' options).
5) [minimap2]( https://github.com/lh3/minimap2). (version 2.11-r797)
6) [Racon](https://github.com/isovic/racon) - (version 0.5.0)
7) [Mash](https://github.com/marbl/Mash).
8) [MUMmer](http://mummer.sourceforge.net/).

All software are packaged as pre-build with OPERA-MS expect MUMmer. Each binary is placed inside of the __utils__ folder.

If the pre-built tools do not work on the user's machine, OPERA-MS will check the user's PATH for the tool. However, the version that the user is using may be different than the one packaged. Alternatively, to specify a different directory for the dependency, edit the configuration file and add the line `(tool)_DIR /path/to/dir` as shown below:

- **(tool)_DIR** : `path/to/tool_directory` - A path to the __directory containing__ the executable file of the specific tool : e.g. blasr, bwa, racon. If commented out the tool within the utils/ directory will be used. 

For example, `BWA_DIR /usr/home/bwa_dir`.
