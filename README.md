# Introduction 
OPERA-MS is first in class hybrid metagenomic assembler which combines the advantages of short and long-read technologies, providing highly contiguous assemblies with low base-pair error rate.

OPERA-MS uses the following strategy: (i) Assembling preliminary contigs using the short-read metagenomic assembler Megahit and overlaying long-read information to construct an assembly graph of all genomes, (ii) disentangling genomes at the species level using a reference-free clustering approach and augmenting clusters using a reference guided approach, (iii) identifying sub-species level clusters and scaffolding and gap-filling them using OPERA-LG. 

# Installation

To install OPERA-MS run the following commands:

~~~~
git clone https://github.com/CSB5/OPERA-MS.git
cd /path/to/OPERA-MS
make
~~~~

### Testing Installation

A set of test files and a sample configuration file is provided to test out the OPERA-MS pipeline. Simply run the following commands: 
~~~~
cd /path/to/OPERA-MS
perl OPERA-MS.pl sample_config.config 2> log.err
diff contig_info.txt sample_output/contig_info.txt
~~~~
This will assemble a low diversity mock community in the folder __OPERA-MS/sample_output__. 

# Running OPERA-MS

### Inputs
OPERA-MS takes short and long reads as inputs.

1) A long read fastq file obtained from either Oxford Nanopore, PacBio or TruSeq Synthetic Long Read sequencing technologies (e.g. sample_files/sample_long_read.fastq).
2) Short read paired end read fastq files (e.g. sample_files/sample_R1.fastq.gz and sample_files/sample_R2.fastq.gz).
3) Optionally, contigs assembled from short reads can be provided (e.g. sample_files/sample_contigs.fasta).

### Executing OPERA-MS

OPERA-MS requires the specfication of a configuration file that indicates the path to the input files and the options used for the assembly. To run OPERA-MS, write `perl /path/to/OPERA-MS/OPERA-MS.pl <config file>`. 

### Required parameters
All paths are relative to the current working directory of your terminal. All paths can be chosen to be either relative or absolute.

- **OUTPUT_DIR** : `path/to/results` - Where the final results of OPERA-MS will be outputted

- **LONG_READ** : `path/to/long-read.fa` - Path to the long read file

- **ILLUMINA_READ_1** : `path/to/illum_read1.fa` - Path to the first illumina read file

- **ILLUMINA_READ_2** : `path/to/illum_read2.fa` - Path to the second illumina read file

### Optional parameters 

- **NUM_PROCESSOR** : `default : 1` - The number of used processors

- **CONTIG_LEN_THR** : `default: 500` - Threshold for contig clustering, smaller contigs will not be considered for clustering

- **STRAIN_CLUSTERING** : `YES/NO` - Indicate if the strain level clustering step is performed (YES) or skipped (NO) 

- **CONTIG_EDGE_LEN** : `default: 80` - When calculating coverage of contigs using SIGMA, this number of bases will not be considered from each end of the contig, to avoid biases due to lower mapping efficiency at contig edges

- **CONTIG_WINDOW_LEN** : `340` - The window size in which the coverage estimation is performed. We recommend using CONTIG_LEN_THR - 2 * CONTIG_EDGE_LEN as the value

- **KMER_SIZE** : `(positive integer)` - The value of kmer used to produce the assembled contigs/scaffolds

- **CONTIGS_FILE** : `path/to/contigs.fa` - Path to the contigs file, if the short reads have already been assembled previously

### Outputs

Outputs can be found in the specified OUTPUT_DIR and will contains all assembled contigs __contig.fasta__.
The file __contig_info.txt__ provides an overview of the assembled contigs according to the following features:
- **CONTIG_ID** : contig identifier. Single strain species contigs are name `opera_contig_X`. Contigs from multi-strain species are named `strainY_opera_contig_X` where `Y` indicate the strain ID
- **LENGTH** : contig length
- **ARRIVAL_RATE** : median contig short read arrival rate
- **SPECIES** : putative species to which the assembled contig belong to
- **NB_STRAIN** : number of strains detected for the species
- **REFERENCE_GENOME** : path to the closest reference genome present in the OPERA-MS database

Strain level scaffold assemblies are provided in the following files: __OUT_DIR/intermediate_files/strain_analysis/*/*/scaffoldSeq.fasta__.

### Dependencies

We require the following software to be functional:
1) [MEGAHIT](https://github.com/voutcn/megahit) - (tesated with version 1.0.4-beta)
2) [samtools](https://github.com/samtools/samtools) - (version 0.1.19 or below).
3) [bwa](https://github.com/lh3/bwa).
4) [blasr](https://github.com/PacificBiosciences/blasr) - (version 5.1 and above which uses '-' options).
5) [minimap2]( https://github.com/lh3/minimap2). (tested with version 2.11-r797)
6) [Racon](https://github.com/isovic/racon) - (version 0.5.0)
7) [Mash](https://github.com/marbl/Mash).
8) [MUMmer](http://mummer.sourceforge.net/).
All softwares are packaged as pre-build with OPERA-MS. Each binary is placed inside of the __utils__ folder.

If a pre-built software does not work on the user's machine, OPERA-MS will check if the tool is present in the user's PATH. However, the version of the software may be different than the one packaged. Alternatively, to specify a different directory for the dependency, a link to the software may be placed in the  __utils__ folder.
