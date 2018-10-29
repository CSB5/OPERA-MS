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
If you encounter any problem during the installation, please see the **Dependencies** Section. 

### Testing Installation

A set of test files and a sample configuration file is provided to test out the OPERA-MS pipeline. To run OPERA-MS on the test dataset, simply use the following commands: 
~~~~
cd /path/to/OPERA-MS
perl OPERA-MS.pl sample_config.config 2> log.err
diff contig_info.txt sample_output/contig_info.txt
~~~~
This will assemble a low diversity mock community in the folder __OPERA-MS/sample_output__. 

# Running OPERA-MS

OPERA-MS requires the specfication of a configuration file that indicates the path to the input files and the options used for the assembly.
The configuration file is formatted as follow:

~~~~
#One space between OPTION and VALUE
<OPTION1> <VALUE1> 
<OPTION2> <VALUE2>
...
<OPTION2> <VALUE3>
~~~~

### Essential parameters

- **OUTPUT_DIR** : `path/to/results` - Directory where OPERA-MS results will be outputted

- **LONG_READ** : `path/to/long-read.fq` - Path to the long-read fastq file obtained from either Oxford Nanopore, PacBio or TruSeq Synthetic Long Read sequencing

- **ILLUMINA_READ_1** : `path/to/illum_read1.fq.gz` - Path to the first Illumina read file

- **ILLUMINA_READ_2** : `path/to/illum_read2.fq.gz` - Path to the second Illumina read file

### Optional parameters 

- **CONTIGS_FILE** : `path/to/contigs.fa` - Path to the contigs file, if the short-reads have been assembled previously

- **NUM_PROCESSOR** : `default : 1` - The number of used processors

- **LONG_READ_MAPPER** `default: blasr` - Software used for long-read mapping blasr/minimap2

- **STRAIN_CLUSTERING** : `default: YES` - Indicate if the strain level clustering step is performed (YES) or skipped (NO)

- **CONTIG_LEN_THR** : `default: 500` - Threshold for contig clustering, smaller contigs will not be considered for clustering

- **CONTIG_EDGE_LEN** : `default: 80` - When calculating contig coverage, number of bases filtered out from each contig ends, to avoid biases due to lower mapping efficiency

- **CONTIG_WINDOW_LEN** : `default: 340` - The window size in which the coverage estimation is performed. We recommend using CONTIG_LEN_THR - 2 * CONTIG_EDGE_LEN as the value

- **KMER_SIZE** : `default: 60` - The kmer value used to produce the assembled contigs


### Outputs

Outputs can be found in the specified OUTPUT_DIR.
The file __contig.fasta__ contains the assembled contigs, and __assembly.stats__ gives overall assembly statistics (e.g. assembly size, N50, longest scaffold ...).
__contig_info.txt__ provides a detailed overview of the assembled contigs according to the following features:
- **CONTIG_ID** : contig identifier. Single strain species contigs are named `opera_contig_X`. Contigs from multi-strain species are named `strainY_opera_contig_X`, where `Y` indicate the strain ID
- **LENGTH** : contig length
- **ARRIVAL_RATE** : median contig short-reads arrival rate
- **SPECIES** : putative species to which the assembled contig belong to
- **NB_STRAIN** : number of strains detected for the species
- **REFERENCE_GENOME** : path to the closest reference genome present in the OPERA-MS database

Finally, strain level scaffold assemblies can be found in the following files: __OUT_DIR/intermediate_files/strain_analysis/\*/\*/scaffoldSeq.fasta__.

# Dependencies

We require the following software to be functional:
1) [MEGAHIT](https://github.com/voutcn/megahit) - (tesated with version 1.0.4-beta)
2) [samtools](https://github.com/samtools/samtools) - (version 0.1.19 or below)
3) [bwa](https://github.com/lh3/bwa)
4) [blasr](https://github.com/PacificBiosciences/blasr) - (version 5.1 and above which uses '-' options)
5) [minimap2]( https://github.com/lh3/minimap2). (tested with version 2.11-r797)
6) [Racon](https://github.com/isovic/racon) - (version 0.5.0)
7) [Mash](https://github.com/marbl/Mash)
8) [MUMmer](http://mummer.sourceforge.net/)

All software are packaged as pre-build with OPERA-MS. Each binary is placed inside of the __utils__ folder.
If a pre-built software does not work on the user's machine, OPERA-MS will check if the tool is present in the user's PATH. However, the version of the software may be different than the one packaged. Alternatively, to specify a different directory for the dependency, a link to the software may be placed in the  __utils__ folder.

OPERA-MS is writen in Python, R and Perl, and makes use of the following Perl modules:
- [Switch](http://search.cpan.org/~chorny/Switch-2.17/Switch.pm)

- [File::Which](https://metacpan.org/pod/File::Which)

- [Statistics::Basic](http://search.cpan.org/~jettero/Statistics-Basic-1.6611/lib/Statistics/Basic.pod)

- [Statistics::R](https://metacpan.org/pod/Statistics::R)

# Contact information
For additional information, help and bug reports please send an email to one of the following: bertrandd@gis.a-star.edu.sg, nagarajann@gis.a-star.edu.sg