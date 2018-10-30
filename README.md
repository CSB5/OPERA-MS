# Introduction 
OPERA-MS is an hybrid metagenomic assembler which combines the advantages of short and long-read technologies providing high quality assemblies by addressing issues of low contiguity of short-read assemblies, and low base-pair quality of long-read assemblies. The method have been extensively tested on mock communities sequenced using different long-read sequencing technology (TruSeq Synthetic Long Read, PacBio and Oxford Nanopore) and provides high quality assemblies for read technologies with high base-pair error rate.

OPERA-MS uses the following strategy:
(1) Assembly of preliminary high base-pair quality contigs using a short-read metagenomic assembler (default: MEGAHIT),
(2) Short and long reads mapping to the assembly to respectively estimates the contig coverage and constructs an assembly graph representing the contig contiguity of the whole metagenome,
(3) Species level contig clustering based on coverage and contiguity information, and cluster augmentation using information provided by a database of reference genomes,
(4) Identification of species for which multiple strains have been assembled, and strain deconvolution using contig coverage information,
(5) Scaffolding and gap-filling of the strain and species clusters using [OPERA-LG](https://sourceforge.net/p/operasf/wiki/The%20OPERA%20wiki/).

OPERA-MS is able to assemble near complete genomes with as little as 9x long-read coverage on mock communities (compared to >30x required by other hybrid assembler). Evaluation on real gut microbiomes reveal that short-read contigs N50 can be improved up to 200x, and that OPERA-MS assembled twice as many draft quality assemblies with N50 >100kbp compared to other hybrid assembler. In comparison to other tools, OPERA-MS provides high quality and contiguous genomes when multiple strain of the same species are present.

# Installation

To install OPERA-MS run the following commands:

~~~~
git clone https://github.com/CSB5/OPERA-MS.git
cd /path/to/OPERA-MS
make
~~~~
If you encounter any problem during the installation, please see the [**Dependencies**](#dependencies) Section. 

A set of test files and a sample configuration file is provided to test out the OPERA-MS pipeline. To run OPERA-MS on the test data-set, simply use the following commands: 
~~~~
cd /path/to/OPERA-MS
perl OPERA-MS.pl sample_config.config 2> log.err
diff sample_files/contig_info.txt sample_output/contig_info.txt
~~~~
This will assemble a low diversity mock community in the folder **OPERA-MS/sample_output**. 

# Running OPERA-MS

OPERA-MS requires the specification of a configuration file that indicates the path to the input files and the options used for the assembly.
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

- **CONTIG_LEN_THR** : `default: 500` - Contig length threshold for clustering, contigs smaller than CONTIG_LEN_THR will be filtered out

- **CONTIG_EDGE_LEN** : `default: 80` - When calculating contig coverage, number of bases filtered out from each contig ends, to avoid biases due to lower mapping efficiency

- **CONTIG_WINDOW_LEN** : `default: 340` - The window length in which the coverage estimation is performed. We recommend using CONTIG_LEN_THR - 2 * CONTIG_EDGE_LEN as the value

- **KMER_SIZE** : `default: 60` - The kmer value used to produce the assembled contigs


### Outputs

Outputs can be found in the specified OUTPUT_DIR.
The file **contig.fasta** contains the assembled contigs, and **assembly.stats** provides overall assembly statistics (e.g. assembly size, N50, longest scaffold ...).
**contig_info.txt** provides a detailed overview of the assembled contigs according to the following features:
- **CONTIG_ID** : contig identifier. Single strain species contigs are named `opera_contig_X`. Contigs from multi-strain species are named `strainY_opera_contig_X`, where `Y` indicate the strain ID
- **LENGTH** : contig length
- **ARRIVAL_RATE** : median contig short-reads arrival rate
- **SPECIES** : putative species to which the assembled contig belong to
- **NB_STRAIN** : number of strains detected for the species
- **REFERENCE_GENOME** : path to the closest reference genome present in the OPERA-MS database

Finally, strain level scaffold assemblies can be found in the following files: **OUT_DIR/intermediate_files/strain_analysis/\*/\*/scaffoldSeq.fasta**.

# Dependencies

The only true dependency is `cpanm`, which is used to automatically install Perl modules. All other required software comes either pre-compiled with OPERA-MS or is build during the installation process. Binaries are placed inside the __utils__
folder:

1) [MEGAHIT](https://github.com/voutcn/megahit) - (tested with version 1.0.4-beta)
2) [samtools](https://github.com/samtools/samtools) - (version 0.1.19 or below)
3) [bwa](https://github.com/lh3/bwa) - (tested with version 0.7.10-r789)
4) [blasr](https://github.com/PacificBiosciences/blasr) - (version 5.1 and above which uses '-' options)
5) [minimap2]( https://github.com/lh3/minimap2) (tested with version 2.11-r797)
6) [Racon](https://github.com/isovic/racon) - (version 0.5.0)
7) [Mash](https://github.com/marbl/Mash) - (tested with version 1.1.1)
8) [MUMmer](http://mummer.sourceforge.net/) (tested with version 3.23)


If a pre-built software does not work on the user's machine, OPERA-MS will check if the tool is present in the user's PATH. However, the version of the software may be different than the one packaged. Alternatively, to specify a different directory for the dependency, a link to the software may be placed in the  **utils** folder.

OPERA-MS is written in C++, Python, R and Perl, and makes use of the following Perl modules (installed using [cpanm](https://metacpan.org/pod/distribution/App-cpanminus/bin/cpanm)):

- [Switch](http://search.cpan.org/~chorny/Switch-2.17/Switch.pm)
- [File::Which](https://metacpan.org/pod/File::Which)
- [Statistics::Basic](http://search.cpan.org/~jettero/Statistics-Basic-1.6611/lib/Statistics/Basic.pod)
- [Statistics::R](https://metacpan.org/pod/Statistics::R)

Note for Mac Users: the system default compiler (`clang`) will likely fail to compile `OPERA-LG`:
please install a recent GNU C++ compiler and point `make` to its path. For example, if you installed
GCC version 8.1.0 via Homebrew to /usr/local, then the following should work:


```
CXX=/usr/local/bin/g++-8 make

```


# Contact information
For additional information, help and bug reports please send an email to one of the following: 

- Denis Bertrand: <bertrandd@gis.a-star.edu.sg>
- Niranjan Nagarajan: <nagarajann@gis.a-star.edu.sg>
