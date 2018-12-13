[![Preprint at bioRxiv](https://img.shields.io/badge/preprint-bioRxiv-yellow.svg)](https://www.biorxiv.org/content/early/2018/10/30/456905)

# Introduction 
OPERA-MS is a hybrid metagenomic assembler which combines the advantages of short and long-read technologies to provide high quality assemblies, addressing issues of low contiguity for short-read only assemblies, and low base-pair quality for long-read only assemblies. OPERA-MS has been extensively tested on mock and real communities sequenced using different long-read technologies, including Oxford Nanopore, PacBio and Illumina Synthetic Long Read, and is particularly robust to noise in the read data.

OPERA-MS employs a staged assembly strategy that is designed to exploit even low coverage long read data to improve genome assembly. It begins by constructing a short-read metagenomic assembly (default: [MEGAHIT](https://github.com/voutcn/megahit)) that provides a good representation of the underlying sequence in the metagenome but may be fragmented. Long and short reads are then mapped to the assembly to identify connectivity between the contigs and to compute read coverage information. This serves as the basis for the core of the OPERA-MS algorithm which is to exploit coverage as well as connectivity information to accurately cluster contigs into genomes using a Bayesian model based approach. Another important feature of OPERA-MS is that it can use information from reference genomes to deconvolute strains in the metagenome. After clustering, individual genomes are further scaffolded and gap-filled using the lightweight and robust scaffolder [OPERA-LG](https://sourceforge.net/p/operasf/wiki/The%20OPERA%20wiki/).

OPERA-MS can assemble near complete genomes from a metagenomic dataset with as little as 9x long-read coverage. Applied to human gut microbiome data it provides hundreds of high quality draft genomes, a majority of which have  N50 >100kbp. We observed the assembly of complete plasmids, many of which were novel and contain previously unseen resistance gene combinations. In addition, OPERA-MS can very accurately assemble genomes even in the presence of multiple strains of a species in a complex metagenome, allowing us to associate plasmids and host genomes using longitudinal data. For further details about these and other results using nanopore sequencing on stool samples from clinical studies see our [bioRxiv paper](https://www.biorxiv.org/content/early/2018/10/30/456905). 

# Installation

To install OPERA-MS on a typical Linux/Unix system run the following commands:

```
git clone https://github.com/CSB5/OPERA-MS.git
cd /path/to/OPERA-MS
make
perl OPERA-MS.pl sample_config.config TEST_INSTALLATION 2> log.err
```
If you encounter any problems during the installation, or if some third party sofware binaries are not functional on your system, please see the [**Dependencies**](#dependencies) section. 

A set of test files and a sample configuration file is provided to test out the OPERA-MS pipeline. To run OPERA-MS on the test data-set, simply use the following commands (please note that the test script runs with 2 cores which is also the minimum.): 
```
cd /path/to/OPERA-MS
perl OPERA-MS.pl sample_config.config 2> log.err
diff sample_output/assembly.stats sample_files/assembly.stats
```
This will assemble a low diversity mock community in the folder **OPERA-MS/sample_output**.

# Running OPERA-MS

OPERA-MS requires the specification of a configuration file that indicates the path to the input files and the options used for the assembly.
The configuration file is formatted as follows:

~~~~
#One space between OPTION and VALUE
<OPTION1> <VALUE1> 
<OPTION2> <VALUE2>
...
<OPTION2> <VALUE3>
~~~~

### Essential parameters

- **OUTPUT_DIR** : `path/to/results` - Directory where OPERA-MS results will be outputted

- **LONG_READ** : `path/to/long-read.fq` - Path to the long-read fastq file obtained from either Oxford Nanopore, PacBio or Illumina Synthetic Long Read sequencing

- **ILLUMINA_READ_1** : `path/to/illum_read1.fq.gz` - Path to the first read for Illumina paired-end read data

- **ILLUMINA_READ_2** : `path/to/illum_read2.fq.gz` - Path to the second read for Illumina paired-end read data

### Optional parameters 

- **CONTIGS_FILE** : `path/to/contigs.fa` - Path to the contig file, if the short-reads have been assembled previously

- **NUM_PROCESSOR** : `default : 2` - Number of processors to use (note that 2 is the minimum)

- **LONG_READ_MAPPER** `default: blasr` - Software used for long-read mapping i.e. blasr or minimap2

- **STRAIN_CLUSTERING** : `default: YES` - Whether strain level clustering should be performed (YES) or skipped (NO)

- **CONTIG_LEN_THR** : `default: 500` - Contig length threshold for clustering; contigs smaller than CONTIG_LEN_THR will be filtered out

- **CONTIG_EDGE_LEN** : `default: 80` - When calculating contig coverage, number of bases filtered out from each contig end, to avoid biases due to lower mapping efficiency

- **CONTIG_WINDOW_LEN** : `default: 340` - Window length in which the coverage estimation is performed. We recommend using CONTIG_LEN_THR - 2 * CONTIG_EDGE_LEN as the value

- **KMER_SIZE** : `default: 60` - Kmer value used to assemble contigs


### Output

The following output files can be found in the specified output directory i.e. OUTPUT_DIR.
The file **contig.fasta** contains the assembled contigs, and **assembly.stats** provides overall assembly statistics (e.g. assembly size, N50, longest contig etc.).
**contig_info.txt** provides a detailed overview of the assembled contigs with the following information:
- **CONTIG_ID** : contig identifier, typically `opera_contig_X`. Contigs from species where OPERA-MS detects multiple genomes are named `strainY_opera_contig_X` to record this information, where `Y` indicates the strain ID
- **LENGTH** : contig length
- **ARRIVAL_RATE** : median short-read arrival rate for the contig
- **SPECIES** : closest reference species in the OPERA-MS database for the assembled contig
- **NB_STRAIN** : number of strain genomes in the metagenome that were detected by OPERA-MS for this species
- **REFERENCE_GENOME** : path to the closest reference genome present in the OPERA-MS database

Finally, strain level contig assemblies can be found in the following directory: **OUT_DIR/strain_assembly/**. For each species, strains are ordered from the most abundant to the less abundant.

# Dependencies

The only true dependency is `cpanm`, which is used to automatically install Perl modules. All other required software comes either pre-compiled with OPERA-MS or is built during the installation process. Binaries are placed inside the __utils__
folder:

1) [MEGAHIT](https://github.com/voutcn/megahit) - (tested with version 1.0.4-beta)
2) [samtools](https://github.com/samtools/samtools) - (version 0.1.19 or below)
3) [bwa](https://github.com/lh3/bwa) - (tested with version 0.7.10-r789)
4) [blasr](https://github.com/PacificBiosciences/blasr) - (version 5.1 and above which uses '-' options)
5) [minimap2]( https://github.com/lh3/minimap2) (tested with version 2.11-r797)
6) [Racon](https://github.com/isovic/racon) - (version 0.5.0)
7) [Mash](https://github.com/marbl/Mash) - (tested with version 1.1.1)
8) [MUMmer](http://mummer.sourceforge.net/) (tested with version 3.23)


If a pre-built software does not work on the user's machine, OPERA-MS will check if the tool is present in the user's PATH. However, the version of the software may be different than the one packaged. Alternatively, to specify a different directory for the dependency, a link to the software may be placed in the **utils** folder.

OPERA-MS is written in C++, Python, R and Perl, and makes use of the following Perl modules (installed using [cpanm](https://metacpan.org/pod/distribution/App-cpanminus/bin/cpanm)):

- [Switch](http://search.cpan.org/~chorny/Switch-2.17/Switch.pm)
- [File::Which](https://metacpan.org/pod/File::Which)
- [File::Spec::Functions](https://perldoc.perl.org/File/Spec/Functions.html)
- [Statistics::Basic](http://search.cpan.org/~jettero/Statistics-Basic-1.6611/lib/Statistics/Basic.pod)
- [Statistics::R](https://metacpan.org/pod/Statistics::R)
- [Getopt::Long](http://perldoc.perl.org/Getopt/Long.html)

Once cpanm is installed, simply run the following command to install all the perl modules:
```
perl utils/install_perl_module.pl
```

# Docker
A simple Dockerfile is provided in the root of the repository. To build the image:
```
docker build -t operams .
```
The generic command to run a OPERA-MS docker container after building:
```
docker run \
    -v /host/path/to/indata/:/indata/ \
    -v /host/path/to/outdata/:/outdata/ \
    operams config.file
```
To process data with the dockerized OPERA-MS, directories for in- and outdata should be mounted into the container. An example is shown below for running the test dataset. In the below example the repo was cloned to /home/myuser/git/OPERA-MS/). The repo is needed only for the `sample_files` directory and the `sample_config.config` file. If Docker is running in a VM, as is the case for Windows or OSX, but also when deployed on a cloud platform such as AWS or Azure, a minimum of 2 available cores is required.  

```
docker run \ 
    -v /home/myuser/git/OPERA-MS/sample_files:/sample_files \
    -v /home/myuser/git/OPERA-MS/sample_output:/sample_output \
     operams sample_config.config
```

# Contact information
For additional information, help and bug reports please send an email to: 

- Denis Bertrand: <bertrandd@gis.a-star.edu.sg>
- Lorenz Gerber (docker related topics): <lorenzottogerber@gmail.com>
- Niranjan Nagarajan: <nagarajann@gis.a-star.edu.sg>
