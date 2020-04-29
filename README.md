[![Preprint at bioRxiv](https://img.shields.io/badge/preprint-bioRxiv-yellow.svg)](https://www.biorxiv.org/content/early/2018/10/30/456905)

# Introduction 
OPERA-MS is a hybrid metagenomic assembler which combines the advantages of short and long-read technologies to provide high quality assemblies, addressing issues of low contiguity for short-read only assemblies, and low base-pair quality for long-read only assemblies. OPERA-MS has been extensively tested on mock and real communities sequenced using different long-read technologies, including Oxford Nanopore, PacBio and Illumina Synthetic Long Read, and is particularly robust to noise in the read data.

OPERA-MS employs a staged assembly strategy that is designed to exploit even low coverage long read data to improve genome assembly. It begins by constructing a short-read metagenomic assembly (default: [MEGAHIT](https://github.com/voutcn/megahit)) that provides a good representation of the underlying sequence in the metagenome but may be fragmented. Long and short reads are then mapped to the assembly to identify connectivity between the contigs and to compute read coverage information. This serves as the basis for the core of the OPERA-MS algorithm which is to exploit coverage as well as connectivity information to accurately cluster contigs into genomes using a Bayesian model based approach. Another important advantage of OPERA-MS is that it can deconvolute strains in the metagenome, optionally using information from reference genomes to support this. After clustering, individual genomes are further scaffolded and gap-filled using the lightweight and robust scaffolder [OPERA-LG](https://sourceforge.net/p/operasf/wiki/The%20OPERA%20wiki/).

OPERA-MS can assemble near complete genomes from a metagenomic dataset with as little as 9x long-read coverage. Applied to human gut microbiome data it provides hundreds of high quality draft genomes, a majority of which have  N50 >100kbp. We observed the assembly of complete plasmids, many of which were novel and contain previously unseen resistance gene combinations. In addition, OPERA-MS can very accurately assemble genomes even in the presence of multiple strains of a species in a complex metagenome, allowing us to associate plasmids and host genomes using longitudinal data. For further details about these and other results using nanopore sequencing on stool samples from clinical studies see our [manuscript](https://www.nature.com/articles/s41587-019-0191-2) or [preprint](https://www.biorxiv.org/content/early/2018/10/30/456905). 

# Installation

To install OPERA-MS on a typical Linux/Unix system run the following commands:

```
git clone https://github.com/CSB5/OPERA-MS.git
cd OPERA-MS
make
perl OPERA-MS.pl check-dependency
```
If you encounter any problems during the installation, or if some third party software binaries are not functional on your system, please see the [**Dependencies**](#dependencies) section. 

A set of test files is provided to test out the OPERA-MS pipeline. To run OPERA-MS on the test dataset, simply use the following commands (please note that the test run requires 2 cores which is also the minimum): 
```
cd test_files
perl ../OPERA-MS.pl \
    --contig-file contigs.fasta \
    --short-read1 R1.fastq.gz \
    --short-read2 R2.fastq.gz \
    --long-read long_read.fastq \
    --no-ref-clustering \
    --out-dir RESULTS 2> log.err
```
This will assemble a low diversity mock community in the folder **RESULTS**.
Note that in the case of interruption during an OPERA-MS run, using the same command-line will re-start the execution after the last completed checkpoint.

# OPERA-MS genome database

To download the precomputed genome database required for the reference based clustering, simply run the following command:
```
perl OPERA-MS.pl install-db
```
The database contains a representative genome of the 23,000 bacteria species from [GTDB](https://gtdb.ecogenomic.org/) and requires 35Gb of free disc space.
Alternativelly a custom database can be generated the using the [**OPERA-MS-UTILS opera-db**](https://github.com/CSB5/OPERA-MS/wiki/Utilities#opera-ms-db) command.

# Usage

### Essential arguments

- **--short-read1** : path to the first read for Illumina paired-end read data (fasta/fastq/fasta.gz/fastq.gz)

- **--short-read2** : path to the second read for Illumina paired-end read data (fasta/fastq/fasta.gz/fastq.gz)

- **--long-read** : path to the long-read file obtained from either Oxford Nanopore, PacBio or Illumina Synthetic Long Read sequencing (fasta/fastq)

- **--out-dir** : directory where OPERA-MS results will be outputted

### Optional arguments 

 - **--genome-db** : path to a custom OPERA-MS genome database use during reference level clustering (defaut OPERA-MS-DB)

- **--no-ref-clustering**  : disable reference based clustering

- **--no-strain-clustering** : disable strain level clustering

- **--polishing** : enable short-read polishing (currently using [Pilon](https://github.com/broadinstitute/pilon/wiki)). The polished contigs can be found in contigs.polished.fasta

- **--long-read-mapper** : software used for long-read mapping i.e. blasr (default) or minimap2

- **--short-read-assembler** : software used for short read assembly i.e. megahit (default) or spades

- **--no-gap-filling** : disable gap-filling stage

- **--kmer-size** : kmer value (default=60) used to assemble contigs

- **--contig-len-thr** : contig length threshold (default=500) for clustering; contigs smaller than contig-len-thr will be filtered out

- **--contig-edge-len** : during contig coverage calculation, number of bases filtered out from each contig end (default=80), to avoid biases due to lower mapping efficiency

- **--contig-window-len** : window length (default=340) in which the coverage estimation is performed. We recommend using contig_len_thr - 2 * contig_edge_len as the value

- **--contig-file** : path to the contig file, if the short-reads have been assembled previously

- **--num-processors** : number of processors to use (default=2; note that 2 is the minimum)

Alternatively, OPERA-MS parameters can be set using a [configuration file](https://github.com/CSB5/OPERA-MS/wiki/Configuration-file).

### Output

The following output files can be found in the specified output directory i.e. **RESULTS**.
The file **contigs.fasta** (and **contigs.polished.fasta** if the assembly has been polished) contains the assembled contigs, **assembly.stats** provides overall assembly statistics (e.g. assembly size, N50, longest contig etc.), and
[**contig_info.txt**](https://github.com/CSB5/OPERA-MS/wiki/Contig-info-file-description) provides detailed overview of the assembled contigs.

Finally, OPERA-MS strain-level clusters (one fasta file per strain) can be found in the directory **RESULTS/opera_ms_clusters/all** and [**cluster_info.txt**](https://github.com/CSB5/OPERA-MS/wiki/Cluster-info-file-description) provides a detailed overview of assembly statistics for these clusters. Note that these clusters are constructed for producing high-quality assemblies and are therefore conservative. They can be binned further using approaches such as [MaxBin2](https://sourceforge.net/projects/maxbin2/) or [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/).

### OPERA-MS-UTILS

Scripts to post-process the assemblies are now availaible using the [**OPERA-MS-UTILS**](https://github.com/CSB5/OPERA-MS/wiki/Utilities) command. We are now providing stream line analysis tools to compute the concordance between the short and long read, the contig binning, bin assessment and annotation of novel species (as described in the paper). A complete desctiption of those tools can be found in the [the OPERA-MS-UTILS wiki section](https://github.com/CSB5/OPERA-MS/wiki/Utilities).

This is a work in progress and additional tools will be available as part of the next release. Please contact us if you would like to add your favorite metagenomic analysis tool.

# Resource Requirements

OPERA-MS's runtime depends on the complexity of the metagenome and the amount of short/long-read data available.
We typically run OPERA-MS with default parameters using 16 threads on an Intel Xeon platinum server with SSD hard drive. With this hardware specification, we obtain the following running time and memory usage characteristics.

| Dataset  | Short-read data (Gbp) | Long-read data (Gbp) | Running time (hours)  | Peak RAM usage (GB) |
|:---:                 |:---: |:---:   |:---:   |:---: |
| CAMI multi-strain mock community (low complexity)  | 3.9  | 2    | 1.4  | 5.5| 
| Human gut microbiome (medium complexity) | 24.4  | 1.6  | 2.7  | 10.2| 
| CAMI environmental mock community (high complexity)  | 9.9  | 4.8  | 4.5    | 12.8|

OPERA-MS is designed to work with deep short-read sequencing, but can work with lower coverage in terms of long-read sequencing. In practice, short-read coverage >15x is recommended, while OPERA-MS can use long-read coverage as low as 9x to boost assembly contiguity. Based on this, we recommend at least 9Gbp of short-read data and 3Gbp of long-read data to allow for assembly of bacterial genomes at 1% relative abundance in the metagenome.

# Dependencies

The only true dependency is `cpanm`, which is used to automatically install Perl modules. All other required programs come either pre-compiled with OPERA-MS or are built during the installation process. Binaries are placed inside the __utils__ folder:

1) [MEGAHIT](https://github.com/voutcn/megahit) - (tested with version 1.0.4-beta)
1) [SPAdes](https://github.com/ablab/spades) - (tested with version 3.13.0)
2) [samtools](https://github.com/samtools/samtools) - (version 0.1.19 or below)
3) [bwa](https://github.com/lh3/bwa) - (tested with version 0.7.10-r789)
4) [blasr](https://github.com/PacificBiosciences/blasr) - (version 5.1 and above which uses '-' options)
5) [minimap2]( https://github.com/lh3/minimap2) (tested with version 2.11-r797)
6) [Racon](https://github.com/isovic/racon) - (version 0.5.0)
7) [Mash](https://github.com/marbl/Mash) - (tested with version 1.1.1)
8) [MUMmer](http://mummer.sourceforge.net/) - (tested with version 3.23)
9) [Pilon](https://github.com/broadinstitute/pilon/wiki) - (tested with version 1.22)

If a pre-built program does not work on the user's machine, OPERA-MS will check if the program is present in the user's PATH. However, the version of the program may be different than the one packaged. Alternatively, to specify a different directory for the dependency, a link to the program may be placed in the __utils__ folder.

OPERA-MS and its dependencies required C++, Java, Python, R and Perl, and use the following Perl modules (installed using [cpanm](https://metacpan.org/pod/distribution/App-cpanminus/bin/cpanm)):

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
If the perl libraries cannot be installed under root, the following line should be added to `.bashrc`:
```
export PERL5LIB="/home/$USER/perl5/lib/perl5${PERL5LIB:+:${PERL5LIB}}";
```

# Docker
A simple Dockerfile is provided in the root of the repository. To build the image:
```
[comment]:docker build -t operams .
```
The generic command to run a OPERA-MS docker container after building:
```
[comment]:docker run \
    -v /host/path/to/indata/:/indata/ \
    -v /host/path/to/outdata/:/outdata/ \
     operams
    --short-read1 /indata/R1.fastq.gz \
    --short-read2 /indata/R2.fastq.gz \
    --long-read /indata/long_read.fastq \
    --out-dir /outdata 
```
To process data with the dockerized OPERA-MS, directories for in- and outdata should be mounted into the container. An example is shown below for running the test dataset. In the below example the repo was cloned to /home/myuser/git/OPERA-MS/). The repo is needed only for the `sample_files` directory and the `sample_config.config` file. If Docker is running in a VM, as is the case for Windows or OSX, but also when deployed on a cloud platform such as AWS or Azure, a minimum of 2 available cores is required.  

```
docker run \
    -v OPERA-MS/test_files/:/sample_files \
    -v OPERA-MS/test_files/RESULTS/:/sample_out \
     operams --contig-file /sample_files/contigs.fasta \
    --short-read1 /sample_files/R1.fastq.gz \
    --short-read2 /sample_files/R2.fastq.gz \
    --long-read /sample_files/long_read.fastq \
    --out-dir /sample_out
```


# Contact information
For additional information, help and bug reports please send an email to: 

- Denis Bertrand <bertrandd@gis.a-star.edu.sg>
- Chengxuan Tong <Tong_Chengxuan@gis.a-star.edu.sg>
- Lorenz Gerber (docker related topics): <lorenzottogerber@gmail.com>


