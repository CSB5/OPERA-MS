Program name: MaxBin
Version: 2.2.4
Developer: Yu-Wei Wu (ywwei@lbl.gov)
Affiliation: Joint BioEnergy Institute, Lawrence Berkeley National Laboratory



===================== Quick installation guideline ===========================================
=== 1. Download MaxBin and unzip it
=== 2. Enter src directory under MaxBin and "make" it.
=== 3. Run "./autobuild_auxiliary" at MaxBin directory to download, compile,
===    and setup the auxiliary software packages
=== 4. MaxBin should be ready to go.
===================== End of Quick installation guideline ====================================



=== What MaxBin can do? ===
MaxBin is a software that is capable of clustering metagenomic contigs into different bins, each consists of contigs from one species. MaxBin uses the nucleotide composition information and contig abundance information to do achieve binning through an Expectation-Maximization algorithm.


=== Support platform ===
MaxBin was developed and maintained on Linux platform. It has been tested on docker images of Ubuntu, Redhat, SUSE, Debian, Fedora, and Mageia. MaxBin core program can now be compiled on Mac platform; however the association of auxiliary software needs to be performed manually on Mac. Please refer to the manual way to install auxiliary software section as described below.


=== Prerequisite ===
1. g++ (c++11 or higher), make
2. Perl5 with CPAN, LWP::Simple, and FindBin modules
   (./autobuild_auxiliary will use CPAN to check and install the other two Perl modules)
3. (for automatic downloading and installing auxiliary software) curl, unzip, which


=== Installation ===
The installation of MaxBin is two-fold. MaxBin needs some auxiliary software packages to run correctly. The installation steps of MaxBin is as follows.

1. Enter MaxBin directory abd run "make" under that directory. This will build the MaxBin executable. The commands are:
   $ cd src
   $ make
   $ cd ..

2. Download and compile auxiliary software packages. There are two ways:

   - The easiest way is execute the script "./autobuild_auxiliary." It will attempt to download all auxiliary software packages from mirror sites, compile them, and set the runtime environment. It is highly recommended that you try this option for less hassle.

   - The manual way. Auxiliary packages include Bowtie2, FragGeneScan, IDBA-UD, and Hmmer3. Please follow the instructions in the software packages to compile and install them.

	* Bowtie2 (tested version: bowtie2-2.2.3)
	  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	* FragGeneScan (tested version: FragGeneScan 1.30)
	  http://omics.informatics.indiana.edu/FragGeneScan/
	* Hmmer3 (tested version: HMMER 3.1b1 64 bit)
	  http://hmmer.janelia.org/
	* IDBA-UD (tested version: IDBA-UD 1.1.1)
	  http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/

	* You also need to Specify the software paths. There are two different ways.

	  a. Export the paths of software executables in system path. For example, in ubuntu linux you can add a line in ~/.bashrc file:
             export PATH=/usr/local/bin/FragGeneScan_1.30:/usr/local/bin/hmmer3:/usr/local/bin/bowtie2:/local/bin/idba-1.1.1:$PATH
             (This is just an example; you can specify your own software locations freely)

	  b. Alternatively, you can specify the software paths in the "setting" file under MaxBin directory. The file is in the follow format:

             [FragGeneScan] /local/bin/FragGeneScan1.30.
             [Bowtie2] /local/bin/bowtie2-2.0.0-beta7
             [HMMER3] /local/bin/hmmer-3.0-linux-intel-x86_64/binaries
             [IDBA-UD] /local/bin/idba-1.1.1/bin

             (This is just an example. Please specify your own software locations.)




---Plotting the number of each marker gene in each contig---
MaxBin provides the functionality to plot the single copy marker genes in each bin using R (with gplots package installed). Please install R beforehand and make sure that R and Rscript can be executed by type the two commands in the command line.


=== Run MaxBin ===
To run MaxBin, please type in "perl run_MaxBin.pl" or just "run_MaxBin.pl". You will see options. Here are the options:

(required) -contig (contig filename)
(required) -out (output file header)

--- at least one of the following parameters is needed
(semi-required) -abund (contig abundance files. To be explained in Abundance session below.)
(semi-required) -reads (reads file in fasta or fastq format. To be explained in Abundance session below.)
(semi-required) -abund_list (a list file of all contig abundance files.)
(semi-required) -reads_list (a list file of all reads file.)
---
(optional) -thread (number of threads; default 1)
(optional) -prob_threshold (minimum probability for EM algorithm; default 0.8)
(optional) -plotmarker (specify this option if you want to plot the markers in each contig. Installing R is a must for this option to work.)
(optional) -verbose (as is. Warning: output log will be LOOOONG.)
(optional) -markerset (choose between 107 marker genes by default or 40 marker genes. see Marker Gene Note for more information.)

Example commands:
run_MaxBin.pl -contig mycontig -abund myabund -out myout
run_MaxBin.pl -contig mycontig -reads myreads -reads2 my2ndreads -reads3 my3rdreads -out myout
run_MaxBin.pl -contig mycontig -abund myabund -abund2 my2ndabund -abund3 my3rdabund -reads myreads -reads2 my2ndreads -out myout -thread 4
run_MaxBin.pl -contig mycontig -abund_list abund_list_file -reads_list reads_list_file -out myout -prob_threshold 0.5 -markerset 40

(Warning: Please do NOT specify abundance and reads that BELONG TO THE SAME METAGENOME together. MaxBin will treat them as two different information and thus will count this metagenome TWICE!)


=== Reassembly Note ===
Reassembly option is still highly experimental. To use this function, you need to feed MaxBin "interleaved paired-end" fastq or fasta file if you were to use this option. An example of interleaved paired-end fasta file is as follows:

>reads1.1
AAAAA
>reads1.2
CCCCC
>reads2.1
TTTTT
>reads2.2
GGGGG
>reads3.1
AAAAA
>reads3.2
TTTTT



=== Marker Gene Note ===
By default MaxBin will look for 107 marker genes present in >95% of bacteria. Alternatively you can also choose 40 marker gene sets that are universal among bacteria and archaea (Wu et al., PLoS ONE 2013). This option may be better suited for environment dominated by archaea; however it tend to split genomes into more bins. You can choose between different marker gene sets and see which one works better.


=== Abundance ===
The contig abundance information can be provided in two ways: user can choose to provide the abundance file or MaxBin will use Bowtie2 to map the sequencing reads against contigs and generate the abundance information.

---if you have the abundance information---
Please make sure that your abundance information is provided in the following format (\t stands for a tab delimiter):

(contig header)\t(abundance)

For example, assume I have three contigs named A0001, A0002, and A0003, then my abundance file will look like

A0001	30.89
A0002	20.02
A0003	78.93

---if you don't have abundance information---
Don't worry, MaxBin will generate this information for you from sequencing reads. Please specify the reads file (in fasta format) in -reads.


---Providing more than one reads or abundance files
The reads and abundance files can be provided via a file consisting of all reads or abundance files. For example, you have 3 abundance files and 5 reads files, and you want to provide them all. You can create two files "abund_list" and "reads_list", which are

$ cat abund_list
(first abundance file) 
(second abundance file) 
(third abundance file) 

$ cat reads_list
(first reads file)
(second reads file)
(third reads file)
(fourth reads file)
(fifth reads file)

Then you can feed all information into MaxBin using -abund_list and -reads_list options. Handy for a large number of reads or abundance files.


===MaxBin Output===
Assume your output file header is (out). MaxBin will generate information using this file header as follows.

(out).0XX.fasta -- the XX bin. XX are numbers, e.g. out.001.fasta
(out).summary -- a summary file describing which contigs are being classified into which bin.
(out).log -- a log file recording the core steps of MaxBin algorithm
(out).marker -- marker gene presence numbers for each bin. This table is ready to be plotted by R or other 3rd-party software.
(out).marker.pdf -- visualization of the marker gene presence numbers using R. Will only appear if -plotmarker is specified.
(out).noclass -- this file stores all sequences that pass the minimum length threshold but are not classified successfully.
(out).tooshort -- this file stores all sequences that do not meet the minimum length threshold.
(out).marker_of_each_gene.tar.gz -- this tarball file stores all markers predicted from the individual bins. Use "tar -zxvf (out).marker_of_each_gene.tar.gz" to extract the markers [(out).0XX.marker.fasta].


===Bug or problems===
Encounter bugs, problems, or have any suggestions? Please contact Yu-Wei Wu (ywwei@lbl.gov).


