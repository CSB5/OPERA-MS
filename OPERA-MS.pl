#!/usr/bin/perl

#use strict;
use warnings;
use Cwd;
use File::Spec;
use Switch;
use Getopt::Std;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use File::Which;

my $STAGE_TO_RUN = "ALL";

#Grab the executing directory; turn all paths into absolute paths.
my $opera_ms_dir = dirname(rel2abs($0)) . "/";

my $main_dir = getcwd;
$main_dir .= "\/";

my $racon_dir;
my $opera_version = "OPERA-LG";
my $runOperaMS_config_name = "runOperaMS.config"; #Generated config file, name can be whatever.
my $output_dir;
my $opera_ms_config_file;
my $long_read_file;
my $illum_read1;
my $illum_read2;
my $contigs_file;
my $kmer_size = 60;
my $OPERAMS_config_name; 
my $opera_ms_cf;
my $config_option;
my $samtools_dir;
my $minimap2_dir;
my $megahit_dir;
my $mash_dir;
my $mummer_dir;
my $blasr_dir;
my $strain_clustering = "YES";
my $short_read_tool_dir;
my $short_read_maptool = "bwa";
my $num_processor = 1;
my $help_line = "To configure OPERA-MS, please look at the example config file inside the OPERA-MS folder.\nUsage: \n\npath/OPERA-MS/OPERA-MS.pl <config_file>\n\nNote that the config file must be inside the path/OPERA-MS directory.\n";
my $incorrect_arguments="Please input the correct arguments. Refer to the documentation for more information or use:\n\n path/OPERA-MS/OPERA-MS.pl -h. \n\n";
my $run_following = 0;
my $contig_edge_len = 80;
my $contig_window_len = 340;
my $contig_len_thr = 500;

if ( @ARGV == 0 ){
    die $incorrect_arguments; 
}

#read from config file
if ( @ARGV >= 1){
    $OPERAMS_config_name = $ARGV[0];
    $STAGE_TO_RUN = $ARGV[1] if(@ARGV == 2);

    if(index($STAGE_TO_RUN, "+") != -1){
	chop $STAGE_TO_RUN;
	$run_following = 1;
    }
    
    die $help_line if($ARGV[0] eq "-h");
    print STDERR "\nReading config file: ".$OPERAMS_config_name."\n";
    die"Config file does not exist. Exiting.\n" if(!(-e $OPERAMS_config_name)); 

    open($opera_ms_cf, "<", $OPERAMS_config_name); 

    while(<$opera_ms_cf>) {
        next if (/^#/);  #skip comments
		 chomp($_);  
		 my @split_line = split('\s+', $_);
		 if (@split_line != 0) {
		     $config_option = $split_line[0];
		     switch ($config_option) {
                
                #Empty cases must be checked so the same config file format
                #can be used for the runOperaMS generated config file. Empty cases indicate
                #that runOperaMS.pl uses those options. 
                  
			 case "CONTIGS_FILE" {
			     $contigs_file = $split_line[1];
			 }

			 case "MAPPING_FILES"{
			 }

			 case "LIB"{
			 }

			 case "EDGE_BUNDLESIZE_THRESHOLD" {
			 }

			 case "OUTPUT_DIR" {
			     $output_dir = $split_line[1];
			 }
			 
			 case "STRAIN_CLUSTERING"{
			     $strain_clustering = $split_line[1];
			 }
			 
			 case "SAMTOOLS_DIR" {
			     $samtools_dir = $split_line[1];
			 }

			 case "BLASR_DIR"{
			     $blasr_dir = $split_line[1];
			 }

			 case "SHORT_READ_TOOL"{
			     $short_read_maptool = $split_line[1];
			 }

			 case "BWA_DIR"{
			     $short_read_tool_dir = $split_line[1];
			 }

			 case "MINIMAP2_DIR" {
			     $minimap2_dir = $split_line[1];
			 }
 
			 case "MEGAHIT_DIR" {
			     $megahit_path = $split_line[1];
			 }
			 
			 case "MASH_DIR" {
			     $mash_dir = $split_line[1];
			 }
			 
			 case "MUMMER_DIR" {
			     $mummer_dir = $split_line[1];
			 }

			 
			 case "SIGMA_CONTIGS_FILE" {
             
			 }


			 case "RACON_DIR"{
			     $racon_dir = $split_line[1];
			 }

			 case "KMER_SIZE" {
			     $kmer_size = $split_line[1];
			     #if (!defined $kmer_size) {
			     #die "KMER_SIZE not provided.\n";
			     #}
			 }

			 #Option uses for runOperaMS.pl
			 case "CONTIGS_FILE_TYPE" {
			 }

			 case "CONTIG_LEN_THR" {
			     $contig_len_thr = $split_line[1];
			 }

			 case "CONTIG_EDGE_LEN" {
			     $contig_edge_len = $split_line[1];
			 }

			 case "CONTIG_WINDOW_LEN" {
			     $contig_window_len = $split_line[1];
			 }

			 case "PDIST_TYPE" {
			 }

			 case "SKIP_OPERA" {
               
			 }

			 case "LONG_READ_MAPPER"{
			     $mapper = $split_line[1];
			     if($mapper ne "blasr" && $mapper ne "minimap2"){
				 die "Unkown LONG_READ_MAPPER $mapper, please use  blasr/minimap2.\n"
			     }
			 }
			 
			 case "LONG_READ"{
			     $long_read_file = $split_line[1];
			     if ($long_read_file =~ /\.gz$/){
				 die "$long_read_file appears to be in gzipped format. Blasr does not accept gzipped files, please unzip.\n";
			     }
			 }

			 case "ILLUMINA_READ_1"{
			     $illum_read1 = $split_line[1];
			 }

			 case "ILLUMINA_READ_2"{
			     $illum_read2 = $split_line[1];
			 }

			 case "NUM_PROCESSOR"{
			     $num_processor = $split_line[1];
			 }

			 else {
			     die "Config option: ".$config_option." unknown, please check the config file. \nExiting. \n";
			 }
		     }

		 }

    }
     
}

else{
    die $incorrect_arguments;
}

#Download the genome reference data base for the first utilasation
if(! -e "$opera_ms_dir/database_updated"){
    print STDOUT " *** First utilization: set-up of reference genome databases. Please wait ...\n";
    #
    print STDOUT " *** (1/3) Download of genomeDB_Sketch.msh\n";
    run_exe("wget --no-check-certificate -O $opera_ms_dir/genomeDB_Sketch.msh https://www.dropbox.com/s/kot086vh26nds6j/genomeDB_Sketch.msh?dl=0");
    #
    print STDOUT " *** (2/3) Download of complete_ref_genomes.tar.gz\n";
    run_exe("wget --no-check-certificate -O $opera_ms_dir/complete_ref_genomes.tar.gz https://www.dropbox.com/s/wcxvy2u0yhp3pw0/complete_ref_genomes.tar.gz?dl=0");
    #
    print STDOUT " *** (3/3) Extraction of complete_ref_genomes.tar.gz\n";
    run_exe("tar -xvzf $opera_ms_dir/complete_ref_genomes.tar.gz --directory $opera_ms_dir");
    #
    run_exe("rm $opera_ms_dir/complete_ref_genomes.tar.gz");
    #
    run_exe("touch $opera_ms_dir/database_updated");
}

my $inter = "intermediate_files"; 
$output_dir = $main_dir . $output_dir."/" if (substr($output_dir, 0, 1) ne "/");
run_exe("mkdir -p $output_dir/$inter") if(! -d "$output_dir/$inter");

#Make paths absolute if they are not already.
$long_read_file = $main_dir . $long_read_file if(substr($long_read_file, 0, 1) ne "/");
$illum_read1 = $main_dir . $illum_read1 if(substr($illum_read1, 0, 1) ne "/");
$illum_read2= $main_dir . $illum_read2 if(substr($illum_read2, 0, 1) ne "/");
if(defined $contigs_file){
    $contigs_file = $main_dir . $contigs_file if(substr($contigs_file, 0, 1) ne "/");
}
$opera_ms_config_file = $output_dir."/".$runOperaMS_config_name;

#
#Check if dependencies are found. Otherwise, die.
#
#
#Short-read bundler
if (! -e "$opera_ms_dir/bin/bundler" or
    ! -e "$opera_ms_dir/OPERA-LG/bin/OPERA-LG"){
    die "bundler and OPERA-LG not found. Please make install before running OPERA-MS.\n";
}

#samtools
if(!defined($samtools_dir)){
    $samtools_dir = "$opera_ms_dir/utils/";
    my $exe_check = "${samtools_dir}samtools 2>&1 | grep Version";
    run_exe($exe_check);
    
    if ($?){
        print STDERR "\nsamtools found in $opera_ms_dir/utils is not functional. Checking path for samtools. \n";
        $samtools_dir = "";
        my $valid_path = which('samtools');
        die "samtools not found in path. Exiting.\n" if (!$valid_path);
	$valid_path = `which samtools`;
	@tmp = split(/\//, $valid_path);
	$samtools_dir = join("\/", @tmp[0..(@tmp-2)])."/";
    }
}
else{
    $samtools_dir = $main_dir . $samtools_dir . "/" if(substr($samtools_dir, 0, 1) ne "/");
    
    if (! -e $samtools_dir . "/samtools") {
        die "Samtools not found at: ".$samtools_dir."\n";
    }
    
}

#blasr
if(!defined($blasr_dir)){
    $blasr_dir = "$opera_ms_dir/utils/";
    my $exe_check = "${blasr_dir}blasr -h";
    run_exe($exe_check);
    if ($?){
	$blasr_dir = "";
	my $valid_path = which("blasr");
	die "blasr not found in path. Exiting.\n" if (!$valid_path);
	$valid_path = `which blasr`;
	@tmp = split(/\//, $valid_path);
	$blasr_dir = join("\/", @tmp[0..(@tmp-2)]);
    }
}
else{
    $blasr_dir = $main_dir . $blasr_dir . "/" if (substr($blasr_dir, 0, 1) ne "/");
    if (! -e $blasr_dir . "/blasr") {
        die "blasr not found at: ".$blasr_dir."\n";
    }
}

#bwa mapping
if(!defined($short_read_tool_dir)){
    $short_read_tool_dir = "$opera_ms_dir/utils/";
    my $exe_check = "${short_read_tool_dir}bwa 2>&1 | grep Version";
    run_exe($exe_check);

    if($?){
        print STDERR "\nbwa found in $opera_ms_dir/utils is not functional. Checking path for bwa. \n";
        $short_read_tool_dir= "";
        my $valid_path = which($short_read_maptool); 
        die "$short_read_maptool not found in path. Exiting.\n" if (!$valid_path);
	$valid_path = `which bwa`;
	@tmp = split(/\//, $valid_path);
	$short_read_tool_dir = join("\/", @tmp[0..(@tmp-2)])."/";
    }
}
else{
    $short_read_tool_dir= $main_dir . $short_read_tool_dir. "/" if (substr($short_read_tool_dir, 0, 1) ne "/");
    if (! -e $short_read_tool_dir . "/$short_read_maptool") {
        die "$short_read_maptool not found at : ".$short_read_tool_dir."\n";
    }
}

#minimap2
if(!defined($minimap2_dir)){
    $minimap2_dir = "$opera_ms_dir/utils/";
    my $exe_check = "${minimap2_dir}minimap2 --version";
    run_exe($exe_check);
    if ($?){
        print STDERR "\nminimap2 found in $opera_ms_dir/utils is not functional. Checking path for minimap2. \n";
        $minimap2_dir = "";
        my $valid_path = which('minimap2');
        die "minimap2 not found in path. Exiting.\n" if (!$valid_path);
	$valid_path = `which minimap2`;
	@tmp = split(/\//, $valid_path);
	$minimap2_dir = join("\/", @tmp[0..(@tmp-2)])."/";
    }
}
else{
    $minimap2_dir = $main_dir . $minimap2_dir . "/" if(substr($minimap2_dir, 0, 1) ne "/");
    if (! -e $minimap2_dir . "/minimap2") {
        die "minimap2 not found at: ".$minimap2_dir."\n";
    }
}

#mash
if(!defined($mash_dir)){
    $mash_dir = "$opera_ms_dir/utils/";
    my $exe_check = "${mash_dir}mash --version";
    run_exe($exe_check);
    if ($?){
        print STDERR "\nmash found in $opera_ms_dir/utils is not functional. Checking path for mash. \n";
        $mash_dir = "";
        my $valid_path = which('mash');
        die "mash not found in path. Exiting.\n" if (!$valid_path);
	$valid_path = `which mash`;
	@tmp = split(/\//, $valid_path);
	$mash_dir = join("\/", @tmp[0..(@tmp-2)])."/";
    }
}
else{
    $mash_dir = $main_dir . $mash_dir . "/" if(substr($mash_dir, 0, 1) ne "/");
    if (! -e $mash_dir . "/mash") {
        die "mash not found at: ".$mash_dir."\n";
    }
}

#mummer
if(!defined($mummer_dir)){
    $mummer_dir = "$opera_ms_dir/utils/MUMmer3.23/";
    if(-d $mummer_dir){
	my $exe_check = "${mummer_dir}nucmer --version";
	run_exe($exe_check);
    }
    if ($? || ! -d $mummer_dir ){
	print STDERR "\nnucmer found in $mummer_dir is not functional. Checking path for mummer. \n";
	print STDERR "\nChecking path for mummer. \n";
	$mummer_dir = " ";
	my $valid_path = which('mummer');
	die "mummer not found in path. Exiting.\n" if (!$valid_path);
	$valid_path = `which mummer`;
	@tmp = split(/\//, $valid_path);
	$mummer_dir = join("\/", @tmp[0..(@tmp-2)])."/";
    }
}
else{
    $mummer_dir = $main_dir . $mummer_dir . "/" if(substr($mummer_dir, 0, 1) ne "/");
    if (! -e $mummer_dir . "/mummer") {
        die "mummer not found at: ".$mummer_dir."\n";
    }
}

if(!defined($megahit_dir)){
    $megahit_dir = "$opera_ms_dir/utils/megahit_v1.0.4-beta_LINUX_CPUONLY_x86_64-bin/";
    my $exe_check = "$megahit_dir/megahit --version";
    run_exe($exe_check);

    if ($?){
        print STDERR "\nMegaHit found in $opera_ms_dir/utils is not functional. Checking path for megahit. \n";
        $megahit_dir = "";
        my $valid_path = which('megahit');
        die "megahit not found in path. Exiting.\n" if (!$valid_path);
	$valid_path = `which megahit`;
	@tmp = split(/\//, $valid_path);
	$megahit_dir = join("\/", @tmp[0..(@tmp-2)])."/";
    }
}
else{
    if (! -e $megahit_dir) {
        die "Megahit not found at: ".$megahit_path."\n";
    }
}
 
#Check for racon.
if(!defined($racon_dir)){
    $racon_dir = "$opera_ms_dir/utils/";
    my $exe_check = "${racon_dir}racon 2>&1 | grep options";
    #my $exe_check = "source activate nanopore;${racon_dir}racon 2>&1 | grep options;source deactivate";
    run_exe($exe_check);

    if($?){
        print STDERR "\nracon found in $opera_ms_dir/utils is not functional. Checking path for racon. \n";
        $racon_dir = "";
        my $valid_path = which("racon"); 
        die "racon not found in path. Exiting.\n" if (!$valid_path);
    }
}
else{
    $racon_dir = $main_dir . $racon_dir . "/" if(substr($racon_dir, 0, 1) ne "/");

    if (! -e $racon_dir . "/racon") {
        die "racon not found at: ".$racon_dir."\n";
    }
}

if( ! -e $illum_read1) {
    die "\nError : $illum_read1 - illumina read 1 file does not exist\n";
}

if( ! -e $illum_read2) {
    die "\nError : $illum_read2 - illumina read 2 file does not exist\n";
}

if( ! -e $long_read_file) {
    die "\nError : $long_read_file - long read file does not exist\n";
}

if( ! -d "$opera_ms_dir/$opera_version"){
    die "\nError : the $opera_version folder is not found. \n";
}
    
#my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -V -l mem_free=20G,h_rt=500:0:0 -pe OpenMP 20";
my $command;

### opera config
#
my ($start_time, $end_time);
my $flag_megahit_assembly = 0;
my $megahit_out_dir = "$output_dir/$inter/megahit_assembly";
if($STAGE_TO_RUN eq "ALL" || $STAGE_TO_RUN eq "SHORT_READ_ASSEMBLY"){
    $STAGE_TO_RUN = "ALL" if($run_following == 1);
    $start_time = time();
    #Stage 0: short assembly using megahit if the contig file is not provided in the config file
    #print STDERR $contigs_file . " " . $main_dir . "\n";<STDIN>;
    if(! defined $contigs_file){
	print STDOUT "\n---  STEP 0: STARTING SHORT READ ASSEMBLY  ---\n";
	$flag_megahit_assembly = 1;
	#
	$contigs_file = "$megahit_out_dir/final.contigs.fa";
	$contigs_file_sed = $contigs_file;
	$contigs_file_sed =~ s/\//\\\//g;      #replace / with \/
	#
	if(! -e $contigs_file){
	    $command = "$megahit_dir/megahit --num-cpu-threads $num_processor -1 $illum_read1 -2 $illum_read2 -o $megahit_out_dir > $output_dir/$inter/megahit.out 2> $output_dir/$inter/megahit.err";
	    run_exe("$command");
	    if($?){
		die "Error Megahit assembly. Please see $output_dir/$inter/megahit.out and $output_dir/$inter/megahit.err for details.\n";
	    }
	}
	$end_time = time();
	print STDOUT "***  Elapsed time: " . ($end_time - $start_time) . "s\n";#<STDIN>;
    }
}

if(! defined $contigs_file){
    $contigs_file = "$megahit_out_dir/final.contigs.fa";
    $contigs_file_sed = $contigs_file;
    $contigs_file_sed =~ s/\//\\\//g;      #replace / with \/
}

#Delete all the intermidate file -_-
#$command = "rm -r $output_dir/$inter;mkdir -p $output_dir";
#run_exe($command);

my $DIR_COV = "$output_dir/$inter/coverage_estimation/";
my $DIR_MAPPING = "$output_dir/$inter/long-read-mapping";
if($STAGE_TO_RUN eq "ALL" || $STAGE_TO_RUN eq "CONTIG_GRAPH"){
    $STAGE_TO_RUN = "ALL" if($run_following == 1);
    $start_time = time;
    #Stage 1: Construction of scaffold graph and contig coverage estimation
    #Use for:
    #mapping of short reads => bwa
    #mapping of long reads => blasr
    #compute long read links
    #
    print STDOUT "\n---  STEP 1: CONSTRUCTION OF THE CONTIG GRAPH AND CONTIG COVERAGE ESTIMATION  ---\n";
    ### Run opera-lr
    #
    $command = "mkdir -p $DIR_MAPPING";
    run_exe($command);
    #print STDERR "\n--- STARTING OPERA-long-read.pl ----\n";
    $start_time_sub = time;
    #$command= "${opera_ms_dir}${opera_version}/bin/OPERA-long-read.pl --contig-file $contigs_file --kmer $kmer_size --long-read-file $long_read_file --output-prefix opera --output-directory $DIR_MAPPING --num-of-processors $num_processor --opera ${opera_ms_dir}${opera_version}/bin/ --illumina-read1 $illum_read1 --illumina-read2 $illum_read2 --samtools-dir $samtools_dir --blasr $blasr_dir --short-read-tooldir $short_read_tool_dir --skip-opera 1  2> $DIR_MAPPING/log.err";

    $long_read_mapper_option = "--minimap2 $minimap2_dir" if($mapper eq "minimap2");
    $long_read_mapper_option = "--blasr $blasr_dir" if($mapper eq "blasr");
    
    $command= "${opera_ms_dir}${opera_version}/bin/OPERA-long-read.pl --contig-file $contigs_file --kmer $kmer_size --long-read-file $long_read_file --output-prefix opera --output-directory $DIR_MAPPING --num-of-processors $num_processor --opera ${opera_ms_dir}${opera_version}/bin/ --illumina-read1 $illum_read1 --illumina-read2 $illum_read2 --samtools-dir $samtools_dir $long_read_mapper_option --short-read-tooldir $short_read_tool_dir --skip-opera 1 > log.out  2> $DIR_MAPPING/log.err";
    
    run_exe($command);
    #
    #Error in STAGE 1
    if($?){
	die "Error in the long read processing. Please see log for details.\n";
    }
    $end_time_sub = time;
    print STDOUT "***  Contig graph generation Elapsed time: " . ($end_time_sub - $start_time_sub) . "s\n";

    ### Softlink opera.bam
    $command="ln -s $DIR_MAPPING/opera.bam $output_dir/contigs.bam";
    run_exe($command);
    
    ### runOperaMS
    #Generate the config file
    $command = "cp $OPERAMS_config_name $opera_ms_config_file";
    run_exe($command);
    open(CONF, '>>',  $opera_ms_config_file);
    if($flag_megahit_assembly){#Add the contig file name to the config file
	print CONF "\nCONTIGS_FILE $contigs_file\n";
    }
    print CONF "\nLIB $output_dir/contigs.bam\n";
    print CONF "\nMAPPING_FILES $output_dir/contigs.bam\n";
    #
    print CONF "SAMTOOLS_DIR $samtools_dir\n";
    
    close(CONF);
    ##
    #SIGMA on short reads used to estimate the contig read coverage 
    #bundling of short reads (can be use if fix the problem with short reads)
    $start_time_sub = time;
    run_exe("mkdir $DIR_COV") if(! -d $DIR_COV);
    $command = "perl ${opera_ms_dir}bin/runOperaMS.pl $opera_ms_config_file > $DIR_COV/log.out 2> $DIR_COV/log.err";
    run_exe($command);
    if($?){
	die "Error in edge bundling and SIGMA execution. Please see $DIR_COV/log.out and $DIR_COV/log.err for details.\n";
    }
    $end_time_sub = time;
    print STDOUT "***  Coverage estimation Elapsed time: " . ($end_time_sub - $start_time_sub) . "s\n";
    #Remove edge with an outlier edge support/contig coverage ratio => good to remove error due to missmapping
    $command = "perl ${opera_ms_dir}bin/filter_edge_coverage_ratio.pl $DIR_MAPPING $DIR_COV/contigs_*";
    run_exe($command);
    $end_time = time;
    print STDOUT "***  Elapsed time: " . ($end_time - $start_time) . "s\n";
}


my $DIR_SIGMA = "$output_dir/$inter/sigma";
print STDERR " *** $DIR_SIGMA\n";
if($STAGE_TO_RUN eq "ALL" || $STAGE_TO_RUN eq "CLUSTERING"){
    $STAGE_TO_RUN = "ALL" if($run_following == 1);
    $start_time = time;
    ### Run sigma for long read (filtered edge obtain using filter_edge_coverage_ratio.pl###
    #
    
    if (-e $DIR_SIGMA ){
	$command = "rm -r $DIR_SIGMA";
	run_exe($command);
    }
    
    $command="mkdir -p $DIR_SIGMA";
    run_exe($command);
    
    generate_sigma_config_file($DIR_COV, $DIR_MAPPING, $DIR_SIGMA);

    #Run sigma
    print STDOUT "\n---  STEP 2: HIERARCHICAL CLUSTERING  ---\n";
    $command="${opera_ms_dir}bin/sigma $DIR_SIGMA/sigma.config > $DIR_SIGMA/log.out 2> $DIR_SIGMA/log.err";
    run_exe($command);

    if($?){
	die "Error in during SIGMA execution. Please see $DIR_SIGMA/log.out and $DIR_SIGMA/log.err for details.\n";
    }
    
    
    #Refine the sigma R parameter
    $command = "${opera_ms_dir}bin/refine_r_estimate.pl $DIR_COV/contigs_340_80 $DIR_SIGMA $opera_ms_dir";
    run_exe($command);
    
    #***
    $end_time = time;
    print STDOUT "***  Elapsed time: " . ($end_time - $start_time) . "s\n";
}


my $DIR_REF_CLUSTERING = "$output_dir/$inter/reference_mapping";
if($STAGE_TO_RUN eq "ALL" || $STAGE_TO_RUN eq "REF_CLUSTERING"){
    $STAGE_TO_RUN = "ALL" if($run_following == 1);
    $start_time = time;
    #Run reference mapping to merge clusters, use both long read and short read.
    #Use the short read and add a link if support > 5, gap estimated < 500 and validated on a reference
    #NEED TO REMOVE LINKS SUPPORTED BY BOTH SHORT AND LONG READS
    print STDOUT "\n---  STEP 3: SPECIES CLUSTERING  ---\n";
    $command="mkdir -p $DIR_REF_CLUSTERING";
    run_exe($command);
    $command="perl ${opera_ms_dir}bin/sequence_similarity_clustering.pl $output_dir/$inter $contigs_file $num_processor $opera_ms_dir $mash_dir $mummer_dir 2> $DIR_REF_CLUSTERING/log.err";
    run_exe($command);
    if($?){
	die "Error in during reference based clustering. Please see $DIR_REF_CLUSTERING/log.err for details.\n";
    }
    $end_time = time;
    print STDOUT "***  Elapsed time: " . ($end_time - $start_time) . "s\n";
}


############################################################## MULTIPLE STRAIN SPECIES ASSEMBLY
my $DIR_STRAIN = "$output_dir/$inter/strain_analysis";
if($strain_clustering eq "YES" && ($STAGE_TO_RUN eq "ALL" || $STAGE_TO_RUN eq "STRAIN_CLUSTERING")){
    $STAGE_TO_RUN = "ALL" if($run_following == 1);
    $start_time = time;
    #Plug-in the strain level assembly
    #->if a species contains the strain each strain will be assembled independently by opera-lg
    #->if a species does not contain the strain they are all pulled together and assembled in a single run by opera-lg
    print STDOUT "\n---  STEP 4: STRAIN LEVEL CLUSTERING AND ASSEMBLY  ---\n";
    run_exe("rm -r $DIR_STRAIN;mkdir -p $DIR_STRAIN");
    $command = "perl ${opera_ms_dir}bin/coverage_clustering.pl $output_dir/$inter $DIR_STRAIN $contigs_file $opera_ms_dir 2>  $DIR_STRAIN/log.err";
    run_exe($command);
    if($?){
	die "Error in during strain clustering. Please see $DIR_STRAIN/log.err for details.\n";
    }
    
    #Assembly of species with multiple strains
    $global_r_value = `tail -n1 $DIR_SIGMA/r_value_evaluation_summary.dat | cut -f1`;chop $global_r_value;
    #
    $r_value_interval = `head -n1 $DIR_SIGMA/r_estimate_value.dat`;chop $r_value_interval;
    @r_value_interval_tab = split(/\s+/, $r_value_interval);
    $r_value_step = ($r_value_interval_tab[1] - $r_value_interval_tab[0]) / 10;
    #
    open(FILE, "$DIR_STRAIN/reference_length.dat");
    my $cmd_strain = "$DIR_STRAIN/assembly.cmd";
    open(CMD_STRAIN, ">$cmd_strain");
    while(<FILE>){
	@line = split(/\t/, $_);
	$species = $line[0];
	$species_length = $line[1];
	$command = "perl ${opera_ms_dir}bin/cluster_strain.pl $DIR_STRAIN/$species $species_length $global_r_value $r_value_step $opera_ms_dir";
	print CMD_STRAIN $command . "\n";
	#run_exe($command);
    }
    close(FILE);
    close(CMD_STRAIN);
    run_exe("cat $cmd_strain | xargs -L 1 -P $num_processor -I COMMAND sh -c \"COMMAND\" 2> $cmd_strain.log");
    if($?){
	die "Error in during strain clustering assembly. Please see $cmd_strain.log for details.\n";
    }
    $end_time = time;
    print STDOUT "***  Elapsed time: " . ($end_time - $start_time) . "s\n";
}
    

############################################################## SINGLE STRAIN SPECIES ASSEMBLY
my $DIR_OPERA_LR = "$output_dir/$inter/opera_long_read";
if($STAGE_TO_RUN eq "ALL" || $STAGE_TO_RUN eq "ASSEMBLY"){
    $STAGE_TO_RUN = "ALL" if($run_following == 1);
    $start_time = time;
    print STDOUT "\n---  STEP 5: ASSEMBLY OF SINGLE STRAIN SPECIES  ---\n";
    $command="mkdir $DIR_OPERA_LR";
    run_exe($command);

    generate_opera_config_file($DIR_MAPPING, $DIR_REF_CLUSTERING, $DIR_OPERA_LR);
    
    #Assembly of all the other contigs were no multiple strain of the same species have inferred
    #Filter that remove contigs with a coverage 1.5 times higher than the mean and contig that are defnied as repeat based on the reference mapping
    #remove the contigs that belong to clusters that are associated to species with multiple strains
    my $reference_cluster_file = "$DIR_REF_CLUSTERING/clusters_seq_similarity";
    $reference_cluster_file = "$DIR_REF_CLUSTERING/clusters_single_strain" if($strain_clustering eq "YES");
    $contig_coverage_file = "$DIR_COV/contigs_$contig_window_len\_$contig_edge_len";
    $command = "${opera_ms_dir}bin/filter_cluster_coverage.pl $contig_coverage_file $reference_cluster_file $DIR_REF_CLUSTERING 1.5 $DIR_REF_CLUSTERING/NO_REPEAT $DIR_REF_CLUSTERING/NUCMER_OUT/ $mummer_dir";
    #$command="outcontig=`ls $DIR_COV/contigs_\*`; ${opera_ms_dir}bin/filter_cluster_coverage.pl `echo \$outcontig` $reference_cluster_file $DIR_REF_CLUSTERING 1.5 $DIR_REF_CLUSTERING/NO_REPEAT $DIR_REF_CLUSTERING/NUCMER_OUT/";
    
    run_exe($command);
    
    
    #finally, run opera
    $command="${opera_ms_dir}${opera_version}/bin/OPERA-LG $DIR_OPERA_LR/opera.config > $DIR_OPERA_LR/log.out 2> $DIR_OPERA_LR/log.err";
    run_exe($command);
    if($?){
	die "Error in during OPERA-LG. Please see $DIR_OPERA_LR/log.out and $DIR_OPERA_LR/log.err for details.\n";
    }
    $end_time = time;
    print STDOUT "***  Elapsed time: " . ($end_time - $start_time) . "s\n";
}

my $DIR_GAPFILLING = "$output_dir/$inter/opera_long_read/GAPFILLING";
if($STAGE_TO_RUN eq "ALL" || $STAGE_TO_RUN eq "GAP_FILLING"){
    $STAGE_TO_RUN = "ALL" if($run_following == 1);
    $start_time = time;
    print STDOUT "\n---  STEP 6: GAP FILLING  ---\n";
    #Add the scaffold obtain during the strain level assembly back in
    
    my $opera_scaff_file = "$DIR_OPERA_LR/scaffolds.scaf";
    if($strain_clustering eq "YES"){
	my $opera_scaff_file_single = "$DIR_OPERA_LR/scaffolds_single.scaf";
	run_exe("mv $opera_scaff_file $opera_scaff_file_single") if(! -e $opera_scaff_file_single);
	run_exe("cp $opera_scaff_file_single $opera_scaff_file");
    }
    
    my $opera_scaff_seq_file = "$DIR_OPERA_LR/scaffoldSeq.fasta";
    if($strain_clustering eq "YES"){
	my $opera_scaff_seq_file_single = "$DIR_OPERA_LR/scaffoldSeq_single.fasta";
	run_exe("mv $opera_scaff_seq_file $opera_scaff_seq_file_single") if(! -e $opera_scaff_seq_file_single);
	run_exe("cp $opera_scaff_seq_file_single $opera_scaff_seq_file");
    }

    my $opera_contig_size_file = "$DIR_OPERA_LR/contigs";
    if($strain_clustering eq "YES"){
	my $opera_contig_size_file_single = "$DIR_OPERA_LR/contigs_single";
	run_exe("mv $opera_contig_size_file $opera_contig_size_file_single") if(! -e $opera_contig_size_file_single);
	run_exe("cp $opera_contig_size_file_single $opera_contig_size_file");
    }

    #Add the multiple strain genome scaffold to the pool
    if($strain_clustering eq "YES"){
	opendir(DIR, "$DIR_STRAIN");
	my @strain_dir = readdir(DIR);
	close(DIR);
	
	my $strain_ID = 1;
	my $scaff_strain_id = 1;
	
	#print STDERR " **** @strain_dir\n";
	
	foreach $strain (@strain_dir){
	    $strain_ID = 1;
	    $strain_scaff_file = "$DIR_STRAIN/$strain/STRAIN_$strain_ID/scaffolds.scaf";
	    $strain_scaff_seq_file = "$DIR_STRAIN/$strain/STRAIN_$strain_ID/scaffoldSeq.fasta";
	    $strain_contig_size_file = "$DIR_STRAIN/$strain/STRAIN_$strain_ID/contigs";
	    #print STDERR " *** $strain_scaff_file\n";
	    while( -e $strain_scaff_file){
		run_exe("sed 's/>opera/>strain$scaff_strain_id\_opera/' $strain_scaff_file >> $opera_scaff_file");
		run_exe("sed 's/>opera/>strain$scaff_strain_id\_opera/' $strain_scaff_seq_file >> $opera_scaff_seq_file");
		run_exe("grep -v Length $strain_contig_size_file >> $opera_contig_size_file");
		#
		$strain_ID++;
		$scaff_strain_id++;
		$strain_scaff_file = "$DIR_STRAIN/$strain/STRAIN_$strain_ID/scaffolds.scaf";
		$strain_scaff_seq_file = "$DIR_STRAIN/$strain/STRAIN_$strain_ID/scaffoldSeq.fasta";
	    }
	}
    }

    #exit(0);

    #perform the gap fill
    $command="rm $output_dir/scaffoldSeq.fasta; rm $output_dir/scaffoldSeq.fasta.stats; ln -s $DIR_OPERA_LR/scaffoldSeq.fasta $output_dir/scaffoldSeq.fasta";
    run_exe($command);
    
    #$command="python ${opera_ms_dir}/SCRIPTS/Pipeline.py --analysis_dir ${output_dir} --main_script_dir ${opera_ms_dir}/SCRIPTS/ --aux_script_dir $opera_ms_dir/$opera_version/bin/ --minimap2_dir $minimap2_dir --pilon_path $pilon_path --racon_dir $racon_dir --bwa_dir $short_read_tool_dir --samtools_dir $samtools_dir --edge_file $lr_output_dir/edge_read_info.dat --scaffolds_file $opera_scaff_file --contig_file $contigs_file --long_read_file $long_read_file --short_read_file_1 $illum_read1  --short_read_file_2 $illum_read2 --num_processor $num_processor";
    $command = "perl $opera_ms_dir/bin/match_scaffolds_clusters.pl $DIR_REF_CLUSTERING $num_processor $opera_ms_dir $DIR_OPERA_LR/scaffoldSeq.fasta $mummer_dir 2> $DIR_REF_CLUSTERING/s_info.err";
    run_exe($command);
    
    $command="$opera_ms_dir/$opera_version/bin/gapfilling.pl $DIR_GAPFILLING $contigs_file $long_read_file $num_processor $opera_ms_dir/$opera_version/bin/ $racon_dir $minimap2_dir $mummer_dir 2> $DIR_OPERA_LR/gapfilling.log";
    run_exe($command);

    if($?){
	die "Error in during gapfilling. Please see $DIR_OPERA_LR/gapfilling.log for details.\n";
    }
    
    $end_time = time;
    print STDOUT "***  Elapsed time: " . ($end_time - $start_time) . "s\n";
}
 
if($STAGE_TO_RUN eq "ALL" || $STAGE_TO_RUN eq "INFO"){
    $STAGE_TO_RUN = "ALL" if($run_following == 1);
    $start_time = time;
    print STDOUT "\n---  STEP 7: ASSEMBLY INFO  ---\n";
    #For the scaffold_info.dat file generation
    $command = "perl $opera_ms_dir/bin/match_scaffolds_clusters.pl $DIR_REF_CLUSTERING $num_processor $opera_ms_dir $DIR_OPERA_LR/scaffoldSeq.fasta.filled $mummer_dir 2> $DIR_REF_CLUSTERING/s_info.err";
    run_exe($command);
    if($?){
	die "Error in during scaffold info genration. Please see $DIR_REF_CLUSTERING/s_info.err for details.\n";
    }
    
    #Remove the contigs that have been used during gapfiling in the .fa file
    #$command = "perl $opera_ms_dir/bin/clean_scaffolds_repeats.pl $output_dir/$inter";
    #run_exe($command);
    
    #$command = "mv $output_dir/lib_* $output_dir/$inter; rm -r $output_dir/$inter/scaffolds; mv $output_dir/runOperaMS.config $output_dir/$inter/";
    #run_exe($command);

    #Rename the scaffold sequence into contig after gap filling
    run_exe("sed 's/_scaffold_/_contig_/' $DIR_OPERA_LR/scaffoldSeq.fasta.filled > $output_dir/contig.fasta");
    run_exe("sed 's/_scaffold_/_contig_/' ${output_dir}/scaffold_info.txt > ${output_dir}/contig_info.txt");
    run_exe("rm ${output_dir}/scaffold_info.txt");
    #$command = "rm $output_dir/contigs.bam";#; ln -s $output_dir/$inter/opera_long_read/scaffoldSeq.fasta.filled $output_dir/scaffoldSeq.fasta.filled";
    #run_exe($command);
    
    #$command="less ${output_dir}/$inter/opera_long_read/statistics";
    #$str = `sort -k2,2nr  ${output_dir}/scaffold_info.txt | head -n10 | cut -f1,2`;
    run_exe($command);
    $end_time = time;
    print STDOUT "***  Elapsed time: " . ($end_time - $start_time) . "s\n";
    #Get the stats
    my @all_size = ();
    open(FILE, "${output_dir}/contig_info.txt");
    my $contig_1mb = 0;
    my $contig_100kb = 0;
    my $assembly_size = 0;
    <FILE>;
    while(<FILE>){
	@line = split(/\t/, $_);
	$size = $line[1];
	$contig_1mb++ if($size > 1000000);
	$contig_100kb++ if($size > 100000);
	$assembly_size += $size;
	push(@all_size, $size);
    }
    close(FILE);
	
    my @sort_tab = sort {$b <=> $a} @all_size;
    my $nb_contig = @all_size+0;
    $nb_seq = compute_Nx(50, $assembly_size, \@sort_tab);$n50 = $sort_tab[$nb_seq];
    my $str_stats =  " *** Assembly stats:\n" . 
	" *** *** " . "Number of contigs: " . $nb_contig . " \n". 
	" *** *** " . "Assembly size: $assembly_size bp\n". 
	#"\t" . "min contig size $sort_tab[$nb_contig-1] bp, ".
	" *** *** " . "Max contig size: $sort_tab[0] bp\n".
	" *** *** " . "Contig(s) longer than 1Mbp: $contig_1mb \n".
	" *** *** " . "Contig(s) longer than 100kb: $contig_100kb contig\n".
	" *** *** " . "Contig N50: $n50 bp\n";
    
    print STDOUT $str_stats . "\n";
    open(OUT, ">${output_dir}/assembly.stats");
    print OUT $str_stats . "\n";
    close(OUT);
}

print STDERR "\n*************OPERA-MS DONE***************\n" if($STAGE_TO_RUN eq "ALL");

sub generate_sigma_config_file{
    my ($dir_cov_estimate, $dir_mapping, $dir_sigma) = @_;
    #
    #long read sigma config
    open(OUT, ">$dir_sigma/sigma.config");
    open(FILE, "$dir_cov_estimate/sigma.config");
    my $max_edge_set = 5;
    my $str_edge = "";
    while(<FILE>){
	@line = split(/\=/, $_);
	
	#no mapping file
	next if($line[0] eq "mapping_files");

	#Change the output
	if($line[0] eq "output_dir"){
	    print OUT "output_dir=$dir_sigma\n";
	    next;
	}
	
	#Get the edge file name
	if($line[0] eq "edges_files"){
	    for(my $i = 0; $i <= $max_edge_set; $i++){
		$str_edge .= "$dir_sigma/pairedEdges_i$i/pairedEdges_i$i,";
		#Creat the directory
		run_exe("mkdir $dir_sigma/pairedEdges_i$i");
		#Link to the correct file
		run_exe("ln -s $dir_mapping/pairedEdges_i$i $dir_sigma/pairedEdges_i$i/pairedEdges_i$i");
		#create the lib file
		open(LIB, ">$dir_sigma/pairedEdges_i$i/lib.txt");
		@val_mean_std = `grep -A2 pairedEdges_i$i $dir_mapping/config | cut -d "=" -f 2`;
		#
		print LIB "For the mapping file: \n";
		print LIB "Mean length of the library is: $val_mean_std[1]";
		print LIB "Standard deviation of the library is: $val_mean_std[2]";
		print LIB "Orientation of paired-reads is: forward ->...->  1st read...2nd read\n";
		print LIB "Number of used paired reads for estimating the mean and the standard deviation of the library is: NA\n";
		print LIB "Percentage of outliers is: NA%";
	    }
	    chop $str_edge;
	    print OUT "edges_files=$str_edge\n";
	    next
	}
	
	print OUT $_;
	    
    }
    #print OUT "samtools_dir=$samtools_dir\n";
    close(OUT);
    close(FILE);
}

sub generate_opera_config_file{
    my ($dir_mapping, $dir_ref_clustering, $dir_opera_lr) = @_;
    #
    #long read sigma config
    open(FILE, "$output_dir/$inter/lib_1_bundles/lib_1.ebconfig");
    open(OUT, ">$dir_opera_lr/opera.config");
    my $max_edge_set = 5;
    while(<FILE>){

	if($_ eq "[LIB]\n"){
	    
	    for(my $i = 0; $i <= $max_edge_set; $i++){
		#create the lib file
		@val_mean_std = `grep -A2 pairedEdges_i$i $dir_mapping/config | cut -d "=" -f 2`;
		print STDERR " *** $dir_mapping/config\n";
		print OUT "[LIB]\n";
		print OUT "map_type=opera\n";
		print OUT "map_file=$dir_ref_clustering/NO_REPEAT/filtered_pairedEdges_i$i\n";
		print OUT "lib_mean=$val_mean_std[1]\n";
		print OUT "lib_std=$val_mean_std[2]\n";
	    }
	    next;
	}
	
	@line = split(/\=/, $_);

	#no mapping file
	next if($line[0] eq "mapping_files" || $line[0] eq "map_file" || $line[0] eq "filter_repeat");

	if($line[0] eq "output_folder"){
	    print OUT "output_folder=$dir_opera_lr\n";
	    next;
	}

	if($strain_clustering eq "YES" && $line[0] eq "contig_file"){
	    print OUT "contig_file=$dir_ref_clustering/single_strain_species.fa\n";
	    next;
	}

	if($line[0] eq "cluster_threshold"){
	    print OUT "cluster_threshold=1\n";
	    next;
	}
	
	print OUT $_;
	
    }
    print OUT "keep_repeat=yes\n";
    #print OUT "samtools_dir=$samtools_dir\n";
    close(FILE);
    close(OUT);
    
}

sub generate_ms_config_file{
    
}

sub compute_Nx{
    my ($x, $length, $ptr_tab) = @_;

    my $nb_x_base = $length * ($x / 100);
    
    my $total = 0;
    
    my $i = 0;
    for($i = 0; $total < $nb_x_base; $i++){
	$total += $ptr_tab->[$i];
    }
    
    #print "$x $pos_x $size\n";<STDIN>;

    #return $ptr_tab->[$i-1];
    return $i -1;
}




sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return . "\n" if($run);
    return $return;
}

    

