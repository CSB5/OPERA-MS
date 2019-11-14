#use strict;
use warnings;
use Cwd;
use Switch;
use File::Which;
use File::Spec;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
#
use Statistics::Basic;
use Statistics::R;
#
use Getopt::Std;
use Getopt::Long;
use Getopt::Long qw(GetOptionsFromArray);
#use Getopt::long::Parser;

my %opera_ms_option = ();
my $opera_ms_full_path = dirname(rel2abs($0)) . "/";
$opera_ms_option{"OPERA_MS_DIR"} = $opera_ms_full_path;
$opera_ms_option{"VERSION"} = "v0.9.0";

my %opera_ms_dependency = ();

require "$opera_ms_full_path/bin/test_time.pl";

#Error message
my $help_line = "To configure OPERA-MS, please look at the example config file inside the OPERA-MS folder.\nUsage: \n\npath/OPERA-MS/OPERA-MS.pl <config_file>\n\n";
my $incorrect_arguments="Please input the correct arguments. Refer to the documentation for more information or use:\n\n path/OPERA-MS/OPERA-MS.pl -h. \n\n";
if ( @ARGV == 0 ){
    print "\n" . print_help();exit(0);
}

#To setup and check if the depndency are functionning
check_dependency(\%opera_ms_option, \%opera_ms_dependency, "CHECK_DEPENDENCY") if($ARGV[0] eq "CHECK_DEPENDENCY");

#read from config file
if ( $ARGV[0] ne "--help" && @ARGV <= 2){
    read_config_file(\%opera_ms_option, @ARGV);
}
else{
    read_argument(\%opera_ms_option,\@ARGV);#Write config file
    #read_config_file(\%opera_ms_option, $opera_ms_option{"CONFIG_PATH"});
    #die $incorrect_arguments;
}
#Make paths absolute if they are not already.
setup_directory_path(\%opera_ms_option);

#
#Check if dependencies are found. Otherwise, die.
check_dependency(\%opera_ms_option, \%opera_ms_dependency, $opera_ms_option{"STAGE"});    

#Check data integrity
check_data(\%opera_ms_option);

setup_opera_ms_database(\%opera_ms_option);

#Starting OPERA-MS pipeline
my ($start_time, $end_time);

print STDERR "version\t" . $opera_ms_option{"VERSION"} . "\n";

short_read_assembly(\%opera_ms_option, \%opera_ms_dependency);

assembly_graph_creation(\%opera_ms_option, \%opera_ms_dependency);

hierarchical_clustering(\%opera_ms_option, \%opera_ms_dependency);

reference_clustering(\%opera_ms_option, \%opera_ms_dependency);#Add check for reference clustering

strain_clustering_and_assembly(\%opera_ms_option, \%opera_ms_dependency);

gap_filling(\%opera_ms_option, \%opera_ms_dependency);
#
generate_assembly_stats(\%opera_ms_option, \%opera_ms_dependency);
write_final_assembly(\%opera_ms_option, \%opera_ms_dependency);
#
#
polishing(\%opera_ms_option, \%opera_ms_dependency);

    
$time = localtime;
print "\n[$time]\tOPERA-MS completed\n";# if($STAGE_TO_RUN eq "ALL");

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
    my ($dir_mapping, $dir_ref_clustering, $dir_opera_lr, $contig_file) = @_;
    #
    #long read sigma config
    open(FILE, "$dir_opera_lr/../lib_1_bundles/lib_1.ebconfig");
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

	if($line[0] eq "contig_file"){
	    print OUT "contig_file=$contig_file\n";
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


sub compute_Nx{
    my ($x, $length, $ptr_tab) = @_;

    my $nb_x_base = $length * ($x / 100);
    
    my $total = 0;
    
    my $i = 0;
    for($i = 0; $total < $nb_x_base; $i++){
	$total += $ptr_tab->[$i];
    }
    
    return $i -1;
}

sub set_default_value{
    my ($opera_ms_option) = @_;
    #GIVE DEFAULT VALUES
    
    #OPERA-MS option
    $opera_ms_option->{"KMER_SIZE"} = 60;
    $opera_ms_option->{"NUM_PROCESSOR"} = 2;
    $opera_ms_option->{"STAGE_FOLLOW"} = 0;
    $opera_ms_option->{"CONTIG_EDGE_LEN"} = 80;
    $opera_ms_option->{"CONTIG_WINDOW_LEN"} = 340;
    $opera_ms_option->{"CONTIG_LEN_THR"} = 500;
    $opera_ms_option->{"REF_CLUSTERING"} = 1;
    $opera_ms_option->{"STRAIN_CLUSTERING"} = 1;
    $opera_ms_option->{"RESTART"} = 0;
    $opera_ms_option->{"STAGE_FOLLOW"} = 0;
    $opera_ms_option->{"MULTI_SAMPLE"} = 0;
    $opera_ms_option->{"REF_REPEAT_DETECTION"} = 0;
    $opera_ms_option->{"POLISHING"} = 0;
    $opera_ms_option->{"GAP_FILLING"} = 1;
    $opera_ms_option->{"LONG_READ_MAPPER"} = "blasr";
    $opera_ms_option->{"SHORT_READ_ASSEMBLER"} = "megahit";
    #my $genomeDB_kraken = "NULL";
    #my $kraken_exe_dir = "";
    $opera_ms_option->{"FILLED_SCAFF_LENGTH"} = 499;
    $opera_ms_option->{"MASH_DB"} = $opera_ms_option->{"OPERA_MS_DIR"}."/genomeDB_Sketch.msh";
    #
    $opera_ms_option->{"KRAKEN_DB"} = "NULL";#$opera_ms_option->{"OPERA_MS_DIR"}."/genomeDB_krakenh.msh";
}

sub read_config_file{

    my ($opera_ms_option, $config_file, $stage_info) = @_;
    
    $opera_ms_option->{"CONFIG_PATH"} = $config_file;
    
    #get the stage
    my $run_following = 0;
    if(defined $stage_info){
	if(index($stage_info, "+") != -1){
	    chop $stage_info;
	    $run_following = 1;
	}
	$opera_ms_option->{"RESTART"} = 1;
    }
    else{
	$stage_info = "ALL";
    }
    $opera_ms_option->{"STAGE"} = $stage_info;
    
        
    set_default_value($opera_ms_option);
    
    $time = localtime;
    print STDOUT "\n[$time]\tReading config file: ".$config_file."\n";
    die"Config file does not exist. Exiting.\n" if(!(-e $config_file)); 
    
    open($opera_ms_cf, "<", $config_file); 

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
		    $opera_ms_option->{$config_option} = $split_line[1];
		}
		
		case "OUTPUT_DIR" {
		    $opera_ms_option->{$config_option} = $split_line[1];
		}
		
		case "KMER_SIZE" {
		    $opera_ms_option->{$config_option} = $split_line[1];
		}
		
		case "CONTIG_LEN_THR" {
		    $opera_ms_option->{$config_option} = $split_line[1];
		}

		case "CONTIG_EDGE_LEN" {
		    $opera_ms_option->{$config_option} = $split_line[1];
		}

		case "CONTIG_WINDOW_LEN" {
		    $opera_ms_option->{$config_option} = $split_line[1];
		}
		
		#Use kraken during sequence similarity clustering in addation to mash and mummer
		case "KRAKEN_DB"{
		    $opera_ms_option->{$config_option} = $split_line[1];
		}
		
		case "FILLED_SCAFF_LENGTH" {
		    $opera_ms_option->{$config_option} = $split_line[1];
		}
		
		case "LONG_READ_MAPPER"{
		    $opera_ms_option->{$config_option} = $split_line[1];
		    if($opera_ms_option->{"LONG_READ_MAPPER"} ne "blasr" && $opera_ms_option->{"LONG_READ_MAPPER"} ne "minimap2"){
			die "Unkown LONG_READ_MAPPER " . $opera_ms_option->{"LONG_READ_MAPPER"} . ", please use  blasr/minimap2.\n"
		    }
		}

		case "SHORT_READ_ASSEMBLER"{
		    $opera_ms_option->{$config_option} = $split_line[1];
		    if($opera_ms_option->{"SHORT_READ_ASSEMBLER"} ne "megahit" && $opera_ms_option->{"SHORT_READ_ASSEMBLER"} ne "spades"){
			die "Unkown SHORT_READ_ASSEMBLER " . $opera_ms_option->{"SHORT_READ_ASSEMBLER"} . ", please use  megahit/spades.\n"
		    }
		}
		
		case "LONG_READ"{
		    $opera_ms_option->{$config_option} = $split_line[1];
		    #if(index($opera_ms_option->{$config_option}, ".gz") != -1){
		    #die "Compressed long-reads are not supported " . $opera_ms_option->{$config_option} . "\n";
		    #}
		}

		case "ILLUMINA_READ_1"{
		    $opera_ms_option->{$config_option} = $split_line[1];
		    if(index($opera_ms_option->{"ILLUMINA_READ_1"}, ",") != -1){#Their is multiple illumina file provided => this is a multiple sample assembly
			$opera_ms_option->{"MULTI_SAMPLE"} = 1;
		    }
		}

		case "ILLUMINA_READ_2"{
		    $opera_ms_option->{$config_option} = $split_line[1];
		}

		case "NUM_PROCESSOR"{
		    $opera_ms_option->{$config_option} = $split_line[1];
		}

		case "REF_REPEAT_DETECTION"{
		    $opera_ms_option->{$config_option} = 1;
		    #$REF_REPEAT_DETECTION = 1;
		    $opera_ms_option->{$config_option} = 0 if( $split_line[1] eq "NO");
		}
		
		case "POLISHING"{
		    $opera_ms_option->{$config_option} = 1;
		    #$REF_REPEAT_DETECTION = 1;
		    $opera_ms_option->{$config_option} = 0 if( $split_line[1] eq "NO");
		}
		
		case "REF_CLUSTERING"{
		    #$genomeDB_msh = $split_line[1];
		    $opera_ms_option->{$config_option} = 1;
		    #$REF_REPEAT_DETECTION = 1;
		    $opera_ms_option->{$config_option} = 0 if( $split_line[1] eq "NO");
		}
		case "STRAIN_CLUSTERING"{
		    $opera_ms_option->{$config_option} = 1;
		    #$REF_REPEAT_DETECTION = 1;
		    $opera_ms_option->{$config_option} = 0 if( $split_line[1] eq "NO");
		}
		case "GAP_FILLING"{
		    $opera_ms_option->{$config_option} = 1;
		    #$REF_REPEAT_DETECTION = 1;
		    $opera_ms_option->{$config_option} = 0 if( $split_line[1] eq "NO");
		}
		else {
		    die "Config option: ".$config_option." unknown, please check the config file. \nExiting. \n";
		}
	    }
	}
    }

    $opera_ms_option->{"CONTIGS_COV_FILE"} = "contigs_" . $opera_ms_option->{"CONTIG_WINDOW_LEN"} . "_" . $opera_ms_option->{"CONTIG_EDGE_LEN"};
    
    #To compute the total number of stage
    $opera_ms_option->{"NB_TOTAL_STAGE"} = 8;
    #$opera_ms_option->{"NB_TOTAL_STAGE"} = 4 if(! $opera_ms_option->{"REF_CLUSTERING"});
    $opera_ms_option->{"NB_TOTAL_STAGE"}++ if($opera_ms_option->{"POLISHING"});
}

sub setup_directory_path{
    my ($opera_ms_option) = @_;
    #Grab the executing directory; turn all paths into absolute paths.
    my $main_dir = getcwd;
    $main_dir .= "\/";
    $opera_ms_option->{"MAIN_DIR"} = $main_dir;

    my $output_dir = $opera_ms_option{"OUTPUT_DIR"};
    my $inter = "intermediate_files";
    #
    $output_dir = $main_dir . $output_dir."/" if (substr($output_dir, 0, 1) ne "/");
    run_exe("mkdir -p $output_dir/$inter") if(! -d "$output_dir/$inter");
    $opera_ms_option->{"OUTPUT_DIR"} = $output_dir;
    $opera_ms_option->{"INTER_DIR"} = "$output_dir/$inter";
    
    #Make paths absolute if they are not already.
    $opera_ms_option->{"LONG_READ"} = $main_dir . $opera_ms_option->{"LONG_READ"} if(substr($opera_ms_option->{"LONG_READ"}, 0, 1) ne "/");
    $opera_ms_option->{"ILLUMINA_READ_1"} = $main_dir . $opera_ms_option->{"ILLUMINA_READ_1"} if(substr($opera_ms_option->{"ILLUMINA_READ_1"}, 0, 1) ne "/");
    $opera_ms_option->{"ILLUMINA_READ_2"} = $main_dir . $opera_ms_option->{"ILLUMINA_READ_2"} if(substr($opera_ms_option->{"ILLUMINA_READ_2"}, 0, 1) ne "/");
    if(defined $opera_ms_option->{"CONTIGS_FILE"}){
	$opera_ms_option->{"CONTIGS_FILE"} = $main_dir . $opera_ms_option->{"CONTIGS_FILE"} if(substr($opera_ms_option->{"CONTIGS_FILE"}, 0, 1) ne "/");
    }
}

sub check_data{
    my ($opera_ms_option) = @_;
    
    data_exists($opera_ms_option->{"ILLUMINA_READ_1"});
    data_exists($opera_ms_option->{"ILLUMINA_READ_2"});
    data_exists($opera_ms_option->{"LONG_READ"});
    #
    if ($opera_ms_option->{"LONG_READ"} =~ /\.gz$/){
	die $opera_ms_option->{"LONG_READ"} . " appears to be in gzipped format. The current version of OPERA-MS does not accept compressed long read file, please unzip.\n";
    }
}

sub data_exists{
    my ($data_file) = @_;
    if(! -e $data_file){
	die "\nError : $data_file does not exist\n";
    }
}
    
sub check_dependency{
    my ($opera_ms_option, $opera_ms_dependency, $stage) = @_;
    #
    my $opera_ms_dir = $opera_ms_option->{"OPERA_MS_DIR"};
    my $utils_dir = "$opera_ms_dir/utils/";
    #
    
    #bundler
    #check_software("$opera_ms_dir/bin", "bundler", $sofware_name, $test_command, $opera_ms_dependency) = @_;
    
    #For the perl verion
    if($stage eq "CHECK_DEPENDENCY"){
	my $perl_path = "$utils_dir/perl";
	`rm $perl_path` if(-e "$perl_path");
	`ln -s $^X $perl_path`;
    }
    
    #OPERA-LG
    check_software("$opera_ms_dir/OPERA-LG/bin/", "OPERA-LG", "OPERA-LG", $opera_ms_dependency, $stage);

    #Check for sigma
    check_software("$opera_ms_dir/bin/", "sigma", "sigma help", $opera_ms_dependency, $stage);
    
    #samtools
    check_software($utils_dir, "samtools", "samtools 2>&1 | grep Version", $opera_ms_dependency, $stage);
    
    #blasr
    check_software($utils_dir, "blasr", "blasr -h", $opera_ms_dependency, $stage);
    
    #bwa
    check_software($utils_dir, "bwa", "bwa 2>&1 | grep Version", $opera_ms_dependency, $stage);
        
    #minimap2
    check_software($utils_dir, "minimap2", "minimap2 --version", $opera_ms_dependency, $stage);
    
    #mash
    check_software($utils_dir, "mash", "mash --version", $opera_ms_dependency, $stage);
    
    #mummer
    check_software("$utils_dir/MUMmer3.23/", "mummer", "nucmer --version", $opera_ms_dependency, $stage);

    #megahit
    check_software("$utils_dir/megahit_v1.0.4-beta_LINUX_CPUONLY_x86_64-bin/", "megahit", "megahit --version", $opera_ms_dependency, $stage);

    #spades
    check_software("$utils_dir/SPAdes-3.13.0-Linux/bin/", "spades", "spades.py --version", $opera_ms_dependency, $stage);
    
    #Check for racon.
    check_software($utils_dir, "racon", "racon 2>&1 | grep options", $opera_ms_dependency, $stage);

    #Check for pilon.
    check_software($utils_dir, "pilon", "java -jar $utils_dir/pilon.jar --version", $opera_ms_dependency, $stage);
    #
    #
    $opera_ms_dependency->{"kraken"} = "NULL";

    if($stage eq "CHECK_DEPENDENCY"){
	print STDOUT "\n *** All compiled software are functional.\n";
	print STDOUT " *** Please try to run OPERA-MS on the test dataset.\n\n";
	exit(0);
    }
}

sub check_software{
    my ($software_dir, $software_name, $test_command, $opera_ms_dependency, $stage) = @_;
    
    if($stage eq "CHECK_DEPENDENCY"){
	$exe_check = "${software_dir}$test_command";
	$exe_check = "$test_command" if(index($test_command, "java") != -1);
	#print STDOUT "\n" . $exe_check . "\n";
	`$exe_check 2> /dev/null > /dev/null`;
	
	if($?){
	    print STDERR "\n$software_name found in $software_dir is not functional. Checking path for $software_name.\n";
	    my $valid_path = which("$software_name"); 
	    die "$software_name not found in path. Exiting.\n" if (!$valid_path);
	    print STDERR "\n$software_name found in path.\n";
	    $valid_path = `which $software_name`;chop $valid_path;
	    run_exe("mv ${software_dir}/$software_name ${software_dir}/$software_name.bk");
	    run_exe("ln -s $valid_path $software_dir");
	}
	print STDOUT " *** $software_name functional\n" if($stage eq "CHECK_DEPENDENCY");
    }
    $opera_ms_dependency->{$software_name} = $software_dir;
}


sub setup_opera_ms_database{
    my ($opera_ms_option) = @_;

    #Download the genome reference data base for the first utilasation
    my $opera_ms_dir = $opera_ms_option->{"OPERA_MS_DIR"};
    if(! -e "$opera_ms_dir/database_updated"){
	print STDOUT " *** First utilization: set-up of reference genome databases. Please wait ...\n";
	#
	print STDOUT " *** (1/3) Download of genomeDB_Sketch.msh\n";	
	run_exe("wget --no-check-certificate -O $opera_ms_dir/genomeDB_Sketch.msh https://ndownloader.figshare.com/files/17290511");
	if($?){
	    die"\nUnfortunately the automated download of the genome data base failed.\nYou can complete the installation process manually by:\n\t(1) Downloading the data base files at: https://figshare.com/articles/genomeDB_Sketch_msh/8425280 and https://figshare.com/articles/complete_ref_genomes_tar_gz/9638207\t(2) Moving the files in $opera_ms_dir\n\t(3) Running the command: tar -xvzf $opera_ms_dir/complete_ref_genomes.tar.gz;touch $opera_ms_dir/database_updated\n\n";
	}
	else{
	    #
	    print STDOUT " *** (2/3) Download of complete_ref_genomes.tar.gz\n";	
	    run_exe("wget --no-check-certificate -O $opera_ms_dir/complete_ref_genomes.tar.gz https://ndownloader.figshare.com/files/17290589");
	    #
	}
	print STDOUT " *** (3/3) Extraction of complete_ref_genomes.tar.gz\n";
	run_exe("tar -xvzf $opera_ms_dir/complete_ref_genomes.tar.gz --directory $opera_ms_dir");
	#
	run_exe("rm $opera_ms_dir/complete_ref_genomes.tar.gz");
	#
	run_exe("touch $opera_ms_dir/database_updated");
    }
    
    #
}

#### OPERA-MS pipeline stage

sub short_read_assembly{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    #
    #
    if($opera_ms_option->{"STAGE"} eq "ALL" || $opera_ms_option->{"STAGE"} eq "SHORT_READ_ASSEMBLY"){
	$opera_ms_option->{"STAGE"} = "ALL" if($opera_ms_option->{"STAGE_FOLLOW"} == 1);

	#Check if the contig file is given in the config file
	if(defined $opera_ms_option->{"CONTIGS_FILE"}){
	    $opera_ms_option->{"ASSEMBLED_CONTIG_FILE"} = 0;
	    $time = localtime;
	    print STDOUT "\n[$time]\tShort read assembly [1/" . $opera_ms_option->{"NB_TOTAL_STAGE"} ."]\n";
	    print STDOUT "[$time]\tSkip [contig file provided as input]\n";
	}
	else{
	    $opera_ms_option->{"ASSEMBLED_CONTIG_FILE"} = 1;
	    $start_time = time();
	    if($opera_ms_option->{"SHORT_READ_ASSEMBLER"} eq "megahit"){
		run_megahit($opera_ms_option, $opera_ms_dependency);
	    }
	    if($opera_ms_option->{"SHORT_READ_ASSEMBLER"} eq "spades"){
		run_spades($opera_ms_option, $opera_ms_dependency);
	    }
	    $end_time = time();
	    write_time($opera_ms_option->{"INTER_DIR"}, "short_read_assembly", ($end_time - $start_time));
	}
    }
}


sub run_megahit{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    my $inter_dir = $opera_ms_option->{"INTER_DIR"};
    my $megahit_out_dir = "$inter_dir/megahit_assembly";
    $opera_ms_option->{"CONTIGS_FILE"} = "$megahit_out_dir/final.contigs.fa";
    #Stage 0: short assembly using megahit if the contig file is not provided in the config file
    if(! check_completed($opera_ms_option->{"CONTIGS_FILE"}, "Short read assembly [MEGAHIT]", 1)){
	run_exe("rm -r $megahit_out_dir") if(-d $megahit_out_dir);
	run_exe($opera_ms_dependency->{"megahit"}."/megahit " . 
		"--num-cpu-threads " . $opera_ms_option->{"NUM_PROCESSOR"} . " " .
		"-1 " . $opera_ms_option->{"ILLUMINA_READ_1"} . " " .
		"-2 " . $opera_ms_option->{"ILLUMINA_READ_2"} . " " .
		"-o " . $megahit_out_dir . "> $inter_dir/megahit.out 2> $inter_dir/megahit.err");
	if($?){
	    die "Error during MEGAHIT assembly. Please see $inter_dir/megahit.out and $inter_dir/megahit.err for details.\n";
	}
	#Delete all the meggahit intermidate files
	run_exe("rm -r $megahit_out_dir/intermediate_contigs");
    }
}

sub run_spades{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    my $inter_dir = $opera_ms_option->{"INTER_DIR"};
    my $assembly_out_dir = "$inter_dir/spades_assembly";
    $opera_ms_option->{"CONTIGS_FILE"} = "$assembly_out_dir/contigs.fasta";
    #Stage 0: short assembly using megahit if the contig file is not provided in the config file
    if(! check_completed($opera_ms_option->{"CONTIGS_FILE"}, "Short read assembly [SPAdes]", 1)){
	run_exe("rm -r $assembly_out_dir") if(-d $assembly_out_dir);
	run_exe($opera_ms_dependency->{"spades"}."/spades.py " .
		"--meta " .
		"--threads " . $opera_ms_option->{"NUM_PROCESSOR"} . " " .
		"-1 " . $opera_ms_option->{"ILLUMINA_READ_1"} . " " .
		"-2 " . $opera_ms_option->{"ILLUMINA_READ_2"} . " " .
		"-k 21,33,55,81 " .
		"-o " . $assembly_out_dir . "> $inter_dir/spades.out 2> $inter_dir/spades.err");
	if($?){
	    die "Error during SPAdes assembly. Please see $inter_dir/spades.out and $inter_dir/spades.err for details.\n";
	}
	
	#Delete all the meggahit intermidate files
	run_exe("rm -r $assembly_out_dir/K* $assembly_out_dir/corrected* $assembly_out_dir/misc");
	my $contig_file = $opera_ms_option->{"CONTIGS_FILE"};
	multi_to_linear_fa($contig_file, "$contig_file.single");
	run_exe("mv $contig_file.single $contig_file");
    }
}

sub multi_to_linear_fa{
    my ($multi_file, $single_file) = @_;
    open(FILE, $multi_file);
    open(OUT, ">$single_file");
    my $first = 1;
    while (<FILE>)
    {
	$line = $_;
	chomp $line;
	if ($line =~ />(.*)/) {
	    print OUT "\n" if(! $first);
	    print OUT $line . "\n";
	    $first = 0
	}
	else { print OUT $line; }
    }
    close(OUT);
    close(FILE);
}
    
sub assembly_graph_creation{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    #
    #
    my $inter_dir = $opera_ms_option->{"INTER_DIR"};
    my $opera_ms_dir = $opera_ms_option->{"OPERA_MS_DIR"};
    
    my $cov_dir = "$inter_dir/coverage_estimation";
    $opera_ms_option->{"COV_DIR"} = $cov_dir;
    my $mapping_dir = "$inter_dir/read_mapping";
    $opera_ms_option->{"MAPPING_DIR"} = $mapping_dir;
    
    #my $DIR_COV = "$output_dir/$inter/coverage_estimation/";
    #my $DIR_MAPPING = "$output_dir/$inter/long-read-mapping";

    if($opera_ms_option->{"STAGE"} eq "ALL" || $opera_ms_option->{"STAGE"} eq "ASSEMBLY_GRAPH"){
	$opera_ms_option->{"STAGE"} = "ALL" if($opera_ms_option->{"STAGE_FOLLOW"} == 1);
	
	$start_time = time;
	#Stage 1: Construction of scaffold graph and contig coverage estimation
	
	#
	
	#Short and long read mapping stage
	if(! check_completed("$mapping_dir/config", "Long-read mapping and assembly graph generation", 2)){
	    init_dir($mapping_dir);
	    #print STDERR "\n--- STARTING OPERA-long-read.pl ----\n";
	    $start_time_sub = time;
	    
	    $long_read_mapper_option = " --minimap2 " . $opera_ms_dependency->{"minimap2"} if($opera_ms_option->{"LONG_READ_MAPPER"} eq "minimap2");
	    $long_read_mapper_option = " --blasr " . $opera_ms_dependency->{"blasr"} if($opera_ms_option->{"LONG_READ_MAPPER"} eq "blasr");
	    
	    run_exe("${opera_ms_dir}utils/perl " . 
		    $opera_ms_dependency->{"OPERA-LG"}."/OPERA-long-read.pl" . 
		    " --contig-file " . $opera_ms_option->{"CONTIGS_FILE"} . 
		    " --kmer " . $opera_ms_option->{"KMER_SIZE"} . 
		    " --long-read-file " . $opera_ms_option->{"LONG_READ"} . 
		    " --output-prefix " . "opera" . 
		    " --output-directory " . $mapping_dir . 
		    " --num-of-processors " . $opera_ms_option->{"NUM_PROCESSOR"} . 
		    " --opera " . $opera_ms_dependency->{"OPERA-LG"} .
		    " --illumina-read1 " . $opera_ms_option->{"ILLUMINA_READ_1"} . 
		    " --illumina-read2 " . $opera_ms_option->{"ILLUMINA_READ_2"} . 
		    " --samtools-dir " . $opera_ms_dependency->{"samtools"} . 
		    $long_read_mapper_option . 
		    " --short-read-tooldir " . $opera_ms_dependency->{"bwa"} .
		    " --skip-opera " . "1 " .
		    " --skip-short-read-mapping " . "1 " .
		    "> $mapping_dir/OPERA-long-read.out  2> $mapping_dir/OPERA-long-read.err"
		);
	    #
	    #Error in STAGE 1
	    if($?){
		die "Error in the long read processing. Please see log for details $mapping_dir/OPERA-long-read.out $mapping_dir/OPERA-long-read.err.\n";
	    }
	    $end_time_sub = time;
	    write_time($inter_dir, "graph_construction", ($end_time_sub - $start_time_sub));
	}
	
	
	### run cov_estimate
	if(! check_completed($cov_dir."/".$opera_ms_option{"CONTIGS_COV_FILE"}, "Short-read mapping and coverage estimation", 3)){
	    run_exe("mkdir $cov_dir") if(! -d $cov_dir);
	    #Generate the config file
	    my $cov_estimate_config = "$cov_dir/cov_estimate.config";
	    run_exe("cp " . $opera_ms_option{"CONFIG_PATH"} . " " . $cov_estimate_config);
	    open(CONF, '>>',  $cov_estimate_config);
	    if($opera_ms_option->{"ASSEMBLED_CONTIG_FILE"}){#Add the contig file name to the config file
		print CONF "\nCONTIGS_FILE " . $opera_ms_option->{"CONTIGS_FILE"} ."\n";
	    }
	    if(index($opera_ms_option->{"ILLUMINA_READ_1"}, ",") == -1){#single sample assembly
		print CONF "\nLIB $mapping_dir/opera.bam\n";
		print CONF "\nMAPPING_FILES $mapping_dir/opera.bam\n";
	    }
	    else{
		@illum_read1_tab = split(/,/, $opera_ms_option->{"ILLUM_READ1"});
		my $str_bam = "";
		for(my $i = 0; $i < @illum_read1_tab; $i++){
		    $str_bam .= "$mapping_dir/opera_$i.bam,";
		}
		chop $str_bam;
		#
		print CONF "\nLIB $mapping_dir/opera_1.bam\n";#NEED TO MERGE THE BAM
		#
		#
		print CONF "\nMAPPING_FILES $str_bam\n";
	    }
	    
	    #
	    print CONF "SAMTOOLS_DIR " . $opera_ms_dependency->{"samtools"} . "\n";
	    
	    close(CONF);
	    ##
	    #SIGMA on short reads used to estimate the contig read coverage 
	    #bundling of short reads (can be use if fix the problem with short reads)
	    $start_time_sub = time;
	    run_exe("${opera_ms_dir}utils/perl ${opera_ms_dir}/bin/cov_estimate.pl $cov_estimate_config > $cov_dir/cov_estimate.out 2> $cov_dir/cov_estimate.err");
	    if($?){
		die "Error in during cov_estimate.pl. Please see $cov_dir/cov_estimate.out and $cov_dir/cov_estimate.err for details.\n";
	    }
	    $end_time_sub = time;
	    write_time($inter_dir, "cov_estimation", ($end_time_sub - $start_time_sub));
	    #Remove edge with an outlier edge support/contig coverage ratio => good to remove error due to missmapping
	    #This step must be skipped in case of multi sample assembly
	    #GGG
	    #
	    if(! $opera_ms_option->{"MULTI_SAMPLE"}){
		$start_time_sub = time;
		run_exe("${opera_ms_dir}utils/perl ${opera_ms_dir}bin/filter_edge_coverage_ratio.pl $mapping_dir " . "$cov_dir/".$opera_ms_option{"CONTIGS_COV_FILE"} . " 2> $mapping_dir/filter_edge_coverage_ratio.err");
		write_time($inter_dir, "cov_filtering", ($end_time_sub - $start_time_sub));
		$end_time_sub = time;
		if($?){
		    die "Error in during filter_edge_coverage_ratio.pl. Please see $mapping_dir/filter_edge_coverage_ratio.err for details.\n";
		}
	    }
	    else{
		#create soft links
	    }
	}
    }
}

sub hierarchical_clustering{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    #
    #
    my $inter_dir = $opera_ms_option->{"INTER_DIR"};
    my $opera_ms_dir = $opera_ms_option->{"OPERA_MS_DIR"};
    my $cov_dir = $opera_ms_option->{"COV_DIR"};
    
    my $sigma_dir = "$inter_dir/hierarchical_clustering";
    $opera_ms_option->{"SIGMA_DIR"} = $sigma_dir;
    

    if($opera_ms_option->{"STAGE"} eq "ALL" || $opera_ms_option->{"STAGE"} eq "CLUSTERING"){
	$opera_ms_option->{"STAGE"} = "ALL" if($opera_ms_option->{"STAGE_FOLLOW"} == 1);
	if(! check_completed("$sigma_dir/r_value_evaluation_summary.dat", "Hierarchical clustering", 4)){
	    init_dir($sigma_dir);
	    $start_time = time;
	    ### Run sigma for long read (filtered edge obtain using filter_edge_coverage_ratio.pl###
	    #
	    generate_sigma_config_file($cov_dir, $opera_ms_option->{"MAPPING_DIR"}, $sigma_dir);
	    
	    #Run sigma
	    run_exe("${opera_ms_dir}bin/sigma $sigma_dir/sigma.config > $sigma_dir/sigma.out 2> $sigma_dir/sigma.err");
	    
	    if($?){
		die "Error during ${opera_ms_dir}bin/sigma execution. Please see $sigma_dir/sigma.out and $sigma_dir/sigma.err for details.\n";
	    }
	    
	    #Refine the sigma R parameter
	    run_exe("${opera_ms_dir}utils/perl ${opera_ms_dir}bin/refine_r_estimate.pl " .
		    $cov_dir."/".$opera_ms_option->{"CONTIGS_COV_FILE"} . " " . 
		    $sigma_dir . " " .
		    $opera_ms_option->{"NUM_PROCESSOR"} . " " .
		    $opera_ms_dir . " 2> $sigma_dir/refine_r.err"
		);
	    
	    $end_time = time;
	    write_time($inter_dir, "sigma", ($end_time - $start_time));
	}
    }
}

sub reference_clustering{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    #
    #
    my $inter_dir = $opera_ms_option->{"INTER_DIR"};
    my $opera_ms_dir = $opera_ms_option->{"OPERA_MS_DIR"};
    
    my $ref_clustering_dir = "$inter_dir/reference_clustering";
    $opera_ms_option->{"REF_CLUSTERING_DIR"} = $ref_clustering_dir;
    
    if($opera_ms_option->{"STAGE"} eq "ALL" || $opera_ms_option->{"STAGE"} eq "REF_CLUSTERING"){
	$opera_ms_option->{"STAGE"} = "ALL" if($opera_ms_option->{"STAGE_FOLLOW"} == 1);
	
	if(! check_completed("$ref_clustering_dir/cluster_species.dat", "Reference based clustering", 5, !$opera_ms_option->{"REF_CLUSTERING"})){
	    init_dir($ref_clustering_dir);
	    $start_time = time;
	    #Run reference mapping to merge clusters, use both long read and short read.
	    #Use the short read and add a link if support > 5, gap estimated < 500 and validated on a reference
	    #NEED TO REMOVE LINKS SUPPORTED BY BOTH SHORT AND LONG READS
	    
	    
	    
	    run_exe("${opera_ms_dir}utils/perl ${opera_ms_dir}bin/sequence_similarity_clustering.pl " . 
		    $inter_dir . " " .
		    $ref_clustering_dir . " " .
		    $opera_ms_option->{"MAPPING_DIR"} . " " .
		    $opera_ms_option->{"SIGMA_DIR"} . " " .
		    $opera_ms_option->{"CONTIGS_FILE"} . " " .
		    $opera_ms_option->{"NUM_PROCESSOR"} . " " .
		    $opera_ms_option->{"MASH_DB"} . " " .
		    $opera_ms_option->{"KRAKEN_DB"} . " " .
		    $opera_ms_dir . " " . 
		    $opera_ms_dependency->{"mash"} . " " .
		    $opera_ms_dependency->{"mummer"} . " " . 
		    $opera_ms_dependency->{"kraken"} . " 2> $ref_clustering_dir/sequence_similarity_clustering.err");
	    if($?){
		die "Error in during reference based clustering. Please see $ref_clustering_dir/sequence_similarity_clustering.err for details.\n";
	    }
	    $end_time = time;
	    write_time($inter_dir, "ref_clustering", ($end_time - $start_time));
	}
    }
}


sub strain_clustering_and_assembly{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    #
    #
    my $inter_dir = $opera_ms_option->{"INTER_DIR"};
    my $opera_ms_dir = $opera_ms_option->{"OPERA_MS_DIR"};

    my $sigma_dir = $opera_ms_option->{"SIGMA_DIR"};
    my $ref_clustering_dir = $opera_ms_option->{"REF_CLUSTERING_DIR"};
    my $strain_dir = "$inter_dir/strain_analysis";
    $opera_ms_option->{"STRAIN_DIR"} = $strain_dir;
    #
    my $opera_lr_dir = "$inter_dir/opera_long_read";
    $opera_ms_option->{"OPERA_LR_DIR"} = $opera_lr_dir;
    
    if(($opera_ms_option->{"STAGE"} eq "ALL" || $opera_ms_option->{"STAGE"} eq "STRAIN_CLUSTERING")){
	$opera_ms_option->{"STAGE"} = "ALL" if($opera_ms_option->{"STAGE_FOLLOW"} == 1);
	if( ! check_completed("$strain_dir/done", "Strain clustering and assembly", 6, !($opera_ms_option->{"REF_CLUSTERING"} && $opera_ms_option->{"STRAIN_CLUSTERING"}))){#Add check for reference clustering
	    init_dir($strain_dir);
	    
	    $start_time = time;
	    #Plug-in the strain level assembly
	    #->if a species contains the strain each strain will be assembled independently by opera-lg
	    #->if a species does not contain the strain they are all pulled together and assembled in a single run by opera-lg
	    
	    run_exe("${opera_ms_dir}utils/perl ${opera_ms_dir}bin/coverage_clustering.pl " . 
		    $opera_ms_option->{"INTER_DIR"} . " " .
		    $strain_dir . " " .
		    $opera_ms_option->{"CONTIGS_FILE"} . " " .
		    $ref_clustering_dir . " "  .
		    $opera_ms_dir . " " .
		    "1  " . "2>  $strain_dir/coverage_clustering.err");
	    if($?){
		die "Error in during strain clustering. Please see $strain_dir/coverage_clustering.err for details.\n";
	    }
	    
	    #Assembly of species with multiple strains
	    my $global_r_value = `tail -n1 $sigma_dir/r_value_evaluation_summary.dat | cut -f1`;chop $global_r_value;
	    #
	    my $r_value_interval = `head -n1 $sigma_dir/r_estimate_value.dat`;chop $r_value_interval;
	    my @r_value_interval_tab = split(/\s+/, $r_value_interval);
	    my $r_value_step = ($r_value_interval_tab[1] - $r_value_interval_tab[0]) / 10;
	    #
	    open(FILE, "$strain_dir/reference_length.dat");
	    my $cmd_strain = "$strain_dir/assembly.cmd";
	    open(CMD_STRAIN, ">$cmd_strain");
	    my ($species, $species_length);
	    while(<FILE>){
		@line = split(/\t/, $_);
		$species = $line[0];
		$species_length = $line[1];
		print CMD_STRAIN "${opera_ms_dir}utils/perl ${opera_ms_dir}bin/cluster_strain.pl $strain_dir/$species $species_length $global_r_value $r_value_step $opera_ms_dir 2> $strain_dir/$species.log" . "\n";
	    }
	    close(FILE);
	    close(CMD_STRAIN);
	    my $num_processor = $opera_ms_option->{"NUM_PROCESSOR"};
	    run_exe("cat $cmd_strain | xargs -L 1 -P $num_processor -I COMMAND sh -c \"COMMAND\"");
	    if($?){
		die "Error in during strain clustering assembly. Please see $strain_dir/*.log for details.\n";
	    }
	    $end_time = time;
	    run_exe("touch $strain_dir/done");
	    write_time($inter_dir, "strain_clustering", ($end_time - $start_time));
	}

	########## SMALL SPECIES CLUSTER AND NO SPECIES CLUSTER ASSEMBLY
	if(! check_completed("$opera_lr_dir/statistics", "Assembly of other clusters",7)){
	    init_dir($opera_lr_dir);
	    
	    $start_time = time;
	    
	    #Add the contigs in cluster strain 0 back to the opera assembly
	    if($opera_ms_option->{"REF_CLUSTERING"} && $opera_ms_option->{"STRAIN_CLUSTERING"}){
		run_exe("rm $opera_lr_dir/*single*") if(-e "$opera_lr_dir/contigs_single");#Deletion required to avoid re-using previous assemblies for the gapfilling
		run_exe("cat $ref_clustering_dir/single_strain_species.fa $strain_dir/*/excluded_contigs.fa > $ref_clustering_dir/single_strain_species_with_excluded.fa");
	    }
		
	    #Assembly of all the other contigs were no multiple strain of the same species have inferred
	    #Filter that remove contigs with a coverage 1.5 times higher than the mean and contig that are defined as repeat based on the reference mapping
	    #remove the contigs that belong to clusters that are associated to species with multiple strains
	    my $contig_mapping_to_reference = "NULL";#No repeat detection based on the reference genome due to high computational time ADD AN OPTION TO TOGGLE IT ON/OFF
	    my $clustering_edge_dir =  $opera_ms_option->{"SIGMA_DIR"};
	    my $reference_cluster_file  = "$clustering_edge_dir/clusters";
	    my $contig_file = $opera_ms_option->{"CONTIGS_FILE"};
	    #
	    if($opera_ms_option->{"REF_CLUSTERING"}){
		$contig_mapping_to_reference = "$ref_clustering_dir/NUCMER_OUT/";#
		$clustering_edge_dir = $ref_clustering_dir;
		$reference_cluster_file = "$ref_clustering_dir/clusters_seq_similarity";
		if($opera_ms_option->{"STRAIN_CLUSTERING"}){
		    $contig_file = "$ref_clustering_dir/single_strain_species_with_excluded.fa";
		    $reference_cluster_file = "$ref_clustering_dir/clusters_single_strain"
		}
	    }
	    
	    $contig_coverage_file = $opera_ms_option->{"COV_DIR"} . "/" . $opera_ms_option->{"CONTIGS_COV_FILE"};
	    
	    run_exe("${opera_ms_dir}utils/perl ${opera_ms_dir}bin/filter_cluster_coverage.pl $contig_coverage_file $reference_cluster_file $clustering_edge_dir 1.5 $clustering_edge_dir/NO_REPEAT $contig_mapping_to_reference " . $opera_ms_dependency->{"mummer"} . "2> $clustering_edge_dir/filter_cluster_coverage.err");
	    if($?){
		die "Error in during filter_cluster_coverage. Please see $clustering_edge_dir/filter_cluster_coverage.err for details.\n";
	    }

	    #Generate the OPERA-long-read config file
	    
	    generate_opera_config_file($opera_ms_option->{"MAPPING_DIR"}, $clustering_edge_dir,  $opera_lr_dir, $contig_file);
	    
	    #finally, run opera
	    my $opera_lg_exe_dir = $opera_ms_dependency->{"OPERA-LG"};
	    run_exe("timeout 5m $opera_lg_exe_dir/OPERA-LG $opera_lr_dir/opera.config > $opera_lr_dir/opera_lg.out 2> $opera_lr_dir/opera_lg.err");
	    if(! -e "$opera_lr_dir/scaffoldSeq.fasta"){#NEED TO IN THE MAKE FILE THE COMMAND TO GET THE FASTS OPERA-MS
		run_exe("$opera_lg_exe_dir/OPERA-LG-fast $opera_lr_dir/opera.config > $opera_lr_dir/opera_lg_fast.out 2> $opera_lr_dir/opera_lg_fast.err");
	    }
	    if($?){
		die "Error in during OPERA-LG assembly. Please see $opera_lr_dir/opera_lg*.out and $opera_lr_dir/opera_lg*.err for details.\n";
	    }
	    $end_time = time;
	    write_time($inter_dir, "other_assembly", ($end_time - $start_time));
	}
    }
}

sub gap_filling{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    #
    #
    my $final_output_dir = $opera_ms_option->{"OUTPUT_DIR"};
    my $inter_dir = $opera_ms_option->{"INTER_DIR"};
    my $opera_ms_dir = $opera_ms_option->{"OPERA_MS_DIR"};
    my $opera_lg_dir = $opera_ms_dependency->{"OPERA-LG"};
    
    my $strain_dir = $opera_ms_option->{"STRAIN_DIR"};
    my $opera_lr_dir = $opera_ms_option->{"OPERA_LR_DIR"};
    my $gap_filling_dir = "$opera_lr_dir/GAPFILLING";
    $opera_ms_option->{"GAPFILLING_DIR"} = $gap_filling_dir;

    
    
    if($opera_ms_option->{"STAGE"} eq "ALL" || $opera_ms_option->{"STAGE"} eq "GAP_FILLING"){
	$opera_ms_option->{"STAGE"} = "ALL" if($opera_ms_option->{"STAGE_FOLLOW"} == 1);
	
	if(! check_completed("$opera_lr_dir/scaffoldSeq.fasta.filled", "Gap filling", 8)){
	    init_dir($gap_filling_dir);
	    
	    $start_time = time;
	    $start_time_sub = time;
	    
	    #Add the scaffold obtain during the strain level assembly back in
	    my $opera_scaff_file = "$opera_lr_dir/scaffolds.scaf";
	    if($opera_ms_option->{"REF_CLUSTERING"} && $opera_ms_option->{"STRAIN_CLUSTERING"}){
		my $opera_scaff_file_single = "$opera_lr_dir/scaffolds_single.scaf";
		run_exe("mv $opera_scaff_file $opera_scaff_file_single") if(! -e $opera_scaff_file_single);
		run_exe("cp $opera_scaff_file_single $opera_scaff_file");
	    }
    
	    my $opera_scaff_seq_file = "$opera_lr_dir/scaffoldSeq.fasta";
	    if($opera_ms_option->{"REF_CLUSTERING"} && $opera_ms_option->{"STRAIN_CLUSTERING"}){
		my $opera_scaff_seq_file_single = "$opera_lr_dir/scaffoldSeq_single.fasta";
		run_exe("mv $opera_scaff_seq_file $opera_scaff_seq_file_single") if(! -e $opera_scaff_seq_file_single);
		run_exe("cp $opera_scaff_seq_file_single $opera_scaff_seq_file");
	    }

	    my $opera_contig_size_file = "$opera_lr_dir/contigs";
	    if($opera_ms_option->{"REF_CLUSTERING"} && $opera_ms_option->{"STRAIN_CLUSTERING"}){
		my $opera_contig_size_file_single = "$opera_lr_dir/contigs_single";
		run_exe("mv $opera_contig_size_file $opera_contig_size_file_single") if(! -e $opera_contig_size_file_single);
		run_exe("cp $opera_contig_size_file_single $opera_contig_size_file");
	    }
    
	    #Add the multiple strain genome scaffold to the pool
	    if($opera_ms_option->{"REF_CLUSTERING"} && $opera_ms_option->{"STRAIN_CLUSTERING"}){
		opendir(DIR, "$strain_dir");
		my @strain_dir = readdir(DIR);
		close(DIR);
	
		my $strain_ID = 1;
		my $scaff_strain_id = 1;
		my ($current_strain_dir, $strain_scaff_file, $strain_scaff_seq_file, $strain_contig_size_file);
		
		foreach $strain (@strain_dir){
		    $strain_ID = 1;
		    $current_strain_dir = "$strain_dir/$strain/STRAIN_$strain_ID/";
		    $strain_scaff_file = "$current_strain_dir/scaffolds.scaf";
		    while( -e $strain_scaff_file){
			#
			$strain_scaff_seq_file = "$current_strain_dir/scaffoldSeq.fasta";
			$strain_contig_size_file = "$current_strain_dir/contigs";
		#
			run_exe("sed 's/>opera/>strain$scaff_strain_id\_opera/' $strain_scaff_file >> $opera_scaff_file");
			run_exe("sed 's/>opera/>strain$scaff_strain_id\_opera/' $strain_scaff_seq_file >> $opera_scaff_seq_file");
			run_exe("grep -v Length $strain_contig_size_file >> $opera_contig_size_file");
			#
			$strain_ID++;
			$scaff_strain_id++;
			$current_strain_dir = "$strain_dir/$strain/STRAIN_$strain_ID/";
			$strain_scaff_file = "$current_strain_dir/scaffolds.scaf";
		    }
		}
	    }

	    $end_time_sub = time;
    
	    #perform the gap fill
	    if($opera_ms_option->{"REF_CLUSTERING"}){
		#run_exe("${opera_ms_dir}utils/perl $opera_ms_dir/bin/match_scaffolds_clusters.pl " . 
		#	$opera_ms_option->{"REF_CLUSTERING_DIR"} . " " .
		#	$opera_ms_option->{"NUM_PROCESSOR"} . " " .
		#	"$opera_lr_dir/scaffoldSeq.fasta" . " " .
		#	"$opera_lr_dir/scaffolds.scaf" . " " .
		#	$opera_ms_dependency->{"mummer"} . " 2> $opera_lr_dir/s_info.err");

		if($?){
		    die "Error in during bin/match_scaffolds_clusters.pl. Please see $opera_lr_dir/s_info.err for details.\n";
		}
	    }

	    #print STDERR " ------ " . $opera_ms_option->{"GAP_FILLING"} . "\n";
	    if($opera_ms_option->{"GAP_FILLING"}){
		
		run_exe("${opera_ms_dir}utils/perl $opera_lg_dir/gapfilling.pl " .
			$gap_filling_dir . " " .
			$opera_ms_option->{"CONTIGS_FILE"} . " " .
			$opera_ms_option->{"LONG_READ"} . " " .
			$opera_ms_option->{"MAPPING_DIR"} . " " .
			$opera_ms_option->{"NUM_PROCESSOR"} . " " .
			$opera_ms_option->{"FILLED_SCAFF_LENGTH"} . " " .
			"$opera_lg_dir" . " " .
			$opera_ms_dependency->{"racon"} . " " .
			$opera_ms_dependency->{"minimap2"} . " " .
			$opera_ms_dependency->{"mummer"} . " " . 
			"2> $gap_filling_dir/gap_filling.err"
		    );
		if($?){
		    die "Error in during gap_filling. Please see $gap_filling_dir/gap_filling.err for details.\n";
		}
	    }
	    else{
		run_exe("ln -s $opera_lr_dir/scaffoldSeq.fasta $opera_lr_dir/scaffoldSeq.fasta.filled");
	    }
	    
	    $end_time = time;
	    write_time($inter_dir, "gap_filling", ($end_time - $start_time));

	    #NEED TO TEST ON THAT
	    #Remove the contigs that have been used during gapfiling in the .fa file
	    #$command = "perl $opera_ms_dir/bin/clean_scaffolds_repeats.pl $output_dir/$inter";
	    #run_exe($command);
	    #$command = "mv $output_dir/lib_* $output_dir/$inter; rm -r $output_dir/$inter/scaffolds; mv $output_dir/runOperaMS.config $output_dir/$inter/";
	    #run_exe($command);
	    
	}
    }
}

sub polishing{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    #
    #
    my $final_output_dir = $opera_ms_option->{"OUTPUT_DIR"};
    my $inter_dir = $opera_ms_option->{"INTER_DIR"};
    my $polishing_dir = "$inter_dir/polished_assembly";
    my $contig = "$final_output_dir/contigs.fasta";
    my $polished_contig = "$final_output_dir/contigs.polished.fasta";
    my $num_processor = $opera_ms_option->{"NUM_PROCESSOR"};
    
    if($opera_ms_option->{"POLISHING"} || $opera_ms_option->{"STAGE"} eq "POLISHING"){
	$opera_ms_option->{"STAGE"} = "ALL" if($opera_ms_option->{"STAGE_FOLLOW"} == 1);
	
	if(! check_completed("$polished_contig", "Polishing", 9)){
	    init_dir($polishing_dir);
	    
	    #get the bwa mapping
	    my $read_1 = $opera_ms_option->{"ILLUMINA_READ_1"};
	    my $read_2 = $opera_ms_option->{"ILLUMINA_READ_2"};
	    my $temp_contig = "$polishing_dir/contigs.fa";
	    run_exe("ln -s $contig $temp_contig");
	    my $bwa_err = "$polishing_dir/bwa_index.err";
	    my $bwa_out = "$polishing_dir/bwa_index.out";
	    #Refence index
	    run_exe($opera_ms_dependency->{"bwa"} . "/bwa index $temp_contig > $bwa_out 2> $bwa_err");
	    if($?){
		die "Error in during bwa indexing. Please see $bwa_out and $bwa_err for details.\n";
	    }
	    #Get a sorted bam file
	    $bwa_err = "$polishing_dir/bwa.err";
	    $bwa_out = "$polishing_dir/bwa.out";
	    run_exe($opera_ms_dependency->{"bwa"} . "/bwa mem -t $num_processor $temp_contig  $read_1 $read_2 2>> $bwa_err| " . 
		    $opera_ms_dependency->{"samtools"} ."/samtools view -Sub - 2>> $bwa_err | " . 
		    $opera_ms_dependency->{"samtools"} . "/samtools sort - $temp_contig > $bwa_out 2>> $bwa_err");
	    if($?){
		die "Error in during bwa mapping. Please see $bwa_out and $bwa_err for details.\n";
	    }
	    #Index the bam file
	    run_exe($opera_ms_dependency->{"samtools"} ."/samtools index $temp_contig.bam >> $bwa_out 2>> $bwa_err");
	    if($?){
		die "Error in during samtools indexing. Please see $bwa_out and $bwa_err for details.\n";
	    }
	    
	    #run pilon [in multi thread]
	    $pilon_path = $opera_ms_dependency->{"pilon"}."/pilon.jar";
	    run_exe("java -jar $pilon_path --fix bases --threads $num_processor --genome $temp_contig --bam $temp_contig.bam --outdir $polishing_dir > $polishing_dir/pilon.out 2> $polishing_dir/pilon.err");
	    if($?){
		die "Error in during pilon. Please see $polishing_dir/pilon.out and $polishing_dir/pilon.err for details.\n";
	    }
	    #Move the final output
	    run_exe("sed 's/_pilon//'  $polishing_dir/pilon.fasta > $polished_contig");
	    run_exe("rm $polishing_dir/pilon.fasta");
	}
	#Generate the assembly stats and write the final assembly after polishing
	write_cluster_assembly(\%opera_ms_option, \%opera_ms_dependency);
    }
   
}

sub generate_assembly_stats{
    my ($opera_ms_option, $opera_ms_dependency) = @_;

    my $final_output_dir = $opera_ms_option->{"OUTPUT_DIR"};
    my $ref_clustering_dir = $opera_ms_option->{"REF_CLUSTERING_DIR"};
    my $opera_ms_dir = $opera_ms_option{"OPERA_MS_DIR"};
    my $opera_lr_dir = $opera_ms_option->{"OPERA_LR_DIR"};
    
    my $cluster_file =  $opera_ms_option->{"SIGMA_DIR"}."/clusters";
    my $species_file = "NULL";
    my $cov_dir = $opera_ms_option->{"COV_DIR"};
    my $coverage_file = "$cov_dir/" . $opera_ms_option{"CONTIGS_COV_FILE"};
    my $read_size = `cat $cov_dir/read_size.dat`;chop $read_size;
    my $long_read_coverage_file = $opera_ms_option->{"MAPPING_DIR"}."/opera.map.cov";
    my $strain_dir = "NULL";
    if($opera_ms_option->{"REF_CLUSTERING"}){
	my $ref_clustering_dir = $opera_ms_option->{"REF_CLUSTERING_DIR"};
	$cluster_file =  "$ref_clustering_dir/clusters_seq_similarity";
	$species_file =  "$ref_clustering_dir/cluster_species.dat";
	$strain_dir = $opera_ms_option->{"STRAIN_DIR"};
    }
    
    $start_time = time;
    $opera_ms_option->{"CONTIG_INFO"} = "$final_output_dir/contig_info.txt";
    
    run_exe("${opera_ms_dir}utils/perl $opera_ms_dir/bin/generate_assembly_stats.pl " .
	    "$opera_lr_dir/scaffolds.scaf" . " " .
	    $cluster_file . " " .
	    $read_size . " " .
	    $coverage_file . " " .
	    $long_read_coverage_file . " " .
	    $species_file . " " .
	    $strain_dir . " " .
	    "$opera_lr_dir/scaffoldSeq.fasta.filled" . " " .
	    $opera_ms_option->{"CONTIG_INFO"} . " " . 
	    $opera_ms_dir . " " .
	    $opera_ms_option->{"NUM_PROCESSOR"} . " " . 
	    $opera_ms_dependency->{"mummer"} . " 2> $opera_lr_dir/g_assembly_stats.err");
    if($?){
	die "Error in during assembly statistics generation. Please see $ref_clustering_dir/s_info.err for details.\n";
    }
    
    #Get the stats
    my @all_size = ();
    open(FILE, "$final_output_dir/contig_info.txt");
    my $contig_1mb = 0;
    my $contig_100kb = 0;
    my $contig_500kb = 0;
    my $assembly_size = 0;
    <FILE>;
    while(<FILE>){
	@line = split(/\t/, $_);
	$size = $line[1];
	$contig_1mb++ if($size >= 1000000);
	$contig_500kb++ if($size >= 500000);
	$contig_100kb++ if($size >= 100000);
	$assembly_size += $size;
	push(@all_size, $size);
    }
    close(FILE);
    
    my @sort_tab = sort {$b <=> $a} @all_size;
    my $nb_contig = @all_size+0;
    my $nb_seq = compute_Nx(50, $assembly_size, \@sort_tab);$n50 = $sort_tab[$nb_seq];
    $time = localtime;
    my $str_stats =  "\n[$time]\tAssembly stats\n" . 
	"Number of contigs: " . $nb_contig . " \n". 
	"Assembly size: $assembly_size bp\n". 
	#"\t" . "min contig size $sort_tab[$nb_contig-1] bp, ".
	"Max contig size: $sort_tab[0] bp\n".
	"Contig(s) longer than 1Mbp: $contig_1mb \n".
	"Contig(s) longer than 500kbp: $contig_500kb\n".
	"Contig(s) longer than 100kbp: $contig_100kb\n".
	"Contig N50: $n50 bp\n";
    
    print STDOUT $str_stats . "\n";
    open(OUT, ">$final_output_dir/assembly.stats");
    print OUT $str_stats . "\n";
    close(OUT);
    
    $end_time = time;
    write_time($opera_ms_option->{"INTER_DIR"}, "stats_generation", ($end_time - $start_time));
}

sub write_final_assembly{
    my ($opera_ms_option, $opera_ms_dependency) = @_;

    my $final_output_dir = $opera_ms_option->{"OUTPUT_DIR"};
    my $opera_lr_dir = $opera_ms_option->{"OPERA_LR_DIR"};
    
    #Write the final contig assembly
    #Get the info from each contig
    my %contig_info = ();
    open(FILE, $opera_ms_option->{"CONTIG_INFO"});
    my ($name, $length, $cov_short, $cov_long, $cluster, $species, $nb_strain, $ref);
    while(<FILE>){
	($name, $length, $cov_short, $cov_long, $cluster, $species, $nb_strain, $ref) = split(/\t/, $_);
	$contig_info{$name} = 
	    "length: $length " .
	    "cov_short: $cov_short " .
	    "cov_long: $cov_long " .
	    "cluster: $cluster ";
    }
    close(FILE);

    #Get the scaff to contig name id
    my %scaff_to_contig = ();
    open(FILE, "$opera_lr_dir/scaffolds.scaf.cname");
    my ($contig_name, $scaff_name);
    while(<FILE>){
	chop $_;
	($scaff_name, $contig_name) = split(/\t/, $_);
	$scaff_to_contig{$scaff_name} = $contig_name;
    }

    #Write the final contig sequence
    open(FILE, "$opera_lr_dir/scaffoldSeq.fasta.filled");
    open(OUT, ">$final_output_dir/contigs.fasta");
    while(<FILE>){
	@line = split(/\s+/, $_);
	$scaff_name = substr($line[0], 1);
	$contig_name = $scaff_to_contig{$scaff_name};
	print OUT ">$contig_name " . $contig_info{$contig_name} . "\n";
	$seq = <FILE>;
	print OUT $seq;
    }
    close(FILE);
    close(OUT);
    
    #Write the cluster assembly
    write_cluster_assembly($opera_ms_option, $opera_ms_dependency);
    
}

sub write_cluster_assembly{
    my ($opera_ms_option, $opera_ms_dependency) = @_;
    
    my $final_output_dir = $opera_ms_option->{"OUTPUT_DIR"};
    
    if($opera_ms_option->{"REF_CLUSTERING"}){
	my $cluster_assembly_dir = "$final_output_dir/opera_ms_clusters";
	init_dir($cluster_assembly_dir);
	$final_contig = "$final_output_dir/contigs.polished.fasta";
	$final_contig = "$final_output_dir/contigs.fasta" if(! -e $final_contig);
	get_cluster_sequence($opera_ms_option->{"CONTIG_INFO"}, $cluster_assembly_dir, $final_contig);
	generate_cluster_stats($opera_ms_option, $opera_ms_dependency, $cluster_assembly_dir);
    }
}

sub get_cluster_sequence{
    my ($info_file, $cluster_dir, $contig_file) = @_;
    
    my %cluster_file = ();
    my %contig_cluster = ();
    open(FILE, $info_file);
    <FILE>;
    while(<FILE>){
	@line = split(/\t/, $_);
	$contig_id = $line[0];
	$cluster_id = $line[4];
	
	if( ! exists $cluster_file{$cluster_id}){
	    #print STDERR " *** Writing $strain_dir/$species\_strain\_$strain_counter.fa => $strain_id\n";
	    open($cluster_file{$cluster_id}, ">$cluster_dir/$cluster_id.fasta");
	}
	$contig_cluster{">".$contig_id} = $cluster_id;
    }
    close(FILE);
    #exit(0);
    open(IN, $contig_file);
    my $print_flag = 0;

    my $OUT_FILE = "";
    my $ID;
    while(<IN>){
	
	if($_ =~ />.*/){
	    #print STDERR "---> $nb_contig\n" if($nb_contig % 500000 == 0);$nb_contig++;
	    @line = split(/\s+/, $_);
	    $ID = $line[0];
	    
	    #print STDERR $ID."\n"; <STDIN>;
	    if(exists $contig_cluster{$ID}){
		$print_flag = 1;
		$OUT_FILE = $cluster_file{$contig_cluster{$ID}};
		if(! defined $OUT_FILE){
		    $print_flag = 0;
		}
		else{
		    #print STDERR " *** $ID => $contig_strain{$ID}\n";
		    #chop $_;
		    #print $_."\n";
		    print $OUT_FILE $ID . "\n";
		}
	    }
	    else{
		$print_flag = 0;
	    }
	}
	else{
	    print $OUT_FILE $_ if($print_flag);
	}
    }
    close(IN);

    foreach $s (keys %strain_file){
	close($strain_file{$s});
    }
    run_exe("mv $cluster_dir/NA.fasta $cluster_dir/unclustered.fasta");
}

sub generate_cluster_stats{
    my ($opera_ms_option, $opera_ms_dependency, $cluster_dir) = @_;

    #Run mash if not runned
    my $mash_dist_file = "$cluster_dir/mash.dist";
    run_exe($opera_ms_dependency->{"mash"} ."/mash". " dist -p 1 -d 0.2 " . $opera_ms_option{"MASH_DB"} . " " .  "$cluster_dir/*fa* > $mash_dist_file") if(! -e $mash_dist_file);
    
    #Get the statistics
    my $opera_ms_dir = $opera_ms_option->{"OPERA_MS_DIR"};
    run_exe("${opera_ms_dir}utils/perl $opera_ms_dir/bin/cluster_info.pl " . $opera_ms_option{"OUTPUT_DIR"} . " $cluster_dir > $cluster_dir/../cluster_info.txt");
}

sub check_completed{
    my ($file, $message, $stage_number, $skip) = @_;
    my $nb_total_stage = $opera_ms_option{"NB_TOTAL_STAGE"};
    $res = 0;
    $time = localtime;
    print STDOUT "\n[$time]\t$message [$stage_number/$nb_total_stage]\n";
    if(defined $skip && $skip){
	print STDOUT "[$time]\tSkip \n";
	$res = 1;
    }
    else{
	if(!$opera_ms_option{"RESTART"} && -e $file){
	    print STDOUT "[$time]\tStage previously completed\n";# [$file detected]\n";
	    $res = 1;
	}
	else{
	    $opera_ms_option{"RESTART"} = 1;
	}
    }
    return $res;
}

sub init_dir{
    my ($dir) = @_;
    run_exe("rm -r $dir")if(-d $dir);
    run_exe("mkdir -p $dir");
}

sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return . "\n" if($run);
    return $return;
}


sub read_argument{

    my ($opera_ms_option, $read_option) = @_;
    my ( $illum_read1, $illum_read2, $long_read, $output_dir, $contig_file,  $strain_cluster, $reference_cluster, $contig_len, $contig_edge_len, $contig_win_len, $kmer_size, $nproc, $long_read_mapping, $flag_help);

    set_default_value($opera_ms_option);
    
    my $help_message = print_help();

    
    #$p = Getopt::Long::Parser->new;
    GetOptionsFromArray(
	$read_option,
	"contig-file=s"       => \$opera_ms_option{"CONTIGS_FILE"},
	"short-read1=s"    => \$opera_ms_option{"ILLUMINA_READ_1"},
	"short-read2=s"    => \$opera_ms_option{"ILLUMINA_READ_2"},
	"long-read=s"    => \$opera_ms_option{"LONG_READ"},
	"out-dir=s"  => \$opera_ms_option{"OUTPUT_DIR"},
	#
	#
	"strain-clustering!"            => \$opera_ms_option{"STRAIN_CLUSTERING"},
	"ref-clustering!"         => \$opera_ms_option{"REF_CLUSTERING"},
	"gap-filling!"         => \$opera_ms_option{"GAP_FILLING"},
	"polishing!"             => \$opera_ms_option{"POLISHING"},
	#
	"contig-len-thr=i"    => \$opera_ms_option{"CONTIG_LEN_THR"},
	"contig-edge-len=i"   => \$opera_ms_option{"CONTIG_EDGE_LEN"},	
	"contig-window-len=i" => \$opera_ms_option{"CONTIG_WINDOW_LEN"},
	"kmer=i"              => \$opera_ms_option{"KMER_SIZE"},	
	#
	#
	
	"long-read-mapper=s" => \$opera_ms_option{"LONG_READ_MAPPER"},
	"short-read-assembler=s" => \$opera_ms_option{"SHORT_READ_ASSEMBLER"},
	"num-processors=i" => \$opera_ms_option{"NUM_PROCESSOR"},
	#	
	"help"                => \$flag_help,
	) or die("Error in command line arguments.\n$help_message");
    
    
    if($flag_help){
	print STDERR $help_message;exit(0);
    }
    
    else{
	#Write config file
	$output_dir = $opera_ms_option{"OUTPUT_DIR"};
	my $config_file = "$output_dir/opera-ms.config";
	run_exe("mkdir -p $output_dir");
	
	$opera_ms_option->{"CONFIG_PATH"} = $config_file;
	open(OUT, ">$config_file");

	my @option_order = (
	    "ILLUMINA_READ_1",
	    "ILLUMINA_READ_2",
	    "LONG_READ",
	    "OUTPUT_DIR",
	    
	    #"POLISHING",
	    #"REF_CLUSTERING",
	    #"STRAIN_CLUSTERING",
	    
	    "CONTIG_EDGE_LEN",
	    "CONTIG_WINDOW_LEN",
	    "CONTIG_LEN_THR",
	    "KMER_SIZE",
	    
	    "CONTIGS_FILE",
	    "LONG_READ_MAPPER",
	    "SHORT_READ_ASSEMBLER",
	    "NUM_PROCESSOR",
	    );
	
	
	foreach $option (@option_order){
	    print OUT $option . " " . $opera_ms_option->{$option} . "\n" if(defined $opera_ms_option->{$option});
	}

	#print STDERR "REF_CLUSTERING" . " " . $opera_ms_option->{"REF_CLUSTERING"} ."\n" . "STRAIN_CLUSTERING" . " "  .$opera_ms_option->{"STRAIN_CLUSTERING"} . "\n POLISHING " . $opera_ms_option{"POLISHING"} . "\n";exit(0);
	
	$opera_ms_option->{"CONTIGS_COV_FILE"} = "contigs_" . $opera_ms_option->{"CONTIG_WINDOW_LEN"} . "_" . $opera_ms_option->{"CONTIG_EDGE_LEN"};
	$opera_ms_option->{"RESTART"} = 0;
	$opera_ms_option->{"STAGE"} = "ALL";
	#To compute the total number of stage
	$opera_ms_option->{"NB_TOTAL_STAGE"} = 8;
	$opera_ms_option->{"NB_TOTAL_STAGE"}++ if($opera_ms_option->{"POLISHING"});
	
	close(OUT);
    }
}

sub print_help{
    "OPERA-MS.pl: OPERA-MS " .	$opera_ms_option{"VERSION"} . "
contacts: Denis Bertrand <bertrandd\@gis.a-star.edu.sg>
          Chengxuan Tong <Tong_Chengxuan\@gis.a-star.edu.sg>

Usage:
  perl OPERA-MS.pl [options] --illumina-read1 <pe1> --illumina-read2 <pe2> --long-read-file <lr> --output-directory <out_dir>

Required arguments:
    
      --short-read1          STR   fasta file of illumina read1 <pe1>
      --short-read2          STR   fasta file of illumina read2 <pe2>
      --long-read            STR   fasta file of long reads <lr>
      --out-dir              STR   output directory for scaffolding results <out_dir>

Optional arguments:
   
    Algorithm options:
      --no-ref-clustering          disable reference level clustering
      --no-strain-clustering       disable strain level clustering
      --no-gap-filling             disable gap-filling stage
      --polishing                  enable assembly polishing (currently using Pilon)
      --long-read-mapper     STR   software used for long-read mapping i.e. blasr or minimap2 [blasr]
      --short-read-assembler STR   software used for short-read assembly i.e. megahit or spades [megahit]
      --kmer-size            INT   kmer value used to assemble contigs [60]
      --contig-len-thr       INT   contig length threshold for clustering; contigs smaller than the threshold will be filtered out [500]
      --contig-edge-len      INT   during contig coverage calculation, number of bases filtered out from each contig end, to avoid biases due to lower mapping efficiency [80]
      --contig-window-len    INT   window length in which the coverage estimation is performed. We recommend using contig-len-thr - 2 * contig-edge-len as the value [340]
		
   Other arguments:
      --contig-file          STR   path to the contig file, if the short-reads have been assembled previously [default assembly using MEGAHIT]
      --num-processors       INT   number of processors to use (note that 2 is the minimum) [2]

";
}



