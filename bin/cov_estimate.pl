#!/usr/bin/perl 

#@authors: M. Senthil Kumar and Denis Bertrand's, Genome Institute of Singapore, 2013., version OperaMS1.0,
#          Janja Paliska, 						  Genome Institute of Singapore, 2014., version OperaMS2.0
#DEPENDENCY: 
#			<GENERIC>: samtools
#			<OPERA>   : Opera's edge bundler, and Opera
#
#
# NOTE: Pass a flag -r 0 before config file for only printing
#       the commands from script without executing them.
#		If r is not defined or is defined as -r 1, commands
#		will be both printed and executed.



use warnings; 
use Cwd;
use File::Spec; 
use Switch; 
use Getopt::Std;
use File::Spec::Functions qw(rel2abs);
use File::Basename;

getopts('r:'); # r=1 -> run cmds in run_exe, else only print


my $main_direct = getcwd;
$main_direct .= "\/";

my $bin_dir = dirname(rel2abs($0)) . "/";

############################## PARSING CONFIG FILE ##############################

my($opera_ms_config_file) = @ARGV;
my($contigs_file, %LIB, $samtools_path, $mapping_files, $sigma_contigs_file,
	$output_dir, $bundle_size_thresh, $kmer_size, $contigs_file_type, 
	$contig_len_thr, $contig_edge_len, 
   $contig_window_len, $pdist_type,
   $short_read_tool_dir, $blasr_dir
    );
my $num_LIB = 0;

$contigs_file_type = "SOAP";
print STDERR "\nReading config file: ".$opera_ms_config_file."\n"; 
open($opera_ms_cf, "<", $opera_ms_config_file); 
while(<$opera_ms_cf>) {
    next if (/^#/);  #skip comments
    chomp($_); 		 
    my @split_line = split('\s+', $_);
    if (@split_line != 0) {
	$config_option = $split_line[0];

	switch ($config_option) {
            
	    case "CONTIGS_FILE" {
		$contigs_file = File::Spec->rel2abs($split_line[1]);
		if (! -e $contigs_file) {
		    die "Contigs file: ".$contigs_file." not found";  
		} 
	    }

	    case "MAPPING_FILES" {
		$mapping_files = $split_line[1]; #comma separated files
		#check if all of them exist
		my @split_mapping_files = split(',', $mapping_files);
		my $mapping_file;
		for $mapping_file (@split_mapping_files) {
		    if (! -e $mapping_file) {
			#die "Library file for computing coverage:"
			#    .$mapping_file." not found";  
		    } 					
		}
	    }

	    
	    case "LIB" {
		my $lib_file_name = $split_line[1]; #filename
		$num_LIB++; 
		my $int_lib_name = "lib_$num_LIB";
		if (@split_line == 3) {
		    $int_lib_name = $split_line[2]; #internal names of library for producing bundles
		}

		if (exists $LIB{$int_lib_name}) {
		    die "Same library names given for different libraries.\n"; 
		}

		$LIB{$int_lib_name} = $lib_file_name;
		print STDERR " *** Add lib $lib_file_name\n";
		
		if (! -e $LIB{$int_lib_name}) {
		    #die "Library (BAM) file: ".$LIB{$int_lib_name}." not found\n";  
		}
	    }

	    case "EDGE_BUNDLESIZE_THRESHOLD" {
		$bundle_size_thresh = $split_line[1]; 
	    }

	    case "OUTPUT_DIR" {
		$output_dir = $split_line[1];
	    }

	    case "SAMTOOLS_DIR" {
		$samtools_path = $split_line[1];
		if (! -e $samtools_path) {
		    die "Samtools not found at: ".$samtools_path."\n";
		}
	    }
	    
	    case "BWA_DIR"{
		$short_read_tool_dir = $split_line[1];
	    }
	    
	    case "SIGMA_CONTIGS_FILE" {
		$sigma_contigs_file = File::Spec->rel2abs($split_line[1]); 
		if(! -e $sigma_contigs_file) {
		    die "Sigma contigs file:"
			.$sigma_contigs_file. " not found.\n"; 
		}	
	    }

	    case "KMER_SIZE" {
		$kmer_size = $split_line[1];
		if (!defined $kmer_size) {
		    die "KMER_SIZE not provided.\n";
		}
	    }
	    
	    case "CONTIGS_FILE_TYPE" {
		$contigs_file_type = $split_line[1];
	    }
	    
	    case "CONTIG_LEN_THR" {
		$contig_len_thr = $split_line[1];
	    }
	    
	    case "CONTIG_EDGE_LEN" {
		$contig_edge_len = $split_line[1];
	    }
	    
	    case "CONTIG_WINDOW_LEN" {
		$contig_window_len = $split_line[1]
	    }

	    case "PDIST_TYPE" {
		$pdist_type = $split_line[1];
	    }
	    
	    case "ILLUMINA_READ_1"{
		$illum_read1 = $split_line[1];
	    }
	    
	    case "ILLUMINA_READ_2"{
		$illum_read2 = $split_line[1];
	    }
	    
	    case "NUM_PROCESSOR"{
		$nproc = $split_line[1];
	    }
	    
	    else {
		#die "Config option: ".$config_option." unknown";
	    }	
	}
    }
}
close($opera_ms_cf);


if(!(defined contigs_file && defined contigs_file_type)) {
	die "Contigs file and contigs_file_type need to be provided. \n";
}
 

if(!defined $bundle_size_thresh) {
	$bundle_size_thresh = 2; 
}

if(!defined $output_dir) {
	$output_dir = "OperaMS"; 
}
unless (-d $output_dir) {
       mkdir $output_dir or die $!;        
}


# Create additional sigma and opera output folders, used mainly 
# for running multiple library combinations simultaneously.
# Output folder name is: $output_dir/ + all libraries names separated by "_"
# NEED TO CHANGE THAT STUFF !!!!
my $results_folder = $output_dir."/"."intermediate_files";
$output_dir = $results_folder;

unless (-d $results_folder) {
       mkdir $results_folder or die $!;        
}


############################## BUNDLING ##############################

# Create comma separated paths to edges files for sigma and
# check which files need to be constructed 
my $edges_files = "";
my @edges_files_not_computed;
foreach $k (keys %LIB) {
	my $lib_name_starting_index;
	my $lib_name_ending_index;
	my $lib_name;
	$lib_name_starting_index = rindex($LIB{$k}, "/") + 1;
	#$lib_name_ending_index = rindex($LIB{$k}, ".");
	$lib_name_ending_index = index($LIB{$k}, ".");#NEED TO UPDATE THAT TO AVOID PROBLEMS WITH DIRECTORY NAMES WITH .
	$lib_name = substr($LIB{$k}, $lib_name_starting_index, $lib_name_ending_index - $lib_name_starting_index);
	
	$edges_file = $output_dir."/".$k."_bundles/clusters_opera";#Use the mapping obtained at previous stage
	#NEED TO CHANGE THAT STUFF !!!!	#$edges_file = $output_dir."/".$k."_bundles/clusters_".$lib_name;
	$edges_files .= $edges_file.",";

	# check which libraries have already been computed
	if(! -e $edges_file) {
	    print("***".$edges_file." needs to be computed\n");
	    push(@edges_files_not_computed, $k);
	}
}
chop($edges_files); # remove last comma


#print STDERR "\nConstructing edge bundles ... \n"; 
# Check which files are not yet computed and compute only those
my $key_lib = "lib_1";
my %LIB_NOT_COMPUTED;
#print STDERR "\nSome edges files need to be computed:\n"; 
#print STDERR "\nComputing ".$LIB{$key_lib}."\n"; 
# call bundling with only those libs that have not yet been computed
$LIB_NOT_COMPUTED{$key_lib} = $LIB{$key_lib};
my $lib_bundles_dir = $output_dir."/".$key_lib."_bundles";
my $lib_bundles_config = $key_lib.".ebconfig"; 
generate__bundler_config_file($contigs_file, \%LIB_NOT_COMPUTED, $bundle_size_thresh, $lib_bundles_dir, $lib_bundles_config); 


############################# SIGMA ##############################

#print STDERR "\nRunning Sigma ... \n"; 
my $sigma_dir = "$results_folder/coverage_estimation"; 
unless (-d $sigma_dir) {
    mkdir $sigma_dir or die $!; 
}

my $sigma_output = "NodePartitions.sigma"; 
generate_sigma_config_file($sigma_dir, $contigs_file_type, $contigs_file, 
	$edges_files, $sigma_contigs_file, $contig_len_thr, 
	$contig_edge_len, $contig_window_len, $pdist_type, $bundle_size_thresh);

get_read_size($illum_read1, "$sigma_dir/read_size.dat");

my $preprocess_read_path = $bin_dir . "../OPERA-LG/bin/preprocess_reads.pl";
my $short_analysis_path = $bin_dir . "short-read-analysis";
my $utils_dir  = $bin_dir . "../utils";
run_exe("$utils_dir/perl $preprocess_read_path --out $sigma_dir/short_read_analysis --tool-dir $utils_dir --nproc $nproc --contig $contigs_file --illumina-read1 $illum_read1 --illumina-read2 $illum_read2 --sigma-conf $sigma_dir/sigma.config --bundler-conf $lib_bundles_dir/$lib_bundles_config 2> $sigma_dir/preprocess_reads.err");
if($? || ! -e "$sigma_dir/assembly_size.dat" || ! -e "$lib_bundles_dir/lib.txt"){
    die "Error during short read pre-preocessing. Please see $sigma_dir/preprocess_reads.err $sigma_dir/short_read_analysis.out $sigma_dir/short_read_analysis.err for details.\n";
}	


############################## OPERA ##############################

# Extract filtered files, mean, and standard deviation for constructing
# opera config file
my %opera_filtered_files = ();
foreach $k (keys %LIB) {
	my $lib_name_starting_index;
	my $lib_name_ending_index;
	my $lib_name;
	
	$lib_name_starting_index = rindex($LIB{$k}, "/") + 1;
	$lib_name_ending_index = index($LIB{$k}, ".");#NEED TO UPDATE THAT TO AVOID PROBLEMS WITH DIRECTORY NAMES WITH .
	$lib_name = substr($LIB{$k}, $lib_name_starting_index, $lib_name_ending_index - $lib_name_starting_index);
	my $filtered_file = $sigma_dir."/filtered_clusters_contigs";
	#NEED TO CHANGE THAT STUFF !!!! #my $filtered_file = $sigma_dir."/filtered_clusters_".$lib_name;
	my $lib_info_file = $output_dir."/".$k."_bundles/lib.txt";
	my $sd;
	my $mean;
	open($lib_info, "<", $lib_info_file);
	while (<$lib_info>) {
		chomp($_);
		if (/^Mean/) {
			my @line = split(":", $_);
			$mean = $line[1];	
			$mean =~ s/^\s+|\s+$//g; # remove leading/trailing white space
		}
		elsif (/^Standard/) {
			my @line = split(":", $_);
			$sd = $line[1];
			$sd =~ s/^\s+|\s+$//g; # remove leading/trailing white space 
		}
 
	}
  
	my @lib_stats = ($mean, $sd, $filtered_file);
	$opera_filtered_files{$k} = \@lib_stats;
}

############################## HELPERS ##############################


sub get_read_size{
    my ($illum_read1, $out_file) = @_;

    my $read_size = 0;
    my $nb_read_read = 1000;
    my $nb_read = 0;
    $cmd_open = $illum_read1;
    $cmd_open = "zcat $illum_read1 |" if(index($illum_read1, ".gz") != -1);
    $read_size = 0;
    open(FILE, $cmd_open);
    while(<FILE>){
	last if($nb_read == $nb_read_read);
	$seq = <FILE>;
	$size = length($seq);
	$read_size = $size if($read_size < $size);
	<FILE>;<FILE>;
	$nb_read++;
    }
    close(FILE);

    open(OUT, ">$out_file");
    print OUT $read_size . "\n";;
    close(OUT);
}

#Constructs config files and runs bundling for all libraries in lib_ref
sub generate__bundler_config_file {
    my ($contigs_file, $lib_ref, $bundle_size_thresh, $lib_bundles_dir, $lib_bundles_config) = @_;
    my %libH = %{$lib_ref};
    foreach $k (keys %libH) {
	run_exe("mkdir $lib_bundles_dir");
        open(my $config, ">", "$lib_bundles_dir/$lib_bundles_config"); #edge bundler config
	print $config "\nsamtools_dir=$samtools_path\n";
        print $config "\noutput_folder=$lib_bundles_dir\n";
        print $config "\ncontig_file=$contigs_file\n";
        print $config "\n[LIB]"; 
        print $config "\nmap_file=$libH{$k}"; #this will cause the bundler to read the read mapping from STDIN
        print $config "\nfilter_repeat=no\n"; 
        print $config "\ncluster_threshold=$bundle_size_thresh\n";
        print $config "\ncluster_increased_step=5\n"; 
        print $config "\nkmer=$kmer_size\n";
	close($config); 
    }
}

# Constructs opera config file 
sub construct_opera_config_file {
    my ($contigs_file, $outputFolder, $LIB_ref, $opera_config_fname) = @_; 
    my %libH = %{$LIB_ref};
    open(my $config, ">", $opera_config_fname); #edge bundler config
    print $config "\noutput_folder=$outputFolder\n";
    print $config "\ncontig_file=$contigs_file\n";
    foreach $k (keys %libH) {
    	print $config "\n[LIB]"; 
    	print $config "\nmap_type=opera"; 
    	print $config "\nmap_file=$libH{$k}->[2]"; 
    	print $config "\nlib_mean=$libH{$k}->[0]";
 		print $config "\nlib_std=$libH{$k}->[1]\n";	
    }
    print $config "\nkeep_repeat=yes\n";
    print $config "\nfilter_repeat=no\n"; 
    print $config "\ncluster_threshold=$bundle_size_thresh\n";
    print $config "\ncluster_increased_step=5\n"; 
    print $config "\nkmer=$kmer_size\n";
    close($config);
}


# Prints out command given as a input parameter AND
# runs it if -r option in command line is not defined
# or is defined as -r 1
sub run_exe {
    my ($exe) = @_;
    
    if (!defined $opt_r || $opt_r==1) {
		print STDERR $exe."\n";
		print STDERR `$exe`;
	} 
    else {
		print STDERR $exe."\n";
	}
}


# Constructs sigma config file
sub generate_sigma_config_file {
    my ($output_dir, $contigs_file_type, $contigs_file, $edges_files, $sigma_contigs_file,
	$contig_len_thr, $contig_edge_len, $contig_window_len, $pdist_type, 
	$bundle_size_thresh) = @_; 
	
	$config_file = $output_dir."/sigma.config";
	open(CONF, ">$config_file");


	####### NOTE ##############################################
	# If coverage is already computed, then it is contained in 
	# sigma_contigs_file so do not specify that file
	#
	# Else:
    # specify sigma_contigs_file as contigs_ + contig_window_len + _ + contig_edge_len,
	# put mapping_files=- (it will be read from stdin)
	# ->coverage will be computed and written in specified sigma_contigs_file
	#############################################################

	if (!defined $sigma_contigs_file) {
	    #print CONF "\nmapping_files=-\n";
	    print CONF "\nmapping_files=$mapping_files\n";
	    print CONF "\nsigma_contigs_file=$output_dir/contigs_".$contig_window_len."_".$contig_edge_len."\n";
	} 
	else {
		print CONF "\nsigma_contigs_file=$sigma_contigs_file\n";
	}

	print CONF "\ncontigs_file_type=$contigs_file_type\n";
	print CONF "\ncontigs_file=$contigs_file\n";
	print CONF "\nedges_files=$edges_files\n";
	print CONF "\noutput_dir=$output_dir\n";
	print CONF "\nedge_bundle_size_thr=$bundle_size_thresh\n";
    print CONF "\nkmer_size=$kmer_size\n";

	# default contig length threshold = 500
	if (defined $contig_len_thr) {
		print CONF "\ncontig_len_thr=$contig_len_thr\n";
	} 
	else {
		print CONF "\ncontig_len_thr=500\n";
	}
	# default contig edge length = 0
	if (defined $contig_edge_len) {
		print CONF "\ncontig_edge_len=$contig_edge_len\n";
	} 
	else {
		print CONF "\ncontig_edge_len=0\n";
	}
	# default contig window length = 0
	if (defined $contig_window_len) {
		print CONF "\ncontig_window_len=$contig_window_len\n";
	} 
	else {
		print CONF "\ncontig_window_len=0\n";
	}
	#default distribution type = Negative binomial
	if (defined $pdist_type) {
		print CONF "\npdist_type=$pdist_type\n";
	} 
	else {
		print CONF "\npdist_type=NegativeBinomial\n";
	}

    print CONF "\nsamtools_dir=$samtools_path\n";
    
	close(CONF);
 }


