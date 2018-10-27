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

my $bundler_path = $bin_dir . "bundler"; 
my $opera_path = $bin_dir . "opera";
my $sigma_path = $bin_dir . "sigma";
my $compute_n50_path = $bin_dir . "scaffold_stats_opt.pl";



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
			die "Library file for computing coverage:"
			    .$mapping_file." not found";  
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
		if (! -e $LIB{$int_lib_name}) {
		    die "Library (BAM) file: ".$LIB{$int_lib_name}." not found\n";  
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
	    case "LONG_READ_MAPPER"{
	    }
	    
	    case "BLASR_DIR"{
		$blasr_dir = $split_line[1];
	    }

	    case "SHORT_READ_TOOL"{
	    }
	    case "STRAIN_CLUSTERING"{
	    }
	    case "BWA_DIR"{
		$short_read_tool_dir = $split_line[1];
	    }

	    
	    case "GRAPHMAP_DIR"{
	    }
	    case "RACON_DIR"{
	    }
	    
	    case "MUMMER_DIR"{
	    }
	    
	    case "MINIMAP2_DIR" {
	    }
	    
	    case "MEGAHIT_DIR" {
	    }
	    
	    case "MASH_DIR"{
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
	    
	    case "LONG_READ"{
	    }
	    
	    case "ILLUMINA_READ_1"{
	    }
	    
	    case "ILLUMINA_READ_2"{
	    }
	    
	    case "NUM_PROCESSOR"{
	    }
	    
	    else {
		die "Config option: ".$config_option." unknown";
	    }	
	}
    }
}
close($opera_ms_cf);


# If coverage has already been computed, it is stored in a sigma_contigs_files, 
# else, mapping_files needs to be provided
if(!(defined $mapping_files || defined $sigma_contigs_file)) {
	die "Calculating coverage information needs a sigma contigs file or
			mapping files.\n"; 
}


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

#foreach $k (sort keys %LIB) {
#	my $lib_size_starting_index;
#	my $lib_size_ending_index;
#	my $lib_size;
#	$lib_size_starting_index = rindex($LIB{$k}, "/") + 1;
#	$lib_size_ending_index = rindex($LIB{$k}, "_");
#	$lib_size = substr($LIB{$k}, $lib_size_starting_index, $lib_size_ending_index - $lib_size_starting_index);

#	$results_folder .= $lib_size."_";
#}
#chop($results_folder);

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
	
	$edges_file = $output_dir."/".$k."_bundles/clusters_contigs";
	#NEED TO CHANGE THAT STUFF !!!!	#$edges_file = $output_dir."/".$k."_bundles/clusters_".$lib_name;
	$edges_files .= $edges_file.",";

	# check which libraries have already been computed
	if(! -e $edges_file) {
	    print("***".$edges_file." needs to be computed\n");
	    push(@edges_files_not_computed, $k);
	}
}
chop($edges_files); # remove last comma


print STDERR "\nConstructing edge bundles ... \n"; 
# Check which files are not yet computed and compute only those
if (scalar(@edges_files_not_computed) == 0) {
	print STDERR "All edges files previously computed\n"; 
}
else {
	my $key_lib;
	my %LIB_NOT_COMPUTED;
	print STDERR "\nSome edges files need to be computed:\n"; 
	for $key_lib (@edges_files_not_computed) {
		print STDERR "\nComputing ".$LIB{$key_lib}."\n"; 
		# call bundling with only those libs that have not yet been computed
      	$LIB_NOT_COMPUTED{$key_lib} = $LIB{$key_lib}
  	}
  	&bundleBAMs($contigs_file, \%LIB_NOT_COMPUTED, $bundle_size_thresh); 
}
print STDERR "Done!\n"; 




############################# SIGMA ##############################

print STDERR "\nRunning Sigma ... \n"; 
my $sigma_dir = "$results_folder/coverage_estimation"; 
unless (-d $sigma_dir) {
    mkdir $sigma_dir or die $!; 
}

my $sigma_output = "NodePartitions.sigma"; 
construct_sigma_config_file($sigma_dir, $contigs_file_type, $contigs_file, 
	$edges_files, $sigma_contigs_file, $contig_len_thr, 
	$contig_edge_len, $contig_window_len, $pdist_type, $bundle_size_thresh);

run_exe("$sigma_path $sigma_dir/sigma.config");

# CHECK IF COVERAGE NEEDS TO BE COMPUTED
#if (0 && !defined $sigma_contigs_file) {			#coverage not computed
#	run_exe("$samtools_path view $mapping_files | $sigma_path $sigma_dir/sigma.config");
#
#} else { 									#coverage computed
#	run_exe("$sigma_path $sigma_dir/sigma.config");
#}


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


print STDERR "\nRunning Opera on partitions ... \n"; 
my $opera_output_folder = "$results_folder/scaffolds";
unless (-d $opera_output_folder){
    mkdir $opera_output_folder  or die $!;  
}
my $opera_config_file = "opera.config"; 
&construct_opera_config_file($contigs_file, $opera_output_folder, \%opera_filtered_files, $opera_output_folder."/".$opera_config_file);

#Uncomment if you want to run OPERA.
#run_exe("$opera_path $opera_output_folder/opera.config > $results_folder/log.txt");
#
##$res_file = "$results_folder/scaffoldSeq.fasta";
#$res_file = "$output_dir/scaffoldSeq.fasta";
#run_exe("rm -f $res_file"); #if(-e $res_file);
#run_exe("ln -s $opera_output_folder/scaffoldSeq.fasta $res_file");
#
#print STDERR "Done!\n"; 
#
##chdir $curDir or die $!; 
#print STDERR "\nALL done.\n"; 
#print STDERR "Result file: $res_file\n\n"; 
#
## COMPUTE N50
#print STDERR "\nComputing N50...\n";
#run_exe("perl $compute_n50_path $res_file > ".$res_file.".stats");
#print STDERR "\nDone!\n";

############################## HELPERS ##############################

#Constructs config files and runs bundling for all libraries in lib_ref
sub bundleBAMs {
    my ($contigs_file, $lib_ref, $bundle_size_thresh) = @_;
    my %libH = %{$lib_ref};
    foreach $k (keys %libH) {
    	my $lib_bundles_dir = $output_dir."/".$k."_bundles";
		my $lib_bundles_config = $k.".ebconfig"; 
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
        run_exe("$bundler_path $lib_bundles_dir/$lib_bundles_config");
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
sub construct_sigma_config_file {
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


