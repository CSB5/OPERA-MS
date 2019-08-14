#!/usr/bin/perl
#use strict;
use warnings;

my ($working_dir, $contig_file, $read_file, $read_mapping_dir, $nb_process, $scaffold_length_threshold, $opera_bin_dir, $racon_dir,  $minimap2_dir, $mummer_dir) = @ARGV;

require "$opera_bin_dir/../../bin/test_time.pl";

run_exe("mkdir $working_dir") if(! -d $working_dir);

my $perl_path = "$opera_bin_dir/../../utils/perl";
my $start_time = time;
my $start_time_begin = time;
my $end_time;
#split the scaffold
my $opera_lr_dir = "$working_dir/../";
my $assembly_scaffold_file = "$opera_lr_dir/scaffolds.scaf";
my $assembly_scaffold_seq_file = "$opera_lr_dir/scaffoldSeq.fasta";

#print STDERR " *** $assembly_scaffold_file\n";
my %scaff_info = ();
my %singleton_contig = ();
my @scaffold_order = ();
split_scaffold_file($working_dir, $assembly_scaffold_file);
$end_time = time;
#print STDERR "***  Construt scaffold dir Elapsed time: " . ($end_time - $start_time) . "s\n";
write_time("$working_dir/../../", "g_dir_constriction", ($end_time - $start_time));
$start_time = time;

################oNLY FOR TESTIN PURPOSE
#merge_contig_from_scaf($working_dir, $contig_file, "$working_dir/../scaffoldSeq.fasta.filled");
#exit(0);
######################################3

#Get the reads that map to a single contig
my $edge_read_info_file = "$read_mapping_dir/edge_read_info.dat";
my $long_read_mapping = "$read_mapping_dir/opera.map.sort.status";

#Actulally never used
#get_read_on_single_contig($working_dir, $long_read_mapping);
#

run_exe("$perl_path $opera_bin_dir/extract_read_sequence.pl --edge-file $edge_read_info_file --contig-file $contig_file --opera-lr-dir $opera_lr_dir --read-file $read_file --output-directory $working_dir 2> $working_dir/extract_read.err");
if($?){
    die "Error in during read sequence extraction. Please see $working_dir/extract_read.err for details.\n";
}

$end_time = time;
write_time("$working_dir/../../", "g_extract_read", ($end_time - $start_time));
#print STDERR "***  Identify read in gap and construt pre_consensus assembly Elapsed time: " . ($end_time - $start_time) . "s\n";
$start_time = time;



#exit(0);

##########################################Get the consensus sequence using minimap2 + racon
#
opendir(DIR, $working_dir);
@all_dir = readdir(DIR);
close(DIR);
my $nb_tiling_run = 0;
my $consensus_cmd = "$working_dir/consensus_cmd.sh";
open(OUT, ">$consensus_cmd");
foreach $scaff_dir (@all_dir){
    #next if($scaff_dir ne "opera_scaffold_1");
    next if($scaff_dir eq ".." || ! -e "$working_dir/$scaff_dir/scaffolds.scaf");
    
    print OUT "$perl_path $opera_bin_dir/run_scaffold_racon.pl $working_dir/$scaff_dir/ $racon_dir $minimap2_dir 2> $working_dir/$scaff_dir/consensus.err\n";
    
    #last if($nb_tiling_run == 10);
    #$nb_tiling_run++;
    #last;
}
close(OUT);
$consensus_log = "$consensus_cmd.log";
run_exe("cat $consensus_cmd | xargs -L 1 -P $nb_process -I COMMAND sh -c \"COMMAND\" 2> $consensus_log");
if($?){
    die "Error during consensus sequence generation. Please see $consensus_log for details.\n";
}
#Combine all the consensus into a single file
my $tilling_dir = "$working_dir/TILLING";
run_exe("mkdir $tilling_dir");
my $consensus_file = "$tilling_dir/consensus.fa";
#run_exe("cat $working_dir/*opera_scaffold_*/racon.fa | sed 's/Consensus_//' > $consensus_file");
run_exe("rm $consensus_file; for f in $working_dir/*opera_scaffold_*/racon.fa; do sed 's/Consensus_//' \$f >> $consensus_file; done");
#
#
#
$end_time = time;
write_time("$working_dir/../../", "g_racon", ($end_time - $start_time));
#print STDERR "***  Consensus sequence creation Elapsed time: " . ($end_time - $start_time) . "s\n";
$start_time = time;

#Perform the scaffold tilling to add back the high bp quality illumina contigs

#########################################First tilling
#
my $first_filling =  "$tilling_dir/first_tilling";
if(1 || ! -e $first_filling){
    $start_time = time;
    run_exe("$perl_path $opera_bin_dir/run_mummer_large_ref.pl $consensus_file $tilling_dir/REF $contig_file $tilling_dir/QUERY $working_dir $first_filling.coords $nb_process $mummer_dir > $tilling_dir/tilling_1.out 2> $tilling_dir/tilling_1.err");
    if($?){
	die "Error during tilling generation. Please see $tilling_dir/tilling_1.out and $tilling_dir/tilling_1.err for details.\n";
    }
    #$end_time = time;
    #print STDERR "***  Mapping (1) completed Elapsed time: " . ($end_time - $start_time) . "s\n";
}

#$start_time = time;
get_paf_file($assembly_scaffold_file, "$first_filling.coords", "$first_filling.paf", "ONLY_SCAFF_CONTIG");
#$end_time = time;
#print STDERR "***  Get paf file Elapsed time: " . ($end_time - $start_time) . "s\n";

get_contig_in_gap("$first_filling.paf", $contig_file, "$first_filling\_contig.fa");

run_exe("$opera_bin_dir/Remap.py $first_filling.paf $consensus_file $first_filling\_contig.fa > $tilling_dir/tilling_1_remap.out");
if($?){
    die "Error during contig remapping. Please see $tilling_dir/tilling_1_remap.out for details.\n";
}
$end_time = time;
write_time("$working_dir/../../", "g_tilling_1", ($end_time - $start_time));
#print STDERR "***  Get first tilling Elapsed time: " . ($end_time - $start_time) . "s\n";

#[minor change with previous version]
#Potential rescue of pre-consensus files if the tilling failed => no contigs mapped to the scaffold may only happen to very short scaffolds
#Rescue pre-consensus contigs if there is no valid tilling available for a given scaffold. 
open(FILE, "$tilling_dir/consensus_scaff_to_rescue.dat");
while(<FILE>){
    chop $_;
    $scaff = $_;
    run_exe("cat $working_dir/$scaff/pre_consensus.fa >> $tilling_dir/consensus_remapped.fasta");
}
close(FILE);

#########################################Second tilling
#
my $second_tilling =  "$tilling_dir/second_tilling";
if(1 || ! -e $second_tilling){
    $start_time = time;
    run_exe("$perl_path $opera_bin_dir/run_mummer_large_ref.pl $tilling_dir/consensus_remapped.fasta $tilling_dir/REF $contig_file $tilling_dir/QUERY $tilling_dir/ $second_tilling.coords $nb_process $mummer_dir  > $tilling_dir/tilling_1.out 2> $tilling_dir/tilling_1.err");
    if($?){
	die "Error during tilling generation. Please see $tilling_dir/tilling_2.out and $tilling_dir/tilling_2.err for details.\n";
    }
}

get_paf_file($assembly_scaffold_file, "$second_tilling.coords", "$second_tilling.paf",  "ALL");

get_contig_in_gap("$second_tilling.paf", $contig_file, "$second_tilling\_contig.fa");
    
run_exe("$opera_bin_dir/Remap.py $second_tilling.paf $tilling_dir/consensus_remapped.fasta $second_tilling\_contig.fa >  $tilling_dir/tilling_2_remap.out");
if($?){
    die "Error during contig remapping. Please see $tilling_dir/tilling_2_remap.out for details.\n";
}
$end_time = time;
#print STDERR "***  Get second tilling Elapsed time: " . ($end_time - $start_time) . "s\n";
write_time("$working_dir/../../", "g_tilling_2", ($end_time - $start_time));

#Rescue pre-consensus contigs if there is no valid tilling available for a given scaffold.
#This can only happen if mummer does not remap the contigs mapped at the previous stage
open(FILE, "$tilling_dir/consensus_remapped_scaff_to_rescue.dat");
while(<FILE>){
    chop $_;
    $scaff = $_;
    run_exe("cat $working_dir/$scaff/pre_consensus.fa >> $tilling_dir/consensus_remapped_remapped.fasta");
}
close(FILE);

run_exe("rm -r $working_dir/*opera_scaffold_*");

##########################################merge the contigs into the gapfilled file
#
$start_time = time;
merge_contig_from_scaf("$tilling_dir/consensus_remapped_remapped.fasta", $assembly_scaffold_seq_file, "$working_dir/../scaffoldSeq.fasta.filled");
$end_time = time;
print STDERR "*** Contig merging Elapsed time: " . ($end_time - $start_time) . "s\n";

print STDERR "***  Get final assembly Elapsed time : " . ($end_time - $start_time_begin) . "s\n";
    
sub merge_contig_from_scaf{
    my ($filled_scaff_file, $scaffold_seq_file, $out_file) = @_;
    
    my %final_scaff_seq = ();
    open(FILE, $scaffold_seq_file);
    while(<FILE>){
	if($_ =~ /^>(.+)/){
	    $scaffold_full_name = $1;
	    @tmp = split(/\s+/, $scaffold_full_name);
	    $scaff_name = $tmp[0];
	    #print STDERR $contig_name . "\n";# <STDIN>;
	    if($scaff_info{$scaff_name} ne "MULTI"){#This scaffold was not filled
		$seq = <FILE>;
		$final_scaff_seq{$scaff_name} = $seq;
		#print STDERR " *** Print singleton contig " . $contig_name . "\n"; #<STDIN>;
	    }
	    else{
		<FILE>;#Scaffold was filled
	    }
	}
	
    }
    close(FILE);    
	
    open(FILE, $filled_scaff_file);
    while(<FILE>){
	if($_ =~ /^>(.+)/){
	    $scaff_name = $1;
	    $seq = <FILE>;
	    print STDERR $scaff_name . "\n";
	    $final_scaff_seq{$scaff_name} = $seq;
	}
    }
    close(FILE);
    
    open(OUT, ">$out_file");
    foreach $scaff_name (@scaffold_order){
	print STDERR $scaff_name . "\n";
	$seq = $final_scaff_seq{$scaff_name};
	$l = length($seq) -1;
	print OUT ">$scaff_name length: $l cov: 0.0\n";
	print OUT $seq;
    }
    close(OUT);
}

sub split_scaffold_file{
    my ($working_dir, $assembly_scaffold_file) = @_;
    open(FILE, "$assembly_scaffold_file");
    my $nb_contig = 0;
    my $scaff_name = "";
    my $scaff_length = 0;
    while(<FILE>){
	if($_ =~ /^>(.+)/){
	    if($nb_contig > 0){
		if($nb_contig > 1 && $scaff_length > $scaffold_length_threshold){
		    $scaff_dir = "$working_dir/$scaff_name";
		    `mkdir $scaff_dir` if(! -d $scaff_dir);
		    open(OUT, ">$scaff_dir/scaffolds.scaf");
		    print OUT $str_scaff;
		    close(OUT);
		    $scaff_info{$scaff_name} = "MULTI";
		}
		else{
		    $scaff_info{$scaff_name} = "NO_FILL";
		}
	    }
	    $str_scaff = $_;
	    @tmp = split(/\s+/, $1);
	    $scaff_name = $tmp[0];
	    $scaff_length = $tmp[2];
	    #print STDERR $scaff_name . "\t" .   $scaff_length . "\n";<STDIN>;
	    push(@scaffold_order, $scaff_name);
	    $nb_contig = 0;
	    #print STDERR " *** $scaff_name\n";#<STDIN>;
	}
	else{
	    $nb_contig++;
	    $str_scaff .= $_;
	}
    }
    close(FILE);

    ####################get the last scaffold !!!
    if($nb_contig > 1 && $scaff_length > $scaffold_length_threshold){
	$scaff_dir = "$working_dir/$scaff_name";
	run_exe("mkdir $scaff_dir") if(! -d $scaff_dir);
	open(OUT, ">$scaff_dir/scaffolds.scaf");
	print OUT $str_scaff;
	close(OUT);
	$scaff_info{$scaff_name} = "MULTI";
    }
    else{
	$scaff_info{$scaff_name} = "NO_FILL";
    }
}


#Extract read that maps to a single contig to avoid adding long-read errors on long contigs during racon
sub get_read_on_single_contig{
    my ($working_dir, $status_maping_file) = @_;
    #print STDERR "get_read_on_single_contig  $working_dir $status_maping_file\n";<STDIN>;
    #print STDERR $status_maping_file . "\n";<STDIN>;
    open(FILE, "grep -v small-alignment $status_maping_file | grep -v small-contig | grep -v partial-match |");

    my $previous_read_name = "";
    my $flag_single_contig = 1;
    my $str_mapping;
    my %contig_read = ();
    my $current_contig = "";
    while(<FILE>){
	chop $_;
	@line = split(/\s+/, $_);
	#$status = $line[12];
	$read_name = $line[0];
	#next if(index($status, "overlapped") != -1 && index($status, "better") == -1);
	
	if($read_name ne $previous_read_name){
	    if($previous_read_name ne ""){
		if($flag_single_contig && $str_mapping ne ""){
		    #print $str_mapping . "\n";
		    $contig_read{$current_contig} = $str_mapping . "\n";
		}
	    }
	    $flag_single_contig = 1;
	    $str_mapping = "";
	    $current_contig = $line[1];
	}
	
	if($flag_single_contig == 1 && $str_mapping eq ""){
	    $str_mapping = $_;
	}
	else{
	    $flag_single_contig = 0;
	}
	$previous_read_name = $read_name;
    }
    close(FILE);
    
    #print STDERR " *** $working_dir\n";<STDIN>; 
    opendir(DIR, $working_dir);
    @all_dir = readdir(DIR);
    close(DIR);
    my $assembly_dir;

    foreach $scaff_dir (@all_dir){
	#next if($scaff_dir ne "opera_scaffold_365");
	next if(! -e "$working_dir/$scaff_dir/scaffolds.scaf");
	$assembly_dir =  "$working_dir/$scaff_dir/";
	`mkdir $assembly_dir` if(! -d $assembly_dir);
	$single_read_map = "$assembly_dir/single_read.map";
	add_single_read("$working_dir/$scaff_dir/scaffolds.scaf", "$single_read_map");
    }
}

sub add_single_read{
    my ($scaff_file, $out_read) = @_;

    #print STDERR " *** single_read_map $scaff_file, $out_read\n";
    
    open(OUT, ">$out_read");
    open(FILE, $scaff_file);
    <FILE>;
    while(<FILE>){
	@line = split(/\t/, $_);
	$contig = $line[0];
	if(exists $contig_read{$contig}){
	    print OUT $contig_read{$contig};
	}
    }
    close(FILE);
    close(OUT);
}

sub get_contig_in_gap{

    my ($mapping_file, $contig_file, $out_contig_file) = @_;


    print STDERR " *** get_contig_and_corrected_read_in_gap $mapping_file  $contig_file  $out_contig_file\n";#<STDIN>;
    
    my $OVERHANG_SEQ = 200;
    my @cut_coordinates = ();

    #Get the list of contig in the mapping file
    open(FILE, "$mapping_file");
    my %contig_mapped = ();
    while(<FILE>){
	@line = split(/\t/, $_);
	$contig_name = $line[0];
	$contig_mapped{">$contig_name"} = 1;
	$nb_contig_to_print++;

    }

    #get the sequence of the contigs in fasta format
    open(IN, $contig_file);
    open(OUT, ">$out_contig_file"); 
    my $print_flag = 0;
    while(<IN>){
	if(index($_, ">") == -1){
	    if($print_flag){
		print OUT $_;
		$nb_contig_to_print--;
		last if($nb_contig_to_print == 0);
	    }
	}
	else{
	    #print STDERR "---> $nb_contig\n" if($nb_contig % 500000 == 0);$nb_contig++;
	    my @line = split(/\s+/, $_);
	    my $ID = substr ($line[0], 0, length($line[0]));
	    #print STDERR $ID."\n"; <STDIN>;
	    if(exists $contig_mapped{$ID}){
		$print_flag = 1;
		#chop $_;
		#print $_."\n";
		print OUT $ID . "\n";
	    }
	    else{
		#print STDERR "**************** IN\n";
		#if($print_flag == 1);
		$print_flag = 0;
	    }
	}
    }
    close(IN);
    close(OUT);
    
}


sub get_paf_file{
    my ($scaff_file, $coord_file, $paf_file, $flag_set_type) = @_;

    print STDERR " ***\n$scaff_file\n$coord_file\n\n";#<STDIN>;

    #This is the set of contigs that should belong to set of contigs used for the tilling
    my %contig_set_to_keep = ();
    #Keep all the contitgs that belong to the current scaffolds
    #No identity threshold is applied to those contigs
    open(FILE, "$scaff_file");
    my $scaff_name = "";
    while(<FILE>){
	chop $_;
	@line = split(/\s+/, $_);
	$name = $line[0];
	if($name =~ m />(.+)/){
	    $scaff_name = $1;
	    $contig_set_to_keep{$scaff_name} = {};
	}
	else{
	    $contig = $name;
	    #print STDERR "ADD" . "\t" . $scaff_name . "\t" . $contig . "\n";
	    $contig_set_to_keep{$scaff_name}->{$contig} = 1; 
	}
    }
    close(FILE);

    #Add the contigs for the same species
    if($flag_set_type eq "ALL"){
	#get_contig_in_species($scaffold_dir, \%contig_set_to_keep);
    }
    
    open(OUT, ">$paf_file");
    open(FILE, "sort -k10,10 -k1,1n $coord_file | grep -v '\\[' | ");### NEED TO FIX
    my $gap_length = 0;
    #my $IDENTITY_THRESHOLD = 0.97;
    #
    my $IDENTITY_THRESHOLD = 0.95;#Low identy threshold to allows mapping of repeat contigs on gaps covered by low number of reads thereshold can be scaffold specific and dependent of the coverage => high threshold due to improve sequence quality
    #
    my $COVERAGE_THRESHOLD = 0.90;
    my $previous_contig_name = "";
    my $prev_scaffold_end = -1;
    my ($contig_name, $scaffold_name, $scaffold_length, $mapping_length, $contig_end, $contig_start, $scaffold_start, $scaffold_end, $tiling_length);
    $scaffold_name = "";
    open(OUT_G, ">$working_dir/S1_gap_size.dat");
    <FILE>;<FILE>;<FILE>;#<FILE>;#skip the header lines
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	$contig_name = $line[10];
	#print STDERR " XXXX $scaffold_name\n";#<STDIN>;
	
	if($contig_name ne $previous_contig_name){
	    if($previous_contig_name ne ""){
		
		$mapping_length = ($contig_end - $contig_start);
		$frac_mapping = ($mapping_length / $contig_len);
		#scaffold length
		#print STDERR " *** MAP " . $scaffold_name . "\t" . $contig_name . "\t" . $previous_contig_name . " " . $mapping_length . " " . $contig_len . " " . $contig_end . " " .$contig_start . "\t" . $frac_mapping . "\t" . $identity . " |" . (exists $contig_set_to_keep{$scaffold_name}) . "|" . (exists $contig_set_to_keep{$scaffold_name}->{$previous_contig_name}) . "|\n";
		
		if(($identity >= $IDENTITY_THRESHOLD || (exists $contig_set_to_keep{$scaffold_name}->{$previous_contig_name})) && #Identity threshold => always try to tile contig that belong to the $contig_set_to_keep
		   ($frac_mapping >= $COVERAGE_THRESHOLD || ($flag_set_type eq  "ONLY_SCAFF_CONTIG" && (exists $contig_set_to_keep{$scaffold_name}->{$previous_contig_name}) && $frac_mapping > 0.55))#More aggressive extention in case of contig that are know to belong to the scaffold, even partial alignements are allowed and will be extended
		    ){
		    $tiled_contig{$previous_contig_name} = 1;
		    
		    print OUT
			$previous_contig_name . "\t" . $contig_len . "\t" . $contig_start . "\t" . $contig_end . "\t" . $contig_orientation . "\t" . 
			$scaffold_name . "\t" . $scaffold_length . "\t" .
			$scaffold_start . "\t" . $scaffold_end . "\t" .
			int(($scaffold_end - $scaffold_start) * $identity) . "\t" . ($scaffold_end - $scaffold_start) . "\t" .
			"60" . "\t" . 
			"tp:A:P" . "\t". "cm:i:175" . "\t" . "s1:i:915"	. "\t" . "s2:i:0" . 
			"\n";

		    #Sum contig length
		    $sum_contig_length += $contig_len;
		    
		    #Tilling length
		    $tiling_length += $mapping_length;

		    if($prev_scaffold_end > 0 && $prev_scaffold_end < $scaffold_start){
			$gap_length += ($scaffold_start - $prev_scaffold_end);
			print  OUT_G ($scaffold_start - $prev_scaffold_end) . "\t" . $previous_contig_name . "\t" .$contig_name . "\n";
		    }
		}
		#
		$prev_scaffold_end = $scaffold_end;
	    }
	    
	    $previous_contig_name = $contig_name;
	    $scaffold_name = $line[9];
	    $scaffold_length = $line[7];
	    #
	    $identity = $line[6]/100;
	    $contig_len = $line[8];
	    $contig_start = $line[2];
	    $contig_end = $line[3];
	    #$identity = $line[6]/100 * ($contig_end -$contig_start);
	    #Get the contig orientation and reverse coordinate in case of -
	    $contig_orientation = "+";
	    if($contig_end < $contig_start){
		$contig_orientation = "-";
		$tmp =  $contig_end; $contig_end = $contig_start;$contig_start = $tmp;
	    }
	    $contig_start--;
	    #
	    $scaffold_start = $line[0] - 1;
	    $scaffold_end = $line[1];
	}
	else{
	    #The alignement can be extended
	    #NEED TO BE CAREFUL IN CASE OF TANDEM REPEAT !!!!

	    #Check if this a correct extention
	    if($contig_orientation eq "+"){
		if($contig_end < $line[3]){#this is an extention
		    $contig_end = $line[3];#the end is extended with next end
		    $scaffold_end = $line[1];
		}
	    }
	    else{
		if($contig_start > $line[3]-1){#this is an extention
		    $contig_start = $line[3] -1;#the end is extended with next start
		    $scaffold_end = $line[1];
		}
	    }
	    #$identity = $line[6]/100 * ($contig_end -$contig_start);
	    
	}
	#$previous_contig_name = $contig_name;
    }
    close(FILE);

    #For the last contig
    $tiled_contig{$previous_contig_name} = 1;
    print OUT
	$previous_contig_name . "\t" . $contig_len . "\t" . $contig_start . "\t" . $contig_end . "\t" . $contig_orientation . "\t" .
	$scaffold_name . "\t" . $scaffold_length . "\t" .
	$scaffold_start . "\t" . $scaffold_end . "\t" .
	int(($scaffold_end - $scaffold_start) * 0.99) . "\t" . ($scaffold_end - $scaffold_start) . "\t" .
	"60" . "\t" . 
	"tp:A:P" . "\t". "cm:i:0" . "\t" . "s1:i:0"	. "\t" . "s2:i:0" . 
	"\n";

    close(OUT);
    $sum_contig_length += $contig_len;
    $tiling_length += ($contig_end - $contig_start);
    if($prev_scaffold_end > 0 && $prev_scaffold_end < $scaffold_start){
	$gap_length += ($scaffold_start - $prev_scaffold_end);
    }

    print STDERR "Gap length: $gap_length\n";
    print STDERR "Tiling length: $tiling_length \t $sum_contig_length\n";
    print STDERR "Fraction gap " . ( $gap_length / ($gap_length + $tiling_length) ) . "\n";
    close(OUT_G);
    #<STDIN>;
    #exit(0);
}

sub get_contig_in_species{
    my ($contig_to_keep, $scaffold_dir) = @_;

    my $scaffold_info = "$scaffold_dir/../../../../scaffold_info.txt";
    my $cluster_file = "$scaffold_dir/../../../reference_mapping/clusters_seq_similarity";
    my @cluster_list = ();
    #Get the species of the scaffold
    #my $genome_path = `grep -w $scaffold_id $scaffold_info | cut -f 6`;chop $genome_path;
    
    #No psecies identified for that cluster
    if($genome_path ne "NA"){
	@tmp = split(/\//, $genome_path);
	$strain = @tmp[@tmp-1];
	@tmp = split(/_/, $strain);
	$species = $tmp[0] . "_" . $tmp[1];
	
	$all_cluster = `grep  $species $scaffold_info | cut -f 4`;
	my @cluster_list = split(/\n/, $all_cluster);
	my %cluster_to_keep = ();
	foreach $c (@cluster_list){
	    $cluster_to_keep{$c} = 1;
	}
	
	#Get the contig from the cluster set
	#print $out_dir . "\n";
	open(FILE_S, $cluster_file);
	while(<FILE_S>){
	    @line = split(/\t/, $_);
	    $cluster = $line[1];
	    $contig = $line[0];
	    if(exists $cluster_to_keep{$cluster}){
		$contig_to_keep->{$contig} = 1;
	    }
	}
	close(FILE_S);
    }
}

sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return . "\n" if($run);
    return $return;
}



