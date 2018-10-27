#!/usr/bin/perl
#use strict;
use warnings;

my ($working_dir, $contig_file, $read_file, $nb_process, $opera_bin_dir, $racon_dir,  $minimap2_dir, $mummer_dir) = @ARGV;


run_exe("mkdir $working_dir") if(! -d $working_dir);

#split the scaffold
my $assembly_scaffold_file = "$working_dir/../scaffolds.scaf";

#print STDERR " *** $assembly_scaffold_file\n";
my %scaff_info = ();
my %singleton_contig = ();
my @scaffold_order = ();
split_scaffold_file($working_dir, $assembly_scaffold_file);

################oNLY FOR TESTIN PURPOSE
#merge_contig_from_scaf($working_dir, $contig_file, "$working_dir/../scaffoldSeq.fasta.filled");
#exit(0);
######################################3

#Get the reads that map to a single contig
my $edge_read_info_file = "$working_dir/../../long-read-mapping/edge_read_info.dat";
my $long_read_mapping = "$working_dir/../../long-read-mapping/opera.map.sort.status";

get_read_on_single_contig($working_dir, $long_read_mapping);


opendir(DIR, $working_dir);
@all_dir = readdir(DIR);
close(DIR);
my $nb_tiling_run = 0;
my $g_cmd = "$working_dir/gapfilling_cmd.sh";
open(OUT, ">$g_cmd");
foreach $scaff_dir (@all_dir){
    #next if($scaff_dir ne "opera_scaffold_1");
    next if($scaff_dir eq ".." || ! -e "$working_dir/$scaff_dir/scaffolds.scaf");
    
    print OUT "perl $opera_bin_dir/gapfill_single_scaffold.pl $contig_file $working_dir/$scaff_dir/ $edge_read_info_file $read_file $opera_bin_dir $racon_dir $minimap2_dir $mummer_dir 2> $working_dir/$scaff_dir/log.err\n";
    
    #last if($nb_tiling_run == 10);
    #$nb_tiling_run++;
    #last;
}
close(OUT);
$g_log = "$g_cmd.log";
run_exe("cat $g_cmd | xargs -L 1 -P $nb_process -I COMMAND sh -c \"COMMAND\" 2> $g_log");


#merge the contigs into the gapfilled file
merge_contig_from_scaf($working_dir, $contig_file, "$working_dir/../scaffoldSeq.fasta.filled");
    
sub merge_contig_from_scaf{
    my ($gap_filling_dir, $contig_file, $out_file) = @_;

    my %singleton_contig_seq = ();
    open(FILE, $contig_file);
    while(<FILE>){
	if($_ =~ /^>(.+)/){
	    $contig_full_name = $1;
	    @tmp = split(/\s+/, $contig_full_name);
	    $contig_name = $tmp[0];
	    #print STDERR $contig_name . "\n";# <STDIN>;
	    if(exists $singleton_contig{$contig_name}){
		$seq = <FILE>;
		$singleton_contig_seq{$contig_name} = $seq;
		#print STDERR " *** Print singleton contig " . $contig_name . "\n"; #<STDIN>;
	    }
	    else{
		<FILE>;#This contig is in a scaffold
	    }
	}
	
    }
    close(FILE);
    
    
    open(OUT, ">$out_file");
    foreach $scaff_name (@scaffold_order){
	#print STDERR " *** Print scaff_name\n";
	if($scaff_info{$scaff_name} eq "MULTI"){
	    open(FILE, "$gap_filling_dir/$scaff_name/racon_remapped_remapped.fasta");
	    @tmp = <FILE>;
	    close(FILE);
	    $seq = $tmp[1];
	    $l = length($seq) -1;
	    print OUT ">$scaff_name length: $l cov: 0.0\n";
	    print OUT $seq;
	}
	else{
	    #print STDERR " *** Print singleton contig " . $scaff_name . "\t" . $scaff_info{$scaff_name} . "\n";<STDIN>;
	    $seq = $singleton_contig_seq{$scaff_info{$scaff_name}};
	    $l = length($seq) -1;
	    print OUT ">$scaff_name length: $l cov: 0.0\n";
	    print OUT $singleton_contig_seq{$scaff_info{$scaff_name}};
	}
    }
    close(OUT);
}

sub split_scaffold_file{
    my ($working_dir, $assembly_scaffold_file) = @_;
    open(FILE, "$assembly_scaffold_file");
    my $nb_contig = 0;
    my $str_scaff = "";
    my $scaff_name = "";
    while(<FILE>){
	if($_ =~ /^>(.+)/){
	    if($nb_contig > 0){
		if($nb_contig > 1){
		    $scaff_dir = "$working_dir/$scaff_name";
		    run_exe("mkdir $scaff_dir") if(! -d $scaff_dir);
		    open(OUT, ">$scaff_dir/scaffolds.scaf");
		    print OUT $str_scaff;
		    close(OUT);
		    $scaff_info{$scaff_name} = "MULTI";
		}
		else{
		    @tmp = split(/\n/, $str_scaff);
		    @c_info = split(/\s+/, $tmp[1]);
		    $contig_name = $c_info[0];
		    #print STDERR $scaff_name . "\t" . $contig_name . "\n";<STDIN>;
		    $singleton_contig{$contig_name} = 1;
		    $scaff_info{$scaff_name} = $contig_name;
		}
	    }
	    $str_scaff = $_;
	    @tmp = split(/\s+/, $1);
	    $scaff_name = $tmp[0];
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
    if($nb_contig > 1){
	run_exe("mkdir $scaff_dir") if(! -d $scaff_dir);
	open(OUT, ">$scaff_dir/scaffolds.scaf");
	print OUT $str_scaff;
	close(OUT);
	$scaff_info{$scaff_name} = "MULTI";
    }
    else{
	@tmp = split(/\n/, $str_scaff);
	@c_info = split(/\s+/, $tmp[1]);
	$contig_name = $c_info[0];
	#print STDERR $scaff_name . "\t" . $contig_name . "\n";<STDIN>;
	$singleton_contig{$contig_name} = 1;
	$scaff_info{$scaff_name} = $contig_name;
    }
}


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




sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return . "\n" if($run);
    return $return;
}



