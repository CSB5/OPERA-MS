#!/usr/bin/perl
#use strict;
use warnings;

my ($contig_file, $scaff_file, $edge_file_dir, $contig_info_file, $out_dir) = @ARGV;

my $MIN_LENGTH_THRESHOLD = 1000;
my $SUPPORT_THRESHOLD = 1;
my $NB_CONTIG_THRESHOLD = 0;

#Read the scaffold file
my %scaffold_info = ();
my %contig_to_scaffold = ();
read_scaffold_file($scaff_file, \%scaffold_info, \%contig_to_scaffold);


#Read the file giving the scaffold to contig corresponding name
my %scaff_name_to_contig_name;
my %contig_name_to_scaff_name;
read_scaff_to_contig_name_file("$scaff_file.cname", \%scaff_name_to_contig_name, \%contig_name_to_scaff_name);

#Update scaffold info the contig info file
update_scaffold_info($contig_info_file, \%contig_name_to_scaff_name, \%scaffold_info);

#Write the list of circular contigs information
write_circular_contig_info($edge_file_dir, \%scaffold_info, \%contig_to_scaffold, \%scaff_name_to_contig_name, "$out_dir/circular_contig_info.dat");

sub update_scaffold_info{
    my ($contig_info_file, $contig_name_to_scaff_name, $scaffold_info) = @_;
    open(FILE, $contig_info_file);
    <FILE>;#Skip the header
    while(<FILE>){
	@line = split(/\t/, $_);
	$scaffold_name = $contig_name_to_scaff_name->{$line[0]};
	#print STDERR $scaffold_name . "\t" . $line[0] . "\n";<STDIN>;
	$length = $line[1];
	$long_read_coverage = $line[3];
	$scaffold_info->{$scaffold_name}->{"LENGTH"} = $length;
	$scaffold_info->{$scaffold_name}->{"COV"} = $long_read_coverage;
    }
    close(FILE);
}

sub write_circular_contig_info{
    my ($edge_file_dir, $scaffold_info, $contig_to_scaffold, $scaff_to_contig, $out_file) = @_;
    
    opendir(DIR, $edge_file_dir);
    my @all_file = readdir(DIR);
    close(DIR);

    open(OUT, ">$out_file");
    my ($c1, $c1_o, $c2, $c2_o, $dist, $std, $support);
    foreach $f (@all_file){
	#if(index($f, "pairedEdges_i") == 0){
	#if(index($f, "pairedEdges_i") != -1){
	if($f eq "pairedEdges"){
	    $f_p =  "$edge_file_dir/$f";
	    print STDERR " *** Read file $f_p\n";
	    open(FILE, "$f_p");
	    while(<FILE>){
		chop $_;
		@line = split(/\t/, $_);
		($c1, $c1_o, $c2, $c2_o, $dist, $std, $support) = @line;
		next if($support <= $SUPPORT_THRESHOLD);

		if(! exists $contig_to_scaffold->{$c1} || ! exists $contig_to_scaffold->{$c2}){
		    #print STDERR " Warning contig c1: $c1 not defined in scaffolds\n" if(! exists $contig_to_scaffold{$c1});
		    #print STDERR " Warning contig c2: $c2 not defined in scaffolds\n" if(! exists $contig_to_scaffold{$c2});#<STDIN>;
		    next;
		}
		
		if($contig_to_scaffold->{$c1} eq $contig_to_scaffold->{$c2}){
		    $scaffold_name = $contig_to_scaffold->{$c1};
		    #Check if the 2 last contig of the scaffold have an edge
		    #print STDERR " **** " . $scaffold_name . "\t" . $scaffold_info{$scaffold_name}->{"CONTIG"}->[0] . "\t" .  $scaffold_info{$scaffold_name}->{"CONTIG"}->[-1] . "\n";<STDIN>;
		    if(
			(($scaffold_info->{$scaffold_name}->{"CONTIG"}->[0] eq $c1 && $scaffold_info{$scaffold_name}->{"CONTIG"}->[-1] eq $c2 &&
			  $scaffold_info->{$scaffold_name}->{"CONTIG_ORIENTATION"}->[0] ne $c1_o && $scaffold_info{$scaffold_name}->{"CONTIG_ORIENTATION"}->[-1] ne $c2_o) 
			 ||
			 ($scaffold_info->{$scaffold_name}->{"CONTIG"}->[0] eq $c2 && $scaffold_info{$scaffold_name}->{"CONTIG"}->[-1] eq $c1 &&
			  $scaffold_info->{$scaffold_name}->{"CONTIG_ORIENTATION"}->[0] ne reverse_ori($c2_o) && $scaffold_info{$scaffold_name}->{"CONTIG_ORIENTATION"}->[-1] ne reverse_ori($c1_o)) 
			) 
			&&
			$scaffold_info->{$scaffold_name}->{"LENGTH"} > $MIN_LENGTH_THRESHOLD
			&&
			$scaffold_info->{$scaffold_name}->{"NB_CONTIG"} > $NB_CONTIG_THRESHOLD
			){
			print OUT $scaff_to_contig->{$scaffold_name} . "\t" . $scaffold_info->{$scaffold_name}->{"LENGTH"} . "\t" . $scaffold_info->{$scaffold_name}->{"COV"} . "\t" . $scaffold_info->{$scaffold_name}->{"NB_CONTIG"} . "\t" . $c1 . "\t" . $c1_o . "\t" . $c2 . "\t" . $c2_o . "\t" . $dist . "\n";
		    }
		}
	    }
	    close(FILE);
	}
    }
    close(OUT);
}


sub read_scaffold_file{
    my ($scaff_file, $scaffold_info, $contig_to_scaffold) = @_;
    open(FILE, $scaff_file);
    my $scaffold_name = "";
    print STDERR " *** Read scaffold file\n";
    while(<FILE>){
	chop $_;
	my @line = split(/\s+/, $_);
	if($line[0] =~ m/^>(.+)/){
	    $scaffold_name = $1;
	    $scaffold_info->{$scaffold_name} = {"LENGTH", 0, "COV", 0, "NB_CONTIG", 0, "CONTIG", [], "CONTIG_ORIENTATION", []};
	}
	else{
	    $contig = $line[0];
	    $contig_to_scaffold->{$contig} = $scaffold_name;
	    #print STDERR " *** " . $scaffold_info{$scaffold_name}->{ . 
	    push(@{$scaffold_info{$scaffold_name}->{"CONTIG"}}, $contig);

	    #For the orientation
	    $ori = "+";
	    $ori = "-" if($line[1] eq "EB");;
	    push(@{$scaffold_info{$scaffold_name}->{"CONTIG_ORIENTATION"}}, $ori);
	    $scaffold_info{$scaffold_name}->{"NB_CONTIG"}++;
	}
	
    }
    close(FILE);
}

sub read_scaff_to_contig_name_file{
    my ($scaff_to_contig_file, $scaff_to_contig, $contig_to_scaff) = @_;
    open(FILE, $scaff_to_contig_file);
    while(<FILE>){
	chop $_;
	my ($scaff, $contig) = split(/\t/, $_);
	$scaff_to_contig->{$scaff} = $contig;
	$contig_to_scaff->{$contig} = $scaff;
    }
    close(FILE);
}

sub reverse_ori{
    my ($o) = @_;
    my $res = "+";
    $res = "-" if($o eq "+");
    return $res;
}
