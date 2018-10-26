#!/usr/bin/perl
use warnings;



open(NUC_MAPPING, "sort -k1,1 -n contig.map | ")
    or die "Parsing problem during read rescue using NUCMER mapping.\n";

#skip the first four lines of the cluster-vs-cluster.txt file.
<NUC_MAPPING>;
<NUC_MAPPING>;
<NUC_MAPPING>;
<NUC_MAPPING>;



my $prev_start = -1;
my $prev_end = -1;
my $prev_length = -1;
#my $prev_percent_map = -1;

while(<NUC_MAPPING>){
    if ($_ eq ""){
	next;
    }
    #print STDERR $_;
    my @line = split(/\t/, $_);
    my $contig = $line[11];
    my $percent_mapped = $line[9];
    
    
    my $index;
