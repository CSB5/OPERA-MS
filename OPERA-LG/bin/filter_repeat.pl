#!/usr/bin/perl
use warnings;
use Getopt::Long;

my ($bam_file, $repeat_file) = @ARGV;

my %repeat = ();
open(FILE, $repeat_file);
while(<FILE>){
    @line = split(/\t/, $_);
    $repeat{$line[0]} = 1;
}
close(FILE);

open(FILE, "samtools view -h $bam_file |");
while(<FILE>){
    @line = split(/\t/, $_);
    $contig = $line[2];
    if(! exists $repeat{$contig}){print $_;}
}
close(FILE);

