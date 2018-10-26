#!/usr/bin/perl
use Data::Dumper;
use warnings;
use Statistics::Basic qw(:all);
use strict;
use List::Util qw(sum);

my $inter_dir = $ARGV[0];
my $matrix = $ARGV[1];
my $out_dir = $ARGV[2];
my $contigs_file = $ARGV[3];
my $thresh1 = $ARGV[4];
my $thresh2 = $ARGV[5];
my %species_contigs = ();

open (FILE, $matrix) or die;
while(<FILE>){
    chomp $_;
    my @line = split(/\t/, $_);
    my $contig = $line[0];

    if($line[1] > $thresh1 and $line[1] < $thresh2){
        $species_contigs{$contig} = 1;
    }
}
close(FILE);

my $cmd = "mkdir $out_dir/pairedEdges_i0 $out_dir/pairedEdges_i1 $out_dir/pairedEdges_i2 $out_dir/pairedEdges_i3 $out_dir/pairedEdges_i4 $out_dir/pairedEdges_i5";
`$cmd`;
open (OUTMAP, ">", "$out_dir/contigs_340_80");
open(OUT_0, , ">", "$out_dir/pairedEdges_i0/pairedEdges_i0");
open(OUT_1, , ">", "$out_dir/pairedEdges_i1/pairedEdges_i1");
open(OUT_2, , ">", "$out_dir/pairedEdges_i2/pairedEdges_i2");
open(OUT_3, , ">", "$out_dir/pairedEdges_i3/pairedEdges_i3");
open(OUT_4, , ">", "$out_dir/pairedEdges_i4/pairedEdges_i4");
open(OUT_5, , ">", "$out_dir/pairedEdges_i5/pairedEdges_i5");
open(CONTIGS, ">", "$out_dir/contigs.fa");
open (WRONG, ">", "$out_dir/wrong");

#my @edge_files = ("$inter_dir/sigma/filtered_pairedEdges_i0","$inter_dir/sigma/filtered_pairedEdges_i1","$inter_dir/sigma/filtered_pairedEdges_i2","$inter_dir/sigma/filtered_pairedEdges_i3","$inter_dir/sigma/filtered_pairedEdges_i4","$inter_dir/sigma/filtered_pairedEdges_i5", "$inter_dir/reference_mapping/edges_between_clusters_good");

my @edge_files = ("$inter_dir/reference_mapping/edges_between_clusters_good", "$inter_dir/long-read-mapping/pairedEdges");

open (MAP, "$inter_dir/coverage_estimation/contigs_340_80");
my $header = <MAP>;
print OUTMAP $header;
while (<MAP>){
    chomp $_;
    #print STDERR $_ . "\n";
    my @line = split(/\t/, $_);
    my $contig = $line[0];
    #print STDERR $contig . "\n";
    if(exists $species_contigs{$contig}){
        print OUTMAP $_ . "\n";
        my $line1 = <MAP>;
        my $line2 = <MAP>;
        print OUTMAP $line1;
        print OUTMAP $line2;
    }
    
    else{
        <MAP>;
        <MAP>;
    }
}
close(MAP);
close(OUTMAP);

foreach my $edge_file (@edge_files){
    #print STDERR $edge_file;
    open(EDGE, $edge_file) or die;
    while(<EDGE>){
        chomp $_;
        my @line = split(/\t/,$_);
        my $contig1 = $line[0];
        my $contig2 = $line[2];
        my $support = $line[6];

        if($edge_file =~ /pairedEdges/){
            if ($support < 2){
                next;
            }
        }

        else{
            if ($support != 1){
                next;
            }
        }

        if (exists $species_contigs{$contig1} and
            exists $species_contigs{$contig2}){
            
            my $length = $line[4];

            if($length < 300){
                print OUT_0 $_ . "\n";
            }

            elsif($length < 1000){
                print OUT_1 $_ . "\n";
            }

            elsif($length < 2000){
                print OUT_2 $_ . "\n";
            }

            elsif($length < 5000){
                print OUT_3 $_ . "\n";
            }

            elsif($length < 15000){
                print OUT_4 $_ . "\n";
            }

            elsif($length < 40000){
                print OUT_5 $_ . "\n";
            }
        }

        else{
            print WRONG $_ . "\n";
        }
    }
    close(EDGE);
}
close(OUT_0);
close(OUT_1);
close(OUT_2);
close(OUT_3);
close(OUT_4);
close(OUT_5);
close(WRONG);

my $prev_scaf;
my $seq = "";
open(FILE, $contigs_file) or die;
while (<FILE>){
    chomp $_;
    if ($_ =~ />/){
        if(defined $prev_scaf){
            my @line = split(/\s/, $prev_scaf);
            my $contigs = substr $line[0], 1;

            if(exists $species_contigs{$contigs}){
                print CONTIGS $prev_scaf . "\n";
                print CONTIGS $seq . "\n";
            }
        }

        $prev_scaf = $_;
        $seq = "";
    }

    else{
        $seq .= $_;
    }
}
close(FILE);

