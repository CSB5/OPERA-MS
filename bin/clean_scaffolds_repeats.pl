#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw(:all);

my $inter_dir = $ARGV[0];

my %contig_cov = ();
my @contigs_file = glob("$inter_dir/coverage_estimation/contigs_*");
print STDERR $contigs_file[0];
my $contigs_cov_file = $contigs_file[0]; 
open (FILE, "$contigs_cov_file") or die;
my $header = <FILE>;
chop $header;
my @line = split (/ /, $header);
my $window_size = $line[3];
while (<FILE>) {
    chomp $_;	
    my @line = split (/\t/, $_);
    my $contig = $line[0];
    my $length = $line[1];
    my $nb_window = $line[4];
    my $read_count = <FILE>;chop $read_count;
    my $str = <FILE>;chop $str;

    if ($length > 500){
        $contig_cov{$contig} = $read_count/($nb_window*$window_size);
    }
}

close FILE;

open (EDGES, "$inter_dir/long-read-mapping/edge_read_info.dat") or die;
my %contigs_contained_within = ();
while (<EDGES>){
    chomp $_;
    #print STDERR "$_\n";
    my @line = split(/\t/, $_);
    $_ =~ /CONTIG_FOR_GAPFILLING:(.*)/;
    my @comma_delim = split(/;/,$1);
    my %contigs_hash = ();
    my $contig1 = $line[0];
    my $contig2 = $line[2];

    foreach my $contigs (@comma_delim){
        my @contig_set = split(/,/, $contigs);
        foreach my $contig (@contig_set){
            if ($contig eq "NONE"){next;}
            my $stripped_contig = substr $contig,1;
            $contigs_hash{$stripped_contig} = 1;
            #print STDERR $stripped_contig . "\n";
        }
    }

    $contigs_contained_within{"$contig1-vs-$contig2"} = \%contigs_hash;
}
close(EDGES);


open (SCAFF, "$inter_dir/opera_long_read/scaffolds.scaf") or die;
my %covered_by_gapfilling = ();
my $contig1;
my $contig2;
while(<SCAFF>){
    chomp $_;
    if ($_ !~ /^>/){
        my @line = split(/\t/, $_);
        my $contig2 = $line[0];
        if ($contig1 ne ""){
            my %contigs_hash = ();

            if (exists $contigs_contained_within{"$contig1-vs-$contig2"} and
                exists $contig_cov{$contig1} and
                exists $contig_cov{$contig2}){
                %contigs_hash = %{$contigs_contained_within{"$contig1-vs-$contig2"}};
            }

            elsif (exists $contigs_contained_within{"$contig2-vs-$contig1"} and
                   exists $contig_cov{$contig1} and
                   exists $contig_cov{$contig2}){
                %contigs_hash = %{$contigs_contained_within{"$contig2-vs-$contig1"}};
            }

            else{
                #print STDERR "not gap filled $contig1 $contig2\n";
                $contig1 = $contig2;
                next;
            }

            my $cov_to_subtract;
            if($contig_cov{$contig1} < $contig_cov{$contig2}){
                $cov_to_subtract = $contig_cov{$contig2};
            }

            else{
                $cov_to_subtract = $contig_cov{$contig1};
            }

            foreach my $contig_in_gap (keys %contigs_hash){

                if (exists $contig_cov{$contig_in_gap}){
                    #Check if the two anchor contigs are repeats.
                    if ($contig_cov{$contig_in_gap} * 1.5 <  $cov_to_subtract){ 
                        next;
                    }
                }
                
                if (!exists $covered_by_gapfilling{$contig_in_gap}){
                    $covered_by_gapfilling{$contig_in_gap} = $cov_to_subtract;
                }
                else{
                    $covered_by_gapfilling{$contig_in_gap} += $cov_to_subtract;
                }
            }
        }
        #print STDERR "$contig1 $contig2\n";
        $contig1 = $contig2;
    }

    else{
        #print STDERR "SCAFFOLD $_\n";
       $contig1 = "";
       $contig2 = "";
    }
}
close (SCAFF);

open (SCAFFOUT, ">", "$inter_dir/opera_long_read/scaffolds_no_repeat.scaf") or die;
open (SCAFF, "$inter_dir/opera_long_read/scaffolds.scaf") or die;
my $current_scaff;
my %repeat_scaffs = ();
while (<SCAFF>){
    chomp $_;
    if ($_ !~ /^>/){
        my @line = split(/\t/, $_);
        my $contig = $line[0];
        my $is_repeat;

        if (exists $covered_by_gapfilling{$contig} and
            exists $contig_cov{$contig}){
            if ($covered_by_gapfilling{$contig} > 0.8 * $contig_cov{$contig}){
                $is_repeat = 1;
            }
        }

        if ($is_repeat){
            print SCAFFOUT "$_ $contig_cov{$contig} $covered_by_gapfilling{$contig} REMOVED\n";
        }

        else{
            if (exists $covered_by_gapfilling{$contig} and
                exists $contig_cov{$contig}){
                 print SCAFFOUT "$_ $contig_cov{$contig} $covered_by_gapfilling{$contig} BELOW THRESHOLD \n" ;
            }

            else{
                print SCAFFOUT "$_\n";
            }
            
            $repeat_scaffs{$current_scaff} = 0;
        }

    }

    else{
        print SCAFFOUT "$_\n";
        $current_scaff = $_;
        $repeat_scaffs{$current_scaff} = 1;
    }
}
close(SCAFF);
close(SCAFFOUT);

open (SCAFF_SEQ, "$inter_dir/opera_long_read/scaffoldSeq.fasta.filled") or die;
open (OUT, ">", "$inter_dir/../scaffoldSeq_no_repeat.fasta.filled") or die;
while(<SCAFF_SEQ>){
    chomp $_;
    if ($_ !~ /^>/){
        next;
    }

    if (!($repeat_scaffs{$_})){
        #print STDERR "$repeat_scaffs{$_}\n";
       print OUT $_ . "\n";
       my $scaffold_sequence = <SCAFF_SEQ>;
       print OUT $scaffold_sequence;
    }

    else{
        #print STDERR "$_ is repeat scaff\n";
    }
}
close(SCAFF_SEQ);

