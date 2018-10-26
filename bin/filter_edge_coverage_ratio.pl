#!/usr/bin/perl
use Statistics::Basic qw(:all);
use warnings;
use strict;

my $output_lr_directory = $ARGV[0];
my $contigs_windows_file = $ARGV[1];

my %contig_cov = ();

open (FILE, $contigs_windows_file) or die;
my $header = <FILE>;chop $header;
my @line = split (/ /, $header);
my $window_size = $line[3];
my $WINDOW_EXCLUSION_THRESHOLD=1.5;
my $edge_std_threshold = 1.5;
while (<FILE>) {
    chomp $_;	
    #The line with the contig ID
    my @line = split (/\t/, $_);
    my $contig = $line[0];
    my $length = $line[1];
    my $nb_window = $line[4];
    #The line with number of arriving reads
    my $read_count = <FILE>;chop $read_count;
    #print STDERR $read_count."\t|".$nb_window."|\t".$window_size."\n";<STDIN>;
    #Skip the next line that contian the windows (need to compute the variance latter)
    my $str = <FILE>;chop $str;
    $contig_cov{$contig} = $read_count/($nb_window*$window_size);
}

for (my $i = 0; $i <=5; $i++){
    my $old_paired_edges_file = "$output_lr_directory/prefiltered_pairedEdges_i$i";
    my $filtered_paired_edges_file = "$output_lr_directory/pairedEdges_i$i"; 
    my $below_thresh_paired_edges_file = "$output_lr_directory/below_thresh_pairedEdges_i$i"; 
    my $rename_cmd = "mv $filtered_paired_edges_file $old_paired_edges_file";
    run_exe($rename_cmd);
    
    my @edges = ();
    my @good_edges = ();
    my @bad_edges = ();
    my @all_edge_ratios = ();

    open (FILE, $old_paired_edges_file) or die;
    while (<FILE>){
       chop $_;
       my @line = split(/\t/, $_);
        
       if(!exists( $contig_cov{$line[0]}) or !exists ($contig_cov{$line[2]})){
           print STDERR $line[0] . " " . $line[2] . " DOESN'T EXIST\n";
           next;
       }
       my $supp = $line[6];

       if($supp == 1){
       #Dont' push into all_edge_ratios, that's
       #for calculating the statistics.
       push (@edges, $_);
       next;
       }
       my $cov1 = $contig_cov{$line[0]};
       my $cov2 = $contig_cov{$line[2]};
       my $cov_to_print;

       if ($cov1 > $cov2){
            $cov_to_print = $cov2;
       }

       else{
            $cov_to_print = $cov1;
       }
       push (@all_edge_ratios, log($supp/($cov_to_print + 0.0001)));
       push (@edges, $_);
    }
    close (FILE);

    my $cutoff;
    open (OUT, ">", $filtered_paired_edges_file) or die;
    open (OUT_FAIL,">", $below_thresh_paired_edges_file) or die;

    my $median_edge_ratio = median(@all_edge_ratios);
    my $std = stddev(@all_edge_ratios);

    #Edges file is empty, move on.
    if(!@edges){
        next;
    }

    $cutoff = $median_edge_ratio - $edge_std_threshold  * $std;
    
    foreach my $edge (@edges){
       my @line = split(/\t/, $edge);
       my $cov1 = $contig_cov{$line[0]};
       my $cov2 = $contig_cov{$line[2]};
       my $cov_to_print;
       my $supp = $line[6];

       if ($cov1 > $cov2){
            $cov_to_print = $cov2;
       }

       else{
            $cov_to_print = $cov1;
       }
       
       if (log($supp/($cov_to_print + 0.0001)) > $cutoff and $supp > 1){
            print OUT "$edge\n";
       }

       else{
            print OUT_FAIL "$edge\n";
       }
    }

    close(OUT);
    close(OUT_FAIL);
}

my $cmd = "cat $output_lr_directory/pairedEdges_i0 $output_lr_directory/pairedEdges_i1 $output_lr_directory/pairedEdges_i2  $output_lr_directory/pairedEdges_i3 $output_lr_directory/pairedEdges_i4 $output_lr_directory/pairedEdges_i5 > $output_lr_directory/pairedEdges_above_thresh";
run_exe ($cmd);

sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return . "\n" if($run);
    return $return;
}


