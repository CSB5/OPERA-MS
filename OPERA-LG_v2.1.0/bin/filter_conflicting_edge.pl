#!/usr/bin/perl
use warnings;
use Getopt::Long;

my ($edge_file, $contig_file, $kmer_size, $cluster_size_threshold) = @ARGV;

my $min_opera_contig_size = 500;
my $std_extention = 6;
my $intersection_threshold = $kmer_size*2;

#Data structure to store the QC data
my %filtered_edge = ();
my %anchor_contig_info = ();
my @filtered_edge_cluser = ();

#get the contig coverage [ONLY FOR PLOTTING]
my %contig_coverage = ();
my %contig_length = ();
open(FILE, $contig_file);
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $contig_coverage{$line[0]} = $line[2];
    $contig_length{$line[0]} = $line[1];
}
close(FILE);

#collect of the contigs from the edge file
my %contig_edges = ();

open(FILE, $edge_file);
my ($c1, $c2, $c1_ID, $c2_ID);
my $nb_edge = 0;
while(<FILE>){
    chop $_;
    $str_edge = $_;
    
    #print " **** $str_edge\n";<STDIN>;

    @line = split(/\t/, $str_edge);
    
    next if(@line < 6);

    #edge distance
    $distance = $line[4];
    $std = $line[5];
    $cluster_size = $line[6];

    #Get the contig information
    $c1 = $line[0];$c1_ori = $line[1];
    $c1_ID = "$c1\_E";
    $c1_ID = "$c1\_B" if($c1_ori eq "-");
    #
    $c2 = $line[2];$c2_ori = $line[3];
    $c2_ID = "$c2\_B";
    $c2_ID = "$c2\_E" if($c2_ori eq "-");
    #

    next if($cluster_size < $cluster_size_threshold ||
	    $contig_length{$c1} < $min_opera_contig_size ||
	    $contig_length{$c2} < $min_opera_contig_size
	);
    
    $nb_edge++;

    
    
    #next if($c2 ne "633463");

    #print STDERR " *** (1) $c1_ID\t$c2_ID\n";#<STDIN>;

    #For c1
    if(!exists $contig_edges{$c1_ID}){
	$contig_edges{$c1_ID} = [];
    }
    detect_all_edge_conflict($c1_ID, $contig_length{$c1}, $contig_coverage{$c1}, $contig_length{$c2}, $contig_coverage{$c2}, $distance, $std, $cluster_size, $str_edge);
    push (@{$contig_edges{$c1_ID}}, $str_edge);
    
    #print STDERR " *** (2) $c1_ID\t$c2_ID\n";#<STDIN>;
    #For c2
    if(!exists $contig_edges{$c2_ID}){
	$contig_edges{$c2_ID} = [];
    }
    detect_all_edge_conflict($c2_ID, $contig_length{$c2}, $contig_coverage{$c2}, $contig_length{$c1}, $contig_coverage{$c1}, $distance, $std, $cluster_size, $str_edge);
    push (@{$contig_edges{$c2_ID}}, $str_edge);
    #<STDIN>;
}
close(FILE);

open(OUT, ">filtered_edges.dat");
open(OUT_C, ">filtered_edges_cov.dat");
$nb_edge_filtered = keys (%filtered_edge);
print STDERR "\n *** ".$nb_edge_filtered."  edges filtered from $nb_edge (".($nb_edge_filtered/$nb_edge).")\n\n";

foreach $e (keys %filtered_edge){
    print OUT $e."\n";
    print OUT_C $filtered_edge{$e}."\n";
}
close(OUT);
close(OUT_C);

open(OUT, ">anchor_contig_info.dat");
foreach $c (keys %anchor_contig_info){
    print OUT $anchor_contig_info{$c}."\n";
}
close(OUT);

open(OUT, ">edge_cluster_info.dat");
foreach $c (@filtered_edge_cluser){
    print OUT $c."\n";
}
close(OUT);



sub detect_all_edge_conflict{
    my ($anchor_contig, $anchor_contig_size, $anchor_contig_coverage, $new_edge_contig_length, $new_edge_contig_coverage, $new_edge_distance, $new_edge_std, $new_edge_cluster_size, $new_edge_str) = @_;
    
    #print STDERR " *** detect_all_edge_conflict $anchor_contig, $anchor_contig_coverage, $new_edge_contig_length, $new_edge_contig_coverage, $new_edge_distance, $new_edge_std, $new_edge_str|\n";

    my ($c1,$c1_ID,$c2,$c2_ID,$res, $other_contig_length, $other_distance, $other_std);

    #print STDERR " *** detect_all_edge_conflict: $anchor_contig \n";#<STDIN>;
    for(my $i = 0; $i < @{$contig_edges{$anchor_contig}}; $i++){

	#print STDERR " \t $i\n";<STDIN>;
	$other_edge_str = $contig_edges{$anchor_contig}->[$i];
	@other_edge = split(/\t/, $other_edge_str);
	
	#edge distance
	$other_cluster_size = $other_edge[6];
	$other_distance = $other_edge[4];
	$other_std = $other_edge[5];
	
	#Get the contig information
	$c1 = $other_edge[0];$c1_ori = $other_edge[1];
	$c1_ID = "$c1\_E";
	$c1_ID = "$c1\_B" if($c1_ori eq "-");
	#
	$c2 = $other_edge[2];$c2_ori = $other_edge[3];
	$c2_ID = "$c2\_B";
	$c2_ID = "$c2\_E" if($c2_ori eq "-");

	$other_contig_length = $contig_length{$c2};
	$other_contig_coverage = $contig_coverage{$c2};
	#
	if($c2_ID eq $anchor_contig){
	    $other_contig_length = $contig_length{$c1};
	    $other_contig_coverage = $contig_coverage{$c1};
	}
	
	$res = 0;
	#print STDERR " *** std |$other_std| |$new_edge_std|\n";<STDIN>;
	$res = (detect_edge_conflict($new_edge_contig_length, $other_contig_length, $new_edge_distance, $other_distance, $new_edge_std, $other_std) && 
		detect_edge_conflict($other_contig_length, $new_edge_contig_length, $other_distance, $new_edge_distance, $other_std, $new_edge_std));
	#print STDERR " *** Res test conflict $res\n";<STDIN>;
	if($res){
	    $filtered_edge{$new_edge_str} = $anchor_contig_coverage."\t".$new_edge_contig_coverage;
	    $filtered_edge{$other_edge_str} = $anchor_contig_coverage."\t".$other_contig_coverage;
	    @c = split(/_/, $anchor_contig);
	    $anchor_contig_info{$c[0]} = $c[0]."\t".$anchor_contig_coverage."\t".$anchor_contig_size;
	    
	    push(@filtered_edge_cluser, $new_edge_cluster_size."\t".$other_cluster_size) if($new_edge_cluster_size <= $other_cluster_size);
	    push(@filtered_edge_cluser, $new_edge_cluster_size."\t".$other_cluster_size) if($new_edge_cluster_size > $other_cluster_size);
	    #last;#No sure if this transitive
	}
    }

}

#c contig with multiple edge
#c1 and c2 extremity of the 2 edges
sub detect_edge_conflict{
    my ($c_e1_length, $c_e2_length, $e1_distance, $e2_distance, $e1_std, $e2_std) = @_;

    #print STDERR "detect_edge_conflict $c_e1_length, $c_e2_length, $e1_distance, $e2_distance, $e1_std, $e2_std\n";#<STDIN>;
    
    my (@c_e1_pos, @c_e2_pos, @tmp);

    #Extend with c1 with min distance
    $c_e1_pos[0] = $e1_distance - ($std_extention * $e1_std);
    $c_e1_pos[1] = $c_e1_pos[0] + $c_e1_length;
    
    #Extend with c2 with max distance
    $c_e2_pos[0] = $e2_distance + ($std_extention * $e2_std);
    $c_e2_pos[1] = $c_e2_pos[0] + $c_e2_length;
    
    #We want to make sure that c_e1 as the smallest begin coordinate
    if($c_e1_pos[0] > $c_e2_pos[0]){
	@tmp = @c_e1_pos;
	@c_e1_pos = @c_e2_pos;
	@c_e2_pos = @tmp;
    }

    my $res = 0;
    
    #print STDERR "[$c_e1_pos[0] $c_e1_pos[1]] [$c_e2_pos[0] $c_e2_pos[1]] => ".($c_e1_pos[1] - $c_e2_pos[0])."\n";

    $res = 1 if($c_e2_pos[0] <  $c_e1_pos[1]  &&  $c_e1_pos[1] - $c_e2_pos[0] > $intersection_threshold); #c_e1_b c_e2_b [c_e2_e] c_e1_e [c_e2_e]

    #print STDERR " *** Res test intersection $res\n\n";

    return $res;
} 

