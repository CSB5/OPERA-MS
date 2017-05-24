#!/usr/bin/perl
use warnings;
use Statistics::Basic qw(:all);
#$use Math::NumberCruncher;

my ($covfile, $cluster_file, $edge_file_dir, $diff_cov_threshold, $dir_res) = @ARGV;
my $command;
my %contig_info = ();

my $WINDOW_EXCLUSION_THRESHOLD = 1.5;

open FILE, "$covfile" or die $!;
my $header = <FILE>;chop $header;
my @line = split (/ /, $header);
my $window_size = $line[3];
my ($mean, $stddev, $median, $nb_exluded_window);
#my (@window_count, @window_count_filtered);
while (<FILE>) {
    chomp $_;	
    #The line with the contig ID
    @line = split (/\t/, $_);
    $contig = $line[0];
    $length = $line[1];
    $nb_window = $line[4];
    #The line with number of arriving reads
    $read_count = <FILE>;chop $read_count;
    #print STDERR $read_count."\t|".$nb_window."|\t".$window_size."\n";<STDIN>;
    #Skip the next line that contian the windows (need to compute the variance latter)
    $str = <FILE>;chop $str;my @window_count = split(/ /, $str);
    
    $mean   = mean(  @window_count); # array refs are ok too
    $median   = median(  @window_count); # array refs are ok too
    $stddev   = stddev( @window_count );
    $variance = variance( @window_count );
    
    #print STDERR " *** $contig\n @window_count\n$median $mean $stddev\n";<STDIN>;
    $contig_info{$contig} = {
	"COV", ($read_count/($nb_window*$window_size)), 
	"LENGTH", $length, 
	#
	#
	"MEAN", $mean,
	"MEDIAN", $median,
	"STDDEV", $stddev,
	"VARIANCE", $variance,
	#
	
    };
}
close(FILE);

open(FILE, $cluster_file);
my %cluster_info;
#
while(<FILE>){
    chomp $_;	
    @line = split(/\t/, $_);
    $contig = $line[0];
    $contig_cov = $contig_info{$contig}->{"COV"};
    $cluster = $line[1];
    $cluster_mean = $line[3];
    
    if(! exists $cluster_info{$cluster}){
	$cluster_info{$cluster} = {"LENGTH", 0, "NB_ELEMENT", 0, "MEAN", $cluster_mean, "MEDIAN", 0, "COV_LIST", [], "CONTIG_LIST", [], "NB_REPEAT", 0} 
    }
    push(@{$cluster_info{$cluster}->{"COV_LIST"}}, $contig_cov);
    push(@{$cluster_info{$cluster}->{"CONTIG_LIST"}}, $contig);
}
close(FILE);

my $median_cov;
my %repeat = ();
foreach $cluster (keys %cluster_info){
    $median_cov = median(@{$cluster_info{$cluster}->{"COV_LIST"}});
    for (my $i = 0; $i < @{$cluster_info{$cluster}->{"COV_LIST"}}; $i++){
	$cov = $cluster_info{$cluster}->{"COV_LIST"}->[$i];
	$contig = $cluster_info{$cluster}->{"CONTIG_LIST"}->[$i];
	if($cov > $median_cov * $diff_cov_threshold){
	    $cluster_info{$cluster}->{"NB_REPEAT"}++;
	    $repeat{$contig} = $cluster;
	    #print STDERR $cluster."\t".$median_cov."\t".$contig."\t".$cov."\n";
	}
	#else{}
    }
}

print STDERR " *** ".(keys(%repeat) +0)." repeats filtered\n"; 

my %lib_to_filter = (
    #"results/clustersInfo_opera-lr", "illumina_edge_filter", 
    "filtered_pairedEdges_i0", "filtered_pairedEdges_i0",
    "filtered_pairedEdges_i1", "filtered_pairedEdges_i1",
    "filtered_pairedEdges_i2", "filtered_pairedEdges_i2",
    "filtered_pairedEdges_i3", "filtered_pairedEdges_i3",
    "filtered_pairedEdges_i4", "filtered_pairedEdges_i4",
    "filtered_pairedEdges_i5",  "filtered_pairedEdges_i5",
    );


print " *** Filtering edge file\n";
`mkdir $dir_res` unless -d $dir_res;
foreach $lib (keys %lib_to_filter){
    open(FILE, "$edge_file_dir/$lib");
    open(OUT, ">$dir_res/$lib_to_filter{$lib}");
    $first = 1;
    print STDERR " *** Processing lib $lib\n";#<STDIN>;
    while(<FILE>){
	if($first){
	    print OUT $_;
	    $first = 0;
	}
	else{
	    @line = split(/\s+/, $_);
	    next if(@line < 3);
	    $c1 = $line[0];
	    $c2 = $line[2];
	    #
	    #print STDERR " *** $c1\t$c2\n";<STDIN>;
	    if(exists $repeat{$c1} || exists $repeat{$c2}){
		#print STDERR " *** Filter edge $c1 $c2\n";
	    }
	    else{
		print OUT $_;
	    }
	}
    }
    close(FILE);
    close(OUT);
}

