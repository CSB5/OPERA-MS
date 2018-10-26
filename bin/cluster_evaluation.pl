#!/usr/bin/perl
use warnings "all";
use Statistics::Basic qw(:all);

my ($coverage_contig_file, $cluster_file, $out_file) = @ARGV;


my %contig_info = ();



#print STDERR " **** $refname $quast_mapping_file\n";

open FILE, "$coverage_contig_file" or die $!;
my $header = <FILE>;chop $header;
my @line = split (/ /, $header);
my $window_size = $line[3];
my ( $mean, $stddev, $median, $nb_exluded_window);
my $WINDOW_EXCLUSION_THRESHOLD = 1.5;
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
    #print STDERR $contig."\t".$length."\t"."@window_count"."\n";<STDIN>;
    $mean   = mean(  @window_count); # array refs are ok too
    $median   = median(  @window_count); # array refs are ok too
    $stddev   = stddev( @window_count );
    $variance = variance (@window_count); 
    
    #print STDERR " *** $contig\n @window_count\n$median $mean $stddev\n";<STDIN>;
    $contig_info{$contig} = {
	"COV", ($read_count/($nb_window*$window_size)), 
	"LENGTH", $length, 
	#
	"WINDOW", \@window_count,
	#
	"MEAN", $mean,
	"MEDIAN", $median,
	"STDDEV", $stddev,
	"VARIANCE", $variance,
	#
	"#EXCLUDED_WINDOW", $nb_exluded_window,
    };
}
close(FILE);

open(FILE, $cluster_file);
my %cluster_info;
my $outlier_length_sum = 0;
my $length_sum = 0;
while(<FILE>){
    chomp $_;	
    @line = split(/\t/, $_);
    $contig = $line[0];
    $cluster = $line[1];
    $cluster_mean = $line[3];
    
    $contig_info{$contig}->{"CLUSTER"} = $cluster;
    $contig_info{$contig}->{"CLUSTER_MEAN"} = $cluster_mean;
    
    $cluster_info{$cluster} = {"LENGTH", 0, "NB_ELEMENT", 0, "NB_OUTLIER", 0,  "OUTLIER_LENGTH", 0, "MEAN", $cluster_mean, "NB_ELEMENT_GENOME", 0, "COV_LIST", [], "COV_LIST_FINAL", []} if(! exists $cluster_info{$cluster});
    
    $cluster_info{$cluster}->{"NB_ELEMENT"}++;
    if( $cluster_mean/2 > $contig_info{$contig}->{"COV"} || $contig_info{$contig}->{"COV"} > $cluster_mean * 2){
	$cluster_info{$cluster}->{"NB_OUTLIER"}++;
	$cluster_info{$cluster}->{"OUTLIER_LENGTH"} += $contig_info{$contig}->{"LENGTH"};
	$outlier_length_sum += $contig_info{$contig}->{"LENGTH"};
    }
    $cluster_info{$cluster}->{"LENGTH"} += $contig_info{$contig}->{"LENGTH"};
    #
    $length_sum += $contig_info{$contig}->{"LENGTH"};
}
close(FILE);

#Identify weird cluster to re-run sigma

#Average
#print STDERR " ####### outlier fraction\n";
#print STDERR $length_sum . "\t" . $outlier_length_sum . "\t" .  ( $outlier_length_sum / $length_sum) . "\n";

#Identify weird cluster to re-run sigma

my $max_length = -1;
my $MIN_LENGTH_THRESHOLD = 200000;
my ($max_cluster_outlier_fraction, $max_cluster_outlier_fraction_length);
$max_cluster_outlier_fraction = 0;
print STDERR "\n\n########## SORTED BY LENGTH\n";
foreach $c (sort { $cluster_info{$b}->{"LENGTH"} <=> $cluster_info{$a}->{"LENGTH"} } keys %cluster_info){
    $cluster_length = $cluster_info{$c}->{"LENGTH"};
    $cluster_outlier_fraction = ($cluster_info{$c}->{"OUTLIER_LENGTH"}/ $cluster_info{$c}->{"LENGTH"});
    $max_length =  $cluster_info{$c}->{"LENGTH"} if($max_length == -1);
    
    print STDERR $c . "\t" . $cluster_length . "\t" . $cluster_info{$c}->{"OUTLIER_LENGTH"}  . "\t" . $cluster_outlier_fraction . "\t" .
	$cluster_info{$c}->{"NB_OUTLIER"} . "\t" . $cluster_info{$c}->{"NB_ELEMENT"}  . "\t" . ($cluster_info{$c}->{"NB_OUTLIER"}/ $cluster_info{$c}->{"NB_ELEMENT"}) . 
	"\n";

    if($cluster_outlier_fraction > $max_cluster_outlier_fraction){
	$max_cluster_outlier_fraction = $cluster_outlier_fraction;
	$max_cluster_outlier_fraction_length = $cluster_length;
    }
    
    last if($cluster_length < $MIN_LENGTH_THRESHOLD);
}
  
my $SIZE_THRESHOLD = 100000;
print STDERR "\n\n########## SORTED BY OUTLIER_LENGTH\n";
my $cmp = 5;
foreach $c (sort { $cluster_info{$b}->{"OUTLIER_LENGTH"} <=> $cluster_info{$a}->{"OUTLIER_LENGTH"} } keys %cluster_info){
    print STDERR $c . "\t" . $cluster_info{$c}->{"LENGTH"} . "\t" . $cluster_info{$c}->{"OUTLIER_LENGTH"}  . "\t" . ($cluster_info{$c}->{"OUTLIER_LENGTH"}/ $cluster_info{$c}->{"LENGTH"}) . "\t" .
	$cluster_info{$c}->{"NB_OUTLIER"} . "\t" . $cluster_info{$c}->{"NB_ELEMENT"}  . "\t" . ($cluster_info{$c}->{"NB_OUTLIER"}/ $cluster_info{$c}->{"NB_ELEMENT"}) . 
	"\n";
    $cmp--;
    last if($cmp == 0);
}

print STDERR "\n" . $length_sum . "\t" . $outlier_length_sum . "\t" . $max_length . "\t" . $max_cluster_outlier_fraction . "\t" . $max_cluster_outlier_fraction_length . "\n";
