#!/usr/bin/perl
use warnings;
use Statistics::Basic qw(:all);
#$use Math::NumberCruncher;

my ($covfile, $cluster_file, $edge_file_dir, $diff_cov_threshold, $dir_res, $nucmer_dir, $mummer_exe_dir) = @ARGV;
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
my %cluster_info = ();
my %contig_in_cluster = ();
#
while(<FILE>){
    chomp $_;	
    @line = split(/\t/, $_);
    $contig = $line[0];
    $contig_cov = $contig_info{$contig}->{"COV"};
    $cluster = $line[1];
    #print STDERR "$cluster\n";
    $cluster_mean = $line[3];
    
    if(! exists $cluster_info{$cluster}){
	$cluster_info{$cluster} = {"LENGTH", 0, "NB_ELEMENT", 0, "MEAN", $cluster_mean, "MEDIAN", 0, "COV_LIST", [], "CONTIG_LIST", [], "NB_REPEAT", 0} 
    }
    push(@{$cluster_info{$cluster}->{"COV_LIST"}}, $contig_cov);
    push(@{$cluster_info{$cluster}->{"CONTIG_LIST"}}, $contig);
    $cluster_info{$cluster} -> {"LENGTH"} += $contig_info{$contig} -> {"LENGTH"};
    $cluster_info{$cluster} -> {"NB_ELEMENT"}++;
    $contig_in_cluster{$contig} = 1;
}
close(FILE);

my $median_cov;
my %repeat = ();
foreach $cluster (keys %cluster_info){
    $median_cov = median(@{$cluster_info{$cluster}->{"COV_LIST"}});
    $cluster_info{$cluster} -> {"MEDIAN"} = $median_cov;
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

opendir(my $nuc_dir, $nucmer_dir) or die "NUCMER directory $nucmer_dir not found.\n";

#Reference genome repeat detection
#A contig is defined as repeat is is mapping length over a genome is higher than 150%
while (my $mapped_file = readdir($nuc_dir)){
    if ($mapped_file =~ /(\d+)-repeat_detection.delta/){
        my $cluster = $1;
       `$mummer_exe_dir/show-coords -lrcT $nucmer_dir/$mapped_file > $nucmer_dir/$1.txt`;
       open (NUC_MAPPING,"$nucmer_dir/$1.txt") or die;
        my %contig_mapped = ();
       <NUC_MAPPING>;
       <NUC_MAPPING>;
       <NUC_MAPPING>;
       <NUC_MAPPING>;
       while(<NUC_MAPPING>){
          if ($_ eq ""){
               next;
          }
          #print STDERR $_;
          my @line = split(/\t/, $_);
          my $contig = $line[11];
          my $percent_mapped = $line[9];

          if(exists $contig_mapped{$contig}){
              $contig_mapped{$contig} += $percent_mapped;
          }

          else{
              $contig_mapped{$contig} = $percent_mapped;
          }
       }

       foreach my $contig( keys %contig_mapped){
            if ($contig_mapped{$contig} > 150){
                $repeat{$contig} = $1;
                $cluster_info{$cluster}->{"NB_REPEAT"}++;
                $contig_info{$contig} -> {"REF_REPEAT"} = 1;
            }
       }
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
	    if(exists $repeat{$c1} || exists $repeat{$c2} || # remove the repeats
	       ! exists $contig_in_cluster{$c1} || ! exists $contig_in_cluster{$c2} #remove the contigs that does not belong to the current clusters (some cluster may have been filter as they assoited to species with multiple strain)
		){
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

open (REPEAT_INFO, ">", "$dir_res/repeat_contigs_info.txt") or die;
open (CLUSTERS_COV, ">", "$dir_res/clusters_cov") or die;
foreach my $cluster (keys %cluster_info){
    #print STDERR $cluster . "\n";

    if (!exists $cluster_info{$cluster} -> {"MEDIAN"}){
        next;
    }

    my $cluster_med = $cluster_info{$cluster}->{"MEDIAN"};
    my $cluster_length = $cluster_info{$cluster}->{"LENGTH"};
    my $cluster_nb_ele = ${cluster_info{$cluster}->{"NB_ELEMENT"}};
    print REPEAT_INFO "$cluster $cluster_med $cluster_length $cluster_nb_ele\n";
    print CLUSTERS_COV "$cluster $cluster_med\n";
    foreach my $contig ( @{$cluster_info{$cluster}->{"CONTIG_LIST"}}){
        if (exists $repeat{$contig}){
            my $contig_cov = ${contig_info{$contig}->{"COV"}};
            my $contig_length = ${contig_info{$contig}->{"LENGTH"}};
            print REPEAT_INFO "\t$contig\t$contig_cov\t$contig_length";
            print REPEAT_INFO "\tREPEAT" if (exists $contig_info{$contig}->{"REF_REPEAT"});#If it detected using the reference 
            print REPEAT_INFO "\n";
        }
    }
}

    
