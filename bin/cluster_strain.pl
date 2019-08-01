#!/usr/bin/perl
use warnings "all";
use Statistics::Basic qw(:all);
use Statistics::R;


my ($end_time, $start_time);
$start_time = time;

my $R = Statistics::R->new();
my @confidance_interval;


my ($species_dir, $reference_size, $dispersion_value, $step_dispersion_value_increase, $opera_ms_dir) = @ARGV;

my $coverage_file = "$species_dir/contigs_window_cov";

my %contig_info = ();

open FILE, "$coverage_file" or die $!;
my $header = <FILE>;chop $header;
#
my @line = split (/ /, $header);
my $window_size = $line[3];
my ( $mean, $stddev, $median, $nb_exluded_window);
while (<FILE>) {
    chomp $_;	
    #The line with the contig ID
    @line = split (/\t/, $_);
    $contig = $line[0];
    $length = $line[1];
    $nb_window = $line[4];
    #The line with number of arriving reads
    $read_count = <FILE>;chop $read_count;
    
    #Skip the next line that contian the windows (need to compute the variance latter)
    $str = <FILE>;chop $str;my @window_count = split(/ /, $str);
    #print STDERR $contig."\t".$length."\t"."@window_count"."\n";<STDIN>;
    $mean   = int(mean(@window_count)); # array refs are ok too
    $median = median(@window_count); # array refs are ok too
    
    #print STDERR " *** $contig\n @window_count\n$median $mean $stddev\n";<STDIN>;
    $contig_info{$contig} = {
	#"READ_COUNT", $read_count, 
	"COV", ($read_count/($nb_window*$window_size)), 
	    "LENGTH", $length, 
	    #
	    "WINDOW", \@window_count,
	    
	    #
	    "MEAN_READ_COUNT", $mean,
	    "MEDIAN_READ_COUNT", $median,
	    #"STDDEV_COV", $stddev,
	    #"VARIANCE_COV", $variance,

	    "STRAIN_ID", 0,
	    "SHARED_CONTIG", 0
    };
}
close(FILE);


#get the initial mean
my $sum_contig_length = 0;
my $sum_read_count = 0;
my @contig_mean_cov = ();
my $window_count = [];
my $min_contig_cov = 0;
#
my $MAX_CONTIG_COV = 100000;
#
foreach $contig (sort {$contig_info{$b}->{"MEAN_READ_COUNT"} <=> $contig_info{$a}->{"MEAN_READ_COUNT"}} keys %contig_info){
    #last if($sum_contig_length > $reference_size);
    $contig_read_count = $contig_info{$contig}->{"MEAN_READ_COUNT"};
    next if($contig_read_count > $MAX_CONTIG_COV);
    
    push(@contig_mean_cov,  $contig_read_count);
    push(@{$window_count}, @{$contig_info{$contig}->{"WINDOW"}});#if($contig_info{$contig}->{"MEAN_READ_COUNT"} <7000);
    $sum_read_count += $contig_read_count;
    $sum_contig_length += $contig_info{$contig}->{"LENGTH"};
    $min_contig_cov = $contig_read_count;
}

#Indicate the contig that have been slected and their assignation
open(OUT_W, ">$species_dir/window_distibution.dat");
foreach $contig (keys %contig_info){
    print OUT_W "" . join("\n", @{$contig_info{$contig}->{"WINDOW"}}) . "\n";
}
print STDERR " *** End window construction\n";
#<STDIN>;

#Identify suitable mode:
#The highest coverage mode should contain large number of contigs as it must contigs of the histes abundance species
#The mode identification can be affected by widows with outlier coverage

my ($mode, $strain_mean_value, $strain_var, $sequence_size_in_highest_mode_distribution);
my $mode_accepted = -1;
my $nb_window_in_dist;
my $coverrage_sun_in_dist;
my $window_count_in_h_mode_distribution;
my @strain_list = ();
while($mode_accepted == -1){
    
    $mode = compute_mode($window_count);

    #Select the higher mode that covers a significant number of contigs
    print STDERR " *** MODE DETECTED @{$mode}\n";

    #if(@{$mode} == 1){
    #update the smoothing parameter
    #next;
    #}

    for(my $i = @{$mode} - 1; $i >= 0; $i--){
	
	
	$curr_mode = $mode->[$i];
	next if(index($curr_mode, "[") != -1);
	
	#get the dispersion from the highest mode
	$strain_mean_value = estimate_mean($curr_mode, $dispersion_value);
	#get the confidence interval
	compute_confidance_interval($strain_mean_value, $dispersion_value);

	print STDERR "MEAN ESTIMATION BASED " . $strain_mean_value . "ON MODE " . $curr_mode . " AND DISPERTION ". $dispersion_value . "\n";
	print STDERR "ESTIMATED CONFIDANCE INTERVAL FOR  " . $curr_mode. " -> [@confidance_interval]" . "\n";#<STDIN>;
	
	#
	@{$window_count_in_h_mode_distribution} = ();
	$nb_window_in_dist = 0;$coverage_sum_in_dist = 0;
	#The coverage estimate is made on contig then windows are extracted
	foreach $contig (sort {$contig_info{$b}->{"MEAN_READ_COUNT"} <=> $contig_info{$a}->{"MEAN_READ_COUNT"}} keys %contig_info){
	    next if($contig_info{$contig}->{"STRAIN_ID"} != 0);
	    $contig_read_count = $contig_info{$contig}->{"MEAN_READ_COUNT"};
	    #$probability = compute_probability( $strain_mean_cov, $dispersion_value, $contig_mean_cov);
	    
	    if($confidance_interval[0] < $contig_read_count && $contig_read_count < $confidance_interval[1]){
		foreach $c (@{$contig_info{$contig}->{"WINDOW"}}){
		    $coverage_sum_in_dist += $c;
		    $nb_window_in_dist ++;
		}
		push(@{$window_count_in_h_mode_distribution}, @{$contig_info{$contig}->{"WINDOW"}});
	    }
	}

	next if($nb_window_in_dist == 0);#Can happen when contigs have widows with very different coverage and the dispertion parameter is very high
	
	$sequence_size_in_h_mode_distribution =  $nb_window_in_dist * $window_size;#@{$window_count_in_h_mode_distribution} * $window_size;
	$mean_cov_in_dist = $coverage_sum_in_dist / $nb_window_in_dist;
	print STDERR " *** Mode distribution evaluation : " . $curr_mode . " sequence_size_in_h_mode_distribution: $sequence_size_in_h_mode_distribution " . "mean_cov_in_dist $mean_cov_in_dist ". " estimated_mean_value $strain_mean_value " . "mean_cov_in_dist/estimated_mean_value " . ($mean_cov_in_dist  / $strain_mean_value) . "\n";#<STDIN>;
	
	if((0.70 < ($mean_cov_in_dist  / $strain_mean_value) && ($mean_cov_in_dist  / $strain_mean_value) < 1.3) && #the mean estimate and the mean obtain by the confidance interval are comparable => good fitting
	   # The sequence is big enough to be a complete genome for the first mode that should be a whole genomes, all the mode will lowest coverage value that pass the pervious test will be taken into account
	   #($mode_accepted != -1 || $sequence_size_in_h_mode_distribution > $reference_size*0.5)){
	   (($mode_accepted != -1 && $sequence_size_in_h_mode_distribution > 100000) || $sequence_size_in_h_mode_distribution > 500000)){

	    $mode_accepted = 0 if($mode_accepted == -1);
	    $mode_accepted++;
	    push(@strain_list, $mode_accepted);
	    print STDERR " *** *** *** MODE ACCEPTED $i $curr_mode\n\n";
	    
	    select_mode_contigs($mode_accepted, $curr_mode, $strain_mean_value, $window_count_in_h_mode_distribution);
	}
    }
    
    if($mode_accepted == -1){
	#All node are pulled into STRAIN_1
	$mode_accepted = 1;
	push(@strain_list, $mode_accepted);
	foreach $contig (keys %contig_info){
	    $contig_info{$contig}->{"STRAIN_ID"} = $mode_accepted;
	}
    }
}
print STDERR " *** End mode detection\n";
#<STDIN>;
#exit(0);

#filter contig base on reference mapping
#contig that share similar sequence are filetred out if they are predicted to belong to the same cluster
#filter_mapping();


#Construct the data for OPERA-LG
my @EDGE_ID_INFO = (300, 1000, 2000, 5000, 15000, 40000);
my ($species_contig_file, $strain_contig_file);
foreach $strain_id (@strain_list){
    $strain_dir = "$species_dir/STRAIN_$strain_id";
    `mkdir $strain_dir` if(! -d $strain_dir);
    
    #Generate the opera_lg config file
    my $OPERA_LG_CONFIG;
    open($OPERA_LG_CONFIG, ">$strain_dir/opera.config");
    print $OPERA_LG_CONFIG "output_folder=$strain_dir\n";
    print $OPERA_LG_CONFIG "contig_file=$strain_dir/contigs.fa\n";
    print $OPERA_LG_CONFIG 
	"keep_repeat=yes" . "\n" . 
	"filter_repeat=no" . "\n" . 
	"cluster_threshold=1" . "\n" . #How to handle the cluster threshold 
	"cluster_increased_step=5" . "\n" . 
	"kmer=60" . "\n";
    
    #get the strain edges
    for(my $edge_id = 0; $edge_id <= 5; $edge_id++){
	extract_edge($edge_id, $strain_id, $OPERA_LG_CONFIG);
    }
    close($OPERA_LG_CONFIG);
}
print STDERR " *** OPERA-LG data generation\n";
#<STDIN>;


#generate the strain contig file
foreach $strain_id (@strain_list){
    $strain_dir = "$species_dir/STRAIN_$strain_id";
    #
    $species_contig_file = "$species_dir/contigs.fa";
    open(FILE, $species_contig_file);
    $strain_contig_file = "$strain_dir/contigs.fa";
    open(OUT, ">$strain_contig_file");
    #Generate the strain contig file
    while(<FILE>){
	if($_ =~ m/>(.*)/){
	    $c = $1;@tmp = split(/\s+/, $c);$c = $tmp[0];
	    #print STDERR $c . "\n";
	    #next if(! exists $contig_info{$c});
	    #print STDERR " *** Contig name $c\n";<STDIN>;
	    $seq = <FILE>;
	    if(exists $contig_info{$c}->{"STRAIN_ID"} && $contig_info{$c}->{"STRAIN_ID"} ne "NA" && ($contig_info{$c}->{"STRAIN_ID"} == $strain_id || $contig_info{$c}->{"SHARED_CONTIG"})){
		print OUT ">$c\n";
		print OUT $seq;
	    }
	}
    }
    close(OUT);
}

#Indicate the contig that have been slected and their assignation
open(OUT, ">$species_dir/strain_cluster.dat");
open(OUT_W, ">$species_dir/window_distibution.dat");
foreach $contig (keys %contig_info){
    if(exists $contig_info{$contig}->{"STRAIN_ID"} && $contig_info{$contig}->{"STRAIN_ID"} ne "NA"){
	print OUT $contig . "\t" . $contig_info{$contig}->{"MEAN_READ_COUNT"} . "\t" . $contig_info{$contig}->{"STRAIN_ID"} . "\t" . $contig_info{$contig}->{"SHARED_CONTIG"} . "\n";
	print OUT_W "" . join("\n", @{$contig_info{$contig}->{"WINDOW"}}) . "\n";
    }
}
$end_time = time;
my @tmp = split(/\//,$species_dir);
my $species_name = $tmp[-1];
print STDOUT "*** ***  Clustering $species_name Elapsed time: " . ($end_time - $start_time) . "\n";

#Run OPERA-LG
foreach $strain_id (@strain_list){
    #next if($strain_id == 2);
    $start_time = time;
    $strain_dir = "$species_dir/STRAIN_$strain_id";
    run_exe("timeout 5m $opera_ms_dir/OPERA-LG/bin/OPERA-LG $strain_dir/opera.config  > $strain_dir/log.txt");
    if(! -e "$strain_dir/scaffoldSeq.fasta"){#NEED TO IN THE MAKE FILE THE COMMAND TO GET THE FASTS OPERA-MS
	run_exe("$opera_ms_dir/OPERA-LG/bin/OPERA-LG-fast $strain_dir/opera.config  > $strain_dir/log_fast.txt");
    }

    if( ! -e "$strain_dir/scaffoldSeq.fasta"){
	die " Error in OPERA-LG. Please see $strain_dir/log.txt and $strain_dir/log_fast.txt for details.\n";
    }
    
    $end_time = time;
    my @tmp = split(/\//,$species_dir);
    my $species_name = $tmp[-1];
    print STDOUT "*** ***  Assembly $species_name $strain_id, Elapsed time: " . ($end_time - $start_time) . "\n";
}

#Write the file with the exluded contigs
get_contig_sequence("$species_dir/excluded_contigs.fa", "$species_dir/contigs.fa", \%contig_info);

sub select_mode_contigs{
    my ($strain_id, $mode_value, $strain_mean_value, $window_count_in_h_mode_distribution) = @_;
    #Restimation of R
    #my $strain_var = int(variance ($window_count_in_h_mode_distribution));
    my $strain_var = variance ($window_count_in_h_mode_distribution);
    #
    my $dispersion_value = ($strain_mean_value * $strain_mean_value) / ($strain_var - $strain_mean_value);
    $dispersion_value = 500 if($strain_var - $strain_mean_value < 0);#Case of distribution very close to poisson => give a high dispertion value
    print STDERR "\n Strain $strain_id\n";
    print STDERR "DISPERTION VALUE " . $dispersion_value . " ESTIMATED USING " . $strain_mean_value . " AND VARIANCE |". $strain_var . "|\n";
    #Restimate the mean using the restimated R value
    $strain_mean_value = estimate_mean($mode_value, $dispersion_value);

    print STDERR "MEAN ESTIMATION BASED " . $strain_mean_value . "\t". $strain_var . "\t|" . $dispersion_value . "|\n";

    #Get the confiance interval for the strain coverage
    compute_confidance_interval($strain_mean_value, $dispersion_value);
    #
    print STDERR $strain_mean_value . "\t". $dispersion_value . "\t" . "[@confidance_interval]" . "\t" . ($sum_read_count / (@contig_mean_cov+0)) . "\t" . $min_contig_cov . "\n";#<STDIN>;

    #Score each contig and assign a strain to it
    foreach $contig (sort {$contig_info{$b}->{"MEAN_READ_COUNT"} <=> $contig_info{$a}->{"MEAN_READ_COUNT"}} keys %contig_info){
	$contig_mean_cov = $contig_info{$contig}->{"MEAN_READ_COUNT"};
	#$probability = compute_probability( $strain_mean_cov, $dispersion_value, $contig_mean_cov);

	if($contig_info{$contig}->{"STRAIN_ID"} == 0 && #Not selected in another strain
	   $confidance_interval[0] < $contig_mean_cov && $contig_mean_cov < $confidance_interval[1]){
	    $contig_info{$contig}->{"STRAIN_ID"} = $strain_id;
	}
	
    }
    
}



sub extract_edge{
    my ($edge_id, $strain_id, $OPERA_LG_CONFIG) = @_;

    my $edge_file = "$species_dir/pairedEdges_i$edge_id/pairedEdges_i$edge_id";
    #my $edge_file = "$species_dir/pairedge_i$edge_id/pairedEdges_i$edge_id";
    print STDERR " *** Analyzing edge file $edge_file\n";
    open(FILE, $edge_file);
    #
    my $strain_edge_dir = "$strain_dir/pairedEdges_i$edge_id";
    run_exe("mkdir $strain_edge_dir") if(!-d $strain_edge_dir);
    my $strain_edge_file = "$strain_edge_dir/pairedEdges_i$edge_id";
    open(OUT, ">$strain_edge_file");
    `touch $strain_edge_dir/lib.txt`;
    #
    print $OPERA_LG_CONFIG "[LIB]\n";
    print $OPERA_LG_CONFIG "map_type=opera\n";
    print $OPERA_LG_CONFIG "map_file=$strain_edge_file\n";
    print $OPERA_LG_CONFIG "lib_mean=" . ($EDGE_ID_INFO[$edge_id]) . "\n";
    print $OPERA_LG_CONFIG "lib_std=" . ($EDGE_ID_INFO[$edge_id] / 10) . "\n";
    #
    while(<FILE>){
	$str = $_;
	@line = split(/\t/, $str);
	#
	$c1 = $line[0];
	$c2 = $line[2];
	#
	$c1_strain = $contig_info{$c1}->{"STRAIN_ID"};
	$c2_strain = $contig_info{$c2}->{"STRAIN_ID"};

	next if(
	    ! (defined $c1_strain) || ! (defined $c2_strain) || 
	    $c1_strain eq "NA" || $c2_strain eq "NA" );
	
	#
	if($c1_strain eq $strain_id && $c2_strain eq $strain_id){
	    print OUT $str;
	}
	else{
	    #to rescue shared region in case of edge between contig from different strain only rescue for strain for lower coverage => higher strain ID
	    if($c1_strain * $c2_strain != 0 &&
	       ($c1_strain == $strain_id || $c2_strain == $strain_id) &&
	       ($c1_strain < $strain_id || $c2_strain < $strain_id)){
		print OUT $str;
		$shared_contig = $c1;
		$shared_contig = $c2 if($c1_strain == $strain_id);
		#$contig_info{$shared_contig}->{"SHARED_CONTIG"} = 1 if($edge_id == 0);#Only rescue local edges
	    }
	}
    }
    close(FILE);
    close(OUT);

    #Then filter rescue edges that gives rise to a local conflict
    
}

#Search for a mode for wich the negative banomial distribution conver a significant fraction of the genome studied
#If the largest contains seems to contains more than 1 genome:
#    * Try to deacrise the smoothing factor
#    * Give a warning indicating that it is not possible to perform strain analysis for that genome
#sub select_mode{
#    my ($mode_set, $assembly_length) = @_;
#    for (my $i = @{$mode_set}-1; $i >= 0; $i--){
#    }
#}

sub compute_confidance_interval{
    my ( $strain_mean_cov, $dispersion_value) = @_;
    $R->set( 'mean', $strain_mean_cov );
    $R->set( 'dispersion', $dispersion_value);

    #my @interval = ()
    
    $a = $R->run(
	q ` s <- qnbinom(c(0.02), size=dispersion, mu=mean)`,
	q `print(s)`
	);
    @tmp = split(/\s+/, $a);
    $confidance_interval[0] = $tmp[1];
    
    $a = $R->run(
	q ` s <- qnbinom(c(0.98), size=dispersion, mu=mean)`,
	q `print(s)`
	);
    @tmp = split(/\s+/, $a);
    $confidance_interval[1] = $tmp[1];
    #return \@interval;

}

sub compute_mode{
    my ($window_distrib) = @_;

    #$window_distrib = [10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10];

    #Resolve weird problem
    my @temp = @{$window_distrib}[0..@{$window_distrib}-2];#NEED TO FIX uninitiallized value push in the array
    #print STDERR $window_distrib[@{$window_distrib}-1];
    #@temp = @{$window_distrib};
    
    #foreach $f (@temp){
    #print STDERR "xx" . $f . "\n";
    #}
    
    #$R->set( 'values', $window_distrib);
    $R->set( 'values', \@temp);
    $R->set( 'span', 11);
    
    $a = $R->run(
	## adpated from EDDA
	#q `print(sessionInfo())`,	
	q `length(values)`,
	q `dens <- density(values)`,
	q `series <-dens$y`,
	q `z <- embed(series, span)`,
	q `s <- span%/%2`,
	q `ind <- apply(z, 1, which.max)`,
	q `v <- ind == (1 + s)`,
	q `result <- c(rep(FALSE, s), v)`,
	q `result <- result[1:(length(result) - s)]`,
	q `print(dens$x[result])`,
	);

    @tmp = split(/\s+/, $a);

    print STDERR " *** $a\n";

    my @res;
    @res = @tmp[1..@tmp-1]; #($tmp[1], $tmp[2]);
    @res = @tmp[2..@tmp-1] if(index($tmp[1], "[") != -1); #($tmp[2], $tmp[3]) if($tmp[1] eq "[1]");
    
    return \@res;
    
}



sub compute_probability{
    my ( $strain_mean_cov, $dispersion_value, $contig_mean_cov) = @_;
    $value = int($contig_mean_cov);
    #
    $R->set( 'mean', $strain_mean_cov );
    $R->set( 'dispersion', $dispersion_value);
    $R->set( 'k', $value);
    my $out1 = $R->run(
	q ` s <- dispersion * log(dispersion / (dispersion + mean)) + lgamma(dispersion + k) + k * log(mean / (dispersion + mean)) - lgamma(k + 1) - lgamma(dispersion)`,
	#  q`a <- $val`,
	q `print(s)`
	);

    print $out1 ."\n";
    
#return dispersion * log(dispersion / (dispersion + mean)) + //Compute once
#	lgammal(dispersion + k) + //Sum of the 2 sons
    #		k * log(mean / (dispersion + mean)) - //k * compute once
#		lgamma(k + 1) -	//Sum of the 2 sons
#		lgamma(dispersion); //Compute once globaly

}


sub estimate_mean{
    my ($mode, $dispersion) = @_;
    print STDERR " estimate_mean $mode $dispersion\n";
    return $mode * ($dispersion / ($dispersion -1));
}

sub filter_mapping{

    #my ($mapping_file);

    my $contig_mapping = "$species_dir/contig.map";
    
    #Mapping to the reference genome
    if(! -e $contig_mapping){
	my @tmp = split(/\//,$species_dir);
	my $species_name = $tmp[-1];
	my $reference = `grep $species_name $species_dir/../reference_length.dat | cut -f4`;chop $reference;

	run_exe("nucmer -p $species_dir/out $reference $species_dir/contigs.fa 2> $species_dir/nucmer.log");
	run_exe("show-coords -lrcT $species_dir/out.delta | sort -k1,1 -n > $species_dir/contig.map 2> $species_dir/coords.log");
    }
    
    #open(NUC_MAPPING, "sort -k1,1 -n $contig_mapping | ")
    open(NUC_MAPPING, "$contig_mapping")
	or die "Parsing problem during read rescue using NUCMER mapping.\n";
    
    #skip the first four lines of the cluster-vs-cluster.txt file.
    <NUC_MAPPING>;
    <NUC_MAPPING>;
    <NUC_MAPPING>;
    <NUC_MAPPING>;
    
    my ($contig, $percent_mapped, $length, $start, $end);

    my $prev_contig = -1;
    my $prev_start = -1;
    my $prev_end = -1;
    my $prev_length = -1;
    #my $prev_percent_map = -1;

    open(OUT, ">$species_dir/filtered_contig.dat");
    
    while(<NUC_MAPPING>){
	if ($_ eq ""){
	    next;
	}
	#print STDERR $_;
	chop $_;
	my @line = split(/\t/, $_);
	$contig = $line[12];
	$percent_mapped = $line[10];
	$length = $line[8];
	$map_length = $line[5];
	$start = $line[0];$end = $line[1];
	$cov = $contig_info{$contig}->{"MEAN_READ_COUNT"};
	#

	next if($map_length < 400 || $percent_mapped < 20);
	
	#
	if(($prev_start <= $start && $end <= $prev_end) ||#current alignement is incuded the previous alignement
	   ($start <= $prev_start && $prev_end <= $end)#previous alignement incuded the current alignement
	    ){
	    if($contig_info{$contig}->{"STRAIN_ID"} ne "NA" && $contig_info{$contig}->{"STRAIN_ID"} != 0 && $contig_info{$contig}->{"STRAIN_ID"} == $contig_info{$prev_contig}->{"STRAIN_ID"} ){
		$filter_contig = $contig;
		$filter_cov = $cov;
		#if($prev_length < $length){
		if($prev_cov < $cov){
		    $filter_contig = $prev_contig;
		    $filter_cov = $prev_cov;
		}
		
		print OUT $filter_contig .  "\t" . $contig_info{$filter_contig}->{"STRAIN_ID"} . "\t" . $filter_cov . "\t" . 
		    $contig . "\t" . $start . "\t" . $end . "\t" . $cov . "\t" .
		    $prev_contig . "\t" . $prev_start . "\t" . $prev_end . "\t" . $prev_cov . "\n";#<STDIN>;
		$contig_info{$filter_contig}->{"STRAIN_ID"} = "NA";
	    }
	}

	#Update the prev values
	$prev_contig = $contig;
	$prev_length = $length;
	$prev_cov = $cov;
	$prev_start = $start;
	$prev_end = $end;
	
    }

    close(OUT);
    
}

sub get_contig_sequence{
    my ($out_file, $contig_file, $contig_info) = @_;
    
    open(IN, $contig_file);
    my $print_flag = 0;

    open(OUT, ">$out_file");
    my $ID;
    while(<IN>){
	
	if($_ =~ />(.*)/){
	    #print STDERR "---> $nb_contig\n" if($nb_contig % 500000 == 0);$nb_contig++;
	    @line = split(/\s+/, $1);
	    $ID = $line[0];
	    	    
	    if(exists $contig_info->{$ID}->{"STRAIN_ID"} && ($contig_info->{$ID}->{"STRAIN_ID"} eq "NA" || $contig_info->{$ID}->{"STRAIN_ID"} == 0 )){
		$print_flag = 1;
		print OUT ">" . $ID . "\n";
	    }
	    else{
		#print STDERR "**************** IN\n";
		$print_flag = 0;
		
	    }
	}
	else{
	    print OUT $_ if($print_flag);
	}
    }
    close(IN);
    close(OUT);
    
}


    
sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR "\n".$exe."\n";#<STDIN>;
    print STDERR `$exe` if($run);
}

#add plot for each species representing the histogram of the window size, the density of the estimated neg-binom and the interval used for the window selection
#perhaps can also add the assembly size
#gg=read.table("ANALYSIS/ASSEMBLY/TLL11/intermediate_files/strain_analysis/Faecalibacterium_prausnitzii/window_distibution.dat");hist(gg[,1], breaks=1000)
#m=10;N <- rnbinom(n = 1000, size = 5, mu = m);g=density(N);plot(g);abline(v=m)
