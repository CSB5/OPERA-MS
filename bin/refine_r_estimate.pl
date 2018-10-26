#!/usr/bin/perl
use warnings "all";
use Statistics::Basic qw(:all);

my ($coverage_contig_file, $sigma_dir, $opera_ms_dir) = @ARGV;

my $NB_STEP = 10;
my $FRACTION_THRESHOLD = 0.01;
my $CLUSTER_FRACTION_THRESHOLD = 0.01;
my $LENGTH_THRESHOLD = 6000000;

my $cluster_file = "$sigma_dir/clusters";

my $tmp = `cat $sigma_dir/r_estimate_value.dat`;chop $tmp;
my @r_estimate_value = split(/\t/, $tmp);
my $step = ($r_estimate_value[1] - $r_estimate_value[0]) / $NB_STEP;

print STDERR " *** R refinement " . "\t" .  $r_estimate_value[0] . "\t" . $r_estimate_value[1] . "\t" . $step . "\n";;

my $tmp_config_file = "$sigma_dir/tmp";
my $sigma_config_file = "$sigma_dir/sigma.config";
open(OUT, ">$sigma_dir/r_value_evaluation_summary.dat");
for (my $i = 0; $i <= $NB_STEP; $i++){
    $r_value = $r_estimate_value[0] + $i * $step;

    if($i != 0 ){
	
	run_exe("grep -v R_VALUE $sigma_config_file > $tmp_config_file; echo \"R_VALUE=$r_value\" >> $tmp_config_file");
	run_exe("mv $tmp_config_file $sigma_config_file");
	
	run_exe("$opera_ms_dir/bin/sigma $sigma_config_file");
	
    }
    run_exe("$opera_ms_dir/bin/cluster_evaluation.pl $coverage_contig_file $cluster_file 2> $sigma_dir/r_estimate_$r_value.dat");

    #get the fraction of the assembly with outlier contigs
    $tmp = `tail -n1 $sigma_dir/r_estimate_$r_value.dat`; chop $tmp;
    @cluster_eval_measure = split(/\t/, $tmp);
    $outlier_fraction = $cluster_eval_measure[1]/$cluster_eval_measure[0];
    $max_cluster_length = $cluster_eval_measure[2];
    #
    $max_cluster_outlier_fraction = $cluster_eval_measure[3];
    $max_cluster_outlier_fraction_length = $cluster_eval_measure[4];
    print OUT $r_value . "\t" . 
	$outlier_fraction . "\t" . $FRACTION_THRESHOLD  . "\t" . 
	$max_cluster_length . "\t" . $LENGTH_THRESHOLD  . "\t" .
	$max_cluster_outlier_fraction . "\t" . $max_cluster_outlier_fraction_length . "\t" .  $CLUSTER_FRACTION_THRESHOLD . "\n";#<STDIN>;

    #r value should be higher than 1 as value below does not fit our model
    last if($r_value > 1 && $outlier_fraction < $FRACTION_THRESHOLD && $max_cluster_length < $LENGTH_THRESHOLD && $max_cluster_outlier_fraction < $CLUSTER_FRACTION_THRESHOLD);
    
    #run_exe("evaluate_cluster.pl
}
close(OUT);


sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR "$return\n" if($run);
    return $return;
}
