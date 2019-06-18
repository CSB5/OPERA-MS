#!/usr/bin/perl
use warnings "all";
use Statistics::Basic qw(:all);

my ($coverage_contig_file, $sigma_dir, $num_processor, $opera_ms_dir) = @ARGV;

my $NB_STEP = 10;
my $FRACTION_THRESHOLD = 0.01;
my $CLUSTER_FRACTION_THRESHOLD = 0.01;
my $LENGTH_THRESHOLD = 6000000;

my $cluster_file = "$sigma_dir/clusters";

my $tmp = `cat $sigma_dir/r_estimate_value.dat`;chop $tmp;
my @r_estimate_value = split(/\t/, $tmp);
my $r_value = $r_estimate_value[0];
my $step = ($r_estimate_value[1] - $r_estimate_value[0]) / $NB_STEP;

print STDERR " *** R refinement " . "\t" .  $r_estimate_value[0] . "\t" . $r_estimate_value[1] . "\t" . $step . "\n";;
open(OUT_SUM, ">$sigma_dir/r_value_evaluation_summary.dat");
#Cher if the first run provide good results
if(! evaluate_cluster($coverage_contig_file, $cluster_file, "$sigma_dir/r_estimate_$r_value.dat")){
    #Run sigma with other R parameters
    my $sigma_config_file = "$sigma_dir/sigma.config";
    my $sigma_cmd = "$sigma_dir/cmd.sh";
    open(OUT, ">$sigma_cmd");
    for (my $i = 1; $i <= $NB_STEP; $i++){
	$r_value = $r_estimate_value[0] + $i * $step;
	
	$new_dir = "$sigma_dir/sigma_$r_value";
	$new_config_file = "$sigma_dir/sigma_$r_value/sigma.config";
	run_exe("mkdir  $new_dir") if(! -d $new_dir);
	run_exe("grep -v R_VALUE $sigma_config_file | grep -v output_dir  > $new_config_file;  echo \"output_dir=$new_dir\" >> $new_config_file; echo \"R_VALUE=$r_value\" >> $new_config_file");
	#run_exe("mv $tmp_config_file $sigma_config_file");
	print OUT "$opera_ms_dir/bin/sigma $new_config_file\n";
    }
    close(OUT);
    #exit(0);
    run_exe("cat $sigma_cmd | xargs -L 1 -P $num_processor -I COMMAND sh -c \"COMMAND\" 2> $sigma_cmd.log");


    my $r_index = 1;
    for ($r_index = 1; $r_index <= $NB_STEP; $r_index++){
	$r_value = $r_estimate_value[0] + $r_index * $step;
	$cluster_file = "$sigma_dir/sigma_$r_value/clusters";
	last if(evaluate_cluster($coverage_contig_file, $cluster_file, "$sigma_dir/r_estimate_$r_value.dat"));
    }
    #Copy the cluster (and edges) used for the next steps
    run_exe("rm $sigma_dir/clusters; ln -s $sigma_dir/sigma_$r_value/clusters $sigma_dir/clusters");
}
close(OUT_SUM);


sub evaluate_cluster{
    my ($coverage_contig_file, $cluster_file, $out_file) = @_;

    $res = 0;
    
    run_exe("${opera_ms_dir}utils/perl $opera_ms_dir/bin/cluster_evaluation.pl $coverage_contig_file $cluster_file 2> $out_file");
    
    #get the fraction of the assembly with outlier contigs
    $tmp = `tail -n1 $out_file`; chop $tmp;
    @cluster_eval_measure = split(/\t/, $tmp);
    $outlier_fraction = $cluster_eval_measure[1]/$cluster_eval_measure[0];
    $max_cluster_length = $cluster_eval_measure[2];
    #
    $max_cluster_outlier_fraction = $cluster_eval_measure[3];
    $max_cluster_outlier_fraction_length = -1;
    $max_cluster_outlier_fraction_length = $cluster_eval_measure[4];
    print OUT_SUM $r_value . "\t" . 
	$outlier_fraction . "\t" . $FRACTION_THRESHOLD  . "\t" . 
	$max_cluster_length . "\t" . $LENGTH_THRESHOLD  . "\t" .
	$max_cluster_outlier_fraction . "\t" . $max_cluster_outlier_fraction_length . "\t" .  $CLUSTER_FRACTION_THRESHOLD . "\n";#<STDIN>;

    #r value should be higher than 1 as value below does not fit our model
    $res = 1 if($r_value > 1 && $outlier_fraction < $FRACTION_THRESHOLD && $max_cluster_length < $LENGTH_THRESHOLD && $max_cluster_outlier_fraction < $CLUSTER_FRACTION_THRESHOLD);    
    
    return $res;
}



sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR "$return\n" if($run);
    return $return;
}
