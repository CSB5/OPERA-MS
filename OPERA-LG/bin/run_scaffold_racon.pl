#!/usr/bin/perl

#use strict;
use warnings;

my ($scaffold_dir, $racon_dir, $minimap2_dir) = @ARGV;

my $start_time_total = time;

my $start_time = time;
my $end_time;

my $re_run = 1;

my @tmp = split(/\//, $scaffold_dir);


my $assembly_dir =  "$scaffold_dir";
run_exe("mkdir $assembly_dir");
    
#Run racon
#Get the initial assembly
my $initial_assembly = "$assembly_dir/pre_consensus.fa";

#Consensus construction
my $mapping_file = "$initial_assembly.filtered.paf";
my $read_file = "$assembly_dir/POOL.fastq";
my $racon_consensus = "$assembly_dir/racon.fa";

if(-e $read_file){
    #run_exe("mkdir $assembly_dir/../CORR/;touch $assembly_dir/../CORR/single_read.map");
    get_read_mapping($read_file, $initial_assembly, $mapping_file);
    #$end_time = time;
    #print STDERR "*** Get read sequence Elapsed time: " . ($end_time - $start_time) . "s\n";
    #$start_time = time;
    
    run_racon($initial_assembly, $read_file, $mapping_file, $racon_consensus, 500, 5, -4, -8, -6, "CONSENSUS");
    #
    $end_time = time;
    print STDERR "*** Racon consensus assembly construction Elapsed time: " . ($end_time - $start_time) . "s\n";
    $start_time = time;
}
else{
    #There no long reads mapping to the scaffold, the scaffolds must be linked by short reads
    run_exe("cp $initial_assembly $racon_consensus");
}


sub run_racon{
    my ($assembly_file, $read_file, $mapping, $out_file, $window_size, $match, $mismatch, $gap_open, $gap_extend, $RACON_TYPE) = @_;
    #my ($window_size, $match, $mismatch, $gap_open, $gap_extend);

    my $pileup_option = "";
    $pileup_option = "--pileup" if(index($out_file, "pileup") != -1);
    
    open(OUT, ">$out_file.sh");
    #print OUT "source activate nanopore\n";
    #print OUT "racon --winlen 2000 --ovl-margin 0.5 --reads $fastq_pool --alnpath $mapping --raw $assembly_file --out $racon_out.fa\n";
    #print OUT "$racon_path $fastq_pool $mapping $assembly_file > $out_file\n";
    print OUT " $racon_dir/racon --winlen $window_size --match $match --mismatch $mismatch --gapopen $gap_open --gapext $gap_extend --reads $read_file --alnpath $mapping --raw $assembly_file --out $out_file\n";
    #print OUT "source deactivate\n";
    close(OUT);
    run_exe("sh $out_file.sh");
}


#Should be done once to avoid reading the read file multiple times
sub get_read_mapping{
    my ($long_read, $assembly, $mapping_file) = @_;
    
    run_exe("$minimap2_dir/minimap2 --cs=short -w5 -m0 $assembly $long_read > $mapping_file");
    $end_time = time;
    print STDERR "*** Map read in gap Elapsed time: " . ($end_time - $start_time) . "s\n";
    $start_time = time;
    
}



sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return . "\n" if($run);
    return $return;
}

    

