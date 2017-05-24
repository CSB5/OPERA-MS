#!/usr/bin/perl
use warnings;
use Switch;
use Getopt::Long;

my ($seq_file, $graphmap_dir, $water_dir, $pbdagcon_dir, $bin_dir) = @ARGV;

#To get the ref and seq
my $ref = $seq_file."_r.fa";
my $queries = $seq_file."_q.fa";
my $mapping = "$seq_file.m5";
run_exe("head -n2 $seq_file.fa > $ref;tail -n+3 $seq_file.fa > $queries");
#For the problem of the missing 1st mapping
#run_exe("head -n4 $seq_file.fa | tail -n2 >> $queries");

#Run the mapping
run_exe("$graphmap_dir/graphmap -t 1 -r $ref -d $queries | python $bin_dir/sam2blasr.py -i - -r $ref > $mapping");
#remove the graphmap indexes
run_exe("rm $ref*gmidx*");
#run_exe("rm $ref*gm* $ref $queries");

#Get the consensu using pbdagcon
my $nb_base_to_trim = 0;
#my $kill_res = "";
my $nb_line_cons = 0;
my $inc = 1;
my $cons= $seq_file."_pbdagcon.fa";
my $pid;
#while($kill_res eq ""){
while($nb_line_cons == 0){
    $cons = $seq_file."_pbdagcon.fa";
    $cmd = "$pbdagcon_dir/pbdagcon -t $nb_base_to_trim -c 1 -m 1 $mapping > $cons";
    #run the command
    print STDERR " *** Running $cmd\n";
    $pid = open(my $ph, "-|", $cmd) or die $!;
    #print STDERR " *** Cmd $pid\n";
    #sleep 5s
    run_exe("sleep 5");
    #try to kill the process
    $nb_line_cons = `wc -l $cons | cut -d \" \" -f1`;chop $nb_line_cons;
    #print $nb_line
    #print STDERR " **** $kill_cmd\n";
    #$kill_res = `$kill_cmd`;chop $kill_res;
    #print STDERR "|$kill_res|\n";
    run_exe("kill -9 $pid");
    $nb_base_to_trim += $inc;
    $inc = 10 if($nb_base_to_trim == 10);
    #last if($inc == 100);
}

#mapping of the contig to the consensus to get the exact position of the bondaries
$queries = $seq_file."_contig.fa";
run_exe("tail -n4 $seq_file.fa > $queries");
#For the problem of the missing 1st mapping
#run_exe("tail -n4 $seq_file.fa | head -n2 >> $queries");
$mapping = "$seq_file\_contig_cons.m5";
#Map using Smith-Waterman alignement method
run_exe("$water_dir/water -gapextend 0.5 -gapopen 10 -asequence $cons -bsequence $queries -outfile $mapping");
#run_exe("graphmap -x illumina -t 1 -r $cons -d $queries | python ~lich/projects_backup/INCSeq/utils/sam2blasr.py -i - -r $cons > $mapping");	

sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";;
    print STDERR `$exe` if($run);
}
