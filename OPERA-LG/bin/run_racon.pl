#!/usr/bin/perl
use warnings;
use Switch;
use Getopt::Long;

my ($seq_file, $bin_dir, $graphmap_dir, $water_dir, $racon_dir);

GetOptions(
    "seq-file=s" => \$seq_file,
    "bin-dir=s" => \$bin_dir,
    "graphmap=s" => \$graphmap_dir,
    "water=s" => \$water_dir,
    "racon=s" => \$racon_dir,
    )
 or die("Error in command line arguments.");

if (defined $graphmap_dir){
    $graphmap_dir .= "/";
}
else{
    $graphmap_dir = "";
}
if (defined $water_dir){
    $water_dir .= "/";
}
else{
    $water_dir = "";
}

#To get the ref and seq
my $ref = $seq_file."_r.fa";
my $queries = $seq_file."_q.fa";
my $mapping = "$seq_file.m5";


run_exe("head -n4 $seq_file.fa > $ref;tail -n+5 $seq_file.fa > $queries");


#For the problem of the missing 1st mapping
#run_exe("head -n4 $seq_file.fa | tail -n2 >> $queries");

#Run the mapping
#OLD GRAPHMAP

$mapping = "$seq_file.sam";
run_exe("${graphmap_dir}graphmap align -t 1 -r $ref -d $queries -o $mapping");

#NEW GRAPHMAP
#run_exe("$graphmap_dir/graphmap align -t 1 -r $ref -d $queries | python $bin_dir/sam2blasr.py -i - -r $ref > $mapping");

#remove the graphmap indexes
run_exe("rm $ref*gmidx*");
#run_exe("rm $ref*gm* $ref $queries");

#Get the consensu using pbdagcon
my $cons;

$cons = $seq_file."_racon.fa";
my $cmd = "$racon_dir/racon --sam $queries $mapping $ref $cons"; 
run_exe($cmd);

#mapping of the contig to the consensus to get the exact position of the bondaries
$queries = $seq_file."_contig.fa";

run_exe("tail -n8 $seq_file.fa | sed -n -e 1p -e 2p -e 5p -e 6p | sed 's/@/>/' > $queries");
#For the problem of the missing 1st mapping
#run_exe("tail -n4 $seq_file.fa | head -n2 >> $queries");
$mapping = "$seq_file\_contig_cons.m5";
#Map using Smith-Waterman alignement method
run_exe("${water_dir}water -gapextend 0.5 -gapopen 10 -asequence $cons -bsequence $queries -outfile $mapping");
#run_exe("graphmap -x illumina -t 1 -r $cons -d $queries | python ~lich/projects_backup/INCSeq/utils/sam2blasr.py -i - -r $cons > $mapping");	

sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";;
    print STDERR `$exe` if($run);
}
