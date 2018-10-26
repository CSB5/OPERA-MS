#!/usr/bin/perl
use warnings;
use Switch;
use Getopt::Long;

#my ($seq_file,  $graphmap_dir, $water_dir, $racon_dir);
my ($seq_file,  $graphmap_dir, $water_dir, $racon_dir, $edge_file_info);

GetOptions(
    "seq-file=s" => \$seq_file,
    "edge-file=s" => \$edge_file_info,
    "graphmap=s" => \$graphmap_dir,
    "water=s" => \$water_dir,
    "racon=s" => \$racon_dir,
    #"edge_file_info=s" => \$edge_file_info,
    )
 or die("Error in command line arguments.");

if (defined $graphmap_dir and $graphmap_dir ne ""){
    $graphmap_dir .= "/";
}
else{
    $graphmap_dir = "";
}
if (defined $water_dir and $water_dir ne ""){
    $water_dir .= "/";
}
else{
    $water_dir = "";
}

if(defined $racon_dir and $racon_dir ne ""){
    $racon_dir .= "/";
}
else{
    $racon_dir = "";
}

#open(F,">>debug.txt");

#To get the ref and seq
my $ref = $seq_file."_r.fa";
#my $dbg = $seq_file."_dbg";
my $queries = $seq_file."_q.fa";
my $mapping = "$seq_file.m5";
my $all_seq = $seq_file.".fa";
#my $edge_file_info = '/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/ANALYSIS/ASSEMBLY/X_OPERA_MS_REF/FULL_SAMPLE/MEGAHIT/NANOPORE/RESULTS/intermediate_files/long-read-mapping/edge_file_info.dat';
my @temp_list = split(/\//, "$seq_file");
#print F "".$#temp_list."\n";
#exit(0);
#open(DBG, ">".$dbg);
my $contigs = $temp_list[-1];
@temp_list = split(/-vs-/,"$contigs");
my $contig_1 = $temp_list[0];
my $contig_2 = $temp_list[1];
#print DBG $contig_1."\t".$contig_2."\n";
#print STDERR "awk '{if(\$1 == $contig_1 && \$3 == $contig_2) print \$0}' $edge_file_info | cut -f 5 | cut -d\";\" -f 1";
my $temp = `awk '{if(\$1 == \"$contig_1\" && \$3 == \"$contig_2\") print \$0}' $edge_file_info | cut -f 5 | cut -d\";\" -f 1`;
#print DBG $temp;
#print STDERR " *** |$temp|\n";<STDIN>;
$temp = `awk '{if(\$1 == \"$contig_2\" && \$3 == \"$contig_1\") print \$0}' $edge_file_info | cut -f 5 | cut -d\";\" -f 1` if($temp eq "");
#print DBG $temp;
#print STDERR " *** |$temp|\n";<STDIN>;
my $read_name = "@".$temp;chop $read_name;
print STDERR " *** read_name |$read_name| => $all_seq \n";

#print STDERR $read_name . "\n";<STDIN>;

#open(F,">>debug.txt");
#print F $contig_1;
#print F $contig_2;
    
#if($contig_1=="k99_701"||$contig_2=="k99_701"){
#    print F "$read_name CATCHSITE";
#}

open(REF,">".$ref);
open(Q,">".$queries);
open(FILE,$all_seq);
#open(DBG, ">".$dbg);
while(<FILE>) {
    #print STDERR " *** $_ |$read_name|\n";<STDIN>;
    if(index($_, $read_name) != -1) {
        #print DBG $read_name."\tb\t".$_;
        print STDERR "Read $read_name\n";
        print REF $_;
        $tmp = <FILE>; print REF $tmp;
        $tmp = <FILE>; print REF $tmp;
        $tmp = <FILE>; print REF $tmp;
    }
    else{
    	print Q $_;
    	$tmp = <FILE>; print Q $tmp;
    	$tmp = <FILE>; print Q $tmp;
    	$tmp = <FILE>; print Q $tmp;
    }
}
close(REF);
close(Q);
close(FILE);
#close(DBG);
#run_exe("head -n4 $seq_file.fa > $ref;tail -n+5 $seq_file.fa > $queries");
#run_exe("tail -n+5 $seq_file.fa > $queries");
#
#Selection based on read_name provided as input 
#

#For the problem of the missing 1st mapping
#run_exe("head -n4 $seq_file.fa | tail -n2 >> $queries");

#Run the mapping
#OLD GRAPHMAP

$mapping = "$seq_file.sam";
$mapping = "$seq_file.paf";

#Test which version of graph map is used
#run_exe("which graphmap > $seq_file.gversion");

#run_exe("${graphmap_dir}graphmap align -t 1 -r $ref -d $queries -o $mapping");
#run_exe("${graphmap_dir}minimap2 -Sw5 -L100 -m0 $ref  $queries > $mapping");

#NEW GRAPHMAP
#run_exe("$graphmap_dir/graphmap align -t 1 -r $ref -d $queries | python $bin_dir/sam2blasr.py -i - -r $ref > $mapping");

#remove the graphmap indexes
#run_exe("rm $ref*gmidx*");
#run_exe("rm $ref*gm* $ref $queries");

#Get the consensus using racon
#my $cons;

#$cons = $seq_file."_racon.fa";
#my $cmd = "${racon_dir}racon --sam $queries $mapping $ref $cons";
#my $cmd = "${racon_dir}racon  $queries $mapping $ref $cons"; 
#run_exe($cmd);

#mapping of the contig to the consensus to get the exact position of the bondaries
#$queries = $seq_file."_contig.fa";

#Convert fastq to fasta
#run_exe("tail -n8 $seq_file.fa | sed -n -e 1p -e 2p -e 5p -e 6p | sed 's/@/>/' > $queries");
#For the problem of the missing 1st mapping
#run_exe("tail -n4 $seq_file.fa | head -n2 >> $queries");
#$mapping = "$seq_file\_contig_cons.m5";
#Map using Smith-Waterman alignement method
#run_exe("${water_dir}water -gapextend 0.5 -gapopen 10 -asequence $cons -bsequence $queries -outfile $mapping");
#run_exe("graphmap -x illumina -t 1 -r $cons -d $queries | python ~lich/projects_backup/INCSeq/utils/sam2blasr.py -i - -r $cons > $mapping");   

sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";;
    print STDERR `$exe` if($run);
}




