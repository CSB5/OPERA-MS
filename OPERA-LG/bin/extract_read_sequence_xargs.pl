#!/usr/bin/perl
use warnings;
use Switch;
use Getopt::Long;
#`source activate nanopore`;



my $edge_file_info;
my $contig_file;
my $scaffold_file;
my $read_file;
my $script_dir;
my $vsearch_dir;
my $water_dig;
my $graphmap_dir;
my $racon_dir;
my $ana_dir;
my $outfile;
my $stage;
my $nb_process;
#my $bin_dir;
my $is_fastq = 0;

GetOptions(
    "edge-file=s"    => \$edge_file_info,
    "contig-file=s"    => \$contig_file,
    "scaffold-file=s" => \$scaffold_file,
    "output-directory=s" => \$ana_dir,
    "num-of-processors=i" => \$nb_process,
    "script=s"     => \$script_dir,
    "read-file=s"     => \$read_file,
    "vsearch=s"       => \$vsearch_dir,
    "water=s"         => \$water_dir,
    "racon=s"         => \$racon_dir,
    "graphmap=s"      => \$graphmap_dir,
    "outfile=s"       => \$outfile,
    "stage=s"         => \$stage,
 #   "bin=s"           => \$bin_dir
    )

 or die("Error in command line arguments.\n");


if(!(-e $contig_file)){
    die "CONTIG FILE NOT FOUND \n";
}	

if(!(-e $scaffold_file)){
    die "SCAFFOLD FILE NOT FOUND \n";
}

if(!(-e $read_file)){
    die "READ FILE NOT FOUND \n";
}

if(!(defined $racon_dir)){
    $racon_dir = "";
}

else{
    $racon_dir = "--racon " . $racon_dir;
}

if(!(defined $water_dir)){
    $water_dir = "";
}

else{
    $water_dir = "--water " . $water_dir;
}

if(!(defined $graphmap_dir)){
    $graphmap_dir = "";
}

else{
    $graphmap_dir = "--graphmap " . $graphmap_dir;
}

if(!(defined $vsearch_dir)){
    $vsearch_dir = "";
}


$stage = "ALL" if(! defined $stage);

#Check if the read file is in fastq format or fasta.
open(READ_FILE, $read_file);
my $first_line = <READ_FILE>;
if($first_line =~ /^\@/){
    $is_fastq = 1;
}

else{
    $is_fastq = 0;
}
close(READ_FILE);

my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -V";

#length of the mapping that include the contig sequence
##?##
my $sequence_overhang = 1000;
#my $additional_contig_sequence = 100;
my $additional_contig_sequence = 0;

#For only getting the consensus
if($stage eq "CONS"){
    print STDERR " *** ConsEtruct the gap consensus sequences\n";
    run_exe("rm -rf $ana_dir/LOG;mkdir $ana_dir/LOG");
    #opendir(DIR, $ana_dir);
    #my @all_cons = readdir(DIR);
    #close(DIR);
    
    #foreach $file (@all_cons){
    #@next if(index($file, ".fa") == -1 || index($file, "con") != -1 || index($file, "_q.fa") != -1 || index($file, "_r.fa") != -1);
    #$edge_ID = $file;chop $edge_ID;chop $edge_ID;chop $edge_ID;
    #run_exe("/home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/run_pbdagcon.pl $ana_dir/$name");
    run_exe("cat GAP_FILLING/cmd.txt | xargs -L 1 -P $nb_process -I COMMAND sh -c \"COMMAND\"");
    #run_exe("$qsub -l mem_free=10G,h_rt=00:30:0 -pe OpenMP 1 -N vsearch -e $ana_dir/LOG/$edge_ID.err -o $ana_dir/LOG/$edge_ID.run -b y vsearch --cluster_fast $ana_dir/$edge_ID.fa -id 0 --msaout $ana_dir/$edge_ID\_cons.fa -uc $ana_dir/$edge_ID\_cons.dat");#<STDIN>;
    #run_exe("$qsub  -l mem_free=10G,h_rt=00:30:0 -pe OpenMP 2 -N pbdagcon -e $ana_dir/LOG/$edge_ID.err -o $ana_dir/LOG/$edge_ID.run -b y /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/run_pbdagcon.pl $ana_dir/$edge_ID");#<STDIN>;
	#run_exe("sleep 0.5");
#}
    exit(0);
}

if($stage eq "FILL"){
    print STDERR " *** Fill the scaffold\n";#<STDIN>;
    #run_exe("rm -r $ana_dir/LOG;mkdir $ana_dir/LOG");
    
    #read the scaffold file to get the edges that should be extarcted
    read_scaffold_file($scaffold_file);
    
    #write the filled sequences
    write_filled_gap_scaffold();

    exit(0);

}



print STDERR " *** Starting the sequence extraction\n";#<STDIN>;
run_exe("rm -r $ana_dir");
run_exe("rm -rf $ana_dir/LOG");
run_exe("mkdir -p $ana_dir/LOG");
my $run_racon_log = "$ana_dir/LOG/run_racon_log.txt";

#Structure that indicate sequence required to fill the gap between the 2 contigs involvolved in the edge
my %edge_info = ();
#indicate which reads belong to which edges
my %read_to_edge = ();
#indicate which contig belong to which edges
my %contig_to_edge = ();
#holds the reference reads
my %reference_reads = ();

#Read the scaffold file to get the set of adjacent edges used to construct the scaffold
my $nb_selected_edges = 0;
my $cmp_multi_edge = 1;my $str_multi = "";

#open(FILE, "head -n5 $scaffold_file |");
read_scaffold_file($scaffold_file);

print STDERR " *** Number of edges selected for gapfilling $nb_selected_edges\n";#<STDIN>;

#Get the read and the coordinates of the gaps
open(FILE, $edge_file_info);
my ($reverse_edge, $reverse_edge_r);
my %contig_in_gap_to_extract = ();#list of contig for which the full sequence need to be stored
my %contig_in_gap_support = ();#To compute the support of the contig in the gap, the contig should be at least supported by 2 reads
my %contig_lengths = ();	# Temporarily holds the contig lengths so as to compute the coordinates for cutting the reads

my @tmp_path = split(/long-read-mapping/,$edge_file_info);
my $contig_len_file = $tmp_path[0]."/opera_long_read/contigs";
open(LEN_FILE, $contig_len_file);

<LEN_FILE>;

while(<LEN_FILE>){
	@fields = split(/\t/, $_);
	$contig_lengths{$fields[0]} = $fields[1];
}

close(LEN_FILE);

# Only for debugging to verify correctness of input
# Verified to match file. Discard this as a concern
#foreach $c (keys %contig_lengths){
#	print STDERR "$c\t$contig_lengths{$c}";
#}
#exit(0);


open(OUT_EXTRA, ">$ana_dir/contig_extention.log");
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $contig_1 = $line[0];$contig_1_ori = $line[1];
    $contig_2 = $line[2];$contig_2_ori = $line[3];

    #Get the edge ID and skip the edge that were not part of scaffold
    #Get the edge_ID and corresponding representative orientation
    $tmp_array = get_edge_ID($contig_1, $contig_1_ori, $contig_2, $contig_2_ori);
    $edge_ID = $tmp_array->[0];
    $edge_ori = $tmp_array->[1];
    
    if(exists $edge_info{$edge_ID}){
#	print STDERR " *** edge ID ".$edge_ID." ".$edge_ori."\n";#<STDIN>;
	$reverse_edge = 0;
	$reverse_edge = 1 if($edge_ori ne $edge_info{$edge_ID}->{"EDGE_ORI"});
    }
    else{next;}

    #Get the read name and coordinates
    @read_name = split(/\;/, $line[4]);

    @tmp = split(/\:/, $line[5]);@coord_contig_1 = split(/\;/, $tmp[1]);
    @tmp = split(/\:/, $line[6]);@coord_contig_2 = split(/\;/, $tmp[1]);
    #
    @tmp = split(/\:/, $line[7]);@coord_contig_1_on_read = split(/\;/, $tmp[1]);
    @tmp = split(/\:/, $line[8]);@coord_contig_2_on_read = split(/\;/, $tmp[1]);
    #
    @tmp = split(/\:/, $line[9]);@contig_in_gap = split(/\;/, $tmp[1]);
    %contig_in_gap_support = ();


    ##?##
    @coord_c1_on_r1 = split(/\_/, $coord_contig_1_on_read[0]);
    @coord_c2_on_r1 = split(/\_/, $coord_contig_2_on_read[0]);
    @coord_c1_for_r1 = split(/\_/, $coord_contig_1[0]);
    @coord_c2_for_r1 = split(/\_/, $coord_contig_2[0]);


    #Calculate the maximum possible overhang based on read 1
    $max_seq_overhang_c1 = abs($coord_c1_on_r1[2] - $coord_c1_on_r1[1]);$max_seq_overhang_c1 = $sequence_overhang if($max_seq_overhang_c1 > $sequence_overhang);
    $max_seq_overhang_c2 = abs($coord_c2_on_r1[2] - $coord_c2_on_r1[1]);$max_seq_overhang_c2 = $sequence_overhang if($max_seq_overhang_c2 > $sequence_overhang);
    $max_seq_overhang_c1 = 0;
    $max_seq_overhang_c2 = 0;
    $edge_info{$edge_ID}->{"CONTIG_OVERHANG_INFO"}->{$contig_1} = $max_seq_overhang_c1;
    $edge_info{$edge_ID}->{"CONTIG_OVERHANG_INFO"}->{$contig_2} = $max_seq_overhang_c2;

    # Stores the reference read name (1st read)
    #$edge_info{$edge_ID}->{"REF_READ"} = $read_name[0];
    $reference_reads{$edge_ID} = $read_name[0];

    # Stores the coordinates of the contigs which map to the reference genome
    # Not required anymore
    #$edge_info{$edge_ID}->{"REF_COORDS"}->{$contig_1} = @coord_c1_on_r1;
    #$edge_info{$edge_ID}->{"REF_COORDS"}->{$contig_2} = @coord_c2_on_r1;
    
    for(my $i = 0; $i < @read_name; $i++){
	$r = $read_name[$i];
	$reverse_edge_r = $reverse_edge;

	#####################################################
	#Initailize the sequence overhang based on the first read
	#Store the read name and the overhang in  $edge_info{$edge_ID} {"REF_READ", "CONTIG1_OVERHANG", aa,  "CONTIG2_OVERHANG", bb}
	#
	#
	#
	#####################################################
	
	#Get the contig orientation to known if the extracted sequences have to be reversed
	@coord_c1 = split(/\_/, $coord_contig_1[$i]);
	@coord_c2 = split(/\_/, $coord_contig_2[$i]);
	#Check the orientation of the contig on the reads
	$ori_c1 = ($coord_c1[0] == 0 ? "+" : "-"); $ori_c2 = ($coord_c2[0] == 0 ? "+" : "-");
	$reverse_edge_r = ($reverse_edge_r + 1) % 2 if($ori_c1 ne $contig_1_ori && $ori_c2 ne $contig_2_ori);
	
	#For the contig that belong to the gap
	@contig_in_gap_read =  split(/\,/, $contig_in_gap[$i]);
	foreach $c (@contig_in_gap_read){
	    if ( $c =~ /(\+|\-)(.+)/ ) {
		$c_o = $1;
		$c_n = $2;
		if($reverse_edge_r){
		    $c_o = reverse_ori($c_o);
		}
	    $contig_in_gap_support{$c_n.$c_o} = 0 if(! exists $contig_in_gap_support{$c_n.$c_o});
	    $contig_in_gap_support{$c_n.$c_o}++;
	    }
	}

	#This read is used by this edge
	if(! exists $read_to_edge{$r}){$read_to_edge{$r}={};}
	#Check for the order of the contigs and compute the coordinate of the read sequences that will be extracted
	@coord_c1_on_read = split(/\_/, $coord_contig_1_on_read[$i]);
	@coord_c2_on_read = split(/\_/, $coord_contig_2_on_read[$i]);


	#PROBLEM WITH READ THAT SUPPORT 2 DISTANCES OF AN EDGE
	$str_multi = "";
	if(exists $read_to_edge{$r}->{$edge_ID}){
	    #This read support twice the distance between the 2 contigs
	    $str_multi = "[MULTI_ID$cmp_multi_edge]";
	    $cmp_multi_edge++;
	}

	#my $seq_overhang_c1 = 0;
	#my $seq_overhang_c2 = 0;

	#Update the coordinates to extend the contig coordinates to a full mapping in the read
	if($coord_c1_on_read[2] < $coord_c2_on_read[2]){#Contig 1 is first on the read
		
	    # Need ends of contigs
	    my $end_contig_1_extra = $contig_lengths{$contig_1} - $coord_c1_for_r1[2];
	    #if($coord_c1_for_r1[0] == 1) {
	    #	$end_contig_1_extra = $contig_lengths{$contig_1} - ($contig_lengths{$contig_1} - $coord_c1_for_r1[1]);
	    #}
	    if($end_contig_1_extra < 0){
		print STDERR "Weird Stuff : Negative lengths $contig_1\t$contig_lengths{$contig_1}\t$end_contig_1_extra\t$coord_c1_for_r1[2]\n";
	    }
	    #$end_contig_1_extra = 0;
	    my $start_contig_2_extra = $coord_c2_for_r1[1];
	    #if($coord_c2_for_r1[0] == 1) {
	    #	$start_contig_2_extra = $contig_lengths{$contig_2} - $coord_c2_for_r1[2];
	    #}
	    if($start_contig_2_extra < 0){
		print STDERR "Weird Stuff : Negative lengths $contig_2\t$contig_lengths{$contig_2}\t$start_contig_2_extra\t$coord_c2_for_r1[1]\n";
	    }
	    #$start_contig_2_extra = 0;
	    
	    
	    print OUT_EXTRA $r . "\t" . $edge_ID . "\t" . $end_contig_1_extra . "\t" . $start_contig_2_extra . "\n";
		
	    $read_to_edge{$r}->{$edge_ID.$str_multi} = [$coord_c1_on_read[2] + $end_contig_1_extra, $coord_c2_on_read[1] - $start_contig_2_extra, $reverse_edge_r];
	    
	    #$edge_info{$edge_ID}->{"READ_OVERHANG_INFO"}->{$readname[$i]} = {"contig_1_overhang", $seq_overhang_c1, "contig_2_overhang", $seq_overhang_c2 }

	}
	else{
	    # Need ends of contigs
	    my $end_contig_2_extra = $contig_lengths{$contig_2} - $coord_c2_for_r1[2];
	    #if($coord_c2_for_r1[0] == 1) {
	    #	$end_contig_2_extra = $contig_lengths{$contig_2} - ($contig_lengths{$contig_2} - $coord_c2_for_r1[1]);
	    #}
	    if($end_contig_2_extra < 0){
		print STDERR "Weird Stuff : Negative lengths $contig_2\t$contig_lengths{$contig_2}\t$end_contig_2_extra\t$coord_c2_for_r1[2]\n";
	    }
	    #$end_contig_2_extra = 0;
	    my $start_contig_1_extra = $coord_c1_for_r1[1];
	    #
	    if($start_contig_1_extra < 0){
		print STDERR "Weird Stuff : Negative lengths $contig_1\t$contig_lengths{$contig_1}\t$start_contig_1_extra\t$coord_c1_for_r1[1]\n";
	    }
	    #$start_contig_1_extra = 0;

	    print OUT_EXTRA $r . "\t" . $edge_ID . "\t" . $end_contig_2_extra . "\t" . $start_contig_1_extra . "\n";
	    
	    $read_to_edge{$r}->{$edge_ID.$str_multi} = [$coord_c2_on_read[2] + $end_contig_2_extra, $coord_c1_on_read[1] - $start_contig_1_extra, $reverse_edge_r];
	    
	    #$edge_info{$edge_ID}->{"READ_OVERHANG_INFO"}->{$readname[$i]} = {"contig_1_overhang", $seq_overhang_c1, "contig_2_overhang", $seq_overhang_c2 }
	    
	}
	
	#Get the coordiate of the gap
	$edge_info{$edge_ID}->{"READ_TO_EXTRACT"}++;
    }
    
    foreach $c (keys %contig_in_gap_support){
	$support = $contig_in_gap_support{$c};
	if($support > 1){
	    #print STDERR " *** add contig in gap $edge_ID\t$c $support\n";<STDIN>; 
	    push(@{$edge_info{$edge_ID}->{"CONTIG_IN_GAP"}}, $c);
	    chop $c;#To remove the oriantation
        #print STDERR "\n$c CONTIG_IN_GAP\n";
	    $contig_in_gap_to_extract{$c} = 1;
	}
    }
}
close(FILE);
close(OUT_EXTRA);
#exit(0);

#read the contigs and extract the sequence
print STDERR " *** Reading contig file\n";
print STDERR " *** Extraction of ".((keys %contig_in_gap_to_extract) +0)." contigs in gaps\n";
open(FILE, $contig_file);
my $contig_seq = "";
my $contig_name = "";
my $extract_seq = 0;
my %all_contig_seq = ();
my $edge_ID_to_investigate = "";
while(<FILE>){
    chop $_;
    if(/^>([\w|\/|\.]*)\s*.*/){#this a new scaffold
	if($extract_seq){
	    #For the contig in gap the full sequenced is collected
	    if(exists $contig_in_gap_to_extract{$contig_name}){
		#print STDERR " *** extract sequence from contig in gap $contig_name ".(length($contig_seq))."\n";#<STDIN>;
		$contig_in_gap_to_extract{$contig_name} = "$contig_seq";
	    }
	    #Contig at at the boundaries of the gap
	    if($contig_to_edge{$contig_name}){
		foreach $edge_ID (keys %{$contig_to_edge{$contig_name}}){
		    $extraction_type = $contig_to_edge{$contig_name}->{$edge_ID};
		    $contig_length = length($contig_seq);
		    my $contig_seq_length_to_extract = 0;
		    if(exists $edge_info{$edge_ID}->{"CONTIG_OVERHANG_INFO"}->{$contig_name}){
		    	$contig_seq_length_to_extract = $edge_info{$edge_ID}->{"CONTIG_OVERHANG_INFO"}->{$contig_name} + $additional_contig_sequence;
		    }
		    else {
		    	print STDERR "**********Weird stuff******** $edge_ID not found or overhang info not found\n";
		    }
		    #the end of the contig
		    if($extraction_type eq ("1:+") || $extraction_type eq ("2:-")){
		    $contig_seq_length_to_extract = $contig_length;
			$start_sub = $contig_length - $contig_seq_length_to_extract;
			$start_sub = 0 if($start_sub < 0);
			$seq = substr $contig_seq, $start_sub, $contig_seq_length_to_extract;
			print STDERR $contig_name . "\t" .  $extraction_type . "\t" . $contig_length . "\t" . ($start_sub) . "\t" .  $contig_seq_length_to_extract . "\n";
		    }
		    #The beginning of the contig
		    if($extraction_type eq ("2:+") || $extraction_type eq ("1:-")){
		    $contig_seq_length_to_extract = $contig_length;
			$seq = substr $contig_seq, 0, $contig_seq_length_to_extract;
			print STDERR $contig_name . "\t" .  $extraction_type . "\t" . $contig_length . "\t" . 0 . "\t" .  $contig_seq_length_to_extract . "\n";
		    }
		    #The reverse complement in case of the contig have been reverse on the scaffold
		    $seq = reverse_complement($seq) if(index($extraction_type, "-") != -1);
		    #$edge_info{$edge_ID}->{"CONTIG_SEQ"} = ">$contig_name\n$seq\n";
		    

		    if($edge_ID eq $edge_ID_to_investigate){
			print STDERR " *** Contig seq extraction info: $contig_name $contig_length $contig_seq_length_to_extract $extraction_type " .(length($seq)). " " . (length($edge_info{$edge_ID}->{"CONTIG_SEQ"}))."\n";
		    }

		    my $indicator = "@";
		    if ($edge_info{$edge_ID}->{"CONTIG_SEQ"} eq ""){
			$edge_info{$edge_ID}->{"CONTIG_SEQ"} = "${indicator}$contig_name\n$seq\n";
		    }

		    else{
			$edge_info{$edge_ID}->{"CONTIG_SEQ"} .= "${indicator}$contig_name\n$seq\n";
		    }

		    $edge_info{$edge_ID}->{"CONTIG_SEQ"} .= "+\n";
		    $edge_info{$edge_ID}->{"CONTIG_SEQ"} .= "~" x length($seq);
		    $edge_info{$edge_ID}->{"CONTIG_SEQ"} .= "\n";
		    
		    if($edge_ID eq $edge_ID_to_investigate){
			print STDERR " *** Contig seq extraction info: $contig_name $contig_length $contig_seq_length_to_extract $extraction_type ".($edge_info{$edge_ID}->{"GAP_LENGTH"})." ".(length($seq)). " " . (length($edge_info{$edge_ID}->{"CONTIG_SEQ"}))."\n";#<STDIN>;
		    }
		}
	    }
	    #$all_contig_seq{$contig_name} = $contig_seq;
	}

	#print STDERR "\n$contig_name    $extract_seq   $contig_seq \n" if($extract_seq);
	$contig_name = $1;
	$contig_seq = "";
	$extract_seq = 0;
	$extract_seq = 1 if(exists $contig_to_edge{$contig_name} || exists $contig_in_gap_to_extract{$contig_name});
    }
    else{
	$contig_seq .= $_ if($extract_seq);
    }
}

#exit(0);

#Extract the sequences requires to fill the gaps in the reads
open(OUT_STATS, ">$ana_dir/gap_stats.dat") or die;
print STDERR " *** Extract gap sequences from read file $read_file\n";
my $cmd_file = "$ana_dir/cmd.txt";
open(OUT_CMD, ">$cmd_file") or die;
open(READ_FILE, $read_file) or die;
my $nb_read = 0;my $nb_copy_of_contig_sequence = 0;
my $read_seq = "";
my $quality_seq = "";
my $keep_read = 0;

my %pure_gap_seq = ();

##?## Debug purposes only
open(ERR_F, ">$ana_dir/edge_err");
open(ERR_F2, ">$ana_dir/edge_err_nf");

while(<READ_FILE>){
    #if ($nb_read > 5000){ last; }
    chop $_;
    #my $indicator;

    #if(!$is_fastq){
    #    $indicator = ">";
    #}

    #else{
    #    $indicator = "@";
    #}

    #if(/^>([\w|\/|\d|\-|\.]+)\s*.*/ or $is_fastq){

    if(/^>(.*)/ or $is_fastq){
	
	if($read_seq ne  ""){
	    #print STDERR " *** read name: |$read_name|\n";<STDIN>
	    #print STDERR "\n$read_name\n";
	    $edge_list = $read_to_edge{$read_name};
	    if(!(defined $edge_list)){
		#print STDERR "$read_name NOT DEFINED \n" if ($nb_read % 1000 == 0);
		#print STDERR "$read_name NOT DEFINED POST 10000\n" if ($nb_read > 10000);
	    }

	    if(defined $edge_list){
		#$read_seq = <FILE>;
		foreach $edge_ID (keys %{$edge_list}){
		    $start = $edge_list->{$edge_ID}->[0];
		    $end = $edge_list->{$edge_ID}->[1];
		    $reverse_complement = $edge_list->{$edge_ID}->[2];
		    #print STDERR " *** Sequence extraction ".$read_name.length($read_seq)."\t".$start."\t".($end-$start)."\n";#<STDIN>;
		    print STDERR " *** nb read processed $nb_read\n"if($nb_read % 1000 == 0);
		    $nb_read++;
		    #Get the read sequence on the rightcomplement
		    #print STDERR " *** substr $read_name\n";
		    if($end<$start){
			$pure_gap_seq{$edge_ID}->{$read_name} = $start - $end;
			#print ERR_F "$edge_ID\t$read_name\t$start\t$end\t".$edge_info{$edge_ID}->{"REF_READ"}."\t$pure_gap_seq{$edge_ID}->{$read_name}\n";
			print ERR_F "$edge_ID\t$read_name\t$start\t$end\t".$reference_reads{$edge_ID}."\t$pure_gap_seq{$edge_ID}->{$read_name}\n";
			$end = 1;
			$start = 1;
		    }


		    #my $read_seqlen_to_extract = $end - $start - 1;
		    if($end > length($read_seq)){
			print STDERR "***** Weird Stuff : Seems to exceed $read_name\t".length($seq)."\t$start\t$end\t$reverse_complement\n";
		    }
		    $seq = substr $read_seq, $start, $end-$start;
		    $seq = reverse_complement($seq) if($reverse_complement);
		    #$seq .= "\n#\n".(reverse_complement($seq)) if($reverse_complement);
		    
		    #this read support an edge twice
		    if(index($edge_ID, "[MULTI_ID") != -1){
			@tmp = split(/\[MULTI_ID/, $edge_ID);
			$edge_ID = $tmp[0];
		    }

		    $edge_info{$edge_ID}->{"COORD"} .= "$read_name\t$start\t$end\t$reverse_complement\n";
		    $edge_info{$edge_ID}->{"READ_EXTRACTED"}++;
		    
		    if($is_fastq){
			$edge_info{$edge_ID}->{"SEQ"} .= "\@$read_name\_$reverse_complement\n$seq\n";
			my $qual_seq_gap = substr $quality_seq, $start, $end-$start; 
			#my $t1 = length($read_seq);
			#my $t2 = length($quality_seq);
			#print STDERR "(SEQ) $t1, (QUAL) $t2 $edge_ID\n"; 
			$qual_seq_gap = reverse($qual_seq_gap) if ($reverse_complement);
			$edge_info{$edge_ID}->{"SEQ"} .= "+\n$qual_seq_gap\n"; 
		    }

		    else{
			$edge_info{$edge_ID}->{"SEQ"} .= "\@$read_name\_$reverse_complement\n$seq\n";
			$edge_info{$edge_ID}->{"SEQ"} .= "+\n";
			#Choose 0 as a filler quality value, seems like a value that was observed 
			#from our nanopore reads.
			$edge_info{$edge_ID}->{"SEQ"} .= "0" x ($end-$start);
			$edge_info{$edge_ID}->{"SEQ"} .= "\n";
		    }
		    
		    #All the reads of that gap have been extracted
		    #Consensus computation
		    if($edge_info{$edge_ID}->{"READ_TO_EXTRACT"} == $edge_info{$edge_ID}->{"READ_EXTRACTED"}){
			$edge_file = "$ana_dir/$edge_ID";
			
			open(OUT, ">$edge_file.fa");
			#Write the reads seuquence that belong to the gap
			print OUT $edge_info{$edge_ID}->{"SEQ"};
			#
			$c_length = 0;
			#Write the contig sequence that belong to the gap
			foreach $c (@{$edge_info{$edge_ID}->{"CONTIG_IN_GAP"}}){
			    $c_o = chop $c;
			    $seq = $contig_in_gap_to_extract{$c};
			    if($seq eq "1"){
				print STDERR "\n$c  $seq\n"; 
			    }
			    $seq = reverse_complement($seq) if($c_o eq "-");
			    $c_length += length($seq);
			    print OUT "\@$c$c_o\n$seq\n"; 
			    print OUT "+\n";
			    print OUT "~" x length($seq);
			    print OUT "\n";
			}
			
			#To collect some stats about the gaps that have to be filled
			print OUT_STATS 
			    $edge_ID."\t".
			    $edge_info{$edge_ID}->{"GAP_LENGTH"}."\t".
			    $edge_info{$edge_ID}->{"READ_TO_EXTRACT"}."\t".
			    (@{$edge_info{$edge_ID}->{"CONTIG_IN_GAP"}}+0)."\t".
			    $c_length.
			    "\n";
			
			#Write the contigs sequence at the boundaries of the gaps
			#Print multiple of the contig sequences to make sure that the consensus will not affect the contig sequences
			$nb_copy_of_contig_sequence = 1;#( $edge_info{$edge_ID}->{"READ_TO_EXTRACT"} + 1 ) / 2;
			for(my $i = 0; $i < $nb_copy_of_contig_sequence; $i++){
			    print OUT $edge_info{$edge_ID}->{"CONTIG_SEQ"};
			}
			close(OUT);
			
			open(OUT, ">$edge_file.coord");
			print OUT $edge_info{$edge_ID}->{"COORD"};
			foreach $k (keys %{$edge_info{$edge_ID}->{"CONTIG_OVERHANG_INFO"}}){
			    print OUT $k . "\t" . $edge_info{$edge_ID}->{"CONTIG_OVERHANG_INFO"}->{$k} . "\n";
			}
			close(OUT);
			
			#print OUT_CMD "$vsearch_dir/vsearch --cluster_fast $edge_file.fa -id 0 --msaout $edge_file\_cons.fa -uc $edge_file\_cons.dat\n";
			print OUT_CMD "$script_dir/print_ref.pl --seq-file $edge_file --edge-file $edge_file_info $racon_dir $graphmap_dir $water_dir \n";
			#run_exe("$qsub -l mem_free=2G,h_rt=00:05:0 -pe OpenMP 1 -N vsearch -e $ana_dir/LOG/$edge_ID.err -o $ana_dir/LOG/$edge_ID.run -b y vsearch --cluster_fast $edge_file.fa -id 0 --msaout $edge_file\_cons.fa -uc $edge_file\_cons.dat");
			#run_exe("$qsub  -l mem_free=10G,h_rt=00:05:0 -pe OpenMP 2 -N pbdagcon -e $ana_dir/LOG/$edge_ID.err -o $ana_dir/LOG/$edge_ID.run -b y $script_dir/run_pbdagcon.pl --seq-file $edge_file --bin-dir $bin_dir --isfastq $is_fastq $graphmap_dir $water_dir $pbdagcon_dir");#<STDIN>;
			delete $edge_info{$edge_ID};
			
		    }
		}
		delete $read_to_edge{$read_name}; 
	    }
	}

	#If the file is in fastq format, we assume that the sequence is on 
	#one line.
	if($is_fastq){
	    @line = split(/\s/, $_);
	    $read_name = substr($line[0],1);
	    $read_seq = <READ_FILE>;
	    <READ_FILE>;
	    $quality_seq = <READ_FILE>;
	    chomp($read_seq);
	    chomp($quality_seq);
	}

	else{
	    my @line = split(/\s/, $1);
	    $read_name = $line[0]; 
	    $read_seq = "";
	}
    }

    else{
        if(/^[ATCG]/){
            $keep_read = 1;
        }

        else{
            $keep_read = 0;
        }

        if($keep_read){
            $read_seq .= $_;
        }

    }
}
close(ERR_F);
close(READ_FILE);
close(OUT_STATS);







#Run the sequence consensus calling
print STDERR " *** Run racon for sequence consensus calling\n";
#run_exe("cat $cmd_file | xargs -L 1 -P $nb_process bash &> /dev/null \n");
#run_exe("cat $cmd_file | xargs -L 1 -P $nb_process bash \n");
#Extract the first read from edge read info can be run in parallel
run_exe("cat $cmd_file | xargs -L 1 -P $nb_process -I COMMAND sh -c \"COMMAND\" 2> $run_racon_log");
#to produce the filled gap scaffolds
write_filled_gap_scaffold();

#Not working right now, TODO fix it.
#my $cmd = "ls $ana_dir | grep \"\\\-vs\\\-\" | xargs rm";
#run_exe($cmd);
close(ERR_F2);

#run_exe("$script_dir/run_racon_4.pl --seq-file $outfile --read-file $read_file $racon_dir $graphmap_dir $water_dir");



sub read_scaffold_file{
    my ($scaffold_file) = @_;
    open(FILE, $scaffold_file);
    my ($contig_1, $contig_1_ori, $contig_2, $contig_2_ori);
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	if(index($line[0], ">") != -1){#this a new scaffold
	    $line_str = <FILE>;#read the first contig
	    chop $line_str;
	    @line = split(/\t/, $line_str);
	    #print "*** @line\n";<STDIN>;
	    $contig_1 = $line[0];
	    $contig_1_ori = "+"; $contig_1_ori = "-" if($line[1] eq "EB");
	    $gap_length = $line[3];
	    $contig_to_edge{$contig_1} = {};
	}
	else{
	    $contig_2 = $line[0];
	    $contig_to_edge{$contig_2} = {};
	    
	    $contig_2_ori = "+"; $contig_2_ori = "-" if($line[1] eq "EB");
	    
	    #Get the edge_ID and corresponding representative orientation
	    $tmp_array = get_edge_ID($contig_1, $contig_1_ori, $contig_2, $contig_2_ori);
	    $edge_ID = $tmp_array->[0];
	    $edge_ori = $tmp_array->[1];
	    
	    #print STDERR " *** EDGE tested $edge_ID $gap_length\n";
	    
	    #print STDERR " *** EDGE selected $edge_ID $gap_length\n";<STDIN>;
	    $edge_info{$edge_ID} = {"GAP_LENGTH", $gap_length, "READ_TO_EXTRACT", 0, "READ_EXTRACTED", 0, "SEQ", "", "COORD", "", "CONTIG_SEQ", "", "EDGE_ORI", $edge_ori, "CONTIG_GAPFILLING", [], "CONTIG_OVERHANG_INFO", {}};#, "REF_READ", ""};# "REF_COORDS", {}};#be carreful with conflicting edges that have the same contig pair
	    $nb_selected_edges++;
	    $contig_to_edge{$contig_1}->{$edge_ID} = "1:$contig_1_ori";
	    $contig_to_edge{$contig_2}->{$edge_ID} = "2:$contig_2_ori";
	    
	    #For the next gap
	    $contig_1 = $contig_2;
	    $contig_1_ori = $contig_2_ori;
	    $gap_length = $line[3];
	}
    }
    close(FILE);
}

sub write_filled_gap_scaffold{
    #need contig sequences, scaf file and smith-waterman mapping of contigs to consensus
    
    #Get all the the contig sequence $contig_to_edge should have been already computed
    my %all_contig_seq = ();
    open(FILE, $contig_file);
    print STDERR " *** Reading contig file\n";
    my ($contig_seq, $contig_name); $contig_name = "";
    while(<FILE>){
	chop $_;
	if(/^>([\w|\/|\.]*)\s*.*/){#this a new scaffold
	    #print STDERR " *** $contig_name\n";<STDIN>;
	    $all_contig_seq{$contig_name} = $contig_seq if($contig_name ne "");
	    #
	    $contig_name = $1;
	    $contig_seq = "";
	}
	else{
	    $contig_seq .= $_;
	}
    }
    $all_contig_seq{$contig_name} = $contig_seq;#the last contig
    close(FILE);

    #read scaffold file
    open(OUT, ">$outfile");    
    #run_exe("rm $ana_dir/gap_substring_report.dat;touch $ana_dir/gap_substring_report.dat");
    open(OUT_S, ">$ana_dir/gap_size.dat");
    open(FILE, $scaffold_file);
    my ($contig_1, $contig_1_ori, $contig_2, $contig_2_ori, $scaff_name, $gap_seq);
    $contig_2 = "";$contig_1 = "";
    my $gap_filled_stats = "";
    my $contig_trimming_length = 0;#This is used to remove starting sequences of a contig when 2 contigs are overlapping
    print STDERR " *** Read the scaffold file and fill gaps\n";
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	$gap_seq = "";
	if(index($line[0], ">") != -1){#this a new scaffold
	    
	    #Write the sequence of the last contig and the singleton scaffolds
	    if($contig_2 ne ""){
		print OUT "".(get_contig_seq($contig_2, $contig_2_ori, \%all_contig_seq))."\n";
	    }
	    else{
		if($contig_1 ne ""){#this is a singleton scaffold
		    print OUT "".(get_contig_seq($contig_1, $contig_1_ori, \%all_contig_seq))."\n";
		}
	    }
	    #Write the scaffold file name (need to update the sequence length)
	    print OUT $_."\n";

	    $line_str = <FILE>;#read the first contig
	    chop $line_str;
	    @line = split(/\t/, $line_str);
	    #print "*** @line\n";<STDIN>;
	    $contig_1 = $line[0];
	    $contig_1_ori = "+"; $contig_1_ori = "-" if($line[1] eq "EB");
	    $gap_length = $line[3];
	    $contig_2 = "";
	}
	else{
	    $contig_2 = $line[0];
	    $contig_2_ori = "+"; $contig_2_ori = "-" if($line[1] eq "EB");
	    @contig_order = sort($contig_1, $contig_2);
	    
	    #Write contig 1
	    $contig_1_seq = get_contig_seq($contig_1, $contig_1_ori, \%all_contig_seq);
	    #Trim the contig in case of contig overlap (always trim on the left side)
	    if($contig_trimming_length != 0){
		print STDERR " ***** Trimming $contig_1 $contig_trimming_length\n";
		$contig_1_seq = substr $contig_1_seq, $contig_trimming_length;
	    }
	    print OUT $contig_1_seq;
	    
	    #print STDERR " *** Fill gap between $contig_1 and $contig_2\n";#<STDIN>;
	    #print STDERR " *** $contig_1 $contig_1_ori\n$contig_1_seq\n";<STDIN>;

	    #Write the gap sequence
	    $edge_ID = join("-vs-", @contig_order);
	    $gap_seq = extract_gap_seq($edge_ID);

	    if($gap_seq eq ""){
		print STDERR "**** REF NOT FOUND ($edge_ID) *****\n\n";
	    }
	    
	    $gap_filled_stats = "NA";#the gap is not filled
	    my $gap_diff_large = "";
	    if($gap_seq ne ""){
		#this is an overlap or a no gap
		
		if ( $gap_seq =~ /^[0-9,.E]+$/ ) { 
		    $contig_trimming_length = $gap_seq;
		    $gap_filled_stats = -($contig_trimming_length);
		    #$contig_trimming_length = 0;
		}
		else{
		    #This a gap sequence
		    $gap_filled_stats = length($gap_seq);
		    $contig_trimming_length = 0;
		    #Formula for taking out edges that are way off from the predicted gap size. 
		    #Formula works out to : skip filling if predicted or actual gap size is greater than 3 times the other by 250.
		    # if((abs($gap_filled_stats - $gap_length)/(abs($gap_length + $gap_filled_stats) + 500)) < 0.5){
		    #     print OUT $gap_seq;
		    # }
		    # 
		    # elsif($gap_length > 0){
		    #     print OUT ( "N" x $gap_length);
		    #     $gap_diff_large = "BIGDIFF";
		    # }
		    # }
		    
		    
		    print OUT $gap_seq;
		}
		
	    }
	    #For some reason the gap was not filled
	    else{
		$contig_trimming_length = 0;
		if($gap_length > 0){print OUT ( "N" x $gap_length);}
		else{print OUT ( "N" x 3);}
	    } 
	    #$contig_trimming_length = 0;
	    #
	    #Get statistics
	    print OUT_S $edge_ID."\t".$gap_length."\t".$gap_filled_stats."\t".$gap_diff_large."\n";

	    #For the next gap
	    $contig_1 = $contig_2;
	    $contig_1_ori = $contig_2_ori;
	    $gap_length = $line[3];
	}
    }
    #The last contig of the last scaffold
    if($contig_2 ne ""){
	print OUT "".(get_contig_seq($contig_2, $contig_2_ori, \%all_contig_seq))."\n";
    }
    else{
	if($contig_1 ne ""){#this is a singleton scaffold
	    print OUT "".(get_contig_seq($contig_1, $contig_1_ori, \%all_contig_seq))."\n";
	}
    }
    

    close(FILE);
    close(OUT);
    close(OUT_S);
    
}



sub extract_gap_seq{
    my ($file_ID) = @_;
    
    #get the cordinate of the sequence that have to be extracted
    my $mapping_file = "$ana_dir/$file_ID\_r.fa";
    return "" if(! -e $mapping_file);
    
    print STDERR " *** get Gap sequence $file_ID $mapping_file\n";#<STDIN>;

    open(REF_FILE, $mapping_file);
    my $seq ="";
    <REF_FILE>;
    chop($seq = <REF_FILE>);

    if($seq ne ""){
    	return $seq;
    }
    else{
    	#my $out = $pure_gap_seq{$file_ID}->{$edge_info{$file_ID}->{"REF_READ"}};
    	my $out = $pure_gap_seq{$file_ID}->{$reference_reads{$file_ID}};

	#print ERR_F2 "$file_ID\t".$edge_info{$file_ID}->{"REF_READ"}."\t".$pure_gap_seq{$file_ID}->{$edge_info{$file_ID}->{"REF_READ"}}."\n";
	print ERR_F2 "$file_ID\t".$reference_reads{$file_ID}."\t".$pure_gap_seq{$file_ID}->{$reference_reads{$file_ID}}."\n";
    	if(! (defined $out) ) {
	    #print STDERR "*************WEIRD STUFF : Size entry absent for overlap*************".$file_ID."\t".$edge_info{$file_ID}->{"REF_READ"}."\n";
	    print STDERR "*************WEIRD STUFF : Size entry absent for overlap*************".$file_ID."\t".$reference_reads{$file_ID}."\n";
    	}
    	return $out;
    }
}


#The convention for the edge ID
#This is not perfect and should be clean up
#Contig are ordered according to their name
#the orientatin of the 1st contig is kept
sub get_edge_ID{
    my ($c1, $o1, $c2, $o2) = @_;
    @contig_order = sort($c1, $c2);
    my @res = (); $res[0] = join("-vs-", @contig_order);
    if($c1 eq $contig_order[0]){
	push(@res, $o1);
    }
    else{
	push(@res, $o2);
    }
    return \@res;
}

sub get_contig_seq{
    my ($contig_name, $contig_ori, $all_contig_seq) = @_;
    $seq = $all_contig_seq->{$contig_name};
    $seq = reverse_complement($seq) if($contig_ori eq "-");
    return $seq;
}


sub reverse_complement{
    my ($seq) = @_;
    #print STDERR " ----- reverse_complement $seq\n";
    my @seq_a = split(//, scalar reverse $seq);
    #print STDERR " ----- reverse_complement @seq_a\n";
    for(my $i = 0; $i < @seq_a; $i++){
 	$seq_a[$i] = opposite_base($seq_a[$i]);
    }
    #print STDERR " ----- reverse_complement @seq_a\n";<STDIN>;
    return (join("", @seq_a));
}

sub opposite_base{
    (my $b) = @_;
    switch($b){
	case "A" {return "T";}
	case "C" {return "G";}
	case "G" {return "C";}
	case "T" {return "A";}
	case "N" {return "N";}
    }
}


sub reverse_ori{
    my ($ori) = @_;
    return "+" if($ori eq "-");
    return "-" if($ori eq "+");
}

sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";;
    print STDERR `$exe` if($run);
}

#vsearch --cluster_fast S288C/FINAL_TEST/GAP_FILLING/171068_171204.fa -id 0 --msaout cons.txt -uc cons.dat

#

#head -n2 S288C/FINAL_TEST/GAP_FILLING/171068_171204.fa > S288C/FINAL_TEST/GAP_FILLING/171068_171204.ref.fa;tail -n82 S288C/FINAL_TEST/GAP_FILLING/171068_171204.fa > S288C/FINAL_TEST/GAP_FILLING/171068_171204.querry.fa
#blasr -nproc 1 -m 1 -minMatch 5 -bestn 10 -noSplitSubreads -advanceExactMatches 1 -nCandidates 1 -maxAnchorsPerPosition 1 -sdpTupleSize  7 S288C/FINAL_TEST/GAP_FILLING/171068_171204.ref.fa S288C/FINAL_TEST/GAP_FILLING/171068_171204.querry.fa

#./OPERA-LG_v2.0.6/bin/extract_read_sequence.pl S288C/FINAL_TEST/edge_read_info.dat S288C/FINAL_TEST/results/scaffolds.scaf /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/test_dataset_long_reads/nanopore.fa S288C/FINAL_TEST/GAP_FILLING

#awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' /home/bertrandd/PROJECT_LINK/OPERA_LG/NANOPORE/DATA/W303_ONT_Raw_reads_filter_rename.fa  > /home/bertrandd/PROJECT_LINK/OPERA_LG/NANOPORE/DATA/W303_ONT_Raw_reads_filter_rename_single_line.fa


#172500	BE	1998	-100
#174658	BE	5735	17

#WHY THIS GAP IS UNFILLED
#more S288C/V6_RAW/GAP_FILLING/172838_174704_contig_cons.m5

#No file produced
#S288C/V6_RAW_400/GAP_FILLING/172440_175582

#CONTIG	opera_scaffold_2_length__921946_cov__799.2	981183
#538147	645501	1	107229	tpg_BK006938.2___organism_Saccharomyces_cerevisiae_S288c___strain_S288c___moltype_genomic___chromosome_IV___note_R64-1-1_	opera_scaffold_2_length__921946_cov__799.2	99.4	
#relocation, inconsistency = 5927
#651430	668288	107231	124073	tpg_BK006938.2___organism_Saccharomyces_cerevisiae_S288c___strain_S288c___moltype_genomic___chromosome_IV___note_R64-1-1_	opera_scaffold_2_length__921946_cov__799.2	99.49	
#relocation, inconsistency = -5129
#668295	677852	129209	138778	tpg_BK006938.2___organism_Saccharomyces_cerevisiae_S288c___strain_S288c___moltype_genomic___chromosome_IV___note_R64-1-1_	opera_scaffold_2_length__921946_cov__799.2	99.53	
#local misassembly

#In /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/MOCK_20/ASSEMBLY/OPERA_MS/FULL_SAMPLE/MEGAHIT/NANOPORE/contigs/scaffolds-long-read-latest
#/mnt/projects/bertrandd/opera_lg/OPERA_LONG_READ/OPERA-LG_v2.1.0/bin/extract_read_sequence_xargs.pl ../OPERA-long-read/edge_read_info.dat /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/MOCK_20/ASSEMBLY/MEGAHIT/final.contigs_soap.fa scaffolds.scaf /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/NANOPORE/LIBRARY/POOL/POOL_all/POOL.fa GAP_FILLING scaffoldSeq.fasta.filled ALL 20246


