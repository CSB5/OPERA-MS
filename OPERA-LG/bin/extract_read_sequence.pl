#!/usr/bin/perl
use warnings;
use Switch;
use Getopt::Long;

my ($edge_file_info, $contig_file, $scaffold_file, $read_file, $ana_dir, $outfile, $stage) = @ARGV;

$stage = "ALL" if(! defined $stage);

my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -V";

#length of the mapping that include the contig sequence
my $sequence_overhang = 400;
my $additional_contig_sequence = 100;

#For only getting the consensus
if($stage eq "CONS"){
    print STDERR " *** Construct the gap consensus sequences\n";
    run_exe("rm -rf $ana_dir/LOG;mkdir $ana_dir/LOG");
    opendir(DIR, $ana_dir);
    my @all_cons = readdir(DIR);
    close(DIR);
    
    foreach $file (@all_cons){
	next if(index($file, ".fa") == -1 || index($file, "con") != -1 || index($file, "_q.fa") != -1 || index($file, "_r.fa") != -1);
	$edge_ID = $file;chop $edge_ID;chop $edge_ID;chop $edge_ID;
	#run_exe("/home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/run_pbdagcon.pl $ana_dir/$name");
	run_exe("$qsub -l mem_free=10G,h_rt=00:30:0 -pe OpenMP 1 -N vsearch -e $ana_dir/LOG/$edge_ID.err -o $ana_dir/LOG/$edge_ID.run -b y vsearch --cluster_fast $ana_dir/$edge_ID.fa -id 0 --msaout $ana_dir/$edge_ID\_cons.fa -uc $ana_dir/$edge_ID\_cons.dat");#<STDIN>;
	#run_exe("$qsub  -l mem_free=10G,h_rt=00:30:0 -pe OpenMP 2 -N pbdagcon -e $ana_dir/LOG/$edge_ID.err -o $ana_dir/LOG/$edge_ID.run -b y /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/run_pbdagcon.pl $ana_dir/$edge_ID");#<STDIN>;
	run_exe("sleep 0.5");
    }

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
#run_exe("rm -r $ana_dir");
run_exe("rm -rf $ana_dir/LOG");
run_exe("mkdir -p $ana_dir/LOG");


#Structure that indicate sequence required to fill the gap between the 2 contigs involvolved in the edge
my %edge_info = ();
#indicate which reads belong to which edges
my %read_to_edge = ();
#indicate which contig belong to which edges
my %contig_to_edge = ();

#Read the scaffold file to get the set of adjacent edges used to construct the scaffold
my $nb_selected_edges = 0;
my $cmp_multi_edge = 1;my $str_multi = "";

#open(FILE, "head -n5 $scaffold_file |");
read_scaffold_file($scaffold_file);

print STDERR " *** nb edges selected for gapfilling $nb_selected_edges\n";#<STDIN>;

#Get the read and the coordinates of the gaps
open(FILE, $edge_file_info);
my ($reverse_edge, $reverse_edge_r);
my %contig_in_gap_to_extract = ();#list of contig for which the full sequence need to be stored
my %contig_in_gap_support = ();#To compute the support of the contig in the gap, the contig should be at least supported by 2 reads
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
	#print STDERR " *** edge ID ".$edge_ID." ".$edge_ori."\n";#<STDIN>;
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

    for(my $i = 0; $i < @read_name; $i++){
	$r = $read_name[$i];
	$reverse_edge_r = $reverse_edge;
	
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
	if($coord_c1_on_read[2] < $coord_c2_on_read[2]){#Contig 1 is first on the read
	    $read_to_edge{$r}->{$edge_ID.$str_multi} = [$coord_c1_on_read[2] - $sequence_overhang, $coord_c2_on_read[1] + $sequence_overhang, $reverse_edge_r];
	}
	else{
	    $read_to_edge{$r}->{$edge_ID.$str_multi} = [$coord_c2_on_read[2] - $sequence_overhang, $coord_c1_on_read[1] + $sequence_overhang, $reverse_edge_r];
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
	    $contig_in_gap_to_extract{$c} = 1;
	}
    }
}
close(FILE);

#read the contigs and extract the sequence
print STDERR " *** Reading contig file\n";
print STDERR " *** Extraction of ".((keys %contig_in_gap_to_extract) +0)." contigs in gaps\n";
open(FILE, $contig_file);
my $contig_seq = "";
my $contig_name = "";
my $extract_seq = 0;
my %all_contig_seq = ();
my $contig_seq_length_to_extract = $sequence_overhang + $additional_contig_sequence;
my $edge_ID_to_investigate = "";
while(<FILE>){
    chop $_;
    if(/^>([\w|\/]*)\s*.*/){#this a new scaffold
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
		    #the end of the contig
		    if($extraction_type eq ("1:+") || $extraction_type eq ("2:-")){
			$seq = substr $contig_seq, $contig_length - $contig_seq_length_to_extract, $contig_seq_length_to_extract;
		    }
		    #The beginning of the contig
		    if($extraction_type eq ("2:+") || $extraction_type eq ("1:-")){
			$seq = substr $contig_seq, 0, $contig_seq_length_to_extract;
		    }
		    #The reverse complement in case of the contig have been reverse on the scaffold
		    $seq = reverse_complement($seq) if(index($extraction_type, "-") != -1);
		    #$edge_info{$edge_ID}->{"CONTIG_SEQ"} = ">$contig_name\n$seq\n";
		    

		    if($edge_ID eq $edge_ID_to_investigate){
			print STDERR " *** Contig seq extraction info: $contig_name $contig_length $contig_seq_length_to_extract $extraction_type " .(length($seq)). " " . (length($edge_info{$edge_ID}->{"CONTIG_SEQ"}))."\n";
		    }

		    if($edge_info{$edge_ID}->{"CONTIG_SEQ"} eq ""){
			#$edge_info{$edge_ID}->{"CONTIG_SEQ"} = "$seq";
			$edge_info{$edge_ID}->{"CONTIG_SEQ"} = ">$contig_name\n$seq\n";
		    }
		    else{
			$edge_info{$edge_ID}->{"CONTIG_SEQ"} .= ">$contig_name\n$seq\n";
			#Add the sequence with the gap
			#if(index($extraction_type, "2") == -1){
			#$edge_info{$edge_ID}->{"CONTIG_SEQ"} .= ("N" x ($edge_info{$edge_ID}->{"GAP_LENGTH"})).$seq;
			#}
			#else{
			#$edge_info{$edge_ID}->{"CONTIG_SEQ"} = $seq.("N" x ($edge_info{$edge_ID}->{"GAP_LENGTH"})).$edge_info{$edge_ID}->{"CONTIG_SEQ"};
			#}
		    }
		    if($edge_ID eq $edge_ID_to_investigate){
			print STDERR " *** Contig seq extraction info: $contig_name $contig_length $contig_seq_length_to_extract $extraction_type ".($edge_info{$edge_ID}->{"GAP_LENGTH"})." ".(length($seq)). " " . (length($edge_info{$edge_ID}->{"CONTIG_SEQ"}))."\n";#<STDIN>;
		    }
		}
	    }
	    #$all_contig_seq{$contig_name} = $contig_seq;
	}
	$contig_name = $1;
	$contig_seq = "";
	$extract_seq = 0;
	$extract_seq = 1 if(exists $contig_to_edge{$contig_name} || exists $contig_in_gap_to_extract{$contig_name});
    }
    else{
	$contig_seq .= $_ if($extract_seq);
    }
}

#Extract the sequences requires to fill the gaps in the reads
open(OUT_STATS, ">$ana_dir/gap_stats.dat");
print STDERR " *** Extract gap sequences from read file $read_file\n";
open(FILE, $read_file);
my $nb_read = 0;my $nb_copy_of_contig_sequence = 0;
while(<FILE>){
    if(/^>([\w|\/|\d|\-|\.]+)\s*.*/){
	$read_name = $1;
	#print STDERR " *** read name: |$read_name|\n";<STDIN>
	$edge_list = $read_to_edge{$read_name};
	next if(!defined $edge_list);
	#The extraction is brute force and expect to have read sequence writen on a single line
	$read_seq = <FILE>;
	foreach $edge_ID (keys %{$edge_list}){
	    $start = $edge_list->{$edge_ID}->[0];
	    $end = $edge_list->{$edge_ID}->[1];
	    $reverse_complement = $edge_list->{$edge_ID}->[2];
	    #print STDERR " *** Sequence extraction ".$read_name.length($read_seq)."\t".$start."\t".($end-$start)."\n";#<STDIN>;
	    print STDERR " *** nb read processed $nb_read\n"if($nb_read % 1000 == 0);
	    $nb_read++;
	    #Get the read sequence on the rightcomplement
	    #print STDERR " *** substr $read_name\n";
	    $seq = substr $read_seq, $start, $end-$start;
	    $seq = reverse_complement($seq) if($reverse_complement);
	    #$seq .= "\n#\n".(reverse_complement($seq)) if($reverse_complement);
	    
	    #this read support an edge twice
	    if(index($edge_ID, "[MULTI_ID") != -1){
		@tmp = split(/\[MULTI_ID/, $edge_ID);
		$edge_ID = $tmp[0];
	    }

	    $edge_info{$edge_ID}->{"SEQ"} .= ">$read_name\_$reverse_complement\n$seq\n";
	    $edge_info{$edge_ID}->{"COORD"} .= "$read_name\t$start\t$end\t$reverse_complement\n";
	    $edge_info{$edge_ID}->{"READ_EXTRACTED"}++;
	    
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
		    $seq = reverse_complement($seq) if($c_o eq "-");
		    $c_length += length($seq);
		    print OUT ">$c$c_o\n$seq\n";
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
		close(OUT);

		run_exe("$qsub -l mem_free=10G,h_rt=00:30:0 -pe OpenMP 1 -N vsearch -e $ana_dir/LOG/$edge_ID.err -o $ana_dir/LOG/$edge_ID.run -b y vsearch --cluster_fast $edge_file.fa -id 0 --msaout $edge_file\_cons.fa -uc $edge_file\_cons.dat");
		run_exe("$qsub  -l mem_free=10G,h_rt=00:30:0 -pe OpenMP 2 -N pbdagcon -e $ana_dir/LOG/$edge_ID.err -o $ana_dir/LOG/$edge_ID.run -b y /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/run_pbdagcon.pl $edge_file");#<STDIN>;
		delete $edge_info{$edge_ID};
	    }
	}
	delete $read_to_edge{$read_name};
    }
}
close(FILE);
close(OUT_STATS);

#to produce the filled gap scaffolds
#write_filled_gap_scaffold();

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
	    
	    if(1 || $gap_length > 9){
		#print STDERR " *** EDGE selected $edge_ID $gap_length\n";<STDIN>;
		$edge_info{$edge_ID} = {"GAP_LENGTH", $gap_length, "READ_TO_EXTRACT", 0, "READ_EXTRACTED", 0, "SEQ", "", "COORD", "", "CONTIG_SEQ", "", "EDGE_ORI", $edge_ori, "CONTIG_GAPFILLING", []};#be carreful with conflicting edges that have the same contig pair
		$nb_selected_edges++;
		$contig_to_edge{$contig_1}->{$edge_ID} = "1:$contig_1_ori";
		$contig_to_edge{$contig_2}->{$edge_ID} = "2:$contig_2_ori";
	    }

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
	if(/^>([\w|\/]*)\s*.*/){#this a new scaffold
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
    run_exe("rm $ana_dir/gap_substring_report.dat;touch $ana_dir/gap_substring_report.dat");
    open(OUT_S, ">$ana_dir/gap_size.dat");
    open(FILE, $scaffold_file);
    my ($contig_1, $contig_1_ori, $contig_2, $contig_2_ori, $scaff_name, $gap_seq);
    $contig_2 = "";$contig_1 = "";
    my $gap_filled_stats = "";
    my $contig_trimming_length = 0;#This is used to remove stating sequences of a contig when 2 contigs are overlapping
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
	    #Trim the contig in case of contig overlap
	    if($contig_trimming_length != 0){
		$contig_1_seq = substr $contig_1_seq, $contig_trimming_length;
	    }
	    print OUT $contig_1_seq;
	    
	    print STDERR " *** Fill gap between $contig_1 and $contig_2\n";#<STDIN>;
	    #print STDERR " *** $contig_1 $contig_1_ori\n$contig_1_seq\n";<STDIN>;

	    #Write the gap sequence
	    $edge_ID = join("_", @contig_order);
	    $gap_seq = extract_gap_seq($edge_ID);
	    $gap_seq = extract_gap_seq_vsearch($edge_ID) if($gap_seq eq "");#padagcon was not able to get the gap we use vsearch concensus to fill it
	    $gap_filled_stats = "NA";#the gap is not filled
	    if($gap_seq ne ""){
		#this is an overlap or a no gap
		if ( $gap_seq =~ /^[0-9,.E]+$/ ) { 
		    $contig_trimming_length = $gap_seq;
		    $gap_filled_stats = -($contig_trimming_length);
		}
		else{
		    #This a gap sequence
		    print OUT $gap_seq;
		    $gap_filled_stats = length($gap_seq);
		}
	    }
	    #For some reason the gap was not filled
	    else{
		if($gap_length > 0){print OUT ( "N" x $gap_length);}
		else{print OUT ( "N" x 3);}
	    }
	    $contig_trimming_length = 0;
	    #
	    #Get statistics
	    print OUT_S $edge_ID."\t".$gap_length."\t".$gap_filled_stats."\n";

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

sub extract_gap_seq_vsearch{
    my ($file_ID) = @_;

    print STDERR " *** get vsearch gap sequence $file_ID\n";#<STDIN>;

    my $mapping_file = "$ana_dir/$file_ID\_cons.fa";
    return "" if(! -e $mapping_file);
    #
    #Extract the contig name from the file ... problem with _ ... -_-'
    my @contig_name_tab = split(/_/, $file_ID);
    #NEED TO CHANGE THAT MESS
    @contig_name_tab = ($contig_name_tab[0]."_".$contig_name_tab[1], $contig_name_tab[2]."_".$contig_name_tab[3]) if(@contig_name_tab > 2);
    #
    #
    #Construct the multiple alignement between the contig and the concensus
    my @multi_alignement = ();
    my $contig_name = "";
    my $flag_extract_seq = -1;
    
    open(FILE_CONS, "$mapping_file");
    while(<FILE_CONS>){
	chop $_;
	if ( $_ =~ /\>(.*)/ ){
	    $contig_name = $1;
	    #print STDERR $contig_name."\n";
	    $flag_extract_seq = -1;
	    if($contig_name eq $contig_name_tab[0] || $contig_name eq $contig_name_tab[1] || $contig_name eq "consensus"){
		push(@multi_alignement, []);
		$flag_extract_seq = @multi_alignement - 1;
		#print STDERR " **** EXT $contig_name $flag_extract_seq\n";
	    }
	}
	else{
	    if($flag_extract_seq != -1){
		@str = split(//, $_);
		push(@{$multi_alignement[$flag_extract_seq]}, @str);
	    }
	}
    }
    close(FILE_CONS);

    #Extract the gap sequence
    my ($left, $right);$left = "";$right = "";
    my $alignement_length = @{$multi_alignement[0]} + 0;
    my $gap_seq = "";
    my $begin_overlap = -1; my $end_overlap = -1;
    my $flag_first_rigth = 0;
    for(my $i = 0; $i < $alignement_length; $i++){
	
	#Identification of the boundary order -_-'
	if($left eq ""){
	    if($multi_alignement[0]->[$i] ne "+" && $multi_alignement[0]->[$i] ne "-"){
		$left = 0;
		$right = 1;
	    }
	    if($multi_alignement[1]->[$i] ne "+" && $multi_alignement[1]->[$i] ne "-"){
		$left = 1;
		$right = 0;
	    }
	}
        
	#print STDERR "$multi_alignement[0]->[$i]$multi_alignement[1]->[$i]$multi_alignement[2]->[$i]\nleft $left right $right\n";#<STDIN>;

	next if($left eq "");

	#It is an overlap
	if($multi_alignement[$left]->[$i] ne "+" && $multi_alignement[$left]->[$i] ne "-"){
	    if($multi_alignement[$right]->[$i] ne "+" && $multi_alignement[$right]->[$i] ne "-"){
		$begin_overlap = $i if($begin_overlap == -1);
		$end_overlap = $i;
	    }
	    else{
		$gap_seq = "";
	    }
	}
	else{
	    next if($begin_overlap != -1);
	    if($multi_alignement[$right]->[$i] ne "+"&& $multi_alignement[$right]->[$i] ne "-"){
		last if($gap_seq ne "");#We have got the gap sequence
	    }
	    else{
		$gap_seq .= $multi_alignement[2]->[$i] if($multi_alignement[2]->[$i] ne "+" && $multi_alignement[2]->[$i] ne "-");#Add the concensus sequence to the gap
	    }
	}
    }


    if($begin_overlap != -1){
	$gap_seq = $end_overlap - $begin_overlap;
    }

    return $gap_seq;
    
}


sub extract_gap_seq{
    my ($file_ID) = @_;
    
    #get the cordinate of the sequence that have to be extracted
    my $mapping_file = "$ana_dir/$file_ID\_contig_cons.m5";
    return "" if(! -e $mapping_file);
    
    print STDERR " *** get Gap sequence $file_ID $mapping_file\n";#<STDIN>;

    open(FILE_CONS, $mapping_file);
    my $tmp = "";
    #
    my @tmp_tab = split(/_/, $file_ID);
    #
    #Extract the contig name from the file ... problem with _ ... -_-'
    my @contig_name_tab = split(/_/, $file_ID);
    #NEED TO CHANGE THAT MESS
    @contig_name_tab = ($contig_name_tab[0]."_".$contig_name_tab[1], $contig_name_tab[2]."_".$contig_name_tab[3]) if(@contig_name_tab > 2);
    #
    my $contig_name = "";
    my $ref_name = "";
    my %contig_coord = ($contig_name_tab[0], [-1,-1,-1,-1], $contig_name_tab[1], [-1,-1,-1,-1]);
    my $flag_begin = 1;
    my $contig_seq_length = $sequence_overhang + $additional_contig_sequence; 
    #Read the alignement of the contig to the consensu
    while(<FILE_CONS>){
	chop $_;
	#this is simply comments
	if(index($_, "#") != -1){
	    $flag_begin = 1;
	    next;
	}
	#Read the 3 lines of the alignement
	@ref_seq = split(/\s+/, $_);    
	next if(@ref_seq < 2);#this an empty line

	$tmp = <FILE_CONS>;
	$tmp = <FILE_CONS>;chop $tmp;@contig_seq =  split(/\s+/, $tmp);
	
	$contig_name = $contig_seq[0];
	
	#print STDERR " **** the contig alignement line $tmp ===> $contig_name\n";

	#extract the 1 and last corrdinates of the alignements
	if($flag_begin){
	    $contig_coord{$contig_name}->[0] = $ref_seq[1];
	    $contig_coord{$contig_name}->[1] = $contig_seq[1];
	    $flag_begin = 0;
	}
	else{
	    $contig_coord{$contig_name}->[2] = $ref_seq[3];
	    $contig_coord{$contig_name}->[3] = $contig_seq[3];
	    $flag_begin = 0
	}
    }
    close(FILE_CONS);

    open(OUT_REP, ">>$ana_dir/gap_substring_report.dat");
    print OUT_REP 
	$file_ID."\t".
	$contig_name_tab[0]."\t".(join(";", @{$contig_coord{$contig_name_tab[0]}}))."\t".
	$contig_name_tab[1]."\t".(join(";", @{$contig_coord{$contig_name_tab[1]}}));
    
    my $contig_gap_1 = $contig_name_tab[0];
    my $contig_gap_2 = $contig_name_tab[1];
    
    #Conare the end mapping of the 2 contig and re-order them
    if($contig_coord{$contig_gap_2}->[2] < $contig_coord{$contig_gap_1}->[2]){
	$contig_gap_1 = $contig_name_tab[1];
	$contig_gap_2 = $contig_name_tab[0];
    }

    #beginning of the gap
    $gap_start = $contig_coord{$contig_gap_1}->[2];#alignemnt end of contig 1
    $gap_end = $contig_coord{$contig_gap_2}->[0];#alignemnt start of contig 2
    
    print OUT_REP "\t"."[ $gap_start $gap_end ]\n";
    close(OUT_REP);
    
    my $gap_seq_res = "";
    #CHECK FOR WEIRD CASE
    #Case were only one contig map ... Need to add more cases here
    #print STDERR " *** gap alignement boundaries ".$gap_start." ".$gap_end."\n";<STDIN>;
    return $gap_seq_res if($gap_start == -1 || $gap_end == -1);

    #Check between contig overlap and proper gaps
    if($gap_end <= $gap_start){
	#This is an over lap or a no gap sequence
	$gap_seq_res = $gap_start - $gap_end;
    }
    else{
	#This is a proper gap and we extract the sequence from the consensus sequence
	my $seq = `tail -n1 $ana_dir/$file_ID\_pbdagcon.fa`;
	my $seq_length = length($seq);
	$gap_seq_res = (substr $seq, $gap_start , $gap_end - $gap_start);
    }

    #print STDERR $gap_seq_res."\n";<STDIN>;

    return $gap_seq_res;
}


#The convention for the edge ID
#This is not perfect and should be clean up
#Contig are ordered according to their name
#the orientatin of the 1st contig is kept
sub get_edge_ID{
    my ($c1, $o1, $c2, $o2) = @_;
    @contig_order = sort($c1, $c2);
    my @res = (); $res[0] = join("_", @contig_order);
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


