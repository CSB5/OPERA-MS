use warnings;
use Switch;
use Getopt::Long;


my $edge_file_info;
my $contig_file;
my $scaffold_file;
my $read_file;
my $ana_dir;
my $is_fastq = 0;

GetOptions(
    "edge-file=s"    => \$edge_file_info,
    "contig-file=s"    => \$contig_file,
    "opera-lr-dir=s" => \$opera_lr_dir,
    #"scaffold-file=s" => \$scaffold_file,
    "read-file=s"     => \$read_file,
    #
    "output-directory=s" => \$ana_dir,
    #
    )

 or die("Error in command line arguments.\n");


if(!(-e $contig_file)){
    die "CONTIG FILE NOT FOUND \n";
}	

#if(!(-e $scaffold_file)){
#    die "SCAFFOLD FILE NOT FOUND \n";
#}

if(!(-e $read_file)){
    die "READ FILE NOT FOUND \n";
}

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

print STDERR " *** Starting the sequence extraction\n";#<STDIN>;
#run_exe("rm -r $ana_dir");
run_exe("rm -rf $ana_dir/LOG");
run_exe("mkdir -p $ana_dir/LOG");
#my $run_racon_log = "$ana_dir/LOG/run_racon_log.txt";

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
my $contig_len_file = "$opera_lr_dir/contigs";
$scaffold_file = "$opera_lr_dir/scaffolds.scaf";
read_scaffold_file($scaffold_file);

print STDERR " *** Number of edges selected for gapfilling $nb_selected_edges\n";#<STDIN>;

#Get the read and the coordinates of the gaps

my ($reverse_edge, $reverse_edge_r);
my %contig_in_gap_to_extract = ();#list of contig for which the full sequence need to be stored
my %contig_in_gap_support = ();#To compute the support of the contig in the gap, the contig should be at least supported by 2 reads

my %contig_lengths = ();	# Temporarily holds the contig lengths so as to compute the coordinates for cutting the reads
#
#my @tmp_path = split(/long-read-mapping/,$edge_file_info);

#

open(LEN_FILE, $contig_len_file);
<LEN_FILE>;
while(<LEN_FILE>){
	@fields = split(/\t/, $_);
	$contig_lengths{$fields[0]} = $fields[1];
}
close(LEN_FILE);

open(FILE, $edge_file_info);
print STDERR " $ana_dir/contig_extention.log\n";
open(OUT_EXTRA, ">$ana_dir/contig_extention.log");
my ($start_contig_2_extra, $end_contig_1_extra);
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
    #@tmp = split(/\:/, $line[9]);@contig_in_gap = split(/\;/, $tmp[1]);
    #%contig_in_gap_support = ();
    
    #@coord_c1_on_r1 = split(/\_/, $coord_contig_1_on_read[0]);
    #@coord_c2_on_r1 = split(/\_/, $coord_contig_2_on_read[0]);
    @coord_c1_for_r1 = split(/\_/, $coord_contig_1[0]);
    @coord_c2_for_r1 = split(/\_/, $coord_contig_2[0]);

    #########################################################################################################################################################
    #Calculate the maximum possible overhang based on read 1
    #$max_seq_overhang_c1 = abs($coord_c1_on_r1[2] - $coord_c1_on_r1[1]);$max_seq_overhang_c1 = $sequence_overhang if($max_seq_overhang_c1 > $sequence_overhang);
    #$max_seq_overhang_c2 = abs($coord_c2_on_r1[2] - $coord_c2_on_r1[1]);$max_seq_overhang_c2 = $sequence_overhang if($max_seq_overhang_c2 > $sequence_overhang);
    #$max_seq_overhang_c1 = 0;
    #$max_seq_overhang_c2 = 0;
    #$edge_info{$edge_ID}->{"CONTIG_OVERHANG_INFO"}->{$contig_1} = $max_seq_overhang_c1;
    #$edge_info{$edge_ID}->{"CONTIG_OVERHANG_INFO"}->{$contig_2} = $max_seq_overhang_c2;
    #####################################################################################################################################################

    
    # Stores the reference read name (1st read)
    #$edge_info{$edge_ID}->{"REF_READ"} = $read_name[0];
    $reference_reads{$edge_ID} = $read_name[0];
        
    for(my $i = 0; $i < @read_name; $i++){
	$r = $read_name[$i];
	$reverse_edge_r = $reverse_edge;
	
	#Get the contig orientation to known if the extracted sequences have to be reversed
	@coord_c1 = split(/\_/, $coord_contig_1[$i]);
	@coord_c2 = split(/\_/, $coord_contig_2[$i]);
	#Check the orientation of the contig on the reads
	$ori_c1 = ($coord_c1[0] == 0 ? "+" : "-"); $ori_c2 = ($coord_c2[0] == 0 ? "+" : "-");
	$reverse_edge_r = ($reverse_edge_r + 1) % 2 if($ori_c1 ne $contig_1_ori && $ori_c2 ne $contig_2_ori);
	
	
	#This read is used by this edge
	if(! exists $read_to_edge{$r}){$read_to_edge{$r}={};}
	#Check for the order of the contigs and compute the coordinate of the read sequences that will be extracted
	@coord_c1_on_read = split(/\_/, $coord_contig_1_on_read[$i]);
	@coord_c2_on_read = split(/\_/, $coord_contig_2_on_read[$i]);


	############################################################
	#PROBLEM WITH READ THAT SUPPORT 2 DISTANCES OF AN EDGE
	$str_multi = "";
	if(exists $read_to_edge{$r}->{$edge_ID}){
	    #This read support twice the distance between the 2 contigs
	    $str_multi = "[MULTI_ID$cmp_multi_edge]";
	    $cmp_multi_edge++;
	}
	$str_multi = "";
	###########################################################
	
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
	    

	    print OUT_EXTRA "XXX" . "\t" . $r . "\t" . $edge_ID . "\t" . $coord_c1_on_read[2] . "\t" . $end_contig_1_extra . "\t" . $coord_c2_on_read[1] . "\t" . $start_contig_2_extra . "\n";
	    #print OUT_EXTRA $r . "\t" . $edge_ID . "\t" . $end_contig_1_extra . "\t" . $start_contig_2_extra . "\n";
		
	    $read_to_edge{$r}->{$edge_ID.$str_multi} = [$coord_c1_on_read[2] + $end_contig_1_extra, $coord_c2_on_read[1] - $start_contig_2_extra, $reverse_edge_r];
	    
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

	    print OUT_EXTRA  "YYY" . "\t" . $r . "\t" . $edge_ID . "\t" . $coord_c2_on_read[2] . "\t" . $end_contig_2_extra . "\t" . $coord_c1_on_read[1] . "\t" . $start_contig_1_extra . "\n";
	    #print OUT_EXTRA $r . "\t" . $edge_ID . "\t" . $end_contig_2_extra . "\t" . $start_contig_1_extra . "\n";
	    
	    $read_to_edge{$r}->{$edge_ID.$str_multi} = [$coord_c2_on_read[2] + $end_contig_2_extra, $coord_c1_on_read[1] - $start_contig_1_extra, $reverse_edge_r];
	    
	}
    }
    
}
close(FILE);
close(OUT_EXTRA);
#exit(0);

#Extract the sequences requires to fill the gaps in the reads
open(OUT_STATS, ">$ana_dir/gap_stats.dat") or die;
print STDERR " *** Extract gap sequences from read file $read_file\n";
open(READ_FILE, $read_file) or die;
my $seq;
my $nb_read = 0;my $nb_copy_of_contig_sequence = 0;
my $read_seq = "";
my $quality_seq = "";
my $keep_read = 0;

my %pure_gap_seq = ();
my %wrote_on_scaff = ();
##?## Debug purposes only
open(ERR_F, ">$ana_dir/edge_err");
#open(ERR_F2, ">$ana_dir/edge_err_nf");
my $wrote_on_scaff = ();

#Identify the scaffold to fill
my %scaff_to_fill = ();
#Only used in the case of gap size threshold
#TO_DO CHECK a single time which scacffold have to be filled
opendir(DIR, "$ana_dir/");
my @all_scaff_dir = readdir(DIR);
foreach $s_dir (@all_scaff_dir){ 
    $scaff_to_fill{$s_dir} = 1;
}
################

while(<READ_FILE>){
    
    chop $_;
    
    if(/^>(.*)/ or $is_fastq){
	
	if($read_seq ne  ""){
	    #print STDERR " *** read name: |$read_name|\n";<STDIN>
	    
	    $edge_list = $read_to_edge{$read_name};
	    if(!(defined $edge_list)){
		#print STDERR "$read_name NOT DEFINED POST 10000\n" if ($nb_read > 10000);
	    }

	    else{
		#$read_seq = <FILE>;
		#
		#next if(! exists $edge_info{$edge_ID}->{"SCAFF"});
		#
		
		%wrote_on_scaff = ();
		foreach $edge_ID (keys %{$edge_list}){
		    
		    #Write read file on the scaffold directory
		    $scaff_name = $edge_info{$edge_ID}->{"SCAFF"};

		    #Only used in the case of gap size threshold
		    #TO_DO CHECK a single time which scacffold have to be filled
		    next if(! exists $scaff_to_fill{$scaff_name});
		    		    
		    #print STDERR " *** Process edge $edge_ID $scaff_name\n";
		    
		    if(! exists $wrote_on_scaff{$scaff_name}){
			#print STDERR " **** Write $read_name to $ana_dir/$scaff_name/POOL.fastq\n";#<STDIN>;
			open(OUT_R, ">>$ana_dir/$scaff_name/POOL.fastq");
			print OUT_R 
			    "\@$read_name" . "\n" .
			    $read_seq . "\n" . 
			    "+" . "\n" . 
			    $quality_seq . "\n";
			close(OUT_R);
			$wrote_on_scaff{$scaff_name} = 1;
		    }
		    
		    #Only keep the first read
		    if($reference_reads{$edge_ID} eq $read_name){# && $edge_info{$edge_ID}->{"SEQ"} eq "NA"){

			$start = $edge_list->{$edge_ID}->[0];
			$end = $edge_list->{$edge_ID}->[1];
			$reverse_complement = $edge_list->{$edge_ID}->[2];
			#
			#Get the read sequence on the rightcomplement
			
			if($end<$start){
			    $pure_gap_seq{$edge_ID}->{$read_name} = $start - $end;
			    $seq = "";
			    #print ERR_F "$edge_ID\t$read_name\t$start\t$end\t".$edge_info{$edge_ID}->{"REF_READ"}."\t$pure_gap_seq{$edge_ID}->{$read_name}\n";
			    print ERR_F "$edge_ID\t$read_name\tNA\t\t$start\t$end\t".($start - $end)."\t"."OVERLAP\n";
			    $end = 1;
			    $start = 1;
			}
			else{
			    $seq = substr $read_seq, $start, $end-$start;
			    $seq = reverse_complement($seq) if($reverse_complement);
			    print ERR_F "$edge_ID\t$read_name\t$reverse_complement\t$start\t$end\t".($end-$start)."\tGAP\n";
			}
			$edge_info{$edge_ID}->{"SEQ"} = $seq;
			
		    }
		}
		delete $read_to_edge{$read_name}; 
	    }
	}

	#If the file is in fastq format, we assume that the sequence is on one line.
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
	    #$quality_seq .= "I" x length($read_seq);
	    #ADD the quality value
        }

    }
}
close(ERR_F);
close(READ_FILE);
close(OUT_STATS);


#to produce the filled gap scaffolds
write_filled_gap_scaffold();


sub read_scaffold_file{
    my ($scaffold_file) = @_;
    open(FILE, $scaffold_file);
    my ($contig_1, $contig_1_ori, $contig_2, $contig_2_ori, $scaff_name);
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	if($line[0] =~ m />(.+)/){
	    #if(index($line[0], ">") != -1){#this a new scaffold
	    $scaff_name = $1;
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
	    $edge_info{$edge_ID} = {"SCAFF", $scaff_name, "GAP_LENGTH", $gap_length, "SEQ", "NA", "EDGE_ORI", $edge_ori};#, "REF_READ", ""};# "REF_COORDS", {}};#be carreful with conflicting edges that have the same contig pair
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
    
    #Get all the contig sequence $contig_to_edge should have been already computed
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
    open(OUT_S, ">$ana_dir/gap_size.dat");
    open(FILE, $scaffold_file);
    my ($contig_1, $contig_1_ori, $contig_2, $contig_2_ori, $scaff_name, $gap_seq);
    $scaff_name = "";$contig_2 = "";$contig_1 = "";
    my $scaff_seq = "";
    my $gap_filled_stats = "";
    my $contig_trimming_length = 0;#This is used to remove starting sequences of a contig when 2 contigs are overlapping
    print STDERR " *** Read the scaffold file and fill gaps\n";
    
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	$gap_seq = "";
	if($line[0] =~ m />(.+)/){
	    
	    #Write the sequence of the last contig, sngleton scaffolds for which only contig_1 is not empty are skipted
	    if( exists $scaff_to_fill{$scaff_name} && $contig_2 ne ""){
		open(OUT, ">$ana_dir/$scaff_name/pre_consensus.fa");
		print OUT ">$scaff_name\n";
		print OUT $scaff_seq."".(get_contig_seq($contig_2, $contig_2_ori, \%all_contig_seq))."\n";
		close(OUT);
	    }
	    
	    #Beginning of a new scaffold
	    $scaff_seq = "";
	    $scaff_name = $1;
	    $contig_trimming_length = 0;
	    
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
	    next if(! exists $scaff_to_fill{$scaff_name});
	    $contig_2 = $line[0];
	    $contig_2_ori = "+"; $contig_2_ori = "-" if($line[1] eq "EB");
	    @contig_order = sort($contig_1, $contig_2);
	    
	    #Write contig 1
	    $contig_1_seq = get_contig_seq($contig_1, $contig_1_ori, \%all_contig_seq);
	    #Trim the contig in case of contig overlap (always trim on the left side)
	    if($contig_trimming_length != 0){
		#print STDERR " ***** Trimming $contig_1 $contig_trimming_length\n";
		#CHECK WHY $contig_trimming_length > $contig_1_seq
		if($contig_trimming_length < length($contig_1_seq)){
		    $contig_1_seq = substr $contig_1_seq, $contig_trimming_length;
		}
	    }
	    $scaff_seq .= $contig_1_seq;
	    
	    #print STDERR " *** Fill gap between $contig_1 and $contig_2\n";#<STDIN>;
	    #print STDERR " *** $contig_1 $contig_1_ori\n$contig_1_seq\n";<STDIN>;

	    #Write the gap sequence
	    $edge_ID = join("-vs-", @contig_order);
	    $gap_seq = extract_gap_seq($edge_ID);

	    #if($gap_seq eq ""){#This can happen when edges are derived from Illumina data
	    #print STDERR "**** REF NOT FOUND ($edge_ID) *****\n\n";
	    #}
	    
	    $gap_filled_stats = "NA";#the gap is not filled
	    my $gap_diff_large = "";
	    if($gap_seq ne ""){
		#this is an overlap or a no gap
		
		if ( $gap_seq =~ /^[0-9,.E]+$/ ) { 
		    $contig_trimming_length = $gap_seq;
		    $gap_filled_stats = -($contig_trimming_length);
		    #print STDERR " *** Overlap: $contig_trimming_length\n";#<STDIN>;
		    #$contig_trimming_length = 0;
		}
		else{
		    #This a gap sequence
		    $gap_filled_stats = length($gap_seq);
		    $contig_trimming_length = 0;
		    
		    #print STDERR " *** Fill gap seq\n$gap_seq\n";#<STDIN>;
		    $scaff_seq .= $gap_seq;
		}
		
	    }
	    #For some reason the gap was not filled
	    else{
		$contig_trimming_length = 0;
		if($gap_length > 0){$scaff_seq .= ( "N" x $gap_length);}
		else{$scaff_seq .= ( "N" x 3);}
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
    #if($contig_2 ne ""){
    if(exists $scaff_to_fill{$scaff_name} && $contig_2 ne ""){
	open(OUT, ">$ana_dir/$scaff_name/pre_consensus.fa");
	print OUT ">$scaff_name\n";
	print OUT $scaff_seq.(get_contig_seq($contig_2, $contig_2_ori, \%all_contig_seq))."\n";
	close(OUT);
    }
    #}
    
    close(FILE);
    close(OUT_S);
    
}

sub extract_gap_seq{
    my ($edge_ID) = @_;
    
    $seq = $edge_info{$edge_ID}->{"SEQ"};
    
    if($seq ne ""){
    	return $seq;
    }
    else{
	my $out = $pure_gap_seq{$edge_ID}->{$reference_reads{$edge_ID}};
	
	#print STDERR $edge_ID . "\t".$reference_reads{$edge_ID}."\t".$pure_gap_seq{$edge_ID}->{$reference_reads{$edge_ID}}."\n";
	
    	if(! (defined $out) ) {
	    print STDERR "Size entry absent from long-read overlap " . $edge_ID . "\t" . $reference_reads{$edge_ID}."\n";
	    $out = "";
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


sub reverse_complement {
    my ($seq) = @_;
    
    # reverse the DNA sequence
    my $revcomp = reverse($seq);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $revcomp;
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
