#!/usr/bin/perl

#use strict;
use warnings;

my ($contig_file, $scaffold_dir, $edge_read_info_file, $seq_read_file, $opera_bin_dir, $racon_dir, $minimap2_dir, $mummer_dir) = @ARGV;

my $re_run = 1;

my @tmp = split(/\//, $scaffold_dir);
my $scaffold_id = @tmp[@tmp-1];

#print STDERR " *** $scaffold_id\n";<STDIN>;

my $assembly_dir =  "$scaffold_dir";
run_exe("mkdir $assembly_dir");
    
#Run racon
#Get the initial assembly
my $initial_assembly = "$scaffold_dir/test_gapfil.fa";
run_exe("perl $opera_bin_dir/extract_read_sequence_xargs.pl --edge-file $edge_read_info_file --contig-file $contig_file --scaffold-file $scaffold_dir/scaffolds.scaf --output-directory  $scaffold_dir/GAP_FILLING --num-of-processors 1 --script $opera_bin_dir --read-file $seq_read_file --outfile $initial_assembly");

#Consensus construction
my $mapping_file = "$assembly_dir/test_gapfil.fa.filtered.paf";
my $read_file = "$assembly_dir/POOL.fastq";
my $racon_consensus = "$assembly_dir/racon.fa";
if($re_run || ! -e $racon_consensus){
    #run_exe("mkdir $assembly_dir/../CORR/;touch $assembly_dir/../CORR/single_read.map");
    get_read_mapping($scaffold_dir, $edge_read_info_file, "$assembly_dir/single_read.map", $seq_read_file, $assembly_dir);
    run_racon($initial_assembly, $read_file, $mapping_file, $racon_consensus, 500, 5, -4, -8, -6, "CONSENSUS");
}

#Check if racon consensus empty due to low read coverage or edge identified by short reads
my $concencus_seq_size = `wc -l $assembly_dir/racon.fa | awk '{print \$1}'`;chop $concencus_seq_size;
if($concencus_seq_size == 0){### 
    #There is not enought coverage to produce a concensus for this gap
    #The pre-consensus will be used instead, it cam be improve througth tiling, even if the gap may be of low quality
    run_exe("cp $initial_assembly $racon_consensus");
}

#Tile the contigs
my $contig_to_tile = $contig_file;

#Idnetification of the corrected read
$prefix_name = "$assembly_dir/consensus";
my $contig_tile_mapping = "$prefix_name\_tiling.paf";
print STDERR " *** Identify corrected reads\n";
#Used minimap2 to insert the contigs
run_nucmmer($contig_to_tile, $racon_consensus, $prefix_name);
#
get_paf_file("$prefix_name\_tiling.coords", $contig_tile_mapping, "ONLY_SCAFF_CONTIG");

#<STDIN>;

#Get the contig to tile sequence and the corrected reads
my $corrected_read = "$assembly_dir/corrected_read";
my $contig_in_tiling = "$assembly_dir/contig_in_tiling";
get_contig_and_corrected_read_in_gap($contig_tile_mapping, $racon_consensus, $contig_to_tile, "$corrected_read.fa", "$contig_in_tiling.fa");

#Intergrate the contigs in concensus
run_exe("$opera_bin_dir/Remap.py $contig_tile_mapping $racon_consensus $contig_in_tiling.fa");
$ncr_assembly = "$assembly_dir/racon_remapped.fasta";
#run_exe("mv $scaffold_dir/test_gapfil_remapped.fasta $ncr_assembly");


#Identify the final set of contigs used in the tiling
$prefix_name = "$assembly_dir/final_contig_set";
$contig_tile_mapping = "$prefix_name\_tiling.paf";
$updated_corrected_read = "$assembly_dir/updated_corrected_read";
$contig_in_tiling = "$assembly_dir/final_contig_in_tiling";
run_nucmmer($contig_to_tile, $ncr_assembly, $prefix_name);
get_paf_file("$prefix_name\_tiling.coords", $contig_tile_mapping, "ALL");
get_contig_and_corrected_read_in_gap($contig_tile_mapping, $ncr_assembly, $contig_to_tile, "$updated_corrected_read.fa", "$contig_in_tiling.fa");

#Tile the contigs based on mapping to reference
run_exe("$opera_bin_dir/Remap.py $contig_tile_mapping $ncr_assembly $contig_in_tiling.fa");

#Delete the intermidate file
#run_exe("rm $assembly_dir/TILING_POOL.fastq $assembly_dir/POOL.fastq");
#run_exe("rm $assembly_dir/TILING_POOL.fastq");

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


sub get_contig_and_corrected_read_in_gap{

    my ($mapping_file, $seq_file, $contig_file, $out_read_file, $out_contig_file) = @_;


    #print STDERR " *** get_contig_and_corrected_read_in_gap $mapping_file  $seq_file  $contig_file  $out_read_file  $out_contig_file\n";<STDIN>;
    
    my $OVERHANG_SEQ = 200;
    my @cut_coordinates = ();

    #Get the list of contig in the mapping file
    open(FILE, "$mapping_file");
    my %contig_mapped = ();
    while(<FILE>){
	@line = split(/\t/, $_);
	$contig_name = $line[0];
	$contig_mapped{">$contig_name"} = 1;
	$nb_contig_to_print++;

    }

    #get the sequence of the contigs in fasta format
    open(IN, $contig_file);
    open(OUT, ">$out_contig_file"); 
    my $print_flag = 0;
    while(<IN>){
	if(index($_, ">") == -1){
	    if($print_flag){
		print OUT $_;
		$nb_contig_to_print--;
		last if($nb_contig_to_print == 0);
	    }
	}
	else{
	    #print STDERR "---> $nb_contig\n" if($nb_contig % 500000 == 0);$nb_contig++;
	    my @line = split(/\s+/, $_);
	    my $ID = substr ($line[0], 0, length($line[0]));
	    #print STDERR $ID."\n"; <STDIN>;
	    if(exists $contig_mapped{$ID}){
		$print_flag = 1;
		#chop $_;
		#print $_."\n";
		print OUT $ID . "\n";
	    }
	    else{
		#print STDERR "**************** IN\n";
		#if($print_flag == 1);
		$print_flag = 0;
	    }
	}
    }
    close(IN);
    close(OUT);
    
    #Map the contig sequence to the corrected assembly
    my $tiling_mapping = "$seq_file.paf";
    run_exe("$minimap2_dir/minimap2 -t1 --cs=short -w5 -m0 $seq_file $out_contig_file > $tiling_mapping");#<STDIN>;

    #Identify the contigs gap sequences in between the contigs or the overlaps
    open(FILE, "sort -k8,8 -n $tiling_mapping | ");
    my $prev_contig_end = -1;
    my $current_contig_start = -1;
    %contig_mapped = ();
    my $nb_contig_to_print = 0;
    my $region_in_contig = 0;
    my $region_in_gap = 0;
    while(<FILE>){
	@line = split(/\t/, $_);

	$contig_name = $line[0];
	$contig_mapped{">$contig_name"} = 1;
	
	$current_contig_start = $line[7];
	$nb_contig_to_print++;

	#For the contig sequence
	$current_contig_end = $line[8];
	$region_in_contig +=  ($current_contig_end - $current_contig_start);
	
	if($prev_contig_end != -1){
	    $gap_s = $prev_contig_end;
	    $gap_e = $current_contig_start;
	    if($gap_e < $gap_s){
		$gap_s = $current_contig_start;
		$gap_e = $prev_contig_end;
	    }
	    else{
		#This is a gap
		$region_in_gap +=  ($gap_e - $gap_s);
	    }
	    print STDERR " *** GAP $gap_s $gap_e \n";
	    
	    push(@cut_coordinates, $gap_s - $OVERHANG_SEQ);
	    push(@cut_coordinates, $gap_e + $OVERHANG_SEQ);
	}

	$prev_contig_end = $line[8];
    }
    
    open(OUT, ">$out_read_file");
    open(FILE, $seq_file);
    <FILE>;
    my $seq = <FILE>;
    my ($start, $end, $prev_start, $prev_end);
    $prev_end = -1; $prev_start = -1;
    my $gap_id = 1;
    #Get the splicing
    for(my $i = 0; $i < @cut_coordinates; $i+=2){

	$start = $cut_coordinates[$i];
	$end = $cut_coordinates[$i+1];

	print STDERR $start . "\t" . $end ."\n";
	
	if($prev_end > $start){
	    #print STDERR " *** WARNNING OVERLAPPING SPLICING $prev_start $prev_end  $start $end\n";
	    #Extend the read sequence to the end
	    $prev_end = $end;
	}
	else{
	    if($prev_start != -1){
		$s = substr $seq, $prev_start, ($prev_end-$prev_start);
		print OUT ">gap_$gap_id\_$prev_start\_$prev_end\n";
		print OUT $s ."\n";
		$gap_id++;
	    }
	    $prev_start = $start;
	    $prev_end = $end;
	}
	
    }
    #Write the last gap
    $s = substr $seq, $prev_start, ($prev_end-$prev_start);
    print OUT ">gap_$gap_id\_$prev_start\_$prev_end\n";
    print OUT $s ."\n";
    close(OUT);

    
    
    #exit(0);
}

#Should be done once to avoid reading the read file multiple times
sub get_read_mapping{
    my ($scaffold_dir, $edge_read_info_file, $single_contig_on_file, $long_read_file, $out_dir) = @_;

    $fastq_pool = "$out_dir/POOL.fastq";
    if($re_run || ! -e $fastq_pool){
	run_exe("python $opera_bin_dir/Filter_read_names_edge.py $edge_read_info_file $single_contig_on_file $scaffold_dir/scaffolds.scaf $long_read_file > $fastq_pool");
    }

    $mapping_file = "$out_dir/test_gapfil.fa.filtered.paf";
    run_exe("$minimap2_dir/minimap2 --cs=short -w5 -m0 $scaffold_dir/test_gapfil.fa $fastq_pool > $mapping_file");

    if(0){
	#Get the read file
	#$fastq_pool = "$current_dir/SECOND_POOL.fastq";
	#$read_line = `nl $scaffold_dir/UPDATED_POOL.fastq | grep  "\@k99_" | head -n1 | awk '{print \$1}'`;chop $read_line;$read_line--;
	#replace_quality_value("head -n $read_line $scaffold_dir/UPDATED_POOL.fastq |", $fastq_pool, "-");
	
	#In this case we are simply filtering the read mapping from the previous run
	my $mapping_file = "$scaffold_dir/test_gapfil.fa.paf";
	my $filtered_paf = "$out_dir/test_gapfil.fa.filtered.paf";
	
	open(FILE, "sort -k8,8 -n $mapping_file | cut -f1-10 | grep k99_ |");
	my @contig_coord = ();
	while(<FILE>){
	    @line = split(/\t/, $_);
	    push(@contig_coord, $line[7]);
	    push(@contig_coord, $line[8]);
	    #print OUT $_;
	}
	close(FILE);

	#print STDERR " *** $filtered_paf\n";<STDIN>;
	open(OUT, ">$filtered_paf");
	open(FILE, "sort -k8,8 -n $mapping_file | grep -v k99_ |");
	my $OVERHANG_SEQ_THRESHOLD = 200;#Read must at least be out of the gap from 200bp => if overhang is to low read filetered
	my $UNMAPPED_SEQ_THRESHOLD = 200;#Fraction of the read that may be unmapped, if read not completly unmapped filtered
	
	my $current_contig_start = $contig_coord[0];
	my $current_contig_end = $contig_coord[1];
	my $cmp_contig_coord = 1;
	while(<FILE>){
	    @line = split(/\t/, $_);
	    #
	    $read_length = $line[1];
	    $read_start_c = $line[2];
	    $read_end_c = $line[3];
	    
	    $read_start = $line[7];
	    $read_end = $line[8];
	    

	    if($read_end_c - $read_start_c < $read_length - ($UNMAPPED_SEQ_THRESHOLD * 2)){
		#print STDERR " *** Read not fully mapped $_ $read_length $read_start $read_end =>  " . ($read_end - $read_start) ."\n";
		next;
	    }
	    
	    #Need to check for the next contig
	    if($read_start > $current_contig_end){
		$current_contig_start = $contig_coord[$cmp_contig_coord];
		$current_contig_end = $contig_coord[$cmp_contig_coord];
		$cmp_contig_coord += 2;
	    }
	    
	    #intersection out of the current contig
	    if($read_end > $current_contig_end + $OVERHANG_SEQ_THRESHOLD){#Read with enought sequence out of the contig
		#print STDERR $_;
		print OUT $_;
	    }
	    
	}
	close(FILE);
	close(OUT);
    }
    #exit(0);
}

sub run_nucmmer{
    my ($contig_to_tile, $assembly_file, $nucmer_file) = @_;
    #print STDERR " *** run_nucmmer $contig_to_tile, $assembly_file, $nucmer_file\n";<STDIN>;
    if($re_run || ! -e "$nucmer_file\_tiling.coords"){
	run_exe("$mummer_dir/nucmer --maxmatch -p $nucmer_file $assembly_file $contig_to_tile");
	run_exe("$mummer_dir/delta-filter -r $nucmer_file.delta > $nucmer_file\_tiling.delta");
	run_exe("$mummer_dir/show-coords -T -l $nucmer_file\_tiling.delta | sort -k1,1 -n  > $nucmer_file\_tiling.coords");
    }
}

sub get_paf_file{
    my ($coord_file, $paf_file, $flag_set_type) = @_;

    #This is the set of contigs that should belong to set of contigs used for the tilling
    my %contig_set_to_keep = ();
    #Keep all the contitgs that belong to the current scaffolds
    #No identity threshold is applied to those contigs
    open(FILE, "$scaffold_dir/scaffolds.scaf");<FILE>;
    while(<FILE>){
	@line = split(/\t/, $_);
	$contig_set_to_keep{$line[0]} = 1;
    }
    close(FILE);

    #Add the contigs for the same species
    if($flag_set_type eq "ALL"){
	get_contig_in_species(\%contig_set_to_keep);
    }

    
    open(OUT, ">$paf_file");
    open(FILE, "$coord_file");
    my $IDENTITY_THRESHOLD = 0.97;
    my $COVERAGE_THRESHOLD = 0.90;
    my $previous_contig_name = "";
    my $prev_scaffold_end = -1;
    my ($contig_name, $scaffold_name, $scaffold_length, $mapping_length, $contig_end, $contig_start, $scaffold_start, $scaffold_end, $tiling_length);
    open(OUT_G, ">S1_gap_size.dat");
    <FILE>;<FILE>;<FILE>;<FILE>;#skip the header lines
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	$contig_name = $line[10];
	$scaffold_name = $line[9];
	$scaffold_length = $line[7];
	$identity = $line[6]/100;
	
	next if(
	    ($flag_set_type eq "ONLY_SCAFF_CONTIG" && !(exists $contig_set_to_keep{$contig_name})) &&
	    ($flag_set_type eq "ALL" && !(exists $contig_set_to_keep{$contig_name})) &&
	     $identity < $IDENTITY_THRESHOLD);
	
	if($contig_name ne $previous_contig_name){
	    if($previous_contig_name ne ""){
		
		$mapping_length = ($contig_end - $contig_start);
		#scaffold length
		#print " *** MAP " . $previous_contig_name . " " . $mapping_length . " " . $contig_end . " " .$contig_start . "\n";
		if(($flag_set_type eq  "ONLY_SCAFF_CONTIG" && (exists $contig_set_to_keep{$contig_name}) && $mapping_length / $contig_len > 0.55) || #Always keep contigs from the contig to keep even if partial alignements
		   $mapping_length / $contig_len > $COVERAGE_THRESHOLD){
		    $tiled_contig{$previous_contig_name} = 1;
		    
		    print OUT
			$previous_contig_name . "\t" . $contig_len . "\t" . $contig_start . "\t" . $contig_end . "\t" . $contig_orientation . "\t" . 
			$scaffold_name . "\t" . $scaffold_length . "\t" .
			$scaffold_start . "\t" . $scaffold_end . "\t" .
			int(($scaffold_end - $scaffold_start) * $identity) . "\t" . ($scaffold_end - $scaffold_start) . "\t" .
			"60" . "\t" . 
			"tp:A:P" . "\t". "cm:i:175" . "\t" . "s1:i:915"	. "\t" . "s2:i:0" . 
			"\n";

		    #Sum contig length
		    $sum_contig_length += $contig_len;
		    
		    #Tilling length
		    $tiling_length += $mapping_length;

		    if($prev_scaffold_end > 0 && $prev_scaffold_end < $scaffold_start){
			$gap_length += ($scaffold_start - $prev_scaffold_end);
			print  OUT_G ($scaffold_start - $prev_scaffold_end) . "\t" . $previous_contig_name . "\t" .$contig_name . "\n";
		    }
		}
		#
		$prev_scaffold_end = $scaffold_end;
	    }
	    #$contig_len = 0;
	    #$contig_start = 0;
	    #$contig_end = 0;
	    #$scaffold_start = 0;
	    #$scaffold_end = 0;
	    $previous_contig_name = $contig_name;
	    $contig_len = $line[8];
	    $contig_start = $line[2];
	    $contig_end = $line[3];
	    #$identity = $line[6]/100 * ($contig_end -$contig_start);
	    #Get the contig orientation and reverse coordinate in case of -
	    $contig_orientation = "+";
	    if($contig_end < $contig_start){
		$contig_orientation = "-";
		$tmp =  $contig_end; $contig_end = $contig_start;$contig_start = $tmp;
	    }
	    $contig_start--;
	    #
	    $scaffold_start = $line[0] - 1;
	    $scaffold_end = $line[1];
	}
	else{
	    #The alignement can be extended
	    #NEED TO BE CAREFUL IN CASE OF TANDEM REPEAT !!!!

	    #Check if this a correct extention
	    if($contig_orientation eq "+"){
		if($contig_end < $line[3]){#this is an extention
		    $contig_end = $line[3];#the end is extended with next end
		    $scaffold_end = $line[1];
		}
	    }
	    else{
		if($contig_start > $line[3]-1){#this is an extention
		    $contig_start = $line[3] -1;#the end is extended with next start
		    $scaffold_end = $line[1];
		}
	    }
	    #$identity = $line[6]/100 * ($contig_end -$contig_start);
	    
	}
	#$previous_contig_name = $contig_name;
    }
    close(FILE);
    
    #For the last contig
    $tiled_contig{$previous_contig_name} = 1;
    print OUT
	$previous_contig_name . "\t" . $contig_len . "\t" . $contig_start . "\t" . $contig_end . "\t" . $contig_orientation . "\t" .
	$scaffold_name . "\t" . $scaffold_length . "\t" .
	$scaffold_start . "\t" . $scaffold_end . "\t" .
	int(($scaffold_end - $scaffold_start) * 0.99) . "\t" . ($scaffold_end - $scaffold_start) . "\t" .
	"60" . "\t" . 
	"tp:A:P" . "\t". "cm:i:0" . "\t" . "s1:i:0"	. "\t" . "s2:i:0" . 
	"\n";
    
    close(OUT);
    $sum_contig_length += $contig_len;
    $tiling_length += ($contig_end - $contig_start);
    if($prev_scaffold_end > 0 && $prev_scaffold_end < $scaffold_start){
	$gap_length += ($scaffold_start - $prev_scaffold_end);
    }
    
    print STDERR "Gap length: $gap_length\n";
    print STDERR "Tiling length: $tiling_length \t $sum_contig_length\n";
    print STDERR "Fraction gap " . ( $gap_length / ($gap_length + $tiling_length) ) . "\n";
    close(OUT_G);
    #<STDIN>;
    #exit(0);
}


sub replace_quality_value{
    my ($fastq_file, $updated_fastq_file, $new_qual) = @_;
    open(OUT, ">$updated_fastq_file");
    open(FILE, $fastq_file);
    my $cmp_line = 1;
    while(<FILE>){
	if($cmp_line != 4){
	    print OUT $_;
	    $cmp_line++;
	}
	else{
	    if($new_qual ne "-"){
		substr($_,0, -1, $new_qual x (length($_)-1));
	    }
	    print OUT $_;
	    $cmp_line =1;
	}
    }
    close(FILE);
    close(OUT);
}
    
sub get_contig_in_species{
    my ($contig_to_keep) = @_;

    my $scaffold_info = "$scaffold_dir/../../../../scaffold_info.txt";
    my $cluster_file = "$scaffold_dir/../../../reference_mapping/clusters_seq_similarity";
    my @cluster_list = ();
    #Get the species of the scaffold
    my $genome_path = `grep -w $scaffold_id $scaffold_info | cut -f 6`;chop $genome_path;
    
    #No psecies identified for that cluster
    if($genome_path ne "NO SPECIES"){
	@tmp = split(/\//, $genome_path);
	$strain = @tmp[@tmp-1];
	@tmp = split(/_/, $strain);
	$species = $tmp[0] . "_" . $tmp[1];
	
	$all_cluster = `grep  $species $scaffold_info | cut -f 4`;
	my @cluster_list = split(/\n/, $all_cluster);
	my %cluster_to_keep = ();
	foreach $c (@cluster_list){
	    $cluster_to_keep{$c} = 1;
	}
	
	#Get the contig from the cluster set
	#print $out_dir . "\n";
	open(FILE, $cluster_file);
	while(<FILE>){
	    @line = split(/\t/, $_);
	    $cluster = $line[1];
	    $contig = $line[0];
	    if(exists $cluster_to_keep{$cluster}){
		$contig_to_keep->{$contig} = 1;
	    }
	}
	close(FILE);
    }
}


sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return . "\n" if($run);
    return $return;
}

#canu -p canu -d  scaffolds/opera_scaffold_69/EDGE_READ/CANU  genomeSize=120000 useGrid=0 maxThreads=1 minThreads=1 stopOnReadQuality=false scaffolds/opera_scaffold_69/EDGE_READ/
