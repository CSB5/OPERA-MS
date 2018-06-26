#!/usr/bin/perl
use warnings;
use Getopt::Long;

my $min_contig_size = 300;
my $min_opera_contig_size = 500;
my $long_read_cluster_threshold = 2;
for($i = 0; $i < 6; $i++){$cluster_threshold_tab[$i] = $long_read_cluster_threshold;}
my $short_read_cluster_threshold = 3;
#my @cluster_threshold_tab = (2,2,2,2,1,1);
#my @cluster_threshold_tab = (2,2,2,2,2,1);
my $FLAG_FILTER_CONFLICTING_ALIGNEMENT = 0;
#my $FLAG_NO_CONFLICTING_EDGE = 1;
my $mapper = "blasr";

#Minimum fraction of contig mapped to be concider fully contained on the read
my $fraction = 0.9;
#Minimum contig alignemnt length
my $min_alignment_length = 500;

#Used for 2 different purpose:
#1) maximum overlap allowed between 2 mapped contigs
#2) non-ovelapping sequences allowed on the contig and read to be concidered as a valide alignement
my $overlap = 200;


#my $mapper_extention 
my $graphmapDir = "/mnt/software/stow/graphmap-0.3.0-1d16f07/bin/";


#Used in case of grapmap mapping
my @contig_id_to_name = ();
#list of repeat contigs
my %repeat_contig = ();

#Init the software directories variables
my $blasrDir = "";
my $operaDir = "";
my $short_read_tooldir = "";
my $samtools_dir = "";
my $short_read_maptool = "bwa";
my $kmer_size = 100;
my $flag_help;


my $help_message = "

Options:
        --contig-file: fasta file of contigs
        --kmer: size of the kmer used to produce the contig (default 100)
        --long-read-file: fasta file of long reads
        --output-prefix: prefix of output mapping file
        --output-directory: output directory for scaffolding results
        --num-of-processors: number of processors used for mapping stages
        --blasr: Folder which contains blasr binary (default PATH)
        --short-read-maptool: Mapping tool can be either bwa (default) or bowtie
	--short-read-tooldir: Directory that contains binaries to the chosen short read mapping tool (default PATH)
        --samtools-dir:	Directory that contains samtools binaries (default PATH)
	--opera: Folder which contains opera binary (default PATH)
	--illumina-read1: fasta file of illumina read1
	--illumina-read2: fasta file of illumina read2
        --help : prints this message
";

if ( @ARGV == 0 ) {
        print $help_message;
        exit 0;
}

GetOptions(
    "contig-file=s"    => \$contigFile,
    "long-read-file=s"    => \$readsFile,
    "output-prefix=s" => \$file_pref,
    "output-directory=s" => \$outputDir,
    "num-of-processors=i" => \$nproc,
    "kmer=i" => \$kmer_size,
    "blasr=s"      => \$blasrDir,
    "graphmap=s"      => \$graphmapDir,
    "opera=s"      => \$operaDir,
    "samtools-dir=s"  => \$samtools_dir,
    "short-read-maptool=s" => \$short_read_maptool,
    "short-read-tooldir=s" => \$short_read_tooldir,
    "illumina-read1=s"      => \$illum_read1,
    "illumina-read2=s"      => \$illum_read2,
    "help"       => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
        print $help_message;
        exit 0;
}
if (!defined($contigFile)) {
        print "contigs fasta file needs to be specified\n";
        exit 0;
}
if (!defined($readsFile)) {
        print "long reads fasta file needs to be specified\n";
        exit 0;
}
if (!defined($file_pref)) {
        print "prefix of output mapping file needs to be specified\n";
        exit 0;
}
if (!defined($outputDir)) {
        print "output directory for scaffolding results needs to be specified\n";
        exit 0;
}
#if (!defined($blasrDir)) {
#        print "Folder which contains blasr binary if it is not in PATH needs to be specified\n";
#        exit 0;
#}
#if (!defined($sspaceDir)) {
#        print "Folder which contains SSPACE perl script if it is not in PATH needs to be specified\n";
#        exit 0;
#}
#if (!defined($operaDir)) {
#        print "older which contains opera binary if it is not in PATH needs to be specified\n";
#        exit 0;
#}
if (!defined($illum_read1) || !defined($illum_read2)) {
    print STDERR " *** WARNING illumina fasta file not fully specified\n";
    $illum_read1 = "NONE";
    $illum_read2 = "NONE";
}


if( $outputDir !~ "/\$" && $outputDir ne "" )
{
    $outputDir .= "/";
}
run_exe("mkdir $outputDir") unless(-d $outputDir);
if( !-d $outputDir ){
    print STDERR "Error: the output directory does not exist, please try again.\n";
    exit( -1 );
}

if( !defined( $nproc ) ){
    $nproc = 1;
}
#To make that those varibles are really directories
if( $blasrDir !~ "/\$" && $blasrDir ne "" )
{
    $blasrDir .= "/";
}
if( $operaDir !~ "/\$" && $operaDir ne "" )
{
    $operaDir .= "/";
}
if( $samtools_dir !~ "/\$" && $samtools_dir ne "" )
{
    $samtools_dir .= "/";
}

chdir( $outputDir );

my $str_full_path = "or please enter the full path";
if ( ! -e $contigFile ) {die "\nError: $contigFile - contig file does not exist $str_full_path\n"};
if ( ! -e $readsFile ) {die "\nError: $readsFile - long read file does not exist $str_full_path\n"};
if ( ! -e $illum_read1 ) {die "\nError: $illum_read1 - illumina read 1 file does not exist $str_full_path\n"};
if ( ! -e $illum_read2 ) {die "\nError: $illum_read2 - illumina read 2 file does not exist $str_full_path\n"};

if ( ! -e "$blasrDir/blasr" && $blasrDir ne "") {die "\nError: $blasrDir - blasr does not exist in the directory $str_full_path\n"};
#if ( ! -e "$graphmapDir/graphmap" ) {die "$! graphmap does not exist in the directory $str_full_path\n"};
if ( ! -e "$operaDir/OPERA-LG" && $operaDir ne "") {die "\nError:$operaDir - OPERA-LG does not exist in the directory $str_full_path\n"};
if ( ! -e "$short_read_tooldir/bowtie" && $short_read_maptool eq "bowtie" && $short_read_tooldir ne "") {die "\nError: $short_read_tooldir - bowtie does not exist in the directory $str_full_path\n"};
if ( ! -e "$short_read_tooldir/bwa" && $short_read_maptool eq "bwa" && $short_read_tooldir ne "") {die "\nError: $short_read_tooldir - bwa does not exist in the directory $str_full_path\n"};


#map illumina reads to the contigs using preprocess_reads.pl
if(! -e "${file_pref}.bam" &&  !($illum_read1 eq "NONE" && $illum_read2 eq "NONE")){
    print "Mapping short-reads using  $short_read_maptool...\n";
    $str_path_dir = "";
    $str_path_dir .= "--tool-dir  $short_read_tooldir" if($short_read_tooldir ne "");
    $str_path_dir .= " --samtools-dir $samtools_dir" if($samtools_dir ne "");
    run_exe("perl $operaDir/preprocess_reads.pl --nproc $nproc --map-tool $short_read_maptool $str_path_dir --contig $contigFile --illumina-read1 $illum_read1 --illumina-read2 $illum_read2 --out ${file_pref}.bam");
}

if(! -e "$file_pref.map.sort"){
    # map using blasr
    print "Mapping long-reads using blasr...\n";
    #run_exe( "${blasrDir}blasr -nproc $nproc -m 1 $readsFile $contigFile  | cut -d ' ' -f1-4,7-13 | sed 's/ /\\t/g' > $file_pref.map");
    if($mapper eq "blasr"){
	run_exe( "${blasrDir}blasr  -nproc $nproc -m 1 -minMatch 5 -bestn 10 -noSplitSubreads -advanceExactMatches 1 -nCandidates 1 -maxAnchorsPerPosition 1 -sdpTupleSize 7 $readsFile $contigFile | cut -d ' ' -f1-5,7-12 | sed 's/ /\\t/g' > $file_pref.map");
	# sort mapping
	print "Sorting mapping results...\n";
	run_exe("sort -k1,1 -k9,9g  $file_pref.map > $file_pref.map.sort") 
    }
    
    if($mapper eq "graphmap"){
	print "Mapping using graphmap...\n";
	run_exe("$graphmapDir/graphmap owler -t 20 -r $contigFile -d $readsFile -o $file_pref.map");
	print "Sorting mapping results...\n";
	run_exe("sort -k1,1 -k6,6g  $file_pref.map > $file_pref.map.sort");
    }
}

#Read the contig file to get an array contig ID -> contig_name as the map file contain only read and contig identifier based on thei line number
if($mapper eq "graphmap"){
    print "Analyse contig file...\n";
    open(FILE, "grep \">\" $contigFile | sed 's/>//' |");
    while(<FILE>){
	@line = split(/\s+/, $_);
	push(@contig_id_to_name, $line[0]);
    }
    close(FILE);
}
# analyze mapping file
print "Analyzing sorted results...\n";
#my $all_edge_file = "pairedEdges";
my $all_edge_file = "pairedEdges";
&checkMapping( "$file_pref.map.sort", $all_edge_file);

# extract edges
print "Extracting linking information...\n";
#extract_edge("pairedEdges");
extract_edge($all_edge_file);

#Write the edge read information
#open(OUT, ">edge_read_info.dat");
#foreach $edge (keys %edge_read_info){
#    print OUT 
#	$edge_read_info{$edge}->{"EDGE"}."\t".
#	join(";", @{$edge_read_info{$edge}->{"READ_LIST"}})."\t".
#	"COORD_CONTIG_1:".join(";", @{$edge_read_info{$edge}->{"COORD_CONTIG_1"}})."\t".
#	"COORD_CONTIG_2:".join(";", @{$edge_read_info{$edge}->{"COORD_CONTIG_2"}})."\t".
#	"COORD_CONTIG_1_ON_READ:".join(";", @{$edge_read_info{$edge}->{"COORD_CONTIG_1_ON_READ"}})."\t".
#	"COORD_CONTIG_2_ON_READ:".join(";", @{$edge_read_info{$edge}->{"COORD_CONTIG_2_ON_READ"}})."\n";
#}
#close(OUT);

my @allEdgeFiles = ();

#Detect repeat and confiltint edges that are filtered out
print "Repeat detection...\n";

open(OUT, ">contig_length.dat");
foreach $c (keys %contig_length){
    print OUT $c."\t".$contig_length{$c}."\t"."20"."\n";
}
close(OUT);

open(OUT, ">repeat.dat");
for (my $i = 0; $i <= 0; $i++){
    $edge_file = $all_edge_file."_i$i";
    #run_exe("$operaDir/filter_conflicting_edge.pl results-short-reads/clustersInfo_opera-lr results-short-reads/contigs 100 $short_read_cluster_threshold");
    run_exe("$operaDir/filter_conflicting_edge.pl pairedEdges_i0 contig_length.dat 100 $cluster_threshold_tab[0]");
    
    #read the anchor_contig_info file and consider all contigs as repeat
    open(FILE, "anchor_contig_info.dat");
    while(<FILE>){
	@line = split(/\t/, $_);
	$contig = $line[0];
	#print STDERR " *** edge_ID $edge_ID\n";<STDIN>;
	if(! exists $repeat_contig{$contig}){
	    $repeat_contig{$contig} = 1;
	    print OUT $_;
	    #print STDERR " *** CONTIG $contig\n";
	}
    }
    close(FILE);
    close(OUT);
}
#delete some intermidate file
run_exe("rm anchor_contig_info.dat contig_length.dat filtered_edges.dat filtered_edges_cov.dat *.sai");

#Filter the bam file to remove repeat contigs
#need to change that and add the repeat module in the OPERA-LG code
if(! -e "$file_pref-with-repeat.bam"){
    run_exe("mv $file_pref.bam $file_pref-with-repeat.bam");
    run_exe("$operaDir/filter_repeat.pl $file_pref-with-repeat.bam repeat.dat $samtools_dir | ${samtools_dir}samtools view - -h -S -b > $file_pref.bam");
    run_exe("rm $file_pref-with-repeat.bam");
}

#Filter the long read edges
for (my $i = 0; $i <= 5; $i++){
    $edge_file = $all_edge_file."_i$i";
    $updated_edge_file = $all_edge_file."_no_repeat_i$i";
    push(@allEdgeFiles, $updated_edge_file);
    #read the anchor_contig_info file and consider all contigs as repeat
    open(FILE, "$edge_file");
    open(OUT, ">$updated_edge_file");
    while(<FILE>){
	@line = split(/\t/, $_);
	$contig1 = $line[0];
	$contig2 = $line[2];
	if(! (exists($repeat_contig{$contig1}) || exists($repeat_contig{$contig2}))){
	    print OUT $_;
	}
    }
    close(FILE);
    close(OUT);
}

#
# create configure file
&CreateConfigFile( $contigFile, "", @allEdgeFiles );

# run opera
&run_exe( "export GLIBCXX_FORCE_NEW; valgrind --tool=memcheck --leak-resolution=high --track-origins=yes --leak-check=full --error-limit=no --show-reachable=yes --num-callers=49 --suppressions=${operaDir}string.supp --track-fds=yes ${operaDir}OPERA-LG config > log 2> errLog" );

#Link to the result file
&run_exe("ln -s results/scaffoldSeq.fasta scaffoldSeq.fasta");

sub extract_edge{
    my ($all_edge_file) = @_;

    my %inter = (
	"i0", [-200, 300],
	"i1", [300, 1000],
	"i2", [1000, 2000],
	"i3", [2000, 5000],
	"i4", [5000, 15000],
	"i5", [15000, 40000]
	);
    
    my %out_edge = ();
    foreach $it (keys %inter){
	print STDERR $it."\t".$inter{$it}->[0]."\t".$inter{$it}->[1]."\n";
	my $OUT;
	open($OUT, ">$all_edge_file\_$it");
	$out_edge{$it} = $OUT;
    }

    open(FILE, $all_edge_file);
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	$dist = $line[4];
  	$support = $line[6];
  	foreach $it (keys %inter){
	    if($support >= 0 && $inter{$it}->[0] < $dist && $dist < $inter{$it}->[1]){
		$OUT = $out_edge{$it};
		print $OUT join("\t", @line)."\n";
		last;
	    }
	}
    }
    

    foreach $it (keys %inter){
	$OUT = $out_edge{$it};
	close($OUT);
    }
}




sub getAlignmentType {

    my ($rn, $cn, $ro, $co, $cs, $ce, $cl, $rs, $re, $rl) = @_;
    
    #the contig is fully contained in the read 
    #the second condition is made for short contig that required to have at leat [ cl - overlap ] based mapped ... NOT TOO RELAXED ... It is not better to use the minimum mapping length here ?
    if($ce-$cs >= $fraction*$cl || $ce-$cs >= $cl-$overlap) {
    #if($ce-$cs >= $fraction*$cl || $ce-$cs >= $cl-$min_alignment_length) {

        return("contig-contained");
    }
    
    #the read is fully contained in the read 
    elsif($re-$rs >= $fraction*$rl || $re-$rs >= $rl-$overlap) {
    #elsif($re-$rs >= $fraction*$rl || $re-$rs >= $rl-$min_alignment_length) {

	return("read-contained");
    }

    
    elsif($ce >= $cl-$overlap && $rs <= $overlap) {

	return("contig-at-start");
    }
    elsif($cs <= $overlap && $re >= $rl-$overlap) {
	
	return("contig-at-end");
    }
    else{            
	return("partial-match");
    }
}

sub printEdges {

    local(*alignments) = @_;
    my $next_allignemnet_to_link;
    my $edge_ID;my @contig_order;
    my $next_contig_found = 0;
    for($i = 0; $i <= $#alignments; $i++ ){
	#Conflicting alignement are filtered out
	next if(index($alignments[$i], "_CONFLICT") != -1);
	my ($rn1, $cn1, $ro1, $co1, $cs1, $ce1, $cl1, $rs1, $re1, $rl1) = split(/ /, $alignments[$i]);
	$alignemnt_to_link = $#alignments;
	$next_contig_found = 0;
	for( $j = $i + 1; $j <= $alignemnt_to_link; $j++ ){
	    #Conflicting alignement are filtered out
	    next if(index($alignments[$j], "_CONFLICT") != -1);
	    my ($rn2, $cn2, $ro2, $co2, $cs2, $ce2, $cl2, $rs2, $re2, $rl2) = split(/ /, $alignments[$j]);
	    
	    $distance = int(($rs2-$re1) + abs($rs2-$re1)*0.09); 
	    $sd = int(abs($distance)*0.1+50);
	    $distance += -$cs2 - ($cl1-$ce1);

	    #print "$alignments[$i]\n$alignments[$j]\n \n" if($cn1 eq $cn2); 
	    next if(
		$cn1 eq $cn2 || #the contig is the same
		($cl1 < $min_opera_contig_size || $cl2 < $min_opera_contig_size) ||#one contig does not pass opera contig size threshold
		exists $repeat_contig{$cn1} ||#one of the contig is a repeat
		exists $repeat_contig{$cn2} 
		#exists $edge_to_filter{$cn1.":".$cn2.":".$distance}#This a conflicting edge
		);
	    
	    @contig_order = sort($cn1, $cn2);
	    $edge_ID = join(" ", @contig_order);

	    #Collect the transitive edge for that read to filter out edge from reads that contain only the transitive edge
	    
	    $ori1 = ($co1 == 0 ? "+" : "-"); $ori2 = ($co2 == 0 ? "+" : "-");
	    
	    if(! defined $edges{$edge_ID}){
		my @tab = (0);
		$edges{$edge_ID} = \@tab;
	    }
	    #Update the stracture where we store the distance and the number of edges that support
	    $edges{join(" ", sort($cn1, $cn2))}->[0]++;
	    push(@{$edges{$edge_ID}}, $distance);
	    
	    if($component{$cn1} == 0 && $component{$cn2} == 0) {

		$component{$cn1} = $component{$cn2} = $component_num; 
		$member{$component_num} = "$cn1 $cn2";
		$length{$component_num} += $cl1 + $cl2; 
		$component_num++;
	    }
	    elsif($component{$cn1} == 0) {

		$component{$cn1} = $component{$cn2};
		$member{$component{$cn1}} .= " $cn1";
		$length{$component{$cn1}} += $cl1;
	    }
	    elsif($component{$cn2} == 0) {

		$component{$cn2} = $component{$cn1};
		$member{$component{$cn2}} .= " $cn2";
		$length{$component{$cn2}} += $cl2;
	    }
	    elsif($component{$cn1} != $component{$cn2}) {

		#print $length{$component{$cn1}}." with ".$length{$component{$cn2}}."\n";
		if($length{$component{$cn1}} >= $length{$component{$cn2}}) {
		    
		    $member{$component{$cn1}} .= " ".$member{$component{$cn2}};
		    $length{$component{$cn1}} += $length{$component{$cn2}}; $length{$component{$cn2}} = 0; 
		    foreach $member (split(/ /, $member{$component{$cn2}})) { $component{$member} = $component{$cn1}; }
		}
		else {

		    $member{$component{$cn2}} .= " ".$member{$component{$cn1}};
		    $length{$component{$cn2}} += $length{$component{$cn1}}; $length{$component{$cn1}} = 0; 
		    foreach $member (split(/ /, $member{$component{$cn1}})) { $component{$member} = $component{$cn2}; }
		}
	    }

	    #print EDGE "$cn1\t$ori1\t$cn2\t$ori2\t$distance\t$sd\t5\n" if($edges{join(" ", sort($cn1, $cn2))} == 1);
	    #print STDERR "NB edges ($cn1, $cn2): ".($edges{join(" ", sort($cn1, $cn2))}->[0])." $cn1\t$ori1\t$cn2\t$ori2\t$distance\t$sd\t";<STDIN>;
	    #There is no bundling 
	    #conflict between different distance etimate between 2 same contigs are NOT handled
	    #The first distance between 2 contigs is used for all the remianing alignemnts 
	    #All the other alignement will be considerd to support the first one
	    
	    if($edges{$edge_ID}->[0] == 1){
		$print{$edge_ID} = "$cn1\t$ori1\t$cn2\t$ori2\t$distance\t$sd\t";
		$str = "$cn1\t$ori1\t$cn2\t$ori2";
		$str = "$cn2\t$ori2\t$cn1\t$ori1" if($contig_order[0] ne $cn1);
		$edge_read_info{$edge_ID} = {"EDGE", $str,
					     "READ_LIST", [],
					     "COORD_CONTIG_1", [],
					     "COORD_CONTIG_2", [],
					     "COORD_CONTIG_1_ON_READ", [],
					     "COORD_CONTIG_2_ON_READ", []
		};
	    }

	    $others{$edge_ID} .= "$rn1,$cn1,$ori1,$cn2,$ori2,$distance,$sd|";

	    #Swap the alignement in case of contig order change
	    if($contig_order[0] ne $cn1){
		my ($rn1, $cn1, $ro1, $co1, $cs1, $ce1, $cl1, $rs1, $re1, $rl1) = split(/ /, $alignments[$j]);
		my ($rn2, $cn2, $ro2, $co2, $cs2, $ce2, $cl2, $rs2, $re2, $rl2) = split(/ /, $alignments[$i]);
	    }
	    push(@{$edge_read_info{$edge_ID}->{"READ_LIST"}}, $rn1);
	    #
	    push(@{$edge_read_info{$edge_ID}->{"COORD_CONTIG_1"}}, $co1."_".$cs1."_".$ce1);
	    push(@{$edge_read_info{$edge_ID}->{"COORD_CONTIG_2"}}, $co2."_".$cs2."_".$ce2);
	    #
	    push(@{$edge_read_info{$edge_ID}->{"COORD_CONTIG_1_ON_READ"}}, $ro1."_".$rs1."_".$re1);
	    push(@{$edge_read_info{$edge_ID}->{"COORD_CONTIG_2_ON_READ"}}, $ro2."_".$rs2."_".$re2);

	}
    }    
}

sub checkMapping{
    my ($mapFile, $all_edge_file) = @_;
    my ($rn, $cn, $score, $unused, $ro, $rs, $re, $rl, $co, $cs, $ce, $cl);
    my $currentScore = 0;my $previousScore = 0;
    %component = (); %length = (); $component_num = 1; %member = ();
    %edges = (); %print = (); %others = ();
    %edge_read_info = ();
    open(EDGE, ">$all_edge_file") or die $!;

    open(MAP, "$mapFile") or die $!;
    #open(MAP, "head -n10000 $mapFile | ") or die $!;
    open(STATUS, ">$mapFile.status") or die $!;

    $prev_rn = ""; @alignments = (); $prev_alignment_end = 0;
    while(<MAP>){
	chomp; 
	@line = split /\s+/; 
	next if(@line < 10);
	#Get the mapping coordinates using blars or graphmap
	#($rn, $cn, $ro, $co, $cs, $ce, $cl, $rs, $re, $rl, $score) = @line if($mapper eq "blasr");
	($rn, $cn, $ro, $co, , $score, $cs, $ce, $cl, $rs, $re, $rl) = @line if($mapper eq "blasr");
	#Conversion from the graphmmap (mhap) format to the blasr format
	if($mapper eq "graphmap"){
	    ($rn, $cn, $score, $unused, $ro, $rs, $re, $rl, $co, $cs, $ce, $cl) = @line;
	    $cn = $contig_id_to_name[$cn];
	}
	@data = ($rn, $cn, $ro, $co, $cs, $ce, $cl, $rs, $re, $rl, $score);

	#Init the contig component
	if(! defined $component{$cn}){
	    $component{$cn} = 0;
	}

	#print STDERR "$mapper -> @data"."\n".$rn."\t".$cl."\n\n";<STDIN>;
	
	#Filter small contigs
	if($cl < $min_contig_size) {
	    print STATUS (join("\t", @data)." | "); print STATUS "small-contig\n";
	    next;
	}
	
	#Minimum alignment length threshold
	if($re-$rs < $min_alignment_length) {
	    print STATUS (join("\t", @data)." | "); print STATUS "small-alignment\n";
	    next;
	}
	
	#Save the contig length for the conflicting edge pipeline
	$contig_length{$cn} = $cl;

	#Do we want to filter base on score as well ? At least for graphmap ...
	$alignmentType = &getAlignmentType(@data);
	
	print STATUS (join("\t", @data)." | "); print STATUS $alignmentType;
	if($alignmentType eq "partial-match") {

	    print STATUS "\n";
	    next;
	}

	if($rn ne $prev_rn) {
	    #Get the edge of the alignement on read prev_rn
	    &printEdges(*alignments) if(@alignments > 1);
	    #Udpdate the variable to strat to colloect information about the alignent on the read $rn
	    $prev_rn = $rn; @alignments = (); 
	    push @alignments, "@data";
	    $prev_alignment_end = $re + $cl-$ce;
	    $previousScore = $currentScore;   
	}
	else {

	    $curr_alignment_start = $rs - $cs;
	    $curr_alignment_end = $re + $cl-$ce;
	    $currentScore = $score;
	    #Non overlapping alignement
	    if($prev_alignment_end-$overlap < $curr_alignment_start) {
		
		push @alignments, "@data";
		$prev_alignment_end = $curr_alignment_end;	
		$previousScore = $currentScore;   
	    }
	    else {
		#WE SIMPLY FLAG THE ALIGNEMENT AS CONFLICTING AND LOOSE THE INFORMATION
		if($FLAG_FILTER_CONFLICTING_ALIGNEMENT){
		    $prev_alignment_end = $curr_alignment_end;
		    #print STDERR " *** ADD CONFLICT TAG ".(($alignments[@alignments-1]))."\n";#<STDIN>;
		    $alignments[@alignments-1] .= "_CONFLICT";
		    print STATUS " overlapped";
		}
		else{
		    #If 2 contig alignements are overlapping, take the one with the highest score
		    #DO WE WANT TO USE THE ALIGNEMENT LENGTH AS WELL TO COMPARE THE ALIGNEMENT ???
		    if( ($mapper eq "blasr" && $previousScore > $currentScore) || #the best score are the more negative
			($mapper eq "graphmap" && $previousScore < $currentScore) #here we use the the fraction of bases covered by seeds to compare the alignements
			){
			# replace mapping
			pop @alignments;
			push @alignments, "@data";
			$prev_alignment_end = $curr_alignment_end;	
			$previousScore = $currentScore; 
			print STATUS " overlapped better";
		    }
		    else{
			print STATUS " overlapped";
		    }
		}
	    }
	}

	print STATUS "\n";
    }
    
    #open(EDATA, ">$mapFile.edge_data") or die $!;
    foreach $key (keys(%print)) {
	print EDGE "$print{$key}$edges{$key}->[0]\n";
    }
    #close EDATA;

    close STATUS;
    close MAP;
    close EDGE;

    $total = 0;
    foreach $component (sort {$length{$b} <=> $length{$a}} keys(%length)) {

	#print "$length{$component}\n";
	$total += $length{$component};
	if($total > 5e6) {

	    print STDERR "N50: $length{$component}\n"; 
	    last;
	}
    }
}


sub CreateConfigFile{
    my( $contigFile, $suffix, @edgeFiles) = @_;
    
    open( CONF, ">config".$suffix ) or die $!;

    print CONF "#\n# Essential Parameters\n#\n\n";

    print CONF "# Output folder for final results\n";
    #print CONF "output_folder=results\n";
    #print CONF "output_folder=results_no_trans\n";
    print CONF "output_folder=results".$suffix."\n";
    #print CONF "output_folder=results_no_trans_no_conflict\n";
    #print CONF "output_folder=results_no_conflict\n";

    print CONF "# Contig file\n";
    print CONF "contig_file=$contigFile\n";

    print CONF "samtools_dir=$samtools_dir\n";

    print CONF "kmer=$kmer_size\n";

    #print CONF "cluster_threshold=2\n";
    print CONF "# Mapped read locations\n";

    #Update the config file to add the illumina mapping
    if(-e "${file_pref}.bam"){
	print CONF "[LIB]\n";
	print CONF "map_file=${file_pref}.bam\n";
	print CONF "cluster_threshold=$short_read_cluster_threshold\n";
    }
    
    $i = 0;
    @means = (300, 1000, 2000, 5000, 15000, 40000);
    @stds = (30, 100, 200, 500, 1500, 4000);
    foreach $edgeFileName ( @edgeFiles ){
	    print CONF "[LIB]\n";
	    #print CONF "cluster_threshold=$cluster_threshold\n";
	    print CONF "cluster_threshold=$cluster_threshold_tab[ $i ]\n";
	    print CONF "map_file=$edgeFileName\n";
	    print CONF "lib_mean=$means[ $i ]\n";
	    print CONF "lib_std=$stds[ $i ]\n";
	    print CONF "map_type=opera\n";
	    $i++;
    }
    
    close CONF;
}



sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";;
    print STDERR `$exe` if($run);
}
