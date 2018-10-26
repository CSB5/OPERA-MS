#!/usr/bin/perl
use Getopt::Long;

#$path = "/mnt/software/stow/bwa-0.7.10/bin/";
# $path = "/mnt/software/stow/bwa-0.7.4/bin/";

# preprocess of reads and map them using bowtie

my $help_message = "

	Options:
    --nproc             Number of processor used during the analysis 
    --contig		Multi-fasta contig file
    --illumina-read1	Fasta/Fastq file of read 1 of the paired-end reads
    --illumina-read2	Fasta/Fastq file of read 2 of the paired-end reads
    --out		Name of the output file
    --map-tool		Mapping tool can be either bwa (default) or bowtie
    --tool-dir		Directory that contains binaries to the chosen mapping tool (default PATH)
    --samtools-dir		Directory that contains samtools binaries (default PATH)
    --help			Help
";

my $path = "";
my $samtoolsDir = "";

GetOptions(
    "nproc=i"    => \$nproc,
    "contig=s"    => \$contigFile,
    "illumina-read1=s"    => \$readFile1,
    "illumina-read2=s" => \$readFile2,
    "out=s" => \$outputFile,
    "map-tool=s"      => \$mapTool,
    "tool-dir=s"      => \$path,
    "samtools-dir=s"  => \$samtoolsDir,
    "help"	=> \$flag_help

) or die("Error in command line arguments.\n");



if ( $flag_help ){
        print $help_message;
        exit 0;
}
if ( !(defined($contigFile)) ){
        print "Multi-fasta contig file missing.\n";
        exit 0;
}
if ( !(defined($readFile1)) ){
        print "Fasta/Fastq file of read 1 of the paired-end reads missing.\n";
        exit 0;
}
if ( !(defined($readFile2)) ){
        print "Fasta/Fastq file of read 2 of the paired-end reads missing.\n";
        exit 0;
}
if ( !(defined($outputFile)) ){
        print "Name of the output file missing.\n";
        exit 0;
}
#if ( ($mapTool eq "bwa" || !(defined($mapTool))) && defined($bowtieDir) ){
#        print "Bowtie directory not needed and ignored for bwa mapping.\n";
#        exit 0;
#}
#if ( $mapTool eq "bowtie" && defined($bwaDir) ){
#        print "Bwa directory not needed and ignored for bowtie mapping.\n";
#        exit 0;
#}

$mapTool = bwa if( !(defined($mapTool)) );

#$path = $bwaDir if( ($mapTool eq "bwa" || !(defined($mapTool))) && defined($bwaDir) );
#$path = $bowtieDir if( $mapTool eq "bowtie" && defined($bowtieDir) );

#To make sure that this is a directory
if( $path !~ "/\$" and $path ne "" )
{
	$path .= "/";
}
#
if( $samtoolsDir !~ "/\$" and $samtoolsDir ne "" )
{
	$samtoolsDir .= "/";
}


# create output folder
&createOutputFolder;

# check read format
&checkReadFormat;

if( $mapTool eq "bowtie" )
{
    # build the index 
    &buildIndexUsingBowtie;

    # map with bowtie
    &mapWithBowtie;

    # read the read file and map
    &readAndMap;

    &clearBowtie;
}
elsif( $mapTool eq "bwa" )
{
    # build the index
    &buildIndexUsingBwa;

    # find the SA coordinates
    &findSAWithBwa;
    &readAndMap;

    # generate alignment
    &generateAlignmentUsingBwa;
    &readAndMap;

    # remove the intermediate file
    `rm ${folder}$outputFile\_read.sai`;
}
else
{
    die "The mapping tool should be either bwa or bowtie.\n";
    
}

# finalize
$time = localtime;
print "[$time]\t";
print "Preprocessing done!\n";



sub createOutputFolder
{
    @dir = split( "/", $outputFile );
    $folder = "";
    for( $i = 0; $i < @dir - 1; $i++ )
    {
	$folder .= "$dir[ $i ]/";
    }
    $outputFile = $dir[ @dir - 1 ];
    if( $folder ne "" )
    { system("mkdir -p $folder"); }
    
    print "PREPROCESS:\n";
}

sub checkReadFormat
{
    if(  $readFile1 =~ /\.gz$/ ){
	# open gzip file
	open( READ, "gunzip -c $readFile1 |" ) or die $!;
    }
    else{
	open( READ, "$readFile1" ) or die $!;
    }

    while( $line = <READ> )
    {
	if( $line =~ /^#/ )
	{ next; }
	if( $line =~ "^>" )
	{ $fasta = "-f"; print "Fasta format is recognized\n"; last; }
	elsif( $line =~ "^@" )
	{ $fasta = "-q"; print "Fastq format is recognized\n"; last; }
	else
	{ print "reads format is not correct\n"; exit( 0 ); }
    }
    
    $type = "fs";
    $line = <READ>;
    if( $line =~ m/(0|1|2|3)/ )
    {
	$type = "cs";
    }
    
    close READ;
}

sub buildIndexUsingBowtie
{
    # build the index 
    $time = localtime;
    print "[$time]\t";
    print "Building bowtie index...\n";
    @contigName = split( "/", $contigFile );
    if( ! -f "$folder$contigName[ -1 ].1.ebwt" )
    {
	if( $type eq "cs" )
	{
	    $command = "${path}bowtie-build -C $contigFile $folder$contigName[ -1 ]";
	    `$command` or die "ERROR! Running bowtie-build error.\n";
	}
	else
	{
	    $command = "${path}bowtie-build $contigFile $folder$contigName[ -1 ]";
	    `$command` or die "ERROR! Running bowtie-build error.\n";
	}
    }
    else
    { # no need to re-build
	print "\t\t\t\tIndex already exist. Skipping building index...\n";
    }
}

sub mapWithBowtie
{
    $time = localtime;
    print "[$time]\t";
    print "Mapping reads using bowtie...\n";
    if( $type eq "cs" )
    {
	$command = "${path}bowtie -v 3 -a -m 1 -S -t -C $fasta -p 15 $folder$contigName[ -1 ] - 2>${folder}bowtie.err | sort -n > $folder$outputFile";
    }
    else
    {
	$command = "${path}bowtie -v 3 -a -m 1 -S -t $fasta -p 15 $folder$contigName[ -1 ] - 2>${folder}bowtie.err | sort -n > $folder$outputFile";
    }
}

sub buildIndexUsingBwa
{
    $time = localtime;
    print "[$time]\t";
    print "Building bwa index...\n";
    @contigName = split( "/", $contigFile );
    if( ! -f "$folder$contigName[ -1 ].amb" )
    {
	if( $type eq "cs" )
	{
	    $command = "${path}bwa index -c -p $folder$contigName[ -1 ] $contigFile";
	    `$command`; #or die "Error during bwa index creation.\n";
	}
	else
	{
	    $command = "${path}bwa index -p $folder$contigName[ -1 ] $contigFile";
	    print $command."\n";
	    `$command`; #or die "Error during bwa index creation.\n";
	}
    }
    else
    { # no need to re-build
	print "\t\t\t\tIndex already exist. Skipping building index...\n";
    }
}

sub readAndMap
{
    if(  $readFile1 =~ /\.gz$/ ){
	# open gzip file
	open( F, "gunzip -c $readFile1 |" ) or die $!;
	open( S, "gunzip -c $readFile2 |" ) or die $!;
    }
    else{
	open( F, "$readFile1" ) or die $!;
	open( S, "$readFile2" ) or die $!;
    }
    
    open( OUTPUT, "|-", $command ) or die $!;
    $num = 1;
    while( $first = <F> )
    {
	if( $first =~ /^#/ )
	{ next; }

	if( $fasta eq "-f" )
	{
	    print OUTPUT ">$num".".1\n";
	    $first = <F>;
	    print OUTPUT $first;
	    
	    while( $second = <S> )
	    {
		if( $second !~ /^#/ )
		{ last; }
	    }
	    
	    # $second = <S>;
	    print OUTPUT ">$num".".2\n";
	    $second = <S>;
	    print OUTPUT $second;
	}
	else
	{ # fastq
	    print OUTPUT "\@$num".".1\n";
	    $first = <F>;
	    print OUTPUT $first;
	    $first = <F>;
	    print OUTPUT "+$num".".1\n";
	    $first = <F>;
	    print OUTPUT $first;
	    
	    while( $second = <S> )
	    {
		if( $second !~ /^#/ )
		{ last; }
	    }
	    # $second = <S>;
	    print OUTPUT "\@$num".".2\n";
	    $second = <S>;
	    print OUTPUT $second;
	    $second = <S>;
	    print OUTPUT "+$num".".2\n";
	    $second = <S>;
	    print OUTPUT $second;
	}
	$num++;
    }
    close S;
    close F;
    close OUTPUT;
}

sub findSAWithBwa
{
    # align with bwa
    $time = localtime;
    print "[$time]\t";
    print "Finding the SA coordinates of the reads using BWA aln...\n";
    $command = "${path}bwa aln -t $nproc $folder$contigName[ -1 ] - > ${folder}$outputFile\_read.sai";
    if($?){
	die "Error during bwa mapping. Please see log for details.\n";
    }
}

sub generateAlignmentUsingBwa
{
    $time = localtime;
    print "[$time]\t";
    print "Generate alignments of reads using bwa samse...\n";
    # $command = "${path}bwa samse -n 1 $folder$contigName[ -1 ] ${folder}read.sai - | awk \'{ if( \$3 != \"*\" ) print \$0 }\' > $folder$outputFile";
    # $command = "${path}bwa samse -n 1 $folder$contigName[ -1 ] ${folder}read.sai - | grep '\\(^@\\|XT:A:U\\)' | samtools view -S -h -b -F 0x4 - | samtools sort -no - ${folder}temporarySam | samtools view -h -b - > $folder$outputFile";
    #$command = "${path}bwa samse -n 1 $folder$contigName[ -1 ] ${folder}read.sai - | grep '\\(^@\\|XT:A:U\\)' | ${samtoolsDir}samtools view -S -h -b -F 0x4 - | ${samtoolsDir}samtools sort -\@ 20 -no - ${folder}temporarySam > $folder$outputFile";
    $command = "${path}bwa samse -n 1 $folder$contigName[ -1 ] ${folder}$outputFile\_read.sai - | grep '\\(^@\\|XT:A:U\\)' | ${samtoolsDir}samtools view -S -h -b -F 0x4 - | ${samtoolsDir}samtools sort -\@ $nproc -no - ${folder}$outputFile\_temp > $folder$outputFile";
    print $command."\n";
    if($?){
	die "Error during bwa mapping. Please see log for details.\n";
    }
}

sub clearBowtie
{
    `rm ${folder}bowtie.err`;
}
