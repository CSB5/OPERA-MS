#!/usr/bin/perl
use strict;
use Switch;
use warnings;
use Cwd;

my $main_direct = getcwd;

my $cmnd = "rm -r TEST_RUN/* ; rm -r long-read";
run_exe($cmnd);

if ( @ARGV < 6 ){
    print "Please input the correct arguments. Refer to the documentation for more information. \n";
    exit 0;
}


#my $cmnd = "rm error.log; touch error.log";

#run_exe($cmnd);

$main_direct .= "\/";
my $output_dir=$ARGV[0];
my $long_read_file=$ARGV[1];  
my $lr_output_dir =$ARGV[2]; 
my $illum_read1=$ARGV[3];   
my $illum_read2=$ARGV[4];  
my $contigs_file=$ARGV[5];



#Make path absolute if user specifies a relative path.
$output_dir = $main_direct . $output_dir if (substr($output_dir, 0, 1) ne "/");
$long_read_file = $main_direct . $long_read_file if(substr($long_read_file, 0, 1) ne "/");
$lr_output_dir = $main_direct . $lr_output_dir if(substr($lr_output_dir, 0, 1) ne  "/");
$illum_read1 = $main_direct . $illum_read1 if(substr($illum_read1, 0, 1) ne "/");
$illum_read2= $main_direct . $illum_read2 if(substr($illum_read2, 0, 1) ne "/");
$contigs_file = $main_direct . $contigs_file if(substr($contigs_file, 0, 1) ne "/");

if( ! -e $contigs_file) {
    die "\nError : $contigs_file - contig file does not exist\n"
};

if( ! -e $illum_read1) {
    die "\nError : $illum_read1 - illumina read 1 file does not exist\n"
}

if( ! -e $illum_read2) {
    die "\nError : $illum_read2 - illumina read 2 file does not exist\n"
}

if( ! -e $long_read_file) {
    die "\nError : $long_read_file - long read file does not exist\n"
}

my $command = "perl bin/OPERA-MS.pl $output_dir $long_read_file $lr_output_dir $illum_read1 $illum_read2 $contigs_file";
run_exe($command);

sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    my $return = `$exe`;
    print STDERR "\n".$exe."\n";
    print STDERR $return if($run);
    return $return;
}


#qsub -terse -m a -M $USER_PRINCIPAL_NAME -cwd -V -l mem_free=20G,h_rt=48:0:0 -pe OpenMP 20 -N OPERA-MS -e error.log -o OPERA-MS.log -b y perl bin/OPERA-MS.pl /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/SOFWARE/OPERA-MS/TEST_RUN /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/SOFWARE/OPERA-MS/test_files/POOL.fa /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/SOFWARE/OPERA-MS/TEST_RUN/long-read/ /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/SOFWARE/OPERA-MS/test_files/mock1000.R1.fastq.gz  /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/SOFWARE/OPERA-MS/test_files/mock1000.R2.fastq.gz /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/SOFWARE/OPERA-MS/test_files/final.contigs_soap.fa  
##
