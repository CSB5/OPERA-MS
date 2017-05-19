
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Spec;
use Switch;
use Getopt::Std;
use File::Spec::Functions qw(rel2abs);
use File::Basename;


#Grab the executing directory; turn all paths into absolute paths.
my $opera_ms_dir = dirname(rel2abs($0)) . "/";
print $opera_ms_dir  . "\n";

my $main_dir = getcwd;
$main_dir .= "\/";

#TODO Remove in release; used for testing.

my $opera_version;
my $runOperaMS_config_name = "runOperaMS.config"; #Generated config file, name can be whatever.
my $output_dir;
my $opera_ms_config_file;
my $long_read_file;
my $lr_output_dir;
my $illum_read1;
my $illum_read2;
my $contigs_file;
my $kmer_size = 60;
my $OPERAMS_config_name; 
my $opera_ms_cf;
my $config_option;
my $samtools_dir;
my $blasr_dir;
my $short_read_tool_dir;
my $short_read_maptool;
my $num_processor = 20;
my $help_line = "To configure OPERA-MS, please look at the example config file inside the OPERA-MS folder.\nUsage: \n\npath/OPERA-MS/OPERA-MS.pl <config_file>\n\nNote that the config file must be inside the path/OPERA-MS directory.\n";
my $incorrect_arguments="Please input the correct arguments. Refer to the documentation for more information or use:\n\n path/OPERA-MS/OPERA-MS.pl -h. \n\n";


if ( @ARGV == 0 ){
    die $incorrect_arguments; 
}

#read from config file
if ( @ARGV == 1){
    $OPERAMS_config_name = $ARGV[0];
    die $help_line if($ARGV[0] eq "-h");
    print STDERR "\nReading config file: ".$OPERAMS_config_name."\n";
    die"Config file does not exist. Exiting." if(!(-e $OPERAMS_config_name)); 

    open($opera_ms_cf, "<", $OPERAMS_config_name); 

    while(<$opera_ms_cf>) {
        next if (/^#/);  #skip comments
        chomp($_); 		 
        my @split_line = split('\s+', $_);
        if (@split_line != 0) {
            $config_option = $split_line[0];
            switch ($config_option) {
                
                #Empty cases must be checked so the same config file format
                #can be used for the runOperaMS generated config file. Empty cases indicate
                #that this specficpscript does not need those specific instances. 
                case "CONTIGS_FILE" {
                    $contigs_file = $split_line[1];
                }

                case "MAPPING_FILES"{
                }

                case "LIB"{
                }

                case "EDGE_BUNDLESIZE_THRESHOLD" {
                }

                case "OUTPUT_DIR" {
                    $output_dir = $split_line[1];
                }

                case "SAMTOOLS" {
                    $samtools_dir = $split_line[1];
                }

                case "BLASR_DIR"{
                    $blasr_dir = $split_line[1];
                }

                case "SHORT_READ_TOOL"{
                    $short_read_maptool = $split_line[1];
                }

                case "SHORT_READ_TOOL_DIR"{
                    $short_read_tool_dir = $split_line[1];
                }

                case "SIGMA_CONTIGS_FILE" {
             
                }

                case "KMER_SIZE" {
                    $kmer_size = $split_line[1];
                    if (!defined $kmer_size) {
                        die "KMER_SIZE not provided.\n";
                    }
                }

                case "CONTIGS_FILE_TYPE" {
                }

                case "CONTIG_LEN_THR" {
                }

                case "CONTIG_EDGE_LEN" {
                }

                case "CONTIG_WINDOW_LEN" {
                }

                case "PDIST_TYPE" {
                }

                case "SKIP_OPERA" {
               
                }

                case "LONG_READ_OUTPUT_DIR"{
                    $lr_output_dir = $split_line[1];
                }

                case "LONG_READ_FILE"{
                    $long_read_file = $split_line[1];
                }

                case "ILLUMINA_READ_1"{
                    $illum_read1 = $split_line[1];
                }

                case "ILLUMINA_READ_2"{
                    $illum_read2 = $split_line[1];
                }

                case "NUM_PROCESSOR"{
                    $num_processor = $split_line[1];
                }

                case "OPERA_VERSION"{
                    $opera_version = $split_line[1];
                }

                else {
                    die "Config option: ".$config_option." unknown, please check the config file. \nExiting. \n";
                }	
            }

        }

    }
     
}

else{
    die $incorrect_arguments;
    
  #  $output_dir=$ARGV[0];
  #      $long_read_file=$ARGV[1];    #/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/NANOPORE/LIBRARY/POOL/POOL_all/POOL.fa
  #  $lr_output_dir =$ARGV[2];    #OPERA_LG/OPERA-long-read/MEGAHIT/NANOPORE_ALL/
  #  $illum_read1=$ARGV[3];       #/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/ILLUMINA/mock20.R1.fastq.gz
  #  $illum_read2=$ARGV[4];       #/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/ILLUMINA/mock20.R2.fastq.gz
  #  $contigs_file=$ARGV[5];
  #  $kmer_size = 60;
  #  $num_processor = 20;
}

#
#Check if dependencies are found. Otherwise, die.
#
if(!defined($samtools_dir)){
    $samtools_dir = "";

    my $is_in_path = `which samtools`;
    die "samtools not found in path. Exiting." if (!$is_in_path);

}

else{
    $samtools_dir = $main_dir . $samtools_dir . "/" if(substr($samtools_dir, 0, 1) ne "/");

    if (! -e $samtools_dir . "/samtools") {
        die "Samtools not found at: ".$samtools_dir."\n";
    }

    $samtools_dir ="--samtools-dir " . $samtools_dir;

}

if(!defined($blasr_dir)){
    $blasr_dir = "";

    my $is_in_path = `which blasr`;
    die "blasr not found in path. Exiting." if (!$is_in_path);
}

else{
    $blasr_dir = $main_dir . $blasr_dir . "/" if (substr($blasr_dir, 0, 1) ne "/");
    if (! -e $blasr_dir . "/blasr") {
        die "blasr not found at: ".$blasr_dir."\n";
    }

    $blasr_dir = "--blasr " . $blasr_dir;
}

if(!defined($short_read_maptool)){
    $short_read_maptool = "bwa";
}

if(!defined($short_read_tool_dir)){
    $short_read_tool_dir= "";

    my $is_in_path = "which $short_read_maptool";
    die "$short_read_maptool not found in path. Exiting." if (!$is_in_path);

}

else{
    $short_read_tool_dir= $main_dir . $short_read_tool_dir. "/" if (substr($short_read_tool_dir, 0, 1) ne "/");
    if (! -e $short_read_tool_dir . "/$short_read_maptool") {
        die "short read mapping tool not found at : ".$short_read_tool_dir."\n";
    }
    
    $short_read_tool_dir = "--short-read-tooldir " . $short_read_tool_dir;
}

#TODO Kill bowtie support?
$short_read_maptool = "--short-read-maptool " . $short_read_maptool;
$output_dir = $main_dir . $output_dir."/" if (substr($output_dir, 0, 1) ne "/");

#Make paths absolute if they are not already.
$long_read_file = $main_dir . $long_read_file if(substr($long_read_file, 0, 1) ne "/");
$lr_output_dir = $main_dir . $lr_output_dir."/" if(substr($lr_output_dir, 0, 1) ne  "/");
$illum_read1 = $main_dir . $illum_read1 if(substr($illum_read1, 0, 1) ne "/");
$illum_read2= $main_dir . $illum_read2 if(substr($illum_read2, 0, 1) ne "/");
$contigs_file = $main_dir . $contigs_file if(substr($contigs_file, 0, 1) ne "/");
$opera_ms_config_file = $output_dir."/".$runOperaMS_config_name;

if( ! -e $contigs_file) {
    die "\nError : $contigs_file - contig file does not exist\n";
}

if( ! -e $illum_read1) {
    die "\nError : $illum_read1 - illumina read 1 file does not exist\n";
}

if( ! -e $illum_read2) {
    die "\nError : $illum_read2 - illumina read 2 file does not exist\n";
}

if( ! -e $long_read_file) {
    die "\nError : $long_read_file - long read file does not exist\n";
}


    
#my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -V -l mem_free=20G,h_rt=500:0:0 -pe OpenMP 20";
my $command;
my $hold;
my $hold2;

my $contigs_file_sed = $contigs_file;
$contigs_file_sed =~ s/\//\\\//g;      #replace / with \/
my $output_dir_sed = $output_dir;
$output_dir_sed =~ s/\//\\\//g;      #replace / with \/


### opera config
#

if (-d $output_dir){
     $command = "rm -r $output_dir; mkdir $output_dir";
}

else {
    $command = "mkdisdfjsdfr $output_dir";
}

$command="cat ${OPERAMS_config_name} | sed 's/CONTIGS_FILE .*/CONTIGS_FILE $contigs_file_sed/' | sed 's/MAPPING_FILES .*/MAPPING_FILES $output_dir_sed\\/contigs\.bam/' | sed 's/LIB .*/LIB $output_dir_sed\\/contigs\.bam/' | sed 's/OUTPUT_DIR .*/OUTPUT_DIR $output_dir_sed/' | sed 's/KMER_SIZE .*/KMER_SIZE $kmer_size/' > $opera_ms_config_file";
run_exe($command);





### Run opera-lr
#
#if(!(-d $lr_output_dir)){
$command = "mkdir -p $lr_output_dir";
run_exe($command);

#Use for:
#mapping of short reads => bwa
#mapping of long reads => blasr
#compute long read links
#run opera NO NEED <= DONE : Jim 
#TODO be careful OPERA-LG_v2.1.0 name
#===> add an option --link that does not run opera <= DONE : Jim. Set --skip-opera = 1
print STDERR "\n-----STARTING LONG READ PROCESSING-----\n";
$command= "${opera_ms_dir}${opera_version}/bin/OPERA-long-read.pl --contig-file $contigs_file --kmer $kmer_size --long-read-file $long_read_file --output-prefix opera --output-directory $lr_output_dir --num-of-processors $num_processor --opera ${opera_ms_dir}${opera_version}/bin/ --illumina-read1 $illum_read1 --illumina-read2 $illum_read2 $samtools_dir $blasr_dir $short_read_tool_dir --skip-opera 1";
#$command= "$qsub -N opera-lr -b y /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/OPERA-long-read.pl --contig-file /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/MOCK_20/ASSEMBLY/MEGAHIT/final.contigs_soap.fa --kmer 100 --long-read-file /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/NANOPORE/LIBRARY/POOL/POOL_all/POOL.fa --output-prefix opera --output-directory OPERA_LG/OPERA-long-read/MEGAHIT/NANOPORE_ALL/ --num-of-processors 20 --opera /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/"
my $long_read_fail = run_exe($command);

if($long_read_fail){
    print "Error in the OPERA-long-read.pl script. Please see log for details.";
}
#}

### Softlink opera.bam
#
#NOT REQUIRED CAN BE REMOVE IF REPLACE BY PATH IN $lr_output_dir/opera.bam
$command="ln -s $lr_output_dir/opera.bam $output_dir/contigs.bam";
run_exe($command);

#GENERATE OperaMS-run.sh CONFIG FILE


### runOperaMS
##
#$command = "perl bin/OperaMS-run.sh $opera_ms_config_file $output_dir";
#SIGMA on short reads
#bundling of short reads (can be use if fix the problem with short reads)
#Run SIGMA
#Run OPERA NO NEED
#===> add an option --no-opera that does not run opera
print STDERR "\n-----STARTING SIGMA AND BUNDLING-----\n";
$command = "perl ${opera_ms_dir}bin/runOperaMS.pl $opera_ms_config_file 2>&1";
run_exe($command);

my $lr_output_dir_sed = $lr_output_dir;
$lr_output_dir_sed =~ s/\//\\\//g;	#replace / with \/
$output_dir_sed = $output_dir;
$output_dir_sed =~ s/\//\\\\\\\\\\\//g; #replace / with \\\\\/ => \\\\ and \/
my $output_dir_sed2 = $output_dir;
$output_dir_sed2 =~ s/\//\\\//g; #replace / with \/



### Run sigma for long read ###
#
# SIGMA1
$command="mkdir -p $output_dir/contigs/sigma-long-read";
run_exe($command);
#
#long read sigma config
# SIGMA2
$command="ls $lr_output_dir/pairedEdges_i* | tr '\\n' ',' | sed \"s/$lr_output_dir_sed\\\//$output_dir_sed\\\\\\\\\\\/contigs\\\\\\\\\\\/sigma-long-read\\\\\\\\\\\//g\" | sed 's/\\(\\\\\\\/pairedEdges_i[0-9]*\\)/\\1\\1/g'";
my $edges=run_exe($command);
chop $edges;
print STDERR "$edges\n";
# SIGMA3
$command="cat $output_dir/contigs/sigma/sigma.config | sed 's/\\(mapping_files=.*\\)/\\#\\1/' | sed 's/output_dir=.*/output_dir=${output_dir_sed2}\\/contigs\\/sigma-long-read\\//' | sed 's/\\(edges_files=\\).*/\\1$edges/' > $output_dir/contigs/sigma-long-read/sigma.config";
run_exe($command);
#
#create pairededges directory
# SIGMA4
$command="for i in `ls $lr_output_dir/pairedEdges_i*`; do base=`basename \$i`; fullpath=`readlink -f \$i` ;mkdir $output_dir/contigs/sigma-long-read/\$base; ln -s \$fullpath $output_dir/contigs/sigma-long-read/\$base/\$base; done";
run_exe($command);
#
#create lib.txt in pairededges directory
# SIGMA5
$command="for i in `ls -d $output_dir/contigs/sigma-long-read/pairedEdges_i*`; do cp $output_dir/*_bundles/lib.txt \$i; done";
run_exe($command);
# SIGMA6
$command="for j in `ls -d $output_dir/contigs/sigma-long-read/pairedEdges_i*`; do i=`echo \$j | sed 's/.*\\(_i.*\\)/\\1/g'`; mean=`grep -A1 \$i $lr_output_dir/config | tail -1 | sed 's/lib_mean=\\(.*\\)/\\1/g'` ; dev=`grep -A2 \$i $lr_output_dir/config | tail -1 | sed 's/lib_std=\\(.*\\)/\\1/g'` ; sed -i -e 's/\\(Mean length of the library is: \\).*/\\1'\"\${mean}\"'/g' -e 's/\\(Standard deviation of the library is: \\).*/\\1'\"\${dev}\"'/g' \${j}/lib.txt; done";
run_exe($command);
#
#finally, run sigma
# SIGMA7
$command="${opera_ms_dir}bin/sigma $output_dir/contigs/sigma-long-read/sigma.config";
run_exe($command);


#Filter that remove contigs with a coverage 1.5 times higher than the mean
$command="outcontig=`ls $output_dir/contigs/sigma/contigs_\*`; ${opera_ms_dir}bin/filter_cluster_coverage.pl `echo \$outcontig` $output_dir/contigs/sigma-long-read/clusters $output_dir/contigs/sigma-long-read/ 1.5 $output_dir/contigs/sigma-long-read/NO_REPEAT";
run_exe($command);


### Run opera after sigma for long read ###
#
# OPERA1
$command="mkdir $output_dir/contigs/scaffolds-long-read";
run_exe($command);
#
#long read opera config
# OPERA2
$command="cp $output_dir/contigs/scaffolds/opera.config $output_dir/contigs/scaffolds-long-read";
run_exe($command);
#
#chg outdir
# OPERA3
$command="sed -i 's/output_folder=.*/output_folder=${output_dir_sed2}\\/contigs\\/scaffolds-long-read\\//' $output_dir/contigs/scaffolds-long-read/opera.config";
run_exe($command);
#
# OPERA4
#comment away short read edge library
$command="sed -i -e 's/\\(\\[LIB\\]\\)/\\#\\1/g' -e 's/\\(map_type\\)/\\#\\1/g' -e 's/\\(map_file\\)/\\#\\1/g' -e 's/\\(lib_mean\\)/\\#\\1/g' -e 's/\\(lib_std\\)/\\#\\1/g' $output_dir/contigs/scaffolds-long-read/opera.config";
run_exe($command);
#
#add new edge libraries
# OPERA5
$command="for j in `ls -d $output_dir/contigs/sigma-long-read/pairedEdges_i*`; do  i=`echo \$j | sed 's/.*\\(_i.*\\)/\\1/g'`; mean=`grep -A1 \$i $lr_output_dir/config | tail -1 | sed 's/lib_mean=\\(.*\\)/\\1/g'` ; dev=`grep -A2 \$i $lr_output_dir/config | tail -1 | sed 's/lib_std=\\(.*\\)/\\1/g'`; base=`basename \$j`; if [[ -s $output_dir/contigs/sigma-long-read/filtered_\${base} ]] ; then sed -i -e '\$a\\[LIB]\\' -e 'map_type=opera\\' -e 'map_file=$output_dir\/contigs\/sigma-long-read\/NO_REPEAT\/filtered_'\"\$base\"'\\' -e 'lib_mean='\"\$mean\"'\\' -e 'lib_std='\"\$dev\"'\\'  $output_dir/contigs/scaffolds-long-read/opera.config; fi; done";

#$command="for j in `ls -d $output_dir/contigs/sigma-long-read/pairedEdges_i*`; do  i=`echo \$j | sed 's/.*\\(_i.*\\)/\\1/g'`; mean=`grep -A1 \$i $lr_output_dir/config | tail -1 | sed 's/lib_mean=\\(.*\\)/\\1/g'` ; dev=`grep -A2 \$i $lr_output_dir/config | tail -1 | sed 's/lib_std=\\(.*\\)/\\1/g'`; base=`basename \$j`; if [[ -s $output_dir/contigs/sigma-long-read/filtered_\${base} ]] ; then sed -i -e '\$a\\[LIB]\\' -e 'map_type=opera\\' -e 'map_file=$output_dir\/contigs\/sigma-long-read\/filtered_'\"\$base\"'\\' -e 'lib_mean='\"\$mean\"'\\' -e 'lib_std='\"\$dev\"'\\'  $output_dir/contigs/scaffolds-long-read/opera.config; fi; done";
run_exe($command);
#
# OPERA 6
#finally, run opera
$command="${opera_ms_dir}bin/OPERA-LG $output_dir/contigs/scaffolds-long-read/opera.config > $output_dir/contigs/scaffolds-long-read/log.txt";
run_exe($command);


$command="rm $output_dir/scaffoldSeq.fasta; rm $output_dir/scaffoldSeq.fasta.stats; ln -s $output_dir/contigs/scaffolds-long-read/scaffoldSeq.fasta $output_dir/scaffoldSeq.fasta";
run_exe($command);

$command="less ${output_dir}/contigs/scaffolds-long-read/log.txt";
run_exe($command);




sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return . "\n" if($run);
    return $return;
}


