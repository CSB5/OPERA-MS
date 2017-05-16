
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Spec;
use Switch;
use Getopt::Std;

###to update code for multi lib

my $output_dir=$ARGV[0];
my $opera_ms_config_file = $output_dir."/OPERA-MS.config";
my $long_read_file=$ARGV[1];    #/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/NANOPORE/LIBRARY/POOL/POOL_all/POOL.fa
my $lr_output_dir =$ARGV[2];    #OPERA_LG/OPERA-long-read/MEGAHIT/NANOPORE_ALL/
my $illum_read1=$ARGV[3];       #/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/ILLUMINA/mock20.R1.fastq.gz
my $illum_read2=$ARGV[4];       #/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/ILLUMINA/mock20.R2.fastq.gz
my $contigs_file=$ARGV[5];
my $kmer_size = 60;

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
$command="cat OPERA-MS.config | sed 's/CONTIGS_FILE .*/CONTIGS_FILE $contigs_file_sed/' | sed 's/MAPPING_FILES .*/MAPPING_FILES $output_dir_sed\\/contigs\.bam/' | sed 's/LIB .*/LIB $output_dir_sed\\/contigs\.bam/' | sed 's/OUTPUT_DIR .*/OUTPUT_DIR $output_dir_sed/' | sed 's/KMER_SIZE .*/KMER_SIZE $kmer_size/' > $opera_ms_config_file";
print STDERR "$command\n";
run_exe($command);



$command = "rm $output_dir/run_log.txt";
run_exe($command);




### Run opera-lr
#
#if(!(-d $lr_output_dir)){
$command = "mkdir -p $lr_output_dir";
run_exe($command);

$command= "OPERA-LG_v2.1.0/bin/OPERA-long-read.pl --contig-file $contigs_file --kmer $kmer_size --long-read-file $long_read_file --output-prefix opera --output-directory $lr_output_dir --num-of-processors 20 --opera /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.1.0/bin/ --illumina-read1 $illum_read1 --illumina-read2 $illum_read2";
#$command= "$qsub -N opera-lr -b y /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/OPERA-long-read.pl --contig-file /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/MOCK_20/ASSEMBLY/MEGAHIT/final.contigs_soap.fa --kmer 100 --long-read-file /home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/DATA/MOCK_20/NANOPORE/LIBRARY/POOL/POOL_all/POOL.fa --output-prefix opera --output-directory OPERA_LG/OPERA-long-read/MEGAHIT/NANOPORE_ALL/ --num-of-processors 20 --opera /home/bertrandd/PROJECT_LINK/OPERA_LG/OPERA_LONG_READ/OPERA-LG_v2.0.6/bin/"
run_exe($command);
#}

### Softlink opera.bam
#
$command="ln -s $lr_output_dir/opera.bam $output_dir/contigs.bam";
run_exe($command);

#GENERATE OperaMS-run.sh CONFIG FILE


### runOperaMS
##
$command = "perl bin/OperaMS-run.sh $opera_ms_config_file $output_dir";
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
$command="bin/sigma $output_dir/contigs/sigma-long-read/sigma.config";
run_exe($command);


#Filter that remove contigs with a coverage 1.5 times higher than the mean
$command="outcontig=`ls $output_dir/contigs/sigma/contigs_\*`; bin/filter_cluster_coverage.pl `echo \$outcontig` $output_dir/contigs/sigma-long-read/clusters $output_dir/contigs/sigma-long-read/ 1.5 $output_dir/contigs/sigma-long-read/NO_REPEAT";
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
$command="bin/opera $output_dir/contigs/scaffolds-long-read/opera.config > $output_dir/contigs/scaffolds-long-read/log.txt";
run_exe($command);


$command="rm $output_dir/scaffoldSeq.fasta; rm $output_dir/scaffoldSeq.fasta.stats; ln -s $output_dir/contigs/scaffolds-long-read/scaffoldSeq.fasta $output_dir/scaffoldSeq.fasta";
run_exe($command);




sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";
    my $return = `$exe` if($run);
    print STDERR $return if($run);
    return $return;
}

