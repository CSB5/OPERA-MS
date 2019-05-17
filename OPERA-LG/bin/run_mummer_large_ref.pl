#!/usr/bin/perl

use warnings; 
use Cwd; 

my ($ref_file, $split_ref_dir, $query_file, $split_query_dir, $dir_res, $final_res_file, $nb_process, $mummer_path) = @ARGV;

#http://sourceforge.net/mailarchive/forum.php?thread_name=6A5AFBFC960692419DB59C142EB027B503260DE5%40ukaprdembx01.rd.astrazeneca.net&forum_name=mummer-help
$ref_file_size = `grep -v ">" $ref_file | wc -c`;chop $ref_file_size;

$max_split_ref_size   = $ref_file_size/$nb_process;#200000000;#400000000 ;#536870908



#$max_split_ref_size   = 4000000 ;#536870908
$max_split_query_size = 3000000000;#4294967295
$max_genome_number = 300000;
run_exe("rm -r $split_query_dir $split_ref_dir");
#SPLIT the file
if(! -d $split_ref_dir){
    run_exe("mkdir $split_ref_dir");
    split_fasta_file($ref_file, $split_ref_dir, $max_split_ref_size);
    #print STDERR " *** Splitting completed !!!\n";<STDIN>;
}
#print STDERR " *** finish !!!\n";<STDIN>;
if(! -d $split_query_dir){
    run_exe("mkdir $split_query_dir");
    split_fasta_file($query_file, $split_query_dir, $max_split_query_size);
}


print STDERR " *** run the mummer mapping\n";
$max_proc = $nb_process;

opendir(DIR, "$split_ref_dir");
@split_ref = readdir(DIR);
close(DIR);

opendir(DIR, "$split_query_dir");
@split_query = readdir(DIR);
close(DIR);


run_exe("mkdir $dir_res") unless -d $dir_res;
my $split_mummer_dir = "$dir_res/mummer";
run_exe("rm -r $split_mummer_dir");
run_exe("mkdir $split_mummer_dir");# unless -d $split_mummer_dir;
open($OUT_delta, ">$split_mummer_dir/cmd_delta.txt");
open($OUT_delta_filter, ">$split_mummer_dir/cmd_delta_filter.txt");
open($OUT_coord, ">$split_mummer_dir/cmd_coord.txt");

foreach $ref (@split_ref){
    if(index($ref, ".fa") != -1){

	foreach $query (@split_query){
	    if(index($query, ".fa") != -1){
		print $OUT_delta "$mummer_path/nucmer --nosimplify --maxmatch -p $split_mummer_dir/$ref\_$query $split_ref_dir/$ref $split_query_dir/$query > /dev/null 2> /dev/null \n";
		#Filtet on identity to avoid alignement comming from contig of other species or repeat
		print $OUT_delta_filter "$mummer_path/delta-filter -i95 -r $split_mummer_dir/$ref\_$query.delta > $split_mummer_dir/$ref\_$query.delta.filter\n";
		print $OUT_coord "$mummer_path/show-coords -T -l  $split_mummer_dir/$ref\_$query.delta.filter > $split_mummer_dir/$ref\_$query.coord\n";
		#
		#print $OUT_coord "$mummer_path/show-coords -I 95 -r -l -c -T $split_mummer_dir/$ref\_$query.delta > $split_mummer_dir/$ref\_$query.coord 2> /dev/null \n";
		#print $OUT_coord "$mummer_path/show-coords -r -l -c -T $split_mummer_dir/$ref\_$query.delta > $split_mummer_dir/$ref\_$query.coord 2> /dev/null \n";
	    }
	}
    }
}

close $OUT_delta;
close $OUT_coord;

#run the stuff
$cmd = "cat $split_mummer_dir/cmd_delta.txt | xargs -I cmd --max-procs=$max_proc bash -c cmd > /dev/null \n";
run_exe($cmd);

$cmd = "cat $split_mummer_dir/cmd_delta_filter.txt | xargs -I cmd --max-procs=$max_proc bash -c cmd > /dev/null \n";
run_exe($cmd);

$cmd = "cat $split_mummer_dir/cmd_coord.txt | xargs -I cmd --max-procs=$max_proc bash -c cmd > /dev/null \n";
run_exe($cmd);


#merging the cood file
print STDERR " *** merge the coord file\n";#<STDIN>;exit(0);
$first = 1;


opendir(DIR,$split_mummer_dir);
@split_coord = readdir(DIR);
close(DIR);


#my $final_res_file = "$final_res_file";
foreach $coord (@split_coord){
    $c_f = "$split_mummer_dir/$coord";
    if(index($c_f, "cmd") == -1 &&
       index($c_f, "coord") != -1){
	if($first){
	    run_exe("cp $c_f $final_res_file");
	    $first = 0;
	}
	else{
	    #to remove the 4 header line
	    run_exe ("sed \'1,4 d\' $c_f >> $final_res_file");
	}
	
	#run_exe ("rm $c_f");
    }
}

#$cmd = "cat $split_mummer_dir/*.coords > $dir_res/contig.coords";
#run_exe($cmd);


sub split_fasta_file{
    my ($fasta_file, $split_dir, $max_split_size) = @_;
    
    open(FILE, $fasta_file);
    my $genome;
    my $genome_number = 0;
    my $genome_size = 0;
    my $file_size = 0;
    my $OUT;
    my $split_file_num = 1;
    open($OUT, ">$split_dir/split_$split_file_num.fa");
    
    print STDERR " *** split $fasta_file\n";
    
    while(<FILE>){
	if(index($_, ">") != -1){
	    if($genome_size != 0){
		$file_size += $genome_size;
		if($file_size > $max_split_size || 
		   $genome_number > $max_genome_number){
		    #print " *** $file_size $split_file_num\n";<STDIN>;
		    $split_file_num++;
		    $file_size = $genome_size;
		    $genome_number = 0;
		    close($OUT);
		    open($OUT, ">$split_dir/split_$split_file_num.fa");
		}
		print $OUT $genome;
	    }
	    $genome = $_;
	    $genome_size = 0;
	    $genome_number++;
	}
	else{
	    $genome .= $_;
	    $genome_size += length($_);
	}
	#last if($split_file_num == 3);
    }
    print $OUT $genome;
}


sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";
    print STDERR `$exe` if($run);
}

#/mnt/software/unstowable/mummer-3.23/nucmer --nosimplify --maxmatch -p /mnt/ScratchPool/bertrandd/OPERA_MS/cow_rumen//ASSEMBLY//SOAP_SCAFFOLD/././200_300_3k_5k/mummer/split_17.fa_split_46.fa GENOME_REF/all_bacterial_genome_split/split_17.fa /mnt/ScratchPool/bertrandd/OPERA_MS/cow_rumen//ASSEMBLY//SOAP_SCAFFOLD/././200_300_3k_5k/split_contig/split_46.fa
