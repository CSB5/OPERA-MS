#!/usr/bin/perl
#use strict;
use warnings;
use Statistics::Basic qw(:all);

#####################################################
#Display information about the constructed scaffolds for downstream analysis.
#####################################################

#First we get the total number of contigs for each cluster
my $ref_directory = $ARGV[0];
my $nb_process = $ARGV[1];
my $opera_ms_dir = $ARGV[2];
my $scaffold_seq_file = $ARGV[3];
my $mummer_dir = $ARGV[4];
#
my $cluster_file = "$ref_directory/clusters_seq_similarity";
my $species_file = "$ref_directory/cluster_species.dat";
my $super_cluster_file = "$ref_directory/super_cluster.dat";
#my $cluster_file = "../sigma/clusters";
my $scaffold_file = "$ref_directory/../opera_long_read/scaffolds.scaf";
#my $scaffold_seq_file = "$ref_directory/../opera_long_read/scaffoldSeq.fasta.filled";
my $outfile = "$ref_directory/../../scaffold_info.txt";
my $contigs_windows_glob = "$ref_directory/../coverage_estimation/contigs_*";
my $clusters_cov_file = "$ref_directory/NO_REPEAT/repeat_contigs_info.txt";
my $strain_directory = "$ref_directory/../strain_analysis/";

my %clusters_to_contigs = ();
my %numb_contigs_cluster = ();
my %cluster_to_species = ();
my %super_to_subclust = ();
my %clusters_cov = ();
my %contig_cov = ();
my @cov_list = [];

my $use_nucmer = 0;

my @contig_window_glob = glob($contigs_windows_glob);
my $contigs_windows_file = $contig_window_glob[0];
open (FILE, $contigs_windows_file) or die;
my $header = <FILE>;
chomp $header;
my @line = split (/ /, $header);
my $window_size = $line[3];
while (<FILE>) {
    chomp $_;	
    #The line with the contig ID
    my @line = split (/\t/, $_);
    my $contig = $line[0];
    my $length = $line[1];
    my $nb_window = $line[4];
    #The line with number of arriving reads
    my $read_count = <FILE>;
    chop $read_count;
    #print STDERR $read_count."\t|".$nb_window."|\t".$window_size."\n";<STDIN>;
    #Skip the next line that contian the windows (need to compute the variance latter)
    my $str = <FILE>;chop $str;
    $contig_cov{$contig} = $read_count/($nb_window*$window_size);
}

open (SUPER_CLUSTER, $super_cluster_file);
print STDERR "mkdir -p $ref_directory/POST_EVAL/\n";
`rm -r $ref_directory/POST_EVAL/;mkdir -p $ref_directory/POST_EVAL/`;
while (<SUPER_CLUSTER>){
    chomp $_;
    my @line = split(/\t/, $_);
    $super_to_subclust{$line[0]} = [@line];
    #print STDERR @{$super_to_subclust{$line[0]}} . "\n";
}
close(SUPER_CLUSTER);

print STDERR " *** Reading cluster coverage file\n";
open (CLUSTERS_COV_FILE, $clusters_cov_file) or die;
while (<CLUSTERS_COV_FILE>){
    chomp $_;
    if ($_ =~ /^\t/){next;}
    my @line = split(/\s/, $_);
    my $cluster = $line[0];
    my $cov = $line[1];
    $clusters_cov{$cluster} = $cov;
    #print STDERR $cluster . "\n";
}
close(CLUSTERS_COV_FILE);

open (CLUSTERS, $cluster_file) or die;
while (<CLUSTERS>){
    chomp $_;
    my @line = split(/\t/, $_);
    my $contig = $line[0];
    my $cluster = $line[1];
    $clusters_to_contigs{$contig} = $cluster;
    
    if (exists $numb_contigs_cluster{$cluster}){
        $numb_contigs_cluster{$cluster}++;
    }

    else{
        $numb_contigs_cluster{$cluster} = 1;
    }
}
close(CLUSTERS);

print STDERR " *** Reading reference clustering file\n";
open (SPECIES, $species_file) or die;
my $cluster_name = "";
my $best_species = "";
while(<SPECIES>){
    chomp $_;
    if ($_ =~ />(.*)/){
        if ($cluster_name ne ""){
	    $cluster_to_species{$cluster_name} = $best_species if ($best_species ne ""); 
	    #print STDERR "$best_species, $cluster_name\n";
        }
        #print STDERR $1 . "\n";
        $cluster_name = $1;
        $best_species = "";
    }

    else{
        if ($best_species eq ""){
            $best_species = $_;
        }
    }
}

#Now we associate each scaffold with a cluster, and find out how many contigs
#from that scaffold are associated with the cluster
open (SCAFFOLDS, $scaffold_file) or die;
open (SCAFFOLDSEQ, $scaffold_seq_file) or die;
open (OUTFILE, '>', $outfile) or die;
#print OUTFILE "SCAFF_ID\tLENGTH\tSCAFF_COVERAGE\tCLUSTER_ID\tCONTIGS/CLUSTER_CONTIGS\tREFERENCE_GENOME\tCLUSTER_INFO\tREF_GENOME_LENGTH\tALIGNED_LENGTH\tFRACTION_GENOME_COVERED\tFRACTION_ALIGNED_TO_GENOME\tIDENTITY\n";
print OUTFILE "SEQ_ID\tLENGTH\tARRIVAL_RATE\tSPECIES\tNB_STRAIN\tREFERENCE_GENOME\tFRACTION_GENOME_COVERED\tIDENTITY\n";
my $scaf_name = "";
my $numb_contigs_in_scaffold = 0;
my $cluster_for_scaffold; 
my $cmd_file = "$ref_directory/POST_EVAL/cmd.txt";

my %filled_scaffold_length = ();

#Preprocessing for nucmer.

#print STDERR " *** Run nucmmer\n";
my $scaffold_name;
my $scaffold_seq;
open (CMD, ">",  $cmd_file) or die;
while (<SCAFFOLDS>){
    chomp $_;
    if ($_ =~ />/){
	my @line = split(/\t/, $_);
	$_ =~ /length:\s(\d*)\s/;
	my $length = $1;

	my $scaffold = $line[0];
	my $scaffname = substr $scaffold,1;

	$scaffold_name = <SCAFFOLDSEQ>;
	$scaffold_seq = <SCAFFOLDSEQ>;
	
	$length = length($scaffold_seq) - 1;
	$filled_scaffold_length{$scaffname} = $length;
	#
	#print STDERR " *** *** $scaffold_name $scaffname $length\n";#<STDIN>;
	
	#Make intermediate files for nucmer only for large enough scaffolds.
	if($use_nucmer && $length > 999){
	    
	    open (SCAF_OUT, ">", "$ref_directory/POST_EVAL/$scaffname") or die;
	    print SCAF_OUT $scaffold_name;
	    print SCAF_OUT $scaffold_seq;
	    close (SCAF_OUT);
	    my @contigline = split(/\t/, <SCAFFOLDS>);
	    my $contig = $contigline[0];
	    print STDERR $contig . "\n";
	    my $cluster_for_scaffold = $clusters_to_contigs{$contig};
	    if (defined $cluster_for_scaffold and
		exists $cluster_to_species{$cluster_for_scaffold}){
		my $command = "${mummer_dir}nucmer --maxmatch -c 400 --banded $ref_directory/POST_EVAL/$scaffname $opera_ms_dir/$cluster_to_species{$cluster_for_scaffold} -p $ref_directory/POST_EVAL/$scaffname; ${mummer_dir}show-coords -lrcT $ref_directory/POST_EVAL/$scaffname.delta > $ref_directory/POST_EVAL/$scaffname-out\n"; 
		print CMD $command;
	    }
	}
    }
}
close(SCAFFOLDS) or die;
close(SCAFFOLDSEQ) or die;
close(CMD);
if ($use_nucmer){
    #exit(0);
    #Go through and parse the scaffolds file and the nucmer.
    my $command = "cat $cmd_file | xargs -L 1 -P $nb_process -I COMMAND sh -c \"COMMAND\" 2> $cmd_file-log.txt";
    print STDERR $command . "\n";
    `$command`;
}

open (SCAFFOLDS, $scaffold_file) or die;
#my $no_species_info_line = "NA\tNA\tNA\tNA\tNA";
my $no_species_info_line = "NA\tNA\tNA";
my $species_info_line = "";
my ($ref_genome, $species_name);
my %nb_species_strain = ();
while(<SCAFFOLDS>){
    chomp $_;
    my @line = split(/\t/, $_);
    my $contig = $line[0];
    my $next_scaffold = $_;

    if ($_ =~ />/){
        if ($scaf_name ne ""){
            #$scaf_name =~ /length:\s(\d*)\s/;
	    my @scaff_delim = split(/\t/, $scaf_name);
            my $scaff_id = $scaff_delim[0];
	    my $scaffold_name = substr($scaff_id,1);
            my $median_cov = median(@cov_list);
            #print STDERR " *** $scaffold_name\n";<STDIN>;# @cov_list . "\n";
            $length= $filled_scaffold_length{$scaffold_name};
	    my $is_plasmid = 0;

            if (defined($cluster_for_scaffold)){
                print OUTFILE $scaff_id . "\t$length\t$median_cov\t";
            }
            else{
                print OUTFILE $scaff_id . "\t$length\tNA\t";
            }
	    #print STDERR " **** $scaf_name $length\n";
            if ($length > 999){
                
                if (exists $cluster_to_species{$cluster_for_scaffold}){
		    $species_info_line = "";
		    
                    #open (REF_GENOME, "$opera_ms_dir/$cluster_to_species{$cluster_for_scaffold}") or die;
		    #my $firstline = <REF_GENOME>;
                    #if ($firstline =~ /plasmid/){
			#$species_name = "Potential_plasmid";
                    #}
		    #else{
			#Get the species name
			$ref_genome = $cluster_to_species{$cluster_for_scaffold};
			@tmp = split(/\//, $ref_genome);
			@tmp_2 = split(/\_/,$tmp[@tmp-2]);
			$species_name = $tmp_2[0] . "_" . $tmp_2[1];
			if(! exists $nb_species_strain{$species_name}){
			    $nb_s = 1;
			    $s_s_dir = "$strain_directory/$species_name";
			    if(-d $s_s_dir){
				#print STDERR "ll $s_s_dir/STRAIN_*/contigs.fa/n";
				$nb_s = `ls -l $s_s_dir/STRAIN_*/contigs.fa | wc -l`;chop $nb_s;
			    }
			    $nb_species_strain{$species_name} = $nb_s;
			}
		    #}
                    #print OUTFILE $cluster_to_species{$cluster_for_scaffold} . "\t";
		    
		    #my @scaff_line = split(/\t/, $scaf_name);
		    #$n_info = get_nucmer_info($length,$scaffold_name);
		    #print OUTFILE $species_name . "\t" . $nb_species_strain{$species_name} . "\t" . "$opera_ms_dir/$ref_genome" . "\t" . $n_info . "\t";
		    print OUTFILE $species_name . "\t" . $nb_species_strain{$species_name} . "\t" . "$opera_ms_dir/$ref_genome";
                }

                else{
		    print OUTFILE $no_species_info_line;
                }
            }
	    
            else{
                print OUTFILE $no_species_info_line;
            }
	    print OUTFILE "\n";
            
        }

        $numb_contigs_in_scaffold = 0; 
        $scaf_name = $next_scaffold;
        $cluster_for_scaffold = undef;
        undef(@cov_list);
    }

    else{
        #if contig is less than 500 bp, it will not exist in the clusters file.
        if (!exists($clusters_to_contigs{$contig})) {  next; }
        $numb_contigs_in_scaffold++;
        $cluster_for_scaffold = $clusters_to_contigs{$contig};
        push @cov_list, $contig_cov{$contig};
    }

}

close(SCAFFOLDS);

sub get_species_metrics{
    my ($cluster, $species, $dir) = @_; 
    my $mashfile = "$dir/MASH/$cluster.dat";
    my $metrics = "";
    open (MASHFILE, $mashfile) or die("$cluster, $species, $dir"); 

    while (<MASHFILE>){
        chomp $_;
        my @line = split(/\t/, $_);
        if ($line[0] eq $species){
            $metrics = "C-$cluster;$line[2];$line[4]";
        }
    }

    close (MASHFILE);

    return $metrics;
}

sub get_nucmer_info{
    my ($length, $scaffold) = @_;
    #print STDERR "$ref_directory/POST_EVAL/$scaffold-out\n";
    my $aligned_length;
    my $ref_length;
    my $percent_scaff_align;
    my $identity_weighted;
    #print STDERR " *** Read file $ref_directory/POST_EVAL/$scaffold-out\n";
    open (NUCMER_INFO, "$ref_directory/POST_EVAL/$scaffold-out") or die " *** File $scaffold-out not found\n";
    <NUCMER_INFO>;
    <NUCMER_INFO>;
    <NUCMER_INFO>;
    <NUCMER_INFO>;

    while(<NUCMER_INFO>){
        chomp $_;
        my @line = split(/\t/, $_);
        $ref_length = $line[8];
        $aligned_length += $line[4];
        $identity_weighted += $line[6] * $line[4]; 
    }

    if(defined $aligned_length){
        $percent_scaff_align = $aligned_length / $length;
    }

    else{
        return "NA\tNA";
    }

    $identity_weighted = $identity_weighted / $aligned_length;

    my $percent_ref_covered = $aligned_length/$ref_length;
    return "$percent_scaff_align\t$identity_weighted";
    #return "$ref_length\t$aligned_length\t$percent_ref_covered\t$percent_scaff_align\t$identity_weighted\t";
}
