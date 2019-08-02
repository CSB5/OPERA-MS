#!/usr/bin/perl
use Statistics::Basic qw(:all);
use warnings;
use strict;


my ($inter_dir, $out_dir, $contigs_file, $reference_mapping_dir, $opera_ms_dir, $FLAG_USE_REF) = @ARGV;
#my $out_dir = $ARGV[1];
#my $contigs_file = $ARGV[2];

my %cluster_lengths = ();
my %cluster_to_analyze = ();
my %species_to_ref_genome = ();
my %species_to_contigs_length = ();
my %species_to_clusters = ();
my %species_to_analyze = ();
my %contig_info = ();

my @windows_file = glob("$inter_dir/coverage_estimation/contigs_*");
my $windows_file = $windows_file[0];
#my $clusters_file = "$inter_dir/reference_mapping/NO_REPEAT/repeat_contigs_info.txt";
my $clusters_species_file = "$reference_mapping_dir/cluster_species.dat";
#
#Identify the clusters used 
my $ref_clusters = "$reference_mapping_dir/clusters_seq_similarity";
$ref_clusters = "$inter_dir/sigma/clusters" if(! $FLAG_USE_REF);

#Get length and coverage of contigs
open (FILE, $windows_file) or die;
my $header = <FILE>;
chop $header;
my @line = split (/ /, $header);
my $window_size = $line[3];
while (<FILE>) {
    chomp $_;   
    my @line = split (/\t/, $_);
    my $contig = $line[0];
    my $length = $line[1];
    my $nb_window = $line[4];
    my $read_count = <FILE>;chop $read_count;
    my $str = <FILE>;chop $str;

    if ($length >= 500){
        $contig_info{$contig} = {"COV", $read_count/($nb_window*$window_size), "LENGTH", $length};;
    }
}
close(FILE);


#Get cluster length
open(REF_CLUS, $ref_clusters) or die;
while(<REF_CLUS>){
    chomp $_;
    @line = split(/\t/,$_); 
    my $contig = $line[0];
    my $cluster = $line[1];

    $cluster_lengths{$cluster} =  0 if(! exists $cluster_lengths{$cluster});
    
    $cluster_lengths{$cluster} += $contig_info{$contig}->{"LENGTH"};
}
close(REF_CLUS);

my $cluster;
my $seq_length;
my %strain_contig = ();

if($FLAG_USE_REF){
    open (STRAINS, $clusters_species_file) or die;
    #Identifiy the reference associated to each species using a vote based method using the most voted reference in each cluster
    #SHOULD WE:
    # * weight based on cluster length
    # * MASH the cluster again
    #Collect the vote and summup the cluster length that belong to the same species
    while (<STRAINS>){
	chomp $_;
	if($_ =~ />(.*)/){
	    $cluster = $1;
	    chomp $cluster;
	    #operate at a species level
	    my $strain = <STRAINS>;
	    $strain = "$opera_ms_dir/$strain";
	    chomp $strain;
	    while($strain =~ />(.*)/){
		$cluster = $1;
		chomp $cluster;
		$strain = <STRAINS>;
		if(defined $strain){
		    chomp $strain;
		    $strain = "$opera_ms_dir/$strain";
		    #print STDERR " *** *** $strain\n";
		}
		else{
		    $strain = "";
		}
	    }

	    if($strain eq ""){
		next;
	    }

	    my @full_path = split(/\//, $strain);
	    my $strain_name = $full_path[@full_path-2];
	    my @strain_delim = split(/_/,$strain_name);
	    my $species = $strain_delim[0] . "_" . $strain_delim[1];
	    
	    #print STDERR " *** $strain $species\n";<STDIN>;
	    
	    #$cluster = $1;
	    #chomp $cluster;
	    #print STDERR $cluster . "\n";
	    
	    if(!exists $cluster_lengths{$cluster}){
		next;
	    }
	    my $seq;
	    my $is_plasmid = 0;
	    if(index($strain, "XXX") == -1){
		print STDERR " *** Open file $strain\n";
		open(CHECK_REF, $strain) or die;
		my $header = <CHECK_REF>;
		
		close(CHECK_REF);

		if($header =~ /plasmid/){
		    $is_plasmid = 1;
		    #print STDERR "PLASMID DETECTED $strain\n";
		}
	    }

	    if(!$is_plasmid){

		if(!exists $species_to_ref_genome{$species}){
		    $species_to_ref_genome{$species} -> {$strain}++;
		}
		else{
		    $species_to_ref_genome{$species} -> {$strain} = 1;
		}

	    }

	    if(!exists $species_to_contigs_length{$species}){
		$species_to_contigs_length{$species} = $cluster_lengths{$cluster};
	    }
	    else{
		$species_to_contigs_length{$species} += $cluster_lengths{$cluster};
	    }

	    $species_to_clusters{$species} -> {$cluster} = 1;

	}

	else{
	    next;
	}

	#print STDERR $cluster . "\n";
	#print STDERR $seq_length . "\n";
	#print STDERR $cluster_lengths{$cluster} . "\n";
    }
    close(STRAINS);

    #Select the best refrence genome for a species
    open(OUT, ">$out_dir/reference_length.dat");
    my ($ref_genome_best, $best_count);
    foreach my $species (keys %species_to_ref_genome){
	$ref_genome_best = "NA";
	$best_count = 0;
	foreach my $ref_genome (keys %{$species_to_ref_genome{$species}}){
	    next if(index($ref_genome, "XXX") != -1);#This is a reference obtain using kraken
	    if (($species_to_ref_genome{$species}->{$ref_genome}) > $best_count){
		$ref_genome_best = $ref_genome;
		$best_count = ($species_to_ref_genome{$species}->{$ref_genome});
	    }
	}
	############### NEED TO UPDATE
	$seq_length = 3500000;
	if($ref_genome_best ne "NA"){
	    my $seq;
	    #print STDERR " *** $opera_ms_dir $ref_genome_best\n";
	    open (REF_GENOME, "$ref_genome_best") or die;
	    <REF_GENOME>;
	    while(<REF_GENOME>){
		chomp $_;
		$seq .= $_;
	    }
	    close(REF_GENOME);
	    $seq_length = length($seq);
	}
	#####################
	
	#print STDERR " *** length -> $seq_length -> $species_to_contigs_length{$species}\n";
	if(#$species_to_contigs_length{$species} > $seq_length + ($seq_length * 0.1) and
	   $species_to_contigs_length{$species} > 1000000
	   #$seq_length > 1000000
	    ){
	    print OUT $species . "\t" . $seq_length . "\t" . $species_to_contigs_length{$species} . "\t" . $ref_genome_best . "\n";
	    $species_to_analyze{$species} = 1;
	}
	else{
	    print STDERR $species . "\t" . $seq_length . "\t" . $species_to_contigs_length{$species} . "\t" . $ref_genome_best . "\n";
	}
    }
    close(OUT);

    #
    #exit();
    #
}
else{
    #Indentify cluster with potnetially multiple strains based on the cluster size
    #mkdir reference_mapping strain_analysis NEED to create those directories need to change that: 1) write all related file in strin_directory and not in reference_mapping 2) create a NO_REF_strain_analysis directory
    my $c_length = 0;my $species;
    my $NO_REF_SIZE_THRESHOLD = 2500000;#This is potentially a completegenome that can contains contigs from another strains
    open(OUT, ">$out_dir/reference_length.dat");
    foreach $cluster (keys %cluster_lengths){
	$c_length = $cluster_lengths{$cluster};
	if($c_length > $NO_REF_SIZE_THRESHOLD){
	    $species = "S_$cluster";
	    $species_to_clusters{$species} -> {$cluster} = 1;
	    print OUT  $species . "\t" . $c_length . "\t" . $c_length . "\t" . "XXX" . "\n";
	    $species_to_analyze{$species} = 1;
	}
    }
}

#Construct the strain level directory
#And write the matrix file: contig to coverage information
foreach my $species (keys %species_to_analyze){
    run_exe("mkdir -p $out_dir/$species");
    
    open(REF_CLUS, $ref_clusters) or die("File $ref_clusters not found\n");
    open (MATRIX, ">$out_dir/$species/matrix");
    while(<REF_CLUS>){
	chomp $_;
	my @line = split(/\t/,$_); 
	my $contig = $line[0];
	my $cluster = $line[1];
	if (exists ($species_to_clusters{$species} -> {$cluster})){
	    print MATRIX  $contig . "\t" . $contig_info{$contig}->{"COV"} . "\n";
	    $strain_contig{$contig} = 1;
	}
    }
    close(MATRIX);
    close(REF_CLUS);
    
    #Extract the edges and contigs from that cluster
    run_exe("${opera_ms_dir}utils/perl $opera_ms_dir/bin/get_edges_from_matrix.pl $inter_dir $out_dir/$species/matrix $out_dir/$species $contigs_file 0 1000 2> $out_dir/$species/get_edge_matrix.log");
    if($?){
	die "Error during get_edges_from_matrix.pl. Please see $out_dir/$species/get_edge_matrix.err for details.\n";
    }
}

    
#Read the orginal cluster file and remove any contigs that belong to species with multiple strain
open(REF_CLUS, $ref_clusters) or die("File $ref_clusters not found\n");
open(OUT, ">$reference_mapping_dir/clusters_single_strain");
while(<REF_CLUS>){
    chomp $_;
    my @line = split(/\t/,$_); 
    my $contig = $line[0];
    my $cluster = $line[1];
    if(! exists $strain_contig{$contig}){
	print OUT $_ . "\n";
    }
}
close(OUT);
close(REF_CLUS);


################3
#print STDERR " **** Exit before creating single_strain_species.fa file\n";
#exit();
####################

#Read the original contig sequence file and remove any contigs that belong to species with multiple strain
my $prev_scaf;
my $seq = "";
my $contigs;
open(FILE, $contigs_file) or die("File $contigs_file not found\n");
open(CONTIGS, ">$reference_mapping_dir/single_strain_species.fa");
while (<FILE>){
    chomp $_;
    if ($_ =~ />/){
        if(defined $prev_scaf){
            my @line = split(/\s/, $prev_scaf);
            $contigs = substr $line[0], 1;

            if(! exists $strain_contig{$contigs}){
                print CONTIGS $prev_scaf . "\n";
                print CONTIGS $seq . "\n";
            }
        }

        $prev_scaf = $_;
        $seq = "";
    }
    else{
        $seq .= $_;
    }
}
close(FILE);

if(! exists $strain_contig{$contigs}){
    print CONTIGS $prev_scaf . "\n";
    print CONTIGS $seq . "\n";
}
close(CONTIGS);


sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";#<STDIN>;
    print STDERR `$exe` if($run);
}

