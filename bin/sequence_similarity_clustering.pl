#!/usr/bin/perl
use warnings;
use Statistics::Basic qw(:all);

#use strict;


#############################################################
#Script for merging sigma clusters based on reference genomes.
#1. Gather all clusters with connections between them
#2. For all clusters with edges between them, check if they have similar reference genomes.
#3. Recover edges between clusters in within the super clusters.
#############################################################
my ($inter_dir, $ref_map_folder, $mapping_folder, $sigma_folder, $contig_seq_file, $nb_process, $opera_ms_db, $opera_ms_dir, $mash_exe_dir, $mummer_dir, $kraken_exe_dir) = @ARGV;

require "$opera_ms_dir/bin/test_time.pl";

my $PARSE_FILES_MASH = 1;
my $PARSE_FILES_NUC = 1;
my $MERGE_CLUSTERS = 1;

my %reference_mappings_ = ();
my %contig_info_ = ();
my %cluster_info_ = ();
my %edge_info = ();
my %cluster_species_set_ = ();
my %super_cluster_species_set_ = ();
my %species_to_ignore_ = ();
my @super_cluster = ();
my %cluster_to_super_cluster_ID;
my %super_cluster_strain_set_ = ();


#my $ref_map_folder = "$inter_dir/reference_mapping";
my $cluster_results_file = "$ref_map_folder/clusters_seq_similarity";
my $cluster_species_file = "$ref_map_folder/cluster_species.dat";
my $super_cluster_file = "$ref_map_folder/super_cluster.dat";
my $edges_between_clusters_file = "$ref_map_folder/edges_between_clusters";
my $good_edges_file = "$ref_map_folder/edges_between_clusters_good";

my $clusters_file = "$sigma_folder/clusters";
my $paired_edges_file = "$mapping_folder/pairedEdges";
my $paired_edges_file_above_thresh = "$mapping_folder/pairedEdges_above_thresh";
#
#Need to update to multi sample assemblies
my $short_edges_file = `ls $inter_dir/lib_1_bundles/clusters_*`;chop $short_edges_file;#"$inter_dir/lib_1_bundles/clusters_opera";
#
my $mash_dir = "$ref_map_folder/MASH/";
my $nucmer_dir = "$ref_map_folder/NUCMER_OUT/"; 

######################################START########################3

my $command;
run_exe("mkdir -p $ref_map_folder") if (!(-d $ref_map_folder));

my ($start_time, $end_time);
get_cluster_contigs_info(\%contig_info_, \%cluster_info_);

$start_time = time;
if ($PARSE_FILES_MASH){
    my %edges_to_cov_diff = ();
    run_mash_on_clusters(\%cluster_info_, \%contig_info_, "$opera_ms_db/genomes.msh");
    generate_edges_between_clusters(\%edges_to_cov_diff, \%contig_info_);
}
$end_time = time;
write_time("$inter_dir", "r_mash", ($end_time - $start_time));
#print STDERR "***  Mash mapping Elapsed time: " . ($end_time - $start_time) . "s\n";#<STDIN>;

generate_clusters_species(\%cluster_species_set_, \%species_to_ignore_, \%cluster_info_);

$start_time = time;
if ($PARSE_FILES_NUC){
    run_nucmer(\%cluster_species_set_, \%contig_info_, \%reference_mappings_, \%cluster_info_, "$opera_ms_db/..");
    get_correct_edges(\%reference_mappings_, \%contig_info_);
}
$end_time = time;
write_time("$inter_dir", "r_mummer", ($end_time - $start_time));
#print STDERR "***  Mummer mapping Elapsed time: " . ($end_time - $start_time) . "s\n";#<STDIN>;

$start_time = time;
if ($MERGE_CLUSTERS){
    merge_clusters(\%contig_info_, \%cluster_species_set_, \%super_cluster_species_set_, \%super_cluster_strain_set_);
    recover_edges(\%contig_info_);
    #
    #ADD annotation from kraken NEED TO UPDATE HE WHOLE CODE TO CHECK FOR KRAKEN DEPENDENCY
    if(0 && $kraken_ref ne "NULL"){
	run_kraken_annotation($contig_seq_file);
	add_kraken_annotation(\%contig_info_, \%cluster_species_set_, "$ref_map_folder/KRAKEN_OUT/kraken.out.annot");
    }
    #
    export_info(\%contig_info_, \%cluster_species_set_, \%super_cluster_strain_set_);
}
$end_time = time;
write_time("$inter_dir", "r_cluster_merging", ($end_time - $start_time));
#print STDERR "***  Cluster merging Elapsed time: " . ($end_time - $start_time) . "s\n";#<STDIN>;
###################################END################################

sub run_kraken_annotation{
    my ($contig_seq_file) = @_;
    run_exe("mkdir $ref_map_folder/KRAKEN_OUT");
    $kraken_out_file = "$ref_map_folder/KRAKEN_OUT/kraken.out";
    #
    run_exe("$kraken_exe_dir/kraken --threads $nb_process -db $kraken_ref --fasta-input --output $kraken_out_file $contig_seq_file 2> $kraken_out_file.log");
    run_exe("$kraken_exe_dir/kraken-report --db $kraken_ref $kraken_out_file > $kraken_out_file.report 2> $kraken_out_file.report.log"); #<STDIN>;
    run_exe("$opera_ms_dir/bin/kraken_bin_file.py -k $kraken_out_file -r $kraken_out_file.report -o $kraken_out_file.annot 2> $kraken_out_file.annot.log");#<STDIN>
}

sub add_kraken_annotation{
    my ($contig_info, $cluster_species_set, $kraken_result_file) = @_;

    open(FILE, $kraken_result_file);
    while(<FILE>){
	($contig, $species) = split(/\s+/, $_);
	$contig_info{$contig}->{"KRAKEN"} = $species;
    }
    close(FILE);
    open(OUT_K, ">$kraken_result_file.log");
    my %annot_cluster = ();
    foreach my $contig (keys %{$contig_info}){
	$cluster = $contig_info->{$contig}->{"CLUSTER"};
	#
	next if(! defined $cluster);
	
	if(! exists $annot_cluster{$cluster}){
	    $annot_cluster{$cluster} = 1;
	    $annot_cluster{$cluster} = 0 if(@{$cluster_species_set->{$cluster}} + 0 == 0);
	}
	
	if(! $annot_cluster{$cluster}){
	    push(@{$cluster_species_set->{$cluster}}, "XXX/" . $contig_info{$contig}->{"KRAKEN"} . "/XXX") if(exists $contig_info{$contig}->{"KRAKEN"});
	}
	else{
	    print OUT_K $contig . "\t" . $contig_info{$contig}->{"KRAKEN"} . "\t" . $cluster . "\t" . join(";", @{$cluster_species_set->{$cluster}}) . "\n";
	}
    }
    close(OUT_K);
}


#Export some info about the clusters. The clusters file, the species file, and a super cluster file.
sub export_info{
    my ($contig_info,$cluster_species_set, $super_cluster_strain_set) = @_;
    open (OUT, ">", $cluster_results_file) or die;
    #print STDOUT "EXPORTING INFO\n";
    my $cluster;
    foreach my $contig (keys %{$contig_info}){
        $cluster = $contig_info->{$contig}->{"CLUSTER"};
        next if (!defined($cluster));
        if(exists $cluster_to_super_cluster_ID{$cluster}){
            print OUT $contig . "\t" . $cluster_to_super_cluster_ID{$cluster}."S". "\t" . "001100" . "\t" . "0.0088888" . "\n";
        }
        else{
            print OUT $contig . "\t" . $cluster . "\t" . "001100" . "\t" . $contig_info->{$contig}->{"CLUSTER_MEAN"} . "\n";
        }
	#print STDOUT $contig . "\n";
    }
    close(OUT);
    
    open (OUT, ">", $super_cluster_file) or die;
    for(my $i = 0; $i < @super_cluster; $i++){
	my $super_c = $super_cluster[$i];
        if ($super_c eq "NAN"){
            #print STDERR "Super cluster $i is 'NAN'";
            next;
        }
        print OUT "${i}S" . "\t" . (join("\t", @{$super_c})) . "\n" if($super_c ne "NAN");
    }
    close(OUT);

    open (OUT, ">", $cluster_species_file) or die;
    foreach my $id (keys %{$super_cluster_strain_set}){
        print OUT ">${id}S\n";
        foreach my $strain (sort {$super_cluster_strain_set->{$id}->{$b} <=> $super_cluster_strain_set->{$id}->{$a}} keys %{$super_cluster_strain_set->{$id}}){
            print OUT "$strain\n";
        }
    }

    foreach my $id (keys %{$cluster_species_set}){
        print OUT ">$id\n";
        foreach my $species (@{$cluster_species_set->{$id}}){
            print OUT "$species\n";
        }
    }
    close(OUT);
}

#Retrieve edges inside of each other as well as rescued edges between clusters inside of the super clusters.
sub recover_edges{
#first recover edges between clusters
    my ($contig_info) = @_;
    open(GOOD, $good_edges_file) or die;
    my @outfiles = ();

    for(my $i = 0; $i < 6; $i++) {
        my $OUT;
        open ($OUT, ">","$ref_map_folder/filtered_pairedEdges_i$i") or die;
        $outfiles[$i] = $OUT;
    }

    #Add edges between cluster. Those edge can be below the support threshold
    my @line;
    my ($contig1, $contig2, $dist, $cluster1, $cluster2);
    while (<GOOD>){
        chomp $_;
	@line = split(/\t/, $_);
	$contig1 = $line[0];
	$contig2 = $line[2];
	$dist = $line[4];
	$cluster1 = $contig_info->{$contig1}->{"CLUSTER"};
	$cluster2 = $contig_info->{$contig2}->{"CLUSTER"};

        if (!exists $cluster_to_super_cluster_ID{$cluster1} or
            !exists $cluster_to_super_cluster_ID{$cluster2}){
            next;
        }

        #print STDERR $cluster_to_super_cluster_ID{$cluster2} . "\t" . $cluster_to_super_cluster_ID{$cluster1} . "\n";
        my $super1 = $cluster_to_super_cluster_ID{$cluster1};
        my $super2 = $cluster_to_super_cluster_ID{$cluster2};

        if($super1 eq $super2){
            my $OUT = "NaN";
            if ($dist < 300){
                $OUT = $outfiles[0];
            }

            elsif ($dist < 1000){
                $OUT = $outfiles[1];
            }
            
            elsif ($dist < 2000){
                $OUT = $outfiles[2];
            }

            elsif ($dist < 5000){
                $OUT = $outfiles[3];
            }

            elsif ($dist < 15000){
                $OUT = $outfiles[4];
            }

            elsif ($dist < 40000){
                $OUT = $outfiles[5];
            }

	    next if($OUT eq "NaN");#NEED TO CHANGE TO TAKE INTO ACCOUNT READ LONGER THAN 40kb
	    
            print $OUT $_ . "\n";
        }
    }
    close(GOOD);

    #Add edges above the threshold, those edge can only be internal to a cluster identified by sigma
    open (FILE, "$paired_edges_file_above_thresh") or die("File $paired_edges_file_above_thresh does not exist\n");
    
    while (<FILE>){
        chomp $_;
        @line = split(/\t/, $_);
        my $contig1 = $line[0];
        my $contig2 = $line[2];
        my $supp = $line[6];
        my $dist = $line[4];
        my $cluster1 = $contig_info->{$contig1}->{"CLUSTER"};
        my $cluster2 = $contig_info->{$contig2}->{"CLUSTER"};

        if($cluster1 eq $cluster2 and $supp > 1){
            my $OUT;
            if ($dist < 300){
                $OUT = $outfiles[0];
            }

            elsif ($dist < 1000){
                $OUT = $outfiles[1];
            }
            
            elsif ($dist < 2000){
                $OUT = $outfiles[2];
            }

            elsif ($dist < 5000){
                $OUT = $outfiles[3];
            }

            elsif ($dist < 15000){
                $OUT = $outfiles[4];
            }

            elsif ($dist < 40000){
                $OUT = $outfiles[5];
            }

            print $OUT $_ . "\n";
        }
    }

    close(FILE);
}

#Merge clusters based on the graph generated in the previous step.
sub merge_clusters{
    print STDERR "***Merging Clusters***\n";
    my ($contig_info, $cluster_species_set, $super_cluster_species_set, $super_cluster_strain_set) = @_;
    open(FILE, "$good_edges_file") or die;
    my ($cluster_1, $cluster_2, $super1, $super2);
    while(<FILE>){
        @line = split(/\t/, $_);
        $cluster_1 = $contig_info->{$line[0]}->{"CLUSTER"};
        $cluster_2 = $contig_info->{$line[2]}->{"CLUSTER"};
        
	$super1 = "";
	$super2 = "";
	
        if (exists $cluster_to_super_cluster_ID{$cluster_1}){
            $super1 = $cluster_to_super_cluster_ID{$cluster_1};
        }

        if (exists $cluster_to_super_cluster_ID{$cluster_2}){
            $super2 = $cluster_to_super_cluster_ID{$cluster_2};
        }

        print STDERR " *** S(${super1}S,${super2}S) " . $cluster_1 . "\t" . $cluster_2 ."\n";

        #Check if two clusters share similar reference genomes and merge.
        if(common_species($cluster_1, $cluster_2, $cluster_species_set, $super_cluster_species_set, \%cluster_to_super_cluster_ID)){
            print STDERR " *** *** *** Clusters to merge $cluster_1 $cluster_2\n";
            #2 clusters not in super cluster
            if(! exists $cluster_to_super_cluster_ID{$cluster_1} && ! exists $cluster_to_super_cluster_ID{$cluster_2}){
                my $super_cluster_ID = @super_cluster+0;
                $cluster_to_super_cluster_ID{$cluster_1} = $super_cluster_ID;
                $cluster_to_super_cluster_ID{$cluster_2} = $super_cluster_ID;
                $super_cluster[$super_cluster_ID] = [$cluster_1, $cluster_2];
                
                my %species_intersection = ();
                undef(%species_intersection); 
                my %strain_set = ();
                undef(%strain_set);
                
                foreach my $s1 ( @{$cluster_species_set->{$cluster_1}}){
		    my $species1 = get_species_name($s1);

		    foreach my $s2 (@{$cluster_species_set->{$cluster_2}}){
			my $species2 = get_species_name($s2);

                        if($species2 eq $species1){
                            if (! exists $species_intersection{$species1}){
                                $species_intersection{$species1} = 1;
                            }
                        }
                    }
                }

                foreach my $s ( @{$cluster_species_set->{$cluster_1}}){
                    $strain_set{$s} = 1;
                }

                foreach my $s2 (@{$cluster_species_set->{$cluster_2}}){

                    if(exists $strain_set{$s2}){
                        $strain_set{$s2}++;
                    }

                    else{
                        $strain_set{$s2} = 1;
                    }

                }
                #my $length = @species_intersection;
                #print STDERR "INTERSECTION @species_intersection, LENGTH $length\n";
                my $strain_to_map_to = $cluster_species_set->{$cluster_1}->[0];
                $super_cluster_strain_set->{$super_cluster_ID} = {%strain_set};
                $super_cluster_species_set->{$super_cluster_ID} = [keys %species_intersection]; 
            }
            else{
            #1 clusters not in super cluster
                if(! exists $cluster_to_super_cluster_ID{$cluster_1} || ! exists $cluster_to_super_cluster_ID{$cluster_2}){
                    if(exists $cluster_to_super_cluster_ID{$cluster_2}){
                        my $tmp = $cluster_2;
                        $cluster_2 = $cluster_1;
                        $cluster_1 = $tmp;
                    }
                    my $super_cluster_ID = $cluster_to_super_cluster_ID{$cluster_1};
                    $cluster_to_super_cluster_ID{$cluster_2} =  $super_cluster_ID;
                    push(@{$super_cluster[$super_cluster_ID]}, $cluster_2);
                    #needed to clear out the array
                    my %strain_set = %{$super_cluster_strain_set->{$super_cluster_ID}};
                    my %species_intersection = ();
                    undef(%species_intersection); 

                    foreach my $s (@{$super_cluster_species_set->{$super_cluster_ID}}){
                        foreach my $s2 (@{$cluster_species_set->{$cluster_2}}){
			    
			    my $species = get_species_name($s2);
				
			    if($species eq $s){
				if (! exists $species_intersection{$species}){
				    $species_intersection{$species} = 1;
				}
                            }
                        }
                    }

                    foreach my $s2 (@{$cluster_species_set->{$cluster_2}}){

                        if(exists $strain_set{$s2}){
                            $super_cluster_strain_set->{$super_cluster_ID}->{$s2}++;
                        }

                        else{
                            $super_cluster_strain_set->{$super_cluster_ID}->{$s2} = 1;
                        }

                    }
                    
                    $super_cluster_species_set->{$super_cluster_ID} = [keys %species_intersection]; 
                    #my $length = @species_intersection;
                    #print STDERR "INTERSECTION @species_intersection, LENGTH $length\n";
                }
                #2clusters in super cluster => need to merge the super clusters
                else{
                    my $super_cluster_ID_1 = $cluster_to_super_cluster_ID{$cluster_1};
                    my $super_cluster_ID_2 = $cluster_to_super_cluster_ID{$cluster_2};
                    if($super_cluster_ID_1 != $super_cluster_ID_2){
                        print STDERR " *** Merge super cluster !!!\n";
                        print STDERR " *** super cluster size ". (@{$super_cluster[$super_cluster_ID_1]}) . " " . (@{$super_cluster[$super_cluster_ID_2]}) ."\n";
                        foreach my $clust (@{$super_cluster[$super_cluster_ID_2]}){
                            #print STDERR " *** $clust " . (@{$super_cluster[$super_cluster_ID_2]}) . "\n";<STDIN>;
                            push(@{$super_cluster[$super_cluster_ID_1]}, $clust);
                            $cluster_to_super_cluster_ID{$clust} = $super_cluster_ID_1;
                        }
                        $super_cluster[$super_cluster_ID_2] = "NAN";

                        my %species_intersection = ();
                        #Needed to clear out the array... some scope/memory issue
                        undef(%species_intersection); 

                        foreach my $s (@{$super_cluster_species_set->{$super_cluster_ID_1}}){
                            foreach my $s2 (@{$super_cluster_species_set->{$super_cluster_ID_2}}){
                                if ($s2 eq $s){
                                    $species_intersection{$s} = 1;
                                }
                            }
                        }
                        #my $length = @species_intersection;
                        #print STDERR "INTERSECTION @species_intersection, LENGTH $length\n";
                        my $strain_set1 = $super_cluster_strain_set->{$super_cluster_ID_1};
                        my $strain_set2 = $super_cluster_strain_set->{$super_cluster_ID_2};

                        foreach my $strain2 (keys %{$strain_set2}){
                            if (exists $strain_set1->{$strain2}){
                                $strain_set1->{$strain2} += $strain_set2->{$strain2};
                            }

                            else{
                                $strain_set1->{$strain2} = $strain_set2->{$strain2};
                            }
                        }

                        $super_cluster_species_set->{$super_cluster_ID_1} = [keys %species_intersection]; 
                        delete $super_cluster_species_set->{$super_cluster_ID_2};
                        delete $super_cluster_strain_set -> {$super_cluster_ID_2};
                    }
                }
            }
        }
    }
    close(FILE);
}

sub common_species{
    my ($cluster_1, $cluster_2, $cluster_species_set, $super_cluster_species_set, $cluster_to_super_cluster_ID) = @_;
    my @common_species = ();#NEVER USED !!!!
    #print STDERR "%$cluster_species_set\n";

    #TEST PROBABLY NOT REQUIRED
    if (!defined $cluster_species_set->{$cluster_1} or !defined $cluster_species_set->{$cluster_2}){
        return;
    }

    #Why copy the arrays ?? point is enought
    my @species1 = @{$cluster_species_set->{$cluster_1}};
    my @species2 = @{$cluster_species_set->{$cluster_2}};
    my $species1_super = 0;
    my $species2_super = 0;

    if (exists $cluster_to_super_cluster_ID->{$cluster_1}){
	my $super_cluster_1 = $cluster_to_super_cluster_ID->{$cluster_1}; 
        @species1 = @{$super_cluster_species_set->{$super_cluster_1}};
        $species1_super = 1;
    }

    if (exists $cluster_to_super_cluster_ID->{$cluster_2}){
       my $super_cluster_2  = $cluster_to_super_cluster_ID->{$cluster_2}; 
       @species2 = @{$super_cluster_species_set->{$super_cluster_2}};
       $species2_super = 1;
    }

    my $res = 0;
    #Species in common
    foreach my $s (@species1){ 
        my $species;

        if(!$species1_super){
	    #To extract the species name from the mash output
            $species = get_species_name($s);
        }

        else{
            $species = $s;
        }

        foreach my $s2 (@species2){
            my $species2;

            if(!$species2_super){
		$species2 = get_species_name($s2);
	    }

            else{
                $species2 = $s2;
            }

            if($species2 eq $species){
                $res = 1;
                push @common_species, $s;
                #print STDERR "$s2, $s\n";
            }
        }
    }

    return $res;
}
    
#Actually verify the edges between clusters. 
#Output a file that is the "true" cluster graph, with verified edges.
sub get_correct_edges{
    my ($reference_mappings, $contig_info) = @_;
    open (EDGES, $edges_between_clusters_file) or die;
    open (GOOD_EDGES, ">", $good_edges_file) or die;
    while(<EDGES>){
        chomp $_;
        @line = split(/\t/, $_);
        my $edge = $_;
        my $contig1 = $line[0];
        my $contig2 = $line[2];
        my $support = $line[6];
        my $stdev = $line[5];
        my $overlap = $line[4];

        #Short edges can be mapped to contigs that are smaller than 500 bp,
        #which are not included in this collection.
        if(!exists($contig_info->{$contig1}->{"CLUSTER"}) or
            !exists($contig_info->{$contig2}->{"CLUSTER"})){
                next;
            }

	#Always checked even for edges for support > 1
        if(check_reference_position($contig1, $contig2, $overlap, $stdev, $reference_mappings, $contig_info)){
            print GOOD_EDGES $edge . "\n";
        }
    }

    close(GOOD_EDGES);
    close(EDGES);
}

#Run nucmer on each cluster with an edge between them. We use
#the nucmer mapping to confirm if an edge is verified by the reference to be
#true.
sub run_nucmer{
    my ($cluster_species_set, $contig_info, $reference_mappings, $cluster_info, $genome_db) = @_;
    
    my %cluster_connection = ();
    run_exe("rm -r $nucmer_dir; mkdir -p $nucmer_dir;mkdir -p $nucmer_dir/temp_genome"); 
    open (my $CMD, ">", "$nucmer_dir/cmd.txt") or die;
    
    foreach my $cluster (keys %{$cluster_species_set}){
        if (@{$cluster_species_set->{$cluster}} != 0){
            my $top_species = @{$cluster_species_set->{$cluster}}[0];
	}
    }

    open (EDGES, $edges_between_clusters_file) or die;
    while(<EDGES>){
        chomp $_;
        @line = split(/\t/, $_);
        my $contig1 = $line[0];
        my $contig2 = $line[2];
        my $cluster1 = $contig_info->{$contig1}->{"CLUSTER"};
        my $cluster2 = $contig_info->{$contig2}->{"CLUSTER"};

        if (!exists $cluster_connection{"$cluster1-vs-$cluster2"} and
            !exists $cluster_connection{"$cluster2-vs-$cluster1"}){
            compare_cluster($cluster1, $cluster2, $cluster_species_set, $CMD, $reference_mappings, $cluster_info, $genome_db);
        }

        $cluster_connection{"$cluster1-vs-$cluster2"} = 1;
    }
    close ($CMD);
    close(EDGES);
    run_exe("cat $nucmer_dir/cmd.txt | xargs -L 1 -P $nb_process -I COMMAND sh -c \"COMMAND\" 2> $nucmer_dir/log.txt");
    if($?){
	die "Error in during mummer mapping. Please see log files $nucmer_dir/log.txt for details.\n";
    }
}


#From the intermediate MASH files, associate each cluster with a set of strains.
sub generate_clusters_species{
    my ($cluster_species_set, $species_to_ignore, $cluster_info) = @_;
    my $top_prediction = 5;
    my $nb_prediction;
    
    #Get the best species for each cluster
    
    print STDERR " *** Read mash files\n";
    my @line;my $species;
    open(FILE, "sort -k3,3 -g $mash_dir/mash_dist.dat |");
    while(<FILE>){
	@line = split(/\t/, $_);
	@tmp = split(/\//, $line[1]);
	$cluster_name = $tmp[@tmp-1];
	$cluster_name =~ s/\.fa//;
	if(! exists $cluster_species_set->{$cluster_name}){
	    $cluster_species_set->{$cluster_name} = [];#Should be an array
	    #$nb_prediction = 0;
	}
	next if(@{$cluster_species_set->{$cluster_name}} == $top_prediction);
	#print STDERR " *** $cluster_name\n";<STDIN>;
	$species = $line[0];
	push @{$cluster_species_set->{$cluster_name}},$species;
	#$nb_prediction++;
	#last if($nb_prediction == $top_prediction);
    }
    close(FILE);
}


#Get information about the contigs, their respective cluster and also informaiton
#about the cluster's coverage.
sub get_cluster_contigs_info{
    my ($contig_info, $cluster_info) = @_;
    print STDERR $clusters_file . "\n";
    open (FILE, "$clusters_file") or die;
    my @line;
    my ($contig, $cluster, $cluster_mean);
    while(<FILE>){
        chomp $_;	
        @line = split(/\t/, $_);
        $contig = $line[0];
        $cluster = $line[1];
        $cluster_mean = $line[3];
        
        #add more info to these collections if we need
        $contig_info->{$contig}->{"CLUSTER"} = $cluster;
        $contig_info->{$contig}->{"CLUSTER_MEAN"} = $cluster_mean;
        $cluster_info->{$cluster}->{"CLUSTER_MEAN"} = $cluster_mean;
	$cluster_info->{$cluster}->{"CLUSTER_SEQ"} = "";
	if (! exists $cluster_info->{$cluster}->{"CLUSTER_COUNT"}){
	    $cluster_info->{$cluster}->{"CLUSTER_COUNT"} = 0;
	}
	$cluster_info->{$cluster}->{"CLUSTER_COUNT"}++;
	    
    }
    close(FILE);
}

#Generate an intermediate file which is a preliminary graph representation
#of the clusters. The edges are not confirmed by the reference to be true yet.
sub generate_edges_between_clusters{
    my ($edges_to_cov_diff, $contig_info) = @_;

    #long read edge file
    open (FILE, "$paired_edges_file") or die;
    my @line;
    my ($contig1, $contig2);
    my ($c1_info, $c2_info, $c1, $c2, $cluster_pair_id);
    while(<FILE>){
        chomp $_;
        @line = split(/\t/, $_);
        $contig1 = $line[0];
        $contig2 = $line[2];

        #we're only interested in edges BETWEEN clusters
	$c1_info = $contig_info->{$contig1};
	$c2_info = $contig_info->{$contig2};
	#
	$c1 = $c1_info->{"CLUSTER"};
	$c2 = $c2_info->{"CLUSTER"};
	
        if ($c1 eq $c2){
            next;
        }

	#Need to change by sorting the cluster by number
	$cluster_pair_id = $c1."-vs-".$c2;
	if($c1 > $c2) {
	    $cluster_pair_id = $c2."-vs-".$c1;
	}
	
	if(! exists $edge_info->{$cluster_pair_id}){
	    $edge_info->{$cluster_pair_id} = {"SEQ", "", "CONTIG", {}};
	}

	if(! exists $edge_info->{$cluster_pair_id}->{"CONTIG"}->{$contig1}){
	    $edge_info->{$cluster_pair_id}->{"SEQ"} .= $contig_info->{$contig1}->{"SEQ"};
	    $edge_info->{$cluster_pair_id}->{"CONTIG"}->{$contig1} = 1;
	}
	if(! exists $edge_info->{$cluster_pair_id}->{"CONTIG"}->{$contig2}){
	    $edge_info->{$cluster_pair_id}->{"SEQ"} .= $contig_info->{$contig2}->{"SEQ"};
	    $edge_info->{$cluster_pair_id}->{"CONTIG"}->{$contig2} = 1;
	}
	
        $edges_to_cov_diff->{$_} = abs($c1_info->{"CLUSTER_MEAN"} - $c2_info->{"CLUSTER_MEAN"});
    }
    close (FILE);

    #short read edge file
    open (FILE, "$short_edges_file") or die;
    my ($overlap, $support);
    while(<FILE>){
        chomp $_;
        @line = split(/\t/, $_);
        $contig1 = $line[0];
        $contig2 = $line[2];
        $overlap = $line[4];
        $support = $line[6];
	#
	$c1_info = $contig_info->{$contig1};
	$c2_info = $contig_info->{$contig2};
	#
	#we're only interested in edges BETWEEN clusters
        if (! defined($c1_info) or ! defined($c2_info)){
            next;
        }
	
	$c1 = $c1_info->{"CLUSTER"};
	$c2 = $c2_info->{"CLUSTER"};

	#Not sure why we have to had the define test again
        if (! defined $c1 || ! defined $c2 || $c1 eq $c2){
            next;
        }

        if (abs($overlap) > 500 or $support < 5){
            next;
        }

	$cluster_pair_id = $c1."-vs-".$c2;
	if($c1 > $c2) {
	    $cluster_pair_id = $c2."-vs-".$c1;
	}
	
	if( ! exists $edge_info->{$cluster_pair_id}) {
	    $edge_info->{$cluster_pair_id} = {"SEQ", "", "CONTIG", {}};
	}
	if(! exists $edge_info->{$cluster_pair_id}->{"CONTIG"}->{$contig1}){
	    $edge_info->{$cluster_pair_id}->{"SEQ"} .= $contig_info->{$contig1}->{"SEQ"};
	    $edge_info->{$cluster_pair_id}->{"CONTIG"}->{$contig1} = 1;
	}
	if(! exists $edge_info->{$cluster_pair_id}->{"CONTIG"}->{$contig2}){
	    $edge_info->{$cluster_pair_id}->{"SEQ"} .= $contig_info->{$contig2}->{"SEQ"};
	    $edge_info->{$cluster_pair_id}->{"CONTIG"}->{$contig2} = 1;
	}
	
	
        $edges_to_cov_diff->{$_} = abs($c1_info->{"CLUSTER_MEAN"} - $c2_info->{"CLUSTER_MEAN"});
    }
    close (FILE);

    
    open (FILE, ">", "$edges_between_clusters_file") or die;
    #Sort collection by coverage difference between contigs. The order matters when clustering since this is a greedy
    #algorithm.
    foreach my $edge (sort {$edges_to_cov_diff->{$a} <=> $edges_to_cov_diff->{$b}} keys %{$edges_to_cov_diff}){
        print FILE $edge . "\n";
        #print STDERR $edges_to_cov_diff -> {$edge} . "\n";
    }
    
    close (FILE);
}

#Run mash on all of the clusters to get their top species.
sub run_mash_on_clusters{
    my ($cluster_info, $contig_info, $mash_ref) = @_;
    my $mash_cmd_file = "$mash_dir/cmd.txt";
    run_exe("rm -r $mash_dir");
    run_exe("mkdir -p $mash_dir");
    my $out;

    my $inter_fa_dir = "$mash_dir/INTER_FILES/";
    run_exe("mkdir $inter_fa_dir");
    my $partial_sketch_dir = "$mash_dir/PARTIAL_SKETCH/";
    run_exe("mkdir $partial_sketch_dir");
    
    
    #extract the sequence of each cluster
    open(FILE, $contig_seq_file) or die "Could not open contigs file $contig_seq_file $!.\n";
    my ($contig_name, $contig_seq);
    $contig_size = 0;
    
    #Generation of the cluster fasta file
    my %cluster_length = ();
    my @tmp;my $cluster_name;my $contig_cluster;
    my $count = 0;
    my $partial_count = 0;
    my $cmd = "";
    my $max_genome_for_sketch = 500;
    while(<FILE>){
	#last if($comp == $max_contig);
	
        chomp $_;
        if($_ ne ""){
            if($_ =~ /^>(.+)/){
		@tmp = split(/\s+/, $1);

		if(defined($contig_name) && defined($contig_cluster)){
		    #print STDERR " Write " . $cluster_info->{$contig_info->{$contig_name}->{"CLUSTER"}}->{"OUT_FILE"} . "\n";<STDIN>;
		    #create multiple folders, each folder contain 50 .fa files
		    $contig_info->{$contig_name}->{"SEQ"} .= $contig_seq;

		    # concat all contig_seq belong to the cluster and reduce the contig count
		    $cluster_info->{$contig_cluster}->{"CLUSTER_SEQ"} .= $contig_seq ;
		    $cluster_info->{$contig_cluster}->{"CLUSTER_COUNT"}-- ;

		    # if all contig in the cluster was added, write to file
		    if ($cluster_info->{$contig_cluster}->{"CLUSTER_COUNT"} == 0){
			open($out,'>', "$inter_fa_dir/$contig_cluster.fa");
			print $out $cluster_info->{$contig_cluster}->{"CLUSTER_SEQ"};
			close($out);
			$count++;
			$cluster_info->{$contig_cluster}->{"CLUSTER_SEQ"} = "";
			# only generate partial sketch if there are 50 fa files
			if ($count == $max_genome_for_sketch){
			    $count = 0;
			    # generate partial sketeches
			    run_exe("$mash_exe_dir/mash sketch -p $nb_process -k 21 -s 1000 -o $partial_sketch_dir/partial_Sketch$partial_count $inter_fa_dir/* > $mash_dir/mash_sketch.out 2> $mash_dir/mash_sketch.err");
			    if($?){
				die "Error in during bin/mash sketch. Please see $mash_dir/mash_sketch.out $mash_dir/mash_sketch.err for details.\n";
			    }
			    
			    run_exe("$mash_exe_dir/mash dist -p $nb_process -d 0.90  $mash_ref $partial_sketch_dir/partial_Sketch$partial_count.msh  > $mash_dir/mash_dist_$partial_count.dat 2> $mash_dir/mash_dist.err");	    
			    if($?){
				die "Error in during bin/mash dist. Please see $mash_dir/mash_dist.err for details.\n";				
			    }

			    
			    # remove intermediate .fa files
			    run_exe("rm $inter_fa_dir/*");
			    $partial_count++;
			}
		    }
		}
		
 		$contig_seq = $_ . "\n";
		$contig_name = $tmp[0];
		$contig_cluster = $contig_info->{$contig_name}->{"CLUSTER"};
		$contig_size = 0;
		#<STDIN>;
		#print STDERR "\n".$contig_seq;
            }
            else{
		$contig_seq .= $_ . "\n";
	    }
        }
    }
    close(FILE);
    #Write the last contig in its cluster
    $contig_info->{$contig_name}->{"SEQ"} .= $contig_seq;
    $cluster_info->{$contig_cluster}->{"CLUSTER_SEQ"} .= $contig_seq ;
    if(defined($contig_name) && defined($contig_cluster)){
	open($out,'>',"$inter_fa_dir/$contig_cluster.fa");
	print $out $cluster_info->{$contig_cluster}->{"CLUSTER_SEQ"};
	close($out);
	$count++;
    }
    if($count != 0){
	run_exe("$mash_exe_dir/mash sketch -p $nb_process -k 21 -s 1000 -o $partial_sketch_dir/partial_Sketch$partial_count $inter_fa_dir/* > $mash_dir/mash_sketch.out 2> $mash_dir/mash_sketch.err");
	if($?){
	    die "Error in during bin/mash sketch. Please see $mash_dir/mash_sketch.out $mash_dir/mash_sketch.err for details.\n";
	}
	run_exe("rm -r $inter_fa_dir");
    }
    
    run_exe("$mash_exe_dir/mash dist -p $nb_process -d 0.90  $mash_ref $partial_sketch_dir/partial_Sketch$partial_count.msh  > $mash_dir/mash_dist_$partial_count.dat 2> $mash_dir/mash_dist.err");
    #run_exe("$mash_exe_dir/mash paste $mash_dir/cluster.msh $partial_sketch_dir/partial_Sketch*");
    if($?){
	die "Error in during bin/mash dist. Please see $mash_dir/mash_dist.err for details.\n";
    }
    run_exe("cat $mash_dir/mash_dist_*.dat > $mash_dir/mash_dist.dat");
    #<STDIN>;
}

sub compare_cluster{
    my ($cluster_1, $cluster_2, $cluster_species_set, $CMD, $reference_mappings, $cluster_info, $genome_db) = @_;
    
    #print STDERR "%$cluster_species_set\n";

    $genome_list_1 = $cluster_species_set->{$cluster_1};
    $genome_list_2 = $cluster_species_set->{$cluster_2};
    if (!defined $genome_list_1 or !defined $genome_list_2){
        return;
    }

    #Creat the edge file
    my $cluster_pair_name = "$cluster_1-vs-$cluster_2";
    open(OUT, ">$nucmer_dir/$cluster_pair_name.fa");
    if( exists $edge_info->{$cluster_pair_name}){
	print OUT $edge_info->{$cluster_pair_name}->{"SEQ"};
    }
    else{
	print OUT $edge_info->{$cluster_2."-vs-".$cluster_1}->{"SEQ"};
    }
    close(OUT);
    
    print STDERR " *** compare_cluster $cluster_1 $cluster_2\n";
    
    #compute the genome in common
    my @common_genome = ();
    foreach my $g1 (@{$genome_list_1}){ 
        foreach my $g2 (@{$genome_list_2}){
            if($g2 eq $g1){
		push @common_genome, $g1;
                #print STDERR "$s2, $s\n";
            }
        }
    }
     
    my $count = 0;
    my $ref_genome_file = "";
    foreach my $genome (@common_genome){
		
	$ref_genome_file = "$nucmer_dir/temp_genome/$cluster_pair_name\_$count.fa";
	$command = "zcat $genome_db/$genome > $ref_genome_file;${mummer_dir}nucmer --maxmatch -c 400 --banded $nucmer_dir/$cluster_pair_name.fa $opera_ms_dir/$ref_genome_file -p $nucmer_dir/$cluster_pair_name\_$count >> $nucmer_dir/LOG.txt;${mummer_dir}show-coords -lrcT $nucmer_dir/$cluster_pair_name\_$count.delta > $nucmer_dir/$cluster_pair_name\_$count.txt;rm $ref_genome_file";
	
        $count++;
        print $CMD $command . "\n";

        if(! exists $reference_mappings->{$cluster_pair_name}->{"NUMB_SPECIES"}){
            $reference_mappings->{$cluster_pair_name}->{"NUMB_SPECIES"} = 0;
        }
	$reference_mappings->{$cluster_pair_name}->{"NUMB_SPECIES"}++;
    }
}


sub run_exe{
    my ($exe) = @_;
    my $run = 1;
    print STDERR "\n".$exe."\n";#<STDIN>;
    print STDERR `$exe` if($run);
}

sub check_reference_position{
    my ($contig1, $contig2, $length, $stdev, $reference_mappings, $contig_info) = @_; 
    
    my $clust1 = $contig_info->{$contig1}->{"CLUSTER"};
    my $clust2 = $contig_info->{$contig2}->{"CLUSTER"};
   my $mapfile_name = "";
   my $result = 0;

   #Only proceed if two clusters have common species containing the two contigs.
   if (exists $reference_mappings->{"$clust1-vs-$clust2"}){
        $mapfile_name = "$clust1-vs-$clust2";
   }

   elsif (exists $reference_mappings->{"$clust2-vs-$clust1"}){
        $mapfile_name = "$clust2-vs-$clust1";
   }

   else {return 0; }

   my $number_of_files = $reference_mappings->{$mapfile_name}->{"NUMB_SPECIES"};
   my @contig_mapping1 = ();
   my @contig_mapping2 = ();
   my @contig_percent_mapped1 = ();
   my @contig_percent_mapped2 = ();
   my $ref_length;

   #print STDERR "CHECKING REF POSITION : $contig1($clust1), $contig2($clust2), $length, $number_of_files\n";
   my @line;
   my @positions_reference;
   my ($contig, $percent_mapped, $index);
   for (my $i=0; $i < $number_of_files; $i++){
       open(NUC_MAPPING, "$nucmer_dir/$mapfile_name\_$i.txt")
	   or die "Parsing problem during read rescue using NUCMER mapping.\n";

        #skip the first four lines of the cluster-vs-cluster.txt file.
        <NUC_MAPPING>;
        <NUC_MAPPING>;
        <NUC_MAPPING>;
        <NUC_MAPPING>;
        while(<NUC_MAPPING>){
           if ($_ eq ""){
                next;
           }
           #print STDERR $_;
           @line = split(/\t/, $_);
           $contig = $line[11];
           $ref_length = $line[8];
           $percent_mapped = $line[9];
	   
           #Scan each line of the cluster-vs-cluster file for the contigs
           #we want and find their positions on the reference.
           if($contig eq $contig1){
                $index = 0;
           }
	   elsif($contig eq $contig2){
	       $index = 1;
           }
	   else {next;}
	   
           @positions_reference = ($line[2],$line[3]);
           my @sorted_positions = sort{ $a <=> $b }(@positions_reference); 
           #print STDERR @sorted_positions . "\n";
            

           #Each contig_mapping file is an array of contigs whose entries are an array
           #of two numbers : their starting and end positions.
           if($index){
               push @contig_mapping2 ,[@sorted_positions]; #Is this a copy or a point !!!!!
               push @contig_percent_mapped2, $percent_mapped;
               #my $toprint = $contig_mapping2[@contig_mapping2-1];
               #print STDERR "CONTIG MAPPING 2 @$toprint\n";
               #print Dumper(@contig_mapping2) . "\n";
           }
           
           else{
               push @contig_mapping1,[ @sorted_positions];
               push @contig_percent_mapped1, $percent_mapped;
           }
        }
        close(NUC_MAPPING);
        if (!@contig_mapping1 or !@contig_mapping2){ next; }

	#Repeat detection based on fraction of the contig mapped
        if (sum(\@contig_percent_mapped1) > 150){ 
            #print STDERR "REPEAT DETECTED $contig1 \n";
            next;
        }

        if(sum(\@contig_percent_mapped2) > 150){
            #print STDERR "REPEAT DETECTED $contig2\n";
            next;
        }


        foreach my $contig_ref_positions1(@contig_mapping1){
            if ($result){
                last;
            }

            foreach my $contig_ref_positions2(@contig_mapping2){
                my $mindist;
                my $maxdist;
                my $maxdist_circ;
                my $mindist_circ;
                if (@$contig_ref_positions2[1] > @$contig_ref_positions1[1]){
                    $mindist = @$contig_ref_positions2[0] - @$contig_ref_positions1[1] - 3*abs($stdev);
                    $maxdist = @$contig_ref_positions2[1] - @$contig_ref_positions1[0] + 3*abs($stdev);

                    $mindist_circ = @$contig_ref_positions1[0] - @$contig_ref_positions2[1] - 3*abs($stdev) + $ref_length;
                    $maxdist_circ = @$contig_ref_positions1[1] - @$contig_ref_positions2[0] + 3*abs($stdev) + $ref_length;
                }
                else{
                    $mindist = @$contig_ref_positions1[0] - @$contig_ref_positions2[1] - 3*abs($stdev);
                    $maxdist = @$contig_ref_positions1[1] - @$contig_ref_positions2[0] + 3*abs($stdev);

                    $mindist_circ = @$contig_ref_positions2[0] - @$contig_ref_positions1[1] - 3*abs($stdev) + $ref_length;
                    $maxdist_circ = @$contig_ref_positions2[1] - @$contig_ref_positions1[0] + 3*abs($stdev) + $ref_length;

                }
                
                #print STDERR "QUERY EDGE : $length = length, $mindist = mindist, $maxdist = maxdist, $contig1($clust1), $contig2($clust2), $mapfile_name\n";
                if (($maxdist > $length && $length > $mindist) or 
                    ($maxdist_circ > $length && $mindist_circ < $length)){
                    print STDERR @$contig_ref_positions2[1] . " " . @$contig_ref_positions2[0]
                    ." " .  @$contig_ref_positions1[1] . " " . @$contig_ref_positions1[0] . " \n";

                    print STDERR @$contig_ref_positions1[1] . " " . @$contig_ref_positions1[0]
                    ." " .  @$contig_ref_positions2[1] . " " . @$contig_ref_positions2[0] . " \n";

                    print STDERR "RESCUED EDGE : $length = length, $mindist = mindist, $maxdist = maxdist, $contig1($clust1), $contig2($clust2), $mapfile_name\n";
                    $result = 1;
                    last;
                }
            }
        }

        if ($result){
            last;
        }
   }

   return $result;
}

sub get_species_name{
    my ($strain) = @_;
    my @full_path = split(/\//, $strain);
    my $strain_name = $full_path[@full_path-1];
    my @strain_delim = split(/__/,$strain_name);
    return $strain_delim[0];
}


sub sum{
    my ($array) = @_;
    my $res = 0;
    foreach $s (@{$array}){
	$res += $s;
    }
    return $res;
}

