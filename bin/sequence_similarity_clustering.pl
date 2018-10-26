#!/usr/bin/perl
use Data::Dumper;
use warnings;
use Statistics::Basic qw(:all);
use strict;
use List::Util qw(sum first);

#############################################################
#Script for merging sigma clusters based on reference genomes.
#1. Gather all clusters with connections between them
#2. For all clusters with edges between them, check if they have similar reference genomes.
#3. Recover edges between clusters in within the super clusters.
#############################################################
my ($inter_dir, $contig_seq_file, $nb_process, $opera_ms_dir, $mash_exe_dir, $mummer_dir) = @ARGV;

my $PARSE_FILES_MASH = 1;
my $PARSE_FILES_NUC = 1;
my $MERGE_CLUSTERS = 1;

my %reference_mappings_ = ();
my %contig_info_ = ();
my %cluster_info_ = ();
my %cluster_species_set_ = ();
my %super_cluster_species_set_ = ();
my %species_to_ignore_ = ();
my @super_cluster = ();
my %cluster_to_super_cluster_ID;
my %super_cluster_strain_set_ = ();

my $ref_map_folder = "$inter_dir/reference_mapping";
my $cluster_results_file = "$inter_dir/reference_mapping/clusters_seq_similarity";
my $cluster_species_file = "$inter_dir/reference_mapping/cluster_species.dat";
my $super_cluster_file = "$inter_dir/reference_mapping/super_cluster.dat";
my $edges_between_clusters_file = "$inter_dir/reference_mapping/edges_between_clusters";
my $good_edges_file = "$inter_dir/reference_mapping/edges_between_clusters_good";
my $clusters_file = "$inter_dir/sigma/clusters";
my $paired_edges_file = "$inter_dir/long-read-mapping/pairedEdges";
my $paired_edges_file_above_thresh = "$inter_dir/long-read-mapping/pairedEdges_above_thresh";
my $short_edges_file = "$inter_dir/lib_1_bundles/clusters_contigs";
my $mash_dir = "$inter_dir/reference_mapping/MASH/";
my $mash_ref = "$opera_ms_dir/genomeDB_Sketch.msh";
my $nucmer_dir = "$inter_dir/reference_mapping/NUCMER_OUT/"; 

######################################START########################3

my $command = "mkdir -p $ref_map_folder"; 
run_exe($command) if (!(-d $ref_map_folder));

get_cluster_contigs_info(\%contig_info_, \%cluster_info_);
if ($PARSE_FILES_MASH){
    my %edges_to_cov_diff = ();
    generate_edges_between_clusters(\%edges_to_cov_diff, \%contig_info_);
    run_mash_on_clusters(\%cluster_info_, \%contig_info_, $contig_seq_file);
}

generate_clusters_species(\%cluster_species_set_, \%species_to_ignore_);

if ($PARSE_FILES_NUC){
    run_nucmer(\%cluster_species_set_, \%contig_info_, \%reference_mappings_);
    get_correct_edges(\%reference_mappings_, \%contig_info_);
}

if ($MERGE_CLUSTERS){
    merge_clusters(\%contig_info_, \%cluster_species_set_, \%super_cluster_species_set_, \%super_cluster_strain_set_);
    recover_edges(\%contig_info_);
    export_info(\%contig_info_, \%cluster_species_set_, \%super_cluster_strain_set_);
}
###################################END################################

#Export some info about the clusters. The clusters file, the species file, and a super cluster file.
sub export_info{
    my ($contig_info,$cluster_species_set, $super_cluster_strain_set) = @_;
    open (OUT, ">", $cluster_results_file) or die;
    #print STDOUT "EXPORTING INFO\n";
    foreach my $contig (keys %{$contig_info}){
        my $cluster = $contig_info->{$contig}->{"CLUSTER"};
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
        foreach my $strain (sort {$super_cluster_strain_set->{$id}->{$b} <=> $super_cluster_strain_set->{$id}->{$a}} 
        keys %{$super_cluster_strain_set->{$id}}){
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
        open ($OUT, ">","$inter_dir/reference_mapping/filtered_pairedEdges_i$i") or die;
        $outfiles[$i] = $OUT;
    }

    #Add edges between cluster. Those edge can be below the support threshold
    while (<GOOD>){
        chomp $_;
        my @line = split(/\t/, $_);
        my $contig1 = $line[0];
        my $contig2 = $line[2];
        my $dist = $line[4];
        my $cluster1 = $contig_info->{$contig1}->{"CLUSTER"};
        my $cluster2 = $contig_info->{$contig2}->{"CLUSTER"};

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
    open (FILE, "$paired_edges_file_above_thresh") or die;
    while (<FILE>){
        chomp $_;
        my @line = split(/\t/, $_);
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
    while(<FILE>){
        my @line = split(/\t/, $_);
        my $cluster_1 = $contig_info->{$line[0]}->{"CLUSTER"};
        my $cluster_2 = $contig_info->{$line[2]}->{"CLUSTER"};
        my $super1;
        my $super2;

        if (exists $cluster_to_super_cluster_ID{$cluster_1}){
            $super1 = $cluster_to_super_cluster_ID{$cluster_1};
        }

        else{
            $super1 = "";
        }

        if (exists $cluster_to_super_cluster_ID{$cluster_2}){
            $super2 = $cluster_to_super_cluster_ID{$cluster_2};
        }

        else{
            $super2 = "";
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
                
                my @species_intersection = ();
                undef(@species_intersection); 
                my %strain_set = ();
                undef(%strain_set);
                
                foreach my $s ( @{$cluster_species_set->{$cluster_1}}){

                    my @full_path = split(/\//,$s);
                    my $strain = $full_path[1];
                    my @strain_delim = split(/_/,$strain);
                    my $species = $strain_delim[0] . "_" . $strain_delim[1];

                    foreach my $s2 (@{$cluster_species_set->{$cluster_2}}){

                        my @full_path2 = split(/\//,$s2);
                        my $strain2 = $full_path[1];
                        #print STDERR "STRAIN " . $strain2 . "\n";
                        my @strain_delim2 = split(/_/,$strain);
                        my $species2 = $strain_delim[0] . "_" . $strain_delim[1];

                        if($species2 eq $species){
                            if (!(first {$_ eq $species} @species_intersection)){
                                push @species_intersection, $species;
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
                $super_cluster_species_set->{$super_cluster_ID} = [@species_intersection]; 
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
                    my @species_intersection = ();
                    undef(@species_intersection); 

                    foreach my $s (@{$super_cluster_species_set->{$super_cluster_ID}}){
                        foreach my $s2 (@{$cluster_species_set->{$cluster_2}}){

                        my @full_path = split(/\//,$s2);
                        my $strain = $full_path[1];
                        my @strain_delim = split(/_/,$strain);
                        my $species = $strain_delim[0] . "_" . $strain_delim[1];

                            if($species eq $s){
                                if (!(first {$_ eq $species} @species_intersection)){
                                    push @species_intersection, $species;
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
                    
                    $super_cluster_species_set->{$super_cluster_ID} = [@species_intersection]; 
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

                        my @species_intersection = ();
                        #Needed to clear out the array... some scope/memory issue
                        undef(@species_intersection); 

                        foreach my $s (@{$super_cluster_species_set->{$super_cluster_ID_1}}){
                            foreach my $s2 (@{$super_cluster_species_set->{$super_cluster_ID_2}}){
                                if ($s2 eq $s){
                                    push @species_intersection, $s;
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

                        $super_cluster_species_set->{$super_cluster_ID_1} = [@species_intersection]; 
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
    my @common_species = ();
    #print STDERR "%$cluster_species_set\n";

    if (!defined $cluster_species_set->{$cluster_1} or !defined $cluster_species_set->{$cluster_2}){
        return;
    }

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
            my @full_path = split(/\//,$s);
            my $strain = $full_path[1];
            my @strain_delim = split(/_/,$strain);
            $species = $strain_delim[0] . "_" . $strain_delim[1];
        }

        else{
            $species = $s;
        }

        foreach my $s2 (@species2){
            my $species2;

            if(!$species2_super){
                my @full_path = split(/\//,$s2);
                my $strain = $full_path[1];
                my @strain_delim = split(/_/,$strain);
                $species2 = $strain_delim[0] . "_" . $strain_delim[1];
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
        my @line = split(/\t/, $_);
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
    my %cluster_connection = ();
    run_exe("rm -r $nucmer_dir; mkdir -p $nucmer_dir"); 
    open (my $CMD, ">", "$nucmer_dir/cmd.txt") or die;
    my ($cluster_species_set, $contig_info, $reference_mappings) = @_;
    foreach my $cluster (keys %{$cluster_species_set}){
        if (@{$cluster_species_set->{$cluster}} != 0){
            my $top_species = @{$cluster_species_set->{$cluster}}[0];
            my $command = "${mummer_dir}nucmer --maxmatch -c 400 --banded $mash_dir/$cluster.fa $opera_ms_dir/$top_species -p $nucmer_dir/$cluster-repeat_detection"; 
            print $CMD $command . "\n";
        }
    }

    open (EDGES, $edges_between_clusters_file) or die;
    while(<EDGES>){
        chomp $_;
        my @line = split(/\t/, $_);
        my $contig1 = $line[0];
        my $contig2 = $line[2];
        my $cluster1 = $contig_info->{$contig1}->{"CLUSTER"};
        my $cluster2 = $contig_info->{$contig2}->{"CLUSTER"};

        if (!exists $cluster_connection{"$cluster1-vs-$cluster2"} and
            !exists $cluster_connection{"$cluster2-vs-$cluster1"}){
            compare_cluster($cluster1, $cluster2, $cluster_species_set, $CMD, $reference_mappings);
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
    my ($cluster_species_set, $species_to_ignore) = @_;
    my $top_prediction = 5;
    my $nb_prediction;
    opendir(DIR, $mash_dir);
    my @cluster_mash = readdir(DIR);

    #Get the best species for each cluster
    print STDERR " *** Read mash files\n";
    foreach my $cluster (@cluster_mash){
        if($cluster =~ /(.+)\.dat/){
            my $cluster_name = $1;
            #print STDERR " *** ".$cluster. "\t" . $cluster_name. "\n";
            open(FILE, "sort -k3,3 -g $mash_dir/$cluster |");
            $cluster_species_set->{$cluster_name} = [];
            $nb_prediction = 0;
            while(<FILE>){
                my @line = split(/\t/, $_);
                #$line[0] looks is the species ID, it looks like
                #/mnt/path/to/#####.fna
                my $species = $line[0];
                my @delim_species = split(/\//, $species);
                #print STDOUT $delim_species[7] . "\n";
                #if (exists $species_to_ignore->{$delim_species[7]}){
                #    next;
                #}

                push @{$cluster_species_set->{$cluster_name}},$species;
                $nb_prediction++;
                last if($nb_prediction == $top_prediction);
            }
            close(FILE);
        }
    }
}

#Get information about the contigs, their respective cluster and also informaiton
#about the cluster's coverage.
sub get_cluster_contigs_info{
    my ($contig_info, $cluster_info) = @_;
    open (FILE, "$clusters_file") or die;
    while(<FILE>){
        chomp $_;	
        my @line = split(/\t/, $_);
        my $contig = $line[0];
        my $cluster = $line[1];
        my $cluster_mean = $line[3];
        
        #add more info to these collections if we need
        $contig_info->{$contig}->{"CLUSTER"} = $cluster;
        $contig_info->{$contig}->{"CLUSTER_MEAN"} = $cluster_mean;
        $cluster_info->{$cluster}->{"CLUSTER_MEAN"} = $cluster_mean;
    }
    close(FILE);
}

#Generate an intermediate file which is a prelimniary graph representation
#of the clusters. The edges are not confirmed by the reference to be true yet.
sub generate_edges_between_clusters{
    my ($edges_to_cov_diff, $contig_info) = @_; 
    open (FILE, "$paired_edges_file") or die;
    while(<FILE>){
        chomp $_;
        my @line = split(/\t/, $_);
        my $contig1 = $line[0];
        my $contig2 = $line[2];

        #we're only interested in edges BETWEEN clusters
        if ($contig_info->{$contig1}->{"CLUSTER"} eq $contig_info->{$contig2}->{"CLUSTER"}){
            next;
        }

        $edges_to_cov_diff->{$_} = abs($contig_info->{$contig1}->{"CLUSTER_MEAN"} - $contig_info->{$contig2}->{"CLUSTER_MEAN"});
    }
    close (FILE);

    open (FILE, "$short_edges_file") or die;
    while(<FILE>){
        chomp $_;
        my @line = split(/\t/, $_);
        my $contig1 = $line[0];
        my $contig2 = $line[2];
        my $overlap = $line[4];
        my $support = $line[6];

        #we're only interested in edges BETWEEN clusters
        if (!exists $contig_info->{$contig1} or !exists $contig_info->{$contig2}){
            next;
        }

        if ($contig_info->{$contig1}->{"CLUSTER"} eq $contig_info->{$contig2}->{"CLUSTER"}){
            next;
        }

        if (abs($overlap) > 500 or $support < 5){
            next;
        }

        $edges_to_cov_diff->{$_} = abs($contig_info->{$contig1}->{"CLUSTER_MEAN"} - $contig_info->{$contig2}->{"CLUSTER_MEAN"});
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
    my ($cluster_info, $contig_info) = @_;
    my $mash_cmd_file = "$mash_dir/cmd.txt";
    run_exe("rm -r $mash_dir");
    run_exe("mkdir -p $mash_dir");
    my $out;
    
    foreach my $clust (keys %{$cluster_info}){
        my $file_writer;

        #May have to had a filtering based on cluster length
        open($file_writer, ">$mash_dir/$clust.fa");
        $cluster_info->{$clust}->{"OUT_FILE"} = "$mash_dir/$clust.fa";#$file_writer;
        close($file_writer)
    }

    #extract the sequence of each cluster
    open(FILE, $contig_seq_file) or die "Could not open contigs file $contig_seq_file $!.\n";
    my ($contig_name, $contig_seq, $contig_size);
    $contig_size = 0;
    
    my %cluster_length = ();
    while(<FILE>){
	#last if($comp == $max_contig);
        chomp $_;
        if($_ ne ""){
            if($_ =~ /^>(.+)/){
            my @tmp = split(/\s+/, $1);

            if(defined($contig_name)){
                my $cluster_name = $contig_info->{$contig_name}->{"CLUSTER"};
                if(defined $cluster_name &&
                   exists $cluster_info->{$cluster_name}->{"OUT_FILE"}){
		    #print STDERR " Write " . $cluster_info->{$contig_info->{$contig_name}->{"CLUSTER"}}->{"OUT_FILE"} . "\n";<STDIN>;
                    #$out = $cluster_info->{$contig_info->{$contig_name}->{"CLUSTER"}}->{"OUT_FILE"};
                    open($out,'>>',$cluster_info->{$contig_info->{$contig_name}->{"CLUSTER"}}->{"OUT_FILE"});
                    if(exists ($cluster_length{$contig_info->{$contig_name}->{"CLUSTER"}})){
                        $cluster_length{$contig_info->{$contig_name}->{"CLUSTER"}}->{"SIZE"} += $contig_size;
                    }

                    else{
                        $cluster_length{$contig_info->{$contig_name}->{"CLUSTER"}} = {"SIZE", $contig_size};
                    }
                    print $out $contig_seq;
                    close($out);
                }
            }

            $contig_seq = $_."\n";
            $contig_name = $tmp[0];
            $contig_size = 0;
            #<STDIN>;
            #print STDERR "\n".$contig_seq;
            }
            else{
            $contig_seq = $contig_seq.$_."\n";
            $contig_size += length($_); 
            }
        }
    }
    close(FILE);
    
    #foreach my $clust (keys %{$cluster_info}){
    #close($cluster_info->{$clust}->{"OUT_FILE"}); 
    #}

    #Run mash
    open (my $CMD, ">", $mash_cmd_file) or die "Could not write to mash command file at $mash_cmd_file\n";
    
    foreach my $clust (keys %{$cluster_info}){
	if(exists $cluster_info->{$clust}->{"OUT_FILE"}){
	    if ($cluster_length{$clust}->{"SIZE"} > 0){
	        # run_exe("$mash_exe_dir/mash dist -d 0.90 $mash_ref $mash_dir/$clust.fa > $mash_dir/$clust.dat");
		#`$qsub -l mem_free=2G,h_rt=00:30:0 -pe OpenMP 1 -N mashing -e $out_dir/LOG/$clust.err -o $mash_dir/$clust.dat -b y $mash_exe_dir/mash dist -d 0.90 $mash_ref $mash_dir/$clust.fa`;
		print $CMD "${mash_exe_dir}mash dist -d 0.90 $mash_ref $mash_dir/$clust.fa > $mash_dir/$clust.dat\n";
	    }
	    
	    #close($cluster_info->{$clust}->{"OUT_FILE"});
	}
    }

    close ($CMD);
    run_exe("cat $mash_cmd_file | xargs -L 1 -P $nb_process -I COMMAND sh -c \"COMMAND\" 2> $mash_dir/log.txt");
    if($?){
	die "Error in during mash distance computation. Please see log files for details.\n";
    }
}

sub compare_cluster{
    my ($cluster_1, $cluster_2, $cluster_species_set, $CMD, $reference_mappings) = @_;
    my @common_species = ();
    #print STDERR "%$cluster_species_set\n";

    if (!defined $cluster_species_set->{$cluster_1} or !defined $cluster_species_set->{$cluster_2}){
        return;
    }

    my @species1 = @{$cluster_species_set->{$cluster_1}};
    my @species2 = @{$cluster_species_set->{$cluster_2}};

    print STDERR " *** compare_cluster $cluster_1 $cluster_2\n";
    my $res = 0;
    #Species in common
    foreach my $s (@species1){ 
        foreach my $s2 (@species2){
            if($s2 eq $s){
                $res = 1;
                push @common_species, $s;
                #print STDERR "$s2, $s\n";
            }
        }
    }

    my $count = 0;
    foreach my $species (@common_species){
        my $command = "cat $mash_dir/$cluster_1.fa $mash_dir/$cluster_2.fa > $nucmer_dir/$cluster_1-vs-$cluster_2.fa";
        run_exe($command);

        $command = "${mummer_dir}nucmer --maxmatch -c 400 --banded $nucmer_dir/$cluster_1-vs-$cluster_2.fa $opera_ms_dir/$species -p $nucmer_dir/$cluster_1-vs-$cluster_2\_$count >> $nucmer_dir/LOG.txt;${mummer_dir}show-coords -lrcT $nucmer_dir/$cluster_1-vs-$cluster_2\_$count.delta > $nucmer_dir/$cluster_1-vs-$cluster_2\_$count.txt"; 
        $count++;
        print $CMD $command . "\n";

        if(! exists $reference_mappings->{"$cluster_1-vs-$cluster_2"}->{"NUMB_SPECIES"}){
            $reference_mappings->{"$cluster_1-vs-$cluster_2"}->{"NUMB_SPECIES"} = 0;
        }
	$reference_mappings->{"$cluster_1-vs-$cluster_2"}->{"NUMB_SPECIES"}++;
        #else{
        #    $reference_mappings->{"$cluster_1-vs-$cluster_2"}->{"NUMB_SPECIES"} = 1;
        #}
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
           my @line = split(/\t/, $_);
           my $contig = $line[11];
           $ref_length = $line[8];
           my $percent_mapped = $line[9];
           
           my $index;
           
           #Scan each line of the cluster-vs-cluster file for the contigs
           #we want and find their positions on the reference.
           if($contig eq $contig1){
                $index = 0;
           }

           elsif($contig eq $contig2){
                $index = 1;
           }

           else {next;}

           my @positions_reference = ($line[2],$line[3]);
           my @sorted_positions = sort{ $a <=> $b }(@positions_reference); 
           #print STDERR @sorted_positions . "\n";
            

           #Each contig_mapping file is an array of contigs whose entries are an array
           #of two numbers : their starting and end positions.
           if($index){
               push @contig_mapping2 ,[@sorted_positions];
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
        if (sum(@contig_percent_mapped1) > 150){ 
            #print STDERR "REPEAT DETECTED $contig1 \n";
            next;
        }

        if(sum(@contig_percent_mapped2) > 150){
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
