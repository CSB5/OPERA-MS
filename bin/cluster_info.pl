#!/usr/bin/perl
use warnings;
use Statistics::Basic qw(:all);
use Switch;


my ($assembly_dir, $bin_dir) = @ARGV;

my %contig_to_bin;
my %contig_to_illumina_assembly;#for each opera-ms contigs identify the contig and their sizes that were used => allows to compute the illumina N50 and assembly size => allows to indentify the size of gaps as well
my %bin_info = ();

analyze_scaffold("$assembly_dir/intermediate_files/opera_long_read/scaffolds.scaf", \%contig_to_illumina_assembly);

#Get the contigs in bin
analyze_bin("$bin_dir", \%contig_to_bin, \%bin_info);

#Get contig info
analyze_contig_info("$assembly_dir/contig_info.txt", \%contig_to_bin, \%bin_info, \%contig_to_illumina_assembly);

#Get checkm info
#analyze_checkm("$bin_dir/eval.dat", \%bin_info);

#Get mash information
analyze_mash("$bin_dir/mash.dist", \%bin_info);

#Get the kraken info
#analyze_kraken("$assembly_dir/kraken/", "species", \%bin_info);
#analyze_kraken("$assembly_dir/kraken/", "genus", \%bin_info);

print_info(\%bin_info);


sub print_info{
    my ($bin_info) = @_;

    my @col_order = (
	"CLUSTER_ID",
	#
	"SHORT_READ_COV",
	"LONG_READ_COV",
	#
	#"COMLETENESS",
	#"CONTAMINATION",
	#
	"SPECIES", 
	"SIMILARITY", 
	"CLUSTER_SIZE",
	"NB_CONTIG",
	"LONGEST_CONTIG", 
	#
	#"ILLUMINA_N90",
	#"CONTIG_SIZE", 
	"N50", 
	"L50", 
	#"NG50",
	"N90",
	"L90",
	#"NG90",
	#
	"SHORT_READ_NB_CONTIG",
	"SHORT_READ_LONGEST_CONTIG", 
	"SHORT_READ_N50",
	"SHORT_READ_L50",
	#"FRAQ_READ",
	#"ILLU_S_ABUND",
	#"ILLU_G_ABUND",
	#"LONG_S_ABUND",
	#"LONG_G_ABUND"
	);

    
    print "" . join("\t", @col_order) . "\n";

    $bin_info->{"unclustered"}->{"SPECIES"} = "NA";
    $bin_info->{"unclustered"}->{"SIMILARITY"} = "NA";
    
    foreach $bin (sort {$bin_info->{$b}->{"N50"} <=> $bin_info->{$a}->{"N50"}} keys %{$bin_info}){
	print $bin;
	foreach $info (@col_order){
	    next if($info eq "CLUSTER_ID");
	    print "\t" .  $bin_info->{$bin}->{$info};
	    #print STDERR $bin . "\t" .$info . "\n";
	}
	print "\n";
    }
    
}

sub analyze_scaffold{
    my ($scaffold_file, $contig_to_illumina_assembly) = @_;

    my %contig_name = ();
    open(FILE, "$scaffold_file.cname");
    while(<FILE>){
	chop $_;
	($scaff_name, $c_name) = split(/\t/, $_);
	$contig_name{$scaff_name} = $c_name;
    }
    close(FILE);
    
    #
    #open(FILE, "sed 's/scaffold/contig/' $scaffold_file |");
    open(FILE, $scaffold_file);
    my $opera_contig = "";
    while(<FILE>){
	chop $_;
	@line = split(/\t/, $_);
	if ($line[0] =~ />(.*)/){
	    $opera_scaff = $1;
	    $opera_contig = $contig_name{$opera_scaff};
	    #print STDERR " *** $opera_contig\n";<STDIN>;
	    #$opera_contig = ">".$1;
	    $contig_to_illumina_assembly->{$opera_contig} = [];
	}
	else{
	    $illumina_contig_length = $line[2];
	    #print STDERR " *** $illumina_contig_length\n";<STDIN>;
	    push(@{$contig_to_illumina_assembly->{$opera_contig}}, $illumina_contig_length);
	    #print STDERR " *** $opera_contig $illumina_contig_length " . "[" . join(" " , @{$contig_to_illumina_assembly->{$opera_contig}}). "]\n";<STDIN>;
	}
    }
    close(FILE);
}

sub analyze_mash{
    my ($mash_file, $bin_info) = @_;
    print STDERR " Analyze mash bin comparison from $mash_file\n";
    open(FILE, "sort -k3,3 -n $mash_file |");
    while(<FILE>){
	@line = split(/\t/, $_);
	#Get the bin
	$str = $line[1];
	@str_tab = split(/\//, $str);
	$bin = $str_tab[@str_tab-1];
	@str_tab = split(/\./, $bin);
	$bin_id = join("\.", @str_tab[0..@str_tab-2]);
	#print " *** " . "@str_tab" . " " . $bin_id . "\n";<STDIN>;

	if($bin_info->{$bin_id}->{"SPECIES"} eq "NA"){
	    
	    $similarity = $line[2];

	    #Get the species
	    $str = $line[0];
	    @str_tab = split(/\//, $str);
	    $file = $str_tab[@str_tab-2];
	    @str_tab = split(/\_/, $file);
	    $species = $str_tab[0] . "_" . $str_tab[1];

	    next if(index($species, "multispecies") != -1);
	    $bin_info->{$bin_id}->{"SPECIES"} = $species;
	    $bin_info->{$bin_id}->{"SIMILARITY"} = $similarity;
	}
    }
    close(FILE);
	
}

sub analyze_kraken{
    my ($kraken_file, $bin_info) = @_;
    print STDERR " Analyze kraken abundance from $kraken_file\n";
    if(-e $kraken_file){
    }
}

sub analyze_checkm{
    my ($checkm_file, $bin_info) = @_;
    print STDERR " Analyze checkm evaluation in $checkm_file\n";
    if(-e $checkm_file){
	open(FILE, $checkm_file);
	<FILE>;
	<FILE>;
	<FILE>;
	while(<FILE>){
	    @line = split(/\s+/, $_);
	    last if(@line < 14);
	    $bin = $line[1];
	    #print STDERR " *** " . $bin . "\n";
	    $bin_info->{$bin}->{"COMLETENESS"}  = $line[7];
	    $bin_info->{$bin}->{"CONTAMINATION"}  = $line[8];
	}
	close(FILE);
    }
}

sub analyze_contig_info{

    my ($contig_info_file, $contig_to_bin, $bin_info, $illumina_contig_length) = @_;

    open(FILE, $contig_info_file);
    <FILE>;
    while(<FILE>){
	chop $_;
	my ($contig, $contig_size, $contig_short_read_cov, $contig_long_read_cov, $contig_cluster, $contig_species, $nb_strain, $ref) = split(/\t/, $_);
	
	#print STDERR " *** $contig\n";<STDIN>;
	$bin = "UN_BINNED";
	$bin = $contig_to_bin->{$contig} if(exists $contig_to_bin->{$contig});

	next if($contig_size < 500);
	#OPERA ASSEMBLY
	push(@{$bin_info->{$bin}->{"CONTIG_SIZE"}}, $contig_size);
	$bin_info->{$bin}->{"CLUSTER_SIZE"} += $contig_size;
	#ORIGINAL ILLUMINA ASSEMBLY
	push(@{$bin_info->{$bin}->{"SHORT_READ_CONTIG_SIZE"}}, @{$illumina_contig_length->{$contig}});
	$bin_info->{$bin}->{"SHORT_READ_ASSEMBLY_SIZE"} += sum($illumina_contig_length->{$contig});
	#print STDERR " *** $contig " . "[" . join(" " , @{$illumina_contig_length->{$contig}}). "]\n";<STDIN>;
	#
	$bin_info->{$bin}->{"SHORT_READ_COV"} +=  $contig_short_read_cov * $contig_size;
	$bin_info->{$bin}->{"LONG_READ_COV"} +=  $contig_long_read_cov * $contig_size;
	
	#print STDERR $bin . "\t" . $contig . "\t" . $contig_size . "\n";
	
    }
    close(FILE);

    my @sort_tab = ();my $assembly_size;
    foreach $bin (keys %{$bin_info}){
	$assembly_size = $bin_info->{$bin}->{"CLUSTER_SIZE"};

	#print STDERR $bin . "\t" . $assembly_size . "\n";<STDIN>;
	#next;
	$bin_info->{$bin}->{"SHORT_READ_COV"} = sprintf("%.2f", $bin_info->{$bin}->{"SHORT_READ_COV"} / $assembly_size);
	$bin_info->{$bin}->{"LONG_READ_COV"} = sprintf("%.2f", $bin_info->{$bin}->{"LONG_READ_COV"} / $assembly_size);
	
	@sort_tab = sort {$b <=> $a} @{$bin_info->{$bin}->{"CONTIG_SIZE"}};
	$n50_pos = compute_Nx(50, $assembly_size, \@sort_tab);
	$n90_pos = compute_Nx(90, $assembly_size, \@sort_tab);
	#
	$bin_info->{$bin}->{"LONGEST_CONTIG"} = $sort_tab[0];
	$bin_info->{$bin}->{"NB_CONTIG"} = @sort_tab+0;
	$bin_info->{$bin}->{"N50"} = $sort_tab[$n50_pos];
	$bin_info->{$bin}->{"L50"} = $n50_pos+1;
	$bin_info->{$bin}->{"N90"} = $sort_tab[$n90_pos];
	$bin_info->{$bin}->{"L90"} = $n90_pos+1;
	#
	#For ILLUMINA ONLY ASSEMBLY
	$assembly_size = $bin_info->{$bin}->{"SHORT_READ_ASSEMBLY_SIZE"};
	#print STDERR " *** assembly size $assembly_size\n";<STDIN>;
	@sort_tab = sort {$b <=> $a} @{$bin_info->{$bin}->{"SHORT_READ_CONTIG_SIZE"}};
	$n50_pos = compute_Nx(50, $assembly_size, \@sort_tab);
	$n90_pos = compute_Nx(90, $assembly_size, \@sort_tab);
	#
	$bin_info->{$bin}->{"SHORT_READ_LONGEST_CONTIG"} = $sort_tab[0];
	$bin_info->{$bin}->{"SHORT_READ_NB_CONTIG"} = @sort_tab+0;
	$bin_info->{$bin}->{"SHORT_READ_N50"} = $sort_tab[$n50_pos];
	$bin_info->{$bin}->{"SHORT_READ_L50"} = $n50_pos+1;
	$bin_info->{$bin}->{"SHORT_READ_N90"} = $sort_tab[$n90_pos];
	#$bin_info->{$bin}->{"ILLUMINA_N90_CONTIG"} = $n90_pos+1;
	#
    }
}


sub analyze_bin{
    my ($bin_dir, $contig_to_bin, $bin_info) = @_;
    opendir(DIR, $bin_dir);
    my @all_bin = readdir(DIR);
    close(DIR);
    
    foreach $bin (@all_bin){
	if($bin =~ m /(.+)\.fa/){
	    $bin_id = $1;
	    #if(index($bin, ".fa") != -1){
	    $bin_info->{$bin_id} = init_bin();
	    open(FILE, "$bin_dir/$bin");
	    while(<FILE>){
		chop $_;
		$name = $_;
		if($name =~ m />(.+)/){
		    #print STDERR $name . " " . $bin . "\n";<STDIN>;
		    $name = $1;
		    #print STDERR " *** $name $bin_id\n";<STDIN>;
		    $contig_to_bin->{$name} = $bin_id;
		}
	    }
	    close(FILE);
	}
    }
}

sub init_bin{
    my %genome_init = ("SPECIES", "NA",
		       "SIMILARITY", "NA",
		       #
		       "CLUSTER_SIZE", 0,
		       "LONGEST_CONTIG", "NA",
		       "CONTIG_SIZE", [],
		       "NB_CONTIG", "NA",
		       "N50", "NA",
		       "L50", "NA",
		       "NG50", "NA",
		       "N90", "NA",
		       "L90", "NA",
		       "NG90", "NA",
		       #
		       "SHORT_READ_ASSEMBLY_SIZE", 0,
		       "SHORT_READ_LONGEST_CONTIG", "NA",
		       "SHORT_READ_CONTIG_SIZE", [],
		       "SHORT_READ_NB_CONTIG", NA,
		       "SHORT_READ_N50", "NA",
		       "SHORT_READ_L50", "NA",
		       #"ILLUMINA_N50_CONTIG", "NA",
		       #"ILLUMINA_NG50", "NA",
		       "SHORT_READ_N90", "NA",
		       #"ILLUMINA_N90_CONTIG", "NA",
		       #"ILLUMINA_NG90", "NA",
		       #
		       "SHORT_READ_COV", 0,
		       "LONG_READ_COV", 0,
		       #
		       "FRAQ_READ", "NA",
		       "COMLETENESS", "NA",
		       "CONTAMINATION", "NA",
		       "KRAKEN_SPECIES_ABUNDANCE", "NA",
		       "KRAKEN_GENUS_ABUNDANCE", "NA"
	);
    return \%genome_init;
}



sub compute_Nx{
    my ($x, $length, $ptr_tab) = @_;

    my $nb_x_base = $length * ($x / 100);
    
    my $total = 0;
    
    my $i = 0;
    for($i = 0; $total < $nb_x_base; $i++){
	$total += $ptr_tab->[$i];
    }
    
    return $i -1;
}

sub sum{
    my ($array) = @_;
    my $res = 0;
    foreach $s (@{$array}){
	$res += $s;
    }
    return $res;
}




