#!/usr/bin/perl
use warnings "all";

($seq_file, $min_contig_size, $min_scaffold_size) = @ARGV;

if(index($seq_file, ".gz") != -1){
    open(FILE, "zcat $seq_file |");
}
else{
    open(FILE, "$seq_file");
}

#open(FILE, "$seq_file");

$min_contig_size = $min_scaffold_size = 500 if(@ARGV == 1);

my @scaffold = (0);
my $nb_scaff = 0;
my $nb_no_N_bp = 0;

my $current_nb_N = 0;
my $total_nb_N = 0;

my @contig = (0);
my $nb_contig = 0;

my @size_stats = (0, 100, 500, 1000, 2000, 10000, 50000, 100000, 500000, 1000000);
my @nb_contig_size;
my @nb_scaff_size;
foreach $s (@size_stats){
    my @tab = (0, 0);
    push(@nb_contig_size, \@tab);
    my @tab1 = (0, 0, 0);
    push(@nb_scaff_size, \@tab1);
}


while(<FILE>){
    
    #print "\t$_ -----> $nb_scaff\n";<STDIN>;

    if(index($_, ">") == -1){

	my $seq = $_; chomp $seq; 
	my $total_bp = length $seq; 

	my @contigs = split(/N+/i, $seq);

	my $i = 0;
	for(; $i <= $#contigs; $i++) {

	    $contig[$nb_contig] += length $contigs[$i];

	    if($#contigs > 0) {

		if($contig[$nb_contig] >= $min_contig_size){
			$nb_contig++;
		}

		$contig[$nb_contig] = 0;		    
	    }
	}

	$seq =~ s/N//gi;
	my $nonN_bp = length $seq;
	$current_nb_N += $total_bp-$nonN_bp;
	$scaffold[$nb_scaff] += $nonN_bp;
	#print "$nb_scaff -> $nb_no_N_bp\n";<STDIN>;
    }
    
    else{

	for($i = 0; $i < @size_stats; $i++){
	    if($scaffold[$nb_scaff] >= $size_stats[$i]){
		$nb_scaff_size[$i]->[0]++;
		$nb_scaff_size[$i]->[1] += $scaffold[$nb_scaff];
		$nb_scaff_size[$i]->[2] += $current_nb_N;
	    }
	}

	#A new scaffold
	if($scaffold[$nb_scaff] >= $min_scaffold_size){
	    $nb_scaff++;
	    $total_nb_N += $current_nb_N;
	}
	$scaffold[$nb_scaff] = 0;
	$current_nb_N = 0;
	
	#a new contig inside at the end of a scaffold
	if($contig[$nb_contig] >= $min_contig_size){
	    $nb_contig++;
	}
	$contig[$nb_contig] = 0;
	
    }
}
$total_nb_N += $current_nb_N;

pop(@scaffold) if($scaffold[@scaffold-1] < $min_scaffold_size);
pop(@contig) if($contig[@contig-1] < $min_contig_size);

my @sort_tab = ();

for($i = 0; $i < 2; $i ++){
    if($i == 0){
	@sort_tab = sort {$b <=> $a} @contig;
	#nb scaff
	print "\nNo. of contig:\t".(@sort_tab)."\n";
    }
    else{
	@sort_tab = sort {$b <=> $a} @scaffold;
	#nb scaff
	print "\nNo. of scaffolds:\t".(@sort_tab)."\n";
	print "gap-size: $total_nb_N\n";
    }

    #Total length
    my $total_length = 0;
    foreach $s (@sort_tab){
	$total_length += $s;
    }
    print "Total length:\t$total_length\n";
    
    #AVG
    print "Average length:\t".int($total_length/(@sort_tab+0))."\n";

    #Max length
    print "Max length:\t$sort_tab[0]\n";

    #Min length
    print "Min length:\t".$sort_tab[@sort_tab-1]."\n";

    #N50
    $nb_seq = compute_Nx(50, $total_length, \@sort_tab);
    print "N50 length:\t".$sort_tab[$nb_seq]."\t".($nb_seq+1)."\n";

    #N90
    $nb_seq = compute_Nx(90, $total_length, \@sort_tab);
    print "N90 length:\t".$sort_tab[$nb_seq]."\t".($nb_seq+1)."\n";
    #print "N90 length:\t".."\n";    
    if($i == 1){
	for(my $j = 0; $j < @size_stats; $j++){
	    $assembly_length_maxsize = $nb_scaff_size[$j]->[1];
	    $a_l_s = $assembly_length_maxsize;
	    $assembly_length_maxsize = $genome_size if(defined $genome_size);
	    #
	    #$nb_seq = compute_Nx(50, $assembly_length_maxsize, \@sort_tab);
	    $N50 = "-";#$sort_tab[$nb_seq];
	    if($N50 >= $size_stats[$j]){
		print "Nb $type >= $size_stats[$j]bp $N50  $nb_scaff_size[$j]->[0] $a_l_s $nb_scaff_size[$j]->[2]\t(".($nb_seq+1).")\n";
	    }
	    else{
		print "Nb $type >= $size_stats[$j]bp NA  $nb_scaff_size[$j]->[0] $a_l_s $nb_scaff_size[$j]->[2]\t(NA)\n";
	    }
	}
    }

}





sub compute_Nx{
    my ($x, $length, $ptr_tab) = @_;

    my $nb_x_base = $length * ($x / 100);
    
    my $total = 0;
    
    my $i = 0;
    for($i = 0; $total < $nb_x_base; $i++){
	$total += $ptr_tab->[$i];
    }
    
    #print "$x $pos_x $size\n";<STDIN>;

    #return $ptr_tab->[$i-1];
    return $i -1;
}

