#use strict;

my $out_file = $ARGV[0];

my @all_config = @ARGV[1..@ARGV-1];

my %processed_sample = ();
open(FILE, $out_file);
while(<FILE>){
    @line = split(/\t/, $_);
    $sample = $line[0];
    $processed_sample{$sample} = 1;
}
close(FILE);
open(OUT_INFO, ">>$out_file");
if(keys(%processed_sample) == 0){
    print OUT_INFO "Sample_path\tshort_read(millon)\tlong_read(millon)\tlong_read_throughput(Gbp)\tlong_read_N50(kbp)\n";
}

#my $head_cmd = "| head -n1000";
#$head_cmd = "";


my @long_read_size = ();
my $total_long_read_size;
foreach $c_file (@all_config){
    $tmp = `grep OUTPUT_DIR $c_file | cut -d ' ' -f2`;chop $tmp;
    @tmp_tab = split(/\//, $tmp);
    $sample = join("\/", @tmp_tab[0..(@tmp_tab-1)]);
    
    if(! exists $processed_sample{$sample}){
	print " *** Compute throughput for sample $sample " . "\n";
	$lib = `grep ILLUMINA_READ_1 $c_file | cut -d ' ' -f2`;chop $lib;
	print " *** Analyse Illumina library " . "\n";
	$nb_illumina_read = `zcat $lib $head_cmd | wc -l`;chop $nb_illumina_read;
	$nb_illumina_read = $nb_illumina_read / 2;

	print " *** Analyse long read library " . "\n";
	$lib = `grep -w LONG_READ $c_file | cut -f2 -d ' ' `;chop $lib;
	print STDERR " *** $lib\n";
	if(! -e $lib){
	    $lib = "$lib.gz";
	}
	
	@long_read_size = ();
	compute_read_size($lib, \@long_read_size);
	$total_long_read_size = sum(\@long_read_size);
	@sort_tab = sort {$b <=> $a} @long_read_size;
	$n50_pos = compute_Nx(50, $total_long_read_size, \@sort_tab);
	$nb_long_read = @long_read_size+0;
	print  OUT_INFO $sample . "\t" . sprintf("%.2f", ($nb_illumina_read/1000000)) . "\t" . sprintf("%.2f", ($nb_long_read/1000000)) . "\t" . sprintf("%.2f", ($total_long_read_size/1000000000)) . "\t" . sprintf("%.1f", ($sort_tab[$n50_pos]/1000)) . "\n";
    }
}
close(OUT_INFO);

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


sub compute_read_size{
    my ($read_file, $long_read_size) = @_;
    my $read_cmd = "$read_file";
    $read_cmd = "zcat $read_file |" if(index($read_file, ".gz") != -1);
    print STDERR $read_cmd . "\n";
    open(IN, "$read_cmd");
    #my $size_threshold = 200;
    my $comp = 0;
    my $n = 0;
    
    my $contig_name = "";
    my $contig_size = 0;
    my $count=1;
    while(<IN>){
	
	chomp $_;
	
	if($_ ne ""){
	    
	    if($count%4 == 1){
		$contig_size = ($contig_size-1)/2;
		if($count !=1){
		    push(@{$long_read_size}, $contig_size);
		}
		$contig_name = $_;
		$contig_size = 0;
	    }
	    else{
		$n = length($_);
		$contig_size += $n; 
	    }
	}
	$count++;
    }
    close(IN);
}

sub sum{
    my ($array) = @_;
    my $res = 0;
    foreach $s (@{$array}){
	$res += $s;
    }
    return $res;
}


