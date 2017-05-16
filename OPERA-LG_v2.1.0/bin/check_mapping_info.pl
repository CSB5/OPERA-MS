#!/usr/bin/perl
use warnings;
use Getopt::Long;

my ($edge_file_info, $c1, $c2, $read) = @ARGV;

open(FILE, $edge_file_info);
while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    $contig_1 = $line[0];
    $contig_2 = $line[2];

    if(($c1 eq $contig_1 || $c1 eq $contig_2) &&
       ($c2 eq $contig_1 || $c2 eq $contig_2)){
	
	if($read eq "ALL"){
	    print "$_"."\n";
	    last;
	}


	@tmp = split(/\:/, $line[5]);@coord_contig_1 = split(/\;/, $tmp[1]);
	@tmp = split(/\:/, $line[6]);@coord_contig_2 = split(/\;/, $tmp[1]);
	#
	@tmp = split(/\:/, $line[7]);@coord_contig_1_on_read = split(/\;/, $tmp[1]);
	@tmp = split(/\:/, $line[8]);@coord_contig_2_on_read = split(/\;/, $tmp[1]);
	#
	@tmp = split(/\:/, $line[9]);@contig_in_gap = split(/\;/, $tmp[1]);
	
	@read_name = split(/\;/, $line[4]);

	if(index($read, "_") == -1){
	    print STDERR 
		$read_name[$read]."\n".
		"COORD_CONTIG_1:".$coord_contig_1[$read]."\t"."COORD_CONTIG_2:".$coord_contig_2[$read]."\t".
		"COORD_CONTIG_1_ON_READ:".$coord_contig_1_on_read[$read]."\t"."COORD_CONTIG_2_ON_READ:".$coord_contig_2_on_read[$read]."\t".
		"CONTIG_FOR_GAPFILLING:".$contig_in_gap[$read].
		"\n";
	}
	else{
	    
	    for(my $i = 0; $i < @read_name; $i++){
		$r = $read_name[$i];
		if($r eq $read){
		    print STDERR 
			"COORD_CONTIG_1:".$coord_contig_1[$i]."\t"."COORD_CONTIG_2:".$coord_contig_2[$i]."\t".
			"COORD_CONTIG_1_ON_READ:".$coord_contig_1_on_read[$i]."\t"."COORD_CONTIG_2_ON_READ:".$coord_contig_2_on_read[$i]."\t".
			"CONTIG_FOR_GAPFILLING:".$contig_in_gap[$i].
			"\n";
		}
	    }
	}
    }
}

