#!/usr/bin/perl
use warnings;



sub write_time{
    my ($out_dir, $step_name, $time) = @_;
    open(OUT, ">>$out_dir/time.log");
    print OUT $step_name . "\t" . $time . "\n";
    close(OUT);
}

1;
