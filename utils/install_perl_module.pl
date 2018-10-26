#!/usr/bin/perl

#use strict;
use warnings;

check_module("Switch");
check_module("Statistics::Basic");
check_module("File::Which");
check_module("Statistics::R");

sub check_module{
    my ($module) = @_;
    print STDERR " *** Check perl module $module\n";
    system("perl -e 'use $module ' 2> tmp_module.str");
    open(F, "tmp_module.str");
    while(<F>){
	if($_ ne ""){
	    print STDERR " *** Install perl module $module\n";
	    print `cpanm --force $module`;
	    `rm tmp_module.str`;
	    last;
	}
    }
    
}
