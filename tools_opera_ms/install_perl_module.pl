#!/usr/bin/perl

#use strict;
use warnings;

#Check if cpanm is installed
#my $path_to_cpanm = `command -v cpanm`;chop $path_to_cpanm;
`cpanm --version 2> /dev/null`;
if ($?){
    die("Please install cpanm for the setup the Perl modules");
}

#List of perl modules to install
check_module("Switch");
check_module("File::Which");
check_module("File::Spec::Functions");
check_module("Statistics::Basic");
check_module("Statistics::R");
check_module("Getopt::Long");


sub check_module{
    my ($module) = @_;
    print STDERR " *** Check perl module $module\n";
    system("/usr/bin/perl -e 'use $module ' 2> tmp_module.str");
    open(F, "tmp_module.str");
    while(<F>){
	if($_ ne ""){
	    print STDERR " *** Install perl module $module\n";
	    print `cpanm --force $module`;
	    last;
	}
    }
    `rm tmp_module.str`;
}
