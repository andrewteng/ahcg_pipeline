#! /usr/bin/perl
use warnings;
use strict;

open (BED_FILE, "breast_ovarian.txt");

foreach (<BED_FILE>){
	my @start;
	my @stop;
	my @t;
	@t = split("\t",$_);
        @start = split(",",$t[9]);
        @stop = split(",",$t[10]);
        my $i = 0;
	my @new_start = map { $_ - 20 } @start;
	my @new_stop = map { $_ + 20} @stop;
	while ($i < @start){
		
        	print "$t[2]\t$new_start[$i]\t$new_stop[$i]\t$t[1]\t$t[8]\t$t[3]\n";
        	$i++;
	}
}


close BED_FILE;
exit;
