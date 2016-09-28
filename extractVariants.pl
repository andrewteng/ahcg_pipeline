#! /usr/bin/perl
use warnings;
use strict;

open (BED_FILE, "exomes_breast_ovarian.bed");
open (VCF_FILE, "variants.vcf");

my @chr;
my @start;
my @stop;
my @column;


foreach (<BED_FILE>){
	@column = split ("\t", $_);
	push @chr, $column[0];
	push @start, $column[1];
	push @stop, $column[2];
}

foreach (<VCF_FILE>){
	next if /^\s*#/;
	my @column1 = split ("\t", $_);

	if ($column1[0] = @chr){
		if ($column1[1] >= @start){
			if ($column1[1] <= @stop){
				print $_;
			}
		}		
	}
} 

close BED_FILE;
close VCF_FILE;

exit;
