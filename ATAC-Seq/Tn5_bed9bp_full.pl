#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Wed Apr  9 14:58:03 2014
###################################
use strict;
use warnings;

my $usage=<<USAGE;
Usage: preShift.pl in.bed out.bed
in.bed  -- Your input alignment in bed format
out.bed -- Your output file
Example: preShift in.bed out_preShift75.bed
USAGE

if(@ARGV !=2){
    print $usage;
    exit 1;
}

my ($in,$out) = @ARGV;

if(! -e "$in") {
    print "input file '$in' doesn't exists\n";
    print $usage;
    exit 1;
}

open(IN,"$in") or die $!;
open(OUT,">$out") or die $!;
while(<IN>){
    chomp;
    my($chr,$start,$end,$name,$score,$strand) = split "\t";
    if($strand eq "+"){
        $start = $start;
        $end = $start + 9;
    }else{
        $end = $end;
        $start = $end - 9;
    }
    $start = 0 if ($start <0);
    print OUT join "\t",($chr,$start,$end,$name,$score,$strand."\n");
}
close IN;
close OUT;