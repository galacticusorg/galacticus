#!/usr/bin/env perl

unless ( $#ARGV == 1 ) {die('Usage: Find_Executable_Size.pl Executable Size_File')};
my $Executable = $ARGV[0];
my $Size_File  = $ARGV[1];

system("size $Executable > $Size_File");

exit
