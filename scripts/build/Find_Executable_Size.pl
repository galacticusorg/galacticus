#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 

die('Usage: Find_Executable_Size.pl Executable Size_File')
    unless ( scalar(@ARGV) == 2 );
my $Executable = $ARGV[0];
my $Size_File  = $ARGV[1];

system("size $Executable > $Size_File");

exit
