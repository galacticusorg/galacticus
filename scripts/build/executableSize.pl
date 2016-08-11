#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 

# Determine the size of an executable file and store the results to a given file.
# Andrew Benson (08-August-2016)

die('Usage: executableSize.pl <executableName> <sizeFileName>')
    unless ( scalar(@ARGV) == 2 );
my $executableName = $ARGV[0];
my $sizeFileName   = $ARGV[1];

# Simply run the "size" command and redirect output to the given file.
system("size ".$executableName." > ".$sizeFileName);

exit;
