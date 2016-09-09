#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/../perl';
use File::Which;

# Export trees from Galacticus and check that they are written correctly.
# Andrew Benson (12-October-2012)

# Run models.
system("cd ..; mkdir -p testSuite/outputs/test-merger-tree-write; scripts/aux/launch.pl testSuite/parameters/test-merger-tree-write.xml");

# Validate the IRATE-format output.
my $validator = &File::Which::which('iratevalidate');
if ( $validator ) {
    system("iratevalidate outputs/test-merger-tree-write/exportedTreesIRATE.hdf5");
    die("FAILED: IRATE-format file ouput by Galacticus did not validate")
	unless ( $? == 0 );
} else {
    print "SKIP: iratevalidate is not installed - validation of IRATE-format file will be skipped\n";
}

exit;
