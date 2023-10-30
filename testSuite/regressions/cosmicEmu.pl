#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;

# Run a power spectrum task using the CosmicEmu nonlinear power spectrum twice. This tests that a pre-computed power spectrum file
# is correctly re-read on subsequent runs.
# Andrew Benson (13-October-2023)

# Run the model and check for successful completion.
for(my $i=0;$i<2;++$i) {
    system("./Galacticus.exe testSuite/regressions/cosmicEmu.xml");
    die("FAILED: cosmicEmu.pl model failed to complete") 
	unless ( $? == 0 );
}

exit;
