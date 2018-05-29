#!/usr/bin/env perl
use strict;
use warnings;

# Run a set of short Galacticus models to explore different cosmological models.
# Andrew Benson (05-Sep-2010)

# Simply run the models.
system("cd ..; scripts/aux/launch.pl testSuite/test-cosmology.xml");

# Check for failed models.
system("grep -q -i fatal outputs/test-cosmology/galacticus_*/galacticus.log");
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i fatal outputs/test-cosmology/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS!\n";
}

exit;
