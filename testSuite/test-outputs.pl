#!/usr/bin/env perl
use strict;
use warnings;

# Run a set of short Galacticus models exploring a range of output options from Galacticus to ensure
# that they at least run to completion.
# Andrew Benson (26-Jan-2011)

# Simply run the models.
system("cd ..; scripts/aux/launch.pl testSuite/test-outputs.xml");

# Check for failed models.
system("grep -q -i fatal outputs/test-outputs/galacticus_*/galacticus.log");
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i fatal outputs/test-outputs/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS!\n";
}

exit;
