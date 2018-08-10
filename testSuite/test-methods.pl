#!/usr/bin/env perl
use strict;
use warnings;

# Run a set of short Galacticus models spanning a full range of method options to ensure
# that they at least run to completion.
# Andrew Benson (04-Sep-2010)

# Remove automatically generated files to force them to be regenerated (and, therefore, to test that the generating code works
# correctly).
#
# FSPS stellar population synthesis code and associated file.
system("rm -f ".$ENV{'GALACTICUS_DATA_PATH'}."/dynamic/stellarPopulations/SSP_Spectra_Conroy-et-al_v2.5_imfSalpeter.hdf5");
system("rm -rf ../aux/FSPS_v2.4");
# Core files (older than 7 days).
system("find ../ -name '".$_."' -ctime +7 -exec rm {} \\;")
    foreach ( "core.*", "vgcore.*" );
# Noninstantaneous recycling files (older than 14 days).
system("find ".$ENV{'GALACTICUS_DATA_PATH'}."/dynamic/stellarPopulations -name '".$_."' -ctime +14 -exec rm {} \\;")
    foreach ( "Stellar_*_Yield_*_*.xml", "Stellar_Recycled_Fraction_*_*.xml", "Stellar_Energy_Input_*_*.xml" );
# CAMB transfer function files (older than 14 days).
system("find ".$ENV{'GALACTICUS_DATA_PATH'}."/dynamic/largeScaleStructure -name 'transfer_function_CAMB_*.xml' -ctime +14 -exec rm {} \\;");

# Simply run the models.
system("cd ..; scripts/aux/launch.pl testSuite/test-methods.xml");

# Check for failed models.
system("grep -q -i fatal outputs/test-methods/galacticus_*/galacticus.log");
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i fatal outputs/test-methods/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS!\n";
}

exit;
