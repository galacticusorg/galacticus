#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;

# Run a single Galacticus model and ensure that outputs are created.
# Andrew Benson (04-April-2014)

# Simply run the models.
system("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/test-output-times.xml");

# Check for outputs.
die("test-output-times.pl: FAILED to run Galacticus model")
    unless ( -e "outputs/test-output-times.hdf5" );
# Find all output times.
my $file    = new PDL::IO::HDF5("outputs/test-output-times.hdf5");
my $outputs = $file->group('Outputs');
my $times = pdl [];
foreach my $outputName ( $outputs->groups() ) {
    my $output = $outputs->group($outputName);
    (my $outputTime) = $output->attrGet('outputTime');
    $times = $times->append($outputTime);
}
# Construct the expected times.
my $HubbleConstant = pdl 70.0;
my $megaParsec     = pdl 3.08528229e22;
my $kilo           = pdl 1.0e3;
my $gigaYear       = pdl 3.15576e16;
my $ageUniverse    = (2.0/3.0)*$megaParsec/$kilo/$gigaYear/$HubbleConstant;
# Fixed times.
my $timesExpected  = pdl [ 8.0, 9.0 ];
# Lookback times.
$timesExpected     = $timesExpected->append($ageUniverse-2.4);
# Redshifts.
my $redshifts      = pdl [ 0.0, 1.0, 2.0 ];
$timesExpected     = $timesExpected->append($ageUniverse/(1.0+$redshifts)**1.5);
# Sort the times.
$timesExpected     = $timesExpected->qsort();
# Compare the times.
if ( nelem($times) != nelem($timesExpected) ) {
    print "test-output-times.pl: FAILED - number of times does not match\n";
} elsif ( any(abs($times-$timesExpected)/$timesExpected > 1.0e-3) ) {
    print "test-output-times.pl: FAILED - times do not match\n";
} else {
    print "test-output-times.pl: SUCCESS\n";
}

exit;
