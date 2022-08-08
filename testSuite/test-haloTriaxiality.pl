#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use Galacticus::Options;

# Test the halo triaxiality model of Menker & Benson (2022) by computing the median axis ratios.
# Andrew Benson (06-June-2022)

# Run the model.
system("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/haloTriaxialityMenkerBenson2022.xml");
unless ( $? == 0 ) {
    print "FAIL: Menker & Benson (2022) halo triaxility failed to run\n";
    exit;
}

# Extract the data.
my $model      = new PDL::IO::HDF5("outputs/haloTriaxialityMenkerBenson2022.hdf5");
my $nodes      = $model->group  ('Outputs/Output1/nodeData'    )       ;
my $axisRatio2 = $nodes->dataset('darkMatterProfileAxisRatiosY')->get();
my $axisRatio3 = $nodes->dataset('darkMatterProfileAxisRatiosZ')->get();

# Get the median axis ratios.
my $order2           = $axisRatio2->qsorti();
my $order3           = $axisRatio3->qsorti();
my $medianAxisRatio2 = $axisRatio2->($order2)->((int(nelem($axisRatio2)/2)));
my $medianAxisRatio3 = $axisRatio3->($order3)->((int(nelem($axisRatio3)/2)));

# Validate the medians.
my $tolerance              = 0.05;
my $medianAxisRatio2Target = 0.73;
my $medianAxisRatio3Target = 0.53;
if ( abs($medianAxisRatio2-$medianAxisRatio2Target) < $tolerance ) {
    print "success: axis ratio 2 median\n";
} else {
    print "FAIL: axis ratio 2 median (".$medianAxisRatio2.")\n";
}
if ( abs($medianAxisRatio3-$medianAxisRatio3Target) < $tolerance ) {
    print "success: axis ratio 3 median\n";
} else {
    print "FAIL: axis ratio 3 median (".$medianAxisRatio3.")\n";
}

exit;
