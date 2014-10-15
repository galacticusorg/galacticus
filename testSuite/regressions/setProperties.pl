#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
require Galacticus::HDF5;
require Galacticus::Inclination;

# Run a test case for setting properties back to the HDF5 output file.
# Andrew Benson (15-October-2014)

# Run the model and check for successful completion.
system("Galacticus.exe testSuite/parameters/setProperties.xml");
die("FAILED: setProperties.pl model failed to complete") 
    unless ( $? == 0 );

# Check that the spheroid size is non-zero whenever the spheroid mass is non-zero.
my $dataSet;
$dataSet->{'file'  } = "testSuite/outputs/test-set-properties.hdf5";
$dataSet->{'store' } = 1;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);
$dataSet->{'tree'} ="all";
eval { &HDF5::Get_Dataset($dataSet,['inclination']) };
if ( $@ ) {
    print "FAILED: setProperties.pl unable to set property in HDF5 file\n";
}

exit;
