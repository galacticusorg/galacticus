#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use Galacticus::HDF5;
use Galacticus::Inclination;

# Run a test case for setting properties back to the HDF5 output file.
# Andrew Benson (15-October-2014)

# Run the model and check for successful completion.
system("./Galacticus.exe testSuite/parameters/setProperties.xml");
die("FAILED: setProperties.pl model failed to complete") 
    unless ( $? == 0 );

# Check that the spheroid size is non-zero whenever the spheroid mass is non-zero.
my $dataSet;
$dataSet->{'file'  } = "testSuite/outputs/test-set-properties.hdf5";
$dataSet->{'store' } = 1;
&Galacticus::HDF5::Get_Parameters($dataSet);
&Galacticus::HDF5::Count_Trees($dataSet);
&Galacticus::HDF5::Select_Output($dataSet,0.0);
$dataSet->{'tree'} ="all";
eval { &Galacticus::HDF5::Get_Dataset($dataSet,['inclination']) };
if ( $@ ) {
    print "FAILED: setProperties.pl unable to set property in HDF5 file\n";
}

exit;
