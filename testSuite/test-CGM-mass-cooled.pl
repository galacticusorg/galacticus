#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use File::Slurp qw(slurp);
use System::Redirect;

# Check calculations of mass cooled out of the CGM.
# Andrew Benson (30-July-2021)

# This model has a cooling CGM hot halo, but no star formation or feedback. So, the mass cooled out of the CGM should equal
# the initial mass minus the final mass.
system("cd ..; ./Galacticus.exe testSuite/parameters/cgmMassCooled.xml");
unless ( $? == 0 ) {
    print "FAILED: model run:\n";
    exit 1;
} else {
    print "SUCCESS: model run\n";
}

# Read the model data and check for consistency.
my $model            = new PDL::IO::HDF5("outputs/cgmMassCooled.hdf5");
my $outputs          = $model   ->group  ('Outputs'                   )       ;
my $output           = $outputs ->group  ('Output1'                   )       ;
my $nodeData         = $output  ->group  ('nodeData'                  )       ;
my $massCGM          = $nodeData->dataset('hotHaloMass'               )->get();
my $massCooled       = $nodeData->dataset('massCooledCGM'             )->get();
my $massCooledTarget = 1.0e12-$massCGM->((0));
if ( abs($massCooledTarget-$massCooled->((0))) < 1.0e6 ) {
    print "SUCCESS: mass cooled out of CGM\n";
} else {
    print "FAILED: mass cooled out of CGM: ".$massCooled->((0))." ≇ ".$massCooledTarget."\n";
}

exit 0;
