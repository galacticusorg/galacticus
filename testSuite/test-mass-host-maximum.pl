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

# Check calculations of maximum host halo mass.
# Andrew Benson (26-July-2021)

# Run the model.
&System::Redirect::tofile("cd ..; ./Galacticus.exe testSuite/parameters/massHostMaximum.xml","outputs/massHostMaximum.log");
unless ( $? == 0 ) {
    print "FAILED:  model run:\n";
    system("cat outputs/massHostMaximum.log");
} else {
    print "SUCCESS: model run\n";
}
# Read the model data and check for consistency.
my $model   = new PDL::IO::HDF5("outputs/massHostMaximum.hdf5");
my $outputs = $model->group('Outputs');
foreach my $outputName ( $outputs->groups() ) {
    my $output            = $outputs ->group  ($outputName            )       ;
    (my $expansionFactor) = $output  ->attrGet('outputExpansionFactor')       ;
    my $nodeData          = $output  ->group  ('nodeData'             )       ;
    my $isIsolated        = $nodeData->dataset('nodeIsIsolated'       )->get();
    my $massHost          = $nodeData->dataset('satelliteHostHaloMass')->get();
    my $massHostMaximum   = $nodeData->dataset('massHostMaximum'      )->get();
    my $satellites        = which($isIsolated == 0);
    my $redshift          = sprintf("%3.1f",1.0/$expansionFactor-1.0);
    my $status = all($massHost->($satellites) == $massHostMaximum->($satellites)) ? "SUCCESS" : "FAILED";
    print $status.": maximum host mass at z=".$redshift."\n";
}

exit 0;
