#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;

# Test that the determinstic spin model matches the constraint.
# Andrew Benson (15-January-2020)

# Run model.
system("cd ..; mkdir -p testSuite/outputs; Galacticus.exe testSuite/parameters/constrainDeterministicSpins.xml");
unless ( $? == 0 ) {
    print "FAILED: Galacticus model failed to run\n";
    exit 0;
}

# Extract the log-likelihood and check it is sufficiently high.
my  $model          = new PDL::IO::HDF5("outputs/constrainDeterministicSpins.hdf5");
my  $analyses       = $model   ->group  ('analyses'                );
my  $analysis       = $analyses->group  ('spinDistributionBett2007');
(my $logLikelihood) = $analysis->attrGet('logLikelihood'           );
if ( $logLikelihood > -303.0 ) {
    print "SUCCESS: deterministic spin model matches constraint [log ℒ = "      .sprintf("%10.2f",$logLikelihood->sclr())."]\n";
} else {
    print "FAILED: deterministic spin model does not match constraint [log ℒ = ".sprintf("%10.2f",$logLikelihood->sclr())."]\n";
}

exit 0;
