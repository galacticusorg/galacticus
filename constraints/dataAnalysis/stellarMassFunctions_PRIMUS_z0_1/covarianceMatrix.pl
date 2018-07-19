#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";

# Find the maximum likelihood estimate of the covariance matrix for the Moustakas et al. (2013) PRIMUS stellar mass functions.
# Andrew Benson (18-May-2014)

# Get the argument specifying which redshift bin to use.
die("Usage: covarianceMatrix.pl <redshiftBin>")
    unless ( scalar(@ARGV) == 1 );
my $redshiftBin = $ARGV[0];
die("covarianceMatrix.pl: redshiftBin must be 0, 1, 2, 3, 4, 5 or 6")
    if ( $redshiftBin < 0 || $redshiftBin > 6 );

# Simply run the generic script with our config file as argument.
system(&galacticusPath()."constraints/dataAnalysis/scripts/covarianceMatrix.pl ".&galacticusPath()."constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/covarianceMatrixControl.xml ".$redshiftBin);

exit;
