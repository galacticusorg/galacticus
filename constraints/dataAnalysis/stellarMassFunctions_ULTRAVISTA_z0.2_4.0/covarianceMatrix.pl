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

# Find the maximum likelihood estimate of the covariance matrix for the Muzzin et al. (2013) ULTRAVISTA stellar mass functions.
# Andrew Benson (12-August-2014)

# Get the argument specifying which redshift bin to use.
die("Usage: covarianceMatrix.pl <redshiftBin>")
    unless ( scalar(@ARGV) == 1 );
my $redshiftBin = $ARGV[0];
die("covarianceMatrix.pl: redshiftBin must be 0, 1, 2, 3, 4, 5, or 6")
    if ( $redshiftBin < 0 || $redshiftBin > 6 );

# Simply run the generic script with our config file as argument.
system($galacticusPath."constraints/dataAnalysis/scripts/covarianceMatrix.pl ".$galacticusPath."constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/covarianceMatrixControl.xml ".$redshiftBin);

exit;
