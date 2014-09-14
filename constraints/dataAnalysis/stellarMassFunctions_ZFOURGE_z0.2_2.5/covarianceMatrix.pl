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

# Find the maximum likelihood estimate of the covariance matrix for the Tomczak et al. (2014) ZFOURGE stellar mass functions.
# Andrew Benson (11-August-2014)

# Get the argument specifying which redshift bin to use.
die("Usage: covarianceMatrix.pl <redshiftBin>")
    unless ( scalar(@ARGV) == 1 );
my $redshiftBin = $ARGV[0];
die("covarianceMatrix.pl: redshiftBin must be 0, 1, 2, 3, 4, 5, 6, or 7")
    if ( $redshiftBin < 0 || $redshiftBin > 7 );

# Simply run the generic script with our config file as argument.
system($galacticusPath."constraints/dataAnalysis/scripts/covarianceMatrix.pl ".$galacticusPath."constraints/dataAnalysis/stellarMassFunctions_ZFOURGE_z0.2_2.5/covarianceMatrixControl.xml ".$redshiftBin);

exit;
