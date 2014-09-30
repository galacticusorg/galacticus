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

# Find the maximum likelihood estimate of the covariance matrix for the Bernardi et al. (2013) SDSS stellar mass function.
# Andrew Benson (12-May-2014)

# Simply run the generic script with our config file as argument.
system($galacticusPath."constraints/dataAnalysis/scripts/covarianceMatrix.pl ".$galacticusPath."constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/covarianceMatrixControl.xml");

exit;
