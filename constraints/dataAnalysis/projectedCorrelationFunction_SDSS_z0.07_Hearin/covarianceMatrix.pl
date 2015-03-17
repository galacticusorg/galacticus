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

# Find the maximum likelihood estimate of the covariance matrix for the Hearin et al. (2014)
# SDSS stellar mass selected projected correlation functions.
# Andrew Benson (22-July-2014)

# Simply run the generic script with our config file as argument.
system($galacticusPath."constraints/dataAnalysis/scripts/covarianceMatrixProjectedCorrelation.pl ".$galacticusPath."constraints/dataAnalysis/projectedCorrelationFunction_SDSS_z0.07_Hearin/covarianceMatrixControl.xml");

exit;
