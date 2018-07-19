#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";

# Find the maximum likelihood estimate of the covariance matrix for the Bernardi et al. (2013) SDSS stellar mass function.
# Andrew Benson (12-May-2014)

# Simply run the generic script with our config file as argument.
system(&galacticusPath()."constraints/dataAnalysis/scripts/covarianceMatrix.pl ".&galacticusPath()."constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/covarianceMatrixControl.xml");

exit;
