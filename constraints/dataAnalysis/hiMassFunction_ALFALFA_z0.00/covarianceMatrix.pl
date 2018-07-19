#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";

# Find the maximum likelihood estimate of the covariance matrix for the Martin et al. (2011) ALFALFA HI mass function.
# Andrew Benson (05-July-2012)

# Simply run the generic script with our config file as argument.
system(&galacticusPath()."constraints/dataAnalysis/scripts/covarianceMatrix.pl ".&galacticusPath()."constraints/dataAnalysis/hiMassFunction_ALFALFA_z0.00/covarianceMatrixControl.xml");

exit;
