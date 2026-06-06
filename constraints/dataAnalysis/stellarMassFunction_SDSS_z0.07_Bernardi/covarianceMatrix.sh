#!/usr/bin/env bash

# Find the maximum likelihood estimate of the covariance matrix for the Bernardi et al. (2013) SDSS stellar mass function.
# Andrew Benson (12-May-2014)

# Simply run the generic script with our config file as argument.
"${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/scripts/covarianceMatrix.pl" "${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/covarianceMatrixControl.xml"
