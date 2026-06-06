#!/usr/bin/env bash

# Find the maximum likelihood estimate of the covariance matrix for the Li & White (2009) SDSS stellar mass function.
# Andrew Benson (05-July-2012)

# Simply run the generic script with our config file as argument.
"${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/scripts/covarianceMatrix.pl" "${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/covarianceMatrixControl.xml"
