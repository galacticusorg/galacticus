#!/usr/bin/env bash

# Find the maximum likelihood estimate of the covariance matrix for the Baldry et al. (2012) GAMA stellar mass function.
# Andrew Benson (11-June-2014)

# Simply run the generic script with our config file as argument.
"${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/scripts/covarianceMatrix.pl" "${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/stellarMassFunction_GAMA_z0.03/covarianceMatrixControl.xml"
