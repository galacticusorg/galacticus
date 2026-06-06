#!/usr/bin/env bash

# Find the maximum likelihood estimate of the covariance matrix for the Muzzin et al. (2013) ULTRAVISTA stellar mass functions.
# Andrew Benson (12-August-2014)

# Get the argument specifying which redshift bin to use.
if [ "$#" -ne 1 ]; then
    echo "Usage: covarianceMatrix.sh <redshiftBin>" >&2
    exit 1
fi
redshiftBin="$1"
if [ "${redshiftBin}" -lt 0 ] || [ "${redshiftBin}" -gt 6 ]; then
    echo "covarianceMatrix.sh: redshiftBin must be 0, 1, 2, 3, 4, 5, or 6" >&2
    exit 1
fi

# Simply run the generic script with our config file as argument.
"${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/scripts/covarianceMatrix.pl" "${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/covarianceMatrixControl.xml" "${redshiftBin}"
