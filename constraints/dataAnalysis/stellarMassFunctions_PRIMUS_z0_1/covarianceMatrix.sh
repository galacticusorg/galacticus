#!/usr/bin/env bash

# Find the maximum likelihood estimate of the covariance matrix for the Moustakas et al. (2013) PRIMUS stellar mass functions.
# Andrew Benson (18-May-2014)

# Get the argument specifying which redshift bin to use.
if [ "$#" -ne 1 ]; then
    echo "Usage: covarianceMatrix.sh <redshiftBin>" >&2
    exit 1
fi
redshiftBin="$1"
if [ "${redshiftBin}" -lt 0 ] || [ "${redshiftBin}" -gt 6 ]; then
    echo "covarianceMatrix.sh: redshiftBin must be 0, 1, 2, 3, 4, 5 or 6" >&2
    exit 1
fi

# Simply run the generic script with our config file as argument.
"${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/scripts/covarianceMatrix.pl" "${GALACTICUS_EXEC_PATH}/constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/covarianceMatrixControl.xml" "${redshiftBin}"
