"""Likelihood evaluation for halo mass function validation models.

This module computes the negative log-likelihood of halo mass function models
(computed by the Galacticus ``haloMassFunction`` task) with respect to N-body
halo mass function data, replicating the calculation performed by the
``posteriorSampleLikelihoodHaloMassFunction`` class
(``source/models/likelihoods/halo_mass_function.F90``) that was used in
calibrating the halo mass function model of Benson et al. (2026;
https://ui.adsabs.harvard.edu/abs/2026arXiv260612137B).

Specifically, for each mass bin of the N-body data within the allowed mass
range, the expected count of halos is

    mu_i = c phi_i,

where phi_i is the model mass function (per log-mass, averaged over the bin),
and c is a conversion factor between mass function and count inferred from the
data. The likelihood of the observed count N_i is then evaluated assuming a
negative binomial distribution - an over-dispersed Poisson distribution in
which the "stopping time" parameter r = 1/varianceFractionalModelDiscrepancy
accounts for fractional variance due to model discrepancy (a pure Poisson
distribution is used if the model discrepancy is zero).

The model mass function is read directly from the ``haloMassFunctionLnMBinAveraged``
dataset of the model output. The mass grid used by the validation models is
constructed to align exactly with the mass bins of the N-body data (the same
convention - bins spaced uniformly in log-mass, averaged between bin edges at
half the logarithmic spacing - is used by both the ``haloMassFunction`` task
and the likelihood class), so no interpolation is needed.
"""

import math

import h5py
import numpy as np

__all__ = ['negativeLogLikelihood', 'readModel', 'readData']

# Do not allow data bins to be matched to model grid points misaligned by more
# than this fraction of the bin width.
toleranceMisalignment = 1.0e-3


def readModel(fileNameModel):
    """Read a halo mass function model output file.

    Returns a dict with the model mass grid, bin-averaged mass function, and
    output redshift.
    """
    with h5py.File(fileNameModel, "r") as model:
        output = model["Outputs/Output1"]
        return {
            "mass"        : output["haloMass"                      ][:],
            "massFunction": output["haloMassFunctionLnMBinAveraged"][:],
            "redshift"    : float(output.attrs["outputRedshift"]),
        }


def readData(fileNameData):
    """Read an N-body halo mass function data file.

    Returns a dict with the bin masses, mass function, halo counts, and any
    environment metadata.
    """
    with h5py.File(fileNameData, "r") as data:
        simulation = data["simulation0001"]
        result = {
            "mass"        : simulation["mass"        ][:],
            "massFunction": simulation["massFunction"][:],
            "count"       : simulation["count"       ][:],
        }
        for attribute in ("massEnvironment", "overdensityEnvironment", "massPrimary"):
            if attribute in simulation.attrs:
                result[attribute] = float(simulation.attrs[attribute])
    return result


def negativeLogLikelihood(model, data, massRangeMinimum, massRangeMaximum=None,
                          varianceFractionalModelDiscrepancy=0.0):
    """Compute the negative log-likelihood of a halo mass function model.

    Parameters
    ----------
    model : dict
        Model mass function as returned by ``readModel()``.
    data : dict
        N-body mass function data as returned by ``readData()``.
    massRangeMinimum : float
        The minimum halo mass to include in the likelihood evaluation.
    massRangeMaximum : float, optional
        The maximum halo mass to include in the likelihood evaluation.
    varianceFractionalModelDiscrepancy : float, optional
        Fractional variance due to model discrepancy. If non-zero, a negative
        binomial distribution is used in place of a Poisson distribution.

    Returns
    -------
    tuple[float, dict]
        A 2-tuple ``(negLogLikelihood, bins)`` where ``bins`` is a dict of
        per-bin arrays (for the bins included in the likelihood): ``mass``,
        ``count`` (observed), ``countMean`` (model expectation),
        ``massFunctionData``, ``massFunctionModel``, and ``massFunctionError``
        (the root-variance of the count distribution, converted to mass
        function units).
    """
    if massRangeMaximum is None:
        massRangeMaximum = np.inf
    # Select data bins within the allowed mass range. Note that (following the
    # likelihood class) bins with zero counts are retained.
    selected = (data["mass"] >= massRangeMinimum) & (data["mass"] <= massRangeMaximum)
    if not np.any(selected):
        raise ValueError("no usable bins in the data mass function")
    mass              = data["mass"        ][selected]
    count             = data["count"       ][selected]
    massFunctionData  = data["massFunction"][selected]
    # Compute the conversion factor between halo count per bin and the mass
    # function, averaged over all selected, non-empty bins.
    nonEmpty = massFunctionData > 0.0
    countConversionFactor = np.mean(count[nonEmpty]/massFunctionData[nonEmpty])
    # Match data bins to model grid points. The grids are constructed to align
    # so matching is by index offset - but verify this.
    intervalLogarithmic = np.log(model["mass"][-1]/model["mass"][0])/(len(model["mass"])-1)
    indices = np.rint(np.log(mass/model["mass"][0])/intervalLogarithmic).astype(int)
    if np.any(indices < 0) or np.any(indices >= len(model["mass"])):
        raise ValueError("data mass bins extend beyond the model mass grid")
    misalignment = np.abs(np.log(mass/model["mass"][indices]))/intervalLogarithmic
    if np.any(misalignment > toleranceMisalignment):
        raise ValueError(
            f"data mass bins are misaligned with the model mass grid "
            f"(maximum misalignment = {np.max(misalignment):.2e} of the bin width)"
        )
    massFunctionModel = model["massFunction"][indices]
    # Expected count of halos per bin.
    countMean = countConversionFactor*massFunctionModel
    # Evaluate the log-likelihood.
    logLikelihood = 0.0
    for i in range(len(mass)):
        N  = int(count[i])
        mu = countMean[i]
        if mu <= 0.0:
            if N > 0:
                raise ValueError(
                    f"model predicts zero halos in bin at mass {mass[i]:.4e} "
                    f"where {N} are observed - likelihood is zero"
                )
            # No halos expected, none observed - no contribution.
            continue
        if varianceFractionalModelDiscrepancy > 0.0:
            # Negative binomial distribution (over-dispersed Poisson).
            r              = 1.0/varianceFractionalModelDiscrepancy
            logLikelihood += (
                +N*math.log(mu)
                -math.lgamma(N+1.0)
                +math.lgamma(r+N)
                -math.lgamma(r)
                -N*math.log(r+mu)
                -r*math.log(mu/r+1.0)
            )
        else:
            # Poisson distribution.
            logLikelihood += N*math.log(mu)-mu-math.lgamma(N+1.0)
    # Per-bin quantities for plotting. The error is the root-variance of the
    # count distribution (Poisson plus any model discrepancy term), converted
    # to mass function units.
    varianceCount = countMean+varianceFractionalModelDiscrepancy*countMean**2
    bins = {
        "mass"             : mass,
        "count"            : count,
        "countMean"        : countMean,
        "massFunctionData" : massFunctionData,
        "massFunctionModel": massFunctionModel,
        "massFunctionError": np.sqrt(varianceCount)/countConversionFactor,
    }
    return -logLikelihood, bins
