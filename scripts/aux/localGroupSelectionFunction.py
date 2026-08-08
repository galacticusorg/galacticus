#!/usr/bin/env python3
"""Build a sky-marginalized detection probability table for Milky Way satellites.

The DELVE Milky Way Census I (Tan, Drlica-Wagner et al. 2025; arXiv:2509.12313) provides an
observational selection function for Milky Way satellite galaxies in the form of a per-survey
gradient-boosted classifier which returns the probability that a satellite would be detected given

  * its heliocentric distance, :math:`D`;
  * its absolute V-band magnitude, :math:`M_\\mathrm{V}`;
  * its azimuthally-averaged, physical half-light radius, :math:`r_{1/2}`; and
  * the local density of foreground stars, :math:`\\rho_\\star`,

together with HEALPix maps of the survey footprints and of :math:`\\rho_\\star`.

Galacticus predicts satellite positions in the frame of the host halo, which carries no meaningful
orientation with respect to the Galactic plane or to the survey footprints. The appropriate
detection probability for such a model is therefore the average of the selection function over all
directions of the satellite from the Galactic center. Since the direction fixes both the sky
position (and so the footprint coverage and stellar density) and, given the Sun's offset from the
Galactic center, the heliocentric distance, this average is

.. math::

   \\langle P_\\mathrm{det}\\rangle(r,M_\\mathrm{V},r_{1/2}) = \\frac{1}{4\\pi} \\int
   P_\\mathrm{det}\\left(D(\\hat{n}),M_\\mathrm{V},r_{1/2},\\rho_\\star(\\hat{n})\\right)
   \\mathrm{d}\\Omega,

where :math:`r` is the galactocentric radius and :math:`D(\\hat{n})=|r\\hat{n}+\\mathbf{R}_\\odot|`.
This script evaluates that integral on a grid and writes the result to an HDF5 file for use by the
``localGroupDetection`` output analysis weight operator, so that neither the classifier nor any
HEALPix machinery is needed at run time.

The average is computed by sampling directions on a HEALPix grid and accumulating, for each
galactocentric radius, the distribution of (heliocentric distance, stellar density) pairs weighted
by the fraction of each pixel which is unmasked. The classifier is then evaluated once on a grid of
its four arguments and contracted with that distribution, rather than being called once per
direction.

Requires the data products released with the DELVE Milky Way Census I, at

   https://github.com/delve-survey/delve_mw_census

(``data/`` directory), and the ``xgboost``, ``healpy``, and ``astropy`` Python packages. These are
needed only to build the table -- not to run Galacticus.

Note that the Pan-STARRS classifier of that release is stored in the deprecated XGBoost binary
format, support for which has since been dropped. Use ``xgboost==3.0.1``, as pinned by the
``conda-env.yml`` of the release itself.

Andrew Benson <abenson@carnegiescience.edu>
"""

import argparse
import os
import sys

import h5py
import numpy as np
import yaml

# Bits of the DELVE census mask which indicate that a region is unusable. Taken from
# `data/read_config.yaml` of the DELVE census data release: regions near bright stars, near
# catalogued objects, of high Galactic extinction, of high foreground stellar density, or with no
# survey coverage. Note that the bits identifying the survey footprints are *not* included here --
# each survey's mask already sets the "no coverage" bit outside of that survey's own footprint.
maskBitsBad = ("ASSOC", "STAR", "EBV", "DENSITY", "FOOT")

# Surveys contributing to the census. The footprints are disjoint after masking, so the total
# detection probability is the sum over surveys.
surveys = ("desy6", "delvedr3", "ps1d")

# Limits imposed by the census, following the reference implementation in Notebook 1 of the DELVE
# census data release: satellites brighter than `magnitudeBright` are taken to be detected
# irrespective of the survey coverage (they are the classical satellites, all of which were known
# long before these surveys), while satellites outside the heliocentric distance range over which
# the classifiers were trained are taken to be undetectable.
magnitudeBright     = -12.5
distanceMinimum     =  16.0
distanceMaximum     = 400.0

# Solar galactocentric radius (kpc; GRAVITY Collaboration 2018; A&A, 615, L15).
radiusSolar         =   8.122


def coverageAndDensity(dataDirectory, config, survey, nSideOutput):
    """Return the unmasked fraction and foreground stellar density of each pixel for a survey.

    The mask is provided at much higher resolution than the stellar density map. It is degraded to
    the resolution of the density map by computing, for each low-resolution pixel, the fraction of
    its sub-pixels which are unmasked. Since the selection function depends on sky position only
    through the mask and the stellar density, and the density is available only at the lower
    resolution, this loses no information.
    """
    import healpy as hp

    bits = config["maskbits"]
    bitsBad = 0
    for name in maskBitsBad:
        bitsBad |= bits[name]
    mask = hp.read_map(os.path.join(dataDirectory, config[survey]["mask"]), dtype=np.int64)
    coverage = hp.ud_grade(
        np.asarray((mask & bitsBad) == 0, dtype=np.float32), nSideOutput, order_in="RING", order_out="RING"
    )
    del mask
    density = hp.read_map(os.path.join(dataDirectory, config[survey]["density"]))
    if hp.get_nside(density) != nSideOutput:
        raise ValueError(f"stellar density map for '{survey}' has unexpected resolution")
    # Regions with no valid density measurement can not be evaluated, so are treated as having no
    # coverage.
    densityValid = np.isfinite(density) & (density > hp.UNSEEN / 2.0) & (density >= 0.0)
    coverage = np.where(densityValid, coverage, 0.0)
    density = np.where(densityValid, density, 0.0)
    return coverage, density


def interpolationWeights(values, grid):
    """Return indices and weights for linear interpolation of `values` onto a monotonic `grid`.

    Values are clamped to the range of the grid. Returned as a pair of (index, weight) pairs, such
    that a quantity tabulated on the grid is interpolated as the weighted sum over the two indices.
    """
    indexUpper = np.clip(np.searchsorted(grid, values, side="left"), 1, grid.size - 1)
    indexLower = indexUpper - 1
    weightUpper = (values - grid[indexLower]) / (grid[indexUpper] - grid[indexLower])
    weightUpper = np.clip(weightUpper, 0.0, 1.0)
    return (indexLower, 1.0 - weightUpper), (indexUpper, weightUpper)


def distributionDistanceDensity(radius, directions, rotation, coverage, density, nSide, gridDistance, gridDensity):
    """Accumulate the distribution of (heliocentric distance, stellar density) over sky directions.

    For a satellite at galactocentric radius `radius`, with direction from the Galactic center
    sampled uniformly by `directions`, this returns the mean coverage-weighted number of directions
    falling in each cell of the (distance, density) grid, plus the fraction of directions for which
    the heliocentric distance lies within the range over which the census is defined.
    """
    import healpy as hp

    # Position of the satellite relative to the Sun, in Galactic Cartesian coordinates with the
    # Galactic center at (R☉,0,0).
    position = radius * directions
    position[0, :] += radiusSolar
    distance = np.sqrt(np.sum(position**2, axis=0))
    inRange = (distance >= distanceMinimum) & (distance <= distanceMaximum)
    # Convert the heliocentric direction to equatorial coordinates -- the frame in which the survey
    # maps are defined -- and find the corresponding pixel.
    directionsHeliocentric = rotation @ (position / distance)
    pixel = hp.vec2pix(nSide, directionsHeliocentric[0, :], directionsHeliocentric[1, :], directionsHeliocentric[2, :])
    weight = np.where(inRange, coverage[pixel], 0.0)
    (indexDistanceLower, weightDistanceLower), (indexDistanceUpper, weightDistanceUpper) = interpolationWeights(
        np.log10(np.maximum(distance, distanceMinimum)), np.log10(gridDistance)
    )
    (indexDensityLower, weightDensityLower), (indexDensityUpper, weightDensityUpper) = interpolationWeights(
        density[pixel], gridDensity
    )
    distribution = np.zeros((gridDistance.size, gridDensity.size))
    for indexDistance, weightDistance in (
        (indexDistanceLower, weightDistanceLower),
        (indexDistanceUpper, weightDistanceUpper),
    ):
        for indexDensity, weightDensity in (
            (indexDensityLower, weightDensityLower),
            (indexDensityUpper, weightDensityUpper),
        ):
            np.add.at(
                distribution,
                (indexDistance, indexDensity),
                weight * weightDistance * weightDensity,
            )
    return distribution / directions.shape[1], np.mean(inRange)


def loadModel(xgboost, fileName):
    """Load a classifier, reporting the version requirement if the file format is unsupported."""
    model = xgboost.Booster()
    try:
        model.load_model(fileName)
    except xgboost.core.XGBoostError as error:
        raise RuntimeError(
            f"failed to load classifier '{fileName}' with xgboost {xgboost.__version__} -- the"
            " Pan-STARRS classifier of the DELVE census release is stored in the deprecated XGBoost"
            " binary format, so `xgboost==3.0.1` (as pinned by that release) is required"
        ) from error
    return model


def probabilityPredict(model, arguments):
    """Return the detection probability for each row of `arguments`.

    The classifier arguments are log₁₀(heliocentric distance/kpc), absolute V-band magnitude,
    log₁₀(half-light radius/kpc), and foreground stellar density. The low-level `Booster` interface
    is used in preference to the `XGBClassifier` interface used by the reference implementation,
    solely to avoid a dependency on `scikit-learn`; the two give identical results.
    """
    import xgboost

    return model.predict(xgboost.DMatrix(arguments))


def probabilityGrid(model, gridDistance, gridMagnitude, gridRadiusHalf, gridDensity):
    """Evaluate a classifier on the outer product of its four arguments."""
    logDistance = np.log10(gridDistance)
    logRadiusHalf = np.log10(gridRadiusHalf / 1.0e3)
    shape = (gridDistance.size, gridMagnitude.size, gridRadiusHalf.size, gridDensity.size)
    arguments = np.empty((np.prod(shape), 4))
    grids = np.meshgrid(logDistance, gridMagnitude, logRadiusHalf, gridDensity, indexing="ij")
    for i, grid in enumerate(grids):
        arguments[:, i] = grid.reshape(-1)
    return probabilityPredict(model, arguments).reshape(shape)


def build(dataDirectory, outputFileName, nSideDirections, validate):
    import astropy.coordinates
    import astropy.units
    import healpy as hp
    import xgboost

    config = yaml.safe_load(open(os.path.join(dataDirectory, "read_config.yaml")))

    # Grids on which the detection probability is tabulated. Deliberately given different lengths so
    # that any transposition of the resulting table is immediately detectable.
    gridRadius     = np.logspace(np.log10(  5.0), np.log10(450.0), 61)  # Galactocentric radius (kpc).
    gridMagnitude  = np.linspace(          -20.0,             2.0, 45)  # Absolute V-band magnitude.
    gridRadiusHalf = np.logspace(np.log10(  1.0), np.log10(5.0e3), 41)  # Projected half-light radius (pc).
    # Grid used internally to marginalize over sky position.
    gridDistance   = np.logspace(np.log10(distanceMinimum), np.log10(distanceMaximum), 100)

    # Rotation from Galactic to equatorial Cartesian coordinates.
    frameGalactic = astropy.coordinates.Galactic(
        u=[1.0, 0.0, 0.0] * astropy.units.kpc,
        v=[0.0, 1.0, 0.0] * astropy.units.kpc,
        w=[0.0, 0.0, 1.0] * astropy.units.kpc,
        representation_type="cartesian",
    )
    frameEquatorial = frameGalactic.transform_to(astropy.coordinates.ICRS())
    frameEquatorial.representation_type = "cartesian"
    rotation = np.vstack(
        [
            frameEquatorial.x.value,
            frameEquatorial.y.value,
            frameEquatorial.z.value,
        ]
    )

    # Directions of the satellite from the Galactic center, sampled uniformly over the sphere. A
    # HEALPix grid is used since its pixels are of equal area, so all directions carry equal weight.
    nPixel = hp.nside2npix(nSideDirections)
    directions = np.vstack(hp.pix2vec(nSideDirections, np.arange(nPixel)))

    print(f"Marginalizing over {nPixel} directions for {gridRadius.size} galactocentric radii...")
    probability = np.zeros((gridRadius.size, gridMagnitude.size, gridRadiusHalf.size))
    fractionInRange = np.zeros(gridRadius.size)
    nSideMaps = None
    for survey in surveys:
        print(f"  survey '{survey}':")
        print("    loading footprint mask and stellar density map...")
        density = hp.read_map(os.path.join(dataDirectory, config[survey]["density"]))
        nSideMaps = hp.get_nside(density)
        coverage, density = coverageAndDensity(dataDirectory, config, survey, nSideMaps)
        # The census masks regions of high foreground stellar density, so only the range of densities
        # which actually survives the mask need be tabulated.
        gridDensity = np.linspace(0.0, 1.001 * np.max(density[coverage > 0.0]), 41)
        print("    evaluating classifier...")
        model = loadModel(xgboost, os.path.join(dataDirectory, config[survey]["model"]))
        probabilitySurvey = probabilityGrid(model, gridDistance, gridMagnitude, gridRadiusHalf, gridDensity)
        print("    marginalizing over sky position...")
        for i, radius in enumerate(gridRadius):
            distribution, inRange = distributionDistanceDensity(
                radius, directions, rotation, coverage, density, nSideMaps, gridDistance, gridDensity
            )
            probability[i, :, :] += np.einsum("dq,dmrq->mr", distribution, probabilitySurvey)
            fractionInRange[i] = inRange
    # Satellites brighter than the census threshold are taken to be detected in any direction for
    # which the heliocentric distance is within range, irrespective of survey coverage.
    isBright = gridMagnitude < magnitudeBright
    probability[:, isBright, :] = fractionInRange[:, np.newaxis, np.newaxis]
    probability = np.clip(probability, 0.0, 1.0)

    if validate:
        validateTable(
            dataDirectory,
            config,
            probability,
            gridRadius,
            gridMagnitude,
            gridRadiusHalf,
            directions,
            rotation,
            nSideMaps,
        )

    print(f"Writing '{outputFileName}'...")
    with h5py.File(outputFileName, "w") as file:
        file.attrs["description"] = (
            "Sky-marginalized probability that a Milky Way satellite galaxy is detected in the "
            "DELVE Milky Way Census I, as a function of galactocentric radius, absolute V-band "
            "magnitude, and projected half-light radius. Computed by averaging the observational "
            "selection function of Tan, Drlica-Wagner et al. (2025) over all directions of the "
            "satellite from the Galactic center, assuming an isotropic satellite distribution."
        )
        file.attrs["source"                 ] = "Tan, Drlica-Wagner et al. (2025); ApJ; 1000; 87"
        file.attrs["sourceURL"              ] = "https://arxiv.org/abs/2509.12313"
        file.attrs["sourceDOI"              ] = "10.5281/zenodo.18383157"
        file.attrs["sourceRepository"       ] = "https://github.com/delve-survey/delve_mw_census"
        file.attrs["surveys"                ] = ", ".join(surveys)
        file.attrs["models"                 ] = ", ".join(config[survey]["model"] for survey in surveys)
        file.attrs["masks"                  ] = ", ".join(config[survey]["mask" ] for survey in surveys)
        file.attrs["radiusSolar"            ] = radiusSolar
        file.attrs["magnitudeBright"        ] = magnitudeBright
        file.attrs["distanceMinimum"        ] = distanceMinimum
        file.attrs["distanceMaximum"        ] = distanceMaximum
        file.attrs["nSideDirections"        ] = nSideDirections
        file.attrs["builtBy"                ] = "scripts/aux/localGroupSelectionFunction.py"
        dataset = file.create_dataset("radiusGalactocentric", data=gridRadius / 1.0e3)
        dataset.attrs["units"] = "Mpc"
        dataset.attrs["description"] = "Galactocentric radius of the satellite."
        dataset = file.create_dataset("magnitudeAbsoluteV", data=gridMagnitude)
        dataset.attrs["units"] = "AB magnitudes"
        dataset.attrs["description"] = "Absolute V-band magnitude of the satellite."
        dataset = file.create_dataset("radiusHalfLight", data=gridRadiusHalf / 1.0e6)
        dataset.attrs["units"] = "Mpc"
        dataset.attrs["description"] = "Projected, azimuthally-averaged half-light radius of the satellite."
        # Transposed on output so that the array is read by Galacticus with the shape
        # (radiusGalactocentric, magnitudeAbsoluteV, radiusHalfLight).
        dataset = file.create_dataset("detectionProbability", data=probability.T)
        dataset.attrs["description"] = (
            "Probability that the satellite is detected in the census, averaged over sky position."
        )
    print("...done")


def validateTable(
    dataDirectory, config, probability, gridRadius, gridMagnitude, gridRadiusHalf, directions, rotation, nSideMaps
):
    """Check the binned marginalization against a direct evaluation at a few points.

    The table is built by binning the distribution of heliocentric distance and stellar density over
    sky directions and contracting it with the classifier evaluated on a grid. This recomputes the
    average by calling the classifier for every direction, at a handful of grid points, and reports
    the largest discrepancy.
    """
    import healpy as hp
    import xgboost

    print("Validating binned marginalization against direct evaluation...")
    models = {}
    maps = {}
    for survey in surveys:
        models[survey] = loadModel(xgboost, os.path.join(dataDirectory, config[survey]["model"]))
        maps[survey] = coverageAndDensity(dataDirectory, config, survey, nSideMaps)
    # Subsample the directions -- the direct evaluation is expensive. Sampling is done by taking the
    # first sub-pixel of each pixel of a coarser HEALPix grid, which gives one direction per equal
    # area cell and so remains an unbiased sample of the sphere.
    nSideDirections = hp.npix2nside(directions.shape[1])
    subsample = hp.nest2ring(nSideDirections, np.arange(0, directions.shape[1], 16))
    directionsSubsampled = directions[:, subsample]
    deviationMaximum = 0.0
    for indexRadius in (10, 25, 40):
        for indexMagnitude in (25, 32, 38):
            for indexRadiusHalf in (10, 20, 30):
                radius = gridRadius[indexRadius]
                magnitude = gridMagnitude[indexMagnitude]
                radiusHalf = gridRadiusHalf[indexRadiusHalf]
                position = radius * directionsSubsampled
                position[0, :] += radiusSolar
                distance = np.sqrt(np.sum(position**2, axis=0))
                inRange = (distance >= distanceMinimum) & (distance <= distanceMaximum)
                directionsHeliocentric = rotation @ (position / distance)
                pixel = hp.vec2pix(nSideMaps, *directionsHeliocentric)
                total = np.zeros(distance.size)
                for survey in surveys:
                    coverage, density = maps[survey]
                    arguments = np.vstack(
                        [
                            np.log10(np.maximum(distance, distanceMinimum)),
                            np.full(distance.size, magnitude),
                            np.full(distance.size, np.log10(radiusHalf / 1.0e3)),
                            density[pixel],
                        ]
                    ).T
                    total += coverage[pixel] * probabilityPredict(models[survey], arguments)
                total = np.where(inRange, total, 0.0)
                if magnitude < magnitudeBright:
                    total = np.where(inRange, 1.0, 0.0)
                direct = np.clip(np.mean(total), 0.0, 1.0)
                binned = probability[indexRadius, indexMagnitude, indexRadiusHalf]
                deviationMaximum = max(deviationMaximum, abs(direct - binned))
                print(
                    f"  r={radius:7.2f} kpc  M_V={magnitude:6.2f}  r½={radiusHalf:8.2f} pc"
                    f"  direct={direct:.5f}  binned={binned:.5f}"
                )
    print(f"  maximum deviation: {deviationMaximum:.5f}")
    if deviationMaximum > 0.01:
        print("  WARNING: binned marginalization deviates from direct evaluation by more than 0.01", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--dataDirectory",
        required=True,
        help="directory containing the `data/` products of the DELVE Milky Way Census I release",
    )
    parser.add_argument("--outputFileName", required=True, help="name of the HDF5 file to write")
    parser.add_argument(
        "--nSideDirections",
        type=int,
        default=256,
        help="HEALPix resolution used to sample directions from the Galactic center (default: 256)",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="check the binned marginalization against a direct evaluation of the selection function",
    )
    arguments = parser.parse_args()
    build(arguments.dataDirectory, arguments.outputFileName, arguments.nSideDirections, arguments.validate)
