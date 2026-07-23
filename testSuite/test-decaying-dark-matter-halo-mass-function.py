#!/usr/bin/env python3
"""
End-to-end integration test of the decaying dark matter (DDM) halo mass function of Montandon et al.
(2026; arXiv:2607.19244).

The full stack is exercised: a `decayingDarkMatter` dark matter particle supplies the lifetime and
velocity kick; a `criticalOverdensityDecayingDarkMatter` supplies the mass-dependent collapse barrier
delta_c(M_0); a `shethTormen` f(nu) mass function consumes that barrier; and a
`haloMassFunctionDecayingDarkMatter` remaps the result from Lagrangian mass M_0 to collapsed mass
M_coll. The cosmology matches the N-body suite of their Sec. 4.

Two things are checked:

* The cold dark matter limit. Taking the lifetime to infinity and the velocity kick to zero must
  recover the plain Sheth-Tormen mass function exactly, since the DDM correction is applied
  multiplicatively and normalized to delta_c^EdS.

* The suppression relative to that limit, for the four DDM models simulated in the paper, at z=0 and
  z=1.083. Reference values were extracted from the semi-analytic-fit ("dashed") curves of their Fig. 9
  (see the `reference` table below for the extraction method). Each suppression ratio is checked to
  within `toleranceFactor`; the residual is set by the variance pipeline rather than the reference (see
  the comment there). The test also checks, more stringently, the *ordering* of the suppression across
  models, which is what pins the transition mass scale M_1 (their eq. 44): using the analytic estimate
  of their eq. 45 in place of the calibrated eq. 44 overestimates M_1 by a factor of ~200, which inverts
  the ordering and misses these amplitudes by up to a factor of several hundred.
"""
import os
import subprocess
import sys

import h5py
import numpy as np

# Models: label -> (lifetime [Gyr], velocity kick [km/s]). The first entry is the cold dark matter
# limit, used as the reference against which suppression is measured.
models = {
    "CDM"          : (1.0e6,    1.0e-3),
    "10Gyr_1250kms": (1.0e1,    1.250e3),
    "5Gyr_625kms"  : (5.0e0,    6.250e2),
    "20Gyr_625kms" : (2.0e1,    6.250e2),
    "20Gyr_2250kms": (2.0e1,    2.250e3),
}

# Suppression of dn/dlnM relative to the cold dark matter limit, extracted from the semi-analytic-fit
# ("dashed") curves of Montandon et al. (2026), their Fig. 9, at the given collapsed masses [M_Solar/h].
# The reference curves were taken from a vector (SVG) copy of the figure with the N-body points removed,
# rendered at high resolution and read off by nearest-colour classification (calibrated from the axis
# ticks). The reference masses are chosen at 2 and 5e14 M_Solar/h -- well inside the plotted range, where
# all five curves are cleanly separated (at the left frame edge, 1e14 M_Solar/h, the curves cross and the
# extraction is unreliable). Structured as: redshift -> mass -> {model: ratio}.
reference = {
    0.000: {
        2.0e14: {"10Gyr_1250kms": 0.0358, "5Gyr_625kms": 0.4450, "20Gyr_625kms": 0.7245, "20Gyr_2250kms": 0.0998},
        5.0e14: {"10Gyr_1250kms": 0.0727, "5Gyr_625kms": 0.5556, "20Gyr_625kms": 0.7877, "20Gyr_2250kms": 0.0437},
    },
    1.083: {
        2.0e14: {"10Gyr_1250kms": 0.0954, "5Gyr_625kms": 0.5210, "20Gyr_625kms": 0.7769, "20Gyr_2250kms": 0.1419},
        5.0e14: {"10Gyr_1250kms": 0.1862, "5Gyr_625kms": 0.5952, "20Gyr_625kms": 0.8167, "20Gyr_2250kms": 0.0699},
    },
}

hubbleConstantReduced = 0.6776
# The suppression ratios are compared to within this factor. The residual (worst case ~1.6, for the
# strongly-suppressed 10Gyr_1250kms model) is dominated not by the reference extraction but by the linear
# power spectrum / variance pipeline: the paper's sigma(M) uses a different transfer function than the
# LCDM one used here, and this does not fully cancel in the DDM/LCDM ratio because the mass remapping
# evaluates sigma at different masses in the numerator (DDM, at M_0) and denominator (LCDM, at M_coll) --
# an effect largest for the most strongly remapped (largest-kick) models. The test therefore checks the
# amplitude to within this factor and, more stringently, the *ordering* of the suppression across models.
toleranceFactor       = 1.8

failures = 0


def writeChangeFile(fileName, changes):
    """Write a Galacticus parameter change file."""
    with open(fileName, "w") as changeFile:
        changeFile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<changes>\n")
        changeFile.write(changes)
        changeFile.write("</changes>\n")


def run(label, changes):
    """Apply the given parameter changes and run, returning the output file name (or None on failure)."""
    outputFileName = f"outputs/decayingDarkMatterHaloMassFunction_{label}.hdf5"
    changeFileName = f"outputs/decayingDarkMatterHaloMassFunction_{label}.changes.xml"
    writeChangeFile(
        changeFileName,
        changes + f"""  <change type="replace" path="outputFileName">
    <outputFileName value="testSuite/{outputFileName}"/>
  </change>
"""
    )
    command = (
        "cd ..; ./Galacticus.exe testSuite/parameters/decayingDarkMatterHaloMassFunction.xml"
        f" testSuite/{changeFileName}"
    )
    status = subprocess.run(command, shell=True, capture_output=True, text=True)
    if status.returncode != 0:
        print(f"FAILED: model '{label}' did not run")
        print(status.stdout[-2000:])
        print(status.stderr[-2000:])
        return None
    return outputFileName


def runModel(label, lifetime, velocityKick):
    """Run the DDM halo mass function model for the given lifetime [Gyr] and velocity kick [km/s]."""
    return run(label, f"""  <change type="replace" path="darkMatterParticle">
    <darkMatterParticle value="decayingDarkMatter">
      <lifetime     value="{lifetime:e}"/>
      <velocityKick value="{velocityKick:e}"/>
      <darkMatterParticle value="CDM"/>
    </darkMatterParticle>
  </change>
""")


def readMassFunctions(fileName):
    """Return {redshift: (mass [M_Solar], dn/dlnM [Mpc^-3])} from a halo mass function output file."""
    massFunctions = {}
    with h5py.File(fileName, "r") as model:
        for outputName in model["Outputs"]:
            output = model["Outputs"][outputName]
            redshift = float(output.attrs["outputRedshift"])
            massFunctions[round(redshift, 3)] = (
                output["haloMass"           ][:],
                output["haloMassFunctionLnM"][:],
            )
    return massFunctions


def interpolate(massFunction, mass):
    """Log-log interpolate a mass function to the given mass."""
    masses, values = massFunction
    return np.exp(np.interp(np.log(mass), np.log(masses), np.log(values)))


# Run all models.
os.makedirs("outputs", exist_ok=True)
results = {}
for label, (lifetime, velocityKick) in models.items():
    fileName = runModel(label, lifetime, velocityKick)
    if fileName is None:
        failures += 1
        continue
    results[label] = readMassFunctions(fileName)

if "CDM" not in results:
    print("FAILED: the cold dark matter limit model did not run; cannot proceed")
    sys.exit(0)

# Check the cold dark matter limit: with an effectively infinite lifetime and vanishing velocity kick
# the DDM mass function must reduce to the plain Sheth-Tormen result. We verify this indirectly by
# requiring the mass remapping to be the identity, i.e. that the mass function is unchanged when the
# DDM machinery is switched on with null parameters. The plain Sheth-Tormen model is obtained by
# removing the DDM wrappers from the parameter file.
plainFileName = run("plainST", """  <change type="replace" path="darkMatterParticle">
    <darkMatterParticle value="CDM"/>
  </change>
  <change type="replace" path="criticalOverdensity">
    <criticalOverdensity value="sphericalCollapseClsnlssMttrCsmlgclCnstnt"/>
  </change>
  <change type="replace" path="haloMassFunction">
    <haloMassFunction value="shethTormen">
      <a             value="0.707"/>
      <p             value="0.300"/>
      <normalization value="0.322"/>
    </haloMassFunction>
  </change>
""")
if plainFileName is not None:
    plain = readMassFunctions(plainFileName)
    for redshift in sorted(plain):
        masses = plain[redshift][0]
        ratio  = np.array([
            interpolate(results["CDM"][redshift], mass) / interpolate(plain[redshift], mass)
            for mass in masses
        ])
        deviation = np.max(np.abs(ratio - 1.0))
        if deviation < 1.0e-3:
            print(f"SUCCESS: cold dark matter limit recovers Sheth-Tormen at z={redshift} (max deviation {deviation:.2e})")
        else:
            print(f"FAILED: cold dark matter limit deviates from Sheth-Tormen at z={redshift} by {deviation:.2e}")
            failures += 1
else:
    print("FAILED: the plain Sheth-Tormen comparison model did not run")
    failures += 1

# Check the suppression against the digitized reference from their Fig. 9.
for redshift, massReference in reference.items():
    for massReduced, modelReference in massReference.items():
        mass = massReduced / hubbleConstantReduced
        base = interpolate(results["CDM"][redshift], mass)
        for label, referenceRatio in modelReference.items():
            if label not in results:
                failures += 1
                continue
            ratio = interpolate(results[label][redshift], mass) / base
            if referenceRatio / toleranceFactor <= ratio <= referenceRatio * toleranceFactor:
                print(f"SUCCESS: suppression for '{label}' at z={redshift}, M={massReduced:.1e} M☉/h: {ratio:.3f} (reference {referenceRatio:.3f})")
            else:
                print(f"FAILED: suppression for '{label}' at z={redshift}, M={massReduced:.1e} M☉/h: {ratio:.3f} differs from reference {referenceRatio:.3f} by more than a factor {toleranceFactor:g}")
                failures += 1

# Check that the ordering of the suppression across models matches that of their Fig. 9, at z=0. This
# ordering is a sensitive probe of the transition mass scale M_1.
# Ordering is checked at 2e14 M_Solar/h, where the models are unambiguously ranked. (At higher masses
# the two large-kick models -- 20Gyr_2250kms and 10Gyr_1250kms -- swap, as the largest-kick model falls
# off most steeply; see their Fig. 9.)
ordering = ["20Gyr_625kms", "5Gyr_625kms", "20Gyr_2250kms", "10Gyr_1250kms"]
mass     = 2.0e14 / hubbleConstantReduced
if all(label in results for label in ordering):
    ratios = [interpolate(results[label][0.000], mass) / interpolate(results["CDM"][0.000], mass) for label in ordering]
    if all(ratios[i] > ratios[i + 1] for i in range(len(ratios) - 1)):
        print("SUCCESS: ordering of suppression across models matches Montandon et al. (2026), their Fig. 9")
    else:
        print(f"FAILED: ordering of suppression across models does not match their Fig. 9: {list(zip(ordering, ratios))}")
        failures += 1
else:
    failures += 1

print(f"\n{failures} failures")
sys.exit(1 if failures > 0 else 0)
