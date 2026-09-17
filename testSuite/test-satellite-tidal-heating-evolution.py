#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Accumulate Gnedin (1999) tidal heating along one satellite's orbit in a static host potential, and compare
# it against values computed independently by `satelliteTidalHeatingEvolution.py` in the galacticusDevTools
# repository.
#
# This is the companion to `tests.satellite_tidal_heating_rate.exe`. That test verifies the heating rate
# pointwise, with the path-integrated tidal tensor set directly; this one verifies that the rate, and the path
# integral it reads, are accumulated correctly along an orbit. A pointwise comparison cannot catch a rate which
# is correct but enters the differential equations with the wrong sign, factor or decay term, and a comparison
# of the accumulated heating alone cannot tell such a mistake from a wrong rate - which is why both exist.
#
# The reference integrates, alongside the orbit and mass loss of `test-satellite-orbit-evolution.py`,
#
#     dG_ij/dt = g_ij - efficiencyDecay G_ij / T_orb,   dQ/dt = (epsilon/3) A(omega T_shock) g_ij G_ij
#
# where g_ij is the host's tidal tensor at the satellite's position and A is the adiabatic correction. The
# satellite's profile is not heated in response - `darkMatterProfile` is `darkMatterOnly` - so Q accumulates
# but nothing consumes it, which keeps this a test of the accumulation rather than of the profile's response;
# the latter is covered analytically by `tests.dark_matter_profiles.heated.exe`.
#
# The orbit and bound mass are asserted here too. They are inputs to the heating, and a drift in either would
# otherwise surface only as an unexplained disagreement in Q.
#
# Andrew Benson, Claude

# The reference trajectory and heating, from `satelliteTidalHeatingEvolution.py --fortran`.
timesReference             = [9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.380000000000000e+01]
radiusReference            = [2.361165239295172e-01, 1.320251014772402e-01, 1.053610869562232e-01, 2.105308038570400e-01, 2.906800920723504e-01, 3.330988989113521e-01]
massBoundReference         = [8.902557651947208e+09, 6.401239927651635e+09, 3.082872403843493e+09, 3.082872403843493e+09, 3.082872403843493e+09, 3.082872403843493e+09]
heatingNormalizedReference = [6.372821594020219e+01, 6.893741740573513e+03, 6.564907938393738e+04, 6.746186256425572e+04, 6.746680471493296e+04, 6.746691671921306e+04]

# Tolerances, set from the measured agreement with a factor of a few in hand rather than tuned to it.
#
# The orbit is reproduced far more sharply than in `test-satellite-orbit-evolution.py` as first written, because
# the reference now starts from the initial conditions the tree states rather than recomputing them: measured
# 3.1e-6 in radius.
#
# The bound mass is limited by neither integrator but by the King (1962) tidal radius, found by
# `radiusEnclosingDensityNumerical` with a relative tolerance of 1e-3; nothing downstream of it can agree
# better. Measured 1.6e-4.
#
# The heating inherits that, amplified: the adiabatic correction depends on the satellite's orbital frequency,
# which is set by its bound mass, and dln A/dln omega approaches 2 gamma = 5 at large omega T_shock. Measured
# 4.9e-4.
toleranceRadius            = 1.0e-4
toleranceMassBound         = 1.0e-3
toleranceHeating           = 2.0e-3

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/satelliteTidalHeatingEvolution.xml", shell=True)
if status.returncode != 0:
    print("FAILED: satellite tidal heating evolution model run")
    sys.exit(0)
print("SUCCESS: satellite tidal heating evolution model run")

# Extract the trajectory and the accumulated heating. The satellite is the single non-isolated node at each output.
try:
    model = h5py.File("outputs/satelliteTidalHeatingEvolution.hdf5", "r")
except OSError as error:
    print(f"FAILED: unable to open the model output: {error}")
    sys.exit(0)

outputs = model["Outputs"]
names   = sorted(outputs.keys(), key=lambda name: int(name[6:]))
if len(names) != len(timesReference):
    print(f"FAILED: expected {len(timesReference)} outputs but found {len(names)}")
    sys.exit(0)

radius    = []
massBound = []
heating   = []
for name in names:
    nodeData   = outputs[name]["nodeData"]
    isIsolated = nodeData["nodeIsIsolated"][:]
    satellites = np.nonzero(isIsolated == 0)[0]
    if len(satellites) != 1:
        print(f"FAILED: expected exactly one satellite in {name} but found {len(satellites)}")
        sys.exit(0)
    i        = satellites[0]
    position = np.array([nodeData[f"satellitePosition{axis}"][i] for axis in ("X", "Y", "Z")])
    radius   .append(np.linalg.norm(position)                          )
    massBound.append(nodeData["satelliteBoundMass"            ][i]     )
    heating  .append(nodeData["satelliteTidalHeatingNormalized"][i]    )
radius    = np.array(radius   )
massBound = np.array(massBound)
heating   = np.array(heating  )

# Compare, reporting the fractional difference of every comparison so that a reader of the log can see how much
# of the tolerance is in use. A test which passes only because its tolerance is loose otherwise looks identical
# to one which passes sharply.
differenceRadius    = np.abs(radius    / np.array(radiusReference           ) - 1.0)
differenceMassBound = np.abs(massBound / np.array(massBoundReference        ) - 1.0)
differenceHeating   = np.abs(heating   / np.array(heatingNormalizedReference) - 1.0)
for time, dRadius, dMass, dHeat in zip(timesReference, differenceRadius, differenceMassBound, differenceHeating):
    print(f"         t = {time:5.2f} Gyr: fractional difference, radius {dRadius:9.3e}, bound mass {dMass:9.3e}, Q {dHeat:9.3e}")
print(
    f"         largest fractional difference: radius {differenceRadius.max():9.3e},"
    f" bound mass {differenceMassBound.max():9.3e}, Q {differenceHeating.max():9.3e}"
)

for label, difference, tolerance, values, reference in (
    ("orbital radius"          , differenceRadius   , toleranceRadius   , radius   , radiusReference           ),
    ("bound mass"              , differenceMassBound, toleranceMassBound, massBound, massBoundReference        ),
    ("accumulated tidal heating", differenceHeating , toleranceHeating  , heating  , heatingNormalizedReference),
):
    if np.all(difference < tolerance):
        print(f"SUCCESS: satellite {label} matches the independent integration")
    else:
        print(f"FAILED: satellite {label} does not match the independent integration")
        for time, got, expected in zip(timesReference, values, reference):
            print(f"\tt = {time:5.2f}: {got} vs {expected}")

# Guards against a degenerate pass. Were the satellite frozen at its initial conditions, or the heating never
# applied, the comparisons above would still be made - and against a reference generated from those same
# initial conditions some of them could even pass. Require that the orbit is resolved, that mass is lost, and
# that the heating actually accumulates over a wide dynamic range rather than staying near its initial value of
# zero.
if radius.max() / radius.min() > 2.0:
    print("SUCCESS: the orbit is resolved, with apocenter more than twice pericenter")
else:
    print(f"FAILED: the orbit is not resolved: apocenter/pericenter = {radius.max() / radius.min()}")

if massBound[-1] < 0.5 * massBound[0]:
    print("SUCCESS: the satellite loses mass over the orbit")
else:
    print(f"FAILED: the satellite loses little mass: final/initial bound mass = {massBound[-1] / massBound[0]}")

if heating[0] > 0.0 and heating[-1] > 100.0 * heating[0]:
    print("SUCCESS: tidal heating accumulates over the orbit, by more than two orders of magnitude")
else:
    print(f"FAILED: tidal heating does not accumulate: Q rises only from {heating[0]} to {heating[-1]}")

sys.exit(0)
