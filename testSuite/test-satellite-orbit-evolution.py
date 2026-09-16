#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Integrate one satellite's orbit and tidal mass loss in a static host potential, and compare the trajectory
# against values computed independently by `satelliteOrbitEvolution.py` in the galacticusDevTools repository.
#
# This is the companion to `tests.satellite_orbit_rates.exe`. That test verifies the dynamical friction
# acceleration and the tidal mass loss rate pointwise; this one verifies that those rates are assembled
# correctly into the orbital differential equations. A pointwise comparison cannot catch a rate which is
# correct but enters the equations of motion with the wrong sign, factor or units, and a comparison of an
# integrated trajectory alone cannot distinguish such a mistake from a wrong rate - which is why both exist.
#
# The reference integrates
#
#     dr/dt = v,   dv/dt = a_host(r) (1 + massRatio) + a_friction,   dM/dt = massLossRate
#
# with an explicit Runge-Kutta method, from the same initial conditions the tree specifies. Galacticus reaches
# the same trajectory through its own evolver, so agreement is a statement about the equations rather than
# about either integrator.
#
# Andrew Benson, Claude

# The reference trajectory, from `satelliteOrbitEvolution.py --fortran`.
timesReference     = [9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.380000000000000e+01]
radiusReference    = [2.361335130769703e-01, 1.320493456127232e-01, 1.053422034374278e-01, 2.105064760854269e-01, 2.906644224382638e-01, 3.330919965940063e-01]
massBoundReference = [8.902780672104675e+09, 6.402069970543102e+09, 3.083010452526562e+09, 3.083010452526562e+09, 3.083010452526562e+09, 3.083010452526562e+09]

# Tolerances.
#
# Both are set from the measured agreement with a factor of a few in hand, not tuned to it. The floor is set by
# the two integrators, not by the physics: Galacticus evolves with a relative ODE tolerance of 1e-8 accumulated
# over 5.8 Gyr and through a pericentre passage, while the reference is converged to 2e-8 in radius and 2e-7 in
# mass (its own `--converge` mode compares tolerances of 1e-11 and 1e-8). Measured differences are 1.9e-4 in
# radius and 1.2e-4 in bound mass.
toleranceRadius    = 1.0e-3
toleranceMassBound = 1.0e-3

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/satelliteOrbitEvolution.xml", shell=True)
if status.returncode != 0:
    print("FAILED: satellite orbit evolution model run")
    sys.exit(0)
print("SUCCESS: satellite orbit evolution model run")

# Extract the satellite's trajectory. The satellite is the single non-isolated node at each output.
try:
    model = h5py.File("outputs/satelliteOrbitEvolution.hdf5", "r")
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
for name in names:
    nodeData   = outputs[name]["nodeData"]
    isIsolated = nodeData["nodeIsIsolated"][:]
    satellites = np.nonzero(isIsolated == 0)[0]
    if len(satellites) != 1:
        print(f"FAILED: expected exactly one satellite in {name} but found {len(satellites)}")
        sys.exit(0)
    i        = satellites[0]
    position = np.array([nodeData[f"satellitePosition{axis}"][i] for axis in ("X", "Y", "Z")])
    radius   .append(np.linalg.norm(position)             )
    massBound.append(nodeData["satelliteBoundMass"][i]    )
radius    = np.array(radius   )
massBound = np.array(massBound)

# Compare, reporting the fractional difference of every comparison so that a reader of the log can see how much
# of the tolerance is in use. A test which passes only because its tolerance is loose otherwise looks identical
# to one which passes sharply.
differenceRadius    = np.abs(radius    / np.array(radiusReference   ) - 1.0)
differenceMassBound = np.abs(massBound / np.array(massBoundReference) - 1.0)
for time, dRadius, dMass in zip(timesReference, differenceRadius, differenceMassBound):
    print(f"         t = {time:5.2f} Gyr: fractional difference, radius {dRadius:9.3e}, bound mass {dMass:9.3e}")
print(f"         largest fractional difference: radius {differenceRadius.max():9.3e}, bound mass {differenceMassBound.max():9.3e}")

if np.all(differenceRadius < toleranceRadius):
    print("SUCCESS: satellite orbital radius matches the independent integration")
else:
    print("FAILED: satellite orbital radius does not match the independent integration")
    for time, got, expected in zip(timesReference, radius, radiusReference):
        print(f"\tt = {time:5.2f}: {got} vs {expected}")

if np.all(differenceMassBound < toleranceMassBound):
    print("SUCCESS: satellite bound mass matches the independent integration")
else:
    print("FAILED: satellite bound mass does not match the independent integration")
    for time, got, expected in zip(timesReference, massBound, massBoundReference):
        print(f"\tt = {time:5.2f}: {got} vs {expected}")

# The trajectory must actually be a trajectory. Were the satellite frozen at its initial conditions - which is
# what happens if the node's time never advances, as when the `cosmicTime` node operator is absent - every
# comparison above would still be made, and against a reference which had been generated from those same
# initial conditions it could even pass. Require that the orbit is resolved and that mass is lost.
if radius.max() / radius.min() > 2.0:
    print("SUCCESS: the orbit is resolved, with apocenter more than twice pericenter")
else:
    print(f"FAILED: the orbit is not resolved: apocenter/pericenter = {radius.max() / radius.min()}")

if massBound[-1] < 0.5 * massBound[0]:
    print("SUCCESS: the satellite loses mass over the orbit")
else:
    print(f"FAILED: the satellite loses little mass: final/initial bound mass = {massBound[-1] / massBound[0]}")

sys.exit(0)
