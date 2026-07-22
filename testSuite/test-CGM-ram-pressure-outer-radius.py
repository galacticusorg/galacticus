#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Regression test: CGM ram pressure stripping must be applied to a satellite whose CGM outer radius
# exceeds its virial radius.
#
# Prior to the CGM refactors the hot halo outer radius was clamped to the virial radius, so the
# guard `radiusOuter <= radiusVirial` in the CGM outer radius ram pressure stripping node operator
# was identically satisfied and stripping always fired for satellites. Once the clamp became
# optional the outer radius tracked the virial radius only to floating point accuracy, and any
# satellite whose outer radius landed fractionally above the virial radius had ram pressure
# stripping silently disabled -- it retained its gas and formed far more stars.
#
# The model sets a satellite's CGM outer radius to 1.5x its virial radius explicitly, so the regime
# is entered deliberately and by a wide margin rather than by floating point accident. With no
# cooling, star formation or feedback, ram pressure stripping is the only process that can change
# the satellite's hot halo mass: without the fix the mass is exactly constant, with it the mass
# decreases.
#
# Andrew Benson

status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/cgmRamPressureOuterRadius.xml", shell=True)
if status.returncode != 0:
    print("FAILED: model run")
    sys.exit(1)
print("SUCCESS: model run")

with h5py.File("outputs/cgmRamPressureOuterRadius.hdf5", "r") as model:
    outputs = model["Outputs"]

    def read(output, name):
        return outputs[f"Output{output}/nodeData/{name}"][:]

    isSatellite  = read(1, "nodeIsIsolated") == 0
    if isSatellite.sum() != 1:
        print(f"FAILED: expected exactly 1 satellite in the first output, found {isSatellite.sum()}")
        sys.exit(1)

    radiusOuterInitial  = read(1, "hotHaloOuterRadius"        )[isSatellite][0]
    radiusVirialInitial = read(1, "darkMatterOnlyRadiusVirial")[isSatellite][0]
    massInitial         = read(1, "hotHaloMass"               )[isSatellite][0]

    isSatelliteFinal = read(2, "nodeIsIsolated") == 0
    if isSatelliteFinal.sum() != 1:
        print(f"FAILED: expected exactly 1 satellite in the final output, found {isSatelliteFinal.sum()}")
        sys.exit(1)
    massFinal = read(2, "hotHaloMass")[isSatelliteFinal][0]

# The virial radius reported above comes from the global darkMatterHaloScale, which is also the one
# the hot halo component's (deferred, clamping) outerRadius getter uses. The ram pressure stripping
# node operator is deliberately given its own darkMatterHaloScale with 8x the virial density
# contrast, and the guard under test used THAT virial radius. Since the virial radius scales as
# densityContrast^(-1/3), the operator sees exactly half the global value.
#
# Keep this factor in step with <densityContrastValue> under the CGMOuterRadiusRamPressureStripping
# node operator in testSuite/parameters/cgmRamPressureOuterRadius.xml.
radiusVirialOperator = radiusVirialInitial * (329.6207717164546 / 2636.9661737316368)**(1.0/3.0)

# Verify the test's own preconditions. Without these the test could silently become vacuous if the
# model setup ever drifted -- it would then pass on both fixed and unfixed code while testing
# nothing.
if abs(radiusOuterInitial - radiusVirialInitial) > 1.0e-6 * radiusVirialInitial:
    print(f"FAILED: test precondition not met -- satellite CGM outer radius ({radiusOuterInitial:.6f} "
          f"Mpc) is not clamped to the virial radius ({radiusVirialInitial:.6f} Mpc) as the hot halo "
          f"component's outerRadius getter should enforce; the model setup has drifted")
    sys.exit(1)

if radiusOuterInitial > radiusVirialOperator:
    print(f"SUCCESS: satellite CGM outer radius exceeds the virial radius seen by the ram pressure "
          f"stripping operator ({radiusOuterInitial:.6f} > {radiusVirialOperator:.6f} Mpc)")
else:
    print(f"FAILED: test precondition not met -- satellite CGM outer radius does not exceed the "
          f"virial radius seen by the ram pressure stripping operator "
          f"({radiusOuterInitial:.6f} <= {radiusVirialOperator:.6f} Mpc); "
          f"the model no longer exercises the regime this test targets")
    sys.exit(1)

# The actual regression check.
if massFinal < massInitial:
    print(f"SUCCESS: satellite CGM was ram pressure stripped "
          f"({massInitial:.4e} -> {massFinal:.4e} M_sun, "
          f"{100.0*(1.0-massFinal/massInitial):.2f}% removed)")
else:
    print(f"FAILED: satellite CGM was not ram pressure stripped "
          f"({massInitial:.4e} -> {massFinal:.4e} M_sun); ram pressure stripping is being skipped "
          f"for a satellite whose CGM outer radius exceeds its virial radius")
    sys.exit(1)
