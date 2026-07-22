#!/usr/bin/env python3
import subprocess
import sys
import h5py

# Regression test: CGM ram pressure stripping must be applied to a satellite whose CGM outer radius
# exceeds its virial radius.
#
# Ram pressure stripping of a satellite's CGM used to be gated on
# `radiusOuter <= darkMatterHaloScale_%radiusVirial(node)`. That comparison was never between like
# quantities: the outer radius returned by the standard hot halo component's deferred getter is
# already clamped to [0.1*rVir, rVir], but against a darkMatterHaloScale that the COMPONENT builds
# from its own subParameters, whereas the guard compared it against a SEPARATE darkMatterHaloScale
# built by the node operator from the top-level parameters. Being distinct instances the two are not
# guaranteed to agree bit-for-bit, and whenever they disagreed (in practice by a relative 1e-6 to
# 1e-4) the clamped outer radius came out fractionally above the operator's virial radius, the guard
# failed, and ram pressure stripping was skipped entirely -- the satellite retained its gas and
# formed far more stars.
#
# Because the component's getter clamps the outer radius, this regime CANNOT be reached by setting a
# large outerRadius in the tree file; such a value is silently clamped away. The model instead
# reproduces the real failure mode, giving the node operator its own darkMatterHaloScale with 8x the
# virial density contrast so that its virial radius is exactly half the component's. The regime is
# then entered deliberately and by a factor of two rather than by floating point accident.
#
# With no cooling, star formation or feedback, ram pressure stripping is the only process that can
# change the satellite's hot halo mass: without the fix the mass is exactly constant, with it the
# mass decreases.
#
# Andrew Benson

status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/cgmRamPressureOuterRadius.xml", shell=True)
if status.returncode != 0:
    print("FAILED: model run")
    sys.exit(0)
print("SUCCESS: model run")

with h5py.File("outputs/cgmRamPressureOuterRadius.hdf5", "r") as model:
    outputs = model["Outputs"]

    def read(output, name):
        return outputs[f"Output{output}/nodeData/{name}"][:]

    isSatellite  = read(1, "nodeIsIsolated") == 0
    if isSatellite.sum() != 1:
        print(f"FAILED: expected exactly 1 satellite in the first output, found {isSatellite.sum()}")
        sys.exit(0)

    radiusOuterInitial  = read(1, "hotHaloOuterRadius"        )[isSatellite][0]
    radiusVirialInitial = read(1, "darkMatterOnlyRadiusVirial")[isSatellite][0]
    massInitial         = read(1, "hotHaloMass"               )[isSatellite][0]

    isSatelliteFinal = read(2, "nodeIsIsolated") == 0
    if isSatelliteFinal.sum() != 1:
        print(f"FAILED: expected exactly 1 satellite in the final output, found {isSatelliteFinal.sum()}")
        sys.exit(0)
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
    sys.exit(0)

if radiusOuterInitial > radiusVirialOperator:
    print(f"SUCCESS: satellite CGM outer radius exceeds the virial radius seen by the ram pressure "
          f"stripping operator ({radiusOuterInitial:.6f} > {radiusVirialOperator:.6f} Mpc)")
else:
    print(f"FAILED: test precondition not met -- satellite CGM outer radius does not exceed the "
          f"virial radius seen by the ram pressure stripping operator "
          f"({radiusOuterInitial:.6f} <= {radiusVirialOperator:.6f} Mpc); "
          f"the model no longer exercises the regime this test targets")
    sys.exit(0)

# The actual regression check.
if massFinal < massInitial:
    print(f"SUCCESS: satellite CGM was ram pressure stripped "
          f"({massInitial:.4e} -> {massFinal:.4e} M_sun, "
          f"{100.0*(1.0-massFinal/massInitial):.2f}% removed)")
else:
    print(f"FAILED: satellite CGM was not ram pressure stripped "
          f"({massInitial:.4e} -> {massFinal:.4e} M_sun); ram pressure stripping is being skipped "
          f"for a satellite whose CGM outer radius exceeds its virial radius")
    sys.exit(0)
