#!/usr/bin/env python3
import subprocess
import sys
import os
import h5py
import numpy as np
import warnings

# Check the splashback radii output by the `haloMassFunction` task. The fitting functions themselves are validated against
# reference values from the Colossus toolkit by the `tests.splashback_radii.exe` unit test - this test checks that they are
# correctly wired into the task, and that the resulting radii and masses behave sensibly as a function of halo mass.
# Andrew Benson with assistance from Claude (07-August-2026).

warnings.filterwarnings('ignore')

print("Running splashback radius model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-splashback-radius.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/splashbackRadius.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("   FAILED: model run:")
    subprocess.run("cat outputs/test-splashback-radius.log",shell=True)
    sys.exit(0)
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"task failed\" -e \"unrecognized parameter\" -e \"Galacticus experienced an error in the GSL library\" outputs/test-splashback-radius.log",shell=True)
if status.returncode == 0:
    print("   FAILED: model run (errors):")
    subprocess.run("cat outputs/test-splashback-radius.log",shell=True)
    sys.exit(0)
print("   SUCCESS: model run")

model = h5py.File('outputs/splashbackRadius.hdf5','r')
outputs = model['Outputs']

for outputName in sorted(outputs.keys()):
    output = outputs[outputName]
    redshift = output.attrs['outputRedshift']
    label = outputName+" (z="+str(round(float(redshift),2))+")"

    # Check that the datasets are present.
    missing = [name for name in ('haloSplashbackRadius','haloSplashbackMass') if name not in output]
    if missing:
        print("   FAILED: "+label+": missing dataset(s): "+", ".join(missing))
        sys.exit(0)
    print("   SUCCESS: "+label+": splashback datasets present")

    massHalo         = output['haloMass'            ][:]
    radiusVirial     = output['haloVirialRadius'    ][:]
    radiusSplashback = output['haloSplashbackRadius'][:]
    massSplashback   = output['haloSplashbackMass'  ][:]

    # Check that all values are finite and positive.
    if not (np.all(np.isfinite(radiusSplashback)) and np.all(radiusSplashback > 0.0)):
        print("   FAILED: "+label+": splashback radii are not all finite and positive")
        sys.exit(0)
    if not (np.all(np.isfinite(massSplashback  )) and np.all(massSplashback   > 0.0)):
        print("   FAILED: "+label+": splashback masses are not all finite and positive")
        sys.exit(0)
    print("   SUCCESS: "+label+": splashback radii and masses are finite and positive")

    # Check that splashback radii and masses are within a physically plausible range of the virial values. The splashback
    # radius lies outside the virial radius, but by a bounded factor.
    ratioRadius = radiusSplashback/radiusVirial
    ratioMass   = massSplashback  /massHalo
    if not np.all((ratioRadius > 0.7) & (ratioRadius < 4.0)):
        print("   FAILED: "+label+": Rsp/Rvir outside plausible range: "+str(ratioRadius.min())+" to "+str(ratioRadius.max()))
        sys.exit(0)
    if not np.all((ratioMass   > 0.5) & (ratioMass   < 3.0)):
        print("   FAILED: "+label+": Msp/Mvir outside plausible range: "+str(ratioMass.min())+" to "+str(ratioMass.max()))
        sys.exit(0)
    print("   SUCCESS: "+label+": Rsp/Rvir and Msp/Mvir are in plausible ranges")

    # The splashback radius must increase with halo mass.
    if not np.all(np.diff(radiusSplashback) > 0.0):
        print("   FAILED: "+label+": splashback radius is not increasing with halo mass")
        sys.exit(0)
    # More massive halos accrete more rapidly, which shrinks the splashback radius relative to the virial radius. We test the
    # overall trend across the mass range rather than requiring strict point-to-point monotonicity: at high redshift the most
    # massive halos here reach peak heights at the edge of the range over which the Diemer (2020) fitting function was
    # calibrated, where that function flattens and can turn over slightly.
    if not ratioRadius[0] > ratioRadius[-1]:
        print("   FAILED: "+label+": Rsp/Rvir does not decrease across the halo mass range")
        sys.exit(0)
    if not np.all(np.diff(ratioRadius) < 0.02*ratioRadius[:-1]):
        print("   FAILED: "+label+": Rsp/Rvir increases significantly with halo mass somewhere in the range")
        sys.exit(0)
    # The mass accretion rate must increase with halo mass, which is what drives the trend above.
    if not np.all(np.diff(output['haloMassAccretionRate'][:]) > 0.0):
        print("   FAILED: "+label+": mass accretion rate is not increasing with halo mass")
        sys.exit(0)
    print("   SUCCESS: "+label+": splashback radius scales correctly with halo mass")

model.close()
