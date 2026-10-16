#!/usr/bin/env python3
"""
Test the satelliteMergingSoliton node operator.

A fuzzy dark matter model is run on a hand-built merger tree in which a single satellite falls into its host. Tidal mass
loss, tidal heating and satellite destruction are omitted so that dynamical friction carries the satellite into the soliton
cores, where this operator must merge it. At the merger the host's solitonic core mass (and its normalization) must become
a fraction 0.7 of the sum of those of the host and satellite.

Three properties of the model make this checkable from the outputs:

* The rate of change of the normalized core mass (solitonMassCoreNormal) does not depend on its current value, and follows
  the analytic core-halo mass relation of Chan et al. (2022; MNRAS; 511; 943; equation 15). So the difference between the
  host's normalized core mass and that relation is constant between merging events, and changes only at a merger.
* The core mass of a satellite is not evolved, so the satellite's values at the last output before the merger are those
  which the merger uses.
* The time of the merger is recorded in the host's merged subhalo properties.

The host's core mass (solitonMassCore) is recomputed only when its density profile is evaluated, so the value which it held
at the moment of the merger is known only approximately, and is checked against a bracket rather than exactly.
"""
import os
import subprocess
import sys
import h5py
import numpy as np

# Parameters of the core-halo mass relation from Chan et al. (2022; MNRAS; 511; 943).
alpha = 0.515
beta  = 8.0e6
gamma = 10.0**(-5.73)
# Fraction of the combined core mass retained by the host at a merger - this must match the value in
# source/nodes/operators/physics/satellite_merging/soliton.F90.
fractionMassRetained = 0.7
# Tolerances. The normalized core mass is integrated with a relative ODE tolerance of 10⁻², and its offset from the analytic
# relation drifts by up to ~2×10⁻⁴ over the ~11 Gyr following the merger. The alternatives which this test is designed to
# reject (a different retained fraction, omitting the satellite's core, or no change at the merger) all differ from the
# expectation by more than 10%.
toleranceConstancy  = 1.0e-3 # Allowed variation in the offset from the analytic relation between mergers.
toleranceMerger     = 1.0e-3 # Allowed error in the offset after the merger.
growthCoreMaximum   = 0.05   # Maximum growth in the host's core mass between the last output before the merger and the merger.

pathParameters      = os.path.abspath("parameters/satelliteMergingSoliton.xml"                                )
pathOutputDirectory = os.path.abspath("outputs/test-soliton-satellite-merging"                                )
pathOutputLog       = os.path.abspath("outputs/test-soliton-satellite-merging/test-soliton-satellite-merging.log" )
pathOutputModel     = os.path.abspath("outputs/test-soliton-satellite-merging/test-soliton-satellite-merging.hdf5")
os.makedirs(pathOutputDirectory,exist_ok=True)

# Run the model and check for completion.
print("Running model...")
with open(pathOutputLog,"w") as log:
    status = subprocess.run(f"cd ..; ./Galacticus.exe {pathParameters}",stdout=log,stderr=log,shell=True)
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run(f"cat {pathOutputLog}",shell=True)
    sys.exit(0)
print("Checking for errors...")
status = subprocess.run(f"grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" {pathOutputLog}",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run(f"cat {pathOutputLog}",shell=True)
    sys.exit(0)
print("SUCCESS: model run")

# Read the host and satellite properties at each output, in order of time.
with h5py.File(pathOutputModel,'r') as model:
    massParticle = model['Parameters/darkMatterParticle'].attrs['mass']*1.0e-22
    outputs      = sorted(model['Outputs'].values(),key=lambda output: output.attrs['outputTime'])
    time                 = np.array([output.attrs['outputTime'] for output in outputs])
    massCoreNormalHost   = np.zeros(len(outputs))
    massCoreHost         = np.zeros(len(outputs))
    massHost             = np.zeros(len(outputs))
    densityContrastHost  = np.zeros(len(outputs))
    expansionFactor      = np.zeros(len(outputs))
    massCoreNormalSatellite = np.full(len(outputs),np.nan)
    massCoreSatellite       = np.full(len(outputs),np.nan)
    timesMerger          = []
    for i, output in enumerate(outputs):
        nodes      = output['nodeData']
        isIsolated = nodes['nodeIsIsolated'][:] == 1
        if np.count_nonzero(isIsolated) != 1 or np.count_nonzero(~isIsolated) > 1:
            print(f"FAILED: expected one host and at most one satellite at t = {time[i]:.3f} Gyr")
            sys.exit(0)
        host = np.flatnonzero(isIsolated)[0]
        massCoreNormalHost [i] = nodes['solitonMassCoreNormal'][host]
        massCoreHost       [i] = nodes['solitonMassCore'      ][host]
        massHost           [i] = nodes['basicMass'            ][host]
        densityContrastHost[i] = nodes['densityContrastVirial'][host]
        expansionFactor    [i] = 1.0/(1.0+nodes['redshift'    ][host])
        timesMerger.append(nodes['mergedSubhaloTimeCurrent'][host])
        if np.any(~isIsolated):
            satellite = np.flatnonzero(~isIsolated)[0]
            massCoreNormalSatellite[i] = nodes['solitonMassCoreNormal'][satellite]
            massCoreSatellite      [i] = nodes['solitonMassCore'      ][satellite]
    # The final output is at z=0, and provides the virial density contrast at the present day needed by the analytic relation.
    if abs(expansionFactor[-1]-1.0) > 1.0e-6:
        print("FAILED: the final output is not at z=0")
        sys.exit(0)
    densityContrastPresent = densityContrastHost[-1]

# Evaluate the analytic core-halo mass relation (Chan et al. 2022; equation 15) for the host.
massCoreNormalAnalytic = (
    +beta
    *(massParticle/8.0e-23)**(-1.5)
    +(np.sqrt(densityContrastHost/densityContrastPresent)*massHost/gamma)**alpha
    *(massParticle/8.0e-23)**(1.5*(alpha-1.0))
)/np.sqrt(expansionFactor)

# Check that exactly one merger occurred, that it was recorded by this operator (the only one in the model which records merged
# subhalo properties), and that the satellite was present before it and absent after.
hasSatellite = ~np.isnan(massCoreNormalSatellite)
countMergers = np.array([len(timesMerger_) for timesMerger_ in timesMerger])
if countMergers[-1] != 1:
    print(f"FAILED: expected exactly one merger to be recorded, found {countMergers[-1]}")
    sys.exit(0)
timeMerger = timesMerger[-1][0]
premerger  = hasSatellite
postmerger = countMergers == 1
if not np.any(premerger) or not np.any(postmerger) or np.any(hasSatellite & postmerger) or time[premerger].max() >= timeMerger or time[postmerger].min() < timeMerger:
    print(f"FAILED: the satellite is not present before, and only before, the recorded merger time of {timeMerger:.4f} Gyr")
    sys.exit(0)
print(f"SUCCESS: the satellite merged at t = {timeMerger:.4f} Gyr, triggered by the satelliteMergingSoliton operator")

# The satellite's core is not evolved, so its properties must be unchanged throughout. (Its core mass is zero until first
# computed, so only computed values are compared.)
massCoreComputed = premerger & (massCoreSatellite > 0.0)
if not np.any(massCoreComputed):
    print("FAILED: the satellite's core mass was never computed")
    sys.exit(0)
massCoreNormalSatelliteMerger = massCoreNormalSatellite[premerger       ][-1]
massCoreSatelliteMerger       = massCoreSatellite      [massCoreComputed][-1]
if np.ptp(massCoreNormalSatellite[premerger]) > 0.0 or np.ptp(massCoreSatellite[massCoreComputed]) > 0.0:
    print("FAILED: the satellite's core properties changed while it was a satellite")
else:
    print("SUCCESS: the satellite's core properties were unchanged while it was a satellite")

# The offset of the host's normalized core mass from the analytic relation must be constant between the satellite's infall and
# the merger, and again after the merger. (At infall the host's mass jumps by that of the satellite while its core mass does
# not, creating an offset. Outputs begin at infall, as before it the host and satellite are both isolated halos.)
offset = massCoreNormalHost-massCoreNormalAnalytic
for label, select in (("between infall and the merger", premerger), ("after the merger", postmerger)):
    variation = np.ptp(offset[select])/massCoreNormalHost[select].min()
    if variation > toleranceConstancy:
        print(f"FAILED: offset of the host's normalized core mass from the analytic relation varies {label} (by {variation:.2e})")
    else:
        print(f"SUCCESS: offset of the host's normalized core mass from the analytic relation is constant {label} (to {variation:.2e})")
offsetPremerger  = np.mean(offset[premerger ])
offsetPostmerger = np.mean(offset[postmerger])

# Evaluate the analytic relation at the time of the merger by cubic interpolation through the four outputs nearest to it
# (drawn only from outputs after infall, across which the relation is smooth), and so find the host's normalized core mass just
# before the merger.
select                        = np.flatnonzero(time >= time[premerger].min())
select                        = select[np.argsort(np.abs(time[select]-timeMerger))[0:4]]
massCoreNormalAnalyticMerger  = np.polyfit(time[select]-timeMerger,massCoreNormalAnalytic[select],3)[-1]
massCoreNormalHostPremerger   = massCoreNormalAnalyticMerger+offsetPremerger
# Predict the offset after the merger.
massCoreNormalHostPostmerger  = fractionMassRetained*(massCoreNormalHostPremerger+massCoreNormalSatelliteMerger)
offsetPostmergerExpected      = massCoreNormalHostPostmerger-massCoreNormalAnalyticMerger
error                         = abs(offsetPostmerger-offsetPostmergerExpected)/massCoreNormalHostPostmerger
if error > toleranceMerger:
    print(f"FAILED: host's normalized core mass after the merger (offset {offsetPostmerger:.6e}) does not match {fractionMassRetained} × (host + satellite) (offset {offsetPostmergerExpected:.6e}; relative error {error:.2e})")
else:
    print(f"SUCCESS: host's normalized core mass after the merger matches {fractionMassRetained} × (host + satellite) (relative error {error:.2e})")

# The host's core mass is not recomputed after the merger (nothing evaluates its density profile), so it retains the value set
# at the merger. The host's core mass implied by that value must be at least its value at the last output before the merger
# (core masses only grow between mergers), and not much greater.
massCoreHostComputed   = premerger & (massCoreHost > 0.0)
if not np.any(massCoreHostComputed):
    print("FAILED: the host's core mass was never computed before the merger")
    sys.exit(0)
massCoreHostPostmerger = massCoreHost[postmerger          ]
massCoreHostPremerger  = massCoreHost[massCoreHostComputed][-1]
massCoreHostImplied    = massCoreHostPostmerger[0]/fractionMassRetained-massCoreSatelliteMerger
growth                 = massCoreHostImplied/massCoreHostPremerger-1.0
if np.ptp(massCoreHostPostmerger) > 0.0:
    print("FAILED: host's core mass changed after the merger although its density profile was not evaluated")
elif growth < 0.0 or growth > growthCoreMaximum:
    print(f"FAILED: host's core mass after the merger ({massCoreHostPostmerger[0]:.6e}) implies a host core mass at the merger {growth:+.2%} relative to its value at the last output before it - expected growth between 0% and {growthCoreMaximum:.0%}")
else:
    print(f"SUCCESS: host's core mass after the merger matches {fractionMassRetained} × (host + satellite) (implied growth of host core mass since the last output of {growth:+.2%})")
sys.exit(0)
