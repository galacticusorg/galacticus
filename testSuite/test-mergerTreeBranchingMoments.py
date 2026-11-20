#!/usr/bin/env python3
import os
import subprocess
import sys
import h5py
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

# Test the mergerTreeBranchingMoments task, which solves the equations for the mean and variance of counts of progenitors and
# subhalos in merger trees, and decomposes the variance into the contributions of splits. The task is run with a white-noise power
# spectrum and the EPS branching rates (the Parkinson, Cole and Helly rates with G₀=1, γ₁=γ₂=0), for which the binary split kernel
# is symmetric and the process is exactly that of extended Press-Schechter theory. The test checks that:
#
#  * the mean count of progenitors above each threshold at the target redshift (statistic A) matches the EPS conditional mass
#    function;
#  * the density of branches matches the EPS conditional mass function per unit ln M, away from the ends of the mass grid;
#  * the identities which relate the forward densities to the solution hold (to rounding error for statistics B and C, and to the
#    accuracy of its scheme for statistic A), and the contributions of splits sum to the variance of the count;
#  * statistic C has contributions only from splits of the main branch;
#  * the output on one thread is identical to that on many;
#  * the moments are unchanged when the decomposition of the variance is switched off.
#
# Tolerances are set from the measured accuracy of the solver on this grid (50 points per decade of mass, 1000 steps in ω):
# statistic A matches EPS to 3e-5, the density to 2e-3 (where the EPS result is itself approximate, as splits below the
# resolution are treated as continuous loss of mass), and the sums for statistic A to 6e-3.
# Andrew Benson, Claude (18-November-2026)

failed = False
def fail(message):
    global failed
    print("FAILED: "+message)
    failed = True

def run(label,changes,environment=None):
    # Run the model from the root of the repository, logging its output, and report whether it succeeded.
    print("Running "+label+"...")
    with open("outputs/mergerTreeBranchingMoments"+label+".log","w") as log:
        status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeBranchingMomentsWhiteNoise.xml "+changes,stdout=log,stderr=log,shell=True,env=environment)
    print("...done ("+str(status.returncode)+")")
    if status.returncode != 0 or subprocess.run("grep -q -i -e fatal -e aborted outputs/mergerTreeBranchingMoments"+label+".log",shell=True).returncode == 0:
        print("FAILED: "+label+" run:")
        subprocess.run("tail -n 50 outputs/mergerTreeBranchingMoments"+label+".log",shell=True)
        sys.exit(0)

subprocess.run("mkdir -p outputs",shell=True)
subprocess.run("rm -f outputs/mergerTreeBranchingMomentsWhiteNoise*.hdf5",shell=True)
# The run without the decomposition goes first. A run which tabulates σ(M) and the critical overdensity afresh differs, by
# rounding error, from one which reads those tables from the files stored by an earlier run, so the two runs compared bit for bit
# (on many threads and on one) must both follow a run which has stored them.
run("WhiteNoiseNoDecomposition","testSuite/parameters/mergerTreeBranchingMomentsWhiteNoiseNoDecomposition.xml")
run("WhiteNoise"               ,""                                                                       )
run("WhiteNoiseSerial"         ,"testSuite/parameters/mergerTreeBranchingMomentsWhiteNoiseSerial.xml"   ,dict(os.environ,OMP_NUM_THREADS="1"))

moments                = h5py.File("outputs/mergerTreeBranchingMomentsWhiteNoise.hdf5"               ,"r")["branchingMoments"]
momentsSerial          = h5py.File("outputs/mergerTreeBranchingMomentsWhiteNoiseSerial.hdf5"         ,"r")["branchingMoments"]
momentsNoDecomposition = h5py.File("outputs/mergerTreeBranchingMomentsWhiteNoiseNoDecomposition.hdf5","r")["branchingMoments"]

# The EPS conditional mass function per unit ln M, dN/dln M = (M₀/M) f(S|S₀,Δω) |dS/dln M|, with
# f(S|S₀,Δω) = Δω exp(-Δω²/2ΔS)/√(2π ΔS³) and ΔS = S-S₀, from the mass variance on the grid of the output.
mass          = moments["mass"][:]
logMass       = np.log(mass)
logVariance   = CubicSpline(logMass,np.log(moments["rootVariance"][:]**2))
massRoot      = moments.attrs["massRoot"]
massResolution= moments.attrs["massResolution"]
varianceRoot  = np.exp(logVariance(np.log(massRoot)))
def massFunctionEPS(logMass_,omegaInterval):
    variance           = np.exp(logVariance(logMass_))
    varianceDifference = variance-varianceRoot
    return (massRoot/np.exp(logMass_))*omegaInterval/np.sqrt(2.0*np.pi)/varianceDifference**1.5*np.exp(-omegaInterval**2/2.0/varianceDifference)*(-variance*logVariance(logMass_,1))

# Statistic A: the mean count at the root is the EPS conditional mass function integrated above the threshold.
toleranceMassFunction = 1.0e-4
omegaInterval = moments.attrs["omegaTarget"]-moments.attrs["omegaRoot"]
for name in sorted(moments["statisticA"]):
    case       = moments["statisticA/"+name]
    threshold  = case.attrs["massThreshold"]
    countEPS   = quad(massFunctionEPS,np.log(threshold),np.log(massRoot),args=(omegaInterval,),limit=200,epsabs=0.0,epsrel=1.0e-10)[0]
    difference = case.attrs["countMeanRoot"]/countEPS-1.0
    print("statistic A, m = {:.1e}: ⟨N⟩ = {:.6f}, EPS {:.6f} (fractional difference {:+.2e})".format(threshold,case.attrs["countMeanRoot"],countEPS,difference))
    if abs(difference) > toleranceMassFunction:
        fail("statistic A, m = {:.1e}: ⟨N⟩ differs from the EPS conditional mass function".format(threshold))

# The density of branches, compared away from the ends of the mass grid (the narrow peak of the root, and the resolution), at ω
# where it is substantial and the treatment of splits below the resolution has had little effect.
toleranceDensity = 5.0e-3
omega            = moments["omega"][:]
density          = moments["densityBranches"][:]
selectMass       = (mass > 4.0*massResolution) & (mass < 0.25*massRoot)
differenceMaximum = 0.0
for level in range(len(omega)):
    omegaInterval = omega[level]-omega[0]
    if omegaInterval < 1.0 or omegaInterval > 6.0:
        continue
    densityEPS        = massFunctionEPS(logMass[selectMass],omegaInterval)
    differenceMaximum = max(differenceMaximum,np.max(np.abs(density[level,selectMass]/densityEPS-1.0)))
print("density of branches: maximum fractional difference from EPS {:.2e}".format(differenceMaximum))
if differenceMaximum > toleranceDensity:
    fail("density of branches differs from the EPS conditional mass function")

# Identities and sums. For statistics B and C the densities are the exact transpose of the scheme, so these hold to rounding error;
# statistic A is solved differently near the boundary of its support, so they hold to the accuracy of its scheme.
toleranceExact       = 1.0e-10
toleranceApproximate = 2.0e-2
for statistic in "ABC":
    tolerance = toleranceApproximate if statistic == "A" else toleranceExact
    for name in sorted(moments["statistic"+statistic]):
        case          = moments["statistic"+statistic+"/"+name]
        countMean     = case.attrs["countMeanRoot"    ]
        countVariance = case.attrs["countVarianceRoot"]
        differences   = {
            "mean identity"            :     case.attrs["identityResidualMean"    ]         /countMean,
            "variance identity"        :     case.attrs["identityResidualVariance"]         /countVariance     ,
            "sum by mass"              : abs(case["varianceContributionsMass"     ][:].sum()/countVariance-1.0),
            "sum by mass ratio"        : abs(case["varianceContributionsMassRatio"][:].sum()/countVariance-1.0),
            "final cumulative variance": abs(case["varianceCumulative"            ][-1]     /countVariance-1.0),
        }
        label = "statistic {:s}, m = {:.1e}".format(statistic,case.attrs["massThreshold"])
        print(label+": "+", ".join("{:s} {:.1e}".format(quantity,difference) for quantity, difference in differences.items()))
        for quantity, difference in differences.items():
            if difference > tolerance:
                fail(label+": "+quantity+" fails")
        if statistic == "C" and np.any(case["varianceContributionsMass"][2:4] != 0.0):
            fail(label+": side branches contribute to the variance")

# The output on one thread must be identical to that on many.
differing = []
def compare(name,object_):
    if isinstance(object_,h5py.Dataset) and not np.array_equal(object_[:],momentsSerial[name][:]):
        differing.append(name)
moments.visititems(compare)
for name in moments.attrs:
    if not np.array_equal(moments.attrs[name],momentsSerial.attrs[name]):
        differing.append("@"+name)
if differing:
    fail("output on one thread differs from that on many: "+", ".join(differing))
else:
    print("output on one thread is identical to that on many")

# The moments must not depend on the decomposition of the variance (which builds the tables of the rate objects in a different
# order, so changes them by rounding error only).
toleranceDecomposition = 1.0e-9
differenceMaximum      = 0.0
for statistic in "ABC":
    for name in sorted(moments["statistic"+statistic]):
        for quantity in ("countMean","countVariance"):
            reference         = moments               ["statistic"+statistic+"/"+name+"/"+quantity][:]
            other             = momentsNoDecomposition["statistic"+statistic+"/"+name+"/"+quantity][:]
            differenceMaximum = max(differenceMaximum,np.max(np.abs(other-reference))/np.max(np.abs(reference)))
        if "varianceContributionsMass" in momentsNoDecomposition["statistic"+statistic+"/"+name]:
            fail("statistic "+statistic+": the decomposition of the variance is output when switched off")
print("moments without the decomposition of the variance: maximum fractional difference {:.1e}".format(differenceMaximum))
if differenceMaximum > toleranceDecomposition:
    fail("moments differ when the decomposition of the variance is switched off")

if not failed:
    print("SUCCESS: mergerTreeBranchingMoments task reproduces the EPS white-noise limit, and its identities hold")
sys.exit(0)
