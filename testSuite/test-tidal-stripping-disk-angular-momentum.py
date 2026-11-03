#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np
from scipy.special import gammainc, gammaincc, gammaln
from scipy.optimize import brentq

# Check that the specific angular momentum of a tidally stripped disk remains consistent with its truncated mass distribution.
#
# When the mass distribution of a disk is truncated (by `galacticComponentTruncationExponential`), the material stripped from it
# carries the specific angular momentum implied by the change in its truncated mass distribution, j/<j> = 1 - d ln<R>/dx, where x
# is the degree of truncation of the disk (such that it retains a fraction exp(-x) of the mass of its untruncated distribution),
# and <R> is its mass-weighted mean radius. Since the degree of truncation increases as dx = -dM/M, if nothing other than tidal
# stripping changes the disk then its mean specific angular momentum evolves as d ln<j> = d ln<R>, and so
#
#     (<j>/<j>_0) / (<R>_x/<R>_0) = 1
#
# throughout its evolution, where the subscript 0 indicates the untruncated disk. This test runs the single-subhalo example model
# (parameters/examples/tidalStrippingBaryonicSubhalo.xml), in which gas, star formation, and feedback are omitted so that only
# tidal effects change the disk, and checks this relation at each output. The mean radius of the truncated disk is computed here
# independently of Galacticus, using incomplete gamma functions. The relation is a consequence of the differential equations for
# the mass, angular momentum, and degree of truncation of the disk, and does not depend on the structure solver or on the rotation
# curve of the disk. It does not depend on heating of the disk either, which conserves its mean radius.
#
# Andrew Benson, Claude

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Run the model.
log    = open("outputs/test-tidal-stripping-disk-angular-momentum.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe parameters/examples/tidalStrippingBaryonicSubhalo.xml testSuite/parameters/tidalStrippingDiskAngularMomentum.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("FAILED: model failed to run")
    sys.exit(0)

# Read the model.
model = h5py.File("outputs/tidalStrippingDiskAngularMomentum.hdf5","r")
if model['Parameters'].attrs['galacticComponentTruncation'].decode() != "exponential":
    print("FAILED: the model must use exponential truncation")
    sys.exit(0)
radiusDecayOverRadiusTruncate = model['Parameters/galacticComponentTruncation'].attrs['radiusDecayOverRadiusTruncate']

# Find the degree of truncation, mass, and angular momentum of the disk of each satellite at each output.
outputs = sorted(model['Outputs'].keys(),key=lambda name: int(name[6:]))
history = {}
for output in outputs:
    nodeData = model['Outputs'][output]['nodeData']
    if 'nodeIndex' not in nodeData:
        continue
    for i in np.flatnonzero(nodeData['nodeIsIsolated'][:] == 0):
        if nodeData['diskMassGas'][i] != 0.0:
            print("FAILED: the disk must contain no gas")
            sys.exit(0)
        mass = nodeData['diskMassStellar'][i]
        if mass <= 0.0:
            continue
        history.setdefault(nodeData['nodeIndex'][i],[]).append((nodeData['diskTidalTruncation'][i],mass,nodeData['diskAngularMomentum'][i]))
model.close()
if len(history) == 0:
    print("FAILED: no satellite disks were found")
    sys.exit(0)

# Moments, I_n = ∫_0^∞ Σ(R) Rⁿ dR, of the surface density of an exponential disk of unit mass and unit scale length truncated at
# radius R_t, beyond which the surface density is multiplied by (R/R_t)^b exp[-b (R/R_t-1)], where b=1/radiusDecayOverRadiusTruncate.
b = 1.0/radiusDecayOverRadiusTruncate
def moment(n,radiusTruncate):
    a     = 1.0+b/radiusTruncate
    s     = n+b+1.0
    inner = np.exp(gammaln(n+1.0))*gammainc(n+1.0,radiusTruncate)
    outer = np.exp(-b*np.log(radiusTruncate)+b-s*np.log(a)+gammaln(s))*gammaincc(s,a*radiusTruncate)
    return inner+outer

# Galacticus holds the shape of the truncated disk fixed for truncation radii below 0.01 scale lengths.
radiusTruncateMinimum = 1.0e-2
truncationMaximum     = -np.log(moment(1,radiusTruncateMinimum))

# Ratio of the mean radius of the truncated disk to that of the untruncated disk (which is 2 scale lengths).
def ratioRadiusMean(truncation):
    if truncation <= 0.0:
        return 1.0
    truncation     = min(truncation,truncationMaximum)
    radiusTruncate = brentq(lambda radius: -np.log(moment(1,radius))-truncation,radiusTruncateMinimum,100.0,xtol=1.0e-14,rtol=1.0e-14)
    return moment(2,radiusTruncate)/moment(1,radiusTruncate)/2.0

# Check the relation for each satellite disk.
toleranceRelative  = 5.0e-3
truncationRequired = 1.0
for nodeIndex, states in history.items():
    truncation, mass, angularMomentum = map(np.array,zip(*states))
    if truncation[0] != 0.0:
        print(f"FAILED: disk of node {nodeIndex} is already truncated when it first becomes a satellite")
        sys.exit(0)
    ratioAngularMomentum = (angularMomentum/mass)/(angularMomentum[0]/mass[0])
    ratioRadius          = np.array([ratioRadiusMean(x) for x in truncation])
    error                = np.max(np.abs(ratioAngularMomentum/ratioRadius-1.0))
    status               = "SUCCESS" if error < toleranceRelative else "FAILED"
    print(f"{status}: disk of node {nodeIndex}, truncated to x={truncation.max():.2f} (<R>_x/<R>_0={ratioRadius.min():.3f}): maximum |(<j>/<j>_0)/(<R>_x/<R>_0)-1| = {error:.2e} (tolerance {toleranceRelative:.0e})")
    # Require substantial truncation, so that the test is sensitive to the specific angular momentum of the stripped material.
    if truncation.max() < truncationRequired:
        print(f"FAILED: disk of node {nodeIndex} was truncated only to x={truncation.max():.2f} - at least x={truncationRequired:.1f} is required for a meaningful test")
