#!/usr/bin/env python3
import numpy as np
import h5py
import subprocess
import os
import sys

# Validate the mass profile produced by the `particulate` merger tree operator.
# Andrew Benson (27-February-2026)

# Create output path.
try:
    os.mkdir("testSuite/outputs")
except FileExistsError:
    pass

# Run the particulate model.
status = subprocess.run("./Galacticus.exe testSuite/regressions/particulate.xml",shell=True)
if status.returncode != 0:
   print("FAILED: particulate test model failed to run")
   sys.exit()

# Read the particle data.
fileParticles = h5py.File('testSuite/outputs/particulateParticles.hdf5','r')
position      = fileParticles['PartType1/Coordinates'][:]
IDs           = fileParticles['PartType1/ParticleIDs'][:]
massTable     = fileParticles['Header'].attrs['MassTable']

# Convert from Gadget to Galacticus units.
position  /= 1.0e03
massTable *= 1.0e10

# Read the mass profile.
model         = h5py.File('testSuite/outputs/particulate.hdf5','r')
nodes         = model['Outputs/Output1/nodeData']
haloIDs       = nodes['nodeIndex'        ][:]
mass          = nodes['massProfile'      ][:]
radius        = nodes['massProfileRadius'][:]

# Test the particle mass profile, iterating over halos.
for i in range(len(haloIDs)):
    # Select all particles in this halo.
    haloID         = haloIDs[i]
    selection      = (IDs >= haloID*10000000) & (IDs < (haloID+1)*10000000)
    # Find the radius of the particles relative to the halo center.
    x              = position[selection,0]-np.mean(position[selection,0])
    y              = position[selection,1]-np.mean(position[selection,1])
    z              = position[selection,2]-np.mean(position[selection,2])
    r              = np.sqrt(x**2+y**2+z**2)
    # Construct the mass profile.
    massParticles  = np.zeros(radius.shape[1])
    countParticles = np.zeros(radius.shape[1])
    for j in range(radius.shape[1]):
        selectRadius     = r < radius[i,j]
        massParticles[j] = np.count_nonzero(selectRadius)*float(massTable[0])
    # Estimate the Poisson noise level.
    noiseLevel = np.sqrt(mass[i,j]/float(massTable[0]))*float(massTable[0])
    # Check that mass profiles agree with 0.1% and 3σ.
    if np.any((countParticles > 0.0) & (np.abs(massParticles-mass[i,:]) > 1.0e-3*mass[i,:]) & (np.abs(massParticles-mass[i,:]) > 3.0*noiseLevel)):
        print(f"FAIL: mass profile does not match for halo {haloID}")
    else:
        print(f"success: mass profile matches for halo {haloID}")
