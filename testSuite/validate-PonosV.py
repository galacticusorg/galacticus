#!/usr/bin/env python3
import sys
import os
import subprocess
import numpy as np
import h5py
import json
import codecs
from git import Repo

# Run models to validate subhalo projected density against the PonosV simulation of Fiacconi et al. (2016; https://ui.adsabs.harvard.edu/abs/2016ApJ...824..144F).
# Andrew Benson (13-February-2024)

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the validate model.
status = subprocess.run("cd ..; export OMP_NUM_THREADS=2; ./Galacticus.exe testSuite/parameters/validate_PonosV.xml",shell=True)
if status.returncode != 0:
   print("FAILED: PonosV validation model failed to run")
   sys.exit()

# Read data.
model       = h5py.File("outputs/validate_PonosV.hdf5","r")
outputs     = model  ['Outputs'         ]
redshift0p0 = outputs['Output2/nodeData']
redshift0p7 = outputs['Output1/nodeData']
data = {
    "0.0": {},
    "0.7": {}
}
for propertyName in 'nodeIsIsolated', 'darkMatterOnlyRadiusVirial':
    data['0.0'][propertyName] = redshift0p0[propertyName][:]
for propertyName in 'nodeIsIsolated', 'hostDarkMatterOnlyRadiusVirial', 'positionOrbitalX', 'positionOrbitalY', 'positionOrbitalZ', 'satelliteBoundMass', 'mergerTreeIndex':
    data['0.7'][propertyName] = redshift0p7[propertyName][:]

# Validate z=0.0 host halo virial radii.
hostsFinal         = np.nonzero(data['0.0']['nodeIsIsolated'] == 1)
radiusVirialTarget = 0.6005 # Virial radius of PonosV from Table 1 of Fiacconi et al. (2016; https://ui.adsabs.harvard.edu/abs/2016ApJ...824..144F).
offsetFractional   = (data['0.0']['darkMatterOnlyRadiusVirial'][hostsFinal]-radiusVirialTarget)/radiusVirialTarget
if np.any(offsetFractional > 0.01):
    print("FAIL: PonosV z=0.0 host virial radii")
    print("   Expected: "+str(radiusVirialTarget))
    print("   Found: "+str(data['0.0']['darkMatterOnlyRadiusVirial'][hostsFinal]))
else:
    print("SUCCESS: PonosV z=0.0 host virial radii")

# Select z=0.7 subhalos.
radiusFractionalMinimum = 0.00e0
radiusFractionalMaximum = 0.04e0
massBoundMinimum8       = 0.50e8
massBoundMaximum8       = 2.00e8
massBoundMinimum9       = 0.50e9
massBoundMaximum9       = 2.00e9
kilo                    = 1.00e3
data['0.7']['radiusOrbital2D'] = np.sqrt(+data['0.7']['positionOrbitalX']**2+data['0.7']['positionOrbitalY']**2                                   )
data['0.7']['radiusOrbital3D'] = np.sqrt(+data['0.7']['positionOrbitalX']**2+data['0.7']['positionOrbitalY']**2+data['0.7']['positionOrbitalZ']**2)
selection = np.nonzero(
    (data['0.7']['nodeIsIsolated'    ]               == 0                                                                    ) # ] Select subhalos.
    &              
    (data['0.7']['radiusOrbital3D'   ]               <=                         data['0.7']['hostDarkMatterOnlyRadiusVirial']) # ] Select subhalos within the host virial radius.
    &              
    (data['0.7']['radiusOrbital2D'   ]               >= radiusFractionalMinimum*data['0.7']['hostDarkMatterOnlyRadiusVirial']) # ⎫ 
    &                                                                                                                          # ⎬ Select subhalos close to projected radius of 0.02 of host virial radius.
    (data['0.7']['radiusOrbital2D'   ]               <= radiusFractionalMaximum*data['0.7']['hostDarkMatterOnlyRadiusVirial']) # ⎭
)
selection8 = np.nonzero(
    (data['0.7']['satelliteBoundMass'][selection] >= massBoundMinimum8                                                       ) # ⎫
    &                                                                                                                          # ⎬ Select subhalos close to a bound mass of 10⁸M☉.
    (data['0.7']['satelliteBoundMass'][selection] <= massBoundMaximum8                                                       ) # ⎭
)
selection9 = np.nonzero(
    (data['0.7']['satelliteBoundMass'][selection] >= massBoundMinimum9                                                       ) # ⎫
    &                                                                                                                          # ⎬ Select subhalos close to a bound mass of 10⁹M☉.
    (data['0.7']['satelliteBoundMass'][selection] <= massBoundMaximum9                                                       ) # ⎭
)

# Compute the number density of selected subhalos in each tree.
treeCount              = model['Parameters/mergerTreeBuildMasses'].attrs['treeCount'][0]
subhaloSurfaceDensity8 = np.zeros(treeCount)
subhaloSurfaceDensity9 = np.zeros(treeCount)
for i in range(treeCount):
    selectTree8 = np.nonzero(data['0.7']['mergerTreeIndex'][selection][selection8] == i+1)
    selectTree9 = np.nonzero(data['0.7']['mergerTreeIndex'][selection][selection9] == i+1)
    if len(selectTree8[0]) > 0:
        countSubhalos8            = len(selectTree8[0])
        subhaloSurfaceDensity8[i] = +(
            +countSubhalos8
            /2.0
            /np.pi
            /data['0.7']['hostDarkMatterOnlyRadiusVirial'][selection][selectTree8][0]**2
            /kilo                                                                    **2
            )  \
        /(
            +radiusFractionalMaximum**2
            -radiusFractionalMinimum**2
        ) \
        /np.log10(
            +massBoundMaximum8
            /massBoundMinimum8
        )
    else:
        subhaloSurfaceDensity8[i] = 0.0
    if len(selectTree9[0]) > 0:
        countSubhalos9            = len(selectTree9[0])
        subhaloSurfaceDensity9[i] = +(
            +countSubhalos9
            /2.0
            /np.pi
            /data['0.7']['hostDarkMatterOnlyRadiusVirial'][selection][selectTree9][0]**2
            /kilo                                                                    **2
            ) \
            /(
                +radiusFractionalMaximum**2
                -radiusFractionalMinimum**2
            ) \
            /np.log10(
                      +massBoundMaximum9
                      /massBoundMinimum9
                     )
    else:
        subhaloSurfaceDensity9[i] = 0.0

# Set target values from the PonosV simulation of Fiacconi et al. (2016; https://ui.adsabs.harvard.edu/abs/2016ApJ...824..144F).
alphaPonosV                  = 0.850
subhaloSurfaceDensity8PonosV = 0.006

# Compute percentage of realizations above/below the PonosV subhalo surface density, and report.
above = np.nonzero(subhaloSurfaceDensity8 > subhaloSurfaceDensity8PonosV)
percentageAbove = 100.0*np.count_nonzero(subhaloSurfaceDensity8 > subhaloSurfaceDensity8PonosV)/treeCount
percentageBelow = 100.0-percentageAbove
statusSurfaceDensity = "SUCCESS" if percentageAbove > 5.0 or percentageBelow > 5.0 else "FAIL"
print(f"{statusSurfaceDensity}: Percentage of realizations above/below the PonosV subhalo surface density: {percentageAbove:5.1f}/{percentageBelow:5.1f}")

# Compute the mean slope of the subhalo mass function, and report.
# We exclude models for which there are no subhalos present in one of the mass cuts. This probably introduces some bias.
nonZero     = np.nonzero((subhaloSurfaceDensity8 > 0.0) & (subhaloSurfaceDensity9 > 0.0))
alphas      = -np.log(        subhaloSurfaceDensity8[nonZero]     /        subhaloSurfaceDensity9[nonZero]     ) \
              /np.log(np.sqrt(massBoundMinimum8*massBoundMaximum8)/np.sqrt(massBoundMinimum9*massBoundMaximum9))
percentageAboveSlope = 100.0*np.count_nonzero(alphas > alphaPonosV)/len(alphas)
percentageBelowSlope = 100.0-percentageAboveSlope
statusSlope = "SUCCESS" if percentageAboveSlope > 5.0 or percentageBelowSlope > 5.0 else "FAIL"
print(f"{statusSlope}: Percentage of realizations above/below the PonosV subhalo mass function slope: {percentageAboveSlope:5.1f}/{percentageBelowSlope:5.1f}")

# Interface with git.
repo         = Repo(os.environ['GALACTICUS_EXEC_PATH'])
actor        = repo.head.commit.author
lastRevision = repo.head.object.hexsha
authorName   = actor.name
authorEmail  = actor.email
authorDate   = str(repo.head.commit.committed_datetime)
message      = repo.head.commit.message

# Generate content for the validation metrics page.
output = {
    "repoUrl"      : "https://github.com/galacticusorg/galacticus",
    "parameterFile": "testSuite/parameters/validate_PonosV.xml",
    "commit"       :
    {
        "author":
        {
            "name" : authorName,
            "email": authorEmail
        },
        "id"       : lastRevision,
        "message"  : message,
        "timestamp": authorDate,
        "url"      : "https://github.com/galacticusorg/galacticus/commit/"+lastRevision
    },
    "surfaceDensity":
    {
        "percentageBelow": percentageBelow             ,
        "percentageAbove": percentageAbove             ,
        "target"         : subhaloSurfaceDensity8PonosV
    },
    "slope":
    {
        "percentageBelow": percentageBelowSlope,
        "percentageAbove": percentageAboveSlope,
        "target"         : alphaPonosV         
    }
}
f = codecs.open("outputs/results_PonosV.json", "w", "utf-8")
f.write("window.PONOSV_DATA = ")
f.write(json.dumps(output,indent=4,ensure_ascii=False))
f.close()
