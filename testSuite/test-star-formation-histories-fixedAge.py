#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check internal self-consistency of fixed age star formation histories.
# Andrew Benson (07-October-2024)

# Run the model and check for completion                                                                                                                                                                                 
print("Running model...")
status = subprocess.run("mkdir -p outputs/test-star-formation-histories-fixedAge",shell=True)
log = open("outputs/test-star-formation-histories-fixedAge/galacticus.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/test-star-formation-histories-fixedAge.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-star-formation-histories-fixedAge/galacticus.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-star-formation-histories-fixedAge/galacticus.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run("cat outputs/test-star-formation-histories-fixedAge/galacticus.log",shell=True)
    sys.exit()
print("SUCCESS: model run")

# Open the model and extract the recycled fraction.
model             = h5py.File('outputs/test-star-formation-histories-fixedAge/galacticus.hdf5','r')
stellarPopulation = model['Parameters/stellarPopulation']
recycledFraction  = stellarPopulation.attrs['recycledFraction']

# Check that the star formation histories integrate to the expected total stellar mass.
# Get the nodeData group.
nodes                           = model["Lightcone/Output1/nodeData"]
# Read stellar masses and star formation histories.
nodesTime                            = nodes['time'                             ][:]
nodesMassStellarDisk                 = nodes['diskMassStellar'                  ][:]
nodesMassStellarSpheroid             = nodes['spheroidMassStellar'              ][:]
sfhDiskStarFormationHistory          = nodes['diskStarFormationHistoryMass'     ][:]
sfhSpheroidStarFormationHistory      = nodes['spheroidStarFormationHistoryMass' ][:]
sfhDiskStarFormationHistoryTimes     = nodes['diskStarFormationHistoryTimes'    ][:]
sfhSpheroidStarFormationHistoryTimes = nodes['spheroidStarFormationHistoryTimes'][:]
# Sum star formation history masses over times and metallicities.
nodesSFHIntegratedDisk     = np.array(list(map(lambda x: np.sum(x)*(1.0-recycledFraction),list(map(lambda x: sum(x),sfhDiskStarFormationHistory    )))))
nodesSFHIntegratedSpheroid = np.array(list(map(lambda x: np.sum(x)*(1.0-recycledFraction),list(map(lambda x: sum(x),sfhSpheroidStarFormationHistory)))))
# Extract final time of each history.
nodesSFHTimeFinalDisk      = np.array(list(map(lambda x,y: x[-1] if len(x) > 0 else y,sfhDiskStarFormationHistoryTimes    ,nodesTime)))
nodesSFHTimeFinalSpheroid  = np.array(list(map(lambda x,y: x[-1] if len(x) > 0 else y,sfhSpheroidStarFormationHistoryTimes,nodesTime)))
# Test accuracy of star formation histories.
tolerance      = 1.0e-3
statusDisk     = "FAILED" if any(abs(nodesSFHIntegratedDisk    -nodesMassStellarDisk    ) > tolerance*nodesMassStellarDisk    ) else "SUCCESS"
statusSpheroid = "FAILED" if any(abs(nodesSFHIntegratedSpheroid-nodesMassStellarSpheroid) > tolerance*nodesMassStellarSpheroid) else "SUCCESS"
print(" -> "+statusDisk    +": disk stellar mass"    )
print(" -> "+statusSpheroid+": spheroid stellar mass")
# Test that the output time matches the final time in star formation histories.
tolerance      = 1.0e-6
statusDisk     = "FAILED" if any(abs(nodesSFHTimeFinalDisk    -nodesTime) > tolerance*nodesTime) else "SUCCESS"
statusSpheroid = "FAILED" if any(abs(nodesSFHTimeFinalSpheroid-nodesTime) > tolerance*nodesTime) else "SUCCESS"
print(" -> "+statusDisk    +": disk final time"    )
print(" -> "+statusSpheroid+": spheroid final time")
