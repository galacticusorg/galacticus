#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check internal self-consistency of in situ star formation histories.
# Andrew Benson (17-July-2024)

# Run the model and check for completion                                                                                                                                                                                 
print("Running model...")
status = subprocess.run("mkdir -p outputs/test-star-formation-histories-inSitu",shell=True)
log = open("outputs/test-star-formation-histories-inSitu/galacticus.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/test-star-formation-histories-inSitu.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-star-formation-histories-inSitu/galacticus.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-star-formation-histories-inSitu/galacticus.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run("cat outputs/test-star-formation-histories-inSitu/galacticus.log",shell=True)
    sys.exit()
print("SUCCESS: model run")

# Open the model and extract the recycled fraction.
model             = h5py.File('outputs/test-star-formation-histories-inSitu/galacticus.hdf5','r')
outputs           = model['Outputs'                     ]
stellarPopulation = model['Parameters/stellarPopulation']
recycledFraction  = stellarPopulation.attrs['recycledFraction']

# Iterate over outputs to check that the star formation histories integrates to the expected total stellar mass.
for output in outputs.keys():
    print(output+":")
    # Get the nodeData group.
    nodes                            = outputs[output+"/nodeData"]
    # Read stellar masses and star formation histories.
    nodesMassStellarDisk             = nodes['diskMassStellar'                 ][:]
    nodesMassStellarSpheroid         = nodes['spheroidMassStellar'             ][:]
    sfhDiskStarFormationHistory      = nodes['diskStarFormationHistoryMass'    ][:]
    sfhSpheroidStarFormationHistory  = nodes['spheroidStarFormationHistoryMass'][:]
    # Sum star formation history masses over times and metallicities. For in-situ star formation histories the first half of the
    # columns correspond to the total star formation history.
    nodesSFHIntegratedDisk           = np.array(list(map(lambda x: np.sum(x)*(1.0-recycledFraction),list(map(lambda x: sum(x[:int(len(x)/2) ]),sfhDiskStarFormationHistory    )))))
    nodesSFHIntegratedSpheroid       = np.array(list(map(lambda x: np.sum(x)*(1.0-recycledFraction),list(map(lambda x: sum(x[:int(len(x)/2) ]),sfhSpheroidStarFormationHistory)))))
    # The second half correspond to the in-situ component.
    nodesSFHIntegratedDiskInSitu     = np.array(list(map(lambda x: np.sum(x)*(1.0-recycledFraction),list(map(lambda x: sum(x[-int(len(x)/2):]),sfhDiskStarFormationHistory    )))))
    nodesSFHIntegratedSpheroidInSitu = np.array(list(map(lambda x: np.sum(x)*(1.0-recycledFraction),list(map(lambda x: sum(x[-int(len(x)/2):]),sfhSpheroidStarFormationHistory)))))
    # Test accuracy of total star formation histories.
    tolerance      = 1.0e-3
    statusDisk     = "FAILED" if any(abs(nodesSFHIntegratedDisk    -nodesMassStellarDisk    ) > tolerance*nodesMassStellarDisk    ) else "SUCCESS"
    statusSpheroid = "FAILED" if any(abs(nodesSFHIntegratedSpheroid-nodesMassStellarSpheroid) > tolerance*nodesMassStellarSpheroid) else "SUCCESS"
    print(" -> "+statusDisk    +": total disk stellar mass"    )
    print(" -> "+statusSpheroid+": total spheroid stellar mass")
    # Test in-situ star formation histories.
    tolerance      = 1.0e-3
    statusDisk     = "FAILED" if any(nodesSFHIntegratedDiskInSitu     > nodesSFHIntegratedDisk    ) else "SUCCESS"
    statusSpheroid = "FAILED" if any(nodesSFHIntegratedSpheroidInSitu > nodesSFHIntegratedSpheroid) else "SUCCESS"
    print(" -> "+statusDisk    +": in-situ disk stellar mass"    )
    print(" -> "+statusSpheroid+": in-situ spheroid stellar mass")
