#!/usr/bin/env python3
import os
import sys
import h5py
import numpy as np

# Check internal self-consistency of adaptive star formation histories.
# Andrew Benson (25-March-2021)

# Run the model and check for completion
if os.system("mkdir -p outputs/test-star-formation-histories-adapative; cd ..; ./Galacticus.exe testSuite/parameters/test-star-formation-histories-adaptive.xml"):
    print("FAILED: model run:")
    sys.exit()
else:
    print("SUCCESS: model run")

# Open the model and extract the recycled fraction.
model             = h5py.File('outputs/test-star-formation-histories-adapative/galacticus.hdf5','r')
outputs           = model['Outputs'                     ]
stellarPopulation = model['Parameters/stellarPopulation']
recycledFraction  = stellarPopulation.attrs['recycledFraction']

# Iterate over outputs to check that the star formation histories integrates to the expected total stellar mass.
for output in outputs.keys():
    print(output+":")
    # Get the nodeData group.
    nodes                           = outputs[output+"/nodeData"]
    # Read stellar masses and star formation histories.
    nodesMassStellarDisk            = nodes['diskMassStellar'                 ][:]
    nodesMassStellarSpheroid        = nodes['spheroidMassStellar'             ][:]
    sfhDiskStarFormationHistory     = nodes['diskStarFormationHistoryMass'    ][:]
    sfhSpheroidStarFormationHistory = nodes['spheroidStarFormationHistoryMass'][:]
    # Sum star formation history masses over times and metallicities.
    nodesSFHIntegratedDisk     = np.array(list(map(lambda x: np.sum(x)*(1.0-recycledFraction),list(map(lambda x: sum(x),sfhDiskStarFormationHistory    )))))
    nodesSFHIntegratedSpheroid = np.array(list(map(lambda x: np.sum(x)*(1.0-recycledFraction),list(map(lambda x: sum(x),sfhSpheroidStarFormationHistory)))))
    # Test accuracy of star formation histories.
    tolerance      = 1.0e-3
    statusDisk     = "FAILED" if any(abs(nodesSFHIntegratedDisk    -nodesMassStellarDisk    ) > tolerance*nodesMassStellarDisk    ) else "SUCCESS"
    statusSpheroid = "FAILED" if any(abs(nodesSFHIntegratedSpheroid-nodesMassStellarSpheroid) > tolerance*nodesMassStellarSpheroid) else "SUCCESS"
    print(" -> "+statusDisk    +": disk stellar mass"    )
    print(" -> "+statusSpheroid+": spheroid stellar mass")
