#!/usr/bin/env python3
import subprocess
import os
import sys
import h5py
import numpy as np
import xml.etree.ElementTree as ET

# Test darkMatterProfileSoliton nodeOperator.
# Yu Zhao (17-September-2025)
pathParameters      = os.path.abspath(f"parameters/darkMatterProfilesChan2022.xml"               )
pathOutputDirectory = os.path.abspath(f"outputs/test-soliton-coreMass"                           )
pathOutputLog       = os.path.abspath(f"outputs/test-soliton-coreMass/test-soliton-coreMass.log" )
pathOutputModel     = os.path.abspath(f"outputs/test-soliton-coreMass/test-soliton-coreMass.hdf5")
status              = subprocess.run(f"mkdir -p {pathOutputDirectory}",shell=True)

# Run the model and check for completion 
print("Running model...")
log    = open(pathOutputLog,"w")
status = subprocess.run(f"cd ..; ./Galacticus.exe {pathParameters}",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run(f"cat {pathOutputLog}",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run(f"grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" {pathOutputLog}",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run(f"cat {pathOutputLog}",shell=True)
    sys.exit()
print("SUCCESS: model run")

# Best-fitting parameters from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
alpha = 0.515            
beta  = 8.0e6
gamma = 10.0**(-5.73)

# Test if the core mass evolution agrees with analytic expectations.
with h5py.File(pathOutputModel, 'r') as file:
    # Extract the fuzzy dark matter particle mass (and convert to units of eV).
    massParticle = file['Parameters/darkMatterParticle'].attrs['mass']*1.0e-22

    # Get z=0 virial density contrast.
    output           = file["/Outputs/Output6/nodeData"]
    densityContrast0 = output["densityContrastVirial"][:]
    
    # Iterate over outputs.
    for i in range(1, 7):  # Output1 ~ Output6
        output           = file[f"/Outputs/Output{i}/nodeData"]
        massCore         = output["solitonMassCoreNormal"][:]
        densityContrastZ = output["densityContrastVirial"][:]
        redshift         = output["redshift"             ][:]
        massHalo         = output["basicMass"            ][:]
        expansionFactor  = 1.0/(1.0+redshift)
        
        # Evaluate the analytic expectation for the core mass using equation (15) of Chan et al. (2022; MNRAS; 551; 943;
        # https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
        massCoreAnalytic = (
            +beta
            *(massParticle/8.0e-23)**(-1.5)
            +(np.sqrt(densityContrastZ/densityContrast0[0])*massHalo/gamma)**alpha
            *(massParticle/8.0e-23)**(1.5*(alpha-1.0))
        )/np.sqrt(expansionFactor)

        # Verify that numerical results agree with the analytic expectation.
        if np.allclose(massCore, massCoreAnalytic, rtol=massCore*1e-4):
            print(f"SUCCESS [output {i}]: results do agree"   )
        else:
            print(f"FAILED  [output {i}]: results do not agree")
