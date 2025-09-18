#!/usr/bin/env python3
import subprocess
import os
import sys
import h5py
import numpy as np
import xml.etree.ElementTree as ET

# Test darkMatterProfileSoliton nodeOperator.
# Yu Zhao (17-September-2025)
path_in_param = os.path.abspath(f"parameters/darkMatterProfilesChan2022.xml"                )
path_out_dir  = os.path.abspath(f"outputs/test-soliton-coreMass"                            )
path_out_log  = os.path.abspath(f"outputs/test-soliton-coreMass/test-soliton-coreMass.log"  )
path_out_hdf  = os.path.abspath(f"outputs/test-soliton-coreMass/test-soliton-coreMass.hdf5" )

status = subprocess.run(f"mkdir -p {path_out_dir}",shell=True)

# Run the model and check for completion 
print("Running model...")

log = open(path_out_log,"w")
status = subprocess.run(f"cd ..; ./Galacticus.exe {path_in_param}",stdout=log,stderr=log,shell=True)
log.close()

print("...done ("+str(status)+")")

if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run(f"cat {path_out_log}",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run(f"grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" {path_out_log}",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run(f"cat {path_out_log}",shell=True)
    sys.exit()
print("SUCCESS: model run")

# Best-fitting parameters from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
alpha = 0.515            
beta  = 8.0e6
gamma = 10.0**(-5.73)
 
with h5py.File(path_out_hdf, 'r') as f:
    for i in range(1, 7):  # Output1 ~ Output6
        g = f[f"/Outputs/Output{i}/nodeData"]

        massCore        = g["solitonMassCore"][:]
        zeta_0          = g["solitonZeta0"][:]
        zeta_z          = g["solitonZetaZ"][:]
        expansionFactor = g["solitonExpansionFactor"][:]
        massParticle    = g["solitonMassParticle"][:]
        massHalo        = g["solitonMassHalo"][:]

        # Equation (15) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
        massCoreAnalytic = (beta*(massParticle/8.0e-23)**(-1.5)
                           +(np.sqrt(zeta_z/zeta_0)*massHalo/gamma)**alpha
                           * (massParticle/8.0e-23)**(1.5*(alpha-1.0))
                           )/np.sqrt(expansionFactor)

        if np.allclose(massCore, massCoreAnalytic, rtol=massCore*1e-4):
            print(f"SUCCESS: results do agree"   )
        else:
            print(f"FAILED: results do not agree")
