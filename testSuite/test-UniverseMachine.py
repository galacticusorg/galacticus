#!/usr/bin/env python3
import subprocess
import os
import sys
import h5py
import numpy as np
import xml.etree.ElementTree as ET

# Check UniverseMachine results.
# Charles Gannon (26-August-2024)

def um_smhm_scaling_param_def():
    # Scaling variables from UniverseMachine table J1 (https://arxiv.org/pdf/1806.07893).
    return dict(
                 e0     = -1.435, 
                 ea     = +1.831,
                 elna   = +1.368,        
                 ez     = -0.217,
                 m0     = +12.035,
                 ma     = +4.556,
                 mlna   = +4.417, 
                 mz     = -0.731, 
                 a0     = +1.963, 
                 aa     = -2.316,
                 alna   = -1.732,
                 az     = +0.178,
                 b0     = +0.482,
                 ba     = -0.841,
                 blna   = +0.000,
                 bz     = -0.471,
                 d0     = +0.411,                
                 da     = +0.000,
                 dlna   = +0.000,
                 dz     = +0.000,
                 g0     = -1.034,
                 ga     = -3.100,
                 glna   = +0.000,
                 gz     = -1.055, 
                )                

def um_smhm_scaling_fromdict(scaling_parameters:dict[str, float], vstr:str, prefix = "y"):
    return {
                prefix + "0"   : scaling_parameters.get(vstr + "0"  ),
                prefix + "a"   : scaling_parameters.get(vstr + "a"  ),
                prefix + "lna" : scaling_parameters.get(vstr + "lna"),
                prefix + "z"   : scaling_parameters.get(vstr + "z"  ),
    }

def um_smhm_scaling(z, y0, ya, ylna, yz):
    # Scaling function from universe machine apprendix equations j3 - j8
    # https://arxiv.org/pdf/1806.07893     
    a = 1.0 / (1.0 + z)
   
    return  (
             +y0 
             +ya
             *(
                 +a
                 -1.0
              ) 
             -ylna
             *np.log(a) 
             +yz
             *z
             ) 

def um_smhm(m, z, umparams = None):
    # Variable M_\star from universe machine, appendix j1
    # https://arxiv.org/pdf/1806.07893     
    
    if umparams is None:
        umparams = um_smhm_scaling_param_def()

    logm1    = um_smhm_scaling(z, **um_smhm_scaling_fromdict(umparams, "m")) 
    loggamma = um_smhm_scaling(z, **um_smhm_scaling_fromdict(umparams, "g")) 

    epsilon  = um_smhm_scaling(z, **um_smhm_scaling_fromdict(umparams, "e")) 
    alpha    = um_smhm_scaling(z, **um_smhm_scaling_fromdict(umparams, "a")) 
    beta     = um_smhm_scaling(z, **um_smhm_scaling_fromdict(umparams, "b")) 
    delta    = um_smhm_scaling(z, **um_smhm_scaling_fromdict(umparams, "d")) 

    m1       = np.power(10.0,logm1   )
    gamma    = np.power(10.0,loggamma)


    x = np.log10(m / m1)

    log_msr = (
               +epsilon
               -np.log10(
                          +np.power(
                                     +10.0,
                                     -alpha 
                                     *x
                                    )
                          +np.power(
                                    +10.0,
                                    -beta
                                    *x
                                   )
                        )
               +gamma
               *np.exp  (
                          -0.5
                          *np.power(
                                     +x
                                     /delta
                                     ,2
                                    ) 
                        )
               )
              
    return (
            +np.power(
                      +10.0,
                      log_msr
                     )
            *m1
           ) 


z_space = np.linspace(0,10,21)
status = subprocess.run("mkdir -p outputs/test-UniverseMachine outputs/test-UniverseMachine/parameters outputs/test-UniverseMachine/logs",shell=True)

for z in z_space:
    # Run the model and check for completion 
    print("Running model...")
    path_param = os.path.abspath(f"outputs/test-UniverseMachine/parameters/test-UniverseMachine-z{z:.1f}.xml" )
    path_log   = os.path.abspath(f"outputs/test-UniverseMachine/logs/test-UniverseMachine-z{z:.1f}.log"       )
    path_out   = os.path.abspath(f"outputs/test-UniverseMachine/test-UniverseMachine-z{z:.1f}.hdf5"           )
    
    log = open(path_log,"w")
    
    # Modify the output redshift of parameters/test-UniverseMachine.xml
    xroot = ET.parse("parameters/test-UniverseMachine.xml")

    element_zbase = xroot.find("mergerTreeConstructor/redshiftBase")
    element_zout  = xroot.find("outputTimes/redshifts"             )
    element_fname = xroot.find("outputFileName"                    )
    
    element_zbase.set("value","{:.2f}".format(z))
    element_zout.set("value","{:.2f}".format(z))
    element_fname.set("value",path_out)   

    xroot.write(path_param,xml_declaration=True, encoding="UTF-8")

    status = subprocess.run(f"cd ..; ./Galacticus.exe {path_param}",stdout=log,stderr=log,shell=True)
    log.close()
    print("...done ("+str(status)+")")
    if status.returncode != 0:
        print("FAILED: model run:")
        subprocess.run(f"cat {path_log}",shell=True)
        sys.exit()
    print("Checking for errors...")
    status = subprocess.run(f"grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" {path_log}",shell=True)
    print("...done ("+str(status)+")")
    if status.returncode == 0:
        print("FAILED: model run (errors):")
        subprocess.run(f"cat {path_log}",shell=True)
        sys.exit()
    print("SUCCESS: model run")
    
    # Open the model and extract the recycled fraction.
    model       = h5py.File(path_out,'r')
    massHalo    = model["Outputs/Output1/nodeData/massHaloEnclosedCurrent"][:]
    massStellar = model["Outputs/Output1/nodeData/spheroidMassStellar"    ][:]
    isIsolated  = model["Outputs/Output1/nodeData/nodeIsIsolated"         ][:].astype(bool)
    redshift    = model["Outputs/Output1/nodeData/redshiftLastIsolated"   ][:]  
    
    massHost, redshift, massStellarHost = massHalo[isIsolated], redshift[isIsolated], massStellar[isIsolated]
    
    massStellarPython = um_smhm(massHost, redshift)
    
    select = massStellarPython / massHost > 1E-4
    
    if np.allclose(np.log10(massStellarHost[select]), np.log10(massStellarPython[select]), rtol=1e-2):
        print(f"SUCCESS: results do agree at redshift z={z:.1f}"   )
    else:
        print(f"FAILED: results do not agree z={z:.1f}")

