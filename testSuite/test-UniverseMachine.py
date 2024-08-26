#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np
from numpy.testing import assert_allclose

# Check internal self-consistency of adaptive star formation histories.
# Andrew Benson (25-March-2021)
def um_smhm_scaling_param_def():
    # Scaling variables from universe machine table j1
    # https://arxiv.org/pdf/1806.07893     
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
    ensurenotnone = lambda y: 0.0 if y is None else y
    return {
                prefix + "0"   : ensurenotnone(scaling_parameters.get(vstr + "0"  )),
                prefix + "a"   : ensurenotnone(scaling_parameters.get(vstr + "a"  )),
                prefix + "lna" : ensurenotnone(scaling_parameters.get(vstr + "lna")),
                prefix + "z"   : ensurenotnone(scaling_parameters.get(vstr + "z"  )),
    }

def um_smhm_scaling(z, y0, ya, ylna, yz):
    # Scaling function from universe machine apprendix equations j3 - j8
    # https://arxiv.org/pdf/1806.07893     
    a = 1 / (1 + z)
   
    return  (
             +y0 
             +ya
             *(
                 +a
                 -1
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

    m1       = np.power(10,logm1   )
    gamma    = np.power(10,loggamma)


    x = np.log10(m / m1)

    log_msr = (
               +epsilon
               -np.log10(
                          +np.power(
                                     +10,
                                     -alpha 
                                     *x
                                    )
                          +np.power(
                                    +10,
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
                    +10,
                    log_msr
                   )
            *m1
           ) 




# Run the model and check for completion 
print("Running model...")
status = subprocess.run("mkdir -p outputs/test-UniverseMachine",shell=True)
log = open("outputs/test-UniverseMachine/galacticus.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/test-UniverseMachine.xml",stdout=log,stderr=log,shell=True)
log.close()
print("...done ("+str(status)+")")
if status.returncode != 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-UniverseMachine/galacticus.log",shell=True)
    sys.exit()
print("Checking for errors...")
status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-UniverseMachine/galacticus.log",shell=True)
print("...done ("+str(status)+")")
if status.returncode == 0:
    print("FAILED: model run (errors):")
    subprocess.run("cat outputs/test-UniverseMachine/galacticus.log",shell=True)
    sys.exit()
print("SUCCESS: model run")

# Open the model and extract the recycled fraction.
model             = h5py.File('outputs/test-UniverseMachine/galacticus.hdf5','r')
massHalo          = model["Outputs/Output1/nodeData/massHaloEnclosedCurrent"][:]
massStellar       = model["Outputs/Output1/nodeData/spheroidMassStellar"    ][:]
isIsolated        = model["Outputs/Output1/nodeData/nodeIsIsolated"         ][:].astype(bool)
redshift          = model["Outputs/Output1/nodeData/redshiftLastIsolated"   ][:]
        


massHost, redshift = massHalo[isIsolated], 0.0

massStellarPython = um_smhm(massHost, redshift)
assert_allclose(massStellar[isIsolated], massStellarPython)

