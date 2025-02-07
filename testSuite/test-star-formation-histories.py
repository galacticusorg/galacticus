#!/usr/bin/env python3
import subprocess
import sys
import h5py
import re
import os
import numpy as np
import xml.etree.ElementTree as ET
import warnings

# Check star formation history calculations. A closed box model - for which the evolution of
# star formation rate and gas metallicity are known analytically - is used. Results from
# Galacticus are compared to analytic solutions.
# Andrew Benson (31-August-2024)

# Iterate over star formation history types.
for sfhType in ( "adaptive", "metallicitySplit" ):
    # Construct the file name suffix.
    suffix = sfhType[0].upper() + sfhType[1:]
    print("Testing '"+sfhType+"' star formation histories...")
    
    # Run the model and check for completion.
    print("   Running model...")
    status = subprocess.run("mkdir -p outputs",shell=True)
    log = open("outputs/test-star-formation-histories.log","w")
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/starFormationHistory"+suffix+".xml",stdout=log,stderr=log,shell=True)
    log.close()
    if status.returncode != 0:
        print("   ...done ("+str(status)+")")
        print("   FAILED: model run:")
        subprocess.run("cat outputs/test-star-formation-histories.log",shell=True)
        sys.exit()
    else:
        print("   ...done")
        print("   Checking for errors...")
        status = subprocess.run("grep -q -i -e fatal -e aborted -e \"task failed\" -e \"Galacticus experienced an error in the GSL library\" outputs/test-star-formation-histories.log",shell=True)
    if status.returncode == 0:
        print("   ...done ("+str(status)+")")
        print("   FAILED: model run (errors):")
        subprocess.run("cat outputs/test-star-formation-histories.log",shell=True)
        sys.exit()
    else:
        print("   ...done")
        print("   SUCCESS: model run")

    # Find Solar metallicity from the Galacticus source file.
    metallicitySolar = None
    for line in open(os.environ['GALACTICUS_EXEC_PATH']+'/source/numerical.constants.astronomical.F90'):
        match = re.match(r'.*variable="metallicitySolar"\s+value="([0-9\.\+\-]+)d([0-9\+\-]+)"',line)
        if match:
            metallicitySolar = float(match.group(1))*10**int(match.group(2))
    if metallicitySolar is None:
        print("   FAILED: unable to find Solar metallicity")
        sys.exit()

    # Extract the initial propeties of the galaxy from the merger tree file.
    tree         = ET.parse('parameters/starFormationHistoryTree.xml')
    timeStart    = float(tree.find('.//basic/time'  ).text)
    massGasStart = float(tree.find('.//disk/massGas').text)
    
    # Open the model and extract the recycled fraction, metal yield, and star formation timescale.
    model                  = h5py.File('outputs/starFormationHistory'+suffix+'.hdf5','r')
    outputs                = model['Outputs'                                                 ]
    stellarPopulation      = model['Parameters/stellarPopulation'                            ]
    starFormation          = model['Parameters/starFormationRateDisks/starFormationTimescale']
    recycledFraction       = stellarPopulation.attrs['recycledFraction']
    metalYield             = stellarPopulation.attrs['metalYield'      ]
    timescaleStarFormation = starFormation    .attrs['timescale'       ]

    # Compute the effective timescale for star formation, accounting for recycling.
    timescaleStarFormationEffective = timescaleStarFormation/(1.0-recycledFraction)

    # Read the star formation history.
    nodes                    = model['Outputs/Output1/nodeData'    ]
    diskStarFormationHistory = nodes['diskStarFormationHistoryMass']
    massFormed               = np.vstack(diskStarFormationHistory[:][0])
    metallicity              = diskStarFormationHistory.attrs['metallicity']
    # Handle the different ways in which the time bin data can be stored.
    if 'time' in diskStarFormationHistory.attrs:
        time = diskStarFormationHistory.attrs['time']
    else:
        time = nodes['diskStarFormationHistoryTimes'][0]
        
    # Find the minimum and maximum time associated with star formation in each bin. 
    warnings.filterwarnings('ignore')
    ## First, find the time at which each metallicity bin boundary is reached in the closed box
    ## chemical evolution model.
    timeMetallicityReached = timescaleStarFormation*metallicity/(metalYield/metallicitySolar)+timeStart
    ## Find the minimum and maximum times associated with the time bin boundaries.
    timeMinimumTime        = np.tile(np.concatenate((np.zeros(1),time[0:-1])),(len(metallicity),1))
    timeMaximumTime        = np.tile(                            time        ,(len(metallicity),1))
    ## Find the minimum and maximum times associated with the metallicity bin boundaries.
    timeMinimumMetallicity = np.transpose(np.tile(np.concatenate((np.zeros(1),timeMetallicityReached[0:-1])),(len(time),1)))
    timeMaximumMetallicity = np.transpose(np.tile(                            timeMetallicityReached        ,(len(time),1)))
    ## Find the limiting factor in setting the minimum and maximum time over which star
    ## formation is accumulated to each bin. Also, do not let the maximum time be les than the
    ## minimum time, or the minimum time be before the start time.
    timeMinimum            = np.maximum(np.maximum(timeMinimumTime,timeMinimumMetallicity),timeStart*np.ones(timeMinimumTime.shape))
    timeMaximum            = np.maximum(np.minimum(timeMaximumTime,timeMaximumMetallicity),                  timeMinimum           )

    # Compute the mass of stars that we expect to have formed in each bin from the closed box evolution model.
    massFormedTarget = (
        +massGasStart/(1.0-recycledFraction)
        *(
            +np.exp(-(timeMinimum-timeStart)/timescaleStarFormationEffective)
            -np.exp(-(timeMaximum-timeStart)/timescaleStarFormationEffective)
        )
    )

    # Report on status.
    status = np.allclose(massFormed,massFormedTarget,rtol=0.01,atol=1.0)
    if status:
        print("   SUCCESS: mass of stars formed in all bins")
    else:
        print("   FAILED: mass of stars formed in all bins")
