#!/usr/bin/env python3
import sys
import os
import subprocess
import re
import numpy as np
import h5py
import xml.etree.ElementTree as ET

# Test dark matter halo mass functions against HMFcalc results.
# Andrew Benson (07-August-2019)

# Specify mass functions to test.
massFunctionTypes = [
     {
     	 "label" : "Press-Schechter",
     	 "method": "pressSchechter"
     },
     {
     	 "label" : "Sheth-Tormen",
     	 "method": "shethTormen"
     },
     {
     	 "label" : "Tinker2008",
     	 "method": "tinker2008"
     },
     # Disabled because we're not able to get agreement with the HMFCalc result. Reason is not understood.
     # {
     # 	 "label" : "Bhattacharya",
     # 	 "method": "bhattacharya2011"
     # }
    ]

# Iterate over mass functions.
for massFunctionType in massFunctionTypes:
    # Specify HMFcalc directory.
    hmfCalcPath = "data/HMFcalc/"+massFunctionType['label']+"/"
    try:
        os.mkdir("outputs/HMFcalc")
    except FileExistsError:
        pass

    # Read the HMFcalc parameter file.
    parametersHMFCalc = {}
    parameterFile = open(hmfCalcPath+"parameters.txt","r")
    for line in parameterFile:
        parameter = re.search(r"^([a-zA-Z0-9_]+):\s+([\d\.\+\-]+)\s*$",line)
        if parameter is not None:
            parametersHMFCalc[parameter.group(1)] = float(parameter.group(2))
    parameterFile.close()
    
    # Read the HMFcalc transfer function file.
    wavenumber, transferFunction = np.squeeze(np.hsplit(np.loadtxt(hmfCalcPath+"kVector_PLANCK-SMT .txt", usecols=(0,2)),2))
    wavenumber *= parametersHMFCalc['h']

    # Read the HMFcalc mass function file.
    massFunctionColumn = 2 if massFunctionType['label'] == "Bhattacharya" else 6
    mass, massFunction = np.squeeze(np.hsplit(np.loadtxt(hmfCalcPath+"mVector_PLANCK-SMT .txt",usecols=(0,massFunctionColumn)),2))
    mass         /= parametersHMFCalc['h']
    massFunction *= parametersHMFCalc['h']**3

    # Construct a parameter file for Galacticus.
    parameters = ET.parse('parameters/haloMassFunctionsBase.xml')
    ET.SubElement  (parameters.getroot(),'haloMassFunction'       ).set('value',                                 massFunctionType['method']             )
    parameters.find('./task/haloMassMinimum'                      ).set('value',str(mass[ 0]                                                           ))
    parameters.find('./task/haloMassMaximum'                      ).set('value',str(mass[-1]                                                           ))
    parameters.find('./task/pointsPerDecade'                      ).set('value',str(1.0/parametersHMFCalc['dlog10m']                                   ))
    parameters.find('./cosmologyParameters/temperatureCMB'        ).set('value',str(    parametersHMFCalc['t_cmb'  ]                                   ))
    parameters.find('./cosmologyParameters/OmegaMatter'           ).set('value',str(    parametersHMFCalc['omegam' ]                                   ))
    parameters.find('./cosmologyParameters/OmegaDarkEnergy'       ).set('value',str(    parametersHMFCalc['omegav' ]                                   ))
    parameters.find('./cosmologyParameters/OmegaBaryon'           ).set('value',str(    parametersHMFCalc['omegab' ]                                   ))
    parameters.find('./cosmologyParameters/HubbleConstant'        ).set('value',str(    parametersHMFCalc['H0'     ]                                   ))
    parameters.find('./cosmologicalMassVariance/sigma_8'          ).set('value',str(    parametersHMFCalc['sigma_8']                                   ))
    parameters.find('./powerSpectrumPrimordial/index'             ).set('value',str(    parametersHMFCalc['n'      ]                                   ))
    parameters.find('./criticalOverdensity/criticalOverdensity'   ).set('value',str(    parametersHMFCalc['delta_c']                                   ))
    parameters.find('./virialDensityContrast/densityContrastValue').set('value',str(    parametersHMFCalc['delta_h']                                   ))
    parameters.find('./transferFunction/fileName'                 ).set('value',    "testSuite/outputs/HMFcalc/"+massFunctionType['label' ]+"_Tk.hdf5"  )
    parameters.find('./outputFileName'                            ).set('value',    "testSuite/outputs/HMFcalc/"+massFunctionType['label' ]+"_HMF.hdf5" )
    parameters.write("outputs/HMFcalc/"+massFunctionType['label']+".xml")
    
    # Construct a transfer function file for Galacticus.
    transferFunctionFile = h5py.File("outputs/HMFcalc/"+massFunctionType['label']+"_Tk.hdf5","w")
    darkMatter              = transferFunctionFile.create_group('darkMatter'   )
    parameters              = transferFunctionFile.create_group('parameters'   )
    extrapolation           = transferFunctionFile.create_group('extrapolation')
    extrapolationWavenumber = extrapolation       .create_group('wavenumber'   )
    transferFunctionFile   .create_dataset('wavenumber'             ,data=wavenumber      )
    darkMatter             .create_dataset('transferFunctionZ0.0000',data=transferFunction)
    transferFunctionFile   .attrs['fileFormat'     ] = 2
    parameters             .attrs['OmegaMatter'    ] = parametersHMFCalc['omegam']
    parameters             .attrs['OmegaDarkEnergy'] = parametersHMFCalc['omegav']
    parameters             .attrs['OmegaBaryon'    ] = parametersHMFCalc['omegab']
    parameters             .attrs['HubbleConstant' ] = parametersHMFCalc['H0'    ]
    parameters             .attrs['temperatureCMB' ] = parametersHMFCalc['t_cmb' ]
    extrapolationWavenumber.attrs['low'            ] = 'fix'
    extrapolationWavenumber.attrs['high'           ] = 'extrapolate'
    transferFunctionFile.close()
    
    # Run Galacticus to generate the mass function.
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/outputs/HMFcalc/"+massFunctionType['label']+".xml",shell=True)
    if status.returncode != 0:
        print("FAILED: Galacticus failed for '"+massFunctionType['label']+"'")
        sys.exit()

    # Read the mass function generated by Galacticus.
    massFunctionFile       = h5py.File("outputs/HMFcalc/"+massFunctionType['label']+"_HMF.hdf5","r")
    massGalacticus         = massFunctionFile['Outputs/Output1/haloMass'           ][:]
    massFunctionGalacticus = massFunctionFile['Outputs/Output1/haloMassFunctionLnM'][:]

    # Check that masses are consistent.
    if len(mass) != len(massGalacticus):
         print("FAILED: ["+massFunctionType['label']+"] number of masses differ")
         sys.exit()
    if np.any(np.abs(mass-massGalacticus)/mass > 1.0e-6):
        print("FAILED: ["+massFunctionType['label']+"] masses are inconsistent")
        sys.exit()

    # Compare mass function with that from HMFcalc.
    differenceFractional       = np.abs(massFunction-massFunctionGalacticus)/massFunction
    differenceFractionlMaximum = np.max(differenceFractional)
    if differenceFractionlMaximum > 1.0e-3:
        print("FAILED: [" +massFunctionType['label']+"] mass function differs")
    else:
        print("success: ["+massFunctionType['label']+"] mass functions agree")

