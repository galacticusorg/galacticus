#!/usr/bin/env python
import numpy as np
import h5py
import argparse
import sys
import os
import shutil
import re
import copy
import lxml.etree as ET
import queueManager
import cloudy
import subprocess
import time
import datetime
import urllib.request
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import json
import codecs
from termcolor import colored
from git import Repo

# Generate tables of Cloudy models for use in emission line calculations.
# Andrew Benson (21-April-2025)

# Replace digits with characters.
def replaceDigits(x):
    string = ""
    for i, c in enumerate(str(x)):
        string += chr(int(c)+65)
    return string

# Convert command line arguments to floats.
def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))
    return x

# Convert command line arguments to ints.
def restricted_int(x):
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a integer literal" % (x,))
    return x

def establishGridSSP(grid,args):
    # Establish the grid of models to use for SSP calculations.
    # Get an SSP model.
    ## Download data if necessary.
    match = re.match(r'^http.*\/([^\?]+).*',args.sspFileName)
    if match:
        fileName = match.group(1)
        if not os.path.isfile(args.workspace+fileName):
            try:
                urllib.request.urlretrieve(args.sspFileName, args.workspace+fileName)
            except urllib.request.HTTPError as e:
                print(f"Failed to download stellar populations file: HTTP Error: {e.code} - {e.reason}")
            except Exception as e:
                print(f"Failed to download stellar populations file: {e}")
        grid['sspURL'] = args.sspFileName
        args.sspFileName = args.workspace+fileName
    ## Read data from file.
    stellarPopulations          = h5py.File(args.sspFileName,'r')
    ages                        = stellarPopulations['ages'         ][:]
    logMetallicities            = stellarPopulations['metallicities'][:]
    wavelength                  = stellarPopulations['wavelengths'  ][:]
    spectra                     = stellarPopulations['spectra'      ][:]
    ## Determine energies in Rydbergs.
    energy                      = plancksConstant*speedOfLight/wavelength/angstroms/electronVolt/rydbergEnergy
    ## Construct wavelength intervals.
    deltaWavelength             = copy.deepcopy(wavelength)
    deltaWavelength[0:-2]       = wavelength[1:-1]-deltaWavelength[0:-2]
    deltaWavelength[  -1]       =                  deltaWavelength[  -2]
    # Refine the range of metallicities, by interpolating the spectra to intermediate points.
    refineMetallicityBy       = 2
    countRefinedMetallicities = refineMetallicityBy*(len(logMetallicities)-1)+1
    logMetallicitiesRefined   = np.zeros([countRefinedMetallicities                                  ])
    spectraRefined            = np.zeros([countRefinedMetallicities,spectra.shape[1],spectra.shape[2]])
    for i in range(len(logMetallicities)-1):
        deltaLogMetallicity = logMetallicities[i+1]-logMetallicities[i]
        for j in range(refineMetallicityBy):
            zero    = (spectra[i,:,:] <= 0.0) | (spectra[i+1,:,:] <= 0.0)
            nonZero = (spectra[i,:,:] >  0.0) & (spectra[i+1,:,:] >  0.0)
            # Spectra are interpolated in log-log space where possible (i.e. where the spectrum is non-zero at both endpoints of the
            # interpolation), and in lin-log space otherwise. Metallicities are always interpolated in log space.
            spectraRefined         [(i*refineMetallicityBy+j),:,:][   zero] =        +       spectra     [(i  ),:,:][   zero] *(1.0-(j/refineMetallicityBy)) \
                                                                                     +       spectra     [(i+1),:,:][   zero] *(0.0+(j/refineMetallicityBy))
            spectraRefined         [(i*refineMetallicityBy+j),:,:][nonZero] = np.exp(
                                                                                     +np.log(spectra     [(i  ),:,:][nonZero])*(1.0-(j/refineMetallicityBy))
                                                                                     +np.log(spectra     [(i+1),:,:][nonZero])*(0.0+(j/refineMetallicityBy))
                                                                               	    )
            logMetallicitiesRefined[(i*refineMetallicityBy+j)    ]          =        +logMetallicities   [(i  )    ]                                         \
 		                                                                     +deltaLogMetallicity                     *     (j/refineMetallicityBy)
    logMetallicitiesRefined[-1    ] = logMetallicities[-1    ]
    spectraRefined         [-1,:,:] = spectra         [-1,:,:]
    grid['ages'            ] = ages
    grid['logMetallicities'] = logMetallicitiesRefined
    grid['spectra'         ] = spectraRefined
    grid['wavelength'      ] = wavelength
    grid['energy'          ] = energy

    # Evaluate the number of Lyman continuum photons emitted per second for each population (age,metallicity).
    ## Construct the integrand. Spectra are in units of L☉ Hz¯¹. We want to evaluate Qₕ = ∫_Eₕ^∞ dν S_ν/hν = ∫₀^λₕ dλ S_ν / hλ.
    integrand                         = grid['spectra']*deltaWavelength/wavelength*luminositySolar/plancksConstant
    ## Find the range to include in each integral.
    hydrogenContinuum                 = wavelength < wavelengthLymanContinuum
    ## Evaluate the integrals.
    grid['ionizingLuminosityPerMass'] = np.sum(integrand[:,:,hydrogenContinuum],axis=2)

    # Determine the normalization method to use.
    if args.normalization == "ionizingLuminosity":
        # Define ionizing luminosities, Qₕ, to tabulate.
        grid['logHydrogenLuminosities'] = np.array([ 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0            ])
    elif args.normalization == "massStellar":
        # Define ionizing luminosities, Qₕ, to tabulate.
        grid['logStellarMasses'       ] = np.array([  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0, 5.5,  6.0 ])
    else:
        sys.exit("Expected `normalization` of `ionizingLuminosity` or `massStellar`")
    # Define hydrogen densities, nₕ, to tabulate.
    grid['logHydrogenDensities'   ] = np.array([ 1.0, 1.5,  2.0,  2.5,  3.0, 3.5, 4.0 ])
    # Specify the iterables in the grid.
    if args.normalization == "ionizingLuminosity":
        grid['iterables'] = ( "ages", "logMetallicities", "logHydrogenLuminosities"   , "logHydrogenDensities" )
        grid['names'    ] = ( "age" , "metallicity"     , "ionizingLuminosityHydrogen", "densityHydrogen"      )
    elif args.normalization == "massStellar":
        grid['iterables'] = ( "ages", "logMetallicities", "logStellarMasses"          , "logHydrogenDensities" )
        grid['names'    ] = ( "age" , "metallicity"     , "massStellar"               , "densityHydrogen"      )
    else:
        sys.exit("Expected `normalization` of `ionizingLuminosity` or `massStellar`")

def establishGridAGN(grid,args):
    # Establish the grid of models to use for AGN calculations.

    # Validate the normalization option.
    if args.normalization != "ionizingLuminosity":
        sys.exit("AGN models require ``--normalization ionizingLuminosity`")

    # Define spectral indices, α, to tabulate.
    grid['spectralIndices'        ] = np.array([ -1.2, -1.4, -1.7, -2.0 ])

    # Define ionization parameters, Uₛ, to tabulate.
    grid['logIonizationParameters'] = np.array([ -4.0, -3.0, -2.0, -1.0 ])
    
    # Define metallicities, Z, to tabulate.
    grid['logMetallicities'       ] = np.array([ -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5 ])
    
    # Define hydrogen densities, nₕ, to tabulate.
    grid['logHydrogenDensities'   ] = np.array([ 1.0,  1.5,  2.0,  2.5,  3.0, 3.5, 4.0 ])
    
    # Specify the iterables in the grid.
    grid['iterables'] = ( "spectralIndices", "logMetallicities", "logIonizationParameters", "logHydrogenDensities" )
    grid['names'    ] = ( "spectralIndex"  , "metallicity"     , "ionizationParameter"    , "densityHydrogen"      )

    # Construct spectra, and their bolometric luminosity normalization factors.
    ## Our spectrum is (Feltre, Charlot & Gutkin; 2016; MNRAS; 456; 3354; https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.3354F):
    ##
    ##          ⎧ A₀ ν⁻ⁿ for 0.001 ≤ λ/μm ≤ 0.25,
    ##  Sᵥ(λ) = ⎨ A₀ A₁ ν⁻⁰·⁵ for 0.25 < λ/μm ≤ 10.0,
    ##          ⎩ A₀ A₁ A₂ ν² for λ/μm > 10.0,
    ##
    ## where A₁ = ν₁ⁿ⁺⁰·⁵, and A₂ = ν₂⁻²·⁵, where ν₂ and ν₃ are the frequencies corresponding to 0.25 and 10μm respectively. This
    ## is integrated to find the bolometric luminosity in terms of A₀, and, from that, the normalization factor for a given
    ## bolometric lumnosity.
    ##
    ## First find the frequencies (in Rydberg units) corresponding to the break points.
    wavelengths   = np.array([ 10.0, 0.25, 0.001 ]) # in μm.
    frequencies   = plancksConstant*speedOfLight/micron/wavelengths/electronVolt/rydbergEnergy
    
    ## Define a grid of frequencies (in Rydberg units) at which to evaluate the spectrum.
    frequencyLow        = 1.0e-8; # Recommended low-energy limit from Hazy documentation.
    frequencyHigh       = frequencies[2]
    countPerDecade      = 100.0
    countFrequencies    = int(np.log10(frequencyHigh/frequencyLow)*countPerDecade)+2
    grid['frequencies'] = np.logspace(np.log10(frequencyLow),np.log10(frequencyHigh),countFrequencies)
    grid['spectra'    ] = []
    normalization       = np.array([])
    # Integrate for each spectral index.
    for i in range(len(grid['spectralIndices'])):
        # Find the continuity factors.
        A1            = frequencies[1]**(grid['spectralIndices'][i]+0.5)
        A2            = frequencies[0]**(                          -2.5)
        # Evaluate the bolometric luminosity normalization factor.
        integral0     = frequencies[0]**(3.0                           )/ 3.0
        integral1     = frequencies[1]**(0.5                           )/ 0.5                            -frequencies[0]**(0.5                           )/ 0.5
        integral2     = frequencies[1]**(1.0+grid['spectralIndices'][i])/(1.0+grid['spectralIndices'][i])-frequencies[0]**(1.0+grid['spectralIndices'][i])/(1.0+grid['spectralIndices'][i])
        integral      = +integral0+integral1+integral2;
        normalization = np.append(normalization,1.0/integral)
        # Evaluate the spectrum.
        range0     =                                           (grid['frequencies'] < frequencies[0])
        range1     = (grid['frequencies'] >= frequencies[0]) & (grid['frequencies'] < frequencies[1])
        range2     = (grid['frequencies'] >= frequencies[1])
        spectrum   = np.zeros(len(grid['frequencies']))
        spectrum[range0] = A1*A2*grid['frequencies'][range0]**(+2.0                       )
        spectrum[range1] = A1   *grid['frequencies'][range1]**(-0.5                       )
        spectrum[range2] =       grid['frequencies'][range2]**(+grid['spectralIndices'][i])
        grid['spectra'].append(spectrum)
    grid['normalization'] = normalization

def adjustAbundances(abundancesReference,metallicity,dustToMetalsRatio):
    # Copy and modify the provided reference abundances to create abundances suitable for the given metallicity (linear, relative
    # to Solar), and dust-to-metals ratio.
    # Make a copy.
    abundances = copy.deepcopy(abundancesReference)
    # Determine metallicity and dust-to-metals ratio for reference abundances.
    elements                 = sorted(abundancesReference,key=lambda x: abundancesReference[x]['atomicNumber'])
    isMetal                  = list(map(lambda x: abundances[x]['atomicNumber'] >  2,elements))
    isNotMetal               = list(map(lambda x: abundances[x]['atomicNumber'] <= 2,elements))
    abundancesByMass         = np.array(list(map(lambda x: abundancesReference[x]['atomicMass']*10.0**abundancesReference[x]['logAbundanceByNumber']                                                   ,elements)))
    dustByMass               = np.array(list(map(lambda x: abundancesReference[x]['atomicMass']*10.0**abundancesReference[x]['logAbundanceByNumber']*(1.0-abundancesReference[x]['undepletedFraction']),elements)))
    metallicityReference     = np.sum(abundancesByMass[isMetal])/np.sum(abundancesByMass         )
    dustMetalsRatioReference = np.sum(dustByMass      [isMetal])/np.sum(abundancesByMass[isMetal])
    # Adjust abundance of all metals with the overall metallicity.
    for element in elements:
        if abundances[element]['atomicNumber'] <= 2:
            continue
        abundances[element]['logAbundanceByNumber'] += np.log10(metallicity)
    # Apply any custom adjustments to individual elements.
    for element in elements:
        if 'adjustAbundance' in abundances[element]:
            abundances[element]['adjustAbundance'](abundances,metallicity*metallicityReference)
    # Renormalize to keep the total metallicity fixed.
    abundancesByMassNew   = np.array(list(map(lambda x: abundances[x]['atomicMass']*10.0**abundances[x]['logAbundanceByNumber'],elements)))
    renormalizationFactor = metallicity*metallicityReference/(1.0-metallicity*metallicityReference)*np.sum(abundancesByMassNew[isNotMetal])/np.sum(abundancesByMassNew[isMetal])
    for element in elements:
        # Skip non-metals.
        if abundances[element]['atomicNumber'] <= 2:
            continue
        # Renormalize.
        abundances[element]['logAbundanceByNumber'] += np.log10(renormalizationFactor)
    # Adjust depletion factors. Here we linearly interpolate between the reference depletion fractions, f☉, and dust-to-metals
    # ratio ξ☉, and the limiting cases of (f,ξ)=(0,0) and (f,ξ)=(1,1) as per Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; §2.3.2).
    for element in elements:
        depletionFractionReference = 1.0-abundancesReference[element]['undepletedFraction']
        if dustToMetalsRatio < dustMetalsRatioReference:
            depletionFraction = depletionFractionReference*dustToMetalsRatio/dustMetalsRatioReference
        else:
            depletionFraction = depletionFractionReference+(dustToMetalsRatio-dustMetalsRatioReference)/(1.0-dustMetalsRatioReference)*(1.0-depletionFractionReference)
        undepletedFraction = 1.0-depletionFraction
        abundances[element]['undepletedFraction'] = undepletedFraction
    # Compute the correction to the grain abundances (relative to Cloudy's default Orion grains) required to get our desired dust-to-metals ratio.
    dustToGasRatioCloudy         = 5.396e-03; # This is the dust-to-gas ratio for Cloudy's default Orion grains.
    abundancesByMassFinal        = np.array(list(map(lambda x: abundances[x]['atomicMass']*10.0**abundances[x]['logAbundanceByNumber'],elements)))
    i                            = -1
    for element in elements:
        i = i+1
        abundances[element]['abundanceByMass'] = abundancesByMassFinal[i]
    metallicityFinal             = np.sum(abundancesByMassFinal[isMetal])/np.sum(abundancesByMassFinal)
    dustToMetalsRatioCloudy      = dustToGasRatioCloudy/metallicityFinal
    dustToMetalsBoostLogarithmic = np.log10(dustToMetalsRatio/dustToMetalsRatioCloudy)
    # Return the modified set of abundances.
    return (dustToMetalsBoostLogarithmic,abundances)

def adjustAbundanceNitrogen(abundances,metallicity):
    # Adjust the abundance of nitrogen following the model of Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; eqn. 11).
    logOH                                   = abundances['O']['logAbundanceByNumber']
    OH                                      = 10.0**logOH
    NH                                      = 0.41*OH*(10.0**(-1.6)+10.0**(2.33+logOH))
    logNH                                   = np.log10(NH)
    abundances['N']['logAbundanceByNumber'] = logNH

def adjustAbundanceCarbon(abundances,metallicity):
    # Adjust the abundance of nitrogen following the model of Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G).
    logOH                                   = abundances['O']['logAbundanceByNumber']
    abundances['C']['logAbundanceByNumber'] = np.log10(0.44)+logOH
    
def adjustAbundanceHelium(abundances,metallicity):
    # Adjust the abundance of nitrogen following the model of Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; eqn. 12).
    massFractionHe                           = 0.2485+1.7756*metallicity
    massFractionH                            = 1.0-massFractionHe-metallicity
    HeH                                      = +(massFractionHe/abundances['He']['atomicMass'])/(massFractionH/abundances['H']['atomicMass'])
    logHeH                                   = np.log10(HeH)
    abundances['He']['logAbundanceByNumber'] = logHeH

def linesParse(job):
    # Parse output from a Cloudy job to extract line data.
    grid    = job['grid'   ]
    indices = job['indices']
    # Check for successful completion.
    label = " ".join(map(lambda x: str(x),indices))
    grid['lineData']['status'][indices] = 0
    checkDisaster = subprocess.run(['grep', '-q', 'DISASTER', job['logOutput']], capture_output=True, text=True)
    if checkDisaster.returncode == 0:
        print("FAIL (Cloudy failed disasterously): "+label+" see "+job['logOutput'])
        if grid['lineData']['status'][indices] == 0:
            grid['lineData']['status'][indices] = 1
        return
    checkFailed   = subprocess.run(['grep', '-q', 'FAILED'  , job['logOutput']], capture_output=True, text=True)
    if checkFailed  .returncode == 0:
        print("FAIL (Cloudy exited with non-zero status): "+label+" see "+job['logOutput'])
        if grid['lineData']['status'][indices] == 0:
            grid['lineData']['status'][indices] = 2
        return
    # Allow multiple attempts to read the lines file in case the file is written over NFS and we must wait for it to catch up.    
    attemptsMaximum = 10
    linesFound      = {}
    for attempt in range(attemptsMaximum):
 	# Read the lines file.
        badFile = True
        if os.path.isfile(job['linesFileName']):
            badFile = False
            linesFile = open(job['linesFileName'],"r")
            for line in linesFile:
                if re.match(r'^#',line):
                    continue
                columns = re.split(r'\t+',line)
                if len(columns) != 5:
                    badFile = True
 		# Extract line properties. Use the "intrinsic" line luminosities. From the Cloudy documentation:
 		#
 		#  | This is the spectrum produced within the cloud and does not include effects of dust that lies outside the
		#  | line-forming region. These intensities do not include the reddening effects of any grains or other opacity
 		#  | sources that lie outside the line-forming region.
                lineEnergy     = columns[0]
                lineLabel      = columns[1]
                lineLuminosity = columns[2]
                if lineLabel in lineList:
                    lineName       = lineList[lineLabel]
                    lineLuminosity = float(lineLuminosity)
                    lineEnergy     = float(lineEnergy    )
                    lineWavelength = plancksConstant*speedOfLight/rydbergEnergy/electronVolt/lineEnergy/angstroms
                    grid['lineData'][lineName]['luminosity'][indices] = 10.0**lineLuminosity
                    grid['lineData'][lineName]['wavelength']          = lineWavelength
                    linesFound[lineLabel] = 1
            linesFile.close()
        if badFile:
            if attempt == attemptsMaximum-1:
                print("FAIL [unable to find Cloudy output lines file]: "+label+" see "+job['logOutput'])
                if grid['lineData']['status'][indices] == 0:
                    grid['lineData']['status'][indices] = 3
            else:
                time.sleep(10)
        else:
            break
    # Check that we found all lines.
    for lineName in lineList:
        if lineName not in linesFound:
            print("FAIL [some emission lines missing]: "+label+" '"+lineName+"' see "+job['logOutput'])
            if grid['lineData']['status'][indices] == 0:
                grid['lineData']['status'][indices] = 3
    # Clean up.
    if grid['lineData']['status'][indices] == 0:
        if not args.noClean:
            os.remove(job['linesFileName'       ])
            os.remove(job['continuumFileName'   ])
            os.remove(job['launchFile'          ])
            os.remove(job['logOutput'           ])
            os.remove(job['cloudyScriptFileName'])

def reprocessSSP(grid,args):
    # Reprocess a single Cloudy job for SSP calculations.
    # Reprocess output files. This can be useful if some previous processing of Cloudy output files failed (we often have tens of
    # thousands of these so some intermittment failures can occur).
    tableFile                  = h5py.File(args.workspace+args.outputFileName,'r')
    lineGroup                  = tableFile['lines' ]
    grid['lineData']['status'] = np.transpose(lineGroup['status'][:])
    for lineIdentifier in lineList:
        lineName = lineList[lineIdentifier]
        grid['lineData'][lineName]['luminosity'] = np.transpose(lineGroup[lineName][:])
        grid['lineData'][lineName]['wavelength'] =              lineGroup[lineName].attrs['wavelength']
    # Determine the normalization method used.
    if args.normalization == "ionizingLuminosity":
        normalization = "logHydrogenLuminosities"
    elif args.normalization == "massStellar":
        normalization = "logStellarMasses"
    else:
        sys.exit("Expected `normalization` of `ionizingLuminosity` or `massStellar`")

    jobNumber = -1
    for iAge in range(grid['ionizingLuminosityPerMass'].shape[1]):
        for iMetallicity in range(grid['ionizingLuminosityPerMass'].shape[0]):
            for iNormalization in range(len(grid[normalization])):
                for iLogHydrogenDensity in range(len(grid['logHydrogenDensities'])):
                    jobNumber += 1
                    if grid['lineData']['status'][iAge,iMetallicity,iNormalization,iLogHydrogenDensity] == 0:
                        continue
 		    # Reset the status before attempting to reprocess.
                    statusOld = grid['lineData']['status'][iAge,iMetallicity,iNormalization,iLogHydrogenDensity]
                    linesParse(
                        {
                            "label":                               "emissionLines"+str(jobNumber)       ,
                            "launchFile":           args.workspace+"emissionLines"+str(jobNumber)+".sh" ,
                            "logOutput":            args.workspace+"emissionLines"+str(jobNumber)+".log",
                            "logError":             args.workspace+"emissionLines"+str(jobNumber)+".log",
                            "jobNumber":            jobNumber                                           ,
                            "grid":                 grid                                                ,
 	                    "cloudyScriptFileName": args.workspace+cloudyScriptFileName                 ,
 	                    "linesFileName":        args.workspace+"lines"        +str(jobNumber)+".out",
 	                    "continuumFileName":    args.workspace+"continuum"    +str(jobNumber)+".out",
 	                    "indices":              ( iAge, iMetallicity, iNormalization, iLogHydrogenDensity ),
                        }
		    )
                    print("Reprocess job number "+str(jobNumber)+" (status = "+str(statusOld)+" ==> "+str(grid['lineData']['status'][iAge,iMetallicity,iNormalization,iLogHydrogenDensity])+")")
    # Set data URL.
    match = re.match(r'^http.*\/([^\?]+).*',args.sspFileName)
    if match:
        grid['sspURL'] = args.sspFileName

def reprocessAGN(grid,args):
    # Reprocess a single Cloudy job for SSP calculations.
    # Reprocess output files. This can be useful if some previous processing of Cloudy output files failed (we often have tens of
    # thousands of these so some intermittment failures can occur).
    tableFile                  = h5py.File(args.workspace+args.outputFileName,'r')
    lineGroup                  = tableFile['lines' ]
    grid['lineData']['status'] = np.transpose(lineGroup['status'][:])
    for lineIdentifier in lineList:
        lineName = lineList[lineIdentifier]
        grid['lineData'][lineName]['luminosity'] = np.transpose(lineGroup[lineName][:])
        grid['lineData'][lineName]['wavelength'] =              lineGroup[lineName].attrs['wavelength']
    jobNumber = -1
    for iSpectralIndex in range(len(grid['spectralIndices'])):
        for iMetallicity in range(len(grid['logMetallicities'])):
            for iIonizationParameter in range(len(grid['logIonizationParameters'])):
                for iLogHydrogenDensity in range(len(grid['logHydrogenDensities'])):
                    jobNumber += 1
                    if grid['lineData']['status'][iSpectralIndex,iMetallicity,iIonizationParameter,iLogHydrogenDensity] == 0:
                        continue
 		    # Reset the status before attempting to reprocess.
                    statusOld = grid['lineData']['status'][iSpectralIndex,iMetallicity,iIonizationParameter,iLogHydrogenDensity]
                    linesParse(
                        {
                            "label":                               "emissionLines"+str(jobNumber)       ,
                            "launchFile":           args.workspace+"emissionLines"+str(jobNumber)+".sh" ,
                            "logOutput":            args.workspace+"emissionLines"+str(jobNumber)+".log",
                            "logError":             args.workspace+"emissionLines"+str(jobNumber)+".log",
                            "jobNumber":            jobNumber                                           ,
                            "grid":                 grid                                                ,
 	                    "cloudyScriptFileName": args.workspace+cloudyScriptFileName                 ,
 	                    "linesFileName":        args.workspace+"lines"        +str(jobNumber)+".out",
 	                    "continuumFileName":    args.workspace+"continuum"    +str(jobNumber)+".out",
 	                    "indices":              ( iSpectralIndex, iMetallicity, iIonizationParameter, iLogHydrogenDensity ),
                        }
		    )
                    print("Reprocess job number "+str(jobNumber)+" (status = "+str(statusOld)+" ==> "+str(grid['lineData']['status'][iSpectralIndex,iMetallicity,iIonizationParameter,iLogHydrogenDensity])+")")

def generateJobSSP(grid,args):
    # Generate a single Cloudy job for SSP calculations.
    # Extract the indices for this job.
    iAge                   = grid['counter'][0]
    iMetallicity           = grid['counter'][1]
    iNormalization         = grid['counter'][2]
    iLogHydrogenDensity    = grid['counter'][3]
    # If this is a rerun, load line data and status.
    if args.rerun:
        if not 'rerunStatusRead' in grid:
            tableFile                  = h5py.File(args.workspace+args.outputFileName,'r')
            lineGroup                  = tableFile['lines']
            grid['lineData']['status'] = np.transpose(lineGroup['status'][:])
            for lineIdentifier in lineList:
                lineName = lineList[lineIdentifier]
                grid['lineData'][lineName]['luminosity'] = np.transpose(lineGroup[lineName][:])
                grid['lineData'][lineName]['wavelength'] =              lineGroup[lineName].attrs['wavelength']
            grid['rerunStatusRead'] = 1
        statusOld = grid['lineData']['status'][iAge,iMetallicity,iNormalization,iLogHydrogenDensity]
        if statusOld == 0:
            return
    # Normalize the spectrum - this is a convenience only as the normalization will be recomputed by Cloudy.
    if not 'normalized' in grid:
        grid['normalized'] = np.zeros([len(grid['logMetallicities']),len(grid['ages'])],dtype=int)
    if grid['normalized'][iMetallicity,iAge] != 1: 
        normalizer = np.max(grid['spectra'][iMetallicity,iAge,:])
        grid['spectra'   ][iMetallicity,iAge,:] /= normalizer
        grid['normalized'][iMetallicity,iAge  ]  = 1
    # Find the logarithmic luminosity in ionizing photons for this model.
    if args.normalization == "ionizingLuminosity":
        logHydrogenLuminosity = grid['logHydrogenLuminosities'][iNormalization]
    elif args.normalization == "massStellar":
        logHydrogenLuminosity = grid['logStellarMasses'       ][iNormalization]+np.log10(grid['ionizingLuminosityPerMass'][iMetallicity,iAge])
    else:
        sys.exit("Expected `normalization` of `ionizingLuminosity` or `massStellar`")
    # Get abundances for this metallicity.
    ## Compute metallicity relative to Solar.
    metallicity                                    = 10.0**grid['logMetallicities'][iMetallicity]
    ## Specify a dust-to-metals ratio. The following corresponds to the dust-to-metals ratio for the reference model of
    # Gutkin, Charlot & Bruzual (2016; https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; table 1). It differs slightly
    # from their stated value, but agrees with our internal calculation of this value from their data.
    dustToMetals                                   = 0.401
    (dustToMetalsBoostLogarithmic, abundances)     = adjustAbundances(abundancesReference,metallicity,dustToMetals)
    # Compute the Strömgren radius for this model.
    radiusStromgrenLogarithmic = (logHydrogenLuminosity-np.log10(4.0*np.pi/3.0*coefficientRecombinationCaseB)-2.0*grid['logHydrogenDensities'][iLogHydrogenDensity])/3.0
    # Determine the inner radius of the cloud.
    radiusCloudInnerLogarithmic = radiusStromgrenLogarithmic+np.log10(args.factorMorphology)
    # If required, determine the cloud outer radius.
    if args.normalization == "massStellar" and args.stopOuterRadius:
        epsilonStar                 = 3.000000000e-02 # Star formation efficiency of clouds (Grudic et al.; 2019; MNRASm 488, 1501).
        massAtomic                  = 1.660539040e-24 # Atomic mass unit in g.
        massSolar                   = 1.989100000e+33 # Solar mass in g.
        elements                    = sorted(abundances,key=lambda x: abundancesReference[x]['atomicNumber'])
        fractionMassHydrogen        = abundances['H']['abundanceByMass']/np.sum(list(map(lambda x: abundances[x]['abundanceByMass'],elements)))
        massCloudGas                = (10.0**grid['logStellarMasses'    ][iNormalization]     )*massSolar                               /epsilonStar
        densityGas                  = (10.0**grid['logHydrogenDensities'][iLogHydrogenDensity])*massAtomic*abundances['H']['atomicMass']/fractionMassHydrogen
        radiusCloudOuter            = (3.0*massCloudGas/4.0/np.pi/densityGas+10.0**(3.0*radiusCloudInnerLogarithmic))**(1.0/3.0)
        radiusCloudOuterLogarithmic = np.log10(radiusCloudOuter)
    # Create a depletion file that can be read by Cloudy. Cloudy does not permit digits in these file names. To work around this,
    # we translate our numerical metallicity index into ASCII characters, by mapping digits (0→A, 1→B, etc.).
    if not args.noGrains:
        encodedMetallicity = replaceDigits(iMetallicity)
        depletionFileName  = "grains_"+encodedMetallicity+".dpl"
        if not 'depletions' in grid:
            grid['depletions'] = [0]*len(grid['logMetallicities'])
        if grid['depletions'][iMetallicity] != 1:
            depletionFile = open(args.workspace+depletionFileName,"w")
            for element in abundances:
                depletionFile.write(abundances[element]['name']+" "+str(abundances[element]['undepletedFraction'])+"\n")
            depletionFile.close()
            grid['depletions'][iMetallicity] = 1
    # Generate a Cloudy parameter file.
    cloudyScript      = "title emission line job number "+str(jobNumber)+"\n"
    cloudyScript     += "# ["+str(iAge               )+"] age     = "+str(grid['ages'                   ][iAge               ])+"\n"
    cloudyScript     += "# ["+str(iMetallicity       )+"] log Z   = "+str(grid['logMetallicities'       ][iMetallicity       ])+"\n"
    cloudyScript     += "# ["+str(iNormalization     )+"] log Q_H = "+str(      logHydrogenLuminosity                         )+"\n"
    if args.normalization == "massStellar":
        cloudyScript += "# ["+str(iNormalization     )+"] log M_* = "+str(grid['logStellarMasses'       ][iNormalization     ])+"\n"
    cloudyScript     += "# ["+str(iLogHydrogenDensity)+"] log n_H = "+str(grid['logHydrogenDensities'   ][iLogHydrogenDensity])+"\n"
    ## Set the input spectrum for Cloudy.
    if not 'cloudySpectrum' in grid:
        grid['cloudySpectrum'] = {}
    if not iAge             in grid['cloudySpectrum']:
        grid['cloudySpectrum'][iAge] = {}
    if not iMetallicity     in grid['cloudySpectrum'][iAge]:
        nonZero        = grid['spectra'][iMetallicity,iAge,:] > 0.0
        energy         =          grid['energy' ][                  :][nonZero]
        logSpectrum    = np.log10(grid['spectra'][iMetallicity,iAge,:][nonZero])
        counter        = -1
        cloudySpectrum = ""
        for iWavelength in reversed(range(len(energy))):
            counter += 1
            if counter % 3 == 0:
                if counter == 0:
                    cloudySpectrum +=  "interpolate"
                else:
                    cloudySpectrum +=  "\ncontinue"
            cloudySpectrum +=  " ("+str(energy[iWavelength])+" "+str(logSpectrum[iWavelength])+")"
        cloudySpectrum += "\n"
        grid['cloudySpectrum'][iAge][iMetallicity] = cloudySpectrum
    cloudyScript += grid['cloudySpectrum'][iAge][iMetallicity]
    ## Set normalization of the spectrum.
    cloudyScript += "q(h) = "+str(logHydrogenLuminosity)+"\n"
    # Set the chemical composition of the HII region.
    # Start with Cloudy's HII region abundances - but no grains - we will add these later.
    cloudyScript += "abundances HII region no grains\n"
    # Set individual element abundances.
    for element in abundances:
        # Skip hydrogen.
        if abundances[element]['atomicNumber'] < 2:
            continue
        # Set the abundance for this element.
        cloudyScript += "element abundances "+abundances[element]['name'].lower()+" "+str(abundances[element]['logAbundanceByNumber'])+"\n"
    # Specify Cloudy's default ISM grains, but with abundance reduced in proportion to the metallicity to retain a fixed
    # dust-to-metals ratio. Include the sublimation suppression function (see section 7.9.5 of Hazy1;
    # https://data.nublado.org/cloudy_releases/c23/c23.01.tar.gz).
    if not args.noGrains:
        cloudyScript += "grains Orion "+str(dustToMetalsBoostLogarithmic)+" _log function sublimation\n"
    # Deplete metals into grains using our custom depletions file. Note that these depletions differ from those assumed by
    # the "grains _ISM" model above. This seems to be at the ~10% level, so we do not worry too much about this.
    if not args.noGrains:
        cloudyScript += "metals deplete \""+depletionFileName+"\"\n"
    # Set HII region density - this is log₁₀(nₕ/cm¯³).
    cloudyScript += "hden "+str(grid['logHydrogenDensities'][iLogHydrogenDensity])+"\n"
    # Set other HII region properties.
    cloudyScript += "sphere expanding\n"
    if args.normalization == "massStellar" and args.stopOuterRadius:
        radiusCloudOuterLogarithmicLabel = f' {radiusCloudOuterLogarithmic:.3f}'
    else:
        radiusCloudOuterLogarithmicLabel = ''
    cloudyScript += f"radius {radiusCloudInnerLogarithmic:.3f}{radiusCloudOuterLogarithmicLabel}\n"
    # Set cosmic rays (needed to avoid problems in Cloudy in neutral gas).
    cloudyScript += "cosmic rays background\n"
    # Set stopping criteria.
    ## Using temperature as a stopping criterion can be problematic as the model can then extend far into the neutral region if
    ## grains are present and cause heating. Instead, we stop based when an electron fraction (default of 1%) is reached, or when
    ## an Lyman limit optical depth (default of 10) is reached which should accurately capture the ionization front.
    cloudyScript += "stop temperature off\n"
    cloudyScript += "stop efrac "                +str(         args.stopElectronFraction  )+"\n"
    cloudyScript += "stop Lyman optical depth = "+str(np.log10(args.stopLymanOpticalDepth))+"\n"
    cloudyScript += "iterate to convergence\n"
    ## Set overview output.
    if args.overview:
        cloudyScript += "save overview \"overview"+str(jobNumber)+".out\"\n"
    ## Output the continuum for reference.
    cloudyScript += "punch continuum \"continuum"+str(jobNumber)+".out\"\n"
    ## Set line output options.
    cloudyScript += "print lines faint _off\n"
    ## WORKAROUND     
    ## Currently disabled as Cloudy has a bug that prevents us from outputting in vacuum wavelengths. See these
    ## two error reports:
    ##   https://cloudyastrophysics.groups.io/g/Main/topic/101921207#5396
    ##   https://cloudyastrophysics.groups.io/g/Main/topic/102424985#5431
    #cloudyScript += "print line vacuum\n"
    cloudyScript += "save lines, array \"lines"+str(jobNumber)+".out\"\n"
    ## Write the Cloudy script to file.
    cloudyScriptFileName = "cloudyInput"+str(jobNumber)+".txt"
    cloudyScriptFile = open(args.workspace+cloudyScriptFileName,"w")
    cloudyScriptFile.write(cloudyScript)
    cloudyScriptFile.close()
    ## Construct a job to run this parameter file.
    job = {
        "label":                               "emissionLines"+str(jobNumber)       ,
        "launchFile":           args.workspace+"emissionLines"+str(jobNumber)+".sh" ,
        "logOutput":            args.workspace+"emissionLines"+str(jobNumber)+".log",
        "logError":             args.workspace+"emissionLines"+str(jobNumber)+".log",
        "nodes":                1                                                   ,
 	"cpusPerTask":          1                                                   ,
 	"walltime":             "10:00:00"                                          ,
        "onCompletion":         linesParse                                          ,
        "jobNumber":            jobNumber                                           ,
        "grid":                 grid                                                ,
 	"cloudyScriptFileName": args.workspace+cloudyScriptFileName                 ,
 	"linesFileName":        args.workspace+"lines"        +str(jobNumber)+".out",
 	"continuumFileName":    args.workspace+"continuum"    +str(jobNumber)+".out",
 	"indices":      	( iAge, iMetallicity, iNormalization, iLogHydrogenDensity ),
        "command":              "cd "+args.workspace+"; ulimit -c 0\n"+cloudyPath+"/source/cloudy.exe < "+cloudyScriptFileName+"\n"+"if [ $? != 0 ]; then\necho CLOUDY FAILED\nfi\n"
    }
    ## Push job to job list.
    if not 'jobs' in grid:
        grid['jobs'] = []
    grid['jobs'].append(job)

def generateJobAGN(grid,args):
    # Generate a single Cloudy job for AGN calculations.
    # Extract the indices for this job.
    iSpectralIndex       = grid['counter'][0]
    iMetallicity         = grid['counter'][1]
    iIonizationParameter = grid['counter'][2]
    iLogHydrogenDensity  = grid['counter'][3]
    # Get abundances for this metallicity.
    ## Compute metallicity relative to Solar.
    metallicity                                    = 10.0**grid['logMetallicities'][iMetallicity]
    ## Specify a dust-to-metals ratio. The following corresponds to the dust-to-metals ratio for the reference model of
    # Gutkin, Charlot & Bruzual (2016; https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G; table 1). It differs slightly
    # from their stated value, but agrees with our internal calculation of this value from their data.
    dustToMetals                                   = 0.401
    (dustToMetalsBoostLogarithmic, abundances)     = adjustAbundances(abundancesReference,metallicity,dustToMetals)
    # Create a depletion file that can be read by Cloudy. Cloudy does not permit digits in these file names. To work around this,
    # we translate our numerical metallicity index into ASCII characters, by mapping digits (0→A, 1→B, etc.).
    if not args.noGrains:
        encodedMetallicity = replaceDigits(iMetallicity)
        depletionFileName  = "grains_"+encodedMetallicity+".dpl"
        if not 'depletions' in grid:
            grid['depletions'] = [0]*len(grid['logMetallicities'])
        if grid['depletions'][iMetallicity] != 1:
            depletionFile = open(args.workspace+depletionFileName,"w")
            for element in abundances:
                depletionFile.write(abundances[element]['name']+" "+str(abundances[element]['undepletedFraction'])+"\n")
            depletionFile.close()
            grid['depletions'][iMetallicity] = 1
    # Generate a Cloudy parameter file.
    cloudyScript  = "title emission line job number "+str(jobNumber)+"\n"
    cloudyScript += "# ["+str(iSpectralIndex      )+"] alpha   = "+str(grid['spectralIndices'        ][iSpectralIndex      ])+"\n"
    cloudyScript += "# ["+str(iMetallicity        )+"] log Z   = "+str(grid['logMetallicities'       ][iMetallicity        ])+"\n"
    cloudyScript += "# ["+str(iIonizationParameter)+"] log U_H = "+str(grid['logIonizationParameters'][iIonizationParameter])+"\n"
    cloudyScript += "# ["+str(iLogHydrogenDensity )+"] log n_H = "+str(grid['logHydrogenDensities'   ][iLogHydrogenDensity ])+"\n"
    ## Set the input spectrum for Cloudy.
    counter = -1
    for iFrequency in range(len(grid['frequencies'])):
        counter += 1
        if counter % 3 == 0:
            if counter == 0:
                cloudyScript +=  "interpolate"
            else:
                cloudyScript +=  "\ncontinue"
        cloudyScript +=  " ("+str(grid['frequencies'][iFrequency])+" "+str(np.log10(grid['spectra'][iSpectralIndex][iFrequency]))+")"
    cloudyScript += "\n"
    ## Set normalization of the spectrum.
    cloudyScript += "ionization parameter = "+str(grid['logIonizationParameters'][iIonizationParameter])+"\n"
    # Set the chemical composition of the HII region.
    # Start with Cloudy's HII region abundances - but no grains - we will add these later.
    cloudyScript += "abundances HII region no grains\n"
    # Set individual element abundances.
    for element in abundances:
        # Skip hydrogen.
        if abundances[element]['atomicNumber'] < 2:
            continue
        # Set the abundance for this element.
        cloudyScript += "element abundances "+abundances[element]['name'].lower()+" "+str(abundances[element]['logAbundanceByNumber'])+"\n"
    # Specify Cloudy's default ISM grains, but with abundance reduced in proportion to the metallicity to retain a fixed
    # dust-to-metals ratio. Include the sublimation suppression function (see section 7.9.5 of Hazy1;
    # https://data.nublado.org/cloudy_releases/c23/c23.01.tar.gz).
    if not args.noGrains:
        cloudyScript += "grains Orion "+str(dustToMetalsBoostLogarithmic)+" _log function sublimation\n"
    # Deplete metals into grains using our custom depletions file. Note that these depletions differ from those assumed by
    # the "grains _ISM" model above. This seems to be at the ~10% level, so we do not worry too much about this.
    if not args.noGrains:
        cloudyScript += "metals deplete \""+depletionFileName+"\"\n"
    # Set HII region density - this is log₁₀(nₕ/cm¯³).
    cloudyScript += "hden "+str(grid['logHydrogenDensities'][iLogHydrogenDensity])+"\n"
    # Set cosmic rays (needed to avoid problems in Cloudy in neutral gas).
    cloudyScript += "cosmic rays background\n"
    # Set stopping criteria.
    ## Using temperature as a stopping criterion can be problematic as the model can then extend far into the neutral region if
    ## grains are present and cause heating. Instead, we stop based when an electron fraction (default of 1%) is reached, or when
    ## an Lyman limit optical depth (default of 10) is reached which should accurately capture the ionization front.
    cloudyScript += "stop temperature off\n"
    cloudyScript += "stop efrac "                +str(         args.stopElectronFraction  )+"\n"
    cloudyScript += "stop Lyman optical depth = "+str(np.log10(args.stopLymanOpticalDepth))+"\n"
    cloudyScript += "iterate to convergence\n"
    ## Output the continuum for reference.
    cloudyScript += "punch continuum \"continuum"+str(jobNumber)+".out\"\n"
    ## Set line output options.
    cloudyScript += "print lines faint _off\n"
    ## WORKAROUND     
    ## Currently disabled as Cloudy has a bug that prevents us from outputting in vacuum wavelengths. See these
    ## two error reports:
    ##   https://cloudyastrophysics.groups.io/g/Main/topic/101921207#5396
    ##   https://cloudyastrophysics.groups.io/g/Main/topic/102424985#5431
    #cloudyScript += "print line vacuum\n"
    cloudyScript += "save lines, array \"lines"+str(jobNumber)+".out\"\n"
    ## Write the Cloudy script to file.
    cloudyScriptFileName = "cloudyInput"+str(jobNumber)+".txt"
    cloudyScriptFile = open(args.workspace+cloudyScriptFileName,"w")
    cloudyScriptFile.write(cloudyScript)
    cloudyScriptFile.close()
    ## Construct a job to run this parameter file.
    job = {
        "label":                               "emissionLines"+str(jobNumber)       ,
        "launchFile":           args.workspace+"emissionLines"+str(jobNumber)+".sh" ,
        "logOutput":            args.workspace+"emissionLines"+str(jobNumber)+".log",
        "logError":             args.workspace+"emissionLines"+str(jobNumber)+".log",
        "nodes":                1                                                   ,
 	"cpusPerTask":          1                                                   ,
 	"walltime":             "10:00:00"                                          ,
        "onCompletion":         linesParse                                          ,
        "jobNumber":            jobNumber                                           ,
        "grid":                 grid                                                ,
 	"cloudyScriptFileName": args.workspace+cloudyScriptFileName                 ,
 	"linesFileName":        args.workspace+"lines"        +str(jobNumber)+".out",
 	"continuumFileName":    args.workspace+"continuum"    +str(jobNumber)+".out",
 	"indices":      	( iSpectralIndex, iMetallicity, iIonizationParameter, iLogHydrogenDensity ),
        "command":              "cd "+args.workspace+"; ulimit -c 0\n"+cloudyPath+"/source/cloudy.exe < "+cloudyScriptFileName+"\n"+"if [ $? != 0 ]; then\necho CLOUDY FAILED\nfi\n"
    }
    ## Push job to job list.
    if not 'jobs' in grid:
        grid['jobs'] = []
    grid['jobs'].append(job)


def validateSSP(grid,args):    
    # Validate the results of the Cloudy calculations for SSPs.
    # Validate the ratio L_Hα/Qₕ. From Osterbrock & Ferland (2006; https://ui.adsabs.harvard.edu/abs/2006agna.book.....O) we
    # expect:
    #
    #    (L_Hα/ergs s¯¹) / (Qₕ/photons s¯¹) = (α^eff_Hα/α_B) hν_Hα = 1.37 x 10¯¹²
    #    
    # where α^eff_Hα is the effective recombination coefficient for Hα, α_B is the case B recombination coefficient, and hν_Hα is
    # the energy of an Hα photon. The numerical value above is for a pure hydrogen/helium nebula, with temperature of 10,000K and
    # density n_e = 100 cm¯³.
    selectAge         = grid['ages'            ]                                            <  0.01                        # Consider young (< 10 Myr) populations.
    selectMetallicity = grid['logMetallicities']                                            == grid['logMetallicities'][0] # Select primordial metallicity.
    selection         = grid['lineData'        ]['status'][selectAge,selectMetallicity,:,:] == 0                           # Select successful models.
    # Compute the ratio.
    ratio = copy.deepcopy(grid['lineData']['balmerAlpha6565']['luminosity'])
    if args.normalization == "ionizingLuminosity":
        for i in range(ratio.shape[2]):
            ratio[:,:,i,:] /= 10.0**grid['logHydrogenLuminosities'][i]
    elif args.normalization == "massStellar":
        for i in range(ratio.shape[2]):
            for iMetallicity in range(ratio.shape[1]):
                for iAge in range(ratio.shape[0]):
                    logHydrogenLuminosity = grid['logStellarMasses'][i]+np.log10(grid['ionizingLuminosityPerMass'][iMetallicity,iAge])
                    ratio[iAge,iMetallicity,i,:] /= 10.0**logHydrogenLuminosity
    else:
        sys.exit("Expected `normalization` of `ionizingLuminosity` or `massStellar`")

    # Find the minimum and maximum ratios.
    ratioMinimum = np.min(ratio[selectAge,selectMetallicity,:,:][selection])
    ratioMaximum = np.max(ratio[selectAge,selectMetallicity,:,:][selection])
    # Report.
    ratioTarget = 1.37e-12;
    print(f"Validation: ratio L_Hα/Qₕ")
    print(f"   Expected: {ratioTarget:8.3e}")
    print(f"   Found range: {ratioMinimum:8.3e} to {ratioMaximum:8.3e}\n")

    # Validate the ionizing photon luminosity per star formation rate. Compute the ratio of star formation rate to ionizing
    # luminosity. This assumes a constant star formation rate of φ=1 M☉/yr. The ratio is then:
    #
    #    (SFR/M☉ yr¯¹)/(Qₕ/photons s¯¹) = φ / ∫₀᪲ φ Qₕ(t) dt
    #
    # a value of 7.4 x 10¯⁵⁴ is expected for a Kroupa IMF (Osterbrock & Ferland, 2006;
    # https://ui.adsabs.harvard.edu/abs/2006agna.book.....O). For a Chabrier IMF a lower value of around 4.7 x 10¯⁵⁴ is expected.
    ageStep       = copy.deepcopy(grid['ages'])
    ageStep[0:-2] = +grid['ages'][1:-1] \
                    -grid['ages'][0:-2]
    ageStep[  -1] = ageStep[-2]
    giga          = 1.0e9
    starFormationRateToIonizingLuminosity = 1.0/np.sum(grid['ionizingLuminosityPerMass'][0,:]*ageStep*giga)
    # Report.
    print( "Validation: ratio SFR/Qₕ")
    print( "   Expected (  Kroupa IMF): 7.54e-54")
    print( "   Expected (Chabrier IMF): 4.74e-54")
    print(f"   Found: {starFormationRateToIonizingLuminosity:8.3e}")

    # Validate BPT diagram location of the "standard" model of Gutkin, Charlot & Bruzual (2016;
    # https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G).
    # Their standard model has:
    #  Z_ISM   = 0.014
    #  ξ_d     = 0.28
    #  log U_S = −3.4
    #  n_H     = 10² cm⁻³
    # and they generate BPT diagrams assuming a constant star formation rate over the past 100 Myr.
    #
    # Set target ratios matching those from Gutkin, Charlot & Bruzual (2016).
    ratioNIIHalphaTarget = -0.60
    ratioOIIIHbetaTarget = +0.00
    ratioSIIHalphaTarget = -0.47
    ratioNIIOIITarget    = -0.61
    ratioOIIOIIITarget   = +0.48
    toleranceSuccess     = +0.20
    if args.normalization == "ionizingLuminosity":
        # Select model based on a distance metric.
        densityHydrogenStandard                =  100.0
        ionizationParameterLogarithmicStandard = -3.4
        metallicityStandard                    =  0.014/0.01524 # Convert to relative to Solar using Solar metallicity from Gutkin, Charlot & Bruzual (2016).
        radiiStromgrenLogarithmic              = (grid['logHydrogenLuminosities']+np.log10(3.0/4.0/np.pi/densityHydrogenStandard**2/coefficientRecombinationCaseB))/3.0
        ionizationParameterLogarithmic         =  grid['logHydrogenLuminosities']-np.log10(4.0*np.pi*(hecto*speedOfLight))-np.log10(densityHydrogenStandard)-2.0*radiiStromgrenLogarithmic
        distanceMetals                         = (grid['logMetallicities'       ]-np.log10(metallicityStandard                   ))**2
        distanceDensity                        = (grid['logHydrogenDensities'   ]-np.log10(densityHydrogenStandard               ))**2
        distanceIonizationParameter            = (ionizationParameterLogarithmic-          ionizationParameterLogarithmicStandard )**2
        selectMetals                           = np.argmin(distanceMetals             )
        selectDensity                          = np.argmin(distanceDensity            )
        selectIonizationParameter              = np.argmin(distanceIonizationParameter)
        selectAge                              = grid['ages'] < 100.0e-3
        deltaAge                               = copy.deepcopy(grid['ages'])
        deltaAge[0:-2]                         = grid['ages'][1:-1]-deltaAge[0:-2]
        deltaAge[  -1]                         =                    deltaAge[  -2]
        # Compute luminosities of all diagnostic lines.
        luminosityHalpha               = np.sum(grid['lineData']['balmerAlpha6565']['luminosity'][selectAge,selectMetals,selectIonizationParameter,selectDensity]*deltaAge[selectAge]*grid['ionizingLuminosityPerMass'][selectMetals,selectAge])
        luminosityHbeta                = np.sum(grid['lineData']['balmerBeta4863' ]['luminosity'][selectAge,selectMetals,selectIonizationParameter,selectDensity]*deltaAge[selectAge]*grid['ionizingLuminosityPerMass'][selectMetals,selectAge])
        luminosityOII                  = np.sum(grid['lineData']['oxygenII3727'   ]['luminosity'][selectAge,selectMetals,selectIonizationParameter,selectDensity]*deltaAge[selectAge]*grid['ionizingLuminosityPerMass'][selectMetals,selectAge])
        luminosityOIII                 = np.sum(grid['lineData']['oxygenIII5008'  ]['luminosity'][selectAge,selectMetals,selectIonizationParameter,selectDensity]*deltaAge[selectAge]*grid['ionizingLuminosityPerMass'][selectMetals,selectAge])
        luminosityNII                  = np.sum(grid['lineData']['nitrogenII6585' ]['luminosity'][selectAge,selectMetals,selectIonizationParameter,selectDensity]*deltaAge[selectAge]*grid['ionizingLuminosityPerMass'][selectMetals,selectAge])
        luminositySIIa                 = np.sum(grid['lineData']['sulfurII6718'   ]['luminosity'][selectAge,selectMetals,selectIonizationParameter,selectDensity]*deltaAge[selectAge]*grid['ionizingLuminosityPerMass'][selectMetals,selectAge])
        luminositySIIb                 = np.sum(grid['lineData']['sulfurII6733'   ]['luminosity'][selectAge,selectMetals,selectIonizationParameter,selectDensity]*deltaAge[selectAge]*grid['ionizingLuminosityPerMass'][selectMetals,selectAge])
        # Construct BPT diagram line ratios.
        ratioNIIHalpha                 = np.log10( luminosityNII                 /luminosityHalpha)
        ratioOIIIHbeta                 = np.log10( luminosityOIII                /luminosityHbeta )
        ratioSIIHalpha                 = np.log10((luminositySIIa+luminositySIIb)/luminosityHalpha)
        ratioNIIOII                    = np.log10( luminosityNII                 /luminosityOII   )
        ratioOIIOIII                   = np.log10( luminosityOII                 /luminosityOIII  )
        # Determine success.
        successNIIHalpha = "within acceptable range" if np.abs(ratioNIIHalpha-ratioNIIHalphaTarget) < toleranceSuccess else "NOT within acceptable range"
        successOIIIHbeta = "within acceptable range" if np.abs(ratioOIIIHbeta-ratioOIIIHbetaTarget) < toleranceSuccess else "NOT within acceptable range"
        successSIIHalpha = "within acceptable range" if np.abs(ratioSIIHalpha-ratioSIIHalphaTarget) < toleranceSuccess else "NOT within acceptable range"
        successNIIOII    = "within acceptable range" if np.abs(ratioNIIOII   -ratioNIIOIITarget   ) < toleranceSuccess else "NOT within acceptable range"
        successOIIOIII   = "within acceptable range" if np.abs(ratioOIIOIII  -ratioOIIOIIITarget  ) < toleranceSuccess else "NOT within acceptable range"
        # Report on agreement.
        print(f'BPT diagram ratio [NII]6584/Hα         (actual : target) = {ratioNIIHalpha:+.2f} : {ratioNIIHalphaTarget:+.2f}; {successNIIHalpha}')
        print(f'BPT diagram ratio [OIII]5007/Hβ        (actual : target) = {ratioOIIIHbeta:+.2f} : {ratioOIIIHbetaTarget:+.2f}; {successOIIIHbeta}')
        print(f'BPT diagram ratio [SII]6717,6731/Hα    (actual : target) = {ratioSIIHalpha:+.2f} : {ratioSIIHalphaTarget:+.2f}; {successSIIHalpha}')
        print(f'BPT diagram ratio [NII]6584/[OII]3727  (actual : target) = {ratioNIIOII   :+.2f} : {ratioNIIOIITarget   :+.2f}; {successNIIOII   }')
        print(f'BPT diagram ratio [OII]3727/[OIII]5007 (actual : target) = {ratioOIIOIII  :+.2f} : {ratioOIIOIIITarget  :+.2f}; {successOIIOIII  }')

    # Make plots of BPT diagrams.
    ## Construct line ratios.
    if args.normalization == "ionizingLuminosity":
        normalizationName = "logHydrogenLuminosities"
    elif args.normalization == "massStellar":
        normalizationName = "logStellarMasses"
    else:
        sys.exit("Expected `normalization` of `ionizingLuminosity` or `massStellar`")        
    with np.errstate(divide='ignore', invalid='ignore'):
        ratioNIIHalphaAll     = np.log10( grid['lineData']['nitrogenII6585' ]['luminosity']                                                /grid['lineData']['balmerAlpha6565']['luminosity'])
        ratioOIIIHbetaAll     = np.log10( grid['lineData']['oxygenIII5008'  ]['luminosity']                                                /grid['lineData']['balmerBeta4863' ]['luminosity'])
        ratioSIIHalphaAll     = np.log10((grid['lineData']['sulfurII6718'   ]['luminosity']+grid['lineData']['sulfurII6733']['luminosity'])/grid['lineData']['balmerAlpha6565']['luminosity'])
        ratioNIIOIIAll        = np.log10( grid['lineData']['nitrogenII6585' ]['luminosity']                                                /grid['lineData']['oxygenII3727'   ]['luminosity'])
        ratioOIIOIIIAll       = np.log10( grid['lineData']['oxygenII3727'   ]['luminosity']                                                /grid['lineData']['oxygenIII5008'  ]['luminosity'])
        ageAll                = np.broadcast_to(grid['ages'                ][   :,None,None,None],(len(grid['ages']),len(grid['logMetallicities']),len(grid[normalizationName]),len(grid['logHydrogenDensities'])))
        metallicityRatioAll   = np.broadcast_to(grid['logMetallicities'    ][None,   :,None,None],(len(grid['ages']),len(grid['logMetallicities']),len(grid[normalizationName]),len(grid['logHydrogenDensities'])))
        ionizingLuminosityAll = np.broadcast_to(grid[normalizationName     ][None,None,   :,None],(len(grid['ages']),len(grid['logMetallicities']),len(grid[normalizationName]),len(grid['logHydrogenDensities'])))
        densityHydrogenAll    = np.broadcast_to(grid['logHydrogenDensities'][None,None,None,   :],(len(grid['ages']),len(grid['logMetallicities']),len(grid[normalizationName]),len(grid['logHydrogenDensities'])))
    rng       = np.random.default_rng()
    selectAll = rng.choice(ratioNIIHalphaAll.ravel().shape[0], 10000)
    ## Read SDSS data.
    sdss = h5py.File(os.environ.get('GALACTICUS_DATA_PATH')+'/static/observations/emissionLines/emissionLineFluxesStarFormingSDSSDR8.hdf5','r')
    flux = {}
    for lineName in sdss:
        flux[lineName] = sdss[lineName][:]
    with np.errstate(divide='ignore', invalid='ignore'):
        ratioNIIHalphaSDSS = np.log10( flux['nitrogenII6585' ]                      /flux['balmerAlpha6565'])
        ratioOIIIHbetaSDSS = np.log10( flux['oxygenIII5008'  ]                      /flux['balmerBeta4863' ])
        ratioSIIHalphaSDSS = np.log10((flux['sulfurII6718'   ]+flux['sulfurII6733'])/flux['balmerAlpha6565'])
        ratioNIIOIISDSS    = np.log10( flux['nitrogenII6585' ]                      /flux['oxygenII3727'   ])
        ratioOIIOIIISDSS   = np.log10( flux['oxygenII3727'   ]                      /flux['oxygenIII5008'  ])
    selectSDSS         = rng.choice(ratioNIIHalphaSDSS.shape[0], 3000)
    ## OIII vs. NII
    figure, axes = plt.subplots(tight_layout=True)
    axes.set_xlabel("log([NII]6584/H$\\alpha$)" )
    axes.set_ylabel("log([OIII]5007/H$\\beta$)" )
    axes.set_xlim([-3.5,1.0])
    axes.set_ylim([-2.5,2.0])
    axes.scatter(ratioNIIHalphaSDSS        [selectSDSS],ratioOIIIHbetaSDSS        [selectSDSS],color="#AAAAAA",marker=".",s=1)
    axes.scatter(ratioNIIHalphaAll .ravel()[selectAll ],ratioOIIIHbetaAll .ravel()[selectAll ],c=metallicityRatioAll.ravel()[selectAll],cmap='viridis',norm='linear',s=1)
    if args.normalization == "ionizingLuminosity":
        axes.scatter(ratioNIIHalpha                        ,ratioOIIIHbeta                        ,s=100,color='black',marker="x")
    target = patches.Rectangle((ratioNIIHalphaTarget-toleranceSuccess, ratioOIIIHbetaTarget-toleranceSuccess), 2.0*toleranceSuccess, 2.0*toleranceSuccess, edgecolor='black', facecolor='none', linewidth=2)
    axes.add_patch(target)
    plt.savefig(args.workspace+'bptDiagramOIIINII.svg')  
    plt.clf()
    ## OIII vs. SII
    figure, axes = plt.subplots(tight_layout=True)
    axes.set_xlabel("log([SII]6716,6731/H$\\alpha$)" )
    axes.set_ylabel("log([OIII]5007/H$\\beta$)" )
    axes.set_xlim([-3.5,1.0])
    axes.set_ylim([-2.5,2.0])
    axes.scatter(ratioSIIHalphaSDSS        [selectSDSS],ratioOIIIHbetaSDSS        [selectSDSS],color="#AAAAAA",marker=".",s=1)
    axes.scatter(ratioSIIHalphaAll .ravel()[selectAll ],ratioOIIIHbetaAll .ravel()[selectAll ],c=metallicityRatioAll.ravel()[selectAll],cmap='viridis',norm='linear',s=1)
    if args.normalization == "ionizingLuminosity":
        axes.scatter(ratioSIIHalpha                        ,ratioOIIIHbeta                        ,s=100,color='black',marker="x")
    target = patches.Rectangle((ratioSIIHalphaTarget-toleranceSuccess, ratioOIIIHbetaTarget-toleranceSuccess), 2.0*toleranceSuccess, 2.0*toleranceSuccess, edgecolor='black', facecolor='none', linewidth=2)
    axes.add_patch(target)
    plt.savefig(args.workspace+'bptDiagramOIIISII.svg')  
    plt.clf()
    ## OIII vs. NIIOII
    figure, axes = plt.subplots(tight_layout=True)
    axes.set_xlabel("log([NII]6584/[OII]3727)" )
    axes.set_ylabel("log([OIII]5007/H$\\beta$)" )
    axes.set_xlim([-3.0,1.5])
    axes.set_ylim([-2.5,2.0])
    axes.scatter(ratioNIIOIISDSS        [selectSDSS],ratioOIIIHbetaSDSS        [selectSDSS],color="#AAAAAA",marker=".",s=1)
    axes.scatter(ratioNIIOIIAll .ravel()[selectAll ],ratioOIIIHbetaAll .ravel()[selectAll ],c=metallicityRatioAll.ravel()[selectAll],cmap='viridis',norm='linear',s=1)
    if args.normalization == "ionizingLuminosity":
        axes.scatter(ratioNIIOII                        ,ratioOIIIHbeta                        ,s=100,color='black',marker="x")
    target = patches.Rectangle((ratioNIIOIITarget-toleranceSuccess, ratioOIIIHbetaTarget-toleranceSuccess), 2.0*toleranceSuccess, 2.0*toleranceSuccess, edgecolor='black', facecolor='none', linewidth=2)
    axes.add_patch(target)
    plt.savefig(args.workspace+'bptDiagramNIIOII.svg')  
    plt.clf()
    ## OIII vs. OIIOII
    figure, axes = plt.subplots(tight_layout=True)
    axes.set_xlabel("log([OII]3727/[OIII]5007)" )
    axes.set_ylabel("log([OIII]5007/H$\\beta$)" )
    axes.set_xlim([-2.5,2.0])
    axes.set_ylim([-2.5,2.0])
    axes.scatter(ratioOIIOIIISDSS        [selectSDSS],ratioOIIIHbetaSDSS        [selectSDSS],color="#AAAAAA",marker=".",s=1)
    axes.scatter(ratioOIIOIIIAll .ravel()[selectAll ],ratioOIIIHbetaAll .ravel()[selectAll ],c=metallicityRatioAll.ravel()[selectAll],cmap='viridis',norm='linear',s=1)
    if args.normalization == "ionizingLuminosity":
        axes.scatter(ratioOIIOIII                        ,ratioOIIIHbeta                        ,s=100,color='black',marker="x")
    target = patches.Rectangle((ratioOIIOIIITarget-toleranceSuccess, ratioOIIIHbetaTarget-toleranceSuccess), 2.0*toleranceSuccess, 2.0*toleranceSuccess, edgecolor='black', facecolor='none', linewidth=2)
    axes.add_patch(target)
    plt.savefig(args.workspace+'bptDiagramOIIOIII.svg')  
    plt.clf()

    # Update results in the datasets repo GitHub pages.
    if args.suffixGitHubPages is not None:
        ## Clone the repo (gh-pages branch).
        if not os.path.isdir(args.workspace+'datasets'):
            try:
                Repo.clone_from("git@github.com:galacticusorg/datasets.git", args.workspace+'datasets', branch='gh-pages')
            except Exception as e:
                print(f"Error cloning repository: {e}")
        ## Parse the JSON definition file if present.
        if os.path.exists(args.workspace+'datasets/hiiRegions/tableDefinitions.json'):
            with open(args.workspace+'datasets/hiiRegions/tableDefinitions.json', 'r') as file:
                next(file)
                data = file.read()
                definition = json.loads(data)
        else:
            definition = {}    
        ## Copy and rename our plots to the repo.
        if not os.path.isdir(args.workspace+'datasets/hiiRegions'):
            os.mkdir(args.workspace+'datasets/hiiRegions')
        shutil.copy(args.workspace+'bptDiagramOIIINII.svg', args.workspace+'datasets/hiiRegions/bptDiagramOIIINII_'+args.suffixGitHubPages+'.svg')
        shutil.copy(args.workspace+'bptDiagramOIIISII.svg', args.workspace+'datasets/hiiRegions/bptDiagramOIIISII_'+args.suffixGitHubPages+'.svg')
        shutil.copy(args.workspace+'bptDiagramNIIOII.svg' , args.workspace+'datasets/hiiRegions/bptDiagramNIIOII_' +args.suffixGitHubPages+'.svg')
        shutil.copy(args.workspace+'bptDiagramOIIOIII.svg', args.workspace+'datasets/hiiRegions/bptDiagramOIIOIII_'+args.suffixGitHubPages+'.svg')
        ## Update and output the JSO1N definition file.
        definition[args.suffixGitHubPages] = {
            "time":          str(datetime.datetime.now()),
            "gitRevision":   grid['gitRevision'  ],
            "cloudyVersion": grid['cloudyVersion'],
            "sspURL":        grid['sspURL'       ],
            "commandLine":   grid['commandLine'  ],
            "fileName":      args.outputFileName
        }
        f = codecs.open(args.workspace+'datasets/hiiRegions/tableDefinitions.json', "w", "utf-8")
        f.write("window.TABLE_DATA =\n")
        f.write(json.dumps(definition,indent=4,ensure_ascii=False))
        f.close()
        ## Display instructions for updating the repo.
        print(colored('\n\nGitHub pages content has been updated.','red', attrs=["bold"])+' To deploy, do:\n')
        print(colored('   cd '+args.workspace+'datasets','green'))
        print(colored('   git add hiiRegions/tableDefinitions.json','green'))
        print(colored('   git add hiiRegions/bptDiagramOIIINII_'+args.suffixGitHubPages+'.svg','green'))
        print(colored('   git add hiiRegions/bptDiagramOIIISII_'+args.suffixGitHubPages+'.svg','green'))
        print(colored('   git add hiiRegions/bptDiagramNIIOII_' +args.suffixGitHubPages+'.svg','green'))
        print(colored('   git add hiiRegions/bptDiagramOIIOIII_'+args.suffixGitHubPages+'.svg','green'))
        print(colored('   git commit -m "feat: Update BPT diagrams"','green'))
        print(colored('   git push\n\n','green'))
    # Validate model success.
    selectSuccess       = grid['lineData']['status'] == 0
    selectDisaster      = grid['lineData']['status'] == 1
    selectExit          = grid['lineData']['status'] == 2
    selectMissingOutput = grid['lineData']['status'] == 3
    selectMissingLine   = grid['lineData']['status'] == 4
    modelsTotal         = grid['lineData']['status'].size
    modelsSuccess       = np.count_nonzero(selectSuccess      )
    modelsDisaster      = np.count_nonzero(selectDisaster     )
    modelsExit          = np.count_nonzero(selectExit         )
    modelsMissingOutput = np.count_nonzero(selectMissingOutput)
    modelsMissingLine   = np.count_nonzero(selectMissingLine  )
    print("Validation: model success")
    print("      Total models: "+str(modelsTotal        ))
    print("        Successful: "+str(modelsSuccess      ))
    print("         Disasters: "+str(modelsDisaster     ))
    print("    Bad exit codes: "+str(modelsExit         ))
    print("   Missing outputs: "+str(modelsMissingOutput))
    print("     Missing lines: "+str(modelsMissingLine  ))

def validateAGN(grid,args):
    # Validate the results of the Cloudy calculations for AGN.
    # Validate the ratio L_Hα/Qₕ. From Osterbrock & Ferland (2006; https://ui.adsabs.harvard.edu/abs/2006agna.book.....O) we
    # expect:
    #
    #    (L_Hα/ergs s¯¹) / (Qₕ/photons s¯¹) = (α^eff_Hα/α_B) hν_Hα = 1.37 x 10¯¹²
    #    
    # where α^eff_Hα is the effective recombination coefficient for Hα, α_B is the case B recombination coefficient, and hν_Hα is
    # the energy of an Hα photon. The numerical value above is for a pure hydrogen/helium nebula, with temperature of 10,000K and
    # density n_e = 100 cm¯³.
    #
    # For this AGN model the source is parameterized in terms of the ionization parameter:
    #
    #  U = Q / 4π r² n c
    #
    # which means that:
    #
    #  Q = 4π U r² n c
    #
    # In this case, the line output from Cloudy is the "energy radiated by a unit area of cloud into 4 π sr (4πJ, erg cm¯²
    # s¯¹). The line luminosity (assuming unit covering factor) is then:
    #
    #  L_Hα = J 4π r²
    #
    # Therefore:
    #
    #  (L_Hα/ergs s¯¹) / (Qₕ/photons s¯¹) = J 4π r² / 4π r² n c = J / U n c
    #
    # where c should be in units of cm s¯¹.
    selectMetallicity = grid['logMetallicities']                                    == grid['logMetallicities'][0] # Select primordial metallicity.
    selection         = grid['lineData'        ]['status'][:,selectMetallicity,:,:] == 0                           # Select successful models.
    # Compute the ratio.
    ratio = copy.deepcopy(grid['lineData']['balmerAlpha6565']['luminosity'])
    for i in range(ratio.shape[2]):
        for j in range(ratio.shape[3]):
            ratio[:,:,i,j] /= 10.0**(grid['logIonizationParameters'][i]+grid['logHydrogenDensities'][j])*speedOfLight*hecto
    # Find the minimum and maximum ratios.
    ratioMinimum = np.min(ratio[:,selectMetallicity,:,:][selection])
    ratioMaximum = np.max(ratio[:,selectMetallicity,:,:][selection])
    # Report.
    ratioTarget = 1.37e-12
    print(f"Validation: ratio L_Hα/Qₕ")
    print(f"   Expected: {ratioTarget:8.3e}")
    print(f"   Found range: {ratioMinimum:8.3e} to {ratioMaximum:8.3e}\n")

    # Validate model success.
    selectSuccess       = grid['lineData']['status'] == 0
    selectDisaster      = grid['lineData']['status'] == 1
    selectExit          = grid['lineData']['status'] == 2
    selectMissingOutput = grid['lineData']['status'] == 3
    selectMissingLine   = grid['lineData']['status'] == 4
    modelsTotal         = grid['lineData']['status'].size
    modelsSuccess       = np.count_nonzero(selectSuccess      )
    modelsDisaster      = np.count_nonzero(selectDisaster     )
    modelsExit          = np.count_nonzero(selectExit         )
    modelsMissingOutput = np.count_nonzero(selectMissingOutput)
    modelsMissingLine   = np.count_nonzero(selectMissingLine  )
    print("Validation: model success")
    print("      Total models: "+str(modelsTotal        ))
    print("        Successful: "+str(modelsSuccess      ))
    print("         Disasters: "+str(modelsDisaster     ))
    print("    Bad exit codes: "+str(modelsExit         ))
    print("   Missing outputs: "+str(modelsMissingOutput))
    print("     Missing lines: "+str(modelsMissingLine  ))

def outputSSP(grid,args):
    # Output the results of the Cloudy calculations for SSPs.
    # Write the line data to file.
    tableFile = h5py.File(args.workspace+args.outputFileName,"w")
    # Add useful metadata.
    tableFile.attrs['time'         ] = str(datetime.datetime.now())
    tableFile.attrs['gitRevision'  ] = grid['gitRevision'  ]
    tableFile.attrs['commandLine'  ] = grid['commandLine'  ]
    tableFile.attrs['cloudyVersion'] = grid['cloudyVersion']
    # Write parameter grid points and attributes.
    datasetAge             = tableFile.create_dataset('age'                       ,data=      grid['ages'                   ])
    datasetMetallicity     = tableFile.create_dataset('metallicity'               ,data=10.0**grid['logMetallicities'       ])
    datasetDensityHydrogen = tableFile.create_dataset('densityHydrogen'           ,data=10.0**grid['logHydrogenDensities'   ])
    datasetAge            .attrs['description'] = "Age of the stellar population."
    datasetAge            .attrs['units'      ] = "Gyr"
    datasetAge            .attrs['unitsInSI'  ] = secondsPerGyr
    datasetMetallicity    .attrs['description'] = "Metallicity relative to Solar."
    datasetDensityHydrogen.attrs['description'] = "Hydrogen density."
    datasetDensityHydrogen.attrs['units'      ] = "cm¯³"
    datasetDensityHydrogen.attrs['unitsInSI'  ] = mega
    if args.normalization == "ionizingLuminosity":
        datasetIonizationLuminosityHydrogen = tableFile.create_dataset('ionizingLuminosityHydrogen',data=10.0**grid['logHydrogenLuminosities'])
        datasetIonizationLuminosityHydrogen.attrs['description'] = "Hydrogen ionizing photon emission rate."
        datasetIonizationLuminosityHydrogen.attrs['units'      ] = "photons s¯¹"
        datasetIonizationLuminosityHydrogen.attrs['unitsInSI'  ] = one
    elif args.normalization == "massStellar":
        datasetMassStellar = tableFile.create_dataset('massStellar',data=10.0**grid['logStellarMasses'])
        datasetMassStellar.attrs['description'] = "Zero-age stellar mass of the HII region."
        datasetMassStellar.attrs['units'      ] = "M☉"
        datasetMassStellar.attrs['unitsInSI'  ] = massSolar
    else:
        sys.exit("Expected `normalization` of `ionizingLuminosity` or `massStellar`")
    # Write index in the tables for each iterable.
    i = 0
    for iterable in grid['names']:
        tableFile[iterable].attrs['index'] = i
        i += 1
    # Write table of ionizing rates per unit mass of stars formed.
    datasetNormalization = tableFile.create_dataset('ionizingLuminosityHydrogenNormalized',data=grid['ionizingLuminosityPerMass'])
    datasetNormalization.attrs['description'] = "Hydrogen ionizing photon emission rate per unit mass of stars."
    datasetNormalization.attrs['units'      ] = "photons s¯¹ M☉¯¹"
    datasetNormalization.attrs['unitsInSI'  ] = 1.0/massSolar    
    # Write line data. Note that the grids are transposed for output in order to match the ordering produced by the original
    # (Perl-based) version of this script.
    lineGroup = tableFile.create_group('lines')
    datasetStatus = lineGroup.create_dataset('status'     ,data=np.transpose(grid['lineData']['status'     ]))
    datasetModel  = lineGroup.create_dataset('modelNumber',data=np.transpose(grid            ['modelNumber']))
    datasetStatus.attrs['description'] = "Cloudy model status: 0 = success; 1 = disaster; 2 = non-zero exit status; 3 = missing output file; 4 = missing emission lines"
    datasetModel .attrs['description'] = "Cloudy model number"
    for lineLabel in lineList:
        lineName    = lineList[lineLabel]
        datasetLine = lineGroup.create_dataset(lineName,data=np.transpose(grid['lineData'][lineName]['luminosity']))
        datasetLine.attrs['description'] = "Energy radiated by a unit area of cloud into 4 π sr."
        datasetLine.attrs['units'      ] = "erg cm¯² s¯¹"
        datasetLine.attrs['unitsInSI'  ] = unitsIntensity
        datasetLine.attrs['wavelength' ] = grid['lineData'][lineName]['wavelength']

def outputAGN(grid,args):
    # Output the results of the Cloudy calculations for AGN.
    # Write the line data to file.
    tableFile = h5py.File(args.workspace+args.outputFileName,'w')
    # Add useful metadata.
    tableFile.attrs['time'         ] = str(datetime.datetime.now())
    tableFile.attrs['gitRevision'  ] = grid['gitRevision'  ]
    tableFile.attrs['commandLine'  ] = grid['commandLine'  ]
    tableFile.attrs['cloudyVersion'] = grid['cloudyVersion']
    # Write parameter grid points and attributes.
    datasetSpectralIndex       = tableFile.create_dataset('spectralIndex'      ,data=      grid['spectralIndices'        ])
    datasetMetallicity         = tableFile.create_dataset('metallicity'        ,data=10.0**grid['logMetallicities'       ])
    datasetIonizationParameter = tableFile.create_dataset('ionizationParameter',data=10.0**grid['logIonizationParameters'])
    datasetDensityHydrogen     = tableFile.create_dataset('densityHydrogen'    ,data=10.0**grid['logHydrogenDensities'   ])
    datasetSpectralIndex      .attrs['description'] = "Spectral index at optical/UV wavelength population."
    datasetMetallicity        .attrs['description'] = "Metallicity relative to Solar."
    datasetIonizationParameter.attrs['description'] = "Ionization parameter."
    datasetDensityHydrogen    .attrs['description'] = "Hydrogen density."
    datasetDensityHydrogen    .attrs['units'      ] = "cm¯³"
    datasetDensityHydrogen    .attrs['unitsInSI'  ] = mega
    # Write index in the tables for each iterable.
    i = 0
    for iterable in grid['names']:
        tableFile[iterable].attrs['index'] = i
        i += 1
    # Write line data. Note that the grids are transposed for output in order to match the ordering produced by the original
    # (Perl-based) version of this script.
    lineGroup = tableFile.create_group('lines')
    datasetStatus = lineGroup.create_dataset('status'     ,data=np.transpose(grid['lineData']['status'     ]))
    datasetModel  = lineGroup.create_dataset('modelNumber',data=np.transpose(grid            ['modelNumber']))
    datasetStatus.attrs['description'] = "Cloudy model status: 0 = success; 1 = disaster; 2 = non-zero exit status; 3 = missing output file; 4 = missing emission lines"
    datasetModel .attrs['description'] = "Cloudy model number"
    for lineLabel in lineList:
        lineName    = lineList[lineLabel]
        datasetLine = lineGroup.create_dataset(lineName,data=np.transpose(grid['lineData'][lineName]['luminosity']))
        datasetLine.attrs['description'] = "Energy radiated by a unit area of cloud into 4 π sr."
        datasetLine.attrs['units'      ] = "erg cm¯² s¯¹"
        datasetLine.attrs['unitsInSI'  ] = unitsIntensity
        datasetLine.attrs['wavelength' ] = grid['lineData'][lineName]['wavelength']


# Parse command line arguments.
parser = argparse.ArgumentParser(prog='analysesPlotcreateEmissionLinesTable.py',description='Generate tables of Cloudy models for use in emission line calculations.')
parser.add_argument('--outputFileName'                                    ,action='store'                            ,help='the file to which the table should be output'                                                                 )
parser.add_argument('--sspFileName'                                       ,action='store'                            ,help='the SSP file for which to compute emission line luminosities'                                                 )
parser.add_argument('--agnModel'                                          ,action='store'                            ,help='the AGN model for which to compute emission line luminosities'                                                )
parser.add_argument('--workspace'            ,default='cloudyTable/'      ,action='store'                            ,help='the path in which temporary files should be created'                                                          )
parser.add_argument('--reprocess'                                         ,action='store_true'                       ,help='reprocess models that failed to be read previously'                                                           )
parser.add_argument('--rerun'                                             ,action='store_true'                       ,help='rerun models that previously failed'                                                                          )
parser.add_argument('--generateOnly'                                      ,action='store_true'                       ,help='only generate model input files, do not run them'                                                             )
parser.add_argument('--overview'                                          ,action='store_true'                       ,help='include the Cloudy overview in the output'                                                                    )
parser.add_argument('--noClean'                                           ,action='store_true'                       ,help='do not clean up temporary files'                                                                              )
parser.add_argument('--noGrains'                                          ,action='store_true'                       ,help='do not include dust grains in the models'                                                                     )
parser.add_argument('--normalization'        ,default='ionizingLuminosity',action='store'                            ,help='specify how the spectrum is to be normalized (`ionizingLuminosity` or `massStellar`)'                         )
parser.add_argument('--factorMorphology'     ,default='1.0'               ,action='store'      ,type=restricted_float,help='set the morphology factor (f=R_{in}/R_{Strömgren}; https://ui.adsabs.harvard.edu/abs/2016A%2526A...594A..37M)')
parser.add_argument('--stopOuterRadius'                                   ,action='store_true'                       ,help='set Cloudy to stop at the cloud outer radius'                                                                 )
parser.add_argument('--stopElectronFraction' ,default='0.01'              ,action='store'      ,type=restricted_float,help='set the elctron fraction at which to stop the Cloudy models'                                                  )
parser.add_argument('--stopLymanOpticalDepth',default='10.0'              ,action='store'      ,type=restricted_float,help='set the Lyman optical depth at which to stop the Cloudy models'                                               )
parser.add_argument('--cloudyVersion'        ,default="23.01"             ,action='store'                            ,help='the version of Cloudy to use'                                                                                 )
parser.add_argument('--suffixGitHubPages'                                 ,action='store'                            ,help='update GitHub pages content using this suffix'                                                                )
parser.add_argument('--model'                                             ,action='store'      ,type=restricted_int  ,help='run only the given model number'                                                                              )
parser.add_argument('--partition'                                         ,action='store'                            ,help='the partition to which to submit jobs'                                                                        )
parser.add_argument('--jobMaximum'                                        ,action='store'      ,type=restricted_int  ,help='the maximum number of active jobs to allow'                                                                   )
parser.add_argument('--waitOnSubmit'                                      ,action='store'      ,type=restricted_int  ,help='the time (in seconds) to wait after submitting each job'                                                      )
parser.add_argument('--waitOnActive'                                      ,action='store'      ,type=restricted_int  ,help='the time (in seconds) to wait after polling active jobs'                                                      )
args = parser.parse_args()

# Validate options.
if args.outputFileName is None:
    sys.exit("An output file name must be specified via the `--outputFileName` option")
if args.sspFileName is not None and args.agnModel is not None:
    sys.exit("Can not specify both  `--sspFileName` and `--agnModel`")
if args.sspFileName is not None:
    print("Computing emission line luminosities for stellar populations using file:\n\t"+args.sspFileName)
    establishGrid = establishGridSSP
    generateJob   = generateJobSSP
    reprocess     = reprocessSSP
    validate      = validateSSP
    output        = outputSSP
elif args.agnModel is not None:
    print("Computing emission line luminosities for AGN using model:\n\t"               +args.agnModel   )
    establishGrid = establishGridAGN
    generateJob   = generateJobAGN
    reprocess     = reprocessAGN 
    validate      = validateAGN
    output        = outputAGN
else:
    sys.exit("Specify either `--sspFileName` or `--agnModel`")

# Move any existing output file to a backup. This avoids attempting to overwrite HDF5 datasets with arrays of different size.
if not args.generateOnly and not args.reprocess and not args.rerun:
    if os.path.isfile(args.workspace+args.outputFileName):
        print("Moving old output file to `"+args.outputFileName+".bak`")
        os.rename(args.workspace+args.outputFileName,args.workspace+args.outputFileName+".bak")

# Create workspace.
if not args.workspace.endswith(os.sep):
    args.workspace += os.sep
os.makedirs(args.workspace, exist_ok=True)

# Initalize Cloudy.
options = {"version": args.cloudyVersion}
(cloudyPath, cloudyVersion) = cloudy.initialize(options)

# Set physical constants.
plancksConstant               = 6.6260700400000e-34 # J s
speedOfLight                  = 2.9979245800000e+08 # m s¯¹
coefficientRecombinationCaseB = 2.6000000000000e-13 # cm³ s⁻¹
electronVolt                  = 1.6021766340000e-19 # J
rydbergEnergy                 = 1.3605693122994e+01 # eV
micron                        = 1.0000000000000e-06 # μm
angstroms                     = 1.0000000000000e-10 # m
luminositySolar               = 3.8390000000000e+26 # W
wavelengthLymanContinuum      = 9.1176000000000e+02 # Å
massSolar                     = 1.9900000000000e+30 # M☉
one                           = 1.0000000000000e+00
hecto                         = 1.0000000000000e+02
mega                          = 1.0000000000000e+06
joulesPerErg                  = 1.0000000000000e-07
secondsPerGyr                 = 3.1557600000000e-16
unitsIntensity                = joulesPerErg*hecto**2

# Specify abundances and depletion model. This is based upon the work by Gutkin, Charlot & Bruzual (2016;
# https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1757G).
abundancesReference = {
     "H": {
	 "atomicNumber"        :   1,
	 "logAbundanceByNumber":   0.00,
	 "undepletedFraction"  :   1.000
     },
     "He": {
	 "atomicNumber"        :   2,
	 "logAbundanceByNumber":  -1.01,
	 "undepletedFraction"  :   1.000,
	 "adjustAbundance"     :   adjustAbundanceHelium
     },
     "Li": {
	 "atomicNumber"        :   3,
	 "logAbundanceByNumber": -10.99,
	 "undepletedFraction"  :   0.160
     },
     "Be": {
	 "atomicNumber"        :   4,
	 "logAbundanceByNumber": -10.63,
	 "undepletedFraction"  :   0.600
     },
     "B": {
	 "atomicNumber"        :   5,
	 "logAbundanceByNumber":  -9.47,
	 "undepletedFraction"  :   0.130
     },
     "C": {
	 "atomicNumber"        :   6,
	 "logAbundanceByNumber":  -3.53,
	 "undepletedFraction"  :   0.500,
	 "adjustAbundance"     :  adjustAbundanceCarbon
     },
     "N": {
	 "atomicNumber"        :   7,
	 "logAbundanceByNumber":  -4.32,
	 "undepletedFraction"  :   1.000,
	 "adjustAbundance"     :   adjustAbundanceNitrogen
     },
     "O": {
	 "atomicNumber"        :   8,
	 "logAbundanceByNumber":  -3.17,
	 "undepletedFraction"  :   0.700
     },
     "F": {
	 "atomicNumber"        :   9,
	 "logAbundanceByNumber":  -7.47,
	 "undepletedFraction"  :   0.300
     },
     "Ne": {
	 "atomicNumber"        :  10,
	 "logAbundanceByNumber":  -4.01,
	 "undepletedFraction"  :   1.000
     },
     "Na": {
	 "atomicNumber"        :  11,
	 "logAbundanceByNumber":  -5.70,
	 "undepletedFraction"  :   0.250
     },
     "Mg": {
	 "atomicNumber"        :  12,
	 "logAbundanceByNumber":  -4.45,
	 "undepletedFraction"  :   0.200
     },
     "Al": {
	 "atomicNumber"        :  13,
	 "logAbundanceByNumber":  -5.56,
	 "undepletedFraction"  :   0.020
     },
     "Si": {
	 "atomicNumber"        :  14,
	 "logAbundanceByNumber":  -4.48,
	 "undepletedFraction"  :   0.100
     },
     "P": {
	 "atomicNumber"        :  15,
	 "logAbundanceByNumber":  -6.57,
	 "undepletedFraction"  :   0.250
     },
     "S": {
	 "atomicNumber"        :  16,
	 "logAbundanceByNumber":  -4.87,
	 "undepletedFraction"  :   1.000
     },
     "Cl": {
	 "atomicNumber"        :  17,
	 "logAbundanceByNumber":  -6.53,
	 "undepletedFraction"  :   0.500
     },
     "Ar": {
	 "atomicNumber"        :  18,
	 "logAbundanceByNumber":  -5.63,
	 "undepletedFraction"  :   1.000
     },
     "K": {
	 "atomicNumber"        :  19,
	 "logAbundanceByNumber":  -6.92,
	 "undepletedFraction"  :   0.300
     },
     "Ca": {
	 "atomicNumber"        :  20,
	 "logAbundanceByNumber":  -5.67,
	 "undepletedFraction"  :   0.003
     },
     "Sc": {
	 "atomicNumber"        :  21,
	 "logAbundanceByNumber":  -8.86,
	 "undepletedFraction"  :   0.005
     },
     "Ti": {
	 "atomicNumber"        :  22,
	 "logAbundanceByNumber":  -7.01,
	 "undepletedFraction"  :   0.008
     },
     "V": {
	 "atomicNumber"        :  23,
	 "logAbundanceByNumber":  -8.03,
	 "undepletedFraction"  :   0.006
     },
     "Cr": {
	 "atomicNumber"        :  24,
	 "logAbundanceByNumber":  -6.36,
	 "undepletedFraction"  :   0.006
     },
     "Mn": {
	 "atomicNumber"        :  25,
	 "logAbundanceByNumber":  -6.64,
	 "undepletedFraction"  :   0.050
     },
     "Fe": {
	 "atomicNumber"        :  26,
	 "logAbundanceByNumber":  -4.51,
	 "undepletedFraction"  :   0.010
     },
     "Co": {
	 "atomicNumber"        :  27,
	 "logAbundanceByNumber":  -7.11,
	 "undepletedFraction"  :   0.010
     },
     "Ni": {
	 "atomicNumber"        :  28,
	 "logAbundanceByNumber":  -5.78,
	 "undepletedFraction"  :   0.040
     },
     "Cu": {
	 "atomicNumber"        :  29,
	 "logAbundanceByNumber":  -7.82,
	 "undepletedFraction"  :   0.100
     },
     "Zn": {
	 "atomicNumber"        :  30,
	 "logAbundanceByNumber":  -7.43,
	 "undepletedFraction"  :   0.250
     }
}

# Augment abundance data with names and atomic masses. These are read from Galacticus' atomic data file.
atomicData = ET.parse(os.environ['GALACTICUS_DATA_PATH']+"/static/abundances/Atomic_Data.xml")
for element in atomicData.findall('./element'):
    elementName = element.find('./name'      )
    shortLabel  = element.find('./shortLabel')
    atomicMass  = element.find('./atomicMass')
    if shortLabel is not None and shortLabel.text in abundancesReference:
 	# Cloudy uses British spelling (presumably because it originated at Cambridge) - so correct for that here so that we can
 	# use the correct names when generating a Cloudy script.
        nameBritish = "Sulphur" if elementName.text == "Sulfur" else  elementName.text
 	# Store the name and atomic mass for this element.
        abundancesReference[shortLabel.text]['atomicMass'] = float(atomicMass.text)
        abundancesReference[shortLabel.text]['name'      ] = nameBritish
        
# Define the list of lines to extract. The following dictionary contains keys which match the line names in the Cloudy emission
# lines output file, and values which are our internal names for these lines.
lineList = {
    ## WORKAROUND     
    ## Currently using in-air wavelengths as Cloudy has a bug that prevents us from outputting in vacuum wavelengths. See these
    ## two error reports:
    ##   https://cloudyastrophysics.groups.io/g/Main/topic/101921207#5396
    ##   https://cloudyastrophysics.groups.io/g/Main/topic/102424985#5431
    # In air wavelengths.
    "H  1                6562.80A": "balmerAlpha6565"  ,
    "H  1                4861.32A": "balmerBeta4863"   ,
    "H  1                4340.46A": "balmerGamma4342"  ,
    "H  1                4101.73A": "balmerDelta4103"  ,
    "H  1                1.87510m": "paschenAlpha18756",
    "H  1                1.28181m": "paschenBeta12822" ,
    "O  2                3726.03A": "oxygenII3727"     ,
    "O  2                3728.81A": "oxygenII3730"     ,
    "O  3                4958.91A": "oxygenIII4960"    ,
    "O  3                5006.84A": "oxygenIII5008"    ,
    "O  3                4931.23A": "oxygenIII4933"    ,
    "N  2                6583.45A": "nitrogenII6585"   ,
    "N  2                6548.05A": "nitrogenII6550"   ,
    "S  2                6730.82A": "sulfurII6733"     ,
    "S  2                6716.44A": "sulfurII6718"     ,
    "S  3                9068.62A": "sulfurIII9071"    ,
    "S  3                9530.62A": "sulfurIII9533"
     # In vacuum wavelengths.
     # "H  1                6564.62A": "balmerAlpha6565"  ,
     # "H  1                4862.69A": "balmerBeta4863"   ,
     # "H  1                4341.68A": "balmerGamma4342"  ,
     # "H  1                4102.89A": "balmerDelta4103"  ,
     # "H  1                1.87561m": "paschenAlpha18756",
     # "H  1                1.28215m": "paschenBeta12822" ,
     # "O  2                3727.09A": "oxygenII3727"     ,
     # "O  2                3729.88A": "oxygenII3730"     ,
     # "O  3                4960.29A": "oxygenIII4960"    ,
     # "O  3                5008.24A": "oxygenIII5008"    ,
     # "O  3                4932.60A": "oxygenIII4933"    ,
     # "N  2                6585.27A": "nitrogenII6585"   ,
     # "N  2                6549.86A": "nitrogenII6550"   ,
     # "S  2                6732.67A": "sulfurII6733"     ,
     # "S  2                6718.29A": "sulfurII6718"     ,
     # "S  3                9071.11A": "sulfurIII9071"    ,
     # "S  3                9533.23A": "sulfurIII9533"
}

# Establish the model grid.
grid = {}
grid['type'] = "SSP" if args.sspFileName is not None else "AGN"
establishGrid(grid,args)

# Find the current git revision hash.
repo                = Repo(os.environ['GALACTICUS_EXEC_PATH'])
lastRevision        = repo.head.object.hexsha
grid['gitRevision'] = lastRevision

# Add the command line arguments for possible output.
grid['commandLine'] = " ".join(sys.argv)

# Store Cloudy version.
grid['cloudyVersion'] = cloudyVersion

# Initialize the line luminosity tables.
dimensions       = list(map(lambda x: grid[x].size,grid['iterables']))
grid['lineData'] = {}
for lineName in lineList.keys():
    grid['lineData'][lineList[lineName]] = {"luminosity": np.zeros(dimensions)}
grid['lineData']['status'] = np.zeros(dimensions,dtype=int)

# Iterate over all iterables to build an array of jobs.
grid['counter'    ] = [0]*len(grid['iterables'])
grid['modelNumber'] = np.zeros(dimensions,dtype=int)
jobNumber = -1
jobCount  =  1
for i in range(len(grid['counter'])):
    jobCount *= len(grid[grid['iterables'][i]])

if args.reprocess:
    # Reprocess output files. This can be useful if some previous processing of Cloudy output files failed (we often have tens of
    # thousands of these so some intermittment failures can occur).
    reprocess(grid,args)
else:
    while True:
        jobNumber += 1
        indices    = ()
        for i in range(len(grid['iterables'])):
            indices = indices+(grid['counter'][i],)
        grid['modelNumber'][indices] = jobNumber
        if args.model is None or args.model == jobNumber:
            if jobNumber % 100 == 0 or args.model is not None:
                print("Generating model "+str(jobNumber)+" of "+str(jobCount))
            generateJob(grid,args)
        for i in range(len(grid['counter'])):
            grid['counter'][i] += 1
            if grid['counter'][i] == dimensions[i]:
                grid['counter'][i] = 0
            else:
                break
        if all(x == 0 for x in grid['counter']):
            break
    # Exit here if we are to only generate models.
    if args.generateOnly:
        sys.exit()
    # Get a job manager to allow us to submit and monitor jobs.
    manager = queueManager.factory(args)
    # Launch all jobs.
    if 'jobs' in grid:
        manager.submitJobs(grid['jobs'])
    
# Validate results.
validate(grid,args)

# Output results.
output(grid,args)

