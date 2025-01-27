#!/usr/bin/env python3
import sys
import os
import subprocess
import h5py
import numpy as np
import re
import xml.etree.ElementTree as ET
import json
import codecs
from git import Repo

# Run models to validate a dark matter only subhalo evolution model for a decaying dark matter scenario.
# Andrew Benson (13-January-2024)

# Create output path.
try:
    os.mkdir("outputs/validateDecayingDarkMatter")
except FileExistsError:
    pass

# Define required constants.
kilo         = 1.000e3
speedOfLight = 2.998e5

# Interface with git.
repo         = Repo(os.environ['GALACTICUS_EXEC_PATH'])
actor        = repo.head.commit.author
lastRevision = repo.head.object.hexsha
authorName   = actor.name
authorEmail  = actor.email
authorDate   = str(repo.head.commit.committed_datetime)
message      = repo.head.commit.message

# Parse the base parameter file.
ET.register_namespace("xi","http://www.w3.org/2001/XInclude")
parameterTree     = ET.parse('parameters/validate_darkMatterOnlySubHalos_decayingDarkMatter.xml')
parameters        = parameterTree.getroot()
nodeLifetime      = parameters.find('darkMatterParticle/lifetime'     )
nodeMassSplitting = parameters.find('darkMatterParticle/massSplitting')
nodeFileName      = parameters.find('outputFileName'                  )

# Open the target dataset file and iterate over models.
dataTarget = h5py.File('data/decayingDarkMatterSubhalosMau2022.hdf5')
for model in dataTarget.keys():
    # Extract combination of lifetime and kick velocity.
    matched      = re.search("lifetime:(\d+)_velocityKick:(\d+)",model)
    lifetime     = float(matched.group(1))
    velocityKick = float(matched.group(2))
    massSplit    = velocityKick/speedOfLight
    # Extract mass and radius CDFs.
    masses           = dataTarget[model]['mass'     ][:]
    radii            = dataTarget[model]['radius'   ][:]
    massCDFTarget    = dataTarget[model]['massCDF'  ][:].astype(np.float64)
    radiusCDFTarget  = dataTarget[model]['radiusCDF'][:].astype(np.float64)
    massCDFNonZero   = massCDFTarget   > 0
    radiusCDFNonZero = radiusCDFTarget > 0
    massesCentral    = np.sqrt(masses[0:-1]*masses[1:])
    radiiCentral     = np.sqrt(radii [0:-1]*radii [1:])
    # Set file names.
    parameterFileName  = f'outputs/validateDecayingDarkMatter/lifetime:{lifetime}_velocityKick:{velocityKick}.xml'
    galacticusFileName = f'outputs/validateDecayingDarkMatter/lifetime:{lifetime}_velocityKick:{velocityKick}.hdf5'
    # Update the parameters.
    nodeLifetime     .attrib['value'] = str(lifetime )
    nodeMassSplitting.attrib['value'] = str(massSplit)
    nodeFileName     .attrib['value'] = 'testSuite/'+galacticusFileName
    # Output a new parameter file.
    with open(parameterFileName, 'wb') as parameterFile:
        parameterTree.write(parameterFile)
    # Run the validation model.
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/"+parameterFileName,shell=True)
    if status.returncode != 0:
        print("FAILED: dark matter-only decaying dark matter subhalos validation model failed to run")
        sys.exit()
    # Initialize likelihoods and results.
    likelihoods         = []
    results             = []
    # Extract data from the Galacticus model.
    galacticus          = h5py.File(galacticusFileName,"r")
    nodes               = galacticus['Outputs/Output1/nodeData']
    isCentral           = nodes['nodeIsIsolated'                     ][:]
    massBound           = nodes['satelliteBoundMass'                 ][:]
    radiusVirial        = nodes['darkMatterOnlyRadiusVirial'         ][:]
    positionX           = nodes['positionOrbitalX'                   ][:]
    positionY           = nodes['positionOrbitalY'                   ][:]
    positionZ           = nodes['positionOrbitalZ'                   ][:]
    velocityMaximum     = nodes['darkMatterProfileDMOvelocityMaximum'][:]
    velocityPeak        = nodes['darkMatterProfileDMOvelocityPeak'   ][:]
    weight              = nodes['nodeSubsamplingWeight'              ][:]
    # Compute orbital radii (and convert to kpc).
    positionR           = np.sqrt(positionX**2+positionY**2+positionZ**2)*kilo
    # Count centrals.
    centrals            = (isCentral == 1)
    countCentrals       = np.count_nonzero(centrals)
    # Find virial radius of centrals.
    radiusVirialCentral = radiusVirial[centrals][0]*kilo
    # Select subhalos.
    subhalosMass        = (isCentral == 0) & (velocityMaximum > 9.0) & (velocityPeak > 10.0) & (positionR < radiusVirialCentral) & (massBound > 2.0e7)
    subhalosRadius      = (isCentral == 0) & (velocityMaximum > 9.0) & (velocityPeak > 10.0) & (positionR < radiusVirialCentral) & (massBound > 1.0e8)
    # Construct CDFs.
    massHistogram  , _ = np.histogram(massBound[subhalosMass  ],weights=weight[subhalosMass  ],bins=masses)
    massCDF            = np.flip(np.cumsum(np.flip(massHistogram  )))/countCentrals
    radiusHistogram, _ = np.histogram(positionR[subhalosRadius],weights=weight[subhalosRadius],bins=radii )
    radiusCDF          =         np.cumsum(        radiusHistogram)  /countCentrals    
    # Evaluate the likelihood.
    logLikelihoodMass   = -0.5*np.sum((  massCDF[  massCDFNonZero]-  massCDFTarget[  massCDFNonZero])**2/  massCDFTarget[  massCDFNonZero])
    logLikelihoodRadius = -0.5*np.sum((radiusCDF[radiusCDFNonZero]-radiusCDFTarget[radiusCDFNonZero])**2/radiusCDFTarget[radiusCDFNonZero])
    # Store CMF data.
    resultMass   = {
        "attributes": {
                "name": "massCumulativeDistribution",
                "xAxisLabel": "M_subhalo [M☉]",
                "yAxisLabel": "N(<M_subhalo)",
                "xAxisIsLog": "1",
                "yAxisIsLog": "1",
                "targetLabel": "Mau et al. (DES), 2022",
                "comment": "Analysis of subhalo bound mass cumulative distribution function",
                "description": f'Subhalo bound mass cumulative distribution function (τ={lifetime} Gyr; vₖ={velocityKick} km/s)',
                "logLikelihood": str(np.abs(logLikelihoodMass  )),
                "type": "function1D"
        },
        "data": {
            "xDataset": list(massesCentral),
            "yDataset": list(massCDF),
            "yDatasetTarget": list(massCDFTarget),
            "yErrorTarget": list(np.sqrt(massCDFTarget))
        }
    }
    resultRadius = {
        "attributes": {
                "name": "radialCumulativeDistribution",
                "xAxisLabel": "r_orbit [kpc]",
                "yAxisLabel": "N(<r_orbit)",
                "xAxisIsLog": "1",
                "yAxisIsLog": "1",
                "targetLabel": "Mau et al. (DES), 2022",
                "comment": "Analysis of subhalo orbital radius cumulative distribution function",
                "description": f'Subhalo orbital radius cumulative distribution function (τ={lifetime} Gyr; vₖ={velocityKick} km/s)',
                "logLikelihood": str(np.abs(logLikelihoodRadius)),
                "type": "function1D"
        },
        "data": {
            "xDataset": list(radiiCentral),
            "yDataset": list(radiusCDF),
            "yDatasetTarget": list(radiusCDFTarget),
            "yErrorTarget": list(np.sqrt(radiusCDFTarget))
        }
    }
    results    .append(resultMass  )
    results    .append(resultRadius)
    likelihoodMass   = {
        "name" : f' (τ={lifetime} Gyr; vₖ={velocityKick} km/s) - Likelihood - mass function',
        "unit" : "-logℒ"                                                                   ,
        "value": str(np.abs(logLikelihoodMass  ))
    }
    likelihoodRadius = {
        "name" : f' (τ={lifetime} Gyr; vₖ={velocityKick} km/s) - Likelihood - orbital radius function',
        "unit" : "-logℒ"                                                                             ,
        "value": str(np.abs(logLikelihoodRadius))
    }
    likelihoods.append(likelihoodMass  )
    likelihoods.append(likelihoodRadius)
    # Generate content for the validation metrics page.
    output = {
        "repoUrl"      : "https://github.com/galacticusorg/galacticus",
        "parameterFile": "testSuite/parameters/validate_darkMatterOnlySubHalos_decayingDarkMatter.xml",
        "commit"       :
        {
            "author":
            {
                "name" : authorName,
                "email": authorEmail
            },
            "id"       : lastRevision,
            "message"  : message,
            "timestamp": authorDate,
            "url"      : "https://github.com/galacticusorg/galacticus/commit/"+lastRevision
        },
        "results": results
    }
    # Output results.
    suffix=f'darkMatterOnlySubhalos_decayingDarkMatter_lifetime{lifetime}_velocityKick{velocityKick}'
    f = codecs.open("outputs/results_"+suffix+".json", "w", "utf-8")
    f.write("window.ANALYSES_DATA = ")
    f.write(json.dumps(output,indent=4,ensure_ascii=False))
    f.close()
    # Write benchmark results.
    f = codecs.open("outputs/validate_"+suffix+".json", "w", "utf-8")
    f.write(json.dumps(likelihoods,indent=4,ensure_ascii=False))
    f.close()

print("SUCCESS: decaying dark matter dark matter-only subhalos validation model")
