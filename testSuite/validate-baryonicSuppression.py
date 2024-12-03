#!/usr/bin/env python3
import sys
import os
import subprocess
import numpy as np
import h5py
import json
import codecs
from git import Repo

# Run models to validate suppression of the halo mass function by baryons.
# Andrew Benson (01-April-2024)

# Define the target data for the mass function suppression. This was read from Figure 2 of Zheng et al. (2024;
# https://ui.adsabs.harvard.edu/abs/2024arXiv240317044Z).
#
# * 'withBaryons' corresponds to the "RI" model of Zheng et al.
# * 'withBaryons_noReionization' corresponds to the "NR" model of Zheng et al.
#
# Array indices correspond to redshift:
#
# * 4 ==> z = 9.27
# * 3 ==> z = 5.72
# * 2 ==> z = 3.06
# * 1 ==> z = 0.00
redshifts = ( "", "9.27", "5.72", "3.06", "0.00" )
target = {"withBaryons": [np.array([])], "withBaryons_noReionization": [np.array([])]}
target['withBaryons'              ].append(
                                           np.array(
                                                    [
                                                     0.6756756756756757,
                                                     0.6540540540540540,
                                                     0.7261261261261260,
                                                     1.0000000000000000,
                                                     1.0072072072072071,
                                                     1.0000000000000000,
                                                     1.0072072072072071,
                                                     1.0000000000000000
                                                    ]
                                                   )
                                         )
target['withBaryons_noReionization'].append(
                                           np.array(
                                                    [
                                                     0.6756756756756757,
                                                     0.6612612612612612,
                                                     0.7333333333333334,
                                                     1.0000000000000000,
                                                     0.9927927927927928,
                                                     1.0000000000000000,
                                                     1.0000000000000000,
                                                     1.0072072072072071
                                                    ]
                                                   )
                                         )
target['withBaryons'              ].append(
                                           np.array(
                                                    [
                                                     0.6685714285714284,
                                                     0.7028571428571427,
                                                     0.6876190476190476,
                                                     0.7371428571428571,
                                                     0.6609523809523808,
                                                     0.9885714285714284,
                                                     0.9885714285714284,
                                                     1.0000000000000000
                                                    ]
                                                   )
                                         )
target['withBaryons_noReionization'].append(
                                           np.array(
                                                    [
                                                     0.6723809523809523,
                                                     0.7104761904761904,
                                                     0.8247619047619046,
                                                     0.9314285714285714,
                                                     0.9466666666666665,
                                                     0.9999999999999999,
                                                     0.9885714285714284,
                                                     1.0000000000000000
                                                    ]
                                                   )
                                         )
target['withBaryons'              ].append(
                                           np.array(
                                                    [
                                                     0.6750000000000000,
                                                     0.7321428571428572,
                                                     0.6821428571428573,
                                                     0.7392857142857143,
                                                     0.5607142857142857,
                                                     0.7464285714285716,
                                                     0.9964285714285716,
                                                     1.0000000000000000
                                                    ]
                                                   )
                                         )
target['withBaryons_noReionization'].append(
                                           np.array(
                                                    [
                                                     0.7321428571428572,
                                                     0.7678571428571429,
                                                     0.8250000000000001,
                                                     0.9464285714285715,
                                                     0.8964285714285715,
                                                     0.9892857142857143,
                                                     0.9964285714285716,
                                                     1.0000000000000000
                                                    ]
                                                   )
                                         )
target['withBaryons'              ].append(
                                           np.array(
                                                    [
                                                     0.7115384615384613,
                                                     0.6961538461538460,
                                                     0.7576923076923074,
                                                     0.6961538461538460,
                                                     0.7730769230769230,
                                                     0.7346153846153844,
                                                     0.6423076923076921,
                                                     0.9500000000000000
                                                    ]
                                                   )
                                         )
target['withBaryons_noReionization'].append(
                                           np.array(
                                                    [
                                                     0.7346153846153844,
                                                     0.7423076923076921,
                                                     0.8269230769230768,
                                                     0.9038461538461537,
                                                     0.9423076923076922,
                                                     1.0346153846153845,
                                                     0.9038461538461537,
                                                     0.9807692307692306
                                                    ]
                                                   )
                                         )

# Define χ² targets for each dataset.
chiSquaredTarget = {"withBaryons": np.array([0.0,6.0,4.5,5.0,3.0]), "withBaryons_noReionization": np.array([0.0,7.0,3.0,5.0,4.0])}

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the validate model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/validate_baryonicSuppression_IGM_evolution.xml",shell=True)
if status.returncode != 0:
    print("FAILED: baryonic suppression (IGM evolution) validation model failed to run")
    sys.exit()

# Read data and repackage into a file suitable for re-reading by other models.
with h5py.File("outputs/validate_baryonicSuppression_IGM_evolution.hdf5","r") as model:
    igmFile       = h5py.File("outputs/validate_baryonicSuppression_IGM.hdf5","w")
    igmProperties = model['igmProperties']
    data          = {}
    for propertyName in ( 'redshift', 'temperature', 'densityHydrogen1', 'densityHydrogen2', 'densityHelium1', 'densityHelium2', 'densityHelium3' ):
        data[propertyName] = igmProperties[propertyName][:]
    data['hIonizedFraction' ] =  data['densityHydrogen2']                                                   /(data['densityHydrogen1']+data['densityHydrogen2']                       )
    data['heIonizedFraction'] = (data['densityHelium2'  ]+data['densityHelium3']                           )/(data['densityHelium1'  ]+data['densityHelium2'  ]+data['densityHelium3'])
    data['electronFraction' ] = (data['densityHydrogen2']+data['densityHelium2']+2.0*data['densityHelium3'])/(data['densityHydrogen1']+data['densityHydrogen2']                       )
    igmFile['redshift'         ] = data['redshift'         ]
    igmFile['matterTemperature'] = data['temperature'      ]
    igmFile['hIonizedFraction' ] = data['hIonizedFraction' ]
    igmFile['heIonizedFraction'] = data['heIonizedFraction']
    igmFile['electronFraction' ] = data['electronFraction' ]
    igmFile.attrs['extrapolationAllowed'] = 1
    igmFile.attrs['fileFormat'          ] = 1
    igmFile.close()

# Establish bins in halo mass.
massHaloLogarithmicBins = np.linspace(4.0,7.5,8)
haloMassFunction        = {}

# Run models with and without baryonic suppression.
for suffix in "withoutBaryons", "withBaryons", "withBaryons_noReionization":
    print("Running model: '"+suffix+"'")
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/validate_baryonicSuppression_evolve_"+suffix+".xml",shell=True)
    if status.returncode != 0:
        print("FAILED: baryonic suppression (evolve: '"+suffix+"') validation model failed to run")
        sys.exit()
    model              = h5py.File("outputs/validate_baryonicSuppression_evolve_"+suffix+".hdf5","r")
    cosmology          = model['Parameters/cosmologyParameters']
    OmegaMatter        = cosmology.attrs['OmegaMatter']
    OmegaBaryon        = cosmology.attrs['OmegaBaryon']
    fractionDarkMatter = (OmegaMatter-OmegaBaryon)/OmegaMatter
    haloMassFunction[suffix] = [None] * 5
    # Iterate over outputs.
    for outputIndex in range(1,5):
        print("\tProcessing output: "+str(outputIndex))
        # Read all required data.
        output = model['Outputs/Output'+str(outputIndex)]
        nodes  = output['nodeData']
        data = {}
        for propertyName in 'nodeIsIsolated', 'mergerTreeIndex', 'nodeIndex', 'parentIndex', 'hotHaloMass', 'diskMassGas', 'massHaloEnclosedCurrent', 'basicMass', 'nodeIsIsolated', 'mergerTreeWeight':
            data[propertyName] = nodes[propertyName][:]
        # Identify isolated and subhalos.
        isolated = np.nonzero(data['nodeIsIsolated'] == 1)[0]
        subhalo  = np.nonzero(data['nodeIsIsolated'] == 0)[0]
        # Build an index mapping each halo to its host halo. We take advantage here of the depth-first ordering of the outputs.
        index = np.zeros(len(data['nodeIsIsolated']),dtype='int32')
        for i in range(len(isolated)):
            indexStart                 = 0 if i == 0 else isolated[i-1]+1
            indexEnd                   =                  isolated[i  ]+1
            index[indexStart:indexEnd] =                  isolated[i  ]
        # Accumulate masses of isolated halos.
        if suffix == "withoutBaryons":
            # In models without baryons, the halo mass is just the dark matter mass.
            massHalo            = +data['massHaloEnclosedCurrent']
        else:
            # In models with baryons we assume that the hot gas of the host is distributed as the dark matter, so compute a
            # correction factor to account for the differing virial density contrast definitions in Galacticus and Zheng et
            # al. Gas in the disk is assumed to always be included within the virial radius.
            correctionFactor = +data['massHaloEnclosedCurrent'] \
                               /data['basicMass'              ]
            massHalo         = +data['massHaloEnclosedCurrent']*fractionDarkMatter \
                               +data['hotHaloMass'            ]*correctionFactor   \
                               +data['diskMassGas'            ]
            # Accumulate masses of subhalos halos. We assume that these are also distributed as the dark matter of the host and so
            # apply the same correction factor.
            massHalo[index][subhalo] += data['hotHaloMass'][subhalo]*correctionFactor[index][subhalo]
            massHalo[index][subhalo] += data['diskMassGas'][subhalo]*correctionFactor[index][subhalo]
        # Construct final quantities needed for the mass function.
        weight              =          data    ['mergerTreeWeight'][isolated]
        massHaloLogarithmic = np.log10(massHalo                    [isolated])
        # Construct the mass function.
        massHaloLogarithmicBinWidth  = massHaloLogarithmicBins[1]-massHaloLogarithmicBins[0]
        massHaloLogarithmicBinsEdges = np.append(massHaloLogarithmicBins-0.5*massHaloLogarithmicBinWidth,massHaloLogarithmicBins[-1]+0.5*massHaloLogarithmicBinWidth)
        haloMassFunction[suffix][outputIndex]                      = {}
        haloMassFunction[suffix][outputIndex]['massFunction'     ] =         np.histogram(massHaloLogarithmic,massHaloLogarithmicBinsEdges,weights=weight   )[0]
        haloMassFunction[suffix][outputIndex]['massFunctionError'] = np.sqrt(np.histogram(massHaloLogarithmic,massHaloLogarithmicBinsEdges,weights=weight**2)[0])

# Compute ratios of mass functions with the dark matter only model mass function.
output     = {"model": {}, "target": {}}
chiSquared = {}
failed     = False
for suffix in "withBaryons", "withBaryons_noReionization":
    chiSquared          [suffix] = [None] * 5
    output    ["model" ][suffix] = [None] * 5
    output    ["target"][suffix] = [None] * 5
    for outputIndex in range(1,5):
        chiSquared          [suffix][outputIndex] = {}
        output    ["model" ][suffix][outputIndex] = {}
        output    ["target"][suffix][outputIndex] = {}
        haloMassFunction[suffix][outputIndex]['ratio'     ] = haloMassFunction[suffix][outputIndex]['massFunction']/haloMassFunction['withoutBaryons'][outputIndex]['massFunction']
        haloMassFunction[suffix][outputIndex]['ratioError'] = +np.sqrt(
                     +(
                       +haloMassFunction[suffix          ][outputIndex]['massFunctionError']
                       /haloMassFunction[suffix          ][outputIndex]['massFunction'     ]
                      )**2 \
                     +(
                       +haloMassFunction['withoutBaryons'][outputIndex]['massFunctionError']
                       /haloMassFunction['withoutBaryons'][outputIndex]['massFunction'     ]
                      )**2
                    ) \
            *           haloMassFunction[suffix          ][outputIndex]['ratio'            ]
        chiSquared[suffix][outputIndex] = +np.sum(
                    +(
                      +haloMassFunction[suffix][outputIndex]['ratio'     ]
                      -target          [suffix][outputIndex]
                     )                                                    **2
                    /  haloMassFunction[suffix][outputIndex]['ratioError']**2
                   ) \
            /   len(   target          [suffix][outputIndex])
        output['model'][suffix][outputIndex]['ratio'     ] = list(haloMassFunction[suffix][outputIndex]['ratio'     ])
        output['model'][suffix][outputIndex]['ratioError'] = list(haloMassFunction[suffix][outputIndex]['ratioError'])
        # Report.
        status     = "SUCCESS" if chiSquared[suffix][outputIndex] < chiSquaredTarget[suffix][outputIndex] else "FAILED"
        if status == "FAILED":
            failed = True
        inequality = "<" if chiSquared[suffix][outputIndex] < chiSquaredTarget[suffix][outputIndex] else "≥"
        padding = " " * (len("withBaryons_noReionization")-len(suffix))
        print(f"{status}: model '{suffix}'{padding} at z={redshifts[outputIndex]} validation (χ² = {chiSquared[suffix][outputIndex]:5.3f} {inequality} {chiSquaredTarget[suffix][outputIndex]:5.3f})")

# Interface with git.
repo         = Repo(os.environ['GALACTICUS_EXEC_PATH'])
actor        = repo.head.commit.author
lastRevision = repo.head.object.hexsha
authorName   = actor.name
authorEmail  = actor.email
authorDate   = str(repo.head.commit.committed_datetime)
message      = repo.head.commit.message

# Generate content for the validation metrics page.
output['repoUrl'      ] = "https://github.com/galacticusorg/galacticus";
output['parameterFile'] = "testSuite/parameters/validate_baryonicSuppression_evolve_withBaryons.xml";
output['commit'       ] = {
    "author":
    {
        "name" : authorName,
        "email": authorEmail
    },
    "id"       : lastRevision,
    "message"  : message,
    "timestamp": authorDate,
    "url"      : "https://github.com/galacticusorg/galacticus/commit/"+lastRevision
}
for suffix in "withBaryons", "withBaryons_noReionization":
    for outputIndex in range(1,5):
        output['target'][suffix][outputIndex] = list(target[suffix][outputIndex])
output['redshift'] = redshifts
output['massHalo'] = list(massHaloLogarithmicBins)
f = codecs.open("outputs/results_baryonicSuppression.json", "w", "utf-8")
f.write("window.BARYONICSUPPRESSION_DATA = ")
f.write(json.dumps(output,indent=4,ensure_ascii=False))
f.close()
if failed:
    print("model failed - results were:\n")
    with open("outputs/results_baryonicSuppression.json", "r") as file:
        for line in file:
            print(line.replace("\n",""))
        
