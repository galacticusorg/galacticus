#!/usr/bin/env python3
import sys
import os
import subprocess
import h5py
import numpy as np
import argparse
import shutil
import re
import lxml.etree as ET
import json
import codecs
import matplotlib.pyplot as plt
import queueManager
import time
from git import Repo
from scipy.optimize import minimize
from num2tex import num2tex

# Run models to measure convergence in halo concentrations determined from merger tree structure with mass resolution.
# Andrew Benson (15-May-2025)

# This validation will:
#  * generate a set of merger trees spanning a range of resolutions and numerical tolerances;
#  * compute concentrations for the halos using the model of Johnson, ......;
#  * create plots showing this concentration distribution;
#  * output JSON files containing the results for inclusion into a web page.
#
# The validation is typically run as:
#
#  ./testSuite/validate_concentrationConvergence.py --outputPath /scratch/concentrationConvergence --jobMaximum 10 --waitOnSubmit 1 --waitOnActive 30
#
# where the first option controls the directory to which results will be output and the other (optional) arguments are passed to
# the queue manager.
#
# Once complete (note that these high-resolution models can take a very long time to finish - ~2 weeks) plots are created and
# output to the `outputPath` showing the convergence with respect to each numerical parameter considered. Additionally, a file
# `results.json` is output to `outputPath`. This contains the raw data on all convergence tests and should be copied to
# `dev/valid/concentrationConvergence/` in the `gh-pages` branch and committed there. It will then appear as a page of convergence
# plots at https://galacticusorg.github.io/galacticus/dev/valid/concentrationConvergence.

# Convert command line arguments to ints.
def restricted_int(x):
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not an integer literal" % (x,))
    return x

# Convert command line arguments to reals.
def restricted_real(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a float literal" % (x,))
    return x

# Convert command line arguments to bools.
def restricted_bool(x):
    if x != "true" and x != "false":
        raise argparse.ArgumentTypeError("%r not a boolean literal"  % (x,))
    return x

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='validate_concentrationConvergence.py',description='Compute convergence in the progenitor mass function.')
parser.add_argument('--outputPath'                 ,action='store'     ,default="."                             ,help='the path to which data should be written'                                                                           )
parser.add_argument('--partition'                  ,action='store'                                              ,help='the partition to which to submit jobs'                                                                              )
parser.add_argument('--jobMaximum'                 ,action='store'                     ,type=restricted_int     ,help='the maximum number of active jobs to allow'                                                                         )
parser.add_argument('--waitOnSubmit'               ,action='store'                     ,type=restricted_int     ,help='the time (in seconds) to wait after submitting each job'                                                            )
parser.add_argument('--waitOnActive'               ,action='store'                     ,type=restricted_int     ,help='the time (in seconds) to wait after polling active jobs'                                                            )
parser.add_argument('--energyBoost'                ,action='store'                                              ,help='the energy boost parameter in the Johnson et al. (2021) model'                                                      )
parser.add_argument('--massExponent'               ,action='store'                                              ,help='the mass exponent parameter in the Johnson et al. (2021) model'                                                     )
parser.add_argument('--peakHeightExponent'         ,action='store'                                              ,help='the peak height exponent parameter in the Johnson et al. (2021) model'                                              )
parser.add_argument('--scatterExcess'              ,action='store'                                              ,help='the excess scatter parameter in the Johnson et al. (2021) model'                                                    )
parser.add_argument('--unresolvedEnergy'           ,action='store'                                              ,help='the multiplier of unresolved energy in the Johnson et al. (2021) model'                                             )
parser.add_argument('--factorMassResolution'       ,action='store'                                              ,help='the factor above the mass resolution at which to begin applying the Johnson et al. (2021) model'                    )
parser.add_argument('--factorMassTrust'            ,action='store'     ,default=100.0      ,type=restricted_real,help='the factor above the mass scale at which the Johnson et al. (2021) model is applied that the results can be trusted')
parser.add_argument('--countSampleEnergyUnresolved',action='store'                                              ,help='the number of samples in Monte Carlo estimates of unresolved energy in the Johnson et al. (2021) model'             )
parser.add_argument('--acceptUnboundOrbits'        ,action='store'     ,default="false"    ,type=restricted_bool,help='if true, allow unbound orbits in the Johnson et al. (2021) model'                                                   )
parser.add_argument('--includeUnresolvedVariance'  ,action='store'     ,default="false"    ,type=restricted_bool,help='if true, include variance in the energy of unresolved orbits in the Johnson et al. (2021) model'                    )
parser.add_argument('--orbitModel'                 ,action='store'     ,default="unchanged"                     ,help='the orbit model to use in the Johnson et al. (2021) model (or "unchanged" to leave it unchanged)'                   )
parser.add_argument('--mainBranchOnly'             ,action='store_true'                                         ,help='if present apply to model to the main branch only'                                                                  )
parser.add_argument('--removeConcentrationLimits'  ,action='store_true'                                         ,help='if present remove any imposed limits on concentration'                                                              )

args = parser.parse_args()

# Find the "time" command.
timeCommand = shutil.which("time")

def extractResults(job):
    # Check for a failed model.
    if "exitStatus" in job and job['exitStatus'] > 0:
        print("FAIL: model '"+job['label']+"' failed to run")
        sys.exit()
    # Wait for the model to appear.
    while not os.path.exists(job['galacticusFileName']):
        time.sleep(1)
    # Call the relevant extractor.
    if job['model']['fullOutput']:
       extractFullOutput(job)
    else:
       extractStandard  (job)

def extractFullOutput(job):
    # Extract and analyze data from a full output model.
    galacticus                     = h5py.File(job['galacticusFileName'],"r")
    outputs                        = galacticus['Outputs']
    factorMassResolution           = galacticus['Parameters/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'].attrs['factorMassResolution']
    nodeMasses                     = np.array([])
    nodeRadii                      = np.array([])
    nodeExpansionFactors           = np.array([])
    nodePeakHeights                = np.array([])
    nodeMassesMainBranch           = np.array([])
    nodeRadiiMainBranch            = np.array([])
    nodeExpansionFactorsMainBranch = np.array([])
    nodePeakHeightsMainBranch      = np.array([])
    nodeTreeIndexMainBranch        = np.array([])
    nodeMassMinimum                = job['model']['massTrees']*job['model']['massResolutionFractional']*factorMassResolution*args.factorMassTrust
    indexMatchMaximum              = 0
    nodePeakHeightFinal            = None
    for outputName, output in outputs.items():
        match = re.match(r'Output(\d+)',outputName)
        if match:
            indexMatch = int(match.group(1))
        else:
            continue
        expansionFactor                = output.attrs['outputExpansionFactor']
        nodes                          = output['nodeData'              ]
        nodeIsIsolated                 = nodes ['nodeIsIsolated'        ][:]
        nodeIsOnMainBranch             = nodes ['nodeIsOnMainBranch'    ][:]
        nodeMass                       = nodes ['basicMass'             ][:]
        nodeRadius                     = nodes ['darkMatterProfileScale'][:]
        nodePeakHeight                 = nodes ['haloPeakHeightNu'      ][:]
        nodeTreeIndex                  = nodes ['mergerTreeIndex'       ][:]
        selection                      = (nodeIsIsolated     == 1) & (nodeMass > nodeMassMinimum)
        mainBranch                     = (nodeIsOnMainBranch == 1) & (nodeMass > nodeMassMinimum)
        nodeMasses                     = np.append(nodeMasses                    ,nodeMass      [selection ])
        nodeMassesMainBranch           = np.append(nodeMassesMainBranch          ,nodeMass      [mainBranch])
        nodeRadii                      = np.append(nodeRadii                     ,nodeRadius    [selection ])
        nodeRadiiMainBranch            = np.append(nodeRadiiMainBranch           ,nodeRadius    [mainBranch])
        nodePeakHeights                = np.append(nodePeakHeights               ,nodePeakHeight[selection ])
        nodePeakHeightsMainBranch      = np.append(nodePeakHeightsMainBranch     ,nodePeakHeight[mainBranch])
        nodeTreeIndexMainBranch        = np.append(nodeTreeIndexMainBranch       ,nodeTreeIndex [mainBranch])
        nodeExpansionFactors           = np.append(nodeExpansionFactors          ,np.ones(np.count_nonzero(selection ))*expansionFactor)
        nodeExpansionFactorsMainBranch = np.append(nodeExpansionFactorsMainBranch,np.ones(np.count_nonzero(mainBranch))*expansionFactor)
        if indexMatch > indexMatchMaximum:
            indexMatchMaximum   = indexMatch
            nodePeakHeightFinal = nodePeakHeight[0]
        # Write the data to a CSV file for easy visualization.
        for i in range(len(nodeMass)):
            if selection[i]:
                fileData.write(str(np.log10(nodeMass[i]))+" , "+str(np.log10(expansionFactor))+" , "+str(np.log10(nodePeakHeight[i]))+" , "+str(np.log10(nodeRadius[i]))+"\n")
    # Store the results.
    job['results'][job['id']] = {
        "massResolutionFractional":       massResolutionFractional         ,
        "factorMassResolution":           factorMassResolution             ,
        "model":                          job['model']                     ,
        "isResolutionModel":              job['model']['isResolutionModel'],
        "fullOutput":                     job['model']['fullOutput'       ],
        "massTrees":                      job['model']['massTrees'        ],
        "nodeMasses":                     nodeMasses                       ,
        "nodeRadii":                      nodeRadii                        ,
        "nodePeakHeights":                nodePeakHeights                  ,
        "nodeExpansionFactors":           nodeExpansionFactors             ,
        "nodeMassesMainBranch":           nodeMassesMainBranch             ,
        "nodeRadiiMainBranch":            nodeRadiiMainBranch              ,
        "nodePeakHeightsMainBranch":      nodePeakHeightsMainBranch        ,
        "nodeExpansionFactorsMainBranch": nodeExpansionFactorsMainBranch   ,
        "nodePeakHeightFinal":            nodePeakHeightFinal              ,
        "nodeTreeIndexMainBranch":        nodeTreeIndexMainBranch
    }
    
def extractStandard(job):
    # Extract data from the Galacticus model.
    galacticus                = h5py.File(job['galacticusFileName'],"r")
    nodes                     = galacticus['Outputs/Output1/nodeData']
    job['results'][job['id']] = {}
    radiusScale               = nodes['darkMatterProfileScale'][:]
    weight                    = nodes['nodeSubsamplingWeight' ][:]
    massResolutionFractional  = galacticus['Parameters/mergerTreeMassResolution'].attrs['massResolutionFractional']
    # Extract parameters of the Johnson model.
    if job['model']['isResolutionModel']:
        energyBoost               = galacticus['Parameters/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'].attrs['energyBoost'         ]
        massExponent              = galacticus['Parameters/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'].attrs['massExponent'        ]
        peakHeightExponent        = galacticus['Parameters/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'].attrs['peakHeightExponent'  ]
        scatterExcess             = galacticus['Parameters/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'].attrs['scatterExcess'       ]
        unresolvedEnergy          = galacticus['Parameters/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'].attrs['unresolvedEnergy'    ]
        factorMassResolution      = galacticus['Parameters/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'].attrs['factorMassResolution']
        output['parametersJohnson'] = {
            "energyBoost":          energyBoost         ,
            "massExponent":         massExponent        ,
            "peakHeightExponent":   peakHeightExponent  ,
            "scatterExcess":        scatterExcess       ,
            "unresolvedEnergy":     unresolvedEnergy    ,
            "factorMassResolution": factorMassResolution
        }
    # Construct scale radius PDF.
    radiusScaleFunctionCount, _ = np.histogram(np.log10(radiusScale*kilo),bins=radiiScale,weights=weight)
    radiusScaleFunction         =         radiusScaleFunctionCount /radiiScaleDelta/float(countTrees)
    radiusScaleFunctionError    = np.sqrt(radiusScaleFunctionCount)/radiiScaleDelta/float(countTrees)
    # Estimate the mean and scatter.
    radiusScaleLogarithmicMean    =         np.average( np.log10(radiusScale*kilo)                               ,weights=weight)
    radiusScaleLogarithmicScatter = np.sqrt(np.average((np.log10(radiusScale*kilo)-radiusScaleLogarithmicMean)**2,weights=weight))
    # Store the results.
    job['results'][job['id']] = {
        "massResolutionFractional":      massResolutionFractional         ,
        "model":                         job['model']                     ,
        "isResolutionModel":             job['model']['isResolutionModel'],
        "fullOutput":                    job['model']['fullOutput'       ],
        "massTrees":                     job['model']['massTrees'        ],
        "redshift":                      job['model']['redshift'         ],
        "radiusScaleFunction":           radiusScaleFunction              ,
        "radiusScaleFunctionError":      radiusScaleFunctionError         ,
        "radiusScaleLogarithmicMean":    radiusScaleLogarithmicMean       ,
        "radiusScaleLogarithmicScatter": radiusScaleLogarithmicScatter    ,
    }

def sigmoid(y0,y1,xt,xw,x):
    # Sigmoid function.
    return y0+(y1-y0)/(1.0+np.exp(-(x-xt)/xw))
    
def radiusMeanPowerLaw(radiusLow,radiusHigh,radiusTransition,radiusWidth,massLow,massHigh,massTransition,massWidth,expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth,mass,expansionFactor,peakHeight):
    # Evaluate the scale radius model.
    radius                  = sigmoid(radiusLow         ,radiusHigh         ,radiusTransition         ,radiusWidth         ,peakHeight)
    exponentMass            = sigmoid(massLow           ,massHigh           ,massTransition           ,massWidth           ,peakHeight)
    exponentExpansionFactor = sigmoid(expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth,peakHeight)
    return radius*(mass/massReference)**exponentMass*expansionFactor**exponentExpansionFactor
    
def radiusScaleFit(parameters,results,computeResidual=False,dataFile=None):
    # Fitting function for the scale radius model.
    ## Extract the parameters of the model.
    radiusLow,radiusHigh,radiusTransition,radiusWidth,massLow,massHigh,massTransition,massWidth,expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth = parameters
    # Iterate over results, accumulating the fit metric.
    fitMetric = 0.0
    count     = 0
    for id, result in results.items():
        # Skip any non-full output model.
        if not result['fullOutput']:
            continue
        radiiModel  = radiusMeanPowerLaw(radiusLow,radiusHigh,radiusTransition,radiusWidth,massLow,massHigh,massTransition,massWidth,expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth,result['nodeMasses'],result['nodeExpansionFactors'],result['nodePeakHeights'])
        residuals   = np.log10(result['nodeRadii']/radiiModel)
        weight      = result['nodeMasses']/result['massTrees']
        count      += np.sum(             weight)
        fitMetric  += np.sum(residuals**2*weight)
        # Write data to file if requested.
        if dataFile is not None:
            for i in range(len(result['nodeMasses'])):
                dataFile.write(str(np.log10(result['nodeMasses'][i]))+" , "+str(np.log10(result['nodeExpansionFactors'][i]))+" , "+str(np.log10(result['nodePeakHeights'][i]))+" , "+str(np.log10(radiiModel[i]))+"\n")
    # Compute the rms residual if requested.
    residualRMS = np.sqrt(fitMetric/count)
    if computeResidual:
        return fitMetric, residualRMS
    else:
        return fitMetric

def correlationFit(parameters,correlationMatrix,massLogarithmic):
    # Fitting function for the scale radius correlation model.
    rateDecay, exponent  = parameters
    correlationMatrixFit = np.zeros(correlationMatrix.shape)
    for i in range(len(massLogarithmic)):
        for j in range(len(massLogarithmic)):
            correlationMatrixFit[i,j] = np.exp(-rateDecay*np.abs(massLogarithmic[i]-massLogarithmic[j])**exponent)
    fitMetric = np.sum((correlationMatrix-correlationMatrixFit)**2)
    return fitMetric
    
# Get a job manager to allow us to submit and monitor jobs.
manager = queueManager.factory(args)

# Create output path.
try:
    os.mkdir(args.outputPath)
except FileExistsError:
    pass

# Interface with git.
repo         = Repo(os.environ['GALACTICUS_EXEC_PATH'])
actor        = repo.head.commit.author
lastRevision = repo.head.object.hexsha
authorName   = actor.name
authorEmail  = actor.email
authorDate   = str(repo.head.commit.committed_datetime)
message      = repo.head.commit.message
output = {
    "repoUrl"      : "https://github.com/galacticusorg/galacticus",
    "commit"       : {
        "author":  {
            "name" : authorName,
            "email": authorEmail
        },
        "id"       : lastRevision,
        "message"  : message,
        "timestamp": authorDate,
        "url"      : "https://github.com/galacticusorg/galacticus/commit/"+lastRevision
    },
    "parameterFile": "testSuite/parameters/validate_concentrationConvergence.xml",
    "results": {}
}
    
# Set default model attributes.
defaults = {
    "mergerTreeBranchingProbability/accuracyFirstOrder":           0.1,
    "mergerTreeBuilder/accretionLimit":                            0.1,
    "mergerTreeBuilder/mergeProbability":                          0.1,
    "mergerTreeBuildController[@value='subsample']/massThreshold": 0.0,
    "mergerTreeBuilder/branchIntervalStep":                        "false"
}

# Initialize output files.
fileData = open(args.outputPath+'/powerLawModelTargetData.csv','w')
fileData.write("# All halos used in constraining the 'simple power law' model of halo scale radii.\n")
fileData.write("log10Mass, log10ExpansionFactor, log10PeakHeight, log10RadiusScale\n")

# Define the reference mass.
massReference = 1.0e12

# Conversion factors.
kilo          = 1.0e3

# Specify the number of trees to run.
countTrees    = 2400
output['countTrees'] = countTrees

# Pre-determined values for the parameters of the power-law and correlation models used to set the initial conditions for
# poorly-resolved halos.
radiusLow                 = +0.0154
radiusHigh                = +0.0962
radiusTransition          = +1.2137
radiusWidth               = +0.5482
massLow                   = +0.3895
massHigh                  = +0.2984
massTransition            = -0.2583
massWidth                 = +1.6605e1
expansionFactorLow        = -0.6977
expansionFactorHigh       = +0.7972
expansionFactorTransition = +0.5395
expansionFactorWidth      = +0.4282
modelScatter              = +0.1513
correlationRateDecay      = +3.8361
correlationExponent       = +1.6198
## Store the values used to the output data structure.
output['parametersInitialConditions'] = {
    "radiusLow"                : radiusLow                ,
    "radiusHigh"               : radiusHigh               ,
    "radiusTransition"         : radiusTransition         ,
    "radiusWidth"              : radiusWidth              ,
    "massLow"                  : massLow                  ,
    "massHigh"                 : massHigh                 ,
    "massTransition"           : massTransition           ,
    "massWidth"                : massWidth                ,
    "expansionFactorLow"       : expansionFactorLow       ,
    "expansionFactorHigh"      : expansionFactorHigh      ,
    "expansionFactorTransition": expansionFactorTransition,
    "expansionFactorWidth"     : expansionFactorWidth     ,
    "modelScatter"             : modelScatter             ,
    "correlationRateDecay"     : correlationRateDecay     ,
    "correlationExponent"      : correlationExponent
    }

# Specify tree masses.
massesTree    = ( 1.0e10, 1.0e12, 1.0e14 )

# Specify redshifts.
redshifts     = ( 0.0, 1.0, 3.0, 6.0 )

# Define radii for the scale radius function.
radiiScale        = np.linspace(-1.0,4.0,51)
radiiScaleDelta   =      radiiScale[1   ]-radiiScale[0 ]
radiiScaleCentral = 0.5*(radiiScale[0:-1]+radiiScale[1:])

# Initialize job and result arrays.
jobs    = []
results = {}

# Iterate over tree masses.
for massTrees in massesTree:

    # Iterate over redshifts.
    for redshift in redshifts:

        # Define models to run.
        models = [
            {
                "massResolutionFractional":                                    1.0e-2,
                "darkMatterProfileScaleRadius":                                "concentration"
            },
            {
                "massResolutionFractional":                                    1.0e-5
            },
            {
                "massResolutionFractional":                                    1.0e-4
            },
            {
                "massResolutionFractional":                                    1.0e-3
            },
            {
                "massResolutionFractional":                                    1.0e-2
            },
            {
                "massResolutionFractional":                                    1.0e-1
            },
            {
                "massResolutionFractional":                                    5.0e-1
            },
            {
                "massResolutionFractional":                                    1.0e-5,
                "fullOutput":                                                  True
            }
        ]
        # For a specific mass and redshift, add other convergence tests.
        if massTrees == 1.0e12 and redshift == 0.0:
            models.extend(
                [
                    {
                        "massResolutionFractional":                                    1.0e-5,
                        "mergerTreeBuilder/accretionLimit":                            0.01
                    },
                    {
                        "massResolutionFractional":                                    1.0e-5,
                        "mergerTreeBuilder/mergeProbability":                          0.01
                    },
                    {
                        "massResolutionFractional":                                    1.0e-5,
                        "mergerTreeBranchingProbability/accuracyFirstOrder":           0.01
                    },
                    {
                        "massResolutionFractional":                                    1.0e-5,
                        "mergerTreeBuilder/branchIntervalStep":                        "true"
                    },
                    {
                        "massResolutionFractional":                                    1.0e-5,
                        "mergerTreeBuildController[@value='subsample']/massThreshold": 1.0e-3*massTrees
                    }
                ]
            )
        # Iterate over models.
        massResolutionFractionalMinimum = min(    map(lambda x: x['massResolutionFractional'],models) )
        countResolutions                = len(set(map(lambda x: x['massResolutionFractional'],models)))
        for model in models:
            # Set default behavior.
            if not 'fullOutput' in model.keys():
                model['fullOutput'] = False
            # Construct the file of parameter changes for this model.
            parameters = ET.Element("changes")
            # Set the task.
            change                     = ET.SubElement(parameters,'change')
            change    .attrib['type' ] = "append"
            change    .attrib['path' ] = ""
            taskEvolve                 = ET.SubElement(change,"task")
            taskEvolve.attrib['value'] = "evolveForests"
            # Set a non-evolving evolver (unless full output is required).
            if not model['fullOutput']:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "replace"
                change.attrib['path' ] = "mergerTreeEvolver"
                mergerTreeEvolver      = ET.SubElement(change           ,"mergerTreeEvolver")
                pruneTree              = ET.SubElement(mergerTreeEvolver,"pruneTree"        )
                mergerTreeEvolver.attrib['value'] = "nonEvolving"
                pruneTree        .attrib['value'] = "true"
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "remove"
                change.attrib['path' ] = "nodeOperator/nodeOperator[@value='darkMatterProfileScaleInterpolate']"
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "remove"
            change.attrib['path' ] = "nodeOperator/nodeOperator[@value='haloAngularMomentumInterpolate']"
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "remove"
            change.attrib['path' ] = "nodeOperator/nodeOperator[@value='subsubhaloPromotion']"
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "remove"
            change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteOrbit']"
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "remove"
            change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteDynamicalFriction']"
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "remove"
            change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteTidalMassLoss']"
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "remove"
            change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteTidalHeating']"
            if not model['fullOutput']:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "remove"
                change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteMergingRadiusTrigger']"
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "remove"
                change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteDestructionMassThreshold']"
            else:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteMergingRadiusTrigger']/radiusVirialFraction"
                change.attrib['value'] = "1.0"
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteDestructionMassThreshold']/massDestructionAbsolute"
                change.attrib['value'] = "0.0"
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "nodeOperator/nodeOperator[@value='satelliteDestructionMassThreshold']/massDestructionMassInfallFraction"
                change.attrib['value'] = "1.0"
            # Set merger tree masses.
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "append"
            change.attrib['path' ] = ""
            mergerTreeBuildMasses  = ET.SubElement(change               ,"mergerTreeBuildMasses")
            massTree               = ET.SubElement(mergerTreeBuildMasses,"massTree"             )
            treeCount              = ET.SubElement(mergerTreeBuildMasses,"treeCount"            )
            mergerTreeBuildMasses.attrib['value'] = "fixedMass"
            massTree             .attrib['value'] = str(massTrees )
            if model['fullOutput']:
                treeCount        .attrib['value'] = "24"
            else:
                treeCount        .attrib['value'] = str(countTrees)
            # Set merger tree resolution.
            change                   = ET.SubElement(parameters,'change')
            change.attrib['type' ]   = "append"
            change.attrib['path' ]   = ""
            mergerTreeMassResolution = ET.SubElement(change                  ,"mergerTreeMassResolution")
            massResolutionMinimum    = ET.SubElement(mergerTreeMassResolution,"massResolutionMinimum"   )
            massResolutionFractional = ET.SubElement(mergerTreeMassResolution,"massResolutionFractional")
            mergerTreeMassResolution.attrib['value'] = "scaled"
            massResolutionMinimum   .attrib['value'] = "0.0"
            massResolutionFractional.attrib['value'] = str(model['massResolutionFractional'])
            # Set tree redshift.
            change                   = ET.SubElement(parameters,'change')
            change.attrib['type' ]   = "update"
            change.attrib['path' ]   = "outputTimes/redshifts"
            change.attrib['value']   = str(redshift)
            change                   = ET.SubElement(parameters,'change')
            change.attrib['type' ]   = "update"
            change.attrib['path' ]   = "mergerTreeConstructor/mergerTreeConstructor/redshiftBase"
            change.attrib['value']   = str(redshift)
            # Set output options.
            if model['fullOutput']:
                # Set many outputs.
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "replace"
                change.attrib['path' ] = "outputTimes"
                outputTimes            = ET.SubElement(change,"outputTimes")
                outputTimes.attrib['value'] = 'logarithmicSpacingInCriticalOverdensity'
                redshiftMinimum             = ET.SubElement(outputTimes,'redshiftMinimum'    )
                redshiftMaximum             = ET.SubElement(outputTimes,'redshiftMaximum'    )
                countTimes                  = ET.SubElement(outputTimes,'countTimes'         )
                criticalOverdensity         = ET.SubElement(outputTimes,'criticalOverdensity')
                redshiftMinimum.attrib['value'] = str(redshift)
                redshiftMaximum.attrib['value'] = '30.0'
                countTimes     .attrib['value'] = '30'
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type'  ] = "replaceWith"
                change.attrib['path'  ] = 'outputTimes/criticalOverdensity'
                change.attrib['target'] = 'criticalOverdensity/criticalOverdensity'
                # Include tree index, peak height and main branch status in output.
                change                           = ET.SubElement(parameters,'change')
                change.attrib['type' ]           = "append"
                change.attrib['path' ]           = "mergerTreeOutputter"
                nodePropertyExtractorMulti       = ET.SubElement(change                    ,"nodePropertyExtractor")
                nodePropertyExtractorIndices     = ET.SubElement(nodePropertyExtractorMulti,"nodePropertyExtractor")
                nodePropertyExtractorMainBranch  = ET.SubElement(nodePropertyExtractorMulti,"nodePropertyExtractor")
                nodePropertyExtractorPeakheight  = ET.SubElement(nodePropertyExtractorMulti,"nodePropertyExtractor")
                nodePropertyExtractorIndicesTree = ET.SubElement(nodePropertyExtractorMulti,"nodePropertyExtractor")
                nodePropertyExtractorMulti      .attrib['value'] = "multi"
                nodePropertyExtractorIndices    .attrib['value'] = "nodeIndices"
                nodePropertyExtractorMainBranch .attrib['value'] = "mainBranchStatus"
                nodePropertyExtractorPeakheight .attrib['value'] = "peakHeight"
                nodePropertyExtractorIndicesTree.attrib['value'] = "indicesTree"
            # Filter out subhalos.
            change                   = ET.SubElement(parameters,'change')
            change.attrib['type' ]   = "append"
            change.attrib['path' ]   = "mergerTreeOutputter"
            galacticFilter           = ET.SubElement(change,"galacticFilter")
            galacticFilter.attrib['value'] = "haloIsolated"
            # Remove limits on concentration
            if args.removeConcentrationLimits:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/concentrationMinimum"
                change.attrib['value'] = "1.0e-30"
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/concentrationMaximum"
                change.attrib['value'] = "1.0e+30"
            # Set parameters of the Johnson et al. (2021) model.
            if args.energyBoost                 is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/energyBoost"
                change.attrib['value'] = args.energyBoost
            if args.massExponent                is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/massExponent"
                change.attrib['value'] = args.massExponent
            if args.peakHeightExponent          is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/peakHeightExponent"
                change.attrib['value'] = args.peakHeightExponent
            if args.scatterExcess               is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/scatterExcess"
                change.attrib['value'] = args.scatterExcess
            if args.unresolvedEnergy            is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/unresolvedEnergy"
                change.attrib['value'] = args.unresolvedEnergy
            if args.factorMassResolution        is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/factorMassResolution"
                change.attrib['value'] = args.factorMassResolution
            if args.acceptUnboundOrbits        is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "componentSatellite/acceptUnboundOrbits"
                change.attrib['value'] = args.acceptUnboundOrbits
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/acceptUnboundOrbits"
                change.attrib['value'] = args.acceptUnboundOrbits
            if args.includeUnresolvedVariance is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/includeUnresolvedVariance"
                change.attrib['value'] = args.includeUnresolvedVariance
            if args.countSampleEnergyUnresolved is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/countSampleEnergyUnresolved"
                change.attrib['value'] = args.countSampleEnergyUnresolved
            if args.countSampleEnergyUnresolved is not None:
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = "darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/countSampleEnergyUnresolved"
                change.attrib['value'] = args.countSampleEnergyUnresolved
            if args.orbitModel                  is not None:
                if args.orbitModel == "unchanged":
                    # Nothing to do.
                    continue
                elif args.orbitModel == "li2020":
                    change                 = ET.SubElement(parameters,'change')
                    change.attrib['type' ] = "replace"
                    change.attrib['path' ] = "virialOrbit/virialOrbit"
                    virialOrbit = ET.SubElement(change,"virialOrbit")
                    virialOrbit.attrib['value'] = "li2020"
            # Set default attributes of the model.
            for attribute in defaults.keys():
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = attribute
                change.attrib['value'] = str(defaults[attribute])
            # Set updated attributes of the model.
            id                = ""
            separator         = ""
            isResolutionModel = True
            for attribute in model.keys():
                if attribute != "massResolutionFractional" and attribute != "fullOutput":
                    isResolutionModel = False
                if attribute in defaults.keys() or attribute == "massResolutionFractional":
                    if isinstance(model[attribute],float):
                        id += f'{separator}{attribute}{model[attribute]:.1e}'
                    else:
                        id += f'{separator}{attribute}{model[attribute].capitalize()}'
                    separator = "_"
                    if attribute != "massResolutionFractional":
                        change                 = ET.SubElement(parameters,'change')
                        change.attrib['type' ] = "update"
                        change.attrib['path' ] = attribute
                        change.attrib['value'] = str(model[attribute])
                elif attribute == "fullOutput":
                    id += f'{separator}{attribute}{model[attribute]}'
            # Switch to Diemer & Joyce model if required.
            if "darkMatterProfileScaleRadius" in model.keys() and model['darkMatterProfileScaleRadius'] == "concentration":
                id += f'{separator}darkMatterProfileScaleRadius{model["darkMatterProfileScaleRadius"].capitalize()}'
                change                          = ET.SubElement(parameters                     ,'change'                         )
                darkMatterProfileScaleRadiusTmp = ET.SubElement(change                         ,"darkMatterProfileScaleRadiusTmp")
                darkMatterProfileScaleRadius    = ET.SubElement(darkMatterProfileScaleRadiusTmp,"darkMatterProfileScaleRadius"   )
                change.attrib['type'  ] = "append"
                change.attrib['path'  ] = ''
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type'  ] = "replaceWith"
                change.attrib['path'  ] = 'darkMatterProfileScaleRadiusTmp/darkMatterProfileScaleRadius'
                change.attrib['target'] = 'darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type'  ] = "replaceWith"
                change.attrib['path'  ] = 'darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'
                change.attrib['target'] = 'darkMatterProfileScaleRadiusTmp/darkMatterProfileScaleRadius'
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type'  ] = "remove"
                change.attrib['path'  ] = 'darkMatterProfileScaleRadiusTmp'
                change                 = ET.SubElement(parameters,'change')
                change.attrib['type' ] = "update"
                change.attrib['path' ] = 'darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/useMeanConcentration'
                change.attrib['value'] = 'false'
            else:
                # Otherwise, switch to the power-law model to set self-consistent initial conditions. Note that we set zero
                # scatter in this model - scatter will be added by the `johnson2021` class (therefore we set the scatter expected
                # for the power law model in the johnson2021 class instead).
                change                                       = ET.SubElement(parameters                  ,'change'                         )
                darkMatterProfileScaleRadius                 = ET.SubElement(change                      ,'darkMatterProfileScaleRadius'   )
                radiusLow_                                   = ET.SubElement(darkMatterProfileScaleRadius,'radiusLow'                      )
                radiusHigh_                                  = ET.SubElement(darkMatterProfileScaleRadius,'radiusHigh'                     )
                radiusTransition_                            = ET.SubElement(darkMatterProfileScaleRadius,'radiusTransition'               )
                radiusWidth_                                 = ET.SubElement(darkMatterProfileScaleRadius,'radiusWidth'                    )
                massLow_                                     = ET.SubElement(darkMatterProfileScaleRadius,'massLow'                        )
                massHigh_                                    = ET.SubElement(darkMatterProfileScaleRadius,'massHigh'                       )
                massTransition_                              = ET.SubElement(darkMatterProfileScaleRadius,'massTransition'                 )
                massWidth_                                   = ET.SubElement(darkMatterProfileScaleRadius,'massWidth'                      )
                expansionFactorLow_                          = ET.SubElement(darkMatterProfileScaleRadius,'expansionFactorLow'             )
                expansionFactorHigh_                         = ET.SubElement(darkMatterProfileScaleRadius,'expansionFactorHigh'            )
                expansionFactorTransition_                   = ET.SubElement(darkMatterProfileScaleRadius,'expansionFactorTransition'      )
                expansionFactorWidth_                        = ET.SubElement(darkMatterProfileScaleRadius,'expansionFactorWidth'           )
                scatter_                                     = ET.SubElement(darkMatterProfileScaleRadius,'scatter'                        )
                correlationRateDecay_                        = ET.SubElement(darkMatterProfileScaleRadius,'correlationRateDecay'           )
                correlationExponent_                         = ET.SubElement(darkMatterProfileScaleRadius,'correlationExponent'            )
                change                      .attrib['path' ] = 'darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/darkMatterProfileScaleRadius'
                change                      .attrib['type' ] = 'replace'
                darkMatterProfileScaleRadius.attrib['value'] = "powerLaw"
                radiusLow_                  .attrib['value'] = str(radiusLow                )
                radiusHigh_                 .attrib['value'] = str(radiusHigh               )
                radiusTransition_           .attrib['value'] = str(radiusTransition         )
                radiusWidth_                .attrib['value'] = str(radiusWidth              )
                massLow_                    .attrib['value'] = str(massLow                  )
                massHigh_                   .attrib['value'] = str(massHigh                 )
                massTransition_             .attrib['value'] = str(massTransition           )
                massWidth_                  .attrib['value'] = str(massWidth                )
                expansionFactorLow_         .attrib['value'] = str(expansionFactorLow       )
                expansionFactorHigh_        .attrib['value'] = str(expansionFactorHigh      )
                expansionFactorTransition_  .attrib['value'] = str(expansionFactorTransition)
                expansionFactorWidth_       .attrib['value'] = str(expansionFactorWidth     )
                scatter_                    .attrib['value'] = "0.000"
                correlationRateDecay_       .attrib['value'] = str(correlationRateDecay     )
                correlationExponent_        .attrib['value'] = str(correlationExponent      )
                change                                       = ET.SubElement(parameters,'change')
                change                      .attrib['path' ] = 'darkMatterProfileScaleRadius/darkMatterProfileScaleRadius/scatter'
                change                      .attrib['type' ] = 'update'
                change                      .attrib['value'] = str(modelScatter)

            id = re.sub("/"         ,":",id)
            id = re.sub("\[@value='","-",id)
            id = re.sub("'\]"       ,"" ,id)
            # Append tree mass and redshift to the id.
            id = id+f'_massTree{massTrees:.1e}_redshift{redshift:.1f}'
            # Store model ID and attributes.
            model['isResolutionModel'] = isResolutionModel
            model['massTrees'        ] = massTrees
            model['redshift'         ] = redshift
            model['id'               ] = id
            # Set file names.
            path               = args.outputPath+'/'
            parameterFileName  = f'{path}/{id}.xml'
            galacticusFileName = f'{path}/{id}.hdf5'
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "append"
            change.attrib['path' ] = ""
            outputFileName         = ET.SubElement(change,"outputFileName")
            outputFileName.attrib['value'] = galacticusFileName
            # Output a parameter changes file.
            with open(parameterFileName, 'wb') as parameterFile:
                parameterTree = ET.ElementTree(parameters)
                parameterTree.write(parameterFile, pretty_print=True)
    
            # Estimate the memory required per tree. This is a simple model calibrated to a few test runs. Results are in megabytes, as
            # required by SLURM.
            memoryBase              = 9.00
            memoryNorm              = 4.00e-2
            memoryExponent          = 0.99
            memoryBuffer            = 1.50
            memoryBoostExponent     = 0.70
            memorySubsampleExponent = 0.70
            memoryBoost             = 1.00
            if "mergerTreeBuilder/branchIntervalStep" in model and model['mergerTreeBuilder/branchIntervalStep'] == "true":
                memoryBoost *= 5.00
            if "mergerTreeBuildController/mergerTreeBuildController/massThreshold" in model:
                memoryBoost *= (massTrees*model['massResolutionFractional']/model['mergerTreeBuildController/massThreshold'])**(memoryExponent*memorySubsampleExponent)
            if "mergerTreeBuilder/accretionLimit"                  in model:
                memoryBoost *= (0.1/model['mergerTreeBuilder/accretionLimit'                 ])**memoryBoostExponent
            if "mergerTreeBuilder/mergeProbability"                in model:
                memoryBoost *= (0.1/model['mergerTreeBuilder/mergeProbability'               ])**memoryBoostExponent
            if "mergerTreeBranchingProbability/accuracyFirstOrder" in model:
                memoryBoost *= (0.1/model['mergerTreeBranchingProbability/accuracyFirstOrder'])**memoryBoostExponent
            memoryPerTree  = int(memoryBuffer*(memoryBase+memoryBoost*memoryNorm*(1.0/model['massResolutionFractional'])**memoryExponent))
            # Construct a job for this validation model.
            job = {
                "label":                   "convergenceConcentration_"+id         ,
                "launchFile":         path+"convergenceConcentration_"+id+".slurm",
                "logOutput":          path+"convergenceConcentration_"+id+".out"  ,
                "logError":           path+"convergenceConcentration_"+id+".err"  ,
                "nodes":               1,
                "memoryPerThread":    memoryPerTree,
                "walltime":           "1-00:00:00",
                "command":            timeCommand+" -v ./Galacticus.exe testSuite/parameters/validate_concentrationConvergence.xml "+parameterFileName,
                "onCompletion":       extractResults,
                "optimizeFor":        "nodes",
                "galacticusFileName": galacticusFileName,
                "results":            results,        
                "id":                 id,
                "model":              model
            }
            if os.path.exists(galacticusFileName):
                extractResults(job)
            else:
                jobs.append(job)

                # Submit all jobs and wait for completion.
manager.submitJobs(jobs)

# Fit a scale radius model to the full output results.
## Perform the fit.
parametersInitial       = np.array([radiusLow,radiusHigh,radiusTransition,radiusWidth,massLow,massHigh,massTransition,massWidth,expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth])
fit                     = minimize(radiusScaleFit, parametersInitial, args=(results), method='nelder-mead',options={'xatol': 1e-4, 'disp': False, 'maxiter': 10000})
fileModel = open(args.outputPath+'/powerLawModelFitData.csv','w')
fileModel.write("# Fits to all halos used in constraining the 'simple power law' model of halo scale radii.\n")
fileModel.write("log10Mass, log10ExpansionFactor, log10PeakHeight, log10RadiusScale\n")
fitMetric, modelScatterFit = radiusScaleFit(fit['x'],results,computeResidual=True,dataFile=fileModel)
fileModel.close()
## Extract the best fit parameter values.
radiusLowFit, radiusHighFit, radiusTransitionFit, radiusWidthFit, massLowFit, massHighFit, massTransitionFit, massWidthFit, expansionFactorLowFit, expansionFactorHighFit, expansionFactorTransitionFit, expansionFactorWidthFit = fit['x']
print( "\n"                                                                                                                                                  )
print( "Simple power-law scale radius model results:"                                                                                                        )
print( "   r' = r(ν) (M/M₁₂)^α(ν) a^β(ν)"                                                                                                                    )
print( " where:"                                                                                                                                             )
print( "  M   = halo mass"                                                                                                                                   )
print( "  a   = expansion factor"                                                                                                                            )
print( "  ν   = peak height"                                                                                                                                 )
print( "  M₁₂ = 10¹²M☉\n"                                                                                                                                    )
print( "r(ν), α(ν), β(ν) are sigmoid functions of the form:\n"                                                                                               )
print( "  y(x) = y₀+(y₁-y₀)/(1+exp[-(x-yᵥ)/Δy]),\n"                                                                                                           )
print( "and"                                                                                                                                                  )
print(f"r₀, r₁, rᵥ, Δr = {radiusLowFit         :+.4f}, {radiusHighFit         :+.4f}, {radiusTransitionFit         :+.4f}, {radiusWidthFit         :+.4f} Mpc")
print(f"α₀, α₁, αᵥ, Δα = {massLowFit           :+.4f}, {massHighFit           :+.4f}, {massTransitionFit           :+.4f}, {massWidthFit           :+.4f}"    )
print(f"β₀, β₁, βᵥ, Δβ = {expansionFactorLowFit:+.4f}, {expansionFactorHighFit:+.4f}, {expansionFactorTransitionFit:+.4f}, {expansionFactorWidthFit:+.4f}"    )
print(f"  Δlog₁₀(r/r') = {modelScatterFit      :+.4f}"                                                                                                        )

# Determine the correlation structure in the scale radius.
## Construct a list of offsets in log₁₀rₛ (relative to the best-fit simple power law model) on a fixed grid of Δ|log₁₀M|
offsetsLogRadius = []
for id, result in results.items():
    # Skip non-full output models.
    if not result['fullOutput']:
        continue
    # Define a grid of logarithmic, relative masses on which we will compute the correlation matrix.
    deltaLogMassGrid = np.linspace(np.log10(result['model']['massResolutionFractional']*result['factorMassResolution']*args.factorMassTrust),0.0,10)
    # Compute mean scale radii based on the power-law model.
    radiusScaleMean = radiusMeanPowerLaw(radiusLow,radiusHigh,radiusTransition,radiusWidth,massLow,massHigh,massTransition,massWidth,expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth,result['nodeMassesMainBranch'],result['nodeExpansionFactorsMainBranch'],result['nodePeakHeightsMainBranch'])
    # Iterate over trees.
    countTrees       = int(max(result['nodeTreeIndexMainBranch'       ]))
    expansionFactor0 =     max(result['nodeExpansionFactorsMainBranch'])
    for indexTree in range(1,countTrees+1):
        select              = (result['nodeTreeIndexMainBranch'] == indexTree)
        select0             = (result['nodeTreeIndexMainBranch'] == indexTree) & (result['nodeExpansionFactorsMainBranch'] == expansionFactor0)
        deltaLogMass        = np.log10(result['nodeMassesMainBranch'][select]/result         ['nodeMassesMainBranch'][select0])
        offsetLogRadius     = np.log10(result['nodeRadiiMainBranch' ][select]/radiusScaleMean                        [select ])
        offsetLogRadiusGrid = np.interp(deltaLogMassGrid,deltaLogMass,offsetLogRadius)
        offsetsLogRadius.append(offsetLogRadiusGrid)
# Stack the offsets and compute their correlation matrix.
offsetsLogRadius  = np.vstack  (offsetsLogRadius             )
correlationMatrix = np.corrcoef(offsetsLogRadius,rowvar=False)
# Fit a simple model for the structure of this correlation matrix.
parametersInitial   = np.array([correlationRateDecay,correlationExponent])
fit                 = minimize(correlationFit, parametersInitial, args=(correlationMatrix,deltaLogMassGrid), method='nelder-mead',options={'xatol': 1e-4, 'disp': False, 'maxiter': 10000})
correlationRateDecay, correlationExponent = fit['x']
# Report on the model.
print( "\n"                                                                                    )
print( "log₁₀rₛ correlation model results:"                                                     )
print( "   C_{Δlog₁₀rₛ}(Δ|log₁₀M|) = exp(-α[Δ|log₁₀M|]^β)"                                      )
print( " where:"                                                                               )
print( "  Δlog₁₀rₛ   = offset in log₁₀ of scale radius from the best-fit simple power law model")
print( "  Δ|log₁₀M| = absolute offset in log₁₀ of halo mass"                                   )
print(f"  α         = {correlationRateDecay:+.4f}"                                             )
print(f"  β         = {correlationExponent :+.4f}"                                             )
# Plot the correlation and fit.
deltaLogMass = np.zeros(correlationMatrix.shape)
for i in range(len(deltaLogMassGrid)):
   for j in range(len(deltaLogMassGrid)):
       deltaLogMass[i,j] = np.abs(deltaLogMassGrid[i]-deltaLogMassGrid[j])
figure, axes = plt.subplots(tight_layout=True)
axes.set_title(f'Scale radius correlation')
axes.set_ylim([0.0,1.1])
axes.set_xlabel('$\Delta |\log_{10} M|$')
axes.set_ylabel('$C(\Delta |\log_{10} M|)$' )
axes.plot(deltaLogMass.flatten(),correlationMatrix.flatten(),marker='.',linestyle='',color='red',label="Measured")
deltaLogMassFit      = np.linspace(0.0,1.2,100)
correlationMatrixFit = np.exp(-correlationRateDecay*np.abs(deltaLogMassFit)**correlationExponent)
axes.plot(deltaLogMassFit,correlationMatrixFit,marker='',linestyle='-',color='blue',label="Fit")
plt.savefig(f'{args.outputPath}/correlation.pdf')  
plt.clf()
# Store the data and fit to the output structure.
output['correlationModel'] = {
    "deltaLogMass":    list(deltaLogMass        .flatten()),
    "correlation":     list(correlationMatrix   .flatten()),
    "deltaLogMassFit": list(deltaLogMassFit               ),
    "correlationFit":  list(correlationMatrixFit          )
    }
        
# Plot the results for the primary resolution convergence.
## Construct radii for plotting.
radii     = 10.0**radiiScaleCentral
## Extract colors for models from a color map.
colorMap  = plt.colormaps['viridis']
colors    = colorMap(np.linspace(0, 1, countResolutions))
## Initialize a fit metric.
fitMetric = 0.0
## Iterate over tree masses.
for massTrees in massesTree:
    # Iterate over redshifts.
    for redshift in redshifts:
        # Find the ID of the corresponding Diemer & Joyce (2019) model.
        idDJ2019     = f'massResolutionFractional1.0e-02_fullOutputFalse_darkMatterProfileScaleRadiusConcentration_massTree{massTrees:.1e}_redshift{redshift:.1f}'
        # Find the ID of the corresponding full output model.
        idFullOutput = f'massResolutionFractional{massResolutionFractionalMinimum:.1e}_fullOutputTrue_massTree{massTrees:.1e}_redshift{redshift:.1f}'
        # Begin plot creation.
        figure, axes = plt.subplots(tight_layout=True)
        axes.set_title(f'Scale radius convergence')
        axes.set_ylim([0.0,5.0])
        axes.set_xscale('log')
        axes.set_xlabel('$r_\mathrm{s}$ [kpc]')
        axes.set_ylabel('$\mathrm{d}p/\mathrm{d}\log r_\mathrm{s}$' )
        ## Iterate over results.
        nonZeroBest        = None
        iColor             = -1
        resultsResolutions = []
        for id, result in results.items():
            # Skip non-resolution models.
            if not result['isResolutionModel']:
                continue
            # Skip full output models.
            if result['fullOutput']:
                continue
            # Skip non-matching tree masses.
            if result['massTrees'] != massTrees:
                continue
            # Skip non-matching redshifts.
            if result['redshift' ] != redshift:
                continue
            # Store result.
            resultsResolutions.append(
                {
                    "massResolutionFractional":      result    ['massResolutionFractional'] ,
                    "radius":                   list(radiiScale                            ),
                    "pdf":                      list(result    ['radiusScaleFunction'     ]),
                    "error":                    list(result    ['radiusScaleFunctionError'])
                }
            )
            # Plot result.
            resolution = '{:.1e}'.format(num2tex(result['massResolutionFractional']))
            iColor     = iColor+1
            nonZero    = result['radiusScaleFunction'] > 0.0
            if result['model']['massResolutionFractional'] == massResolutionFractionalMinimum:
                nonZeroBest = nonZero
            axes.plot    (radii[nonZero],result['radiusScaleFunction'][nonZero],marker='o',markerfacecolor=colors[iColor],markeredgecolor=colors[iColor],linestyle='',label=f"$f_\mathrm{{res}}={resolution}$; $(\mu,\sigma) = ({result['radiusScaleLogarithmicMean']:.2f},{result['radiusScaleLogarithmicScatter']:.2f})$")
            axes.errorbar(radii[nonZero],result['radiusScaleFunction'][nonZero],yerr=result['radiusScaleFunctionError'][nonZero],ecolor=colors[iColor],fmt='none')
            if result['model']['massResolutionFractional'] == massResolutionFractionalMinimum:
                axes.axvline (10.0**result['radiusScaleLogarithmicMean'])
            # Accumulate the fit metric.
            fitMetric += np.sum((result['radiusScaleFunction']-results[idDJ2019]['radiusScaleFunction'])**2)
        # Add the reference Diemer & Joyce (2019) model.
        nonZero = results[idDJ2019]['radiusScaleFunction'] > 0.0
        axes.plot    (radii[nonZero],results[idDJ2019]['radiusScaleFunction'][nonZero],marker='',linestyle='-',color='black',label=f"D&J2019; $(\mu,\sigma) = ({results[idDJ2019]['radiusScaleLogarithmicMean']:.2f},{results[idDJ2019]['radiusScaleLogarithmicScatter']:.2f})$")
        axes.errorbar(radii[nonZero],results[idDJ2019]['radiusScaleFunction'][nonZero],yerr=results[idDJ2019]['radiusScaleFunctionError'][nonZero],ecolor='black',fmt='none')
        # Add simple power-law model.
        expansionFactor = 1.0/(1.0+redshift)
        radiusMeanModel = 1.0e3*radiusMeanPowerLaw(radiusLow,radiusHigh,radiusTransition,radiusWidth,massLow,massHigh,massTransition,massWidth,expansionFactorLow,expansionFactorHigh,expansionFactorTransition,expansionFactorWidth,massTrees,expansionFactor,results[idFullOutput]['nodePeakHeightFinal'])
        resultsModel    = np.exp(-0.5*(np.log(radii[nonZeroBest]/radiusMeanModel)/modelScatter/np.log(10.0))**2)/np.sqrt(2.0*np.pi)/modelScatter
        axes.plot    (radii[nonZeroBest],resultsModel,marker='',linestyle='-',color='gray',label=f"power-law; $(\mu,\sigma) = ({np.log10(radiusMeanModel):.2f},{modelScatter:.2f})$")
        radiusMeanModel = 1.0e3*radiusMeanPowerLaw(radiusLowFit,radiusHighFit,radiusTransitionFit,radiusWidthFit,massLowFit,massHighFit,massTransitionFit,massWidthFit,expansionFactorLowFit,expansionFactorHighFit,expansionFactorTransitionFit,expansionFactorWidthFit,massTrees,expansionFactor,results[idFullOutput]['nodePeakHeightFinal'])
        resultsModel    = np.exp(-0.5*(np.log(radii[nonZeroBest]/radiusMeanModel)/modelScatterFit/np.log(10.0))**2)/np.sqrt(2.0*np.pi)/modelScatterFit
        axes.plot    (radii[nonZeroBest],resultsModel,marker='',linestyle='dotted',color='gray')
        axes.axvline (radiusMeanModel,linestyle="dotted")
        # Store our results.
        resultsMassResolutionFractional = {
            "massTrees": massTrees         ,
            "redshift":  redshift          ,
            "results":   resultsResolutions
            }
        if not 'massResolutionFractional' in output['results']:
            output['results']['massResolutionFractional'] = []
        output['results']['massResolutionFractional'].append(resultsMassResolutionFractional)
        # Finalize the plot.
        axes.legend(loc='upper center',bbox_to_anchor=(0.40,0.95),fancybox=True,shadow=True,ncol=2)
        plt.savefig(f'{args.outputPath}/convergence_massTree{massTrees:.1e}_redshift{redshift:.1f}.pdf')  
        plt.clf()

# Plot the results for other convergence tests.
## Iterate over tree masses.
for massTrees in massesTree:
    # Iterate over redshifts.
    for redshift in redshifts:
        ## Iterate over results.
        for id, result in results.items():
            # Skip resolution models.
            if result['isResolutionModel']:
                continue
            # Skip full output models.
            if result['fullOutput']:
                continue
            # Skip concentration models.
            if 'darkMatterProfileScaleRadius' in result['model']:
                continue
            # Skip non-matching tree masses.
            if result['massTrees'] != massTrees:
                continue
            # Skip non-matching redshifts.
            if result['redshift' ] != redshift:
                continue
            # Find the reference model.
            referenceModelFound = False
            for idReference, resultReference in results.items():
                # Skip non-resolution models.
                if not resultReference['isResolutionModel']:
                    continue
                # Skip full output models.
                if resultReference['fullOutput']:
                    continue
                # Skip concentration models.
                if 'darkMatterProfileScaleRadius' in result['model']:
                    continue
                # Skip non-matching tree masses.
                if resultReference['massTrees'] != massTrees:
                    continue
                # Skip non-matching redshifts.
                if resultReference['redshift' ] != redshift:
                    continue
                # Reference model is found.
                referenceModelFound = True
                break
            if not referenceModelFound:
                print("FAIL: reference model not found")
                sys.exit()

            # Find the attribute of this model.
            attribute       = [x for x in result['model'].keys() if x not in ["id", "isResolutionModel", "massTrees", "redshift", "massResolutionFractional", "darkMatterProfileScaleRadius", "fullOutput"]][0]
            attributeSuffix = re.sub(r'.*/(.+)',r'\1',attribute)
            ## Begin plot creation.
            figure, axes = plt.subplots(tight_layout=True)
            axes.set_title(f'Scale radius convergence')
            axes.set_ylim([0.0,4.0])
            axes.set_xscale('log')
            axes.set_xlabel('$r_\mathrm{s}$ [kpc]')
            axes.set_ylabel('$\mathrm{d}p/\mathrm{d}\log r_\mathrm{s}$' )
            # Store results.
            output['results'][attributeSuffix] = {
                "massTrees":                massTrees                                        ,
                "redshift":                 redshift                                         ,
                "massResolutionFractional":      result         ['massResolutionFractional'] ,
                "radius":                   list(radiiScale                                 ),
                "pdf":                      list(result         ['radiusScaleFunction'     ]),
                "error":                    list(result         ['radiusScaleFunctionError']),
                "pdfDefault":               list(resultReference['radiusScaleFunction'     ]),
                "errorDefault":             list(resultReference['radiusScaleFunctionError']),
                "default":                  defaults         [attribute],
                "updated":                  result  ['model'][attribute]
            }
            # Plot result.
            nonZero = (result['radiusScaleFunction'] > 0.0) | (resultReference['radiusScaleFunction'] > 0.0)
            axes.plot    (radii[nonZero],result         ['radiusScaleFunction'][nonZero],marker='o',markerfacecolor=colors[0],markeredgecolor=colors[0],linestyle='',label='test'     )
            axes.errorbar(radii[nonZero],result         ['radiusScaleFunction'][nonZero],yerr=result         ['radiusScaleFunctionError'][nonZero],ecolor=colors[0],fmt='none')
            axes.plot    (radii[nonZero],resultReference['radiusScaleFunction'][nonZero],marker='o',markerfacecolor=colors[1],markeredgecolor=colors[1],linestyle='',label='reference')
            axes.errorbar(radii[nonZero],resultReference['radiusScaleFunction'][nonZero],yerr=resultReference['radiusScaleFunctionError'][nonZero],ecolor=colors[1],fmt='none')
            # Finalize the plot.
            axes.legend(loc='upper center',bbox_to_anchor=(0.40,0.95),fancybox=True,shadow=True,ncol=2)
            plt.savefig(f'{args.outputPath}/convergence_{attributeSuffix}.pdf')  
            plt.clf()

# Output a file containing the fit metric.
fitMetricOutput = open(f'{args.outputPath}/fitMetric.txt', "w")
fitMetricOutput.write(str(fitMetric)+"\n")
fitMetricOutput.close()

# Add the command line to the output structure.
output['commandLine'] = " ".join(sys.argv)

# Output JSON file containing all results.
f = codecs.open(f'{args.outputPath}/results.json', "w", "utf-8")
f.write("window.VALIDATION_DATA = ")
f.write(json.dumps(output,indent=4,ensure_ascii=False))
f.close()
    
print("SUCCESS: convergence model")
