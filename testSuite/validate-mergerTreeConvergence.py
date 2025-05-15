#!/usr/bin/env python3
import sys
import os
import subprocess
import h5py
import numpy as np
import argparse
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

# Run models to measure convergence in merger tree construction with mass resolution.
# Andrew Benson (11-February-2025)

# This validation will:
#  * generate a set of merger trees spanning a range of resolutions and numerical tolerances;
#  * fit a model for the completeness of the progenitor mass function as a function of progenitor mass;
#  * create plots showing this completeness model;
#  * create plots showing convergence with other numerical tolerances;
#  * output JSON files containing the results for inclusion into a web page.
#
# The validation is typically run as:
#
#  ./convergence.py --outputPath /scratch/mergerTreeConvergence --jobMaximum 10 --waitOnSubmit 1 --waitOnActive 30
#
# where the first option controls the directory to which results will be output and the other (optional) arguments are passed to
# the queue manager.
#
# Once complete (note that these high-resolution models can take a very long time to finish - ~2 weeks) plots are created and
# output to the `outputPath` showing the convergence with respect to each numerical parameter considered. Additionally, a file
# `results.json` is output to `outputPath`. This contains the raw data on all convergence tests and should be copied to
# `dev/valid/mergerTreeConvergence/` in the `gh-pages` branch and committed there. It will then appear as a page of convergence
# plots at https://galacticusorg.github.io/galacticus/dev/valid/mergerTreeConvergence.

# Convert command line arguments to ints.
def restricted_int(x):
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not an integer literal" % (x,))
    return x

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='convergence.py',description='Compute convergence in the progenitor mass function.')
parser.add_argument('--outputPath'  ,action='store',default="."        ,help='the path to which data should be written'               )
parser.add_argument('--partition'   ,action='store'                    ,help='the partition to which to submit jobs'                  )
parser.add_argument('--jobMaximum'  ,action='store',type=restricted_int,help='the maximum number of active jobs to allow'             )
parser.add_argument('--waitOnSubmit',action='store',type=restricted_int,help='the time (in seconds) to wait after submitting each job')
parser.add_argument('--waitOnActive',action='store',type=restricted_int,help='the time (in seconds) to wait after polling active jobs')
args = parser.parse_args()

def extractResults(job):
    # Wait for the model to appear.
    while not os.path.exists(job['galacticusFileName']):
        time.sleep(1)
    # Extract data from the Galacticus model.
    galacticus                = h5py.File(job['galacticusFileName'],"r")
    outputs                   = galacticus['Outputs']
    job['results'][job['id']] = {}
    for outputName, output in outputs.items():
        matched = re.match("Output(\d+)",outputName)
        if not matched:
            continue
        indexOutput         = matched.group(1)
        nodes               = output['nodeData']
        expansionFactor     = output.attrs['outputExpansionFactor']
        if expansionFactor >= 1.0:
            continue
        criticalOverdensity = output.attrs['criticalOverdensity']
        barrierHeight       = criticalOverdensity/expansionFactor
        redshift            = 1.0/expansionFactor-1.0
        massHalo            = nodes['basicMass'            ][:]
        weight              = nodes['nodeSubsamplingWeight'][:]
        massResolution      = galacticus['Parameters/mergerTreeMassResolution'].attrs['massResolution']
        # Construct halo mass CDF.
        massFunctionCount, _ = np.histogram(np.log10(massHalo),bins=masses,weights=weight)
        massFunction         = massFunctionCount/massesDelta/float(countTrees)    
        # Store the results.
        job['results'][job['id']][indexOutput] = {
            "massResolution":    massResolution                         ,
            "massFunctionCount": massFunctionCount                      , 
            "redshift":          redshift                               ,
            "barrierHeight":     barrierHeight                          ,
            "isResolutionModel": job['model']['isResolutionModel']
        }

def completeness(massesScaleFree,barrierHeight,a,b):
    # Evaluate the completeness function.
    resolved = massesScaleFree > 1.0
    factor   = np.zeros([len(massesScaleFree)])
    factor[resolved] = np.exp(-a*barrierHeight/(massesScaleFree[resolved]-1.0)**b)
    return factor
    
def massFunctionFit(parameters,massesLogarithmic,results,resultBestID):
    # Fitting function for the progenitor mass function.
    ## Extract the parameters of the completeness model.
    a, b   = parameters
    ## Construct masses.
    masses = 10.0**massesLogarithmic
    # Iterate over results, accumulating the fit metric.
    fitMetric = 0.0
    for id, result in results.items():
        # Skip the "best" model as we use it as our reference, and non-resolution models.
        if id == resultBestID or not result['1']['isResolutionModel']:
            continue
        # Iterate over outputs.
        for indexOutput in result.keys():
            # Get masses scaled to the mass resolution.
            massesScaleFree     = masses/result               [indexOutput]['massResolution']
            massesScaleFreeBest = masses/results[resultBestID][indexOutput]['massResolution']
            # Get the barrier height.
            barrierHeight       =        result               [indexOutput]['barrierHeight' ]
            # Compute completeness factors.
            factorCompleteness     = completeness(massesScaleFree    ,barrierHeight,a,b)
            factorCompletenessBest = completeness(massesScaleFreeBest,barrierHeight,a,b)
            # Determine which masses are resolved.
            resolved = (massesScaleFree >= 1.0)
            # Find the fractional difference between the count of halos in this result, and that in the best model after applying our
            # completeness model to it.
            difference = (
                          +result                              [indexOutput]['massFunctionCount'][resolved]
                          -factorCompleteness                                                    [resolved]
                          /factorCompletenessBest                                                [resolved]
                          *results               [resultBestID][indexOutput]['massFunctionCount'][resolved]
                         )
            variance   = (
                          +1.0
                          /result                              [indexOutput]['massFunctionCount'][resolved]
                          +factorCompletenessBest                                                [resolved]
                          /factorCompleteness                                                    [resolved]
                          /results               [resultBestID][indexOutput]['massFunctionCount'][resolved]
                         )
            # Accumulate our fit metric.
            fitMetric += np.sum(difference**2/variance)
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
    "parameterFile": "testSuite/parameters/validate_mergerTreeConvergence.xml",
    "results": {}
}
    
# Set default model attributes.
defaults = {
    "mergerTreeBranchingProbability/accuracyFirstOrder":                                     0.1,
    "mergerTreeBuilder/accretionLimit":                                                      0.1,
    "mergerTreeBuilder/mergeProbability":                                                    0.1,
    "mergerTreeBuildController/mergerTreeBuildController[@value='subsample']/massThreshold": 0.0,
    "mergerTreeBuilder/branchIntervalStep":                                                  "false"
}

# Specify the number of trees to run.
countTrees    = 240

# Specify the tree mass.
massTrees     = 1.0e12

# Define masses for the mass function.
masses        = np.linspace(4.0,9.0,41)
massesDelta   = masses[1]-masses[0]
massesCentral = 0.5*(masses[0:-1]+masses[1:])

# Define models to run.
models = [
    {
        "massResolution":                                                                        1.0e4
    },
    {
        "massResolution":                                                                        1.0e5
    },
    {
        "massResolution":                                                                        1.0e6
    },
    {
        "massResolution":                                                                        1.0e7
    },
    {
        "massResolution":                                                                        1.0e4,
        "mergerTreeBuilder/accretionLimit":                                                      0.01
    },
    {
        "massResolution":                                                                        1.0e4,
        "mergerTreeBuilder/mergeProbability":                                                    0.01
    },
    {
        "massResolution":                                                                        1.0e4,
        "mergerTreeBranchingProbability/accuracyFirstOrder":                                     0.01
    },
    {
        "massResolution":                                                                        1.0e4,
        "mergerTreeBuilder/branchIntervalStep":                                                  "true"
    },
    {
        "massResolution":                                                                        1.0e4,
        "mergerTreeBuildController/mergerTreeBuildController[@value='subsample']/massThreshold": 1.0e6
    }
]

# Iterate over mass resolutions.
massResolutionMinimum = min(    map(lambda x: x['massResolution'],models) )
countResolutions      = len(set(map(lambda x: x['massResolution'],models)))
jobs                  = []
results               = {}
resultBestID          = None
for model in models:
    # Construct the file of parameter changes for this model.
    parameters = ET.Element("changes")
    # Set the task.
    change                 = ET.SubElement(parameters,'change')
    change.attrib['type' ] = "append"
    change.attrib['path' ] = ""
    taskMulti                           = ET.SubElement(change   ,"task"                               )
    taskHMF                             = ET.SubElement(taskMulti,"task"                               )
    taskEvolve                          = ET.SubElement(taskMulti,"task"                               )
    includeUnevolvedSubhaloMassFunction = ET.SubElement(taskHMF  ,"includeUnevolvedSubhaloMassFunction")
    includeMassAccretionRate            = ET.SubElement(taskHMF  ,"includeMassAccretionRate"           )
    criticalOverdensity                 = ET.SubElement(taskHMF  ,"criticalOverdensity"                )
    darkMatterHaloBias                  = ET.SubElement(taskHMF  ,"darkMatterHaloBias"                 )
    taskMulti                          .attrib['value'] = "multi"
    taskHMF                            .attrib['value'] = "haloMassFunction"
    taskEvolve                         .attrib['value'] = "evolveForests"
    includeUnevolvedSubhaloMassFunction.attrib['value'] = "false"
    includeMassAccretionRate           .attrib['value'] = "false"
    criticalOverdensity                .attrib['value'] = "sphericalCollapseClsnlssMttrCsmlgclCnstnt"
    darkMatterHaloBias                 .attrib['value'] = "pressSchechter"
    # Nullify unused components.
    change                 = ET.SubElement(parameters,'change')
    change.attrib['type' ] = "update"
    change.attrib['path' ] = "componentSatellite"
    change.attrib['value'] = "null"
    change                 = ET.SubElement(parameters,'change')
    change.attrib['type' ] = "update"
    change.attrib['path' ] = "componentSpin"
    change.attrib['value'] = "null"
    # Set simple dark matter profile.
    change                       = ET.SubElement(parameters,'change')
    change.attrib['type' ]       = "append"
    change.attrib['path' ]       = ""
    darkMatterProfileDMO         = ET.SubElement(change,"darkMatterProfileDMO"        )
    darkMatterProfileScaleRadius = ET.SubElement(change,"darkMatterProfileScaleRadius")
    darkMatterProfileDMO        .attrib['value'] = "isothermal"
    darkMatterProfileScaleRadius.attrib['value'] = "zero"
    # Set a non-evolving evolver.
    change                 = ET.SubElement(parameters,'change')
    change.attrib['type' ] = "replace"
    change.attrib['path' ] = "mergerTreeEvolver"
    mergerTreeEvolver      = ET.SubElement(change           ,"mergerTreeEvolver")
    pruneTree              = ET.SubElement(mergerTreeEvolver,"pruneTree"        )
    mergerTreeEvolver.attrib['value'] = "nonEvolving"
    pruneTree        .attrib['value'] = "true"
    change                 = ET.SubElement(parameters,'change')
    change.attrib['type' ] = "remove"
    change.attrib['path' ] = "nodeOperator"
    # Set merger tree masses.
    change                 = ET.SubElement(parameters,'change')
    change.attrib['type' ] = "append"
    change.attrib['path' ] = ""
    mergerTreeBuildMasses  = ET.SubElement(change               ,"mergerTreeBuildMasses")
    massTree               = ET.SubElement(mergerTreeBuildMasses,"massTree"             )
    treeCount              = ET.SubElement(mergerTreeBuildMasses,"treeCount"            )
    mergerTreeBuildMasses.attrib['value'] = "fixedMass"
    massTree             .attrib['value'] = str(massTrees )
    treeCount            .attrib['value'] = str(countTrees)
    # Set merger tree resolution.
    change                   = ET.SubElement(parameters,'change')
    change.attrib['type' ]   = "append"
    change.attrib['path' ]   = ""
    mergerTreeMassResolution = ET.SubElement(change                  ,"mergerTreeMassResolution")
    massResolution           = ET.SubElement(mergerTreeMassResolution,"massResolution"          )
    mergerTreeMassResolution.attrib['value'] = "fixed"
    massResolution          .attrib['value'] = str(model['massResolution'])
    # Set merger tree maximum redshift.
    change                 = ET.SubElement(parameters,'change')
    change.attrib['type' ] = "update"
    change.attrib['path' ] = "mergerTreeBuilder/redshiftMaximum"
    change.attrib['value'] = "1.0"
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
        if isinstance(model[attribute],float):
            id += f'{separator}{attribute}{model[attribute]:.1e}'
        else:
            id += f'{separator}{attribute}{model[attribute].capitalize()}'
        separator = "_"
        if attribute != "massResolution":
            isResolutionModel = False
            change                 = ET.SubElement(parameters,'change')
            change.attrib['type' ] = "update"
            change.attrib['path' ] = attribute
            change.attrib['value'] = str(model[attribute])
    id = re.sub("/"         ,":",id)
    id = re.sub("\[@value='","-",id)
    id = re.sub("'\]"       ,"" ,id)
    model['isResolutionModel'] = isResolutionModel
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
    # Indicate the best model.
    if model['isResolutionModel'] and model['massResolution'] == massResolutionMinimum:
        resultBestID = id
    # Output a parameter changes file.
    with open(parameterFileName, 'wb') as parameterFile:
        parameterTree = ET.ElementTree(parameters)
        parameterTree.write(parameterFile, pretty_print=True)

    # Estimate the memory required per tree. This is a simple model calibrated to a few test runs. Results are in megabytes, as
    # required by SLURM.
    memoryBase              = 65.00
    memoryNorm              =  1.10e-4
    memoryExponent          =  0.99
    memoryBuffer            =  1.50
    memoryBoostExponent     =  0.70
    memorySubsampleExponent =  0.70
    memoryBoost             =  1.00
    if "mergerTreeBuilder/branchIntervalStep" in model and model['mergerTreeBuilder/branchIntervalStep'] == "true":
        memoryBoost *= 5.00
    if "mergerTreeBuildController/mergerTreeBuildController/massThreshold" in model:
        memoryBoost *= (model['massResolution']/model['mergerTreeBuildController/mergerTreeBuildController/massThreshold'])**(memoryExponent*memorySubsampleExponent)
    if "mergerTreeBuilder/accretionLimit"                  in model:
        memoryBoost *= (0.1/model['mergerTreeBuilder/accretionLimit'                 ])**memoryBoostExponent
    if "mergerTreeBuilder/mergeProbability"                in model:
        memoryBoost *= (0.1/model['mergerTreeBuilder/mergeProbability'               ])**memoryBoostExponent
    if "mergerTreeBranchingProbability/accuracyFirstOrder" in model:
        memoryBoost *= (0.1/model['mergerTreeBranchingProbability/accuracyFirstOrder'])**memoryBoostExponent
    memoryPerTree  = int(memoryBuffer*(memoryBase+memoryBoost*memoryNorm*(massTrees/model['massResolution'])**memoryExponent))
    
    # Construct a job for this validation model.
    job = {
        "label":                   "convergencePMF_"+id         ,
        "launchFile":         path+"convergencePMF_"+id+".slurm",
        "logOutput":          path+"convergencePMF_"+id+".out"  ,
        "logError":           path+"convergencePMF_"+id+".err"  ,
        "nodes":               1,
        "memoryPerThread":    memoryPerTree,
        "command":            "/usr/bin/time -v ./Galacticus.exe testSuite/parameters/validate_mergerTreeConvergence.xml "+parameterFileName,
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

# Fit a completeness model to the accumulated results.
## Initial guesses for the parameters of the model.
a = 0.01
b = 1.16
## Perform the fit.
parametersInitial = np.array([a,b])
fit               = minimize(massFunctionFit, parametersInitial, args=(massesCentral, results, resultBestID), method='nelder-mead',options={'xatol': 1e-4, 'disp': False, 'maxiter': 10000})
## Extract the best fit parameter values.
a, b              = fit['x']

# Plot the results for the primary resolution convergence.
## Construct masses for plotting.
massesCentralFine = np.linspace(5.0,9.0,201)
masses            = 10.0**massesCentral
massesFine        = 10.0**massesCentralFine
## Extract colors for models from a color map.
colorMap = plt.colormaps['viridis']
colors   = colorMap(np.linspace(0, 1, countResolutions-1))
## Iterate over redshifts.
for indexOutput in results[resultBestID]:
    ## Begin plot creation.
    figure, axes = plt.subplots(tight_layout=True)
    axes.set_title(f'Progenitor mass function convergence at z={results[resultBestID][indexOutput]["redshift"]:.2f}')
    axes.set_ylim([-0.2,+0.1])
    axes.set_xscale('log')
    axes.annotate(f"$\\frac{{N(m_\mathrm{{prog}})}}{{N_0(m_\mathrm{{prog}})}}=\exp\left(-\\frac{{{a:.3f}}}{{[m_\mathrm{{prog}}/m_\mathrm{{res}}-1]^{{{b:.3f}}}}}\\right)$", xy=(1.7e7, -0.15))
    axes.set_xlabel('$m_\mathrm{prog}$')
    axes.set_ylabel('$N(m)/N_0(m)$' )
    ## Iterate over results.
    iColor             = -1
    resultsResolutions = []
    for id, result in results.items():
        # Skip the best model as we use it as our reference or it not a resolution model.
        if result[indexOutput]['massResolution'] <= results[resultBestID][indexOutput]['massResolution'] or not result[indexOutput]['isResolutionModel']:
            continue
        # Construct masses for model results.
        massesScaleFreeFine     = massesFine/result               [indexOutput]['massResolution']
        massesScaleFreeFineBest = massesFine/results[resultBestID][indexOutput]['massResolution']
        # Get the barrier height.
        barrierHeight           =            result               [indexOutput]['barrierHeight' ]
        # Evaluate the completeness function.
        factorCompleteness     = completeness(massesScaleFreeFine    ,barrierHeight,a,b)
        factorCompletenessBest = completeness(massesScaleFreeFineBest,barrierHeight,a,b)
        # Find offsets relative to the best model.
        resolved                         = masses     > result[indexOutput]['massResolution']
        resolvedFine                     = massesFine > result[indexOutput]['massResolution']
        result[indexOutput]['offset'     ] = result            [indexOutput]['massFunctionCount'][resolved    ]/results               [resultBestID][indexOutput]['massFunctionCount'][resolved    ]-1.0
        result[indexOutput]['offsetFit'  ] = factorCompleteness                                  [resolvedFine]/factorCompletenessBest                                                [resolvedFine]-1.0
        result[indexOutput]['offsetError'] = np.sqrt(
                                                     1.0/result[indexOutput]['massFunctionCount'][resolved]+1.0/results[resultBestID][indexOutput]['massFunctionCount'][resolved]
                                                    )*  (result[indexOutput]['offset'           ]+1.0)
        # Store result.
        if not 'massResolution' in output['results']:
            output['results']['massResolution'] = []
        resultsResolutions.append(
            {
                "massResolution": result[indexOutput]['massResolution'],
                "mass":           list(masses    [resolved    ]               ),
                "massFit":        list(massesFine[resolvedFine]               ),
                "ratio":          list(result    [indexOutput ]['offset'     ]),
                "error":          list(result    [indexOutput ]['offsetError']),
                "ratioFit":       list(result    [indexOutput ]['offsetFit'  ]),
            }
        )
        # Plot result.
        resolution = '{:.1e}'.format(num2tex(result[indexOutput]['massResolution']))
        iColor     = iColor+1
        axes.plot(masses    [resolved    ],result[indexOutput]['offset'   ],marker='o',markerfacecolor=colors[iColor],markeredgecolor=colors[iColor],linestyle='',label=f"$m_\mathrm{{res}}={resolution}\mathrm{{M}}_\odot$")
        axes.plot(massesFine[resolvedFine],result[indexOutput]['offsetFit'],marker='',linestyle='-',color=colors[iColor])        
    # Create our results section if necessary.
    if not 'massResolution' in output['results']:
        output['results']['massResolution'] = []
    # Append a section for this redshift.
    output['results']['massResolution'].append(
        {
            "redshift":    results[resultBestID][indexOutput]["redshift"],
            "resolutions": resultsResolutions
        }
    )
    # Finalize the plot.
    axes.legend(loc='upper center',bbox_to_anchor=(0.40,0.95),fancybox=True,shadow=True,ncol=2)
    plt.savefig(f'{args.outputPath}/convergence_z{results[resultBestID][indexOutput]["redshift"]:.2f}.pdf')  
    plt.clf()

# Add meta-data on best model and fit.
output['massResolutionFit'] = {
    "massResolution": results[resultBestID]['1']['massResolution'],
    "a":              a                                           ,
    "b":              b                                           ,
    "massTrees":      massTrees                                   ,
    "countTrees":     countTrees
}
    
# Plot the results for other convergence tests.
for model in models:
    if model['isResolutionModel']:
        continue
    modelReference = None
    for modelSeek in models:
        if not modelSeek['isResolutionModel']:
            continue
        if modelSeek['massResolution'] == model['massResolution']:
            modelReference = modelSeek
            break
    result          = results[model         ['id']]
    resultReference = results[modelReference['id']]
    # Find the attribute of this model.
    attribute       = [x for x in model.keys() if x not in ["id", "massResolution"]][0]
    attributeSuffix = re.sub(r'.*/(.+)',r'\1',attribute)
    # Iterate over outputs.
    for indexOutput in result:
        ## Begin plot creation.
        figure, axes = plt.subplots(tight_layout=True)
        axes.set_title(f'Progenitor mass function convergence at z={results[resultBestID][indexOutput]["redshift"]:.2f}')
        axes.set_ylim([-0.2,+0.1])
        axes.set_xscale('log')
        axes.set_xlabel('$m_\mathrm{prog}$')
        axes.set_ylabel('$N(m)/N_0(m)$' )
        axes.axhline(y=0.0, color='black', linestyle='dotted')
        axes.annotate(f"{attributeSuffix} = {defaults[attribute]} $\\rightarrow$ {model[attribute]}", xy=(1.7e7, -0.15))
        # Find offsets relative to the reference model.
        resolved = masses > result[indexOutput]['massResolution']
        result[indexOutput]['offset'     ] =             result[indexOutput]['massFunctionCount'][resolved]/    resultReference[indexOutput]['massFunctionCount'][resolved]-1.0
        result[indexOutput]['offsetError'] = np.sqrt(
                                                     1.0/result[indexOutput]['massFunctionCount'][resolved]+1.0/resultReference[indexOutput]['massFunctionCount'][resolved]
                                                    )*  (result[indexOutput]['offset'           ]+1.0)
        # Store result.
        if not attributeSuffix in output['results']:
            output['results'][attributeSuffix] = []
        output['results'][attributeSuffix].append(
            {
                "redshift": results[resultBestID][indexOutput]["redshift"],
                "mass":     list(masses[resolved]),
                "ratio":    list(result[indexOutput]['offset'     ]),
                "error":    list(result[indexOutput]['offsetError']),
                "default":  defaults[attribute],
                "updated":  model   [attribute]
            }
        )
        # Plot result.
        axes.plot(masses[resolved],result[indexOutput]['offset'],marker='o',markerfacecolor=colors[0],markeredgecolor=colors[0],linestyle='')
        # Finalize the plot.
        plt.savefig(f'{args.outputPath}/convergence_{attributeSuffix}_z{result[indexOutput]["redshift"]:.2f}.pdf')  
        plt.clf()

# Output JSON fle containing all results.
f = codecs.open(f'{args.outputPath}/results.json', "w", "utf-8")
f.write("window.VALIDATION_DATA = ")
f.write(json.dumps(output,indent=4,ensure_ascii=False))
f.close()
    
print("SUCCESS: convergence model")
