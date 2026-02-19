#!/usr/bin/env python3
import numpy as np
import h5py
import argparse
import re
import sys
import matplotlib.pyplot as plt
import lxml.etree as ET
from queue import LifoQueue
import os

# Make plots of all analyses that were output to a Galacticus model file(s).
# Andrew Benson (21-June-2024)

# Convert command line arguments to floats.
def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))
    return x

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='analysesPlot.py',description='Make plots of all analyses that were output to a Galacticus model file(s).')
parser.add_argument('filename',nargs='+')
parser.add_argument('--outputDirectory'          , action='store'     ,default='.'                      ,help='the directory to which plots should be output'                     )
parser.add_argument('--showFullRange'            , action='store_true'                                  ,help='show the full range of the y-axis, not just the target range'      )
parser.add_argument('--factorLogLikelihoodAccept', action='store'     ,default=0.5,type=restricted_float,help='the fractional shift in logℒ to accept when updating prior ranges')
parser.add_argument('--configFileName'           , action='store'                                       ,help='the name of the posterior sampling config file to modify'          )
args = parser.parse_args()

# Open a file for reports.
reports = open(args.outputDirectory+'/report.txt','w')

# Open all files, and find all available analyses.
models = LifoQueue()
for filename in args.filename:
    models.put({'filename': filename, 'group': 'analyses'})
# Find all available analyses.
analysisNames = {}
analyses      = []
while not models.empty():
    model = models.get()
    file  = h5py.File(model['filename'],"r")
    group = file[model['group']]
    # Look for subgroups of analyses.
    subgroupRE       = re.compile(r"^step\d+:chain\d+$")
    subgroupsPresent = False
    for analysisName in group.keys():
         if subgroupRE.match(analysisName):
            subgroupsPresent = True
            models.put({'filename': model['filename'], 'group': "analyses/"+analysisName})
    # Create a dictionary of all analyses that we find.
    if not subgroupsPresent:
        analyses.append({'filename': model['filename'], 'file': file, 'groupname': model['group'], 'group': group})
        for analysisName in group.keys():
            if isinstance(group[analysisName], h5py.Group):
                analysisNames[analysisName] = 1

# Record if we are plotting a single model.
singleModel = len(analyses) == 1

# Iterate over analyses finding the best model for each analysis, and overall.
analysisBest      = {"global": None}
logLikelihoodBest = {"global": -np.inf}
for analysisName in analysisNames.keys():
    analysisBest     .update({analysisName: None   })
    logLikelihoodBest.update({analysisName: -np.inf})
for analysis in analyses:
    analysis.update({'isBest': {'global': False}})
    for analysisName in analysisNames.keys():
        analysis['isBest'].update({analysisName: False})
    logLikelihoodGlobal = 0.0
    for analysisName in analysisNames.keys():
        groupAnalysis = analysis['group'][analysisName]
        if "logLikelihood" in groupAnalysis.attrs.keys():
            logLikelihood = groupAnalysis.attrs['logLikelihood']
            analysis['logLikelihood'] = logLikelihood
            logLikelihoodGlobal += logLikelihood
            if logLikelihood > logLikelihoodBest[analysisName]:
                logLikelihoodBest[analysisName] = logLikelihood
                analysisBest[analysisName] = analysis
        else:
            analysis['logLikelihood'] = -np.inf
    if "logLikelihood" in analysis['group'].attrs.keys():
        logLikelihoodGlobal = analysis['group'].attrs['logLikelihood']
        analysis['logLikelihoodGlobal'] = logLikelihoodGlobal
    else:
        analysis['logLikelihoodGlobal'] = -np.inf
    if logLikelihoodGlobal > logLikelihoodBest['global']:
        logLikelihoodBest['global'] = logLikelihoodGlobal
        analysisBest['global'] = analysis
analysisBest['global']['isBest'].update({'global': True})
for analysisName in analysisNames.keys():
    analysisBest[analysisName]['isBest'].update({analysisName: True})

# Report on likelihoods and models.
if not singleModel:
    reports.write("-> Report on best fit models\n")
    names = list(analysisNames.keys())
    names.insert(0,'global')
    for analysis in analyses:
        bestFor = []
        for name in names:
            if analysis['isBest'][name]:
                bestFor.append(name)
        if len(bestFor) > 0:
            stateVector    =      analysis['group']['simulationState'][:]
            parameterNames = list(analysis['group']['parameterNames' ][:])
            for i in range(len(parameterNames)):
                parameterNames[i] = parameterNames[i].decode('utf-8')
            lengthMaximum = len(max(parameterNames,key=len))
            reports.write(" -> "+analysis['filename']+" : "+analysis['groupname']+" - is best for:\n")
            reports.write("  -> "+", ".join(bestFor)+"\n")
            reports.write("  -> State:\n")
            for i in range(len(stateVector)):
                reports.write("      "+(" "*(lengthMaximum-len(parameterNames[i])))+parameterNames[i]+" = %+12.6e\n" % stateVector[i])
        
# Report on update parameter ranges.
logLikelihoodAccept = -np.inf
if not singleModel:
    reports.write("-> Report on updated parameter ranges\n")
    # Get parameter names.
    parameterNames = list(analysisBest['global']['group']['parameterNames'][:])
    for i in range(len(parameterNames)):
        parameterNames[i] = parameterNames[i].decode('utf-8')
    lengthMaximum = len(max(parameterNames,key=len))
    # Load the config file if available.
    if args.configFileName:
        tree             = ET.parse(args.configFileName)
        activeParameters = tree.findall(".//modelParameter[@value='active']")
    # Find accepted ranges.
    names = [ 'global' ]
    for name in names:
        # Find the best likelihood, and offset by our acceptance factor.
        logLikelihoodAccept = logLikelihoodBest[name]*(1.0+args.factorLogLikelihoodAccept)
        # Gather states and likelihoods across all analyses.
        logLikelihoods = []
        states = np.column_stack([analysis['group']['simulationState'][:] for analysis in analyses])
        for analysis in analyses:
            logLikelihood = -np.inf
            if name == "global":
                if "logLikelihood" in analysis['group'].attrs.keys():
                    logLikelihood = analysis['group'].attrs['logLikelihood']
            else:
                groupAnalysis = analysis['group'][name]
                if "logLikelihood" in groupAnalysis.attrs.keys():
                    logLikelihood = groupAnalysis.attrs['logLikelihood']
            logLikelihoods = np.append(logLikelihoods,logLikelihood)
        # Iterate over parameters, find the change in likelihood one step before and after the maximum likelihood.
        for i in range(states.shape[0]):
            order = np.argsort(states[i])
            j     = np.argmax(logLikelihoods[order])
            if j >            0:
                logLikelihoodAccept = np.min([logLikelihoodAccept,logLikelihoods[order][j-1]])
            if j < len(order)-1:
                logLikelihoodAccept = np.min([logLikelihoodAccept,logLikelihoods[order][j+1]])
        # Find accepted states.
        accepted = logLikelihoods >= logLikelihoodAccept
        # Check if any best fit models were not accepted.
        for analysis, analysisAccepted in zip(analyses,accepted):
            if not analysisAccepted:
                bestFor = []
                for comparison, isBest in analysis['isBest'].items():
                    if isBest:
                        bestFor.append(comparison)
                if len(bestFor) > 0:
                    reports.write("WARNING: "+analysis['filename']+" : "+analysis['groupname']+" was not accepted but is best for\n")
                    reports.write("  -> "+", ".join(bestFor)+"\n")
        # Iterate over parameters, finding the range which are above the acceptance likelihood.
        for i in range(states.shape[0]):
            statesAccepted       = np.sort(states[i][accepted])
            stateAcceptedMinimum = np.min(statesAccepted)
            stateAcceptedMaximum = np.max(statesAccepted)
            reports.write(" -> "+(" "*(lengthMaximum-len(parameterNames[i])))+parameterNames[i]+": %+12.6e — %+12.6e\n" % ( stateAcceptedMinimum, stateAcceptedMaximum ))
            if args.configFileName:
                for activeParameter in activeParameters:
                    activeParameterName = activeParameter.find('name').attrib['value']
                    if activeParameterName == parameterNames[i]:
                        prior        = activeParameter.find('distributionFunction1DPrior')
                        priorMinimum = prior          .find('limitLower')
                        priorMaximum = prior          .find('limitUpper')
                        priorMinimum.set('value',str(stateAcceptedMinimum))
                        priorMaximum.set('value',str(stateAcceptedMaximum))
    if args.configFileName:
        with open(args.outputDirectory+"/"+os.path.basename(args.configFileName), 'wb') as f:
            tree.write(f)
                    
# Iterate over analyses making plots.
reports.write("-> Generating plots\n")
for analysisName in analysisNames.keys():
    reports.write(" -> "+analysisName+"\n")
    targetPlotted = False
    axesCreated   = False
    for analysis in analyses:
        if singleModel and "logLikelihood" in analysis['group'].attrs.keys():
            reports.write("   -> logℒ = "+str(analysis['group'].attrs['logLikelihood'])+"\n")
        groupAnalysis = analysis['group'][analysisName]
        if "type" in groupAnalysis.attrs.keys():
            if groupAnalysis.attrs['type'] == b"function1D":
                datasets = {
                    'xDataset'         : {'required': True },
	            'yDataset'         : {'required': True },
	            'yDatasetTarget'   : {'required': False},
	            'yCovariance'      : {'required': False},
	            'yCovarianceTarget': {'required': False},
	            'yErrorLower'      : {'required': False},
	            'yErrorUpper'      : {'required': False},
	            'yErrorLowerTarget': {'required': False},
	            'yErrorUpperTarget': {'required': False}
                }
                for datasetName, dataset in datasets.items():
                    if datasetName in groupAnalysis.attrs.keys():
                        dataset['data'] = groupAnalysis[groupAnalysis.attrs[datasetName]][:]
                    else:
                        if dataset['required']:
                            sys.exit("Error: attribute '"+datasetName+"' is missing from analysis '"+analysisName+"' but is required.")
                for phase in range(5):
                    if not axesCreated:
                        figure, axes = plt.subplots(tight_layout=True)
                        if groupAnalysis.attrs['yAxisIsLog'] == 1:
                            axes.set_yscale('log')
                        if groupAnalysis.attrs['xAxisIsLog'] == 1:
                            axes.set_xscale('log')
                        # Replace double backslash with single. Note that we have to escape "\" twice here - once for Python and once for the re engine.
                        xAxisLabel  = re.sub("\\\\\\\\","\\\\",groupAnalysis.attrs['xAxisLabel' ].decode('utf-8'))
                        yAxisLabel  = re.sub("\\\\\\\\","\\\\",groupAnalysis.attrs['yAxisLabel' ].decode('utf-8'))
                        description = re.sub("\\\\\\\\","\\\\",groupAnalysis.attrs['description'].decode('utf-8'))
                        # matplotlib seems to not know "\hbox", so replace it with #\mathrm here.
                        xAxisLabel  = re.sub("\\\\hbox","\\\\mathrm",xAxisLabel )
                        yAxisLabel  = re.sub("\\\\hbox","\\\\mathrm",yAxisLabel )
                        description = re.sub("\\\\hbox","\\\\mathrm",description)
                        # matplotlib seems to not know "\le", so replace it with the corresponding Unicode symbol here.
                        description = re.sub("\\\\le","≤",description)
                        axes.set_xlabel(xAxisLabel )
                        axes.set_ylabel(yAxisLabel )
                        axes.set_title (description)
                        axesCreated = True
                    if not targetPlotted and phase == 2:
                        # Plot the target dataset.
                        targetPlotted = True
                        if 'data' in datasets['yDatasetTarget']:
                            haveErrorBars = False
                            y = datasets['yDatasetTarget']['data']
                            if 'data' in datasets['yErrorLowerTarget']:
                                errors = np.concatenate((datasets['yErrorLowerTarget']['data'], datasets['yErrorUpperTarget']['data']), axis=1)
                                yLower = y-datasets['yErrorLowerTarget']['data']
                                yUpper = y+datasets['yErrorUpperTarget']['data']
                                haveErrorBars = True
                            elif 'data' in datasets['yCovarianceTarget']:
                                errors = np.sqrt(np.diagonal(datasets['yCovarianceTarget']['data']))
                                yLower = y-errors
                                yUpper = y+errors
                                haveErrorBars = True
                            else:
                                yLower = y
                                yUpper = y
                            nonZero = np.nonzero(datasets['yDatasetTarget']['data'])
                            if haveErrorBars:
                                axes.errorbar(datasets['xDataset']['data'][nonZero],datasets['yDatasetTarget']['data'][nonZero],yerr=errors[nonZero],fmt='none',ecolor='#4c3af2',zorder=phase*10)
                                haveErrorBars = False
                            axes    .plot    (datasets['xDataset']['data'][nonZero],datasets['yDatasetTarget']['data'][nonZero],marker='o',linestyle='',markerfacecolor='#4c3af2',markeredgecolor='#2a0d83',zorder=phase*10,label=groupAnalysis.attrs['targetLabel'].decode('utf-8'))
                            # If only the target data range is to be shown, compute and set a suitable y-axis range here.
                            if not args.showFullRange:
                                if groupAnalysis.attrs['yAxisIsLog'] == 1:
                                    yMinimum  = np.min(yLower,initial=+np.inf,where=yLower > 0.0)
                                    yMaximum  = np.max(yUpper,initial=-np.inf,where=yUpper > 0.0)
                                    yRange    = yMaximum/yMinimum
                                    yMinimum /= np.power(yRange,0.05)
                                    yMaximum *= np.power(yRange,0.05)                              
                                else:
                                    yMinimum  = np.min(yLower,initial=+np.inf,where=y != 0.0)
                                    yMaximum  = np.max(yUpper,initial=-np.inf,where=y != 0.0)
                                    yRange    = yMaximum-yMinimum
                                    yMinimum -= 0.05*yRange
                                    yMaximum += 0.05*yRange
                                axes.set_ylim([yMinimum,yMaximum])
                            # Set a suitable x-axis range.
                            x        = datasets['xDataset']['data']
                            xMinimum = np.min(x)
                            xMaximum = np.max(x)
                            if groupAnalysis.attrs['xAxisIsLog'] == 1:
                                xRange    = xMaximum/xMinimum
                                xMinimum /= np.power(xRange,0.05)
                                xMaximum *= np.power(xRange,0.05)
                            else:
                                xRange    = xMaximum-xMinimum
                                xMinimum -= 0.05*xRange
                                xMaximum += 0.05*xRange
                            axes.set_xlim([xMinimum,xMaximum])
                    # Extract error bar information.
                    if 'data' in datasets['yErrorLower']:
                        errors = np.concatenate((datasets['yErrorLower']['data'], datasets['yErrorUpper']['data']), axis=1)
                        haveErrorBars = True
                    elif 'data' in datasets['yCovariance']:
                        errors = np.sqrt(np.diagonal(datasets['yCovariance']['data']))
                        haveErrorBars = True
                    # Detect best-fit analyses and set colors and styles appropriately.
                    if analysis['isBest'][analysisName]:
                        faceColor = '#ffc70f'
                        edgeColor = '#ff7d00'
                        marker    = 'o'
                        lineStyle = ''
                        label     = "Galacticus (best)"
                    elif analysis['isBest']['global']:
                        faceColor = '#24ab24'
                        edgeColor = '#246224'
                        marker    = 'o'
                        lineStyle = ''
                        label     = "Galacticus (global best)"
                    elif analysis['logLikelihoodGlobal'] >= logLikelihoodAccept:
                        faceColor = '#cccccc'
                        edgeColor = '#bbbbbb'
                        marker    = ''
                        lineStyle = 'solid'
                        label     = None
                    else:
                        faceColor = '#eeeeee'
                        edgeColor = '#dddddd'
                        marker    = ''
                        lineStyle = 'solid'
                        label     = None
                    # Determine which datasets to plot:
                    #  phase 0: analyses that are not a the best fit, and are not close to the best fit
                    #  phase 1: analyses that are not a the best fit, but are not close to the best fit
                    #  phase 2: nothing (this is for the target dataset)
                    #  phase 3: analyses that are the global (but not local) best fit
                    #  phase 4: the local best fit
                    # This ensures that z-ordering is as desired.
                    if   phase == 0:
                        doPlot = not analysis['isBest']['global'] and not analysis['isBest'][analysisName] and not analysis['logLikelihoodGlobal'] >= logLikelihoodAccept
                    elif phase == 1:
                        doPlot = not analysis['isBest']['global'] and not analysis['isBest'][analysisName] and     analysis['logLikelihoodGlobal'] >= logLikelihoodAccept
                    elif phase == 2:
                        doPlot = False
                    elif phase == 3:
                        doPlot =     analysis['isBest']['global'] and not analysis['isBest'][analysisName]
                    elif phase == 4:
                        doPlot =                                          analysis['isBest'][analysisName]
                    if doPlot:
                        nonZero = np.nonzero(datasets['yDataset']['data'])
                        if haveErrorBars and analysis['isBest'][analysisName]:
                            axes.errorbar(datasets['xDataset']['data'][nonZero],datasets['yDataset']['data'][nonZero],yerr=errors[nonZero],fmt='none',ecolor=faceColor,zorder=phase*10)
                        axes    .plot    (datasets['xDataset']['data'][nonZero],datasets['yDataset']['data'][nonZero],marker=marker,linestyle=lineStyle,color=faceColor,markerfacecolor=faceColor,markeredgecolor=edgeColor,zorder=phase*10,label=label)
        else:
            print("Warning: analysis '"+analysisName+"' has no 'type' attribute, so it can not be processed.")
    axes.legend(loc='upper center',bbox_to_anchor=(0.5,-0.05),fancybox=True,shadow=True,ncol=3)
    plt.savefig(args.outputDirectory+"/"+analysisName+'.pdf')  
    plt.clf()

# Close the report file.
reports.close()
