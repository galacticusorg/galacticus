#!/usr/bin/env python3
# Plot expected evolution of disk/spheroid gas/stellar masses for a model with non-negative evolution.
import numpy as np
import h5py
import subprocess
import sys
import os
import argparse

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='test-enforceNonNegativity.py',description='Test non-negative evolution.')
parser.add_argument('--makeplots', action='store_true',help='make plots showing mass evolution')

# Create output path.
try:
    os.mkdir("outputs/enforceNonNegativity")
except FileExistsError:
    pass

# Inintialize status.
testStatus = "SUCCESS: no negative values found when enforcing"

# Iterate over model variations.
for variation in ( "ODETol0.1", "ODETol0.01", "ODETol0.001", "ODETol0.0001", "ODETol1e-6", "ODETol1e-9", "EnforceNonNegativityODETol0.1", "EnforceNonNegativityODETol0.01", "EnforceNonNegativityODETol0.001", "EnforceNonNegativityODETol0.0001", "EnforceNonNegativityODETol1e-6", "EnforceNonNegativityODETol1e-9" ):

    # Run the model.
    status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/enforceNonNegativity/rapidDepletion.xml testSuite/parameters/enforceNonNegativity/rapidDepletion"+variation+".xml",shell=True)
    if status.returncode != 0:
        print("FAILED: model failed to run")
        sys.exit()
        
    # Begin a figure showing evolution of masses with time.
    if args.makeplots:
        import matplotlib.pyplot as plt
        figure, axes = plt.subplots(tight_layout=True)
        axes.set_yscale('log'                         )
        axes.set_xlabel('$t/$Gyr'                     )
        axes.set_ylabel('$M/\mathrm{M}_\odot$'        )
        axes.set_title (f'ODE evolution ({variation})')

    # Initialize arrays that will hold times and masses.
    time                      = np.array([])
    massGasSpheroid           = np.array([])
    massStellarSpheroid       = np.array([])
    massMetalsGasSpheroid     = np.array([])
    massMetalsStellarSpheroid = np.array([])
    
    # Access the model.
    model = h5py.File("outputs/enforceNonNegativity/rapidDepletion"+variation+".hdf5","r")
    ## Extract parameters needed for analytic solution.
    parameters             = model     ['Parameters'                                           ]
    recycledFraction       = parameters['stellarPopulation'                                    ].attrs['recycledFraction'    ]
    metalYield             = parameters['stellarPopulation'                                    ].attrs['metalYield'          ]
    timescaleStarFormation = parameters['starFormationRateSpheroids/starFormationTimescale'    ].attrs['timescale'           ]
    massLoading            = parameters['nodeOperator/nodeOperator[10]/stellarFeedbackOutflows'].attrs['fraction'            ]
    enforcing              = parameters['mergerTreeNodeEvolver'                                ].attrs['enforceNonNegativity'].decode("utf-8")
    ## Extract meta-data on number of evaluations required.
    countEvaluations       = np.sum(model['metaData']['evolverProfiler']['evaluationCount'][:])
    ## Iterate over all outputs.
    outputs = model['Outputs']
    for outputName in outputs:
        output                    = outputs[outputName]
        nodes                     = output ['nodeData']
        time                      = np.append(time                     ,output.attrs['outputTime'                     ]   )
        massGasSpheroid           = np.append(massGasSpheroid          ,nodes       ['spheroidMassGas'                ][0])
        massStellarSpheroid       = np.append(massStellarSpheroid      ,nodes       ['spheroidMassStellar'            ][0])
        massMetalsGasSpheroid     = np.append(massMetalsGasSpheroid    ,nodes       ['spheroidAbundancesGasMetals'    ][0])
        massMetalsStellarSpheroid = np.append(massMetalsStellarSpheroid,nodes       ['spheroidAbundancesStellarMetals'][0])
        
    # Determine the correct time-ordering of model points.
    order = np.argsort(time)

    # Test for negative values if enforcing.
    if enforcing == "true":
        if np.any(massGasSpheroid           < 0.0):
            testStatus = "FAIL: negative gas mass detected"
        if np.any(massStellarSpheroid       < 0.0):
            testStatus = "FAIL: negative stellar mass detected"
        if np.any(massMetalsGasSpheroid     < 0.0):
            testStatus = "FAIL: negative gas metal mass detected"
        if np.any(massMetalsStellarSpheroid < 0.0):
            testStatus = "FAIL: negative stellar metal mass detected"
    
    # Evaluate analytic solutions (Cole et al.; 2000; https://ui.adsabs.harvard.edu/abs/2000MNRAS.319..168C; equations B2-B8).
    timescaleStarFormationEffective   = timescaleStarFormation/(1.0-recycledFraction+massLoading)
    timeFine                          = np.arange(time[order][0],time[order][-1],0.001)
    timeDelta                         = timeFine-timeFine[0]
    massGasSpheroidAnalytic           =                                                                         massGasSpheroid[order][0]*                                                     np.exp(-timeDelta/timescaleStarFormationEffective)
    massStellarSpheroidAnalytic       =            (1.0-recycledFraction)/(1.0-recycledFraction+massLoading)   *massGasSpheroid[order][0]*(1.0-                                                np.exp(-timeDelta/timescaleStarFormationEffective))
    massMetalsStellarSpheroidAnalytic = metalYield*(1.0-recycledFraction)/(1.0-recycledFraction+massLoading)**2*massGasSpheroid[order][0]*(1.0-(1.0+timeDelta/timescaleStarFormationEffective)*np.exp(-timeDelta/timescaleStarFormationEffective))
    massMetalsGasSpheroidAnalytic     = metalYield/(1.0-recycledFraction)*massStellarSpheroidAnalytic-(1.0-recycledFraction+massLoading)/(1.0-recycledFraction)*massMetalsStellarSpheroidAnalytic
    
    # Plot model results.
    if args.makeplots:
        axes.text(12.5,0.01,f'$N_\mathrm{{eval}}={countEvaluations}$')
        ## Gas mass.
        positive = massGasSpheroid           >= 0.0
        negative = massGasSpheroid           <  0.0
        axes.plot(time    [order][positive],+massGasSpheroid                  [order][positive],marker='o',linestyle='' ,label="$M_\mathrm{gas}$"    ,markerfacecolor="burlywood" ,markeredgecolor="burlywood" ,markersize=8)
        axes.plot(time    [order][negative],-massGasSpheroid                  [order][negative],marker='o',linestyle='' ,label=""                    ,markerfacecolor="burlywood" ,markeredgecolor="burlywood" ,markersize=2)
        axes.plot(timeFine                 , massGasSpheroidAnalytic                           ,marker='' ,linestyle='-',label=""                    ,color          ="burlywood"                                           )
        ## Stellar mass.
        positive = massStellarSpheroid       >= 0.0
        negative = massStellarSpheroid       <  0.0
        axes.plot(time    [order][positive],+massStellarSpheroid              [order][positive],marker='*',linestyle='' ,label="$M_\star$"           ,markerfacecolor="dodgerblue",markeredgecolor="dodgerblue",markersize=8)
        axes.plot(time    [order][negative],+massStellarSpheroid              [order][negative],marker='*',linestyle='' ,label=""                    ,markerfacecolor="dodgerblue",markeredgecolor="dodgerblue",markersize=2)
        axes.plot(timeFine                 , massStellarSpheroidAnalytic                       ,marker='' ,linestyle='-',label=""                    ,color          ="dodgerblue"                                          )
        ## Gas metal mass.
        positive = massMetalsGasSpheroid     >= 0.0
        negative = massMetalsGasSpheroid     <  0.0
        axes.plot(time    [order][positive],+massMetalsGasSpheroid            [order][positive],marker='d',linestyle='' ,label="$M_{Z,\mathrm{gas}}$",markerfacecolor="lightcoral",markeredgecolor="lightcoral",markersize=8)
        axes.plot(time    [order][negative],+massMetalsGasSpheroid            [order][negative],marker='d',linestyle='' ,label=""                    ,markerfacecolor="lightcoral",markeredgecolor="lightcoral",markersize=2)
        axes.plot(timeFine                 , massMetalsGasSpheroidAnalytic                     ,marker='' ,linestyle='-',label=""                    ,color          ="lightcoral"                                          )
        ## Stellar metal mass.
        positive = massMetalsStellarSpheroid >= 0.0
        negative = massMetalsStellarSpheroid <  0.0
        axes.plot(time    [order][positive],+massMetalsStellarSpheroid        [order][positive],marker='H',linestyle='' ,label="$M_{Z,\star}$"       ,markerfacecolor="plum"      ,markeredgecolor="plum"      ,markersize=8)
        axes.plot(time    [order][negative],+massMetalsStellarSpheroid        [order][negative],marker='H',linestyle='' ,label=""                    ,markerfacecolor="plum"      ,markeredgecolor="plum"      ,markersize=2)
        axes.plot(timeFine                 , massMetalsStellarSpheroidAnalytic                 ,marker='' ,linestyle='-',label=""                    ,color          ="plum"                                                )    
        ## Finalize the plot.
        axes.legend(loc='upper center',bbox_to_anchor=(0.4,0.12),fancybox=True,shadow=True,ncol=4)
        plt.savefig('outputs/enforceNonNegativity/massEvolution'+variation+'.pdf')
        plt.clf()

# Report status.
print(testStatus)
