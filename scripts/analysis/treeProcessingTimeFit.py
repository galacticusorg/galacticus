#!/usr/bin/env python3
import numpy as np
import h5py
import sys
import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom

# Fit a quadratic log-log relation of tree processing time versus tree mass (and versus tree node count) using the raw per-tree
# timing data recorded by the `mergerTreeOperatorTreeProcessingTimer` operator (in the `metaData/treeTiming` group). Timing data
# from multiple Galacticus output files may be combined into a single fit. The resulting coefficients can be written to an XML file
# (consumable by the `file` metaTreeProcessingTime class) and/or written back into a Galacticus HDF5 file (in the same format that
# the timer operator itself writes, so that the file can subsequently be used directly by the `file` class).
#
# This is the "middle step" that turns recorded timing data into a usable cost model. The Fortran operator performs the same fit
# in-code at the end of a run; this script exists to combine data across multiple runs, or to (re)fit after the fact.
#
# Andrew Benson (10-July-2026)

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='treeProcessingTimeFit.py',description='Fit tree processing time versus mass and node count from recorded timing data.')
parser.add_argument('timingFiles',nargs='+',help='One or more Galacticus output files containing recorded tree timing data in the "metaData/treeTiming" group.')
parser.add_argument('--outputXML' ,default=None,help='If given, write the mass-based fit coefficients to this XML file (consumable by the "file" metaTreeProcessingTime class).')
parser.add_argument('--outputHDF5',default=None,help='If given, write the fit results (coefficients, residual, range) back into this Galacticus HDF5 file, in the "metaData/treeTiming" group.')
parser.add_argument('--degree'    ,type=int,default=2,help='The maximum polynomial degree to fit (default: 2).')
args = parser.parse_args()

# Accumulate the raw timing data across all supplied files.
mass       = np.array([],dtype='float64')
countNodes = np.array([],dtype='float64')
timeProcess= np.array([],dtype='float64')
for fileName in args.timingFiles:
    with h5py.File(fileName,'r') as file:
        if not "metaData" in file or not "treeTiming" in file["metaData"]:
            sys.exit('file "'+fileName+'" contains no "metaData/treeTiming" group - was the "mergerTreeOperatorTreeProcessingTimer" operator active in the run that produced it?')
        timing        = file["metaData/treeTiming"]
        massFile      = timing["treeMass"     ][:].astype('float64')
        countNodesFile= timing["countNodes"   ][:].astype('float64')
        # The total processing time is the sum of construction and evolution times.
        timeFile      = timing["timeConstruct"][:].astype('float64')+timing["timeEvolve"][:].astype('float64')
        mass          = np.append(mass       ,massFile      )
        countNodes    = np.append(countNodes ,countNodesFile)
        timeProcess   = np.append(timeProcess,timeFile      )

# Perform a least-squares polynomial fit of log10(time) versus log10(predictor), reducing the degree automatically if the data do
# not support the requested degree. Returns the three coefficients (padded with zeros), the RMS residual (in dex), and the range of
# the predictor.
def fit(predictor,time):
    select = (predictor > 0.0) & (time > 0.0)
    x      = np.log10(predictor[select])
    y      = np.log10(time     [select])
    if x.size < 1:
        return None
    degree       = min(args.degree,x.size-1,2)
    coefficients = np.zeros(3)
    while degree >= 0:
        try:
            fitCoefficients = np.polyfit(x,y,degree)
            break
        except (np.linalg.LinAlgError,ValueError):
            degree -= 1
    else:
        return None
    # np.polyfit returns highest-order coefficient first; reverse to give [C0, C1, C2].
    for i in range(degree+1):
        coefficients[i] = fitCoefficients[degree-i]
    model    = np.zeros_like(x)
    for i in range(3):
        model += coefficients[i]*x**i
    if x.size > degree+1:
        residual = np.sqrt(np.sum((y-model)**2)/(x.size-degree-1))
    else:
        residual = 0.0
    return coefficients,residual,np.array([np.min(predictor[select]),np.max(predictor[select])])

fitMass       = fit(mass      ,timeProcess)
fitCountNodes = fit(countNodes,timeProcess)
if fitMass is None:
    sys.exit('unable to fit - no valid timing data found')
coefficientsMass      ,residualMass      ,rangeMass       = fitMass
coefficientsCountNodes,residualCountNodes,rangeCountNodes = fitCountNodes if fitCountNodes is not None else (np.zeros(3),0.0,np.array([0.0,0.0]))

# Report the fit to standard output.
print("Fitted %d trees from %d file(s)." % (timeProcess.size,len(args.timingFiles)))
print("Mass-based fit       : log10(t/s) = %+.4e %+.4e log10(M) %+.4e log10(M)^2 ; residual = %.3f dex ; mass range = [%.3e, %.3e] Msun" % (coefficientsMass[0],coefficientsMass[1],coefficientsMass[2],residualMass,rangeMass[0],rangeMass[1]))
print("Node-count-based fit : log10(t/s) = %+.4e %+.4e log10(N) %+.4e log10(N)^2 ; residual = %.3f dex ; node count range = [%d, %d]" % (coefficientsCountNodes[0],coefficientsCountNodes[1],coefficientsCountNodes[2],residualCountNodes,int(rangeCountNodes[0]),int(rangeCountNodes[1])))

# Optionally write the mass-based fit coefficients to an XML file.
if args.outputXML is not None:
    timing = ET.Element("timing")
    fitXML = ET.SubElement(timing,"fit")
    for i in range(3):
        coefficient      = ET.SubElement(fitXML,"coefficient")
        coefficient.text = "%.8e" % coefficientsMass[i]
    xmlString = minidom.parseString(ET.tostring(timing)).toprettyxml(indent=" ")
    with open(args.outputXML,"w") as xmlFile:
        xmlFile.write(xmlString)
    print('Wrote mass-based fit coefficients to XML file "'+args.outputXML+'".')

# Optionally write the fit results back into a Galacticus HDF5 file.
if args.outputHDF5 is not None:
    with h5py.File(args.outputHDF5,"a") as file:
        timing = file.require_group("metaData/treeTiming")
        for name,value in (
            ("fitCoefficientMass"      ,coefficientsMass      ),
            ("fitResidualMass"         ,np.array([residualMass])),
            ("fitRangeMass"            ,rangeMass             ),
            ("fitCoefficientCountNodes",coefficientsCountNodes),
            ("fitResidualCountNodes"   ,np.array([residualCountNodes])),
            ("fitRangeCountNodes"      ,rangeCountNodes       ),
        ):
            if name in timing:
                del timing[name]
            timing.create_dataset(name,data=value)
    print('Wrote fit results to HDF5 file "'+args.outputHDF5+'".')
