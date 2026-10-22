#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np
import lxml.etree as ET

# Check that list property extractors (which output variable-length datasets) with more than one element are handled
# correctly by the standard merger tree outputter when used directly, rather than inside a "multi" property extractor. The
# outputter must reserve, name, and write one output buffer for each element of the list. Previously it reserved only one
# buffer per list extractor, so that a list extractor with more than one element overflowed the buffers and corrupted
# memory.
#
# The model is run three times: once with a double list extractor (radiusBlackHoles, with two elements) used directly, once
# with an integer list extractor (galaxyMergerTreeIndices, with two elements) used directly, and once with both inside a
# "multi" extractor. The datasets written by each direct extractor must match those written for it by the "multi" extractor.
# Claude (21-October-2026)

subprocess.run("mkdir -p outputs",shell=True)

# Construct the parameter files for each run.
parametersBase = ET.parse("parameters/outputListExtractorTopLevel.xml")
runs           = {
    "radiusBlackHoles"       : [ "radiusBlackHoles"                            ],
    "galaxyMergerTreeIndices": [ "galaxyMergerTreeIndices"                     ],
    "multi"                  : [ "radiusBlackHoles", "galaxyMergerTreeIndices" ]
}
datasets       = {
    "radiusBlackHoles"       : [ "radiusBlackHoles"         , "radialMigrationRateBlackHoles" ],
    "galaxyMergerTreeIndices": [ "galaxyMergerTreeNodeIndex", "galaxyMergerTreeCount"         ]
}
for run in runs:
    parameters = ET.ElementTree(ET.fromstring(ET.tostring(parametersBase.getroot())))
    root       = parameters.getroot()
    root.find("outputFileName").set("value","testSuite/outputs/outputListExtractorTopLevel_"+run+".hdf5")
    root.remove(root.find("nodePropertyExtractor"))
    extractor  = ET.SubElement(root,"nodePropertyExtractor")
    if run == "multi":
        extractor.set("value","multi")
        for name in runs[run]:
            ET.SubElement(extractor,"nodePropertyExtractor").set("value",name)
    else:
        extractor.set("value",runs[run][0])
    parameters.write("outputs/outputListExtractorTopLevel_"+run+".xml")

# Run the models and check for completion.
for run in runs:
    print("Running model ("+run+")...")
    logFileName = "outputs/test-output-list-extractor-top-level_"+run+".log"
    with open(logFileName,"w") as logFile:
        status = subprocess.run(
            "cd ..; ./Galacticus.exe testSuite/outputs/outputListExtractorTopLevel_"+run+".xml",
            stdout=logFile, stderr=logFile, shell=True
        )
    if status.returncode != 0:
        print("FAILED: model run ("+run+")")
        subprocess.run("cat "+logFileName,shell=True)
        sys.exit(0)
    errors = subprocess.run(
        "grep -q -i -e fatal -e aborted -e \"task failed\" "+logFileName,
        shell=True
    )
    if errors.returncode == 0:
        print("FAILED: model run ("+run+") (errors in log)")
        subprocess.run("cat "+logFileName,shell=True)
        sys.exit(0)
    print("SUCCESS: model run ("+run+")")

# Compare the datasets written by each direct extractor with those written by the "multi" extractor.
def nodeOrder(output):
    # Return the order of nodes in an output sorted by tree index, so that outputs can be compared regardless of the order in
    # which trees were written.
    treeIndex = output["mergerTreeIndex"     ][:]
    treeStart = output["mergerTreeStartIndex"][:]
    treeCount = output["mergerTreeCount"     ][:]
    if len(treeIndex) == 0:
        return np.array([],dtype=int)
    return np.concatenate([np.arange(treeStart[i],treeStart[i]+treeCount[i]) for i in np.argsort(treeIndex)]).astype(int)

modelMulti = h5py.File("outputs/outputListExtractorTopLevel_multi.hdf5","r")
for run in datasets:
    model = h5py.File("outputs/outputListExtractorTopLevel_"+run+".hdf5","r")
    for outputName in modelMulti["Outputs"]:
        if outputName not in model["Outputs"]:
            print("FAILED: output '"+outputName+"' is absent from the '"+run+"' model")
            sys.exit(0)
        output      = model     ["Outputs"][outputName]
        outputMulti = modelMulti["Outputs"][outputName]
        order       = nodeOrder(output     )
        orderMulti  = nodeOrder(outputMulti)
        for name in datasets[run]:
            if name not in output["nodeData"]:
                print("FAILED: dataset '"+name+"' is absent from output '"+outputName+"' of the '"+run+"' model")
                sys.exit(0)
            data      = output     ["nodeData"][name]
            dataMulti = outputMulti["nodeData"][name]
            if h5py.check_vlen_dtype(data.dtype) is None:
                print("FAILED: dataset '"+name+"' in output '"+outputName+"' is not a variable-length dataset")
                sys.exit(0)
            if data.shape != dataMulti.shape:
                print("FAILED: dataset '"+name+"' in output '"+outputName+"' has shape "+str(data.shape)+", but "+str(dataMulti.shape)+" when output by the 'multi' extractor")
                sys.exit(0)
            values      = data     [:][order     ]
            valuesMulti = dataMulti[:][orderMulti]
            if not all(np.array_equal(value,valueMulti) for value, valueMulti in zip(values,valuesMulti)):
                print("FAILED: dataset '"+name+"' in output '"+outputName+"' differs from that output by the 'multi' extractor")
                sys.exit(0)
            attributes      = dict(data     .attrs)
            attributesMulti = dict(dataMulti.attrs)
            if set(attributes) != set(attributesMulti) or any(not np.array_equal(np.asarray(attributes[key]),np.asarray(attributesMulti[key])) for key in attributes):
                print("FAILED: attributes of dataset '"+name+"' in output '"+outputName+"' differ from those output by the 'multi' extractor")
                sys.exit(0)
    print("SUCCESS: datasets output directly by '"+run+"' match those output by the 'multi' extractor")
sys.exit(0)
