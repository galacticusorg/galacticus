#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test postprocessing of galaxies read back from a Galacticus output file, using the `postprocessOutputs` task. A lightcone model
# is run which outputs star formation histories and random inclinations, along with the SEDs (unattenuated, and attenuated by
# both an inclination-independent dust model and by the inclination-dependent dust compendium), emission line luminosities
# (unattenuated, and attenuated by the dust compendium), and stellar masses which postprocessing must then reproduce. The galaxies are postprocessed twice: once in full, writing each galaxy as it is output and
# including the properties of the rebuilt node components; and once for a selection of galaxies (unordered, with repeats, in
# small batches) writing only the extracted properties, with output buffered.
# Andrew Benson (16-November-2026)

outputPath = "outputs/test-postprocess-outputs"

def runModel(label,parameterFile,changesFile=None):
    # Run a model, and check for completion and for errors. Returns true on success.
    print("Running model ("+label+")...")
    logFileName = outputPath+"/"+label+".log"
    log         = open(logFileName,"w")
    command     = "cd ..; ./Galacticus.exe testSuite/parameters/"+parameterFile
    if changesFile is not None:
        command += " testSuite/parameters/"+changesFile
    status      = subprocess.run(command,stdout=log,stderr=log,shell=True)
    log.close()
    print("...done ("+str(status)+")")
    if status.returncode != 0:
        print("FAILED: model run ("+label+"):")
        subprocess.run("cat "+logFileName,shell=True)
        return False
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" "+logFileName,shell=True)
    if status.returncode == 0:
        print("FAILED: model run ("+label+") (errors):")
        subprocess.run("cat "+logFileName,shell=True)
        return False
    print("SUCCESS: model run ("+label+")")
    return True

def compare(label,nodesOriginal,nodesPostprocessed,names):
    # Compare datasets in the postprocessed output with those of the original model, matching galaxies through the index of each
    # in the original output. Returns true if all agree.
    indices = nodesPostprocessed['nodeIndexPostprocessSource'][:]
    success = True
    for name in names:
        if name not in nodesPostprocessed:
            print("FAILED: "+label+": dataset '"+name+"' is missing from the postprocessed output")
            success = False
            continue
        valuesOriginal      = nodesOriginal     [name][:][indices]
        valuesPostprocessed = nodesPostprocessed[name][:]
        if valuesOriginal.shape != valuesPostprocessed.shape:
            print("FAILED: "+label+": dataset '"+name+"' has shape "+str(valuesPostprocessed.shape)+" but expected "+str(valuesOriginal.shape))
            success = False
        elif not np.allclose(valuesPostprocessed,valuesOriginal,rtol=1.0e-10,atol=0.0):
            difference = np.max(np.abs(valuesPostprocessed-valuesOriginal)/np.maximum(np.abs(valuesOriginal),np.finfo(float).tiny))
            print("FAILED: "+label+": dataset '"+name+"' differs from the original (maximum relative difference "+str(difference)+")")
            success = False
    if success:
        print("SUCCESS: "+label+": "+str(len(names))+" datasets agree with the original model")
    return success

def compareAttributes(label,fileName):
    # Check that the output group written has the type and time attributes of the output from which its galaxies were read.
    groupOriginal      = model                       ["Lightcone/Output1"]
    groupPostprocessed = h5py.File(fileName,"r")["Lightcone/Output1"]
    success            = True
    for name in ( "outputType", "outputTime", "outputExpansionFactor", "outputComovingDistance" ):
        if name not in groupPostprocessed.attrs or groupPostprocessed.attrs[name] != groupOriginal.attrs[name]:
            print("FAILED: "+label+": attribute '"+name+"' of the output group is "+str(groupPostprocessed.attrs.get(name))+" but expected "+str(groupOriginal.attrs[name]))
            success = False
    if success:
        print("SUCCESS: "+label+": output group attributes match the original model")

# Run the model whose galaxies are to be postprocessed.
subprocess.run("mkdir -p "+outputPath,shell=True)
if not runModel("model","test-postprocess-outputs-model.xml","test-postprocess-outputs-modelChanges.xml"):
    sys.exit(0)
model         = h5py.File(outputPath+"/model.hdf5","r")
nodesOriginal = model["Lightcone/Output1/nodeData"]
countGalaxies = nodesOriginal["time"].size
# Names of the properties which postprocessing must reproduce.
namesTarget = [
    "massStellarTotal",
    "inclination",
    "diskStellarSED:identity",
    "spheroidStellarSED:identity",
    "diskStellarSED:identity:dustAttenuated:stellarMassRedshift",
    "spheroidStellarSED:identity:dustAttenuated:stellarMassRedshift",
    "diskStellarSED:identity:dustAttenuated:atlasCompendium",
    "spheroidStellarSED:identity:dustAttenuated:atlasCompendium"
]
for component in ( "Disk", "Spheroid" ):
    for line in ( "balmerAlpha6565", "oxygenIII5008" ):
        for suffix in ( "", ":dustAttenuated:atlasCompendium" ):
            namesTarget.append("luminosityEmissionLine"+component+":"+line+suffix)
# Names of the properties of the rebuilt node components. (The time of last isolation is not included: nodes are rebuilt in
# isolation, so that it is reported as the current time.)
namesComponent = [
    "basicMass"
]
for component in ( "disk", "spheroid" ):
    for property in ( "MassStellar", "MassGas", "AbundancesStellarMetals", "AbundancesGasMetals", "Radius", "AngularMomentum", "Velocity" ):
        namesComponent.append(component+property)
# Check that the model is a meaningful test: it must contain galaxies with non-zero SEDs, galaxies with no star formation
# history in one component (which exercises the reading of empty variable-length entries), and galaxies whose disks are
# significantly attenuated by the dust compendium (so that an error in the inclination, or in the gas properties, used in
# postprocessing would be seen), and galaxies with emission lines.
countNonZero      = np.count_nonzero(np.sum(nodesOriginal["diskStellarSED:identity"][:],axis=1))
lengthHistories   = np.array([entry.size for entry in nodesOriginal["diskStarFormationHistoryMass"][:]])
indicesEmpty      = np.nonzero(lengthHistories == 0)[0]
sedUnattenuated   = np.sum(nodesOriginal["diskStellarSED:identity"                                       ][:],axis=1)
sedAttenuated     = np.sum(nodesOriginal["diskStellarSED:identity:dustAttenuated:atlasCompendium"][:],axis=1)
countAttenuated   = np.count_nonzero(sedAttenuated < 0.9*sedUnattenuated)
countLines        = np.count_nonzero(nodesOriginal["luminosityEmissionLineDisk:balmerAlpha6565"][:] > 0.0)
if countGalaxies == 0 or countNonZero == 0 or indicesEmpty.size == 0 or countAttenuated == 0 or countLines == 0:
    print("FAILED: model does not provide a meaningful test (galaxies: "+str(countGalaxies)+"; with non-zero disk SEDs: "+str(countNonZero)+"; with empty disk star formation histories: "+str(indicesEmpty.size)+"; with disks attenuated by the dust compendium: "+str(countAttenuated)+"; with disk H-alpha emission: "+str(countLines)+")")
    sys.exit(0)

# Write the selection of galaxies to postprocess: an unordered set, with repeats, which includes a galaxy with an empty star
# formation history.
randomGenerator = np.random.default_rng(seed=1066)
selection       = randomGenerator.integers(0,countGalaxies,size=countGalaxies//2)
selection       = np.concatenate((selection,selection[0:5],[indicesEmpty[0]]))
randomGenerator.shuffle(selection)
with h5py.File(outputPath+"/selection.hdf5","w") as selectionFile:
    selectionFile.create_dataset("Lightcone/Output1/selection",data=selection)

# Postprocess all galaxies.
if runModel("postprocessAll","test-postprocess-outputs-all.xml"):
    nodesPostprocessed = h5py.File(outputPath+"/postprocessedAll.hdf5","r")["Lightcone/Output1/nodeData"]
    indices            = nodesPostprocessed["nodeIndexPostprocessSource"][:]
    if np.array_equal(np.sort(indices),np.arange(countGalaxies)):
        print("SUCCESS: all: each galaxy was postprocessed once")
    else:
        print("FAILED: all: galaxies postprocessed do not match those in the original model")
    compare("all: target properties"   ,nodesOriginal,nodesPostprocessed,namesTarget   )
    compareAttributes("all",outputPath+"/postprocessedAll.hdf5")
    compare("all: component properties",nodesOriginal,nodesPostprocessed,namesComponent)

# Postprocess the selection of galaxies.
if runModel("postprocessSelection","test-postprocess-outputs-selection.xml"):
    nodesPostprocessed = h5py.File(outputPath+"/postprocessedSelection.hdf5","r")["Lightcone/Output1/nodeData"]
    indices            = nodesPostprocessed["nodeIndexPostprocessSource"][:]
    if np.array_equal(np.sort(indices),np.sort(selection)):
        print("SUCCESS: selection: the selected galaxies were postprocessed")
    else:
        print("FAILED: selection: galaxies postprocessed do not match those selected")
    compare("selection: target properties",nodesOriginal,nodesPostprocessed,namesTarget)
    compareAttributes("selection",outputPath+"/postprocessedSelection.hdf5")
    # With output of component properties switched off, only the extracted properties (and the descriptors of their columns)
    # should be present.
    namesExpected = set(namesTarget+["nodeIndexPostprocessSource"])
    namesUnexpected = sorted([name for name in nodesPostprocessed.keys() if name not in namesExpected and not name.endswith(("Columns","ColumnValues"))])
    if len(namesUnexpected) == 0:
        print("SUCCESS: selection: only extracted properties were output")
    else:
        print("FAILED: selection: unexpected datasets were output: "+", ".join(namesUnexpected))

# Postprocess the selection of galaxies through the Python interface. The parameters are reconstructed from the model's own
# output file, and the property extractors are those of the model itself, taken from that same record.
try:
    from Galacticus.Parameters.extract import parameters_from_hdf5
    from galacticus_launcher import postprocess as launcherPostprocess
    from galacticus_launcher import paths as launcherPaths
    haveLauncher = True
except ImportError as error:
    # A missing launcher dependency means only that this sub-check can not be run - the marker word is deliberately not used.
    print("The Python interface was not tested, as it could not be imported: "+str(error))
    haveLauncher = False
if haveLauncher:
    # The launcher must use the executable built in this tree. Without GALACTICUS_EXEC_PATH and GALACTICUS_DATA_PATH it would
    # instead resolve to a managed install, and download a released executable to test.
    install = launcherPaths.resolve()
    if install.managed:
        print("The Python interface was not tested, as GALACTICUS_EXEC_PATH and GALACTICUS_DATA_PATH are not both set")
        haveLauncher = False
if haveLauncher:
    print("Running model (postprocessPython)...")
    extractorsModel = parameters_from_hdf5(outputPath+"/model.hdf5").find("nodePropertyExtractor")
    extractors      = [extractor for extractor in extractorsModel.findall("nodePropertyExtractor") if extractor.get("value") in ( "massStellar", "inclination", "dustAttenuation" )]
    try:
        launcherPostprocess.postprocess(
            outputPath+"/model.hdf5",
            extractors,
            outputPath+"/postprocessedPython.hdf5",
            selection      ={1: selection},
            outputs_group  ="Lightcone",
            count_batch    =16,
            parameter_file =outputPath+"/postprocessPython.xml",
            log_file       =outputPath+"/postprocessPython.log",
            install        =install
        )
        print("SUCCESS: model run (postprocessPython)")
        nodesPostprocessed = h5py.File(outputPath+"/postprocessedPython.hdf5","r")["Lightcone/Output1/nodeData"]
        indices            = nodesPostprocessed["nodeIndexPostprocessSource"][:]
        if np.array_equal(np.sort(indices),np.sort(selection)):
            print("SUCCESS: Python: the selected galaxies were postprocessed")
        else:
            print("FAILED: Python: galaxies postprocessed do not match those selected")
        compare("Python: target properties",nodesOriginal,nodesPostprocessed,namesTarget)
        compareAttributes("Python",outputPath+"/postprocessedPython.hdf5")
    except (RuntimeError, ValueError) as error:
        print("FAILED: model run (postprocessPython): "+str(error))
        subprocess.run("cat "+outputPath+"/postprocessPython.log",shell=True)
    # An invalid selection must be rejected before Galacticus is run.
    try:
        launcherPostprocess.write_selection(outputPath+"/selectionInvalid.hdf5",np.array([countGalaxies]),outputs_group="Lightcone",catalog=outputPath+"/model.hdf5")
        print("FAILED: Python: an out-of-range selection was not rejected")
    except ValueError:
        print("SUCCESS: Python: an out-of-range selection was rejected")

sys.exit(0)
