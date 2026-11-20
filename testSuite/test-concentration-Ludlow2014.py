#!/usr/bin/env python3
import subprocess
import sys
import xml.etree.ElementTree as ET
import h5py
import numpy as np

# Check that the Ludlow et al. (2014) scale radius class uses its own default for the parameter C (853, from Fig. 3 of Ludlow
# et al. 2014), rather than inheriting the default of the Ludlow et al. (2016) class (400) from which it derives. The same small
# model is run three times: with C left at its default, with C=853 given explicitly, and with C=400 given explicitly. The first
# two must give the same concentrations, and the third must differ from them (which shows that the concentrations do depend
# on C in this model, so that the first comparison can fail). Only isolated halos (the roots of the trees) are compared, as the
# concentration affects the evolution of satellites, and so which satellites survive to be output. Models are run
# single-threaded, so that the trees are output in the same order, and the masses of the roots are required to agree between
# the runs (to a relative tolerance of 1e-6, far smaller than the differences between the masses of different roots), which
# confirms that the halos are paired correctly. Concentrations are compared with a small relative tolerance, as a model that
# tabulates sigma(M) afresh can differ at the 1e-7 level from one that reads the tabulation back.
# Andrew Benson, Claude (19-November-2026)

subprocess.run("mkdir -p outputs/", shell=True)

def makeModel(label,C):
    # Build a parameter file for the model from the Ludlow et al. (2016) concentration distribution test.
    tree       = ET.parse("parameters/concentrationDistributionLudlow2016.xml")
    parameters = tree.getroot()
    # Reduce the number of trees.
    parameters.find("mergerTreeBuildMasses/treesPerDecade").set("value","6")
    # Switch to the Ludlow et al. (2014) class, and set C (if required). The parameter f is left at its default.
    scaleRadius = parameters.find("darkMatterProfileScaleRadius/darkMatterProfileScaleRadius")
    scaleRadius.set("value","ludlow2014")
    for name in ( "C", "f" ):
        element = scaleRadius.find(name)
        if element is not None:
            scaleRadius.remove(element)
    if C is not None:
        element = ET.Element("C")
        element.set("value",C)
        scaleRadius.insert(0,element)
    parameters.find("outputFileName").set("value","testSuite/outputs/concentrationLudlow2014"+label+".hdf5")
    fileName = "outputs/concentrationLudlow2014"+label+".xml"
    tree.write(fileName)
    # Run the model.
    status = subprocess.run(f"cd ..; OMP_NUM_THREADS=1 ./Galacticus.exe testSuite/{fileName}",shell=True)
    if status.returncode != 0:
        print(f"FAILED: Ludlow2014 concentration model ({label}) failed to run")
        sys.exit(0)
    with h5py.File(f"outputs/concentrationLudlow2014{label}.hdf5","r") as model:
        nodeData      = model["Outputs/Output1/nodeData"]
        isolated      = nodeData["nodeIsIsolated"][:] == 1
        concentration = nodeData["concentration" ][:][isolated]
        massRoot      = nodeData["basicMass"     ][:][isolated]
    return massRoot, concentration

massDefault    , concentrationDefault  = makeModel("Default" ,None )
mass853        , concentration853      = makeModel("C853"    ,"853")
mass400        , concentration400      = makeModel("C400"    ,"400")

if len(concentrationDefault) == 0:
    print("FAILED: Ludlow2014 concentration models produced no isolated halos")
    sys.exit(0)
if not (len(massDefault) == len(mass853) == len(mass400) and np.allclose(mass853,massDefault,rtol=1.0e-6,atol=0.0) and np.allclose(mass400,massDefault,rtol=1.0e-6,atol=0.0)):
    print("FAILED: Ludlow2014 concentration models do not contain the same isolated halos")
    sys.exit(0)
toleranceRelative = 1.0e-4
same853           = np.isclose(concentration853,concentrationDefault,rtol=toleranceRelative,atol=0.0)
same400           = np.isclose(concentration400,concentrationDefault,rtol=toleranceRelative,atol=0.0)
if np.all(same853) and not np.all(same400):
    print(f"SUCCESS: Ludlow2014 default C matches C=853 for all {len(concentrationDefault)} halos (and differs from C=400 for {np.count_nonzero(~same400)})")
elif np.all(same400):
    print("FAILED: Ludlow2014 concentrations with the default C are those for C=400 (the Ludlow et al. 2016 default)")
else:
    print(f"FAILED: Ludlow2014 concentrations with the default C differ from those with C=853 for {np.count_nonzero(~same853)} halos")
