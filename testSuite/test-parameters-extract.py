#!/usr/bin/env python3
import subprocess
import sys
import h5py
import re
import filecmp
import lxml.etree as ET

# Check that parameter extraction produces consistent results. A simple parameter file is run with Galacticus. A new parameter
# file is then extracted from the HDF5 output. That new parameter file is then run with Galacticus, and another parameter file is
# extracted from this second HDF5 file. The two extracted parameter files are then compared - they should be identical. Finally,
# the first extracted parameter file is compared against the original parameter file - every parameter given in the original must
# be present, with the same value, in the extracted file.
# Andrew Benson (06-September-2024)

# Run the model and check for completion.
print("Running model...")
status = subprocess.run("mkdir -p outputs",shell=True)
log = open("outputs/test-parameters-extract.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/parametersExtract.xml; ./scripts/parameters/parametersExtract.py testSuite/outputs/parametersExtract.hdf5 testSuite/outputs/parametersExtract.xml; ./Galacticus.exe testSuite/outputs/parametersExtract.xml; ./scripts/parameters/parametersExtract.py testSuite/outputs/parametersExtract.hdf5 testSuite/outputs/parametersExtractSecond.xml; ./Galacticus.exe testSuite/parameters/parametersExtract.xml --dry-run;  ./scripts/parameters/parametersExtract.py testSuite/outputs/parametersExtract.hdf5 testSuite/outputs/parametersExtractDryRun.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0:
    print("...done ("+str(status)+")")
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-parameters-extract.log",shell=True)
    sys.exit(0)
else:
    print("...done")
    print("Checking for errors...")
    status = subprocess.run("grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" outputs/test-parameters-extract.log",shell=True)
    if status.returncode == 0:
        print("...done ("+str(status)+")")
        print("FAILED: model run (errors):")
        subprocess.run("cat outputs/test-parameters-extract.log",shell=True)
        sys.exit(0)
    else:
        print("...done")
        print("SUCCESS: model run")


# Compare the two extracted parameter files - they should be identical.
status = "SUCCESS" if filecmp.cmp("outputs/parametersExtract.xml","outputs/parametersExtractSecond.xml") else "FAILED"
print(status+": consistent extracted parameters")

# Parse the output parameter file.
parameters = ET.parse("outputs/parametersExtractDryRun.xml")

# Check that default parameters have correct values.
linearGrowth  = parameters.findall("./linearGrowth[@value='collisionlessMatter']")
if len(linearGrowth) == 1:
    print('SUCCESS: default `<linearGrowth value="collisionlessMatter"/>` is present')
else:
    print('FAILED: default `<linearGrowth value="collisionlessMatter"/>` is not present')

# Check that references and targets are correctly included.
cosmologyTarget  = parameters.findall("./cosmologyParameters[@id='refCosmo']")
cosmologyPointer = parameters.findall("./transferFunction/cosmologyParameters[@idRef='refCosmo']")
if len(cosmologyTarget ) == 1:
    print('SUCCESS: target `<cosmologyParameters/>` is present')
else:
    print('FAILED: target `<cosmologyParameters/>` is not present')
if len(cosmologyPointer) == 1:
    print('SUCCESS: reference `<cosmologyParameters/>` is present')
else:
    print('FAILED: reference `<cosmologyParameters/>` is not present')

# Check that the extracted parameter file is equivalent to the parameter file that was run. Comparing two *extracted* files (as
# above) can not detect a parameter which never reaches the output file's `Parameters` group at all, since both extractions are
# then equally lossy and so still agree. Comparing against the original catches that: a parameter which the user supplied but
# which is absent from the extracted file means the extracted file would run different physics.
def normalizedValue(value):
    # Normalize a parameter value for comparison. Values are written to the output file verbatim, but may be reformatted (e.g.
    # by dropping insignificant zeros) when read back and written out again, so numerical values are compared numerically.
    value = " ".join(value.split())
    try:
        return " ".join(map(lambda x: "%.10g" % float(x),value.split(" ")))
    except ValueError:
        return value

def parameterCounts(fileName):
    # Return a count, for each (path, value) pair, of the parameters in the given parameter file. Paths are used (rather than
    # elements) so that a parameter is matched only where it appears in the same place in the parameter hierarchy, and counts are
    # used so that repeated instances of the same parameter must be repeated in both files.
    counts = {}
    def accumulate(element,path):
        for child in element:
            if not isinstance(child.tag,str):
                continue
            childPath = path+"/"+child.tag
            if   'idRef' in child.attrib:
                value = "{idRef:"+child.attrib['idRef']+"}"
            elif 'value' in child.attrib:
                value = normalizedValue(child.attrib['value'])
            else:
                value = None
            counts[(childPath,value)] = counts.get((childPath,value),0)+1
            accumulate(child,childPath)
    accumulate(ET.parse(fileName).getroot(),"")
    return counts

# Parameters which are expected to be absent from the extracted file. The parameter file deliberately includes a duplicated
# `mergerTreeOutputter` (to ensure that such duplicates do not cause errors in parameter extraction); the duplicate is never
# used, and so is never recorded in the output file.
ignoredParameters = (
    ("/mergerTreeOutputter","null"),
)
countsOriginal  = parameterCounts("parameters/parametersExtract.xml")
countsExtracted = parameterCounts("outputs/parametersExtract.xml"   )
missing         = []
for (path, value), count in countsOriginal.items():
    if (path,value) in ignoredParameters:
        continue
    if countsExtracted.get((path,value),0) < count:
        missing.append(path+(" = "+value if value is not None else ""))
if len(missing) == 0:
    print("SUCCESS: extracted parameters match the original parameter file")
else:
    print("FAILED: parameters present in the original parameter file are absent from the extracted parameter file:")
    for parameter in sorted(missing):
        print("   "+parameter)
