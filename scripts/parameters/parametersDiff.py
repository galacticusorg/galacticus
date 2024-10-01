#!/usr/bin/env python3
import os
import xml.etree.ElementTree as ET
import argparse
from pathlib import Path
import urllib.request
import tarfile
import subprocess
import sys
import tempfile
import re

# Show differences between two Galacticus parameter files. 
# Andrew Benson (07-September-2024)

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='parametersDiff.py',description='Show differences between two Galacticus parameter files.')
parser.add_argument('parameterFile1')
parser.add_argument('parameterFile2')
parser.add_argument('--respectOrder', action='store_true',help='respect the order of elements when comparing files')
parser.add_argument('--canonicalizeValues', action='store',help='canonicalize all numerical values before comparing')
args = parser.parse_args()

# Install xdiff if necessary.
dynamicPath = os.environ['GALACTICUS_DATA_PATH']+"/dynamic"
xdiffPath   = dynamicPath+"/xdiff-2.4"
xdiff       = Path(xdiffPath+"/xdiff.py")
if not xdiff.is_file():
    tarFile = xdiffPath+".tar.gz"
    urllib.request.urlretrieve("https://hg.sr.ht/~nolda/xdiff/archive/2.4.tar.gz", tarFile)
    tarball = tarfile.open(tarFile)
    tarball.extractall(dynamicPath)
    tarball.close() 

# Create list of filwnames to compare.
fileNames    = [ args.parameterFile1, args.parameterFile2 ]
fileNamesTmp = [ ]

# If parameter order is not to be respected, created copies of our files with parameters sorted by name.
if not args.respectOrder:
    for i in range(2):
        parametersDoc = ET.parse(fileNames[i])
        parameters    = parametersDoc.getroot()
        for parent in parameters.iter():
            parent[:] = sorted(parent,key=lambda x: x.tag)
        ET.indent(parametersDoc, space="  ", level=0)
        fileOut = tempfile.NamedTemporaryFile(mode="w",encoding="utf8",delete=False)
        fileOut.write(ET.tostring(parameters, encoding="unicode"))
        fileOut.close()
        # Replace the file name with that of our temporary file.
        fileNames[i] = fileOut.name
        fileNamesTmp.append(fileOut.name)

# If values are to be canonicalized, do so.
if args.canonicalizeValues:
    formatCanonical = "{:"+args.canonicalizeValues+"}"
    for i in range(2):
        parametersDoc = ET.parse(fileNames[i])
        parameters    = parametersDoc.getroot()
        for parent in parameters.iter():
            if "value" not in parent.attrib:
                continue
            values          = re.split(r'\s+',parent.attrib['value'].strip())
            valuesCanonical = []
            for value in values:
                try:
                    valueNumeric = float(value)
                    valuesCanonical.append(formatCanonical.format(valueNumeric))
                except ValueError:
                    valuesCanonical.append(value)
            parent.attrib['value'] = " ".join(valuesCanonical)
                
        ET.indent(parametersDoc, space="  ", level=0)
        fileOut = tempfile.NamedTemporaryFile(mode="w",encoding="utf8",delete=False)
        fileOut.write(ET.tostring(parameters, encoding="unicode"))
        fileOut.close()
        # Replace the file name with that of our temporary file.
        fileNames[i] = fileOut.name
        fileNamesTmp.append(fileOut.name)    
            
# Run `xdiff` to compare the files.
status = subprocess.run("python3 "+str(xdiff)+" "+" ".join(fileNames),shell=True)

# Remove any temporary files.
if not args.respectOrder:
    for i in range(len(fileNamesTmp)):
        os.unlink(fileNamesTmp[i])

# Return diff status.
sys.exit(status.returncode)

