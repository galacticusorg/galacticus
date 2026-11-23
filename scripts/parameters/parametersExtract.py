#!/usr/bin/env python3
import argparse
import sys
import lxml.etree as ET
from Galacticus.Parameters.extract import parameters_from_hdf5

# Extract parameters from a Galacticus HDF5 file and output as an Galacticus XML parameter file.
# Andrew Benson (06-September-2024)

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='parametersExtract.py',description='Extract parameters from a Galacticus HDF5 file and output as an Galacticus XML parameter file.')
parser.add_argument('hdf5FileName')
parser.add_argument('xmlFileName' )
parser.add_argument('--dereference', action='store_true',help='dereference any `idRef` pointers')
parser.add_argument('--omit-defaulted', action='store_true',help='omit parameters which took their default values')
args = parser.parse_args()

# Extract the parameters.
try:
    parameters = parameters_from_hdf5(args.hdf5FileName,dereference=args.dereference,include_defaulted=not args.omit_defaulted)
except ValueError as error:
    print(error)
    sys.exit(1)

# Write out the parameter file.
tree = ET.ElementTree(parameters)
ET.indent(tree, space="  ", level=0)
fileOut = open(args.xmlFileName,"w")
fileOut.write(str(ET.tostring(parameters, encoding="unicode")))
fileOut.close()
