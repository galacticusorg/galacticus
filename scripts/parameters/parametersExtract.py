#!/usr/bin/env python3
import h5py
import numpy as np
import xml.etree.ElementTree as ET
import sys
import re
import argparse

# Extract parameters from a Galacticus HDF5 file and output as an Galacticus XML parameter file.
# Andrew Benson (06-September-2024)

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='parametersExtract.py',description='Extract parameters from a Galacticus HDF5 file and output as an Galacticus XML parameter file.')
parser.add_argument('hdf5FileName')
parser.add_argument('xmlFileName' )
args = parser.parse_args()

# Create the root `parameters` element.
parameters = ET.Element("parameters")

# Define a function to handle creation of the parameter file structure.
def createStructure(name, node):
    # Ignore anything other than groups.
    if isinstance(node, h5py.Group):
        # Split the parameter name into individual elements.
        elementNames = name.split("/")
        # Find the parent element for this new parameter.
        if len(elementNames) > 1:
            parent = parameters.find("/".join(elementNames[0:-1]))
        else:
            parent = parameters
        # Array parameters (i.e. multiple copies of the same named parameter) are defined by a `[N]` suffix. Detect these.
        isArray = re.match(r'(.*)\[(\d+)\]',elementNames[-1])
        child = None
        if isArray:
            # For array parameters, construct a list of elements of length equal to the current array index.
            elementName  =     isArray.group(1)
            elementIndex = int(isArray.group(2))
            members      = parent.findall(elementName)
            if len(members) <= elementIndex:
                for i in range(len(members),elementIndex):
                    child = ET.SubElement(parent, elementName)
        else:
            # For non-array parameters simply add an appropriately-named child element.
            child = ET.SubElement(parent, elementNames[-1])
            
# Create a function to assign parameter values.
def assignValues(name, node):
    # Ignore everything other than groups.
    if isinstance(node, h5py.Group):
        # Find the parent element for this group.
        if name == "parameters":
            parent = parameters
        else:
            parent = parameters.find(name)
        # Iterate over all attributes in this group.
        for attributeName in node.attrs.keys():
            isReference = None
            idReference = None
            # Extract the attribute value and convert to a string.
            value = node.attrs[attributeName]
            if isinstance(value,bytes):
                value = value.decode()
            elif isinstance(value,list):
                value = " ".join(map(lambda x: x.decode(),value))
            elif isinstance(value,np.ndarray):
                if isinstance(value[0],bytes):
                    value = " ".join(map(lambda x: x.decode(),value))
                else:
                    value = " ".join(map(lambda x: str(x),value))
            value = str(value)
            isReference = re.match(r'\{idRef:([a-zA-Z0-9_]+)\}',value)
            if isReference:
                idReference = isReference.group(1)
            isTarget = re.match(r'(.*)\{id:([a-zA-Z0-9_]+)\}',attributeName)
            if isTarget:
                attributeName = isTarget.group(1)
                idTarget      = isTarget.group(2)
            # Set the value of the parameter (creating the element first if needed).
            if attributeName in node:
                child = parent.find(attributeName)
                if isReference:
                    child.set('idRef', idReference)
                else:
                    child.set('value', value      )
                    if isTarget:
                        child.set('id',idTarget)
            else:
                newChild = ET.SubElement(parent, attributeName)
                if isReference:
                    newChild.set('idRef', idReference)
                else:
                    newChild.set('value', value      )
                    if isTarget:
                        newChild.set('id',idTarget)

# Open the HDF5 file and get the `Parameters` group.
fileIn = h5py.File(args.hdf5FileName,"r")
parametersGroup = fileIn['Parameters']

# Add a `formatVersion` element.
formatVersion      = ET.SubElement(parameters, 'formatVersion')
formatVersion.text = "2"

# Add a `lastModified` element.
version            = fileIn['Version']
if 'gitHash' in version.attrs.keys():
    lastModified      = ET.SubElement(parameters, 'lastModified')
    lastModified.set('revision',version.attrs['gitHash'].decode())

# Construct the structure of the parameter file by visiting all groups in `Parameters`.
parametersGroup.visititems(createStructure)

# Assign top-level parameter values.
assignValues('parameters',parametersGroup)

# Assign parameter values in all sub-parameters.
parametersGroup.visititems(assignValues)

# Write out the parameter file.
tree = ET.ElementTree(parameters)
ET.indent(tree, space="  ", level=0)
fileOut = open(args.xmlFileName,"w")
fileOut.write(str(ET.tostring(parameters, encoding="unicode")))
fileOut.close()
