#!/usr/bin/env python3
# Generate sets of cosmological parameters drawn at random from the WMAP-9 constraints using the full covariance matrix.
# Andrew Benson (15-September-2010)

import os
import sys
import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom

# Read the parameters and their covariances.
dataPath = os.environ.get('GALACTICUS_DATA_PATH', '')
xmlFile = os.path.join(dataPath, 'static', 'cosmology', 'Cosmological_Parameters_WMAP-9.xml')
tree = ET.parse(xmlFile)
root = tree.getroot()

parameters = root.findall('parameter')
parameterCount = len(parameters)
parameterMap = {p.find('label').text: i for i, p in enumerate(parameters)}

# Create the covariance matrix and means vector.
mean       = np.zeros(parameterCount)
covariance = np.zeros((parameterCount, parameterCount))

for paramA in parameters:
    labelA = paramA.find('label').text
    indexA = parameterMap[labelA]
    mean[indexA] = float(paramA.find('mean').text)
    subParams = paramA.findall('parameter')
    for paramB in subParams:
        labelB = paramB.find('label').text
        indexB = parameterMap[labelB]
        cov = float(paramB.find('covariance').text)
        covariance[indexA, indexB] = cov
        covariance[indexB, indexA] = cov

# Perform a Cholesky decomposition on the covariance matrix.
choleskyDecomposed = np.linalg.cholesky(covariance)

# Generate Gaussian random numbers.
deviates = np.random.normal(size=parameterCount)

# Generate a set of parameters.
params = mean + choleskyDecomposed @ deviates

# Compute required parameters.
h = params[parameterMap['H_0']] / 100.0
Omega_M            = params[parameterMap['omega_M']] / h**2
Omega_DE           = 1.0 - Omega_M
Omega_b            = params[parameterMap['omega_B']] / h**2
sigma_8            = params[parameterMap['sigma_8']]
H_0                = params[parameterMap['H_0']]
powerSpectrumIndex = params[parameterMap['n_s']]

# Build XML output.
root_out = ET.Element('parameters')
cosmo = ET.SubElement(root_out, 'cosmologyParameters')
cosmo.set('value', 'simple')
for name, value in [
    ('Omega_M'           , Omega_M           ),
    ('Omega_DE'          , Omega_DE          ),
    ('Omega_b'           , Omega_b           ),
    ('sigma_8'           , sigma_8           ),
    ('H_0'               , H_0              ),
    ('powerSpectrumIndex', powerSpectrumIndex),
]:
    elem = ET.SubElement(cosmo, name)
    elem.set('value', str(value))

xmlStr = minidom.parseString(ET.tostring(root_out, encoding='unicode')).toprettyxml(indent='  ')
print(xmlStr)
