#!/usr/bin/env python3
# Compute the cosmological parameter covariance matrix from the WMAP-9 year results.
# Andrew Benson (28-December-2012)

import os
import sys
import subprocess
import numpy as np
from datetime import datetime, timezone
import xml.etree.ElementTree as ET
from xml.dom import minidom

# Create a working directory.
execPath = os.environ.get('GALACTICUS_EXEC_PATH', '')
dataPath = os.environ.get('GALACTICUS_DATA_PATH', '')
workDirectory = os.path.join(dataPath, 'dynamic', 'WMAP-9')
os.makedirs(workDirectory, exist_ok=True)

# Download the Monte Carlo Markov Chains.
chainsURL = "http://lambda.gsfc.nasa.gov/data/map/dr5/dcp/chains/wmap_lcdm_wmap9_spt_act_snls3_chains_v5.tar.gz"
chainsFile = os.path.join(workDirectory, 'chains.tar.gz')
if not os.path.exists(chainsFile):
    subprocess.run(['wget', chainsURL, '-O', chainsFile], check=True)

# Unpack the chains.
if not os.path.exists(os.path.join(workDirectory, 'description.txt')):
    subprocess.run(['tar', 'xvfz', 'chains.tar.gz'], cwd=workDirectory, check=True)

# Read the weights.
weightsFile = os.path.join(workDirectory, 'weight_including_2012bao_h0')
weights = np.loadtxt(weightsFile, usecols=(1,))
weights /= weights.sum()

# Define datasets.
datasets = [
    {
        'label'      : 'omega_B',
        'file'       : 'omegabh2',
        'units'      : 'none',
        'description': 'Baryon density parameter, Omega_b*h^2.'
    },
    {
        'label'      : 'omega_M',
        'file'       : 'omegamh2',
        'units'      : 'none',
        'description': 'Matter density parameter, Omega_M*h^2.'
    },
    {
        'label'      : 'tau',
        'file'       : 'tau',
        'units'      : 'none',
        'description': 'Optical depth to reionization.'
    },
    {
        'label'      : 'H_0',
        'file'       : 'H0',
        'units'      : 'km/s/Mpc',
        'description': 'Hubble parameter.'
    },
    {
        'label'      : 'n_s',
        'file'       : 'ns002',
        'units'      : 'none',
        'description': 'Primordial power spectrum spectral index.'
    },
    {
        'label'      : 'sigma_8',
        'file'       : 'sigma8',
        'units'      : 'none',
        'description': 'Root-variance of mass fluctuations in 8Mpc/h radius top hat spheres.'
    }
]

for dataset in datasets:
    filePath = os.path.join(workDirectory, dataset['file'])
    data = np.loadtxt(filePath, usecols=(1,))
    dataset['data'] = data
    dataset['mean'] = np.sum(data * weights)

# Compute covariances.
for i, datasetI in enumerate(datasets):
    datasetI['covariance'] = {}
    for j in range(i + 1):
        datasetJ = datasets[j]
        datasetI['covariance'][datasetJ['label']] = np.sum(
            (datasetI['data'] - datasetI['mean']) *
            (datasetJ['data'] - datasetJ['mean']) *
            weights
        )

# Compute correlations.
for i, datasetI in enumerate(datasets):
    datasetI['correlation'] = {}
    for j in range(i + 1):
        datasetJ = datasets[j]
        datasetI['correlation'][datasetJ['label']] = (
            datasetI['covariance'][datasetJ['label']] /
            np.sqrt(
                datasetI['covariance'][datasetI['label']] *
                datasetJ['covariance'][datasetJ['label']]
            )
        )

# Build XML output.
root = ET.Element('parameters')

descriptions = [
    'WMAP-9 Lambda CDM cosmological parameter best fit values with covariance and correlation matrices.',
    'Uses: wmap9+spt+act+snls3+bao+h0 constraints.'
]
for desc in descriptions:
    ET.SubElement(root, 'description').text = desc

urls = [
    'http://lambda.gsfc.nasa.gov/product/map/dr5/params/lcdm_wmap9_spt_act_snls3.cfm',
    'http://lambda.gsfc.nasa.gov/data/map/dr5/dcp/chains/wmap_lcdm_wmap9_spt_act_snls3_chains_v5.tar.gz'
]
for url in urls:
    ET.SubElement(root, 'url').text = url

for creator in ['Galacticus', 'scripts/aux/parameterCovarianceWMAP9.py']:
    ET.SubElement(root, 'createdBy').text = creator

ET.SubElement(root, 'source').text = 'Computed from Monte Carlo Markov Chains.'

now = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3] + '+00:00'
ET.SubElement(root, 'timestamp').text = now

for i, datasetI in enumerate(datasets):
    paramElem = ET.SubElement(root, 'parameter')
    ET.SubElement(paramElem, 'label'            ).text = datasetI['label']
    ET.SubElement(paramElem, 'description'      ).text = datasetI['description']
    ET.SubElement(paramElem, 'units'            ).text = datasetI['units']
    ET.SubElement(paramElem, 'mean'             ).text = str(datasetI['mean'])
    ET.SubElement(paramElem, 'standardDeviation').text = str(np.sqrt(datasetI['covariance'][datasetI['label']]))
    for j in range(i + 1):
        datasetJ = datasets[j]
        subParam = ET.SubElement(paramElem, 'parameter')
        ET.SubElement(subParam, 'label'      ).text = datasetJ['label']
        ET.SubElement(subParam, 'covariance' ).text = str(datasetI['covariance' ][datasetJ['label']])
        ET.SubElement(subParam, 'correlation').text = str(datasetI['correlation'][datasetJ['label']])

# Pretty-print and write.
xmlStr = minidom.parseString(ET.tostring(root, encoding='unicode')).toprettyxml(indent='  ')
outputFile = os.path.join(dataPath, 'static', 'cosmology', 'Cosmological_Parameters_WMAP-9.xml')
with open(outputFile, 'w') as f:
    f.write(xmlStr)
