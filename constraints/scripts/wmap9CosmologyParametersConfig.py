#!/usr/bin/env python3
# Generates a configuration for cosmological parameters based on the WMAP-9 covariance matrix and suitable for use with the
# Galacticus MCMC constraints architecture. Parameters are expressed as linear combinations of six independent normal deviates
# (labelled "cosmology0" through "cosmology5"). The output of this script can be included into the base parameter file for a
# Galacticus MCMC run, and actual cosmological parameters set by referencing those created in this file.
# Andrew Benson (25-May-2012)

import os
import sys
import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom

# Mapping of WMAP-9 parameter labels to Galacticus parameter names.
parameter_name_mapping = {
    "H_0":     "wmap9HubbleConstant",
    "sigma_8": "wmap9Sigma8",
    "omega_M": "wmap9OmegaMatter",
    "omega_B": "wmap9OmegaBaryon",
    "n_s":     "wmap9PowerSpectrumIndex",
    "tau":     "wmap9ReionizationSuppressionOpticalDepth",
}

# Read the parameters and their covariances from the pre-computed XML file.
data_path = os.environ.get("GALACTICUS_DATA_PATH", "")
xml_file = os.path.join(data_path, "static", "cosmology", "Cosmological_Parameters_WMAP-9.xml")
tree = ET.parse(xml_file)
xml_root = tree.getroot()

# Build ordered list of parameter labels and assign sequential indices.
top_level_params = xml_root.findall("parameter")
param_labels = [p.find("label").text for p in top_level_params]
param_count = len(param_labels)
param_index = {label: i for i, label in enumerate(param_labels)}

# Populate mean vector and covariance matrix.
mean = np.zeros(param_count)
covariance = np.zeros((param_count, param_count))
for param_a in top_level_params:
    label_a = param_a.find("label").text
    idx_a = param_index[label_a]
    mean[idx_a] = float(param_a.find("mean").text)
    for param_b in param_a.findall("parameter"):
        label_b = param_b.find("label").text
        idx_b = param_index[label_b]
        cov_val = float(param_b.find("covariance").text)
        covariance[idx_a, idx_b] = cov_val
        covariance[idx_b, idx_a] = cov_val

# Perform Cholesky decomposition: covariance = L @ L.T
L = np.linalg.cholesky(covariance)

# Build XML output.
root = ET.Element("parameters")
for parameter in param_labels:
    index = param_index[parameter]
    value = f"{mean[index]}"
    for i in range(param_count):
        coeff = L[index, i]
        if coeff != 0.0:
            value += f"+[cosmology{i}]*{coeff}"
    param_name = parameter_name_mapping.get(parameter, parameter)
    # Density parameters are stored as ω=Ωh²; convert to Omega.
    if "omega" in parameter:
        value = f"({value})/([wmap9OmegaMatterHubbleConstant]/100.0)**2"
    # Add an initial "=" so that Galacticus evaluates this parameter.
    value = "=" + value
    ET.SubElement(root, param_name, attrib={"value": value})

print(minidom.parseString(ET.tostring(root)).toprettyxml(indent="  "))
