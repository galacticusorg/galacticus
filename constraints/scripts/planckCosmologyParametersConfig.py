#!/usr/bin/env python3
# Generates a configuration for cosmological parameters based on the Planck covariance matrix and suitable for use with the
# Galacticus MCMC constraints architecture. Parameters are expressed as linear combinations of six independent normal deviates
# (labelled "cosmology0" through "cosmology5"). The output of this script can be included into the base parameter file for a
# Galacticus MCMC run, and actual cosmological parameters set by referencing those created in this file.
# Andrew Benson (16-October-2014)

import os
import re
import sys
import tempfile
import subprocess
import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom

# Specify parameter names of interest and their Galacticus equivalents.
parameter_map = {
    "omegabh2":  "planckCosmologyOmegaBaryon",
    "omegamh2*": "planckCosmologyOmegaMatter",
    "tau":       "planckCosmologyReionizationSuppressionOpticalDepth",
    "ns":        "planckCosmologyPowerSpectrumIndex",
    "H0*":       "planckCosmologyHubbleConstant",
    "sigma8*":   "planckCosmologySigma8",
}

with tempfile.TemporaryDirectory() as tmpdir:
    # Download Planck MCMC chains data set.
    tar_path = os.path.join(tmpdir, "COM_CosmoParams_base-plikHM-TT-lowTEB_R2.00.tar.gz")
    if not os.path.exists(tar_path):
        subprocess.run([
            "wget",
            "http://pla.esac.esa.int/pla-sl/data-action?COSMOLOGY.COSMOLOGY_OID=1200",
            "-O", tar_path,
        ], check=True)

    # Unpack Planck MCMC chains data set.
    planck_dir = os.path.join(tmpdir, "base", "plikHM_TT_lowTEB_lensing")
    if not os.path.exists(planck_dir):
        subprocess.run(["tar", "xvfz", tar_path], cwd=tmpdir, check=True)

    planck_dir = planck_dir + "/"

    # Read Planck parameter names to find the column index of each parameter in the chain files.
    paramnames_file = os.path.join(
        planck_dir,
        "base_plikHM_TT_lowTEB_lensing_post_BAO_H070p6_JLA.paramnames",
    )
    parameters = {}
    with open(paramnames_file) as f:
        for i, line in enumerate(f):
            cols = line.strip().split()
            if cols and cols[0] in parameter_map:
                # Chain files have weight and log-likelihood in columns 0-1; parameters start at column 2.
                parameters[cols[0]] = {"column": i + 2, "chain": []}

    # Check all required parameters were found.
    for key in parameter_map:
        if key not in parameters:
            sys.exit(f"planckCosmologyParametersConfig.py: parameter '{key}' not found")

    # Read Planck chain files, accumulating samples for each parameter.
    sorted_keys = sorted(parameter_map.keys())
    usecols = [parameters[k]["column"] for k in sorted_keys]
    chain_pattern = re.compile(
        r"base_plikHM_TT_lowTEB_lensing_post_BAO_H070p6_JLA_\d+\.txt$"
    )
    for fname in sorted(os.listdir(planck_dir)):
        if chain_pattern.match(fname):
            data = np.loadtxt(os.path.join(planck_dir, fname), usecols=usecols)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            for i, key in enumerate(sorted_keys):
                parameters[key]["chain"].append(data[:, i])

    for key in sorted_keys:
        parameters[key]["chain"] = np.concatenate(parameters[key]["chain"])

    # Compute means and covariance matrix.
    param_count = len(parameter_map)
    mean = np.array([np.mean(parameters[k]["chain"]) for k in sorted_keys])
    covariance = np.zeros((param_count, param_count))
    for i, ki in enumerate(sorted_keys):
        for j, kj in enumerate(sorted_keys):
            covariance[i, j] = np.mean(
                (parameters[ki]["chain"] - mean[i])
                * (parameters[kj]["chain"] - mean[j])
            )

    # Perform Cholesky decomposition: covariance = L @ L.T
    L = np.linalg.cholesky(covariance)

    # Build XML output.
    root = ET.Element("parameters")
    for index, parameter in enumerate(sorted_keys):
        value = f"{mean[index]}"
        for i in range(param_count):
            coeff = L[index, i]
            if coeff != 0.0:
                value += f"+[cosmology{i}]*{coeff}"
        param_name = parameter_map[parameter]
        # Density parameters are stored as ω=Ωh² in Planck chains; convert to Omega.
        if "omega" in parameter:
            value = f"({value})/([planckCosmologyHubbleConstant]/100.0)**2"
        # Apply physical bounds.
        if parameter == "omegamh2*":
            value = f"min({value},1.0)"
        if parameter == "tau":
            value = f"max({value},0.0)"
        # Add an initial "=" so that Galacticus evaluates this parameter.
        value = "=" + value
        ET.SubElement(root, param_name, attrib={"value": value})

print(minidom.parseString(ET.tostring(root)).toprettyxml(indent="  "))
