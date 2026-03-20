#!/usr/bin/env python3
import subprocess
import sys
import re
import h5py
import numpy as np

# Run models that test reproducibility.
# Andrew Benson (ported to Python)

# Make output directory.
subprocess.run("mkdir -p outputs/reproducibility", shell=True)

# Define constants.
gravitationalConstant = 4.3011827419096073e-9

# Define tests.
tests = [
    {
        "name":           "cooling",
        "parameters":     "testSuite/parameters/reproducibility/cooling.xml",
        "outputFileName": "testSuite/outputs/reproducibility/cooling.hdf5",
        "assertions": [
            {
                "name":              "hot halo mass",
                "output":            1,
                "property":          "hotHaloMass",
                "values":            np.array([7.99799857410098e10]),
                "toleranceRelative": 4.0e-6,
            }
        ],
    },
    {
        "name":           "closedBox",
        "parameters":     "testSuite/parameters/reproducibility/closedBox.xml",
        "outputFileName": "testSuite/outputs/reproducibility/closedBox.hdf5",
        "assertions": [
            {"name": "gas mass",        "output": 1, "property": "diskMassGas",                 "values": np.array([9.0717953e9]),   "toleranceRelative": 1.0e-2},
            {"name": "stellar mass",    "output": 1, "property": "diskMassStellar",              "values": np.array([9.0928205e10]), "toleranceRelative": 1.0e-2},
            {"name": "gas metals",      "output": 1, "property": "diskAbundancesGasMetals",      "values": np.array([9.0717953e8]),  "toleranceRelative": 1.0e-2},
            {"name": "stellar metals",  "output": 1, "property": "diskAbundancesStellarMetals",  "values": np.array([2.8814957e9]),  "toleranceRelative": 1.0e-2},
        ],
    },
    {
        "name":           "leakyBox",
        "parameters":     "testSuite/parameters/reproducibility/leakyBox.xml",
        "outputFileName": "testSuite/outputs/reproducibility/leakyBox.hdf5",
        "assertions": [
            {"name": "gas mass",       "output": 1, "property": "diskMassGas",                 "values": np.array([4.0762204e9]),   "toleranceRelative": 1.1e-2},
            {"name": "stellar mass",   "output": 1, "property": "diskMassStellar",              "values": np.array([3.5971417e10]), "toleranceRelative": 1.0e-2},
            {"name": "gas metals",     "output": 1, "property": "diskAbundancesGasMetals",      "values": np.array([2.03811e8]),    "toleranceRelative": 1.0e-2},
            {"name": "stellar metals", "output": 1, "property": "diskAbundancesStellarMetals",  "values": np.array([4.85624e8]),    "toleranceRelative": 1.0e-2},
        ],
    },
    {
        "name":           "adiabaticContraction",
        "parameters":     "testSuite/parameters/reproducibility/adiabaticContraction.xml",
        "outputFileName": "testSuite/outputs/reproducibility/adiabaticContraction.hdf5",
        "assertions": [
            {
                "name":              "spheroid radius",
                "output":            1,
                "property":          "spheroidRadius",
                "values":            np.array([0.00360702918954165]),
                "toleranceRelative": 2.0e-4,
            },
            {
                "name":       "spheroid angular momentum",
                "output":     1,
                "expression": "(%[spheroidRadius]*%[spheroidVelocity]*%[spheroidMassStellar])/%[spheroidAngularMomentum]",
                "values":     np.array([0.5]),
                "toleranceRelative": 2.0e-4,
            },
            {
                "name":       "rotation curve",
                "output":     1,
                "expression": "%[rotationCurve{0}]/%[spheroidVelocity]",
                "values":     np.array([1.0]),
                "toleranceRelative": 2.0e-4,
            },
            {
                "name":       "initial specific angular momentum",
                "output":     1,
                "expression": (
                    "+np.sqrt(0.84333333)"
                    "*%[darkMatterOnlyVelocityVirial]"
                    "/%[spheroidVelocity]"
                    "*("
                    "  +%[darkMatterOnlyRadiusVirial]"
                    "  *("
                    "    +%[rotationCurve{1}]**2"
                    "    *%[spheroidRadius]"
                    "    /gravitationalConstant"
                    "    /0.83333333"
                    "   )"
                    "  /%[basicMass]"
                    " )"
                    "/%[spheroidRadius]"
                ),
                "values":     np.array([1.0]),
                "toleranceRelative": 3.0e-3,
            },
        ],
    },
]


def evaluate_expression(expression, nodeData, gravitationalConstant):
    """Evaluate an expression with %[propertyName] and %[propertyName{index}] substitutions."""
    # Find all property references.
    propertyNames = set()
    for m in re.finditer(r'%\[([a-zA-Z0-9_]+)(?:\{(\d+)\})?\]', expression):
        propertyNames.add(m.group(1))

    # Load all needed properties.
    properties = {}
    for name in propertyNames:
        properties[name] = nodeData[name][:]

    # Build the substituted expression.
    expr = expression
    # Replace indexed references first (more specific pattern).
    expr = re.sub(
        r'%\[([a-zA-Z0-9_]+)\{(\d+)\}\]',
        lambda m: f"properties['{m.group(1)}'][:,{m.group(2)}]",
        expr
    )
    # Replace plain references.
    expr = re.sub(
        r'%\[([a-zA-Z0-9_]+)\]',
        lambda m: f"properties['{m.group(1)}']",
        expr
    )
    return eval(expr)


# Run tests.
for test in tests:
    print(f"Running test \"{test['name']}\": Galacticus.exe {test['parameters']}")
    status = subprocess.run(f"cd ..; ./Galacticus.exe {test['parameters']}", shell=True)
    if status.returncode == 0:
        with h5py.File(f"../{test['outputFileName']}", "r") as model:
            outputs = model["Outputs"]
            for assertion in test["assertions"]:
                output   = outputs[f"Output{assertion['output']}"]
                nodeData = output["nodeData"]
                if "property" in assertion:
                    values = nodeData[assertion["property"]][:]
                elif "expression" in assertion:
                    values = evaluate_expression(assertion["expression"], nodeData, gravitationalConstant)
                else:
                    raise ValueError("assertion gives no property or expression to test")

                difference   = np.abs(values - assertion["values"])
                allowedError = assertion["toleranceRelative"] * np.abs(values)
                if np.all(difference < allowedError):
                    print(f"SUCCESS: assertion '{assertion['name']}' of reproducibility test '{test['name']}' passed")
                else:
                    print(f"FAIL: assertion '{assertion['name']}' of reproducibility test '{test['name']}' failed")
                    for i in range(len(difference)):
                        print(f"\t{values[i]}\t{assertion['values'][i]}\t{difference[i]}\t{allowedError[i]}")
    else:
        print(f"FAIL: reproducibility test '{test['name']}' model failed to run")
