#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Run a test case which checks that volume weights of merger trees read from file are assigned
# correctly when they are computed from the simulation box size.
# Andrew Benson (16-January-2014; ported to Python)

# Run the model and check for successful completion.
status = subprocess.run("./Galacticus.exe testSuite/regressions/mergerAtFinalTimeInTree.xml", shell=True)
if status.returncode != 0:
    print("FAILED: mergerTreeBoxSizeWeight model failed to complete")
    sys.exit(1)

# Extract required information from files.
with h5py.File("testSuite/data/mergerTrees/mergerAtFinalTimeInTree.hdf5", "r") as trees:
    box_size                   = trees["simulation" ].attrs["boxSize"               ]
    length_hubble_exponent      = trees["units"      ].attrs["lengthHubbleExponent"  ]
    length_scale_factor_exponent = trees["units"      ].attrs["lengthScaleFactorExponent"]
    length_units_in_SI          = trees["units"      ].attrs["lengthUnitsInSI"       ]
    hubble                     = trees["cosmology"   ].attrs["HubbleParam"           ]

with h5py.File("testSuite/outputs/regressions/mergerAtFinalTimeInTree.hdf5", "r") as model:
    tree_weight = model["Outputs/Output1/mergerTreeWeight"][:]

# Compute expected tree weight.
scale_factor = 1.0
mega_parsec  = 3.08567758135e+22
box_size     = (box_size
                * hubble**length_hubble_exponent
                * scale_factor**length_scale_factor_exponent
                * length_units_in_SI / mega_parsec)
weight_expected  = 1.0 / box_size**3

# Compute fractional difference.
fractional_error = np.abs(tree_weight - weight_expected) / weight_expected

# Check for consistency.
if np.any(fractional_error >= 1.0e-6):
    print("FAILED: mergerTreeBoxSizeWeight tree weights do not equal expected values")
    sys.exit(1)

print("SUCCESS: mergerTreeBoxSizeWeight")
