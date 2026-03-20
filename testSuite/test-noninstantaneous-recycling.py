#!/usr/bin/env python3
import subprocess
import sys
import re
import h5py
import numpy as np

# Check calculations of noninstantaneous recycling.
# Andrew Benson (ported to Python)

# Run the model.
subprocess.run("mkdir -p outputs", shell=True)
with open("outputs/noninstantaneous_recycling.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/noninstantaneous_recycling.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED:  model run:")
    with open("outputs/noninstantaneous_recycling.log") as f:
        print(f.read())
else:
    print("SUCCESS: model run")

yields = {}
with open("outputs/noninstantaneous_recycling.log") as logFile:
    for line in logFile:
        m = re.search(r'(\S+/yield(Metals|Fe)_[a-z0-9]+\.hdf5)', line)
        if m:
            fileName = m.group(1)
            ytype    = m.group(2)
            yields.setdefault('file', {})[ytype] = fileName
            with h5py.File(fileName, "r") as f:
                yields[ytype] = f["yield" + ytype][:]

if "Metals" not in yields or "Fe" not in yields:
    for ytype in ("Metals", "Fe"):
        if ytype in yields.get("file", {}):
            print(f"found {ytype} yield in '{yields['file'][ytype]}'")
        else:
            print(f"failed to find file for {ytype} yield")
            print("log file is:")
            with open("outputs/noninstantaneous_recycling.log") as f:
                print(f.read())
    print("FAILED: could not read yield tables")
    sys.exit(0)

haveMetals   = np.where(yields["Metals"].flat[:] > 0.0)[0]
ratioMaximum = np.max(yields["Fe"].flat[haveMetals] / yields["Metals"].flat[haveMetals])

# Read the model data and check for consistency.
with h5py.File("outputs/noninstantaneous_recycling.hdf5", "r") as model:
    nodeData  = model["Outputs/Output10/nodeData"]
    dataNames = [
        "spheroidAbundancesStellarMetals",
        "spheroidAbundancesStellarFe",
        "diskAbundancesStellarMetals",
        "diskAbundancesStellarFe",
    ]
    data = {name: nodeData[name][:] for name in dataNames}

massMetals = data["spheroidAbundancesStellarMetals"] + data["diskAbundancesStellarMetals"]
massFe     = data["spheroidAbundancesStellarFe"]     + data["diskAbundancesStellarFe"]
nonZero    = np.where(massMetals > 1.0)[0]
ratio      = massFe[nonZero] / massMetals[nonZero]
status_str = "SUCCESS" if np.all(ratio < ratioMaximum) else "FAILED"
print(f"{status_str}: Fe/Z ratio")
