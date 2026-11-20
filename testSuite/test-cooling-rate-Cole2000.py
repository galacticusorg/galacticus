#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check that the Cole et al. (2000) cooling rate continues to evolve between halo formation events. In that model the density
# profile of the hot gas is frozen at the time of halo formation, but the time available for cooling (the time since formation)
# continues to grow, so the cooling radius grows and the cooling rate declines. The main branch of a single tree is followed through
# a sequence of outputs, and between consecutive outputs with the same formation time (i.e. no formation event between them) the
# cooling rate must change (and fall).
# Andrew Benson, Claude (19-November-2026)

subprocess.run("mkdir -p outputs",shell=True)
log    = open("outputs/test-cooling-rate-Cole2000.log","w")
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/coolingRateCole2000.xml",stdout=log,stderr=log,shell=True)
log.close()
if status.returncode != 0 or subprocess.run("grep -q -i -e fatal -e aborted outputs/test-cooling-rate-Cole2000.log",shell=True).returncode == 0:
    print("FAILED: model run:")
    subprocess.run("cat outputs/test-cooling-rate-Cole2000.log",shell=True)
    sys.exit(0)
print("SUCCESS: model run")

# Follow the most massive node of the (single) tree through the outputs.
rows = []
with h5py.File("outputs/coolingRateCole2000.hdf5","r") as model:
    outputs = model["Outputs"]
    for name in sorted(outputs.keys(),key=lambda name: int(name.replace("Output",""))):
        nodes = outputs[name]["nodeData"]
        i     = np.argmax(nodes["basicMass"][:])
        rows.append((outputs[name].attrs["outputTime"],nodes["nodeFormationTime"][i],nodes["hotHaloRateCooling"][i]))
pairs     = 0
unchanged = 0
rising    = 0
for (time1,timeFormation1,rate1), (time2,timeFormation2,rate2) in zip(rows[:-1],rows[1:]):
    # Consider only consecutive outputs with no intervening formation event, and with cooling occurring.
    if timeFormation1 != timeFormation2 or rate1 <= 0.0 or rate2 <= 0.0:
        continue
    pairs += 1
    if rate2 == rate1:
        unchanged += 1
    elif rate2 > rate1:
        rising    += 1
if pairs == 0:
    print("FAILED: no pairs of outputs without an intervening formation event")
elif unchanged > 0:
    print(f"FAILED: cooling rate unchanged between {unchanged} of {pairs} pairs of outputs with no intervening formation event")
elif rising > 0:
    print(f"FAILED: cooling rate rises between {rising} of {pairs} pairs of outputs with no intervening formation event")
else:
    print(f"SUCCESS: cooling rate falls between all {pairs} pairs of outputs with no intervening formation event")
