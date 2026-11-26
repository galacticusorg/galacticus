#!/usr/bin/env python3
import subprocess
import sys
import xml.etree.ElementTree as ET
import h5py
import numpy as np

# Check that the Cole et al. (2000) cooling rate continues to evolve between halo formation events, with the time available for
# cooling measured from the most recent formation event. In that model the density profile of the hot gas is frozen at the time of
# halo formation, but the time available for cooling (the time since formation) continues to grow. The main branch of a single tree
# is followed through a sequence of outputs.
#
# 1. In the standard model, the cooling rate must change between consecutive outputs with no formation event between them.
#
# 2. With the cooling radius of an isothermal profile, r_c = r_v (t_a/t_c)^{1/2}, where t_a = t - t_f is the time available and
#    t_c the cooling time at the virial radius r_v (both of the frozen profile, and so fixed between formation events, for a
#    cooling function which does not depend on epoch). With the beta-profile (beta=2/3) density of the hot gas, rho = rho_0/(1+r^2/r_0^2),
#    the cooling rate, 4 pi r_c^2 rho(r_c) dr_c/dt, is then K t_a^{1/2}/(1+B t_a) for constants K and B. So t_a^{1/2}/rate must be
#    exactly linear in t_a between formation events - which holds only if t_a is measured from the most recent formation event.
# Andrew Benson, Claude (19-November-2026)

def mainBranch(fileName):
    # Follow the most massive node of the (single) tree through the outputs.
    rows = []
    with h5py.File(fileName,"r") as model:
        outputs = model["Outputs"]
        for name in sorted(outputs.keys(),key=lambda name: int(name.replace("Output",""))):
            nodes = outputs[name]["nodeData"]
            i     = np.argmax(nodes["basicMass"][:])
            rows.append((outputs[name].attrs["outputTime"],nodes["nodeFormationTime"][i],nodes["hotHaloRateCooling"][i]))
    return np.array(rows)

def runModel(parameterFile,label):
    log    = open(f"outputs/test-cooling-rate-Cole2000-{label}.log","w")
    status = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFile}",stdout=log,stderr=log,shell=True)
    log.close()
    if status.returncode != 0 or subprocess.run(f"grep -q -i -e fatal -e aborted outputs/test-cooling-rate-Cole2000-{label}.log",shell=True).returncode == 0:
        print(f"FAILED: model run [{label}]:")
        subprocess.run(f"cat outputs/test-cooling-rate-Cole2000-{label}.log",shell=True)
        sys.exit(0)
    print(f"SUCCESS: model run [{label}]")

subprocess.run("mkdir -p outputs",shell=True)

# Standard model.
runModel("testSuite/parameters/coolingRateCole2000.xml","standard")
rows      = mainBranch("outputs/coolingRateCole2000.hdf5")
pairs     = 0
unchanged = 0
for (time1,timeFormation1,rate1), (time2,timeFormation2,rate2) in zip(rows[:-1],rows[1:]):
    # Consider only consecutive outputs with no intervening formation event, and with cooling occurring.
    if timeFormation1 != timeFormation2 or rate1 <= 0.0 or rate2 <= 0.0:
        continue
    pairs += 1
    if rate2 == rate1:
        unchanged += 1
if pairs == 0:
    print("FAILED: no pairs of outputs without an intervening formation event")
elif unchanged > 0:
    print(f"FAILED: cooling rate unchanged between {unchanged} of {pairs} pairs of outputs with no intervening formation event")
else:
    print(f"SUCCESS: cooling rate changes between all {pairs} pairs of outputs with no intervening formation event")

# Model with an isothermal cooling radius, and a cooling function which does not depend on epoch.
tree       = ET.parse("parameters/coolingRateCole2000.xml")
parameters = tree.getroot()
for name, value in (("coolingInfallRadius","coolingRadius"),("coolingRadius","isothermal"),("coolingFunction","atomicCIECloudy")):
    element = parameters.find(name)
    if element is None:
        element = ET.SubElement(parameters,name)
    for child in list(element):
        element.remove(child)
    element.set("value",value)
parameters.find("outputFileName").set("value","testSuite/outputs/coolingRateCole2000Isothermal.hdf5")
tree.write("outputs/coolingRateCole2000Isothermal.xml")
runModel("testSuite/outputs/coolingRateCole2000Isothermal.xml","isothermal")
rows      = mainBranch("outputs/coolingRateCole2000Isothermal.hdf5")
intervals = 0
failed    = False
for timeFormation in np.unique(rows[:,1]):
    select = (rows[:,1] == timeFormation) & (rows[:,2] > 0.0)
    if np.count_nonzero(select) < 3:
        continue
    intervals    += 1
    timeAvailable = rows[select,0]-timeFormation
    y             = np.sqrt(timeAvailable)/rows[select,2]
    residual      = np.max(np.abs(np.polyval(np.polyfit(timeAvailable,y,1),timeAvailable)/y-1.0))
    if residual > 1.0e-6:
        print(f"FAILED: for formation time {timeFormation:.4f} Gyr, (time available)^1/2/(cooling rate) is not linear in the time available since formation (maximum fractional residual {residual:.3e})")
        failed = True
if intervals == 0:
    print("FAILED: no intervals between formation events with at least three outputs with cooling")
elif not failed:
    print(f"SUCCESS: (time available)^1/2/(cooling rate) is linear in the time available since formation in all {intervals} intervals between formation events")
