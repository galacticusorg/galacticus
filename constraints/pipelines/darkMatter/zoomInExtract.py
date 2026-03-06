#!/usr/bin/env python3
import numpy as np
import h5py
import sys
import os
import re
import time
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
    
# Locate and extract properties of the primary halo in the zoom-in simulation.
# Andrew Benson (15-March-2022)

# Get arguments.
if len(sys.argv) != 9:
    sys.exit("Usage: zoomInExtract.py <pathName> <primaryHaloFileName> <expansionFactor> <hubbleConstant> <massParticle> <massHostLogMin> <massHostLogMax> <hostHaloID>")
pathName            =       sys.argv[1]
primaryHaloFileName =       sys.argv[2]
expansionFactor     = float(sys.argv[3])
hubbleConstant      = float(sys.argv[4])
massParticle        = float(sys.argv[5])
massHostLogMin      = float(sys.argv[6])
massHostLogMax      = float(sys.argv[7])
hostHaloID          = int  (sys.argv[8])

# Extract box and region size from the Music config file.
print("\tExtracting high-resolution region")
musicFileName = pathName+"music.conf"
musicFile = open(musicFileName,'r')
for line in musicFile:
    match = re.match(r'^boxlength\s*=\s*([\d\.]+)',line)
    if match:
        boxSizeOriginal = float(match.group(1))
    match = re.match(r'^ref_extent\s*=\s*([\d\.\,\s]+)',line)
    if match:
        regionExtent = np.array(list(map(lambda x: float(x),re.split(r'\s*,\s*',match.group(1)))))
musicFile.close()
boxSizeOriginal /= hubbleConstant
boxSize          = boxSizeOriginal*np.product(regionExtent)**(1.0/3.0)
# Sleep while joint processing is underway.
while os.path.exists(pathName+"tree.lock"):
    time.sleep(5)
file_path = Path(pathName+"tree.lock")
file_path.touch()
# Extract the name of the tree file being processed.
treeFileName = pathName+"tree_?_?_?.dat"
# Read the merger tree file looking for the position and virial radius of the most massive halo.
print("\tPre-processing tree file")
if not os.path.exists(pathName+"treeReduced.dat"):
    status = subprocess.run("awk '{if (substr($1,1,1) != \"#\" && NF > 1 && $15 == 1) print $1,$2,$4,$11,$12,$18,$19,$20}' "+treeFileName+" > "+pathName+"treeReduced.dat",shell=True)
if not os.path.exists(pathName+"treeCount.dat"):
    status = subprocess.run("wc -l "+pathName+"treeReduced.dat > "+pathName+"treeCount.dat",shell=True)
treeMetaData = open(pathName+"treeCount.dat",'r')
treeCount = int(re.split(r'\s+',treeMetaData.readline())[0])
os.remove(pathName+"tree.lock")
print(f'\t\tPre-processed tree contains {treeCount} halos')
print('\tReading tree file')
# Sleep while joint processing is underway.
while os.path.exists(pathName+"treeRaw.lock"):
    time.sleep(5)
file_path = Path(pathName+"treeRaw.lock")
file_path.touch()
if os.path.exists(pathName+"treeReduced.hdf5"):
    treeRaw = h5py.File(pathName+"treeReduced.hdf5",'r')
    haloExpansionFactor = treeRaw['haloExpansionFactor'][:]
    haloID              = treeRaw['haloID'             ][:]
    descID              = treeRaw['descID'             ][:]
    massHalo            = treeRaw['massHalo'           ][:]
    radiusHalo          = treeRaw['radiusHalo'         ][:]
    x                   = treeRaw['x'                  ][:]
    y                   = treeRaw['y'                  ][:]
    z                   = treeRaw['z'                  ][:]
else:
    haloExpansionFactor = np.zeros(treeCount)
    haloID              = np.zeros(treeCount)
    descID              = np.zeros(treeCount)
    massHalo            = np.zeros(treeCount)
    radiusHalo          = np.zeros(treeCount)
    x                   = np.zeros(treeCount)
    y                   = np.zeros(treeCount)
    z                   = np.zeros(treeCount)
    tree = open(pathName+"treeReduced.dat",'r')
    for i in range(treeCount):
        line = tree.readline()
        columns = re.split(r'\s+',line)
        haloExpansionFactor[i] = columns[0]
        haloID             [i] = columns[1]
        descID             [i] = columns[2]
        massHalo           [i] = columns[3]
        radiusHalo         [i] = columns[4]
        x                  [i] = columns[5]
        y                  [i] = columns[6]
        z                  [i] = columns[7]
    tree.close()
    treeRaw = h5py.File(pathName+"treeReduced.hdf5",'w')
    treeRaw.create_dataset('haloExpansionFactor',data=haloExpansionFactor)
    treeRaw.create_dataset('haloID'             ,data=haloID             )
    treeRaw.create_dataset('descID'             ,data=descID             )
    treeRaw.create_dataset('massHalo'           ,data=massHalo           )
    treeRaw.create_dataset('radiusHalo'         ,data=radiusHalo         )
    treeRaw.create_dataset('x'                  ,data=x                  )
    treeRaw.create_dataset('y'                  ,data=y                  )
    treeRaw.create_dataset('z'                  ,data=z                  )
    treeRaw.flush()
    treeRaw.close()
os.remove(pathName+"treeRaw.lock")

# Convert units to Galacticus standards.
print("\tConverting units")
massHalo       /=        hubbleConstant # Convert Msun/h to Msun.
radiusHalo     *= 1.0e-3/hubbleConstant # Convert kpc/h to Mpc.
x              /=        hubbleConstant # Convert Mpc/h to Mpc.
y              /=        hubbleConstant # Convert Mpc/h to Mpc.
z              /=        hubbleConstant # Convert Mpc/h to Mpc.
# Find region center.
selection = (haloExpansionFactor > 0.9999)
xCenter   = np.sum(x[selection]*massHalo[selection])/np.sum(massHalo[selection])
yCenter   = np.sum(y[selection]*massHalo[selection])/np.sum(massHalo[selection])
zCenter   = np.sum(z[selection]*massHalo[selection])/np.sum(massHalo[selection])

# Find focus halo unless we have already been given a host halo ID.
print("\tLocating focus halo")
indexHaloFocus = None
index = np.arange(len(haloID))
if hostHaloID < 0:
    halosCurrent = (haloExpansionFactor > 0.9999) & (np.log10(massHalo) > massHostLogMin) & (np.log10(massHalo) < massHostLogMax)
    if np.count_nonzero(halosCurrent) <= 0:
	# No halo meets the criterion. Simply find the most massive halo at z=0.
        halosPresent = (haloExpansionFactor > 0.9999)
        if np.count_nonzero(halosPresent) <= 0:
            sys.exit("ERROR: unable to locate any viable halos")
        mostMassive = np.argmax(massHalo[halosPresent])
        indexHaloFocus = index[halosPresent][mostMassive]
    else:
        distanceFromCenter     = np.sqrt(+(x-xCenter)**2+(y-yCenter)**2+(z-zCenter)**2)
        indexHalosCurrentFocus = np.argmin(distanceFromCenter[halosCurrent])
        indexHaloFocus         = index[halosCurrent][indexHalosCurrentFocus]
else:
    halosCurrent = (haloID == hostHaloID)
    if np.count_nonzero(halosCurrent) <= 0:
        sys.exit('unable to locate host halo')
    if np.count_nonzero(halosCurrent) >  1:
        sys.exit('ambiguous host halo')
    indexHaloFocus = index[halosCurrent][0]

# Determine the upper mass limit.
massCentral = massHalo[indexHaloFocus]
# Find the primary progenitor at this expansion factor.
print("\tLocating primary progenitor")
indexPrimaryProgenitor = indexHaloFocus
while abs(haloExpansionFactor[indexPrimaryProgenitor]-expansionFactor) > 1.0e-3:
    selection = (descID == haloID[indexPrimaryProgenitor])
    if np.count_nonzero(selection) < 1:
        sys.exit('no progenitor found' )
    elif np.count_nonzero(selection) > 1:
        sys.exit('ambiguous progenitor')
    else:
        indexPrimaryProgenitor = index[selection][0]

# Store the relevant data.
primaryHaloData = ET.Element('primaryHalo')
primaryHaloData.set('i' ,str(int(haloID    [indexPrimaryProgenitor])))
primaryHaloData.set('x' ,str(    x         [indexPrimaryProgenitor]) )
primaryHaloData.set('y' ,str(    y         [indexPrimaryProgenitor]) )
primaryHaloData.set('z' ,str(    z         [indexPrimaryProgenitor]) )
primaryHaloData.set('r' ,str(    radiusHalo[indexPrimaryProgenitor]) )
primaryHaloData.set('m' ,str(    massHalo  [indexPrimaryProgenitor]) )
primaryHaloData.set('l' ,str(    boxSize                           ) )
primaryHaloData.set('xc',str(    xCenter                           ) )
primaryHaloData.set('yc',str(    yCenter                           ) )
primaryHaloData.set('zc',str(    zCenter                           ) )
primaryHaloData.set('mc',str(    massCentral                       ) )
primaryHaloFile = open(primaryHaloFileName,'w')
primaryHaloFile.write(ET.tostring(primaryHaloData, encoding='unicode'))
primaryHaloFile.close()
