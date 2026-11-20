#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys
import h5py
import numpy as np

# Check construction of constrained merger trees.
# Andrew Benson (ported to Python)

# Statistics of constrained branches. A branch is described by the logarithm of its mass at a set of redshifts below the
# constraint, and by the number of halos merging directly into it below the constraint, in bins of their mass. For constrained
# trees the branch is the constrained branch. For unconstrained trees (used to make the reference) every progenitor within
# ±0.1 dex of the constrained mass at the constrained redshift is selected, and traced to the root: this is the rejection
# method, which samples branches with the correct distribution.
redshiftsBranch        = np.array([0.25,0.5,1.0,1.5,2.0,3.0])
binsMassSecondary      = np.array([10.0,10.5,11.0,11.5,12.0,15.0])
halfWidthSelection     = np.log(10.0)*0.1

def treeLinks(tree):
    # Return masses, redshifts, the index of each node's parent (-1 for none), and a dictionary of each node's children.
    index    = tree["nodeIndex"  ][:]
    mass     = tree["massBasic"  ][:]
    redshift = tree["redshift"   ][:]
    position = {int(nodeIndex): i for i, nodeIndex in enumerate(index)}
    parent   = np.array([position.get(int(parentIndex),-1) for parentIndex in tree["parentIndex"][:]])
    children = {}
    for i, j in enumerate(parent):
        if j >= 0:
            children.setdefault(int(j),[]).append(i)
    return mass, redshift, parent, children

def branchStatistics(redshifts, masses, secondaries, mass):
    order          = np.argsort(redshifts)
    massBranch     = np.interp(redshiftsBranch,np.array(redshifts)[order],np.log10(np.array(masses)[order]))
    countSecondary = np.histogram(np.log10(mass[secondaries]),binsMassSecondary)[0] if len(secondaries) > 0 else np.zeros(len(binsMassSecondary)-1)
    return massBranch, countSecondary

def branchesConstrained(fileName, massConstraint, redshiftConstraint):
    # Statistics of the constrained branch of each tree. Also returns the number of trees in which the branch does not pass
    # through the constrained mass. (The node of the constrained mass is placed slightly before the constrained redshift, at
    # z=4.005 for z=4, so the redshift is matched only to 1%.)
    branches = []; secondaries = []; missed = 0
    with h5py.File(fileName,"r") as model:
        trees = model["mergerTreeStructures"]
        for treeName in trees.keys():
            tree                             = trees[treeName]
            mass, redshift, parent, children = treeLinks(tree)
            constrained = np.where(tree["nodeIsConstrained"][:] == 1)[0]
            below       = constrained[redshift[constrained] < redshiftConstraint*(1.0-1.0e-6)]
            atConstraint= constrained[(np.abs(redshift[constrained]-redshiftConstraint) < 1.0e-2*redshiftConstraint) & (np.abs(mass[constrained]/massConstraint-1.0) < 1.0e-3)]
            if len(atConstraint) == 0:
                missed += 1
                continue
            onBranch = set(constrained.tolist())
            merging  = [k for j in below for k in children.get(int(j),[]) if k not in onBranch]
            b, s     = branchStatistics(list(redshift[below])+[redshiftConstraint],list(mass[below])+[massConstraint],merging,mass)
            branches.append(b); secondaries.append(s)
    return np.array(branches), np.array(secondaries), missed

def branchesRejection(fileName, massConstraint, redshiftConstraint, offset):
    # Statistics of all branches through the constrained mass and redshift, found by rejection in unconstrained trees. Returns
    # also the tree (numbered from the given offset) in which each branch was found.
    branches = []; secondaries = []; treeOf = []
    with h5py.File(fileName,"r") as model:
        trees = model["mergerTreeStructures"]
        for treeNumber, treeName in enumerate(trees.keys()):
            mass, redshift, parent, children = treeLinks(trees[treeName])
            primary = {j: max(c,key=lambda k: mass[k]) for j, c in children.items()}
            for i in np.where((parent >= 0) & (redshift >= redshiftConstraint))[0]:
                j = parent[i]
                if redshift[j] >= redshiftConstraint:
                    continue
                # Mass at the constrained redshift: interpolated along a primary progenitor, otherwise that of the merging halo.
                if primary[int(j)] == i:
                    weight      = (redshift[i]-redshiftConstraint)/(redshift[i]-redshift[j])
                    massCrossing= np.exp((1.0-weight)*np.log(mass[i])+weight*np.log(mass[j]))
                else:
                    massCrossing= mass[i]
                if abs(np.log(massCrossing/massConstraint)) > halfWidthSelection:
                    continue
                path = []
                while j >= 0:
                    path.append(j); j = parent[j]
                onBranch = set(path+[i])
                merging  = [k for j in path for k in children.get(int(j),[]) if k not in onBranch]
                b, s     = branchStatistics([redshiftConstraint]+list(redshift[path]),[massCrossing]+list(mass[path]),merging,mass)
                branches.append(b); secondaries.append(s); treeOf.append(offset+treeNumber)
    return branches, secondaries, treeOf

def compareBranches(fileName, reference):
    # Compare the statistics of constrained branches with those found by rejection in unconstrained trees, printing the
    # outcome of each check. The allowed offsets combine a fixed margin with three times the combined statistical error. Over
    # 640 trees the constraint-conditioned branching rates match rejection to within about 0.01 dex in branch mass and 2% in
    # the number of merging halos, so the fixed margin is a safeguard, not an allowance for a known bias.
    toleranceMassBranch = 0.1
    toleranceSecondary  = 0.1
    significance        = 3.0
    branches, secondaries, missed = branchesConstrained(fileName,reference["massConstraint"],reference["redshiftConstraint"])
    statusStr = "SUCCESS" if missed == 0 and len(branches) > 0 else "FAILED"
    print(f"{statusStr}: constrained branches pass through the constraint ({missed} of {len(branches)+missed} do not)")
    if len(branches) < 2:
        return
    # Mass of the branch.
    offset    = branches.mean(axis=0)-np.array(reference["massBranchMean"])
    error     = np.sqrt(branches.std(axis=0,ddof=1)**2/len(branches)+np.array(reference["massBranchError"])**2)
    allowed   = toleranceMassBranch+significance*error
    statusStr = "SUCCESS" if np.all(np.abs(offset) <= allowed) else "FAILED"
    print(f"{statusStr}: mass of constrained branches matches rejection")
    for z, d, a in zip(reference["redshifts"],offset,allowed):
        print(f"    z = {z:4.2f}: Δlog₁₀M = {d:+6.3f} (allowed ±{a:5.3f})")
    # Halos merging into the branch.
    meanReference = np.array(reference["secondaryMean"])
    ratio         = secondaries.mean(axis=0)/meanReference
    error         = ratio*np.sqrt((secondaries.std(axis=0,ddof=1)/np.sqrt(len(secondaries))/secondaries.mean(axis=0))**2+(np.array(reference["secondaryError"])/meanReference)**2)
    allowed       = toleranceSecondary+significance*error
    statusStr     = "SUCCESS" if np.all(np.abs(ratio-1.0) <= allowed) else "FAILED"
    print(f"{statusStr}: halos merging into constrained branches match rejection")
    for low, high, r, a in zip(reference["binsMass"][:-1],reference["binsMass"][1:],ratio,allowed):
        print(f"    log₁₀M = {low:4.1f}-{high:4.1f}: ratio = {r:5.3f} (allowed 1±{a:5.3f})")

# Make the reference from rejection runs if requested.
parser = argparse.ArgumentParser()
parser.add_argument("--reference",nargs="+",help="Make the reference summary from the given models of constrainedMergerTreesRejection.xml, then exit.")
arguments = parser.parse_args()
if arguments.reference is not None:
    branches = []; secondaries = []; treeOf = []; countTrees = 0
    for fileName in arguments.reference:
        b, s, t = branchesRejection(fileName,1.0e12,4.0,countTrees)
        branches += b; secondaries += s; treeOf += t
        with h5py.File(fileName,"r") as model:
            countTrees += len(model["mergerTreeStructures"].keys())
    branches = np.array(branches); secondaries = np.array(secondaries); treeOf = np.array(treeOf)
    # Errors are found by bootstrap resampling of trees, as branches found in the same tree are correlated.
    generator      = np.random.default_rng(4321)
    branchesOfTree = [np.where(treeOf == tree)[0] for tree in range(countTrees)]
    resampled      = [np.concatenate([branchesOfTree[k] for k in generator.integers(0,countTrees,countTrees)]) for _ in range(500)]
    reference      = {
        "description"       : "Statistics of branches passing through a progenitor of mass 10^12 Msun at z=4 in trees of 10^14 Msun at z=0, found by rejection in unconstrained trees (testSuite/parameters/constrainedMergerTreesRejection.xml) and made by test-constrained-merger-trees.py.",
        "massConstraint"    : 1.0e12,
        "redshiftConstraint": 4.0,
        "countTrees"        : countTrees,
        "countBranches"     : int(len(branches)),
        "redshifts"         : redshiftsBranch.tolist(),
        "binsMass"          : binsMassSecondary.tolist(),
        "massBranchMean"    : branches.mean(axis=0).tolist(),
        "massBranchError"   : np.std([branches[r].mean(axis=0) for r in resampled],axis=0).tolist(),
        "secondaryMean"     : secondaries.mean(axis=0).tolist(),
        "secondaryError"    : np.std([secondaries[r].mean(axis=0) for r in resampled],axis=0).tolist()
    }
    with open("data/constrainedMergerTreesRejection.json","w") as referenceFile:
        json.dump(reference,referenceFile,indent=1)
    sys.exit(0)

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Run the model.
with open("outputs/constrainedMergerTrees.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/constrainedMergerTrees.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED:  model run:")
    with open("outputs/constrainedMergerTrees.log") as f:
        print(f.read())
else:
    print("SUCCESS: model run")

# Define the constraint.
redshiftConstraint = 8.0
massConstraint     = 1.0e12

# Read the model data and check constraint is met.
with h5py.File("outputs/constrainedMergerTrees.hdf5", "r") as model:
    trees = model["mergerTreeStructures"]
    massesClosest = []
    for treeName in trees.keys():
        tree                  = trees[treeName]
        isOnConstrainedBranch = tree["nodeIsConstrained"][:]
        redshift              = tree["redshift"][:]
        mass                  = tree["massBasic"][:]
        constrainedBranch     = np.where((isOnConstrainedBranch == 1) & (redshift >= redshiftConstraint))[0]
        if len(constrainedBranch) > 0:
            redshiftDelta   = np.abs(redshift[constrainedBranch] - redshiftConstraint)
            redshiftClosest = np.argmin(redshiftDelta)
            massesClosest.append(mass[constrainedBranch[redshiftClosest]])

    massesClosest = np.array(massesClosest)
    tolerance     = 1.0e-3
    status_str    = "SUCCESS" if np.all(np.abs(massesClosest - massConstraint) < massConstraint * tolerance) else "FAILED"
    print(f"{status_str}: constrained merger trees")

    # Check propagation of constrained branch indicator.
    output        = model["Outputs/Output1"]
    nodes         = output["nodeData"]
    treeIndices   = output["mergerTreeIndex"][:]
    treeIndex     = nodes["mergerTreeIndex"][:]
    isConstrained = nodes["nodeIsConstrained"][:]
    statusIndicator = "SUCCESS"
    for i in range(len(treeIndices)):
        selection = np.where((treeIndex == treeIndices[i]) & (isConstrained == 1))[0]
        if len(selection) != 1:
            statusIndicator = "FAILED"
            break
    print(f"{statusIndicator}: constrained branch indicator")

# Check the statistics of constrained branches against those found by rejection in unconstrained trees.
with open("outputs/constrainedMergerTreesRegression.log", "w") as logFile:
    status = subprocess.run(
        "cd ..; ./Galacticus.exe testSuite/parameters/constrainedMergerTreesRegression.xml",
        shell=True, stdout=logFile, stderr=subprocess.STDOUT
    )
if status.returncode != 0:
    print("FAILED: model run (constrained branch statistics):")
    with open("outputs/constrainedMergerTreesRegression.log") as f:
        print(f.read())
    sys.exit(0)
print("SUCCESS: model run (constrained branch statistics)")
with open("data/constrainedMergerTreesRejection.json") as referenceFile:
    reference = json.load(referenceFile)
compareBranches("outputs/constrainedMergerTreesRegression.hdf5",reference)
