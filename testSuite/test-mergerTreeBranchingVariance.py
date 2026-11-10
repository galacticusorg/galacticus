#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test the branchingVariance merger tree operator, which decomposes the variance of counts in merger trees into contributions from
# individual splits. A small set of trees is processed in "table" mode (which also exports the trees) and then in "jumps" mode,
# and the analysis script, scripts/analysis/mergerTreeVariance.py, is run on the exported trees. The test checks that:
#
#  * the operator reproduces the analysis script (its reference implementation) on the same trees, for every statistic and
#    threshold;
#  * the identity N = n₁(root) + ΣΔ + Σdrift holds for every tree, and Var(N) = <ΣΔ²> holds within its error;
#  * splits contribute only where they should: on the main branch for statistic C, before the target epoch for statistic A;
#  * the analysis script reproduces the exact result for a Yule (pure-birth) process, Var(N) = <N>(<N>-1), so that the chain of
#    comparisons is anchored to an analytic result.
# Andrew Benson, Claude (03-November-2026)

failed = False
def fail(message):
    global failed
    print("FAILED: "+message)
    failed = True

def run(label,command):
    # Run a command from the root of the repository, logging its output, and report whether it succeeded.
    print("Running "+label+"...")
    with open("outputs/mergerTreeBranchingVariance"+label+".log","w") as log:
        status = subprocess.run("cd ..; "+command,stdout=log,stderr=log,shell=True)
    print("...done ("+str(status.returncode)+")")
    if status.returncode != 0 or subprocess.run("grep -q -i -e fatal -e aborted outputs/mergerTreeBranchingVariance"+label+".log",shell=True).returncode == 0:
        print("FAILED: "+label+" run:")
        subprocess.run("tail -n 50 outputs/mergerTreeBranchingVariance"+label+".log",shell=True)
        sys.exit(0)

subprocess.run("mkdir -p outputs",shell=True)
subprocess.run("rm -f outputs/mergerTreeBranchingVariance*.hdf5",shell=True)

# Run the operator in both modes, and the analysis script on the exported trees. The script must assign trees to the same folds
# as the operator (by seed, which is the base seed, 220, plus the tree index), use the same binning and ω, and evaluate the
# cumulative variance on the same grid. (Within one pair of runs a uniform shift of the fold labels changes nothing, since each
# tree is evaluated with the tables of the other two folds whatever their labels, so this test does not check that folds follow
# the seed rather than the tree index. That matters only when the two runs use different base seeds.)
run("Table","./Galacticus.exe testSuite/parameters/mergerTreeBranchingVarianceTable.xml")
run("Jumps","./Galacticus.exe testSuite/parameters/mergerTreeBranchingVarianceJumps.xml")
run("Script","python3 scripts/analysis/mergerTreeVariance.py testSuite/outputs/mergerTreeBranchingVarianceTrees.hdf5 --statistics ABC --massThresholds 1.0e10,3.0e10 --redshiftTarget 1.0 --binWidthMass 0.0125 --binWidthOmega 0.025 --tableMethod linear --crossProduct --bootstrapCount 200 --foldOffset 220 --cumulativeOmegaTable testSuite/outputs/mergerTreeBranchingVarianceJumps.hdf5 --outputFileName testSuite/outputs/mergerTreeBranchingVarianceScript.hdf5")

operator    = h5py.File("outputs/mergerTreeBranchingVarianceJumps.hdf5" ,"r")["branchingVariance"]
script      = h5py.File("outputs/mergerTreeBranchingVarianceScript.hdf5","r")
forestIndex = h5py.File("outputs/mergerTreeBranchingVarianceTrees.hdf5" ,"r")["forestIndex/forestIndex"][:].astype(np.int64)
# Trees are written by the operator in the order in which they were processed, so sort them by index.
order       = np.argsort(operator["trees/treeIndex"][:])
fold        = operator["trees/treeSeed"][:][order] % 3
omegaTarget = operator.attrs["omegaTarget"]
omegaGrid   = operator["omegaCumulative"][:]

def covarianceWithinFolds(x,y):
    # The cumulative variance is the covariance of the expectations evaluated with the tables of the two other folds, taken within
    # each fold and combined weighted by the number of trees (as in the analysis script).
    counts = np.array([np.count_nonzero(fold == h) for h in range(3)])
    return np.dot(counts,[np.cov(x[fold == h],y[fold == h])[0,1] for h in range(3)])/counts.sum()

# Tolerances. Per-tree counts are integers, and must agree exactly. The script interpolates ω(z) from a tabulation, while the
# operator evaluates it, so jumps differ slightly (by ~1e-4 at most, with typical values of order unity); aggregates over all
# trees and splits agree far more closely.
toleranceAggregate = 1.0e-5
toleranceSplit     = 2.0e-3
toleranceTelescope = 1.0e-10
for statistic in "ABC":
    for indexThreshold in range(2):
        name = "statistic"+statistic+"/threshold"+str(indexThreshold)
        o    = operator[name]
        p    = script  [name]
        count                         = o["countTotal"                   ][:][order]
        sumDeltaSquared               = o["sumDeltaSquared"              ][:][order]
        sumDeltaSquaredAfterFormation = o["sumDeltaSquaredAfterFormation"][:][order]
        print(name+":")
        # Comparison with the script. The sum of Δ² after formation includes the split at which the main branch falls below half
        # of the root mass, which often carries a substantial part of it, so the comparison checks that both place the cut at the
        # same node.
        differences = {
            "<N>"                  : abs(count                        .mean()      -p.attrs["countMean"          ]),
            "Var(N)"               : abs(count                        .var(ddof=1) -p.attrs["countVariance"      ]),
            "<ΣΔ²>"                : abs(sumDeltaSquared              .mean()      -p.attrs["sumDeltaSquaredMean"])/p.attrs["sumDeltaSquaredMean"],
            "<ΣΔ² after formation>": abs(sumDeltaSquaredAfterFormation.mean()      -p["trees/sumDeltaSquaredAfterFormation"][:].mean())/p.attrs["sumDeltaSquaredMean"],
        }
        expectation          = o["cumulativeExpectation"         ][:][order,:].T
        expectationAlternate = o["cumulativeExpectationAlternate"][:][order,:].T
        sumBefore            = o["cumulativeSumDeltaSquared"     ][:][order,:].T
        differences["cumulative Var(E[N|ω])"] = np.max(np.abs(np.array([covarianceWithinFolds(expectation[k],expectationAlternate[k]) for k in range(len(omegaGrid))])-p["cumulative/varianceConditionalMean"][:]))
        differences["cumulative <ΣΔ²>"      ] = np.max(np.abs(sumBefore.mean(axis=1)-p["cumulative/sumDeltaSquared"][:]))
        for quantity, difference in differences.items():
            print("   {:24s} difference {:.2e}".format(quantity,difference))
        if differences["<N>"] != 0.0:
            fail(name+": <N> differs from that of the analysis script")
        if differences["Var(N)"] > toleranceAggregate*p.attrs["countVariance"]:
            fail(name+": Var(N) differs from that of the analysis script")
        for quantity in ("<ΣΔ²>","<ΣΔ² after formation>","cumulative Var(E[N|ω])","cumulative <ΣΔ²>"):
            if differences[quantity] > toleranceAggregate*(1.0 if quantity in ("<ΣΔ²>","<ΣΔ² after formation>") else p.attrs["sumDeltaSquaredMean"]):
                fail(name+": "+quantity+" differs from that of the analysis script")
        # Per-split values, matched by tree, ω, and mass.
        splitsOperator = o["splits"]
        splitsScript   = p["splits"]
        treeScript     = forestIndex[splitsScript["tree"][:]]
        keyOperator    = np.lexsort((splitsOperator["mass"][:],splitsOperator["omega"][:],splitsOperator["treeIndex"][:]))
        keyScript      = np.lexsort((splitsScript  ["mass"][:],splitsScript  ["omega"][:],treeScript                  ))
        if len(keyOperator) != len(keyScript) or not np.array_equal(splitsOperator["treeIndex"][:][keyOperator],treeScript[keyScript]):
            fail(name+": splits differ from those of the analysis script ({:d} vs. {:d})".format(len(keyOperator),len(keyScript)))
        else:
            for quantity in ("delta","deltaSquared","n1Parent","countAdded","isMainBranch","hierarchyDepth"):
                difference = np.max(np.abs(splitsOperator[quantity][:][keyOperator]-splitsScript[quantity][:][keyScript]),initial=0.0)
                print("   split {:18s} difference {:.2e}".format(quantity,difference))
                if difference > toleranceSplit:
                    fail(name+": per-split "+quantity+" differs from that of the analysis script")
        # Identities.
        telescope = np.max(np.abs(o["telescopeResidual"][:]))
        if telescope > toleranceTelescope:
            fail(name+": N = n₁(root)+ΣΔ+Σdrift fails by {:.2e}".format(telescope))
        varianceMinusSum = p.attrs["countVariance"]-p.attrs["sumDeltaSquaredMean"]
        error            = p.attrs["varianceMinusSumDeltaSquaredError"]
        print("   Var(N)-<ΣΔ²> = {:+.4f} ± {:.4f}".format(varianceMinusSum,error))
        if abs(varianceMinusSum) > 4.0*error:
            fail(name+": Var(N) = <ΣΔ²> fails by more than 4σ")
        if statistic == "C" and np.any(splitsOperator["isMainBranch"][:] != 1):
            fail(name+": statistic C has jumps from side-branch splits")
        if statistic == "A":
            if np.any(splitsOperator["omega"][:] >= omegaTarget):
                fail(name+": statistic A has jumps from splits after the target epoch")
            beyond = omegaGrid >= omegaTarget
            if not np.all(expectation[beyond] == count):
                fail(name+": statistic A expectation differs from N at or beyond the target epoch")
if not failed:
    print("SUCCESS: branchingVariance operator agrees with the analysis script, and its identities hold")

# The Yule (pure-birth) limit: every branch splits into two branches of the same mass at rate β per unit ω, with no accretion.
# The number of branches alive at ω₁ is then geometric, with <N> = exp(β(ω₁-ω₀)) and Var(N) = <N>(<N>-1). Trees are generated in
# the export format with the relation between ω and redshift used by the analysis script (δc0/D(z)), and analyzed as statistic A.
countTrees     = 20000
beta           = 2.0
redshiftTarget = 1.0
OmegaMatter    = 0.3153
OmegaLambda    = 0.6847
expansionFactor = np.logspace(-4.0,0.0,4001)
hubble          = np.sqrt(OmegaMatter/expansionFactor**3+OmegaLambda)
integrand       = 1.0/(expansionFactor*hubble)**3
growth          = hubble*(np.concatenate(([0.0],np.cumsum(0.5*(integrand[1:]+integrand[:-1])*np.diff(expansionFactor))))+expansionFactor[0]**2.5/2.5/OmegaMatter**1.5)
growth         /= growth[-1]
def omegaOfRedshift(redshift):
    return 1.686/np.interp(np.log(1.0/(1.0+redshift)),np.log(expansionFactor),growth)
redshifts   = np.linspace(0.0,20.0,200001)
omegas      = omegaOfRedshift(redshifts)
omegaStart  = omegaOfRedshift(0.0)
omegaFinal  = omegaOfRedshift(redshiftTarget)+0.5
rng         = np.random.default_rng(14)
nodeIndex, descendantIndex, omegaNode, firstNode, numberOfNodes = [], [], [], [], []
for tree in range(countTrees):
    index, descendant, omega = [1], [-1], [omegaStart]
    stack     = [(1,omegaStart)]
    indexNext = 2
    while stack:
        node, omegaNode_ = stack.pop()
        omegaSplit = omegaNode_+rng.exponential(1.0/beta)
        if omegaSplit >= omegaFinal:
            # End the branch with a final node beyond the target epoch.
            index.append(indexNext); descendant.append(node); omega.append(omegaFinal); indexNext += 1
            continue
        for i in range(2):
            index.append(indexNext); descendant.append(node); omega.append(omegaSplit)
            stack.append((indexNext,omegaSplit)); indexNext += 1
    firstNode    .append(len(nodeIndex))
    numberOfNodes.append(len(index    ))
    nodeIndex      .extend(index     )
    descendantIndex.extend(descendant)
    omegaNode      .extend(omega     )
with h5py.File("outputs/mergerTreeBranchingVarianceYuleTrees.hdf5","w") as treeFile:
    treeFile.create_dataset("forestHalos/nodeIndex"      ,data=np.array(nodeIndex      ))
    treeFile.create_dataset("forestHalos/descendantIndex",data=np.array(descendantIndex))
    treeFile.create_dataset("forestHalos/nodeMass"       ,data=np.full(len(nodeIndex),1.0e12))
    treeFile.create_dataset("forestHalos/redshift"       ,data=np.interp(np.array(omegaNode),omegas,redshifts))
    treeFile.create_dataset("forestIndex/firstNode"      ,data=np.array(firstNode      ))
    treeFile.create_dataset("forestIndex/numberOfNodes"  ,data=np.array(numberOfNodes  ))
    treeFile.create_group("cosmology")
    treeFile["cosmology"].attrs["OmegaMatter"] = OmegaMatter
    treeFile["cosmology"].attrs["OmegaLambda"] = OmegaLambda
run("Yule","python3 scripts/analysis/mergerTreeVariance.py testSuite/outputs/mergerTreeBranchingVarianceYuleTrees.hdf5 --statistics A --massThresholds 1.0e11 --redshiftTarget 1.0 --binWidthMass 0.0125 --binWidthOmega 0.025 --tableMethod linear --crossProduct --bootstrapCount 200 --segmentEndMass node --outputFileName testSuite/outputs/mergerTreeBranchingVarianceYule.hdf5")
yule         = h5py.File("outputs/mergerTreeBranchingVarianceYule.hdf5","r")["statisticA/threshold0"]
countMean    = np.exp(beta*(yule.file.attrs["omegaTarget"]-omegaStart))
varianceYule = countMean*(countMean-1.0)
# The difference of <ΣΔ²> from the exact variance is its difference from the sample variance, whose error the script gives, plus
# the difference of the sample variance from the exact variance. For a geometric distribution with large <N> (nearly
# exponential, with kurtosis 9) the error of the sample variance is close to √(8/countTrees) times the variance.
errorSum = np.sqrt(yule.attrs["varianceMinusSumDeltaSquaredError"]**2+8.0*varianceYule**2/countTrees)
print("Yule limit: <N> = {:.3f} (exact {:.3f}), Var(N) = {:.3f}, <ΣΔ²> = {:.3f} ± {:.3f} (exact {:.3f})".format(yule.attrs["countMean"],countMean,yule.attrs["countVariance"],yule.attrs["sumDeltaSquaredMean"],errorSum,varianceYule))
if abs(yule.attrs["countMean"]-countMean) > 4.0*np.sqrt(varianceYule/countTrees):
    fail("Yule limit: <N> differs from the exact result")
elif abs(yule.attrs["sumDeltaSquaredMean"]-varianceYule) > 4.0*errorSum:
    fail("Yule limit: <ΣΔ²> differs from the exact Var(N)")
else:
    print("SUCCESS: analysis script reproduces the Yule limit")
sys.exit(0)
