#!/usr/bin/env python3
import os
import sys
import argparse
# The only linear algebra is small matrix products for bootstrapping, for which multithreaded BLAS spends far more time spinning
# than computing.
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import h5py

# Measure the split-by-split decomposition of the halo-to-halo variance of progenitor and subhalo counts in merger trees.
#
# A merger tree built by the Cole et al. (2000) algorithm is a Markov branching process in the variable ω=δc/D. For a count
# statistic N, let n₁(M,ω) be the expected contribution to N of a single branch of mass M at ω. Then E[N|tree built to ω] is a
# martingale whose only jumps occur at splits, so that
#
#   Var(N) = E[ Σ_splits Δ² ],  with  Δ = Σ_progenitors n₁(Mᵢ,ω) [+ any count added by the split] - n₁(Σᵢ Mᵢ,ω).
#
# This script reads trees exported (before evolution) by the "export" merger tree operator in "galacticus" format, measures n₁
# empirically from all nodes of all trees (by the Markov property every node is a valid starting point), evaluates Δ at every
# split, and checks the identity above. Three count statistics are supported:
#
#   A: the number of progenitors with mass ≥ m at a fixed earlier redshift z₁;
#   B: the number of splits anywhere in the tree producing a non-primary progenitor of mass ≥ m (the unevolved subhalo mass
#      function summed over all hierarchy depths);
#   C: as B, but counting only splits of the main branch (direct infall onto the main progenitor).
#
# Along each branch segment n₁ also changes deterministically (the "drift", from smooth accretion and the passage of ω). For each
# tree N = n₁(root) + ΣΔ + Σdrift holds exactly, which checks the bookkeeping; E[ΣΔ + Σdrift] = 0 and Var(N) = E[ΣΔ²] check the
# framework and the n₁ table.
#
# Andrew Benson (20-October-2026)

statisticNames = {
    'A': 'progenitors with mass above threshold at a fixed redshift',
    'B': 'splits producing a non-primary progenitor above threshold, all hierarchy depths',
    'C': 'splits producing a non-primary progenitor above threshold, main branch only'
}

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='mergerTreeVariance.py',description='Decompose the variance of progenitor/subhalo counts in merger trees into contributions from individual splits.')
parser.add_argument('treeFileName'                                                                                                                                                                   )
parser.add_argument('--outputFileName'  ,action='store',default='mergerTreeVariance.hdf5'                    ,help='the HDF5 file to which results are written'                                       )
parser.add_argument('--statistics'      ,action='store',default='ABC'                                        ,help='the statistics to compute, any of "A", "B", "C"'                                  )
parser.add_argument('--massThresholds'  ,action='store',default='1.0e10,3.0e10,1.0e11'                       ,help='comma-separated list of threshold masses, m'                                      )
parser.add_argument('--redshiftTarget'  ,action='store',default=1.0                   ,type=float            ,help='the redshift, z₁, at which progenitors are counted for statistic A'               )
parser.add_argument('--binWidthMass'    ,action='store',default=0.05                  ,type=float            ,help='width (in dex) of mass bins in the n₁ table'                                      )
parser.add_argument('--binWidthOmega'   ,action='store',default=0.05                  ,type=float            ,help='width of ω bins in the n₁ table'                                                  )
parser.add_argument('--crossProduct'    ,action='store_true'                                                 ,help='split the trees into three folds; evaluate Δ for each tree twice, with n₁ tables built from each of the two other folds, and estimate Δ² by the product of the two (whose table noise is independent, so that noise does not inflate the estimate)')
parser.add_argument('--tableMethod'     ,action='store',default='mean'                ,choices=['mean','linear'],help='estimator of n₁ in each table bin: the mean count of its nodes, or a local linear fit evaluated at the bin center')
parser.add_argument('--crossFit'        ,action='store_true'                                                 ,help='evaluate Δ for each half of the trees using an n₁ table built from the other half')
parser.add_argument('--bootstrapCount'  ,action='store',default=1000                  ,type=int              ,help='number of bootstrap resamples used for uncertainties'                             )
parser.add_argument('--deltaCritical0'  ,action='store',default=1.686                 ,type=float            ,help='critical overdensity used to define ω=δc/D(z)'                                    )
parser.add_argument('--segmentEndMass',action='store',default='progenitorSum'   ,choices=['progenitorSum','node'],help='mass of a branch at the end of a segment: the summed progenitor mass (default; accounts for accretion during the segment), or the node mass (for test trees without accretion that do not conserve mass)')
parser.add_argument('--cumulativeCount' ,action='store',default=12                    ,type=int              ,help='number of ω values at which the cumulative variance identity is checked'          )
parser.add_argument('--seed'            ,action='store',default=8783                  ,type=int              ,help='seed for the bootstrap random number generator'                                   )
args = parser.parse_args()

def linearGrowthFactor(redshift,OmegaMatter,OmegaLambda):
    # Compute the linear growth factor, normalized to unity at z=0, for a matter plus cosmological constant universe (the
    # integral solution of Heath 1977, valid for any curvature).
    OmegaCurvature = 1.0-OmegaMatter-OmegaLambda
    expansionFactorGrid = np.logspace(-4.0,0.0,4001)
    def hubble(a):
        return np.sqrt(OmegaMatter/a**3+OmegaCurvature/a**2+OmegaLambda)
    integrand = 1.0/(expansionFactorGrid*hubble(expansionFactorGrid))**3
    integral  = np.concatenate(([0.0],np.cumsum(0.5*(integrand[1:]+integrand[:-1])*np.diff(expansionFactorGrid))))
    # Add the (matter-dominated) contribution from a < 10⁻⁴.
    integral += expansionFactorGrid[0]**2.5/2.5/OmegaMatter**1.5
    growth    = hubble(expansionFactorGrid)*integral
    growth   /= growth[-1]
    expansionFactor = 1.0/(1.0+redshift)
    return np.interp(np.log(expansionFactor),np.log(expansionFactorGrid),growth)

class nodeMeanTable:
    # A table of the mean count per node, binned in log₁₀(mass) and ω, interpolated bilinearly between bin centers. With
    # method="mean" each bin holds the mean count of its nodes, which estimates n₁ at the nodes' centroid rather than at the bin
    # center (a first-order bias where node density varies across the bin). With method="linear" a plane, count ≈ a + b δlog₁₀M +
    # c δω (offsets from the bin center), is fitted by least squares to the nodes in each bin and a is used, giving n₁ at the bin
    # center; bins with too few nodes, or a degenerate spread of nodes, fall back to the mean.
    def __init__(self,logMass,omega,count,logMassRange,omegaRange,binWidthMass,binWidthOmega,method='mean'):
        self.logMassEdges = logMassRange[0]+binWidthMass *np.arange(int(round((logMassRange[1]-logMassRange[0])/binWidthMass ))+1)
        self.omegaEdges   = omegaRange  [0]+binWidthOmega*np.arange(int(round((omegaRange  [1]-omegaRange  [0])/binWidthOmega))+1)
        self.logMassCenters = 0.5*(self.logMassEdges[1:]+self.logMassEdges[:-1])
        self.omegaCenters   = 0.5*(self.omegaEdges  [1:]+self.omegaEdges  [:-1])
        countNodes, _, _ = np.histogram2d(omega,logMass,bins=[self.omegaEdges,self.logMassEdges]                )
        sumCount  , _, _ = np.histogram2d(omega,logMass,bins=[self.omegaEdges,self.logMassEdges],weights=count  )
        sumCount2 , _, _ = np.histogram2d(omega,logMass,bins=[self.omegaEdges,self.logMassEdges],weights=count**2)
        self.countNodes = countNodes
        with np.errstate(invalid='ignore',divide='ignore'):
            self.mean  = sumCount/countNodes
            self.error = np.sqrt(np.maximum(sumCount2/countNodes-self.mean**2,0.0)/countNodes)
        if method == 'linear':
            # Offsets of each node from its bin center, in units of the bin width.
            i  = np.clip(np.searchsorted(self.omegaEdges  ,omega  ,side='right')-1,0,len(self.omegaCenters  )-1)
            j  = np.clip(np.searchsorted(self.logMassEdges,logMass,side='right')-1,0,len(self.logMassCenters)-1)
            dx = (logMass-self.logMassCenters[j])/binWidthMass
            dy = (omega  -self.omegaCenters  [i])/binWidthOmega
            def binSum(weights):
                return np.histogram2d(omega,logMass,bins=[self.omegaEdges,self.logMassEdges],weights=weights)[0]
            Sx , Sy , Sxx, Syy, Sxy = binSum(dx), binSum(dy), binSum(dx*dx), binSum(dy*dy), binSum(dx*dy)
            Sc , Sxc, Syc           = sumCount  , binSum(dx*count), binSum(dy*count)
            matrix = np.stack([np.stack([countNodes,Sx ,Sy ],axis=-1),
                               np.stack([Sx        ,Sxx,Sxy],axis=-1),
                               np.stack([Sy        ,Sxy,Syy],axis=-1)],axis=-2)
            vector = np.stack([Sc,Sxc,Syc],axis=-1)
            # Require enough nodes and a well-conditioned fit (nodes spread across the bin in both directions).
            determinant = np.linalg.det(matrix)
            usable      = (countNodes >= 10) & (determinant > 1.0e-3*np.maximum(countNodes,1.0)**3)
            solution    = np.zeros_like(vector)
            solution[usable] = np.linalg.solve(matrix[usable],vector[usable][...,None])[...,0]
            self.mean   = np.where(usable,solution[...,0],self.mean)
        # Fill empty bins by interpolating in mass along each row of ω (extending end values beyond the populated range). Empty
        # rows are filled from the nearest populated row.
        self.empty  = countNodes == 0
        self.filled = self.mean.copy()
        rowsPopulated = np.nonzero(~np.all(self.empty,axis=1))[0]
        for i in rowsPopulated:
            populated = ~self.empty[i,:]
            self.filled[i,~populated] = np.interp(self.logMassCenters[~populated],self.logMassCenters[populated],self.mean[i,populated])
        for i in np.nonzero(np.all(self.empty,axis=1))[0]:
            self.filled[i,:] = self.filled[rowsPopulated[np.argmin(np.abs(rowsPopulated-i))],:]
        self.evaluationsEmpty = 0

    def __call__(self,logMass,omega):
        x = (logMass-self.logMassCenters[0])/(self.logMassCenters[1]-self.logMassCenters[0])
        y = (omega  -self.omegaCenters  [0])/(self.omegaCenters  [1]-self.omegaCenters  [0])
        x = np.clip(x,0.0,len(self.logMassCenters)-1.0)
        y = np.clip(y,0.0,len(self.omegaCenters  )-1.0)
        i = np.minimum(np.floor(y).astype(np.int64),len(self.omegaCenters  )-2)
        j = np.minimum(np.floor(x).astype(np.int64),len(self.logMassCenters)-2)
        u = y-i
        v = x-j
        self.evaluationsEmpty += np.count_nonzero(self.empty[i,j]|self.empty[i+1,j]|self.empty[i,j+1]|self.empty[i+1,j+1])
        return +(1.0-u)*(1.0-v)*self.filled[i  ,j  ] \
               +     u *(1.0-v)*self.filled[i+1,j  ] \
               +(1.0-u)*     v *self.filled[i  ,j+1] \
               +     u *     v *self.filled[i+1,j+1]

# Read the exported trees.
print('Reading trees from "'+args.treeFileName+'"')
treeFile        = h5py.File(args.treeFileName,'r')
nodeIndex       = treeFile['forestHalos/nodeIndex'      ][:].astype(np.int64)
descendantIndex = treeFile['forestHalos/descendantIndex'][:].astype(np.int64)
mass            = treeFile['forestHalos/nodeMass'       ][:]
redshift        = treeFile['forestHalos/redshift'       ][:]
firstNode       = treeFile['forestIndex/firstNode'      ][:].astype(np.int64)
numberOfNodes   = treeFile['forestIndex/numberOfNodes'  ][:].astype(np.int64)
OmegaMatter     = treeFile['cosmology'].attrs['OmegaMatter']
OmegaLambda     = treeFile['cosmology'].attrs['OmegaLambda']
treeFile.close()
countNodesTotal = len(nodeIndex)
forest          = np.empty(countNodesTotal,dtype=np.int64)
for i in range(len(firstNode)):
    forest[firstNode[i]:firstNode[i]+numberOfNodes[i]] = i

# Compute ω for each node.
omega    = args.deltaCritical0/linearGrowthFactor(redshift,OmegaMatter,OmegaLambda)
logMass  = np.log10(mass)
omega1   = args.deltaCritical0/linearGrowthFactor(args.redshiftTarget,OmegaMatter,OmegaLambda)

# Build the tree structure. Node indices are unique only within a forest, so key on (forest, node index).
keyStride   = nodeIndex.max()+1
key         = forest*keyStride+nodeIndex
order       = np.argsort(key)
parent      = np.full(countNodesTotal,-1,dtype=np.int64)
hasParent   = descendantIndex >= 0
position    = np.searchsorted(key[order],forest[hasParent]*keyStride+descendantIndex[hasParent])
parent[hasParent] = order[position]
if np.any(key[parent[hasParent]] != forest[hasParent]*keyStride+descendantIndex[hasParent]):
    sys.exit('mergerTreeVariance.py: descendant node not found')
children = np.nonzero(hasParent)[0]
# Depth of each node below its root, and the root of each node, by pointer jumping.
depth   = hasParent.astype(np.int64)
jump    = parent.copy()
while True:
    active = np.nonzero(jump >= 0)[0]
    if len(active) == 0:
        break
    depthJump     = depth[jump[active]]
    jumpJump      = jump [jump[active]]
    depth[active] += depthJump
    jump [active]  = jumpJump
root = np.arange(countNodesTotal)
root[hasParent] = parent[hasParent]
while True:
    rootNext = root[root]
    if np.all(rootNext == root):
        break
    root = rootNext
roots      = np.nonzero(~hasParent)[0]
countTrees = len(roots)
treeOfNode = np.searchsorted(roots,root)
# Group nodes by depth.
depthOrder      = np.argsort(depth,kind='stable')
depthMaximum    = depth.max()
depthBoundaries = np.searchsorted(depth[depthOrder],np.arange(depthMaximum+2))
def nodesAtDepth(d):
    return depthOrder[depthBoundaries[d]:depthBoundaries[d+1]]
# Identify the primary (most massive) progenitor of each node.
orderChildren    = children[np.lexsort((mass[children],parent[children]))]
lastOfGroup      = np.append(parent[orderChildren][1:] != parent[orderChildren][:-1],True)
primaryChild     = np.full(countNodesTotal,-1,dtype=np.int64)
primaryChild[parent[orderChildren[lastOfGroup]]] = orderChildren[lastOfGroup]
isPrimary        = np.zeros(countNodesTotal,dtype=bool)
isPrimary[primaryChild[primaryChild >= 0]] = True
countChildren    = np.bincount(parent[children],minlength=countNodesTotal)
massChildren     = np.bincount(parent[children],weights=mass[children],minlength=countNodesTotal)
if args.segmentEndMass == 'node':
    massChildren = mass.copy()
secondary        = children[~isPrimary[children]]
massSecondary    = np.bincount(parent[secondary],weights=mass[secondary],minlength=countNodesTotal)
omegaSplit       = np.where(primaryChild >= 0,omega[np.maximum(primaryChild,0)],np.nan)
# Main-branch membership and hierarchy depth (number of non-primary links between a node and its root).
isMain         = ~hasParent.copy()
hierarchyDepth = np.zeros(countNodesTotal,dtype=np.int64)
for d in range(1,depthMaximum+1):
    nodes = nodesAtDepth(d)
    isMain        [nodes] = isMain        [parent[nodes]] & isPrimary[nodes]
    hierarchyDepth[nodes] = hierarchyDepth[parent[nodes]]+(~isPrimary[nodes])
# Splits are nodes with two or more progenitors.
splits = np.nonzero(countChildren >= 2)[0]
print('  {:d} trees, {:d} nodes, {:d} splits, maximum depth {:d}'.format(countTrees,countNodesTotal,len(splits),depthMaximum))
if np.any(countChildren > 2):
    print('  note: {:d} nodes have more than two progenitors'.format(np.count_nonzero(countChildren > 2)))
# Check that the progenitors of each split share a single ω.
omegaSpread = np.zeros(countNodesTotal)
np.maximum.at(omegaSpread,parent[children],np.abs(omega[children]-omegaSplit[parent[children]]))
print('  maximum ω difference between progenitors of a split: {:.3e}'.format(omegaSpread.max()))
# Formation redshift (main branch first drops below half of the root mass) of each tree.
massRoot          = mass[roots]
redshiftFormation = np.full(countTrees,np.nan)
omegaFormation    = np.full(countTrees,np.nan)
mainNodes         = np.nonzero(isMain)[0]
mainNodes         = mainNodes[np.lexsort((omega[mainNodes],treeOfNode[mainNodes]))]
belowHalf         = mass[mainNodes] < 0.5*massRoot[treeOfNode[mainNodes]]
firstBelow        = mainNodes[belowHalf][np.diff(treeOfNode[mainNodes[belowHalf]],prepend=-1) != 0]
for node in firstBelow:
    tree   = treeOfNode[node]
    before = parent[node]
    factor = (np.log(0.5*massRoot[tree])-np.log(mass[before]))/(np.log(mass[node])-np.log(mass[before]))
    redshiftFormation[tree] = redshift[before]+factor*(redshift[node]-redshift[before])
    omegaFormation   [tree] = omega   [before]+factor*(omega   [node]-omega   [before])

def countSubtree(own,primaryOnly):
    # Sum a per-node count over each node's subtree (or, if primaryOnly, over its own primary-progenitor chain).
    total = own.astype(np.float64).copy()
    for d in range(depthMaximum,0,-1):
        nodes = nodesAtDepth(d)
        if primaryOnly:
            nodes = nodes[isPrimary[nodes]]
        np.add.at(total,parent[nodes],total[nodes])
    return total

def bootstrap(values,function,rng):
    # Return the bootstrap standard deviation of a statistic of per-tree values.
    weights = rng.multinomial(countTrees,np.full(countTrees,1.0/countTrees),size=args.bootstrapCount).astype(np.float64)
    return np.std(function(weights,values))

def weightedVarianceMinusMean(weights,values):
    # Var(N) - <ΣΔ²> for each set of bootstrap weights; values is (N, ΣΔ²).
    total = weights.sum(axis=1)
    mean  = (weights@values[0])/total
    var   = (weights@values[0]**2)/total-mean**2
    var  *= total/(total-1.0)
    return var-(weights@values[1])/total

def weightedVariance(weights,values):
    total = weights.sum(axis=1)
    mean  = (weights@values)/total
    return ((weights@values**2)/total-mean**2)*total/(total-1.0)

def weightedMean(weights,values):
    return (weights@values)/weights.sum(axis=1)

rng = np.random.default_rng(args.seed)
countFolds = 3 if args.crossProduct else (2 if args.crossFit else 1)
massThresholds = [float(m) for m in args.massThresholds.split(',')]
massResolution = mass.min()
print('  minimum node mass {:.3e}; mean root mass {:.3e}'.format(massResolution,massRoot.mean()))
print('  statistic A target: z₁={:.3f}, ω₁={:.4f} (ω₀={:.4f})'.format(args.redshiftTarget,omega1,args.deltaCritical0))
# Ranges of the n₁ table. Bins are aligned so that the largest node mass and the smallest ω (i.e. the roots of trees of fixed
# mass) lie at bin centers, avoiding the bias of clamping at the table edge. At least two bins are used in each direction, as
# needed for interpolation.
countBinsMass  = max(int(np.ceil((logMass.max()-logMass.min())/args.binWidthMass +0.5)),2)
countBinsOmega = max(int(np.ceil((omega  .max()-omega  .min())/args.binWidthOmega+0.5)),2)
logMassRange   = (logMass.max()+0.5*args.binWidthMass -countBinsMass *args.binWidthMass ,logMass.max()+0.5*args.binWidthMass )
omegaRange     = (omega  .min()-0.5*args.binWidthOmega,omega  .min()-0.5*args.binWidthOmega+countBinsOmega*args.binWidthOmega)

output = h5py.File(args.outputFileName,'w')
output.attrs['treeFileName'      ] = args.treeFileName
output.attrs['redshiftTarget'    ] = args.redshiftTarget
output.attrs['omegaTarget'       ] = omega1
output.attrs['crossFit'          ] = args.crossFit
output.attrs['crossProduct'      ] = args.crossProduct
output.attrs['binWidthMass'      ] = args.binWidthMass
output.attrs['binWidthOmega'     ] = args.binWidthOmega
output.attrs['tableMethod'       ] = args.tableMethod
output.attrs['massResolution'    ] = massResolution
treesGroup = output.create_group('trees')
treesGroup.create_dataset('massRoot'         ,data=massRoot         )
treesGroup.create_dataset('redshiftFormation',data=redshiftFormation)
treesGroup.create_dataset('omegaFormation'   ,data=omegaFormation   )

for statistic in args.statistics:
    if statistic not in statisticNames:
        sys.exit('mergerTreeVariance.py: unknown statistic "'+statistic+'"')
    print('\nStatistic '+statistic+': '+statisticNames[statistic])
    print('  {:>10s} {:>9s} {:>9s} {:>9s} {:>19s} {:>9s} {:>9s} {:>8s} {:>8s} {:>9s} {:>9s} {:>9s} {:>9s}'.format('m','<N>','n₁(root)','Var(N)','<ΣΔ²>','Var/<N>','σ_I²','<ΣΔ+Σd>','telescope','M′≥m/Var','M′<m/Var','main/Var','side/Var'))
    for indexThreshold, massThreshold in enumerate(massThresholds):
        # Count of objects added by each split (non-primary progenitors above threshold).
        countAdded = np.bincount(parent[secondary],weights=(mass[secondary] >= massThreshold).astype(np.float64),minlength=countNodesTotal)
        # Per-node counts, and the set of nodes which are valid starting points for the n₁ table.
        if statistic == 'A':
            # A branch segment from node p crosses ω₁ if ω_p < ω₁ ≤ ω_split. Its mass at ω₁ is found by interpolating log mass
            # linearly in ω between the node and the summed progenitor mass at the split.
            crosses = (omega < omega1) & (countChildren > 0) & (omegaSplit >= omega1)
            logMassTarget = np.full(countNodesTotal,-np.inf)
            p = np.nonzero(crosses)[0]
            logMassTarget[p] = logMass[p]+(np.log10(massChildren[p])-logMass[p])*(omega1-omega[p])/(omegaSplit[p]-omega[p])
            own        = (logMassTarget >= np.log10(massThreshold)).astype(np.float64)
            count      = countSubtree(own,primaryOnly=False)
            tableNodes = omega < omega1
            activeSplit = (countChildren >= 2) & (omegaSplit < omega1)
        else:
            own        = countAdded
            count      = countSubtree(own,primaryOnly=(statistic == 'C'))
            tableNodes = np.ones(countNodesTotal,dtype=bool)
            activeSplit = countChildren >= 2
            if statistic == 'C':
                activeSplit &= isMain
        countTree = count[roots]
        # Build the n₁ table(s).
        tables = []
        for h in range(countFolds):
            select = tableNodes if countFolds == 1 else tableNodes & (treeOfNode % countFolds == h)
            tables.append(nodeMeanTable(logMass[select],omega[select],count[select],logMassRange,omegaRange,args.binWidthMass,args.binWidthOmega,method=args.tableMethod))
        def n1(logMass_,omega_,tree_,variant=0):
            # Evaluate n₁ using the appropriate table: with folds, a table built from another fold of the trees (variant 0 or 1
            # selects which of the other folds, for the cross-product estimator).
            if countFolds == 1:
                return tables[0](logMass_,omega_)
            result = np.empty(len(logMass_))
            for h in range(countFolds):
                s = (tree_ % countFolds) == h
                result[s] = tables[(h+1+variant) % countFolds](logMass_[s],omega_[s])
            return result
        # n₁ at each node, and at the end state of each segment (the summed progenitor mass at the split, or the branch end).
        segments = np.nonzero(countChildren > 0)[0]
        def n1NodeEnd(variant):
            n1Node_ = n1(logMass,omega,treeOfNode,variant)
            n1End_  = np.zeros(countNodesTotal)
            n1End_[segments] = n1(np.log10(massChildren[segments]),omegaSplit[segments],treeOfNode[segments],variant)
            return n1Node_, n1End_
        n1Node, n1End = n1NodeEnd(0)
        # Markov check: by the Markov property the mean of (count - n₁) must vanish for every kind of node (roots, the primary and
        # non-primary progenitors created at a split, and nodes created by accretion alone).
        nodeTypes = (('root'              ,~hasParent                                                       ),
                     ('primary at split'  , hasParent &  isPrimary & (countChildren[np.maximum(parent,0)] >= 2)),
                     ('secondary at split', hasParent & ~isPrimary                                          ),
                     ('accretion'         , hasParent & (countChildren[np.maximum(parent,0)] == 1)          ))
        markovCheck = []
        for name, select in nodeTypes:
            select   = select & tableNodes & (n1Node > 0.0)
            residual = count[select]-n1Node[select]
            if len(residual) > 1:
                markovCheck.append((name,len(residual),residual.mean(),residual.std()/np.sqrt(len(residual))))
            else:
                markovCheck.append((name,len(residual),np.nan       ,np.nan                                   ))
        # Jumps at splits.
        s = np.nonzero(activeSplit)[0]
        def jumps(n1Node_,n1End_):
            delta_ = np.zeros(countNodesTotal)
            if statistic == 'C':
                # Role-dependent n₁: only the primary progenitor carries the main branch forward.
                n1Children = np.zeros(countNodesTotal)
                n1Children[segments] = n1Node_[primaryChild[segments]]
            else:
                n1Children = np.bincount(parent[children],weights=n1Node_[children],minlength=countNodesTotal)
            delta_[s] = n1Children[s]+(countAdded[s] if statistic != 'A' else 0.0)-n1End_[s]
            return delta_
        delta = jumps(n1Node,n1End)
        # Deterministic drift of n₁ along each segment.
        drift = np.zeros(countNodesTotal)
        if statistic == 'A':
            live = omega < omega1
            # Segments ending at a split before ω₁, or terminating (leaf) before ω₁.
            ends = live & ~crosses
            drift[ends] = np.where(countChildren[ends] > 0,n1End[ends],0.0)-n1Node[ends]
            # Segments crossing ω₁ end with the count of the surviving branch.
            drift[crosses] = own[crosses]-n1Node[crosses]
        elif statistic == 'B':
            drift = np.where(countChildren > 0,n1End,0.0)-n1Node
        else:
            drift[isMain] = np.where(countChildren[isMain] > 0,n1End[isMain],0.0)-n1Node[isMain]
        # Per-tree sums.
        terminal     = massSecondary >= massThreshold
        if args.crossProduct:
            # Product of jumps evaluated with tables from the two other folds: the table noise in the two is independent, so its
            # variance does not contribute to the expectation.
            deltaSquared = delta*jumps(*n1NodeEnd(1))
        else:
            deltaSquared = delta**2
        def treeSum(values,select=None):
            if select is None:
                return np.bincount(treeOfNode,weights=values,minlength=countTrees)
            return np.bincount(treeOfNode[select],weights=values[select],minlength=countTrees)
        sumDelta                  = treeSum(delta                           )
        sumDeltaSquared           = treeSum(deltaSquared                    )
        sumDeltaSquaredTerminal   = treeSum(deltaSquared, terminal          )
        sumDeltaSquaredStructural = treeSum(deltaSquared,~terminal          )
        sumDeltaSquaredMain       = treeSum(deltaSquared, isMain            )
        sumDeltaSquaredSide       = treeSum(deltaSquared,~isMain            )
        sumDrift                  = treeSum(drift                           )
        beforeFormation           = omegaSplit < omegaFormation[treeOfNode]
        sumDeltaSquaredBeforeFormation = treeSum(deltaSquared,beforeFormation)
        n1Root                    = n1Node[roots]
        telescopeResidual         = countTree-n1Root-sumDelta-sumDrift
        meanCount = countTree.mean()
        # Cumulative variance: at each ω on a grid, Var(E[N|tree built to ω]) should equal <Σ_{ω_split<ω} Δ²>. E[N|tree built to
        # ω] is the count already realized (for statistics B, C) plus n₁ summed over the branches alive at ω.
        omegaGrid           = np.linspace(omega.min(),omega1 if statistic == 'A' else np.quantile(omegaSplit[s],0.999),args.cumulativeCount)
        varianceConditional = np.empty(args.cumulativeCount)
        meanConditional     = np.empty(args.cumulativeCount)
        meanConditionalError= np.empty(args.cumulativeCount)
        sumDeltaSquaredCumulative = np.empty(args.cumulativeCount)
        errorCumulative     = np.empty(args.cumulativeCount)
        for i, omegaG in enumerate(omegaGrid):
            alive = (countChildren > 0) & (omega <= omegaG) & (omegaSplit > omegaG)
            if statistic == 'C':
                alive &= isMain
            p = np.nonzero(alive)[0]
            logMassG = logMass[p]+(np.log10(massChildren[p])-logMass[p])*(omegaG-omega[p])/(omegaSplit[p]-omega[p])
            expectation = np.bincount(treeOfNode[p],weights=n1(logMassG,np.full(len(p),omegaG),treeOfNode[p]),minlength=countTrees)
            before      = activeSplit & (omegaSplit <= omegaG)
            if statistic != 'A':
                expectation += treeSum(countAdded,before)
            sumBefore = treeSum(deltaSquared,before)
            varianceConditional      [i] = expectation.var(ddof=1)
            meanConditional          [i] = expectation.mean()
            meanConditionalError     [i] = np.std(expectation-countTree)/np.sqrt(countTrees)
            sumDeltaSquaredCumulative[i] = sumBefore.mean()
            errorCumulative          [i] = bootstrap(np.array([expectation,sumBefore]),weightedVarianceMinusMean,rng)
        # Summary statistics.
        meanCount    = countTree.mean()
        varianceCount= countTree.var(ddof=1)
        meanSDS      = sumDeltaSquared.mean()
        errorDiff    = bootstrap(np.array([countTree,sumDeltaSquared]),weightedVarianceMinusMean,rng)
        errorVar     = bootstrap(countTree,weightedVariance,rng)
        errorMean    = bootstrap(sumDelta+sumDrift,weightedMean,rng)
        sigmaI2      = (varianceCount-meanCount)/meanCount**2 if meanCount > 0.0 else np.nan
        print('  {:10.3e} {:9.3f} {:9.3f} {:9.3f} {:9.3f} ± {:7.3f} {:9.3f} {:9.4f} {:8.4f} {:8.1e} {:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
            massThreshold,meanCount,n1Root.mean(),varianceCount,meanSDS,errorDiff,varianceCount/meanCount if meanCount > 0.0 else np.nan,sigmaI2,
            (sumDelta+sumDrift).mean(),np.abs(telescopeResidual).max(),
            sumDeltaSquaredTerminal.mean()/varianceCount,sumDeltaSquaredStructural.mean()/varianceCount,
            sumDeltaSquaredMain    .mean()/varianceCount,sumDeltaSquaredSide      .mean()/varianceCount))
        # Write results.
        group = output.create_group('statistic'+statistic+'/threshold{:d}'.format(indexThreshold))
        group.attrs['statistic'               ] = statistic
        group.attrs['description'             ] = statisticNames[statistic]
        group.attrs['massThreshold'           ] = massThreshold
        group.attrs['countMean'               ] = meanCount
        group.attrs['countVariance'           ] = varianceCount
        group.attrs['countVarianceError'      ] = errorVar
        group.attrs['sumDeltaSquaredMean'     ] = meanSDS
        group.attrs['varianceMinusSumDeltaSquaredError'] = errorDiff
        group.attrs['meanJumpPlusDrift'       ] = (sumDelta+sumDrift).mean()
        group.attrs['meanJumpPlusDriftError'  ] = errorMean
        group.attrs['telescopeResidualMaximum'] = np.abs(telescopeResidual).max()
        group.attrs['tableEvaluationsEmpty'   ] = sum(t.evaluationsEmpty for t in tables)
        cumulativeGroup = group.create_group('cumulative')
        cumulativeGroup.create_dataset('omega'                       ,data=omegaGrid                )
        cumulativeGroup.create_dataset('varianceConditionalMean'     ,data=varianceConditional      )
        cumulativeGroup.create_dataset('sumDeltaSquared'             ,data=sumDeltaSquaredCumulative)
        cumulativeGroup.create_dataset('differenceError'             ,data=errorCumulative          )
        cumulativeGroup.create_dataset('meanConditionalMean'         ,data=meanConditional          )
        cumulativeGroup.create_dataset('meanConditionalMeanError'    ,data=meanConditionalError     )
        print('      Markov check: node type, count, <N-n₁>, error')
        for name, number, mean, error in markovCheck:
            print('        {:20s} {:9d} {:+9.4f} {:8.4f}'.format(name,number,mean,error))
        print('      cumulative check: ω, Var(E[N|ω]), <ΣΔ²(<ω)>, error on difference, <E[N|ω]>-<N>, error')
        for i in range(len(omegaGrid)):
            print('        {:7.4f} {:9.4f} {:9.4f} {:7.4f} {:+8.4f} {:7.4f}'.format(omegaGrid[i],varianceConditional[i],sumDeltaSquaredCumulative[i],errorCumulative[i],meanConditional[i]-meanCount,meanConditionalError[i]))
        treeGroup = group.create_group('trees')
        for name, values in (('count',countTree),('n1Root',n1Root),('sumDelta',sumDelta),('sumDrift',sumDrift),('sumDeltaSquared',sumDeltaSquared),
                             ('sumDeltaSquaredTerminal',sumDeltaSquaredTerminal),('sumDeltaSquaredStructural',sumDeltaSquaredStructural),
                             ('sumDeltaSquaredMain',sumDeltaSquaredMain),('sumDeltaSquaredSide',sumDeltaSquaredSide),
                             ('sumDeltaSquaredBeforeFormation',sumDeltaSquaredBeforeFormation)):
            treeGroup.create_dataset(name,data=values)
        splitGroup = group.create_group('splits')
        s = np.nonzero(activeSplit)[0]
        for name, values in (('tree',treeOfNode[s]),('delta',delta[s]),('omega',omegaSplit[s]),
                             ('redshift',redshift[primaryChild[s]]),('mass',massChildren[s]),('massPrimary',mass[primaryChild[s]]),
                             ('massSecondary',massSecondary[s]),('massRoot',massRoot[treeOfNode[s]]),('n1Parent',n1End[s]),
                             ('isMainBranch',isMain[s].astype(np.int8)),('hierarchyDepth',hierarchyDepth[s]),('countAdded',countAdded[s])):
            splitGroup.create_dataset(name,data=values,compression='gzip')
        tableGroup = group.create_group('n1Table')
        tableGroup.create_dataset('logMassCenters',data=tables[0].logMassCenters)
        tableGroup.create_dataset('omegaCenters'  ,data=tables[0].omegaCenters  )
        for i, table in enumerate(tables):
            tableGroup.create_dataset('mean{:d}'      .format(i),data=table.mean      ,compression='gzip')
            tableGroup.create_dataset('error{:d}'     .format(i),data=table.error     ,compression='gzip')
            tableGroup.create_dataset('countNodes{:d}'.format(i),data=table.countNodes,compression='gzip')
output.close()
print('\nColumns: Var(N) is the sample variance of the count; <ΣΔ²> its split-sum estimate (± bootstrap error on the difference);')
print('σ_I² = (Var-<N>)/<N>²; <ΣΔ+Σd> should be zero; "telescope" is the largest per-tree bookkeeping residual (should be ~0);')
print('the last four columns are the fractions of <ΣΔ²> from splits with non-primary mass M′≥m and M′<m, and from main/side-branch splits.')
