#!/usr/bin/env python3
import os
import sys
import argparse
# The only linear algebra is small matrix products for bootstrapping, for which multithreaded BLAS spends far more time spinning
# than computing.
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import multiprocessing
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
# Trees are processed in chunks of whole forests, so that memory use is set by the chunk size rather than by the number of trees.
# Everything evaluated at a node depends only on its own tree, except for the n₁ table (a sum over bins, accumulated over chunks)
# and a few global ranges, so the trees are read three times: to find the ranges of the table, to accumulate the table, and to
# evaluate the jumps. Chunks can be processed in parallel by several worker processes.
#
# Andrew Benson (20-October-2026)

statisticNames = {
    'A': 'progenitors with mass above threshold at a fixed redshift',
    'B': 'splits producing a non-primary progenitor above threshold, all hierarchy depths',
    'C': 'splits producing a non-primary progenitor above threshold, main branch only'
}

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='mergerTreeVariance.py',description='Decompose the variance of progenitor/subhalo counts in merger trees into contributions from individual splits.')
parser.add_argument('treeFileNames'     ,nargs='+'                                                          ,help='one or more files of exported trees, analyzed together as a single sample (they must share a cosmology)')
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
parser.add_argument('--splitTreeStride' ,action='store',default=1                     ,type=int              ,help='write per-split values only for trees whose index is a multiple of this (per-split output is large for high-resolution trees; all trees are still analyzed)')
parser.add_argument('--processes'       ,action='store',default=1                     ,type=int              ,help='number of worker processes over which chunks are distributed; peak memory scales with this times the chunk size')
parser.add_argument('--chunkNodes'      ,action='store',default=20000000              ,type=int              ,help='approximate number of nodes processed at once (whole forests are always kept together); peak memory scales with this'                                                       )
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
    # center; bins with too few nodes, or a degenerate spread of nodes, fall back to the mean. Nodes are added in any number of
    # batches with accumulate() (every quantity needed is a sum over the nodes in a bin), and the table is completed by finalize().
    def __init__(self,logMassRange,omegaRange,binWidthMass,binWidthOmega,method='mean'):
        self.logMassEdges = logMassRange[0]+binWidthMass *np.arange(int(round((logMassRange[1]-logMassRange[0])/binWidthMass ))+1)
        self.omegaEdges   = omegaRange  [0]+binWidthOmega*np.arange(int(round((omegaRange  [1]-omegaRange  [0])/binWidthOmega))+1)
        self.logMassCenters = 0.5*(self.logMassEdges[1:]+self.logMassEdges[:-1])
        self.omegaCenters   = 0.5*(self.omegaEdges  [1:]+self.omegaEdges  [:-1])
        self.binWidthMass   = binWidthMass
        self.binWidthOmega  = binWidthOmega
        self.method         = method
        self.sums           = {}
        self.evaluationsEmpty = 0

    def accumulate(self,logMass,omega,count):
        def binSum(name,weights=None):
            values = np.histogram2d(omega,logMass,bins=[self.omegaEdges,self.logMassEdges],weights=weights)[0]
            self.sums[name] = self.sums[name]+values if name in self.sums else values
        binSum('countNodes'           )
        binSum('sumCount'  ,count     )
        binSum('sumCount2' ,count**2  )
        if self.method == 'linear':
            # Offsets of each node from its bin center, in units of the bin width.
            i  = np.clip(np.searchsorted(self.omegaEdges  ,omega  ,side='right')-1,0,len(self.omegaCenters  )-1)
            j  = np.clip(np.searchsorted(self.logMassEdges,logMass,side='right')-1,0,len(self.logMassCenters)-1)
            dx = (logMass-self.logMassCenters[j])/self.binWidthMass
            dy = (omega  -self.omegaCenters  [i])/self.binWidthOmega
            for name, weights in (('Sx',dx),('Sy',dy),('Sxx',dx*dx),('Syy',dy*dy),('Sxy',dx*dy),('Sxc',dx*count),('Syc',dy*count)):
                binSum(name,weights)

    def addSums(self,sums):
        # Add bin sums accumulated elsewhere (e.g. by another process).
        for name, values in sums.items():
            self.sums[name] = self.sums[name]+values if name in self.sums else values.copy()

    def finalize(self):
        countNodes = self.sums['countNodes']
        sumCount   = self.sums['sumCount'  ]
        sumCount2  = self.sums['sumCount2' ]
        self.countNodes = countNodes
        with np.errstate(invalid='ignore',divide='ignore'):
            self.mean  = sumCount/countNodes
            self.error = np.sqrt(np.maximum(sumCount2/countNodes-self.mean**2,0.0)/countNodes)
        if self.method == 'linear':
            Sx , Sy , Sxx, Syy, Sxy = (self.sums[name] for name in ('Sx','Sy','Sxx','Syy','Sxy'))
            Sc , Sxc, Syc           = sumCount, self.sums['Sxc'], self.sums['Syc']
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
        self.sums = None
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

class chunkOfTrees:
    # The structure of the trees in a contiguous range of forests: tree topology, primary progenitors, main-branch membership,
    # hierarchy depth, and formation times. Tree indices (treeOfNode, roots) are global: offset by the number of trees in earlier
    # chunks.
    def __init__(self,treeFile,forestStart,forestEnd,treeOffset):
        firstNode     = treeFile['forestIndex/firstNode'    ][forestStart:forestEnd].astype(np.int64)
        numberOfNodes = treeFile['forestIndex/numberOfNodes'][forestStart:forestEnd].astype(np.int64)
        nodeStart     = firstNode[0]
        nodeEnd       = firstNode[-1]+numberOfNodes[-1]
        firstNode    -= nodeStart
        nodeIndex       = treeFile['forestHalos/nodeIndex'      ][nodeStart:nodeEnd].astype(np.int64)
        descendantIndex = treeFile['forestHalos/descendantIndex'][nodeStart:nodeEnd].astype(np.int64)
        self.mass       = treeFile['forestHalos/nodeMass'       ][nodeStart:nodeEnd]
        self.redshift   = treeFile['forestHalos/redshift'       ][nodeStart:nodeEnd]
        mass            = self.mass
        countNodesTotal = len(nodeIndex)
        self.countNodes = countNodesTotal
        forest          = np.repeat(np.arange(len(firstNode),dtype=np.int64),numberOfNodes)
        # Compute ω for each node.
        self.omega   = args.deltaCritical0/linearGrowthFactor(self.redshift,OmegaMatter,OmegaLambda)
        self.logMass = np.log10(mass)
        omega        = self.omega
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
        del key, order, position, forest, nodeIndex, descendantIndex
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
        del jump, active
        root = np.arange(countNodesTotal)
        root[hasParent] = parent[hasParent]
        while True:
            rootNext = root[root]
            if np.all(rootNext == root):
                break
            root = rootNext
        del rootNext
        roots           = np.nonzero(~hasParent)[0]
        self.countTrees = len(roots)
        self.treeOffset = treeOffset
        treeOfNode      = np.searchsorted(roots,root)
        del root
        # Group nodes by depth.
        self.depthOrder      = np.argsort(depth,kind='stable')
        self.depthMaximum    = depth.max()
        self.depthBoundaries = np.searchsorted(depth[self.depthOrder],np.arange(self.depthMaximum+2))
        del depth
        # Identify the primary (most massive) progenitor of each node.
        orderChildren    = children[np.lexsort((mass[children],parent[children]))]
        lastOfGroup      = np.append(parent[orderChildren][1:] != parent[orderChildren][:-1],True)
        primaryChild     = np.full(countNodesTotal,-1,dtype=np.int64)
        primaryChild[parent[orderChildren[lastOfGroup]]] = orderChildren[lastOfGroup]
        del orderChildren, lastOfGroup
        isPrimary        = np.zeros(countNodesTotal,dtype=bool)
        isPrimary[primaryChild[primaryChild >= 0]] = True
        countChildren    = np.bincount(parent[children],minlength=countNodesTotal)
        massChildren     = np.bincount(parent[children],weights=mass[children],minlength=countNodesTotal)
        if args.segmentEndMass == 'node':
            massChildren = mass.copy()
        secondary        = children[~isPrimary[children]]
        massSecondary    = np.bincount(parent[secondary],weights=mass[secondary],minlength=countNodesTotal)
        omegaSplit       = np.where(primaryChild >= 0,omega[np.maximum(primaryChild,0)],np.nan)
        self.parent, self.hasParent, self.children, self.roots, self.primaryChild, self.isPrimary = parent, hasParent, children, roots, primaryChild, isPrimary
        self.countChildren, self.massChildren, self.secondary, self.massSecondary, self.omegaSplit = countChildren, massChildren, secondary, massSecondary, omegaSplit
        # Main-branch membership and hierarchy depth (number of non-primary links between a node and its root).
        isMain         = ~hasParent.copy()
        hierarchyDepth = np.zeros(countNodesTotal,dtype=np.int64)
        for d in range(1,self.depthMaximum+1):
            nodes = self.nodesAtDepth(d)
            isMain        [nodes] = isMain        [parent[nodes]] & isPrimary[nodes]
            hierarchyDepth[nodes] = hierarchyDepth[parent[nodes]]+(~isPrimary[nodes])
        self.isMain, self.hierarchyDepth = isMain, hierarchyDepth
        # Splits are nodes with two or more progenitors.
        self.countSplits      = np.count_nonzero(countChildren >= 2)
        self.countMultiple    = np.count_nonzero(countChildren >  2)
        # Check that the progenitors of each split share a single ω.
        omegaSpread = np.zeros(countNodesTotal)
        np.maximum.at(omegaSpread,parent[children],np.abs(omega[children]-omegaSplit[parent[children]]))
        self.omegaSpreadMaximum = omegaSpread.max()
        del omegaSpread
        # Formation redshift (main branch first drops below half of the root mass) of each tree.
        massRoot          = mass[roots]
        redshiftFormation = np.full(self.countTrees,np.nan)
        omegaFormation    = np.full(self.countTrees,np.nan)
        mainNodes         = np.nonzero(isMain)[0]
        mainNodes         = mainNodes[np.lexsort((omega[mainNodes],treeOfNode[mainNodes]))]
        belowHalf         = mass[mainNodes] < 0.5*massRoot[treeOfNode[mainNodes]]
        firstBelow        = mainNodes[belowHalf][np.diff(treeOfNode[mainNodes[belowHalf]],prepend=-1) != 0]
        before            = parent[firstBelow]
        tree              = treeOfNode[firstBelow]
        factor            = (np.log(0.5*massRoot[tree])-np.log(mass[before]))/(np.log(mass[firstBelow])-np.log(mass[before]))
        redshiftFormation[tree] = self.redshift[before]+factor*(self.redshift[firstBelow]-self.redshift[before])
        omegaFormation   [tree] = omega        [before]+factor*(omega        [firstBelow]-omega        [before])
        self.massRoot, self.redshiftFormation, self.omegaFormation = massRoot, redshiftFormation, omegaFormation
        # Local tree index of each node (for per-tree sums within the chunk), and global tree index (for fold assignment and
        # output).
        self.treeLocal  = treeOfNode
        self.treeOfNode = treeOfNode+treeOffset

    def nodesAtDepth(self,d):
        return self.depthOrder[self.depthBoundaries[d]:self.depthBoundaries[d+1]]

    def countSubtree(self,own,primaryOnly):
        # Sum a per-node count over each node's subtree (or, if primaryOnly, over its own primary-progenitor chain).
        total = own.astype(np.float64).copy()
        for d in range(self.depthMaximum,0,-1):
            nodes = self.nodesAtDepth(d)
            if primaryOnly:
                nodes = nodes[self.isPrimary[nodes]]
            np.add.at(total,self.parent[nodes],total[nodes])
        return total

    def treeSum(self,values,select=None):
        if select is None:
            return np.bincount(self.treeLocal,weights=values,minlength=self.countTrees)
        return np.bincount(self.treeLocal[select],weights=values[select],minlength=self.countTrees)

    def statisticCounts(self,statistic,massThreshold):
        # Per-node counts for a statistic, the set of nodes which are valid starting points for the n₁ table, and the active splits.
        c = {}
        # Count of objects added by each split (non-primary progenitors above threshold).
        c['countAdded'] = np.bincount(self.parent[self.secondary],weights=(self.mass[self.secondary] >= massThreshold).astype(np.float64),minlength=self.countNodes)
        if statistic == 'A':
            # A branch segment from node p crosses ω₁ if ω_p < ω₁ ≤ ω_split. Its mass at ω₁ is found by interpolating log mass
            # linearly in ω between the node and the summed progenitor mass at the split.
            crosses       = (self.omega < omega1) & (self.countChildren > 0) & (self.omegaSplit >= omega1)
            logMassTarget = np.full(self.countNodes,-np.inf)
            p = np.nonzero(crosses)[0]
            logMassTarget[p] = self.logMass[p]+(np.log10(self.massChildren[p])-self.logMass[p])*(omega1-self.omega[p])/(self.omegaSplit[p]-self.omega[p])
            c['crosses'    ] = crosses
            c['own'        ] = (logMassTarget >= np.log10(massThreshold)).astype(np.float64)
            c['count'      ] = self.countSubtree(c['own'],primaryOnly=False)
            c['tableNodes' ] = self.omega < omega1
            c['activeSplit'] = (self.countChildren >= 2) & (self.omegaSplit < omega1)
        else:
            c['own'        ] = c['countAdded']
            c['count'      ] = self.countSubtree(c['own'],primaryOnly=(statistic == 'C'))
            c['tableNodes' ] = np.ones(self.countNodes,dtype=bool)
            c['activeSplit'] = self.countChildren >= 2
            if statistic == 'C':
                c['activeSplit'] &= self.isMain
        return c

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

treeFilesWorker = {}
def treeFileOpen(indexFile):
    # Return this process' handle to a tree file (HDF5 file handles must not be shared between processes).
    if indexFile not in treeFilesWorker:
        treeFilesWorker[indexFile] = h5py.File(args.treeFileNames[indexFile],'r')
    return treeFilesWorker[indexFile]

def mapChunks(function,items):
    # Apply a function to each item, yielding results in order, using a pool of worker processes if requested. Workers are
    # forked when iteration begins, so they inherit everything (including the finalized n₁ tables) defined before then.
    if args.processes <= 1:
        yield from map(function,items)
        return
    with multiprocessing.get_context('fork').Pool(args.processes) as pool:
        yield from pool.imap(function,items)

# Open the exported trees, and divide the forests of each file into chunks of roughly the requested number of nodes. A chunk is
# (index of file, first forest, last forest+1); trees are numbered globally in the order of files, then of chunks.
nodesBefore = []
chunks      = []
countForests, countNodesAll = 0, 0
for indexFile, treeFileName in enumerate(args.treeFileNames):
    print('Reading trees from "'+treeFileName+'"')
    treeFile      = h5py.File(treeFileName,'r')
    cosmology     = (treeFile['cosmology'].attrs['OmegaMatter'],treeFile['cosmology'].attrs['OmegaLambda'])
    if indexFile == 0:
        OmegaMatter, OmegaLambda = cosmology
    elif cosmology != (OmegaMatter,OmegaLambda):
        sys.exit('mergerTreeVariance.py: tree files do not share a cosmology')
    numberOfNodes = treeFile['forestIndex/numberOfNodes'][:].astype(np.int64)
    nodesBefore.append(np.concatenate(([0],np.cumsum(numberOfNodes))))
    # Chunks are read as contiguous ranges of nodes, so forests must be stored contiguously and in order.
    if np.any(treeFile['forestIndex/firstNode'][:].astype(np.int64) != nodesBefore[-1][:-1]):
        sys.exit('mergerTreeVariance.py: forests are not stored contiguously in order')
    chunkBoundaries = np.unique(np.concatenate((np.searchsorted(nodesBefore[-1],np.arange(0,nodesBefore[-1][-1],args.chunkNodes),side='right')-1,[len(numberOfNodes)])))
    chunks.extend((indexFile,forestStart,forestEnd) for forestStart, forestEnd in zip(chunkBoundaries[:-1],chunkBoundaries[1:]))
    countForests  += len(numberOfNodes)
    countNodesAll += nodesBefore[-1][-1]
    treeFile.close()
del numberOfNodes
omega1        = args.deltaCritical0/linearGrowthFactor(args.redshiftTarget,OmegaMatter,OmegaLambda)
print('  {:d} forests, {:d} nodes, in {:d} file(s) and {:d} chunk(s)'.format(countForests,countNodesAll,len(args.treeFileNames),len(chunks)))

# Pass 1: global ranges of mass and redshift (which set the n₁ table ranges), and the number of trees in each chunk (which sets
# the global index of each tree, used to assign trees to folds).
def passRanges(chunk):
    indexFile, forestStart, forestEnd = chunk
    nodeStart, nodeEnd = nodesBefore[indexFile][forestStart], nodesBefore[indexFile][forestEnd]
    treeFile_     = treeFileOpen(indexFile)
    massChunk     = treeFile_['forestHalos/nodeMass'       ][nodeStart:nodeEnd]
    redshiftChunk = treeFile_['forestHalos/redshift'       ][nodeStart:nodeEnd]
    countRoots    = np.count_nonzero(treeFile_['forestHalos/descendantIndex'][nodeStart:nodeEnd] < 0)
    return massChunk.min(), massChunk.max(), redshiftChunk.min(), redshiftChunk.max(), countRoots
ranges = list(mapChunks(passRanges,chunks))
massMinimum     = min(r[0] for r in ranges)
massMaximum     = max(r[1] for r in ranges)
redshiftMinimum = min(r[2] for r in ranges)
redshiftMaximum = max(r[3] for r in ranges)
treeCounts      = np.array([r[4] for r in ranges],dtype=np.int64)
treeOffsets     = np.concatenate(([0],np.cumsum(treeCounts)))[:-1]
countTrees      = int(treeCounts.sum())
del ranges
# ω increases monotonically with redshift.
omegaMinimum, omegaMaximum = args.deltaCritical0/linearGrowthFactor(np.array([redshiftMinimum,redshiftMaximum]),OmegaMatter,OmegaLambda)
logMassMinimum, logMassMaximum = np.log10(massMinimum), np.log10(massMaximum)
massResolution = massMinimum

rng = np.random.default_rng(args.seed)
countFolds = 3 if args.crossProduct else (2 if args.crossFit else 1)
massThresholds = [float(m) for m in args.massThresholds.split(',')]
for statistic in args.statistics:
    if statistic not in statisticNames:
        sys.exit('mergerTreeVariance.py: unknown statistic "'+statistic+'"')
# Ranges of the n₁ table. Bins are aligned so that the largest node mass and the smallest ω (i.e. the roots of trees of fixed
# mass) lie at bin centers, avoiding the bias of clamping at the table edge. At least two bins are used in each direction, as
# needed for interpolation.
countBinsMass  = max(int(np.ceil((logMassMaximum-logMassMinimum)/args.binWidthMass +0.5)),2)
countBinsOmega = max(int(np.ceil((omegaMaximum  -omegaMinimum  )/args.binWidthOmega+0.5)),2)
logMassRange   = (logMassMaximum+0.5*args.binWidthMass -countBinsMass *args.binWidthMass ,logMassMaximum+0.5*args.binWidthMass )
omegaRange     = (omegaMinimum  -0.5*args.binWidthOmega,omegaMinimum  -0.5*args.binWidthOmega+countBinsOmega*args.binWidthOmega)
cases          = [(statistic,indexThreshold,massThreshold) for statistic in args.statistics for indexThreshold, massThreshold in enumerate(massThresholds)]
def newTables():
    return {(statistic,indexThreshold): [nodeMeanTable(logMassRange,omegaRange,args.binWidthMass,args.binWidthOmega,method=args.tableMethod) for h in range(countFolds)] for statistic, indexThreshold, massThreshold in cases}

# Pass 2: build the tree structure of each chunk, and accumulate the n₁ tables, per-tree quantities, and the ω of active splits
# (needed for the range of the cumulative check). Each worker is given a contiguous group of chunks, and returns the bin sums of
# its tables (so that only one set of sums per worker is passed back).
def passTables(group):
    tablesGroup = newTables()
    summary     = {'perTree': {name: [] for name in ('massRoot','redshiftFormation','omegaFormation')},
                   'omegaSplitsActive': {statistic: [] for statistic in args.statistics},
                   'countSplits': 0, 'countMultiple': 0, 'depthMaximum': 0, 'omegaSpreadMaximum': 0.0}
    for indexChunk in group:
        indexFile, forestStart, forestEnd = chunks[indexChunk]
        chunk = chunkOfTrees(treeFileOpen(indexFile),forestStart,forestEnd,treeOffsets[indexChunk])
        if chunk.countTrees != treeCounts[indexChunk]:
            sys.exit('mergerTreeVariance.py: inconsistent count of trees in chunk')
        summary['countSplits'       ] += chunk.countSplits
        summary['countMultiple'     ] += chunk.countMultiple
        summary['depthMaximum'      ]  = max(summary['depthMaximum'      ],chunk.depthMaximum      )
        summary['omegaSpreadMaximum']  = max(summary['omegaSpreadMaximum'],chunk.omegaSpreadMaximum)
        for name in summary['perTree']:
            summary['perTree'][name].append(getattr(chunk,name))
        for statistic, indexThreshold, massThreshold in cases:
            c      = chunk.statisticCounts(statistic,massThreshold)
            select = c['tableNodes']
            for h in range(countFolds):
                selectFold = select if countFolds == 1 else select & (chunk.treeOfNode % countFolds == h)
                tablesGroup[(statistic,indexThreshold)][h].accumulate(chunk.logMass[selectFold],chunk.omega[selectFold],c['count'][selectFold])
            if indexThreshold == 0:
                summary['omegaSplitsActive'][statistic].append(chunk.omegaSplit[c['activeSplit']])
        del chunk, c
    summary['sums'] = {key: [table.sums for table in tablesGroup[key]] for key in tablesGroup}
    return summary
countGroups = min(max(args.processes,1),len(chunks))
groups      = [list(g) for g in np.array_split(np.arange(len(chunks)),countGroups)]
tables      = newTables()
perTree     = {name: [] for name in ('massRoot','redshiftFormation','omegaFormation')}
omegaSplitsActive = {statistic: [] for statistic in args.statistics}
countSplits, countMultiple, depthMaximum, omegaSpreadMaximum = 0, 0, 0, 0.0
for summary in mapChunks(passTables,groups):
    countSplits       += summary['countSplits'  ]
    countMultiple     += summary['countMultiple']
    depthMaximum       = max(depthMaximum      ,summary['depthMaximum'      ])
    omegaSpreadMaximum = max(omegaSpreadMaximum,summary['omegaSpreadMaximum'])
    for name in perTree:
        perTree[name].extend(summary['perTree'][name])
    for statistic in args.statistics:
        omegaSplitsActive[statistic].extend(summary['omegaSplitsActive'][statistic])
    for key in tables:
        for table, sums in zip(tables[key],summary['sums'][key]):
            table.addSums(sums)
    del summary
for name in perTree:
    perTree[name] = np.concatenate(perTree[name])
for key in tables:
    for table in tables[key]:
        table.finalize()
massRoot, redshiftFormation, omegaFormation = perTree['massRoot'], perTree['redshiftFormation'], perTree['omegaFormation']
print('  {:d} trees, {:d} nodes, {:d} splits, maximum depth {:d}'.format(countTrees,countNodesAll,countSplits,depthMaximum))
if countMultiple > 0:
    print('  note: {:d} nodes have more than two progenitors'.format(countMultiple))
print('  maximum ω difference between progenitors of a split: {:.3e}'.format(omegaSpreadMaximum))
print('  minimum node mass {:.3e}; mean root mass {:.3e}'.format(massResolution,massRoot.mean()))
print('  statistic A target: z₁={:.3f}, ω₁={:.4f} (ω₀={:.4f})'.format(args.redshiftTarget,omega1,args.deltaCritical0))
# Grid of ω for the cumulative check.
omegaGrids = {}
for statistic in args.statistics:
    omegaGrids[statistic] = np.linspace(omegaMinimum,omega1 if statistic == 'A' else np.quantile(np.concatenate(omegaSplitsActive[statistic]),0.999),args.cumulativeCount)
del omegaSplitsActive

output = h5py.File(args.outputFileName,'w')
output.attrs['treeFileName'      ] = ','.join(args.treeFileNames)
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

# Pass 3: evaluate the jumps, drifts, and their per-tree sums, chunk by chunk. Each chunk returns its per-tree sums, Markov
# check sums, cumulative-check expectations, per-split values, and the number of table evaluations touching empty bins; these
# are gathered in chunk order.
treeQuantities = ('count','n1Root','sumDelta','sumDrift','sumDeltaSquared','sumDeltaSquaredTerminal','sumDeltaSquaredStructural',
                  'sumDeltaSquaredMain','sumDeltaSquaredSide','sumDeltaSquaredBeforeFormation','telescopeResidual')
nodeTypeNames  = ('root','primary at split','secondary at split','accretion')
def passJumps(indexChunk):
    indexFile, forestStart, forestEnd = chunks[indexChunk]
    chunk   = chunkOfTrees(treeFileOpen(indexFile),forestStart,forestEnd,treeOffsets[indexChunk])
    returns = {}
    for statistic, indexThreshold, massThreshold in cases:
        tablesCase  = tables[(statistic,indexThreshold)]
        emptyBefore = sum(table.evaluationsEmpty for table in tablesCase)
        result      = {'trees': {}, 'markov': {}, 'expectation': np.empty((args.cumulativeCount,chunk.countTrees)), 'sumBefore': np.empty((args.cumulativeCount,chunk.countTrees))}
        c           = chunk.statisticCounts(statistic,massThreshold)
        countAdded, own, count, tableNodes, activeSplit = c['countAdded'], c['own'], c['count'], c['tableNodes'], c['activeSplit']
        countTree   = count[chunk.roots]
        def n1(logMass_,omega_,tree_,variant=0):
            # Evaluate n₁ using the appropriate table: with folds, a table built from another fold of the trees (variant 0 or 1
            # selects which of the other folds, for the cross-product estimator).
            if countFolds == 1:
                return tablesCase[0](logMass_,omega_)
            result_ = np.empty(len(logMass_))
            for h in range(countFolds):
                s_ = (tree_ % countFolds) == h
                result_[s_] = tablesCase[(h+1+variant) % countFolds](logMass_[s_],omega_[s_])
            return result_
        # n₁ at each node, and at the end state of each segment (the summed progenitor mass at the split, or the branch end).
        segments = np.nonzero(chunk.countChildren > 0)[0]
        def n1NodeEnd(variant):
            n1Node_ = n1(chunk.logMass,chunk.omega,chunk.treeOfNode,variant)
            n1End_  = np.zeros(chunk.countNodes)
            n1End_[segments] = n1(np.log10(chunk.massChildren[segments]),chunk.omegaSplit[segments],chunk.treeOfNode[segments],variant)
            return n1Node_, n1End_
        n1Node, n1End = n1NodeEnd(0)
        # Markov check: by the Markov property the mean of (count - n₁) must vanish for every kind of node (roots, the primary and
        # non-primary progenitors created at a split, and nodes created by accretion alone). Accumulate the count, sum, and sum of
        # squares of the residual for each node type.
        parentChildren = chunk.countChildren[np.maximum(chunk.parent,0)]
        nodeTypes = (~chunk.hasParent                                                          ,
                      chunk.hasParent &  chunk.isPrimary & (parentChildren >= 2)               ,
                      chunk.hasParent & ~chunk.isPrimary                                       ,
                      chunk.hasParent & (parentChildren == 1)                                  )
        for name, select in zip(nodeTypeNames,nodeTypes):
            select   = select & tableNodes & (n1Node > 0.0)
            residual = count[select]-n1Node[select]
            result['markov'][name] = np.array([len(residual),residual.sum(),(residual**2).sum()])
        del parentChildren, nodeTypes
        # Jumps at splits.
        s = np.nonzero(activeSplit)[0]
        def jumps(n1Node_,n1End_):
            delta_ = np.zeros(chunk.countNodes)
            if statistic == 'C':
                # Role-dependent n₁: only the primary progenitor carries the main branch forward.
                n1Children = np.zeros(chunk.countNodes)
                n1Children[segments] = n1Node_[chunk.primaryChild[segments]]
            else:
                n1Children = np.bincount(chunk.parent[chunk.children],weights=n1Node_[chunk.children],minlength=chunk.countNodes)
            delta_[s] = n1Children[s]+(countAdded[s] if statistic != 'A' else 0.0)-n1End_[s]
            return delta_
        delta = jumps(n1Node,n1End)
        # Deterministic drift of n₁ along each segment.
        drift = np.zeros(chunk.countNodes)
        if statistic == 'A':
            crosses = c['crosses']
            live = chunk.omega < omega1
            # Segments ending at a split before ω₁, or terminating (leaf) before ω₁.
            ends = live & ~crosses
            drift[ends] = np.where(chunk.countChildren[ends] > 0,n1End[ends],0.0)-n1Node[ends]
            # Segments crossing ω₁ end with the count of the surviving branch.
            drift[crosses] = own[crosses]-n1Node[crosses]
        elif statistic == 'B':
            drift = np.where(chunk.countChildren > 0,n1End,0.0)-n1Node
        else:
            isMain = chunk.isMain
            drift[isMain] = np.where(chunk.countChildren[isMain] > 0,n1End[isMain],0.0)-n1Node[isMain]
        # Per-tree sums.
        terminal     = chunk.massSecondary >= massThreshold
        if args.crossProduct:
            # Product of jumps evaluated with tables from the two other folds: the table noise in the two is independent, so its
            # variance does not contribute to the expectation.
            deltaSquared = delta*jumps(*n1NodeEnd(1))
        else:
            deltaSquared = delta**2
        t = result['trees']
        t['count'                         ] = countTree
        t['n1Root'                        ] = n1Node[chunk.roots]
        t['sumDelta'                      ] = chunk.treeSum(delta                       )
        t['sumDeltaSquared'               ] = chunk.treeSum(deltaSquared                )
        t['sumDeltaSquaredTerminal'       ] = chunk.treeSum(deltaSquared, terminal      )
        t['sumDeltaSquaredStructural'     ] = chunk.treeSum(deltaSquared,~terminal      )
        t['sumDeltaSquaredMain'           ] = chunk.treeSum(deltaSquared, chunk.isMain  )
        t['sumDeltaSquaredSide'           ] = chunk.treeSum(deltaSquared,~chunk.isMain  )
        t['sumDrift'                      ] = chunk.treeSum(drift                       )
        beforeFormation = chunk.omegaSplit < chunk.omegaFormation[chunk.treeLocal]
        t['sumDeltaSquaredBeforeFormation'] = chunk.treeSum(deltaSquared,beforeFormation)
        t['telescopeResidual'             ] = countTree-t['n1Root']-t['sumDelta']-t['sumDrift']
        # Cumulative variance: at each ω on a grid, Var(E[N|tree built to ω]) should equal <Σ_{ω_split<ω} Δ²>. E[N|tree built to
        # ω] is the count already realized (for statistics B, C) plus n₁ summed over the branches alive at ω.
        for i, omegaG in enumerate(omegaGrids[statistic]):
            alive = (chunk.countChildren > 0) & (chunk.omega <= omegaG) & (chunk.omegaSplit > omegaG)
            if statistic == 'C':
                alive &= chunk.isMain
            p = np.nonzero(alive)[0]
            logMassG = chunk.logMass[p]+(np.log10(chunk.massChildren[p])-chunk.logMass[p])*(omegaG-chunk.omega[p])/(chunk.omegaSplit[p]-chunk.omega[p])
            expectation = np.bincount(chunk.treeLocal[p],weights=n1(logMassG,np.full(len(p),omegaG),chunk.treeOfNode[p]),minlength=chunk.countTrees)
            before      = activeSplit & (chunk.omegaSplit <= omegaG)
            if statistic != 'A':
                expectation += chunk.treeSum(countAdded,before)
            result['expectation'][i] = expectation
            result['sumBefore'  ][i] = chunk.treeSum(deltaSquared,before)
        # Per-split values, for the trees selected for output.
        s = s[chunk.treeOfNode[s] % args.splitTreeStride == 0]
        result['splits'] = (('tree',chunk.treeOfNode[s]),('delta',delta[s]),('deltaSquared',deltaSquared[s]),('omega',chunk.omegaSplit[s]),
                            ('redshift',chunk.redshift[chunk.primaryChild[s]]),('mass',chunk.massChildren[s]),('massPrimary',chunk.mass[chunk.primaryChild[s]]),
                            ('massSecondary',chunk.massSecondary[s]),('massRoot',chunk.massRoot[chunk.treeLocal[s]]),('n1Parent',n1End[s]),
                            ('isMainBranch',chunk.isMain[s].astype(np.int8)),('hierarchyDepth',chunk.hierarchyDepth[s]),('countAdded',countAdded[s]))
        result['evaluationsEmpty'] = sum(table.evaluationsEmpty for table in tablesCase)-emptyBefore
        returns[(statistic,indexThreshold)] = result
        del c, n1Node, n1End, delta, drift, deltaSquared, terminal, beforeFormation
    return returns
results = {}
for statistic, indexThreshold, massThreshold in cases:
    results[(statistic,indexThreshold)] = {
        'trees'           : {name: np.zeros(countTrees) for name in treeQuantities},
        'markov'          : {name: np.zeros(3) for name in nodeTypeNames},
        'expectation'     : np.zeros((args.cumulativeCount,countTrees)),
        'sumBefore'       : np.zeros((args.cumulativeCount,countTrees)),
        'splits'          : None,
        'evaluationsEmpty': 0
    }
for indexChunk, returns in enumerate(mapChunks(passJumps,range(len(chunks)))):
    trees = slice(treeOffsets[indexChunk],treeOffsets[indexChunk]+treeCounts[indexChunk])
    for statistic, indexThreshold, massThreshold in cases:
        result = results[(statistic,indexThreshold)]
        chunkResult = returns[(statistic,indexThreshold)]
        for name in treeQuantities:
            result['trees'][name][trees] = chunkResult['trees'][name]
        for name in nodeTypeNames:
            result['markov'][name] += chunkResult['markov'][name]
        result['expectation'][:,trees] = chunkResult['expectation']
        result['sumBefore'  ][:,trees] = chunkResult['sumBefore'  ]
        result['evaluationsEmpty']    += chunkResult['evaluationsEmpty']
        # Per-split output, appended chunk by chunk.
        if result['splits'] is None:
            splitGroup = output.create_group('statistic'+statistic+'/threshold{:d}/splits'.format(indexThreshold))
            result['splits'] = {name: splitGroup.create_dataset(name,data=values,maxshape=(None,),chunks=True,compression='gzip',compression_opts=1) for name, values in chunkResult['splits']}
            # Record how many trees contribute per-split values, so that sums over splits can be normalized per tree.
            splitGroup.attrs['treeStride'] = args.splitTreeStride
            splitGroup.attrs['treeCount' ] = len(range(0,countTrees,args.splitTreeStride))
        else:
            for name, values in chunkResult['splits']:
                dataset = result['splits'][name]
                countBefore = dataset.shape[0]
                dataset.resize((countBefore+len(values),))
                dataset[countBefore:] = values
    del returns

# Summarize each statistic and threshold.
for statistic in args.statistics:
    print('\nStatistic '+statistic+': '+statisticNames[statistic])
    print('  {:>10s} {:>9s} {:>9s} {:>9s} {:>19s} {:>9s} {:>9s} {:>8s} {:>8s} {:>9s} {:>9s} {:>9s} {:>9s}'.format('m','<N>','n₁(root)','Var(N)','<ΣΔ²>','Var/<N>','σ_I²','<ΣΔ+Σd>','telescope','M′≥m/Var','M′<m/Var','main/Var','side/Var'))
    for indexThreshold, massThreshold in enumerate(massThresholds):
        result            = results[(statistic,indexThreshold)]
        t                 = result['trees']
        countTree         = t['count'          ]
        sumDelta          = t['sumDelta'       ]
        sumDrift          = t['sumDrift'       ]
        sumDeltaSquared   = t['sumDeltaSquared']
        n1Root            = t['n1Root'         ]
        telescopeResidual = t['telescopeResidual']
        meanCount         = countTree.mean()
        omegaGrid         = omegaGrids[statistic]
        varianceConditional       = np.empty(args.cumulativeCount)
        meanConditional           = np.empty(args.cumulativeCount)
        meanConditionalError      = np.empty(args.cumulativeCount)
        sumDeltaSquaredCumulative = np.empty(args.cumulativeCount)
        errorCumulative           = np.empty(args.cumulativeCount)
        for i in range(args.cumulativeCount):
            expectation = result['expectation'][i]
            sumBefore   = result['sumBefore'  ][i]
            varianceConditional      [i] = expectation.var(ddof=1)
            meanConditional          [i] = expectation.mean()
            meanConditionalError     [i] = np.std(expectation-countTree)/np.sqrt(countTrees)
            sumDeltaSquaredCumulative[i] = sumBefore.mean()
            errorCumulative          [i] = bootstrap(np.array([expectation,sumBefore]),weightedVarianceMinusMean,rng)
        # Summary statistics.
        varianceCount= countTree.var(ddof=1)
        meanSDS      = sumDeltaSquared.mean()
        errorDiff    = bootstrap(np.array([countTree,sumDeltaSquared]),weightedVarianceMinusMean,rng)
        errorVar     = bootstrap(countTree,weightedVariance,rng)
        errorMean    = bootstrap(sumDelta+sumDrift,weightedMean,rng)
        sigmaI2      = (varianceCount-meanCount)/meanCount**2 if meanCount > 0.0 else np.nan
        print('  {:10.3e} {:9.3f} {:9.3f} {:9.3f} {:9.3f} ± {:7.3f} {:9.3f} {:9.4f} {:8.4f} {:8.1e} {:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
            massThreshold,meanCount,n1Root.mean(),varianceCount,meanSDS,errorDiff,varianceCount/meanCount if meanCount > 0.0 else np.nan,sigmaI2,
            (sumDelta+sumDrift).mean(),np.abs(telescopeResidual).max(),
            t['sumDeltaSquaredTerminal'].mean()/varianceCount,t['sumDeltaSquaredStructural'].mean()/varianceCount,
            t['sumDeltaSquaredMain'    ].mean()/varianceCount,t['sumDeltaSquaredSide'      ].mean()/varianceCount))
        # Write results.
        group = output.require_group('statistic'+statistic+'/threshold{:d}'.format(indexThreshold))
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
        group.attrs['tableEvaluationsEmpty'   ] = result['evaluationsEmpty']
        cumulativeGroup = group.create_group('cumulative')
        cumulativeGroup.create_dataset('omega'                       ,data=omegaGrid                )
        cumulativeGroup.create_dataset('varianceConditionalMean'     ,data=varianceConditional      )
        cumulativeGroup.create_dataset('sumDeltaSquared'             ,data=sumDeltaSquaredCumulative)
        cumulativeGroup.create_dataset('differenceError'             ,data=errorCumulative          )
        cumulativeGroup.create_dataset('meanConditionalMean'         ,data=meanConditional          )
        cumulativeGroup.create_dataset('meanConditionalMeanError'    ,data=meanConditionalError     )
        print('      Markov check: node type, count, <N-n₁>, error')
        for name in nodeTypeNames:
            number, total, totalSquared = result['markov'][name]
            if number > 1:
                mean  = total/number
                error = np.sqrt(max(totalSquared/number-mean**2,0.0))/np.sqrt(number)
            else:
                mean, error = np.nan, np.nan
            print('        {:20s} {:9d} {:+9.4f} {:8.4f}'.format(name,int(number),mean,error))
        print('      cumulative check: ω, Var(E[N|ω]), <ΣΔ²(<ω)>, error on difference, <E[N|ω]>-<N>, error')
        for i in range(len(omegaGrid)):
            print('        {:7.4f} {:9.4f} {:9.4f} {:7.4f} {:+8.4f} {:7.4f}'.format(omegaGrid[i],varianceConditional[i],sumDeltaSquaredCumulative[i],errorCumulative[i],meanConditional[i]-meanCount,meanConditionalError[i]))
        treeGroup = group.create_group('trees')
        for name in treeQuantities:
            if name != 'telescopeResidual':
                treeGroup.create_dataset(name,data=t[name])
        tableGroup = group.create_group('n1Table')
        tablesCase = tables[(statistic,indexThreshold)]
        tableGroup.create_dataset('logMassCenters',data=tablesCase[0].logMassCenters)
        tableGroup.create_dataset('omegaCenters'  ,data=tablesCase[0].omegaCenters  )
        for i, table in enumerate(tablesCase):
            tableGroup.create_dataset('mean{:d}'      .format(i),data=table.mean      ,compression='gzip')
            tableGroup.create_dataset('error{:d}'     .format(i),data=table.error     ,compression='gzip')
            tableGroup.create_dataset('countNodes{:d}'.format(i),data=table.countNodes,compression='gzip')
output.close()
print('\nColumns: Var(N) is the sample variance of the count; <ΣΔ²> its split-sum estimate (± bootstrap error on the difference);')
print('σ_I² = (Var-<N>)/<N>²; <ΣΔ+Σd> should be zero; "telescope" is the largest per-tree bookkeeping residual (should be ~0);')
print('the last four columns are the fractions of <ΣΔ²> from splits with non-primary mass M′≥m and M′<m, and from main/side-branch splits.')
