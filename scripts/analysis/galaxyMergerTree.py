#!/usr/bin/env python3
import sys
import numpy as np
import h5py
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology

# Make a simple diagram showing a galaxy merger tree.
# Andrew Benson (18-October-2024)

# Convert command line arguments to integers.
def restricted_integer(x):
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a integer literal" % (x,))
    return x

# Parse command line arguments.
parser = argparse.ArgumentParser(prog='galaxyMergerTree.py',description='Make a plot of a galaxy merger tree.')
parser.add_argument('filename')
parser.add_argument('--mergerTreeIndex', action='store',default='1'                   ,type=restricted_integer,help='the index of the merger tree to plot')
parser.add_argument('--outputFileName' , action='store',default='galaxyMergerTree.pdf'                        ,help='the name of the output file'         )
args = parser.parse_args()

# Open the model output and extract merger tree datasets.
model         = h5py.File(args.filename,"r")
output        = model ['Outputs/Output1']
nodes         = output['nodeData']
treePropertyNames = [ 'mergerTreeIndex', 'mergerTreeStartIndex', 'mergerTreeCount']
treeProperties = {}
for property in treePropertyNames:
    treeProperties[property] = output[property][:]
treeSelect = treeProperties['mergerTreeIndex'][:] == args.mergerTreeIndex
if np.count_nonzero(treeSelect) != 1:
    print("unable to locate merger tree")
    sys.exit(1)
nodesStart =            treeProperties['mergerTreeStartIndex'][treeSelect][0]
nodesEnd   = nodesStart+treeProperties['mergerTreeCount'     ][treeSelect][0]
propertyNames = [ 'galaxyMergerTreeNodeIndex'          , 'galaxyMergerTreeCount'    , 'galaxyMergerTreeTotalStarFormationRate',
                  'galaxyMergerTreeMassStellarTotal'   , 'galaxyMergerTreeTime'     , 'galaxyMergerTreeMergeCentralIndex'     ,
                  'galaxyMergerTreeMergeSatelliteIndex', 'galaxyMergerTreeMergeTime', 'nodeIsIsolated'                        ,
                  'mergerTreeIndex'                    , 'nodeIndex'                                                            ]
properties = {}
for property in propertyNames:
    properties[property] = nodes[property][nodesStart:nodesEnd]

# Identify the central galaxy, and extract the relevant tree.
indexCentral = np.nonzero(properties['nodeIsIsolated'] == 1)[0][0]
for property in propertyNames:
    properties[property] = properties[property][indexCentral]
    
# Generate a color map to color nodes by their star formation rate.
colorMap    = mpl.colormaps['coolwarm_r']
countColors = 100
colors      = colorMap(np.linspace(0,1,countColors))

# Generate a range spanning the number of nodes in tree.
indexNodes = np.array(list(range(len(properties['galaxyMergerTreeTime']))))

# Find the start and indexEnd index for each branch of the tree. This is found by simply cumulating the number of nodes in each
# branch.
indexEnd   = np.cumsum(properties['galaxyMergerTreeCount'])
indexStart = np.concatenate((np.array([0]),indexEnd[0:len(indexEnd)-1]))

# Initialize arrays of nodes (x,y positions, size, and color) and edges (tuples giving pairs of nodes that are connected by an
# edge).
edges      =          []
nodesX     = np.array([]                 )
nodesY     = np.array([]                 )
nodesIndex = np.array([], dtype=np.uint32)
nodesSize  = np.array([]                 )
nodesColor =          []
countNodes = 0

# Find maximum mass and star formation rate so that we can scale positions and colors to these.
massMaximum = np.max(properties['galaxyMergerTreeMassStellarTotal'      ])
sfrMaximum  = np.max(properties['galaxyMergerTreeTotalStarFormationRate'])

# Set size scale.
sizeScale = 300.0

# Iterate over branches.
for indexBranch in range(len(properties['galaxyMergerTreeCount'])):
    for indexNode in range(indexStart[indexBranch],indexEnd[indexBranch]):
        # Determine a color for this node based no its star formation rate.
        nodeColor = colors[int((countColors-1)*properties['galaxyMergerTreeTotalStarFormationRate'][indexNode]/sfrMaximum)]
        # Determine a size for the node based on its mass.
        nodeSize  = sizeScale*(properties['galaxyMergerTreeMassStellarTotal'][indexNode]/massMaximum)**(1.0/3.0)
        # Append node properties to our lists. The x positions are initially all set to the (arbitary) value of 0.
        nodesIndex = np.append(nodesIndex,countNodes                                   )
        nodesX     = np.append(nodesX    ,0.0                                          )
        nodesY     = np.append(nodesY    ,properties['galaxyMergerTreeTime'][indexNode])
        nodesSize  = np.append(nodesSize ,nodeSize                                     )
        nodesColor.append(nodeColor)
        # Add an edge (except for the final node in the branch) connecting to the next node.
        if indexNode < indexEnd[indexBranch]-1:
            edges.append((countNodes,countNodes+1))
        # Increment the count of nodes.
        countNodes += 1
        
# Set the fractional x-offset between nodes about to merge.
xOffsetMerger = 0.05

# Iterate over mergers in order of decreasing time.
order = np.flip(np.argsort(properties['galaxyMergerTreeMergeTime']))
for indexMerger in range(len(order)):
    # Find the nodes participating in this merger.
    indexSatellite = np.nonzero(properties['galaxyMergerTreeNodeIndex'] == properties['galaxyMergerTreeMergeSatelliteIndex'][order[indexMerger]])[0][0]
    indexCentral   = np.nonzero(properties['galaxyMergerTreeNodeIndex'] == properties['galaxyMergerTreeMergeCentralIndex'  ][order[indexMerger]])[0][0]
    # Find the corresponding branches for these nodes. For the central we find just the final node in the branch.
    branchSatellite  = (indexNodes < indexEnd[indexSatellite]) & (indexNodes >= indexStart[indexSatellite])
    branchCentralEnd = (indexNodes < indexEnd[indexCentral  ]) & (indexNodes >= indexStart[indexCentral  ]) & (properties['galaxyMergerTreeTime'] == properties['galaxyMergerTreeMergeTime'][order[indexMerger]])
    # Determine the amount by which to shift the satellite branch in the x-direction.
    xShift = (properties['galaxyMergerTreeMassStellarTotal'][branchSatellite][-1]/massMaximum)**(1.0/2.0)
    # Find all othr nodes that must be shift to make space for this branch.
    toShift = (nodesX > nodesX[branchCentralEnd][0]) & (nodesY < nodesY[branchSatellite][-1])
    # Shift all nodes (both on our satellite branch and others that must move to make space for it). The shift increases with time
    # before the merger to create a track that merges onto the central.
    nodesX[toShift        ] +=                             xShift*(1.0-(1.0-xOffsetMerger)*(nodesY[toShift        ]/nodesY[branchSatellite][-1])**4)   
    nodesX[branchSatellite]  = nodesX[branchCentralEnd][0]+xShift*(1.0-(1.0-xOffsetMerger)*(nodesY[branchSatellite]/nodesY[branchSatellite][-1])**4)   
    # Add an edge connecting the two nodes that merge.
    edges.append((nodesIndex[branchSatellite][-1],nodesIndex[branchCentralEnd][0]))

# Create a figure.
plt.rcParams['font.size'] = 18
fig, ax = plt.subplots(figsize=(10,15))
# Iterate over edges, drawing lines between the connected nodes.
for indexEdge in range(len(edges)):
    edgeX = (nodesX[edges[indexEdge][0]],nodesX[edges[indexEdge][1]])
    edgeY = (nodesY[edges[indexEdge][0]],nodesY[edges[indexEdge][1]])
    ax.plot(edgeX,edgeY,marker='',linestyle='-',color='#aaaaaa',zorder=10)
# Add points showing the nodes.
ax.scatter(nodesX,nodesY,nodesSize,zorder=20,color=nodesColor)
# Add a color bar.
sm = mpl.cm.ScalarMappable(cmap=colorMap, norm=mpl.colors.Normalize(0, sfrMaximum/1.0e9))
fig.colorbar(sm,ax=ax,label='SFR [$\mathrm{M}_\odot/$yr]',pad=0.15)
# Add redshift axis.
cosmo = cosmology.setCosmology('planck18')
timeMinimum = np.min(nodesY)*0.95
timeMaximum = np.max(nodesY)*1.05
redshiftMinimum = cosmo.age(timeMinimum,inverse=True)
redshiftMaximum = cosmo.age(timeMaximum,inverse=True)
ax2 = ax.twinx()
ax.set_ylim([timeMinimum,timeMaximum])
ax2.set_ylim([redshiftMinimum,redshiftMaximum])
ax2.set_ylabel('$z$')
# Remove ticks and labels from the x-axis - position on the x-axis is arbitrary.
ax.tick_params(left = True, right = False , labelleft = True , 
                labelbottom = False, bottom = False)
# Label the y-axis.
ax.set_ylabel('time [Gyr]')
# Plot points to indicate masses.
logMass = int(np.log10(massMaximum))+1
xMaximum  = np.max(nodesX)
xPosition = np.min(nodesX)
yPosition = timeMinimum*0.8
ax.annotate("$\log_{10}(M_\star/\mathrm{M}_\odot)$",xy=((xMaximum-xPosition)/2.0,yPosition+0.02),annotation_clip=False,horizontalalignment='center')
for i in range(3):
    logMass -= 1
    xPosition += xMaximum/4.0
    mass = 10.0**logMass
    pointSize = sizeScale*(mass/massMaximum)**(1.0/3.0)
    ax.scatter(xPosition, yPosition, pointSize, color='b', clip_on=False)
    ax.annotate(str(logMass),xy=(xPosition,yPosition-0.03),annotation_clip=False,horizontalalignment='center')
    
# Save our figure.
fig.savefig(args.outputFileName)
