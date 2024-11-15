#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt

# Make a simple diagram showing a galaxy merger tree.
# Andrew Benson (18-October-2024)

# Open the model output and extract merger tree datasets.
model         = h5py.File("galaxyMergerTree.hdf5","r")
nodes         = model['Outputs/Output1/nodeData']
propertyNames = [ 'galaxyMergerTreeNodeIndex'          , 'galaxyMergerTreeCount'    , 'galaxyMergerTreeTotalStarFormationRate',
                  'galaxyMergerTreeMassStellarTotal'   , 'galaxyMergerTreeTime'     , 'galaxyMergerTreeMergeCentralIndex'     ,
                  'galaxyMergerTreeMergeSatelliteIndex', 'galaxyMergerTreeMergeTime', 'nodeIsIsolated'                          ]
properties = {}
for property in propertyNames:
    properties[property] = nodes[property][:]

# Identify the central galaxy, and extract the relevant tree.
indexCentral = np.nonzero(properties['nodeIsIsolated'] == 1)[0][0]
for property in propertyNames:
    properties[property] = properties[property][indexCentral]

# Generate a color map to color nodes by their star formation rate.
colorMap    = mpl.colormaps['coolwarm']
countColors = 100
colors      = colorMap(np.linspace(1,0,countColors))

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

# Iterate over branches.
for indexBranch in range(len(properties['galaxyMergerTreeCount'])):
    for indexNode in range(indexStart[indexBranch],indexEnd[indexBranch]):
        # Determine a color for this node based no its star formation rate.
        nodeColor = colors[int((countColors-1)*properties['galaxyMergerTreeTotalStarFormationRate'][indexNode]/sfrMaximum)]
        # Determine a size for the node based on its mass.
        nodeSize  = 100.0*(properties['galaxyMergerTreeMassStellarTotal'][indexNode]/massMaximum)**(1.0/3.0)
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
fig, ax = plt.subplots(figsize=(5,15))
# Iterate over edges, drawing lines between the connected nodes.
for indexEdge in range(len(edges)):
    edgeX = (nodesX[edges[indexEdge][0]],nodesX[edges[indexEdge][1]])
    edgeY = (nodesY[edges[indexEdge][0]],nodesY[edges[indexEdge][1]])
    ax.plot(edgeX,edgeY,marker='',linestyle='-',color='#aaaaaa',zorder=10)
# Add points showing the nodes.
ax.scatter(nodesX,nodesY,nodesSize,zorder=20,color=nodesColor)
# Remove ticks and labels from the x-axis - position on the x-axis is arbitrary.
ax.tick_params(left = True, right = False , labelleft = True , 
                labelbottom = False, bottom = False)
# Label the y-axis.
ax.set_ylabel('time [Gyr]')
# Save our figure.
fig.savefig('galaxyMergerTree.pdf')
