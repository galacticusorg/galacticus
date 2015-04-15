# Kernel density estimation on BIE statelog files.
# Andrew Benson (10-June-2012)
from optparse import OptionParser
from scipy import stats, mgrid, c_, reshape, sqrt
import numpy as np
import sys
from string import join

# Parse arguments.
usage = "usage: %prog [options] statelog column1 [column2 [column3]]"
parser = OptionParser(usage)
parser.add_option("-o", "--output", dest="output",
                  help="write density to FILE", metavar="FILE")
parser.add_option("--ngood", dest="ngood", type="int",
                  help="number of states to use from end of statelog")
parser.add_option("--ngrid", dest="ngrid", type="int", default=100,
                  help="number of points to use in output grid")
(options, args) = parser.parse_args()
if len(args) != 2 and len(args) != 3 and len(args) != 4:
    parser.error("incorrect number of arguments")
statelog = args[0]
column1  = args[1]
if len(args) == 4:
    dimensions = 3
    columns    = tuple([int(args[1]),int(args[2]),int(args[3])])
elif len(args) == 3:
    dimensions = 2
    columns    = tuple([int(args[1]),int(args[2])])
else:
    dimensions = 1
    columns    = tuple([int(args[1])])

# Read data.
state = np.loadtxt(statelog,skiprows=0,usecols=columns)

# Find range of states to use.
rangeStart = 0
rangeEnd   = state.shape[0]
if options.ngood:
    rangeStart = -options.ngood
    rangeEnd   = -1

# Extract required states and find ranges.
if dimensions == 1:
    stateX = state[rangeStart:rangeEnd]
else:
    stateX = state[rangeStart:rangeEnd,0]
xRange = stateX.max()-stateX.min()
xMin   = stateX.min()-0.05*xRange
xMax   = stateX.max()+0.05*xRange

if dimensions >= 2:
    stateY = state[rangeStart:rangeEnd,1]
    yRange = stateY.max()-stateY.min()
    yMin   = stateY.min()-0.05*yRange
    yMax   = stateY.max()+0.05*yRange

if dimensions >= 3:
    stateZ = state[rangeStart:rangeEnd,2]
    zRange = stateZ.max()-stateZ.min()
    zMin   = stateZ.min()-0.05*zRange
    zMax   = stateZ.max()+0.05*zRange

# Create grid, positions and values.
ngrid = options.ngrid*sqrt(-1)
if ( dimensions == 1 ):
    X         = mgrid[xMin:xMax:ngrid                                  ]
    positions = c_   [X.ravel()                                        ]
    values    = c_   [stateX                                           ]
elif ( dimensions == 2 ):
    X, Y      = mgrid[xMin:xMax:ngrid, yMin:yMax:ngrid                 ]
    positions = c_   [X.ravel()      , Y.ravel()                       ]
    values    = c_   [stateX         , stateY                          ]
else:
    X, Y, Z   = mgrid[xMin:xMax:ngrid, yMin:yMax:ngrid, zMin:zMax:ngrid]
    positions = c_   [X.ravel()      , Y.ravel()      , Z.ravel()      ]
    values    = c_   [stateX         , stateY         , stateZ         ]

# Perform a kernel density estimator on the results.
kernel = stats.kde.gaussian_kde(values.T)
rho    = reshape(kernel(positions.T).T, X.T.shape)

# Output data.
if dimensions == 1:
    D = np.column_stack((X.flat,rho.flat))
elif ( dimensions == 2 ):
    D = np.column_stack((X.flat,Y.flat,rho.flat))
else:
    D = np.column_stack((X.flat,Y.flat,Z.flat,rho.flat))
if options.output:
    f = open(options.output, 'w')
else:
    f = sys.stdout
# For a 3D mesh, we output the header required for reading by the Ifrit package.
if dimensions == 3:
    f.write(str(options.ngrid)+" "+str(options.ngrid)+" "+str(options.ngrid)+"\n")
for d in D:
        f.write(join(str(x) for x in d)+"\n")
f.close()
