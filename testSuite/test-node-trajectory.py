#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Check recording of node evolutionary trajectories by the `outputTrajectory` node operator.
# Andrew Benson

# Event identifiers, as enumerated in source/nodes/operators/output/trajectory_events.F90.
eventNodeInitialize                = 0
eventDifferentialEvolutionPostStep = 1
eventDifferentialEvolutionPost     = 2
eventNodePromote                   = 3

# Run the base model (no trajectory recording), and the identical model with trajectory recording enabled.
for name in ( "nodeTrajectoryBase", "nodeTrajectory" ):
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/parameters/{name}.xml", shell=True)
    if status.returncode != 0:
        print(f"FAILED: {name} model run")
        sys.exit(0)
print("SUCCESS: model runs")

modelBase = h5py.File("outputs/nodeTrajectoryBase.hdf5", "r")
model     = h5py.File("outputs/nodeTrajectory.hdf5"    , "r")

# Recording trajectories must not perturb the model in any way. In particular, the `mergerTreeExtraOutput` event hook (which
# modifies the star formation histories of nodes) must not be triggered by trajectory output. Verify that every regular output
# dataset is unchanged.
nodeDataBase = modelBase["Outputs/Output1/nodeData"]
nodeData     = model    ["Outputs/Output1/nodeData"]
if set(nodeDataBase.keys()) != set(nodeData.keys()):
    print("FAILED: regular output datasets differ between models")
    sys.exit(0)
perturbed = [
    name for name in nodeDataBase.keys()
    if not np.array_equal(nodeDataBase[name][:], nodeData[name][:])
]
if perturbed:
    print(f"FAILED: trajectory recording perturbed the model; datasets differ: {perturbed}")
else:
    print("SUCCESS: trajectory recording does not perturb the model")

# The trajectory itself must have been recorded.
if "nodeTrajectories/Trajectory1/nodeData" not in model:
    print("FAILED: trajectory group was not created")
    sys.exit(0)
print("SUCCESS: trajectory group created")

trajectory = model["nodeTrajectories/Trajectory1/nodeData"]
time       = trajectory["time"           ][:]
event      = trajectory["trajectoryEvent"][:]
nodeIndex  = trajectory["nodeIndex"      ][:]
treeIndex  = trajectory["mergerTreeIndex"][:]

# Records must be ordered in time, and span the full evolution of the node (from its own time of 1.0 Gyr to the final output
# time of 13.8 Gyr).
if np.all(np.diff(time) >= 0.0):
    print("SUCCESS: trajectory records are ordered in time")
else:
    print("FAILED: trajectory records are not ordered in time")
if time[0] <= 1.0 + 1.0e-6 and time[-1] >= 13.8 - 1.0e-6:
    print("SUCCESS: trajectory spans the full evolution")
else:
    print(f"FAILED: trajectory spans only {time[0]} to {time[-1]} Gyr")

# Records should have been made at ODE solver steps - which is the whole point of this class, and is the only way to obtain many
# more records than there are output times.
if np.count_nonzero(event == eventDifferentialEvolutionPostStep) > 1:
    print("SUCCESS: trajectory records were made at ODE solver steps")
else:
    print("FAILED: no trajectory records were made at ODE solver steps")
if np.any(event == eventNodeInitialize) and np.any(event == eventDifferentialEvolutionPost):
    print("SUCCESS: trajectory records were made at initialization and post-evolution")
else:
    print("FAILED: trajectory records missing for initialization and/or post-evolution")

# All records must refer to the tree that we asked to follow.
if np.all(treeIndex == 1):
    print("SUCCESS: all trajectory records belong to the requested tree")
else:
    print("FAILED: trajectory records belong to more than one tree")

# The node is promoted twice during evolution (node 1 → 2 at 6 Gyr, node 2 → 3 at 13.8 Gyr). The trajectory must follow the
# galaxy across those promotions - i.e. the node index recorded must change, and records must continue to be made after the
# first promotion.
if np.count_nonzero(event == eventNodePromote) == 2:
    print("SUCCESS: both promotions were recorded")
else:
    print(f"FAILED: {np.count_nonzero(event == eventNodePromote)} promotions recorded, expected 2")
timePromotionFirst = time[event == eventNodePromote][0] if np.any(event == eventNodePromote) else None
if timePromotionFirst is not None and np.count_nonzero(time > timePromotionFirst) > 1:
    print("SUCCESS: trajectory continues after promotion")
else:
    print("FAILED: trajectory does not continue after promotion")
# The node index recorded must change across the first promotion: before it the galaxy lives in node 1, after it in node 2. Note
# that a promotion record is made *before* the promotion takes place, so it carries the index of the node being promoted - and
# that node 3 never appears, because the second promotion (node 2 -> node 3) occurs at the final time, after which there is no
# further evolution to record.
if timePromotionFirst is not None:
    indicesBefore = set(np.unique(nodeIndex[time <= timePromotionFirst]).tolist())
    indicesAfter  = set(np.unique(nodeIndex[time >  timePromotionFirst]).tolist())
    if indicesBefore == {1} and indicesAfter == {2}:
        print("SUCCESS: trajectory follows the galaxy across the promotion into its parent node")
    else:
        print(f"FAILED: node indices before/after first promotion are {sorted(indicesBefore)}/{sorted(indicesAfter)}, expected [1]/[2]")

# The trajectory must be far more finely sampled than the regular output - which is the whole point of this class.
countRegular = len(nodeDataBase["basicMass"][:])
if len(time) > 10 * countRegular:
    print(f"SUCCESS: trajectory is finely sampled ({len(time)} records vs {countRegular} regular output row(s))")
else:
    print(f"FAILED: trajectory has only {len(time)} records vs {countRegular} regular output row(s)")

# The trajectory group must not carry the time-related attributes of a regular output, which would be meaningless for a group
# spanning many times.
attributes = set(model["nodeTrajectories/Trajectory1"].attrs.keys())
if attributes.isdisjoint({"outputTime", "outputExpansionFactor", "outputComovingDistance"}):
    print("SUCCESS: trajectory group omits the regular-output time attributes")
else:
    print("FAILED: trajectory group carries misleading regular-output time attributes")

modelBase.close()
model    .close()
