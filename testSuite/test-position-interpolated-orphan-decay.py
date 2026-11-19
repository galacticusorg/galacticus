#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np

# Test the "decay" orphan trajectory of the `positionInterpolated` node operator. Once the position history of a satellite is
# exhausted the orphan should follow a decaying orbit about its host, on which the square of the orbital radius decreases linearly
# in time to zero at the time of merging, and the tangential speed is constant.
# Andrew Benson, Claude (22-October-2026)

# The conversion from Mpc/(km/s) to Gyr, matching `MpcPerKmPerSToGyr` in Galacticus (the GSL parsec, and the sidereal year).
megaParsecPerKmPerSToGyr = 3.08567758135e22/1.0e3/(3.15581497635456e7*1.0e9)

def run(label, changes):
    """Run the model with the given parameter changes, returning the name of the output file (or None on failure)."""
    outputFileName = f"outputs/positionInterpolatedOrphanDecay_{label}.hdf5"
    changeFileName = f"outputs/positionInterpolatedOrphanDecay_{label}.changes.xml"
    with open(changeFileName, "w") as changeFile:
        changeFile.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<changes>
{changes}  <change type="replace" path="outputFileName">
    <outputFileName value="testSuite/{outputFileName}"/>
  </change>
</changes>
""")
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/parameters/positionInterpolatedOrphanDecay.xml testSuite/{changeFileName}", shell=True)
    if status.returncode != 0:
        print(f"FAILED: model run ({label} orphan trajectory)")
        return None
    return outputFileName

def extract(fileName):
    """Extract the orbital position, velocity, and satellite status of every satellite at each output."""
    model  = h5py.File(fileName, "r")
    states = {}
    for output in model["Outputs"].values():
        time     = float(output.attrs["outputTime"])
        nodes    = output["nodeData"]
        position = np.column_stack([nodes["positionOrbital"+axis][:] for axis in "XYZ"])
        velocity = np.column_stack([nodes["velocityOrbital"+axis][:] for axis in "XYZ"])
        for i in np.where(nodes["nodeIsIsolated"][:] == 0)[0]:
            states[(int(nodes["nodeIndex"][i]), time)] = (int(nodes["satelliteStatus"][i]), position[i], velocity[i])
    model.close()
    return states

# Run the model with "decay" orphan trajectories, and, for comparison, with "host" orphan trajectories.
subprocess.run("mkdir -p outputs", shell=True)
fileNameDecay = run("decay", "")
fileNameHost  = run("host" , """  <change type="update" path="nodeOperator/nodeOperator[@value='positionInterpolated']/orphanTrajectory" value="host"/>
""")
if fileNameDecay is None or fileNameHost is None:
    sys.exit(0)
statesDecay = extract(fileNameDecay)
statesHost  = extract(fileNameHost )

# Find the states of orphans on decaying orbits. These are orphans (satellite status 2, indicating a satellite whose position
# history is exhausted) whose positions differ from those found with the "host" trajectory. Orphans following their host (e.g.
# those with no finite time of merging) are thereby excluded. Satellites must be identical in the two models until they become
# orphans.
orphans   = {}
timesLast = {}
failures  = 0
for (index, time), (status, position, velocity) in statesDecay.items():
    timesLast[index] = max(time, timesLast.get(index, -1.0))
    if (index, time) not in statesHost:
        print(f"FAILED: satellite {index} at t={time:.4f} Gyr is present only with the \"decay\" orphan trajectory")
        failures += 1
        continue
    separation = np.linalg.norm(position-statesHost[(index, time)][1])
    if status == 2:
        if separation > 1.0e-8:
            orphans.setdefault(index, []).append((time, position, velocity))
    elif separation > 1.0e-9:
        print(f"FAILED: satellite {index} at t={time:.4f} Gyr differs between orphan trajectories before becoming an orphan ({separation:.3e} Mpc)")
        failures += 1

# Check the decaying orbit of each orphan.
countChecked = 0
for index, states in sorted(orphans.items()):
    states = sorted(states, key=lambda state: state[0])
    if len(states) < 3:
        continue
    time     = np.array([state[0] for state in states])
    radius   = np.array([np.linalg.norm(state[1]) for state in states])
    speedTangential = np.array([np.linalg.norm(np.cross(state[1], state[2]))/np.linalg.norm(state[1]) for state in states])
    # The square of the radius must be linear in time, reaching zero no earlier than the final orphan output, and no later than
    # the final output at which the node remains a satellite.
    slope, intercept = np.polyfit(time, radius**2, 1)
    residual         = np.max(np.abs(np.polyval([slope, intercept], time)-radius**2))/np.max(radius**2)
    timeZero         = -intercept/slope
    if residual > 1.0e-6:
        print(f"FAILED: orphan {index}: square of orbital radius is not linear in time (maximum relative residual {residual:.3e})")
        failures += 1
    if slope >= 0.0 or timeZero < time[-1] or timeZero > timesLast[index]*(1.0+1.0e-6):
        print(f"FAILED: orphan {index}: orbit decays to the host center at t={timeZero:.4f} Gyr, outside the range {time[-1]:.4f}-{timesLast[index]:.4f} Gyr")
        failures += 1
    # The tangential speed must be constant.
    variation = (speedTangential.max()-speedTangential.min())/speedTangential.max() if speedTangential.max() > 0.0 else 0.0
    if variation > 1.0e-6:
        print(f"FAILED: orphan {index}: tangential speed is not constant (relative variation {variation:.3e})")
        failures += 1
    # The radial speed must be consistent with the decrease in radius, dr/dt = (dr²/dt)/(2r).
    speedRadial          = np.array([np.dot(state[1], state[2])/np.linalg.norm(state[1]) for state in states])
    speedRadialPredicted = slope/(2.0*radius)*megaParsecPerKmPerSToGyr
    errorSpeedRadial     = np.max(np.abs(speedRadial/speedRadialPredicted-1.0))
    if errorSpeedRadial > 1.0e-6:
        print(f"FAILED: orphan {index}: radial speed is inconsistent with the decay of the orbital radius (maximum relative error {errorSpeedRadial:.3e})")
        failures += 1
    # The orbit must remain in the plane defined by the initial position and velocity.
    normal      = np.cross(states[0][1], states[0][2])
    normal     /= np.linalg.norm(normal)
    errorPlane  = max(abs(np.dot(state[1], normal))/np.linalg.norm(state[1]) for state in states)
    if errorPlane > 1.0e-6:
        print(f"FAILED: orphan {index}: orbit does not remain in its initial plane (maximum fractional offset {errorPlane:.3e})")
        failures += 1
    # The orbital phase must advance as the integral of the angular speed, (tangential speed)/(radius), giving
    # θ(t) = −2 v_t √(t_m−t)/√(−dr²/dt) + constant. Signed angles in the orbital plane are compared, modulo 2π, between successive
    # outputs.
    phase      = -2.0*speedTangential[0]/megaParsecPerKmPerSToGyr*np.sqrt(np.maximum(timeZero-time, 0.0))/np.sqrt(-slope) if slope < 0.0 else np.zeros(len(time))
    errorPhase = 0.0
    for i in range(len(states)-1):
        anglePredicted = phase[i+1]-phase[i]
        angleMeasured  = np.arctan2(np.dot(np.cross(states[i][1], states[i+1][1]), normal), np.dot(states[i][1], states[i+1][1]))
        errorPhase     = max(errorPhase, abs((anglePredicted-angleMeasured+np.pi) % (2.0*np.pi)-np.pi))
    if errorPhase > 1.0e-6:
        print(f"FAILED: orphan {index}: orbital phase does not advance at the angular speed of the orbit (maximum error {errorPhase:.3e} radians)")
        failures += 1
    countChecked += 1

# Require that a reasonable number of orphans were checked, so that the test can not pass vacuously.
countRequired = 5
if countChecked < countRequired:
    print(f"FAILED: only {countChecked} orphans on decaying orbits were found (at least {countRequired} expected)")
    failures += 1
if failures == 0:
    print(f"SUCCESS: decaying orbits of {countChecked} orphans")
sys.exit(0)
