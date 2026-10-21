#!/usr/bin/env python3
import subprocess
import sys
import h5py
import numpy as np
from scipy.integrate   import solve_ivp
from scipy.interpolate import CubicSpline

# Check the ram pressure history accumulated by the `ramPressureHistory` node operator - the time integral, I, and running
# maximum, P, of the ram pressure experienced by a satellite - against values computed independently from the satellite's
# recorded trajectory.
#
# A single satellite orbits on an eccentric orbit in a static host which carries a hot halo (see
# testSuite/parameters/ramPressureHistory.xml). The `outputTrajectory` node operator records, after every accepted ODE step, the
# satellite's position and velocity, the instantaneous ram pressure p, and the evolved I and P. The model is run twice:
#
#  * with efficiencyDecay=0, I must be the plain time integral of p, which is computed here by integrating a cubic spline through
#    the recorded values of p; and P must be the running maximum of p over the recorded steps - exactly, since P is updated from
#    the same state that is then recorded;
#
#  * with efficiencyDecay=1, I and P obey
#
#        dI/dt = p - I/T_orb,   dP/dt = -P/T_orb  (with P raised to p whenever p > P),
#
#    where T_orb = 2 pi / max(|r x v|/r^2, |r . v|/r^2) is computed here from the recorded position and velocity. The equation for I
#    is integrated here with SciPy's ODE solver, using cubic splines through the recorded p and 1/T_orb. P changes discontinuously
#    only at recorded steps, so it is stepped from one record to the next, decaying by the integral of 1/T_orb between them.
#
# Andrew Benson, Claude

# Event identifier for records made after each accepted ODE step, as enumerated in
# source/nodes/operators/output/trajectory_events.F90.
eventDifferentialEvolutionPostStep = 1

# Conversion from Mpc/(km/s) to Gyr, matching Galacticus' definitions of the parsec (GSL) and the (sidereal) year.
MpcPerKmPerSToGyr = 3.08567758135e22 / 1.0e3 / 3.15581497635456e16

# Tolerances, set from the measured agreement with a factor of ten or more in hand. The ODE steps, and so the records, are up to
# 0.1 Gyr apart, which limits the accuracy of the interpolated references: a trapezoidal rule over the records was in error by
# 7e-4, whereas the cubic splines used here agree to 4.5e-7 without decay. With decay the agreement is 1.0e-5 for I, limited by
# the spline through 1/T_orb, which has kinks where the larger of the two orbital frequencies changes.
toleranceIntegral      = 1.0e-5
toleranceIntegralDecay = 1.0e-4
toleranceMaximum       = 1.0e-4

def run(efficiencyDecay):
    """Run the model with the given decay efficiency, and return the recorded trajectory."""
    with open("parameters/ramPressureHistory.xml") as file:
        parameters = file.read()
    parameters = parameters.replace('<efficiencyDecay value="0.0"/>', f'<efficiencyDecay value="{efficiencyDecay}"/>')
    parameters = parameters.replace("testSuite/outputs/ramPressureHistory.hdf5", f"testSuite/outputs/ramPressureHistory{efficiencyDecay}.hdf5")
    with open(f"outputs/ramPressureHistory{efficiencyDecay}.xml", "w") as file:
        file.write(parameters)
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/outputs/ramPressureHistory{efficiencyDecay}.xml", shell=True)
    if status.returncode != 0:
        print(f"FAILED: model run with efficiencyDecay={efficiencyDecay}")
        return None
    print(f"SUCCESS: model run with efficiencyDecay={efficiencyDecay}")
    try:
        model = h5py.File(f"outputs/ramPressureHistory{efficiencyDecay}.hdf5", "r")
    except OSError as error:
        print(f"FAILED: unable to open the model output: {error}")
        return None
    if "nodeTrajectories/Trajectory1/nodeData" not in model:
        print("FAILED: trajectory was not recorded")
        return None
    trajectory = model["nodeTrajectories/Trajectory1/nodeData"]
    record     = trajectory["trajectoryEvent"][:] == eventDifferentialEvolutionPostStep
    data       = {
        "time"      : trajectory["time"                              ][:][record],
        "pressure"  : trajectory["satelliteRamPressure"              ][:][record],
        "maximum"   : trajectory["satelliteRamPressureMaximum"       ][:][record],
        "integral"  : trajectory["satelliteRamPressureTimeIntegrated"][:][record],
        "position"  : np.column_stack([trajectory[f"positionOrbital{axis}"][:][record] for axis in "XYZ"]),
        "velocity"  : np.column_stack([trajectory[f"velocityOrbital{axis}"][:][record] for axis in "XYZ"]),
    }
    model.close()
    # Verify the test's own preconditions: enough records to integrate over, and a ram pressure which varies substantially along
    # the orbit. Without these the test could pass while testing nothing.
    if len(data["time"]) < 100:
        print(f"FAILED: too few trajectory records ({len(data['time'])}) to test the time integral")
        return None
    if np.any(np.diff(data["time"]) <= 0.0):
        print("FAILED: trajectory records are not strictly ordered in time")
        return None
    if data["pressure"].min() <= 0.0 or data["pressure"].max() < 10.0 * data["pressure"].min():
        print(f"FAILED: test precondition not met -- ram pressure does not vary substantially along the orbit "
              f"({data['pressure'].min():.4e} to {data['pressure'].max():.4e})")
        return None
    return data

def check(name, model, reference, tolerance):
    """Compare a model quantity against its reference, relative to the largest reference value."""
    error = np.max(np.abs(model - reference)) / np.max(np.abs(reference))
    if error < tolerance:
        print(f"SUCCESS: {name} (maximum relative error {error:.2e})")
    else:
        print(f"FAILED: {name} (maximum relative error {error:.2e} exceeds tolerance {tolerance:.1e})")

# No decay: the plain time integral and running maximum.
data = run(0.0)
if data is not None:
    time      = data["time"    ]
    pressure  = data["pressure"]
    # The integral is zero at infall, before the first recorded step, so integrate from the first record using its value of I.
    spline    = CubicSpline(time, pressure)
    integral  = data["integral"][0] + np.array([spline.integrate(time[0], timeEnd) for timeEnd in time])
    check("time integral of ram pressure, without decay", data["integral"], integral, toleranceIntegral)
    # The maximum includes the initial state, which is not recorded; so allow the first record's maximum as a floor.
    maximum   = np.maximum.accumulate(np.maximum(pressure, data["maximum"][0]))
    check("running maximum of ram pressure, without decay", data["maximum"], maximum, toleranceMaximum)
    if np.all(data["maximum"] >= pressure * (1.0 - 1.0e-12)):
        print("SUCCESS: running maximum of ram pressure is never below the instantaneous ram pressure")
    else:
        print("FAILED: running maximum of ram pressure falls below the instantaneous ram pressure")

# Decay with efficiency 1: integrate the decaying equations independently.
data = run(1.0)
if data is not None:
    time      = data["time"    ]
    pressure  = data["pressure"]
    position  = data["position"]
    velocity  = data["velocity"]
    radius    = np.linalg.norm(position, axis=1)
    frequency = np.maximum(
        np.linalg.norm(np.cross(position, velocity), axis=1),
        np.abs(np.sum(position * velocity, axis=1))
    ) / radius**2 / MpcPerKmPerSToGyr
    splinePressure  = CubicSpline(time, pressure          )
    splineRateDecay = CubicSpline(time, frequency/2.0/np.pi)
    solution        = solve_ivp(
        lambda timeCurrent, integral: splinePressure(timeCurrent) - splineRateDecay(timeCurrent) * integral,
        (time[0], time[-1]),
        [data["integral"][0]],
        t_eval=time,
        rtol=1.0e-10,
        atol=1.0e-12 * np.max(data["integral"])
    )
    if not solution.success:
        print(f"FAILED: reference integration failed: {solution.message}")
        sys.exit(0)
    integral    = solution.y[0]
    maximum     = np.zeros_like(time)
    maximum [0] = data["maximum"][0]
    for i in range(1, len(time)):
        maximum[i] = max(pressure[i], maximum[i-1] * np.exp(-splineRateDecay.integrate(time[i-1], time[i])))
    check("time integral of ram pressure, with decay", data["integral"], integral, toleranceIntegralDecay)
    check("running maximum of ram pressure, with decay", data["maximum"], maximum, toleranceMaximum)
    # The decay must have had an effect, or the comparison above would not distinguish the two cases.
    if np.any(data["maximum"] < 0.9 * np.maximum.accumulate(pressure)):
        print("SUCCESS: decay reduces the running maximum of ram pressure below the peak ram pressure")
    else:
        print("FAILED: test precondition not met -- the running maximum of ram pressure has not decayed below the peak")
