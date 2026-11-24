#!/usr/bin/env python3
"""Exercise both sides of a resolved change of host in a 19-node forest."""
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

import h5py
import numpy as np
from scipy.interpolate import CubicHermiteSpline

# Contributors: Andrew Robertson, Codex.


def main():
    suite = Path(__file__).resolve().parent
    model_root = Path.cwd()
    if model_root.name == "testSuite":
        model_root = model_root.parent
    work = model_root / "testSuite/outputs/positionInterpolatedBranchJump"
    work.mkdir(parents=True, exist_ok=True)
    forest = work / "forest.hdf5"
    output = work / "model.hdf5"

    # A and B are stationary centrals. S is isolated at step 0, a satellite of A at
    # steps 1 and 2, and of B at steps 3 through 5. Only B remains at step 6.
    host_a, host_b, satellite = 1, 2, 3
    a = 0.40 + 0.01 * np.arange(7)
    origin = np.array([20., 30., 40.])
    satellite_offset = np.array([[x, 0.02, 0.01] for x in (0.06, 0.04, 0.03, 0.18, 0.19, 0.20)])
    velocity = np.array([[0., 50. + step, 10. + 3 * step] for step in range(6)])
    initial_mass = {host_a: 1.e12, host_b: 2.e12, satellite: 8.e10}
    records = []
    for step in reversed(range(7)):
        for halo in ([host_b] if step == 6 else [host_a, host_b, satellite]):
            node_id = 10 * step + halo
            host = halo
            offset = np.array([0.25, 0., 0.]) if halo == host_b else np.zeros(3)
            peculiar_velocity = np.zeros(3)
            if halo == satellite:
                if step > 0:
                    host = host_a if step <= 2 else host_b
                offset = satellite_offset[step]
                peculiar_velocity = velocity[step]
            descendant = node_id + 10
            if step == 5:
                descendant = 62
            if step == 6:
                descendant = -1
            records.append((node_id, 10 * step + host, descendant, a[step], origin + offset,
                            peculiar_velocity, initial_mass[halo] * (1 + 0.02 * step)))

    with h5py.File(suite / "data/mergerTrees/orbitOnZAxis.hdf5") as template, h5py.File(forest, "w") as f:
        f.attrs.update(template.attrs)
        for name in ("units", "cosmology", "groupFinder", "simulation"):
            template.copy(name, f)
        index = f.create_group("forestIndex")
        index["firstNode"] = [0]
        index["forestIndex"] = [62]
        index["numberOfNodes"] = [len(records)]
        halos = f.create_group("forestHalos")
        halos.attrs.update(template["forestHalos"].attrs)
        columns = ("nodeIndex", "hostIndex", "descendantIndex", "expansionFactor", "position", "velocity", "nodeMass")
        for name, values in zip(columns, zip(*records)):
            halos[name] = np.asarray(values)
        halos["scaleRadius"] = np.full(len(records), 0.01)

    # Use the same reader and ordinary snapshot output configuration as the orphan trajectory test.
    params = ET.parse(suite / "parameters/positionInterpolatedOrphanDecay.xml")
    params.find("mergerTreeConstructor/fileNames").set("value", str(forest))
    params.find("outputFileName").set("value", str(output))
    params.find('nodeOperator/nodeOperator[@value="positionInterpolated"]/orphanTrajectory').set("value", "host")
    samples = np.linspace(a[1], a[4], 7)  # four knots and three intervening outputs
    params.find("outputTimes/redshifts").set("value", " ".join(str(1 / scale - 1) for scale in samples))
    config = work / "parameters.xml"
    params.write(config)
    subprocess.run(["./Galacticus.exe", str(config)], cwd=model_root, check=True)

    with h5py.File(output) as f, h5py.File(forest) as inp:
        if "runEndTime" not in f["Version"].attrs:
            print("FAILED: Galacticus output is missing the end-of-run marker")
            return
        outputs = sorted(f["Outputs"].values(), key=lambda g: g.attrs["outputTime"])
        np.testing.assert_allclose([g.attrs["outputExpansionFactor"] for g in outputs], samples, rtol=1e-6)
        times = np.array([g.attrs["outputTime"] for g in outputs])
        units = outputs[0].attrs
        hubble = inp["cosmology"].attrs["HubbleParam"]
        length = inp["units"].attrs["lengthUnitsInSI"] / units["distanceUnitInSI"] / hubble
        rate = 1000 * units["timeUnitInSI"] / units["distanceUnitInSI"]
        positions, velocities = [], []
        for group in outputs:
            data = group["nodeData"]
            i, = np.flatnonzero(data["nodeIsIsolated"][:] == 0)
            positions.append([data["positionPosition" + c][i] for c in "XYZ"])
            velocities.append([data["positionVelocity" + c][i] for c in "XYZ"])
        positions, velocities = np.asarray(positions), np.asarray(velocities)
        measured = (origin + satellite_offset[1:5]) * length
        np.testing.assert_allclose(positions[::2], measured * a[1:5, None], rtol=0, atol=1e-5)
        np.testing.assert_allclose(velocities[::2], velocity[1:5], rtol=0, atol=1e-4)

        # Both stationary hosts and every measured satellite velocity have zero x
        # velocity. A spurious displacement of the host frame violates this invariant.
        np.testing.assert_allclose(velocities[[1, 5], 0], 0., rtol=0, atol=1e-4)

        # Across the transfer, test S in absolute comoving coordinates. The cubic's
        # tangent at each endpoint is the peculiar velocity divided by the scale factor.
        cross = CubicHermiteSpline(times[[2, 4]], measured[1:3], rate * velocity[2:4] / a[[2, 3], None])
        np.testing.assert_allclose(positions[3], samples[3] * cross(times[3]), rtol=0, atol=1e-5)
        np.testing.assert_allclose(velocities[3], samples[3] * cross(times[3], 1) / rate, rtol=0, atol=.02)

    print("SUCCESS: resolved satellite trajectory on both sides of a branch jump")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, OSError, KeyError, ValueError, subprocess.CalledProcessError) as error:
        print(f"FAILED: resolved satellite branch-jump test: {error}")
