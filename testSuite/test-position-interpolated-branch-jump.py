#!/usr/bin/env python3
"""Exercise the incoming trajectory and both resolved and orphaned host transfers."""
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

import h5py
import numpy as np
from scipy.interpolate import CubicHermiteSpline

# Contributors: Andrew Robertson, Codex.


def run_tree(suite, model_root, case, records, samples, trajectory):
    work = model_root / "testSuite/outputs/positionInterpolatedBranchJump" / case
    work.mkdir(parents=True, exist_ok=True)
    forest = work / "forest.hdf5"
    output = work / "model.hdf5"
    with h5py.File(suite / "data/mergerTrees/orbitOnZAxis.hdf5") as template, h5py.File(forest, "w") as f:
        f.attrs.update(template.attrs)
        for name in ("units", "cosmology", "groupFinder", "simulation"):
            template.copy(name, f)
        index = f.create_group("forestIndex")
        index["firstNode"] = [0]
        index["forestIndex"] = [max(row[0] for row in records)]
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
    params.find('nodeOperator/nodeOperator[@value="positionInterpolated"]/orphanTrajectory').set("value", trajectory)
    params.find("outputTimes/redshifts").set("value", " ".join(str(1 / scale - 1) for scale in samples))
    config = work / "parameters.xml"
    params.write(config)
    subprocess.run(["./Galacticus.exe", str(config)], cwd=model_root, check=True)

    with h5py.File(output) as model:
        np.testing.assert_equal("runEndTime" in model["Version"].attrs, True, err_msg=f"Incomplete output for {case}")
    return forest, output


def check_resolved(suite, model_root):
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

    samples = np.linspace(a[1], a[4], 7)  # four knots and three intervening outputs
    forest, output = run_tree(suite, model_root, "resolved", records, samples, "host")

    with h5py.File(output) as f, h5py.File(forest) as inp:
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

def check_orphan(suite, model_root, transfer):
    # Give B a measurable displacement and velocity, so using its later snapshot
    # at the handoff is detectable even if the orbit stays centered on B.
    origin = np.array([20., 30., 40.])
    records = []
    for step in reversed(range(9)):
        for halo in (1, 2, 3):
            if step > 2 and (halo == 3 or (transfer and halo == 1)):
                continue
            node = 10 * step + halo
            host, descendant = node, node + 10
            if step == 8 or (transfer and step == 2 and halo == 1):
                descendant = -1
            pos, vel = origin.copy(), np.zeros(3)
            if halo == 2:
                pos[0] += .25 + .012 * (step - 2)
                vel[0] = 60. + 2 * step
            if halo == 3:
                if step > 0:
                    host = 10 * step + 1
                if step == 2:
                    descendant = 32 if transfer else 31
                pos += [(0.06, 0.04, 0.03)[step], .02, .01]
                vel = np.array([0., 50. + step, 10. + 3 * step])
            mass = {1: 1.e12, 2: 2.e12, 3: 8.e10}[halo] * (1 + .02 * step)
            records.append((node, host, descendant, .40 + .01 * step, pos, vel, mass))

    # Include the last measurement, both sides of the handoff, and the endpoint
    # of the finite decay clock (the last snapshot in this forest).
    samples = np.array([.415, .42, .4225, .425, .4275, .43, .48])
    name = "orphan-transfer" if transfer else "ordinary-disappearance"
    forest, output = run_tree(suite, model_root, name, records, samples, "decay")
    with h5py.File(output) as f, h5py.File(forest) as inp:
        groups = sorted(f["Outputs"].values(), key=lambda g: g.attrs["outputTime"])
        np.testing.assert_allclose([g.attrs["outputExpansionFactor"] for g in groups], samples)
        times = np.array([g.attrs["outputTime"] for g in groups])
        units = groups[0].attrs
        rate = 1000 * units["timeUnitInSI"] / units["distanceUnitInSI"]
        length = inp["units"].attrs["lengthUnitsInSI"] / units["distanceUnitInSI"] / inp["cosmology"].attrs["HubbleParam"]
        states = {row[0]: row for row in records}
        host_id = 22 if transfer else 21
        h0, h1, sat = (states[key] for key in (host_id, host_id + 10, 23))
        center = CubicHermiteSpline(times[[1, 5]], np.array([h0[4], h1[4]]) * length,
                                    rate * np.array([h0[5], h1[5]]) / samples[[1, 5], None])
        offset = (sat[4] - h0[4]) * length * samples[1]
        radius = np.linalg.norm(offset)
        radial = offset / radius
        relative_v = sat[5] - h0[5]
        tangent_v = relative_v - radial * np.dot(relative_v, radial)
        speed = np.linalg.norm(tangent_v)
        tangent = tangent_v / speed
        duration = times[-1] - times[1]
        for k, group in enumerate(groups[:-1]):
            data = group["nodeData"]
            i, = np.flatnonzero(data["nodeIndex"][:] == 3)
            x, v = [np.array([data[p + axis][i] for axis in "XYZ"])
                   for p in ("positionPosition", "positionVelocity")]
            if k == 0:
                np.testing.assert_allclose(v[0], 0., atol=1e-4)
                continue
            if k == 1:
                np.testing.assert_allclose(x, sat[4] * length * samples[k], rtol=0, atol=2e-6)
                continue
            np.testing.assert_equal(data["satelliteStatus"][i], 2, err_msg=f"{name} at a={samples[k]}")
            fraction = np.sqrt((times[-1] - times[k]) / duration)
            angle = 2 * speed * rate * duration / radius * (1 - fraction)
            rhat = radial * np.cos(angle) + tangent * np.sin(angle)
            that = tangent * np.cos(angle) - radial * np.sin(angle)
            expected_x = samples[k] * center(times[k]) + radius * fraction * rhat
            expected_v = (samples[k] * center(times[k], 1) / rate + speed * that
                          - radius / (2 * duration * fraction * rate) * rhat)
            np.testing.assert_allclose(x, expected_x, rtol=0, atol=2e-6, err_msg=f"{name} position at a={samples[k]}")
            np.testing.assert_allclose(v, expected_v, rtol=0, atol=.01, err_msg=f"{name} velocity at a={samples[k]}")


def main():
    suite = Path(__file__).resolve().parent
    model_root = Path.cwd()
    if model_root.name == "testSuite":
        model_root = model_root.parent
    check_resolved(suite, model_root)
    for transfer in (False, True):
        check_orphan(suite, model_root, transfer)
    print("SUCCESS: resolved host transfer, orphan host transfer and ordinary disappearance")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, OSError, KeyError, ValueError, subprocess.CalledProcessError) as error:
        print(f"FAILED: satellite branch-jump test: {error}")
