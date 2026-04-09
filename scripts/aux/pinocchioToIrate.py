#!/usr/bin/env python3
import sys
import os
import re
import numpy as np
import h5py

# Convert Pinocchio output catalogs to IRATE format.
# Andrew Benson (28-August-2014) [Python port]

if len(sys.argv) != 3:
    print("Usage: pinocchioToIrate.py <pinocchioDirectoryName> <irateFileName>", file=sys.stderr)
    sys.exit(1)

simulation_directory_name = sys.argv[1]
irate_file_name = sys.argv[2]

if not simulation_directory_name.endswith(os.sep):
    simulation_directory_name += os.sep

# Read the parameter file.
parameters = {}
parameter_file_path = os.path.join(simulation_directory_name, "parameter_file")
with open(parameter_file_path, 'r') as f:
    for line in f:
        match = re.match(r'^([^\#\s]+)\s+([^\s]+)', line)
        if match:
            parameter_name = match.group(1)
            parameter_value = match.group(2)
            parameters[parameter_name] = parameter_value

# Read the list of outputs.
outputs = []
output_list_path = os.path.join(simulation_directory_name, parameters['OutputList'])
with open(output_list_path, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        outputs.append(float(line))

outputs.sort()

# Open IRATE file.
with h5py.File(irate_file_name, 'w') as irate_file:
    # Iterate over outputs.
    for i, z_val in enumerate(outputs):
        snapshot_idx = i + 1
        redshift_label = f"{z_val:6.4f}"
        catalog_file_name = os.path.join(
            simulation_directory_name, 
            f"pinocchio.{redshift_label}.{parameters['RunFlag']}.catalog.out"
        )
        
        print(f"Reading {catalog_file_name}")
        
        try:
            # Column mapping (0-indexed): 1=Mass, 5=x, 6=y, 7=z, 8=vx, 9=vy, 10=vz
            data = np.loadtxt(catalog_file_name, usecols=(1, 5, 6, 7, 8, 9, 10), unpack=True)
            if data.ndim == 1:
                # Handle single-row case
                mass, x, y, z, vx, vy, vz = data.reshape(-1, 1)
            else:
                mass, x, y, z, vx, vy, vz = data
        except (OSError, ValueError) as e:
            print(f"Warning: Could not read {catalog_file_name}: {e}")
            continue

        # Construct 3D datasets.
        # PDL zeros(3, N) corresponds to (3, N) shape in numpy
        center = np.array([x, y, z])
        velocity = np.array([vx, vy, vz])

        # Snapshot group.
        snapshot_label = f"Snapshot{snapshot_idx:05d}"
        snap_grp = irate_file.create_group(snapshot_label)
        snap_grp.attrs['Redshift'] = float(z_val)

        # Halo catalog.
        halo_catalog = snap_grp.create_group('HaloCatalog')
        
        ds_center = halo_catalog.create_dataset('Center', data=center)
        ds_vel = halo_catalog.create_dataset('Velocity', data=velocity)
        ds_mass = halo_catalog.create_dataset('Mass', data=mass)

        # Set attributes.
        ds_center.attrs['unitname'] = "Mpc"
        ds_center.attrs['unitscgs'] = np.array([3.08568e+24, 0.0, -1.0])
        
        ds_vel.attrs['unitname'] = "km/s"
        ds_vel.attrs['unitscgs'] = np.array([1.00000e+05, 0.0, 0.0])
        
        ds_mass.attrs['unitname'] = "Msolar"
        ds_mass.attrs['unitscgs'] = np.array([1.98892e+33, 0.0, 0.0])

    # Cosmology.
    cosmology = irate_file.create_group('Cosmology')
    cosmology.attrs['HubbleParam'] = float(parameters['Hubble100'])
    cosmology.attrs['OmegaBaryon'] = float(parameters['OmegaBaryon'])
    cosmology.attrs['OmegaLambda'] = float(parameters['OmegaLambda'])
    cosmology.attrs['OmegaMatter'] = float(parameters['Omega0'])
    cosmology.attrs['sigma_8']     = float(parameters['Sigma8'])

    # Simulation properties.
    simulation = irate_file.create_group('SimulationProperties')
    simulation.attrs['boxSize'] = float(parameters['BoxSize'])
