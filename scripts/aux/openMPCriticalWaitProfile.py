#!/usr/bin/env python3
import sys
import re
import argparse

try:
    import h5py
    import numpy as np
except ImportError:
    print("This script requires the 'h5py' and 'numpy' packages. Install them with: pip install h5py numpy", file=sys.stderr)
    sys.exit(1)

# Extract metadata on OpenMP critical section wait times from a Galacticus run.
# Andrew Benson (30-November-2016)

parser = argparse.ArgumentParser(prog='openMPCriticalWaitProfile.py', description='Extract metadata on OpenMP critical section wait times from a Galacticus run.')
parser.add_argument('modelFileName', help='Galacticus model HDF5 output file')
args = parser.parse_args()

# Extract meta-data from file.
with h5py.File(args.modelFileName, 'r') as model_file:
    open_mp_group              = model_file['metaData']['openMP']
    critical_section_names     = open_mp_group['criticalSectionNames'][:]
    critical_section_wait_times = open_mp_group['criticalSectionWaitTimes'][:]

# Find total critical section wait time.
wait_time_total = critical_section_wait_times.sum()

# If there is no total wait time (or it is non-positive), there is nothing meaningful to report.
if wait_time_total <= 0:
    print("Total wait time for critical section across all threads is non-positive; no significant contributors to report.")
    sys.exit(0)

# Report.
print(f"Total wait time for critical section across all threads = {wait_time_total} s")
print("Significant contributing critical sections are:")

wait_time_ranks = np.argsort(critical_section_wait_times)
n_sections = len(critical_section_wait_times)
i = 0
wait_time_percentage = 100.0
max_sections = min(10, n_sections)
while (i < max_sections or wait_time_percentage > 10.0) and i < n_sections:
    i += 1
    wait_time_percentage = 100.0 * critical_section_wait_times[wait_time_ranks[-i]] / wait_time_total
    # Decode byte string if needed and strip surrounding quotes/whitespace.
    name = critical_section_names[wait_time_ranks[-i]]
    if isinstance(name, (bytes, np.bytes_)):
        name = name.decode('utf-8', errors='replace')
    name = re.sub(r"'(\S+)\s*'", r'\1', name.strip()) if hasattr(name, 'strip') else str(name)
    print(f" -> {wait_time_percentage:8.3f}% : {name}")
