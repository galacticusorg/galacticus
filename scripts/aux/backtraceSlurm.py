#!/usr/bin/env python3
import sys
import os
import re
import argparse
import subprocess

# Obtain stack traces from all threads of all running Galacticus processes that are part of a given Slurm job and output to file.
# Andrew Benson (14-May-2025)

parser = argparse.ArgumentParser(prog='backtraceSlurm.py', description='Obtain stack traces from all threads of all running Galacticus processes in a Slurm job.')
parser.add_argument('slurmJobID', help='Slurm job ID')
parser.add_argument('outputFile', help='Output file for stack traces')
args = parser.parse_args()

# Parse the node list.
print("Parsing node list....")
node_list = None
result = subprocess.run(['scontrol', 'show', 'job', args.slurmJobID], capture_output=True, text=True)
for line in result.stdout.splitlines():
    match = re.match(r'^\s*NodeList=(.+)', line)
    if match:
        node_list = match.group(1)
        break

if node_list is None:
    print("unable to find NodeList", file=sys.stderr)
    sys.exit(1)

# Extract canonical node names.
print("Extracting canonical node names....")
node_names = []
result = subprocess.run(['scontrol', 'show', 'hostnames', node_list], capture_output=True, text=True)
for line in result.stdout.splitlines():
    line = line.strip()
    if line:
        node_names.append(line)

# Open output file and iterate over nodes.
user = os.environ.get('USER', '')
with open(args.outputFile, 'w') as output:
    for node_name in node_names:
        print(f"Connecting to node '{node_name}'....")
        output.write(f"Node: {node_name}\n")
        cmd = (
            f"srun --jobid={args.slurmJobID} -w {node_name} --pty bash -c "
            f"'pgrep -u {user} \"Galacticus.exe\" | xargs -i{{}} gdb --pid {{}} -batch "
            f"-ex \"info threads\" -ex \"thread apply all where\"'"
        )
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        output.write(result.stdout)
        # Reset the terminal as, for some reason, the above srun command seems to leave it in a weird state.
        subprocess.run(['reset'])
        output.write("\n\n")
