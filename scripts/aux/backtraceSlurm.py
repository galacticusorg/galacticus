#!/usr/bin/env python3
import sys
import re
import argparse
import subprocess

# Obtain stack traces from all threads of all running Galacticus processes that are part of a given Slurm job and output to file.
# Andrew Benson (14-May-2025)
# Job-based process selection and process identity reporting added with assistance from Claude (23-July-2026).

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

# Script to be run on each node. Processes are taken from the process list of the job itself (via `scontrol listpids`), rather than
# from a simple match on process name for the user - matching on process name alone will also find Galacticus processes belonging to
# any other job of the same user that happens to share the node, resulting in stack traces that do not correspond to the job being
# debugged. The identity of each process (command line, working directory, start time, and state) is recorded alongside its stack
# trace so that any such misattribution is immediately apparent.
remoteScript = r"""
found=0
while read -r pid rest; do
  case "$pid" in
    ""|PID) continue ;;
  esac
  comm=$(cat /proc/$pid/comm 2>/dev/null)
  if [ "$comm" != "Galacticus.exe" ]; then continue; fi
  found=$((found+1))
  echo "===================================================================="
  echo "PID              : $pid"
  echo "Command          : $(tr "\0" " " < /proc/$pid/cmdline)"
  echo "Working directory: $(readlink /proc/$pid/cwd)"
  echo "Started          : $(ps -o lstart= -p $pid)"
  echo "Elapsed          : $(ps -o etime= -p $pid)"
  echo "State            : $(ps -o stat= -p $pid)"
  echo "===================================================================="
  gdb --pid $pid -batch -ex "info threads" -ex "thread apply all where"
done < <(scontrol listpids {JOBID})
if [ $found -eq 0 ]; then
  echo "WARNING: no Galacticus.exe processes belonging to job {JOBID} were found on this node"
fi
""".replace("{JOBID}", args.slurmJobID)

# Open output file and iterate over nodes.
with open(args.outputFile, 'w') as output:
    for node_name in node_names:
        print(f"Connecting to node '{node_name}'....")
        output.write(f"Node: {node_name}\n")
        # The `--overlap` option (Slurm ≥ 20.11) allows this step to share the resources already allocated to the job - without it
        # the step will block until resources become free, which will never happen for a job that is fully occupying its allocation.
        cmd = (
            f"srun --jobid={args.slurmJobID} --overlap -w {node_name} --pty bash -c "
            f"'{remoteScript}'"
        )
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        output.write(result.stdout)
        # Reset the terminal as, for some reason, the above srun command seems to leave it in a weird state.
        subprocess.run(['reset'])
        output.write("\n\n")
