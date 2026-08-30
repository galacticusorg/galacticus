#!/bin/bash
# Run one Galacticus model for the scaling experiment.
# Usage: run_one.sh <runName> <massResolution> <treeCount> [key=value ...]
set -u
name="$1"; mres="$2"; ntree="$3"; shift 3

export GALACTICUS_EXEC_PATH=/home/user/galacticus
export GALACTICUS_DATA_PATH=/home/user/datasets
export GALACTICUS_TOOLS_PATH=/home/user/tools
export OMP_NUM_THREADS=1

cd /home/user/perf
param="params/${name}.xml"
hdf5="out/${name}.hdf5"
log="logs/${name}.log"

python3 make_params.py "$param" "$mres" "$ntree" "$hdf5" "$@" || exit 1
rm -f "$hdf5"
/usr/bin/time -v /home/user/galacticus/Galacticus.exe "$param" > "$log" 2>&1
status=$?
echo "exit=$status" >> "$log"
if [ $status -ne 0 ] || grep -q -e FATAL -e "Aborted" "$log"; then
    echo "RUN ${name}: FAILED (exit=$status)"
    tail -5 "$log"
else
    summary=$(grep 'Run-time summary' "$log" | tail -1)
    walltime=$(grep 'Elapsed (wall clock)' "$log" | awk '{print $NF}')
    maxrss=$(grep 'Maximum resident set' "$log" | awk '{print $NF}')
    echo "RUN ${name}: OK wall=$walltime maxRSS=${maxrss}kB :: $summary"
fi
