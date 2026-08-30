#!/bin/bash
# Profile one Galacticus run with perf sampling (self-contained: generates
# params, runs under perf record, emits flat + call-graph reports).
# Usage: profile_one.sh <runName> <massResolution> <treeCount> [key=value ...]
set -u
name="$1"; mres="$2"; ntree="$3"; shift 3
PERF=/usr/lib/linux-tools/6.8.0-138-generic/perf

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
$PERF record -F 297 -g -o "out/${name}.perf.data" -- \
    /home/user/galacticus/Galacticus.exe "$param" > "$log" 2>&1
echo "exit=$?" >> "$log"
$PERF report -i "out/${name}.perf.data" --stdio --no-children --percent-limit 0.5 \
    > "out/${name}.perf.flat.txt" 2>/dev/null
$PERF report -i "out/${name}.perf.data" --stdio --children --percent-limit 1.0 -g none \
    > "out/${name}.perf.cum.txt" 2>/dev/null
echo "PROFILE ${name}: done"
grep 'Run-time summary' "$log" | tail -1
