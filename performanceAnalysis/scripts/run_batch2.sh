#!/bin/bash
# Batch 2: extend sweep deeper, then diagnostics and decomposition variants.
set -u
cd /home/user/perf

# Extend the base sweep to higher resolution (the break was not yet reached at 1e8).
./run_one.sh base_m3.16e7 3.16e7 4
python3 collect.py base_m3.16e7
./run_one.sh base_m1.0e7 1.0e7 4
python3 collect.py base_m1.0e7

# Step-size / limiter diagnostics via the evolve profiler (not for timing).
./run_one.sh diag_m1.0e9  1.0e9  4 profiler=simple profileOdeEvolver=true
./run_one.sh diag_m1.0e8  1.0e8  4 profiler=simple profileOdeEvolver=true
./run_one.sh diag_m3.16e7 3.16e7 4 profiler=simple profileOdeEvolver=true

# Census runs (multi-redshift standard output; not for timing).
./run_one.sh census_m1.0e9 1.0e9 4 outputRedshifts=4.0:3.0:2.0:1.5:1.0:0.7:0.5:0.3:0.15:0.0
./run_one.sh census_m1.0e8 1.0e8 4 outputRedshifts=4.0:3.0:2.0:1.5:1.0:0.7:0.5:0.3:0.15:0.0

# Mechanism-decomposition variants at m_res=1e8.
./run_one.sh varNoDefer_m1.0e8     1.0e8 4 fractionTimestepSatelliteMinimum=0.0
./run_one.sh varNoBacktrack_m1.0e8 1.0e8 4 backtrackToSatellites=false
./run_one.sh varReuseStep_m1.0e8   1.0e8 4 reuseODEStepSize=true
./run_one.sh varSyncLoose_m1.0e8   1.0e8 4 timeOffsetMaximumAbsolute=0.1 timeOffsetMaximumRelative=0.01
./run_one.sh varSyncTight_m1.0e8   1.0e8 4 timeOffsetMaximumAbsolute=0.001 timeOffsetMaximumRelative=0.0001

python3 collect.py diag_m1.0e9 diag_m1.0e8 diag_m3.16e7 census_m1.0e9 census_m1.0e8 \
    varNoDefer_m1.0e8 varNoBacktrack_m1.0e8 varReuseStep_m1.0e8 varSyncLoose_m1.0e8 varSyncTight_m1.0e8
echo BATCH2_DONE
