#!/bin/bash
# Batch 3: (a) deeper M/m_res probe using a 1e12 host (10x less memory per M/m),
#          (b) perf sampling profiles across resolution.
set -u
cd /home/user/perf

# 1e12-host series: M/m from 1e4 to 3.16e6 (checks host-mass invariance and
# extends the reach in M/m by 0.5 dex over the 1e13 series).
./run_one.sh base12_m1.0e8  1.0e8  4 massTree=1.0e12
./run_one.sh base12_m1.0e7  1.0e7  4 massTree=1.0e12
python3 collect.py base12_m1.0e8 base12_m1.0e7
./run_one.sh base12_m3.16e6 3.16e6 4 massTree=1.0e12
python3 collect.py base12_m3.16e6
./run_one.sh base12_m1.0e6  1.0e6  2 massTree=1.0e12
python3 collect.py base12_m1.0e6
./run_one.sh base12_m3.16e5 3.16e5 2 massTree=1.0e12
python3 collect.py base12_m3.16e5

# Defaults-series: evolver settings as the CODE defaults them
# (fractionTimestepSatelliteMinimum=0, backtrackToSatellites=false) rather than
# the reference-file mitigations - measures the exponent users get by default.
./run_one.sh def_m1.0e9  1.0e9  4 fractionTimestepSatelliteMinimum=0.0 backtrackToSatellites=false
./run_one.sh def_m3.16e8 3.16e8 4 fractionTimestepSatelliteMinimum=0.0 backtrackToSatellites=false
./run_one.sh def_m1.0e8  1.0e8  4 fractionTimestepSatelliteMinimum=0.0 backtrackToSatellites=false
python3 collect.py def_m1.0e9 def_m3.16e8 def_m1.0e8
./run_one.sh def_m3.16e7 3.16e7 4 fractionTimestepSatelliteMinimum=0.0 backtrackToSatellites=false
python3 collect.py def_m3.16e7

# Perf profiles (shares, not timing): resolutions spanning the transition.
./profile_one.sh prof_m1.0e9 1.0e9 4
./profile_one.sh prof_m1.0e8 1.0e8 2
./profile_one.sh prof12_m1.0e6 1.0e6 1 massTree=1.0e12

echo BATCH3_DONE
