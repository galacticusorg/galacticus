#!/usr/bin/env python3
import subprocess
import os
import json
import numpy as np

# Run models to benchmark performance of a Milky Way model.
# Andrew Benson (10-August-2022; ported to Python)

# Make output directory.
os.makedirs("outputs", exist_ok=True)

# Run the benchmark model multiple times.
run_times = []
for i in range(12):
    # Run the model.
    result = subprocess.run(
        "cd ..; /usr/bin/time --format=\"%e\" --output=testSuite/outputs/benchmark_milkyWay.log ./Galacticus.exe testSuite/parameters/benchmark_milkyWay.xml",
        shell=True
    )
    if result.returncode != 0:
        print("FAIL: Milky Way benchmark model failed to run")
        raise SystemExit(1)
    # Extract timing data (skip first two runs as warm-up).
    if i >= 2:
        with open("outputs/benchmark_milkyWay.log") as f:
            run_times.append(float(f.read().strip()))

# Find average and standard deviation of run times.
run_times = np.array(run_times)
run_time_average        = run_times.mean()
run_time_standard_error = run_times.std() / np.sqrt(len(run_times))
print(f"Timings: {run_times}")
print(f"Benchmark results: {run_time_average} ± {run_time_standard_error} s")

# Generate JSON report.
output = [
    {
        "name":  "Milky Way model - Wall Time",
        "unit":  "seconds",
        "value": float(run_time_average),
        "range": float(run_time_standard_error)
    }
]
with open("outputs/benchmark_milkyWayModel.json", "w") as f:
    json.dump(output, f, indent=4)

print("SUCCESS: Milky Way benchmark model")
