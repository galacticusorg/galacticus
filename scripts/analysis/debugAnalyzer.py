#!/usr/bin/env python3
# Extract information from debug logs and export to CSV for further analysis.
# Andrew Benson (22-August-2023)

import re
import sys
from collections import defaultdict

# Get the name of the debug log file.
if len(sys.argv) != 2:
    sys.exit("Usage: debugAnalyzer.py <debugFile>")
debug_file_name = sys.argv[1]

# Initialize step counter and data.
i = -1
times = []
# properties[name]['value'][step] = value
# properties[name]['rate'][function][step] = rate
properties = defaultdict(lambda: {'value': {}, 'rate': defaultdict(dict)})

# Open and parse the debug log file.
re_step  = re.compile(r'^\s*step:\s+([0-9\.eE\+\-]+)')
re_value = re.compile(r'^\s*value:\s+([a-zA-Z:]+)\s+([\d\.eE\+\-]+)')
re_rate  = re.compile(r'^\s*rate:\s+\(([a-zA-Z0-9_]+)\)\s+([a-zA-Z:\d\[\]]+)\s+([\d\.eE\+\-]+)')

try:
    with open(debug_file_name) as debug_log:
        for line in debug_log:
            m = re_step.match(line)
            if m:
                i += 1
                times.append(m.group(1))
            else:
                m = re_value.match(line)
                if m:
                    properties[m.group(1)]['value'][i] = m.group(2)
                else:
                    m = re_rate.match(line)
                    if m:
                        properties[m.group(2)]['rate'][m.group(1)][i] = m.group(3)
                    else:
                        print(line, end='')
                        sys.exit("failed to parse line")
except OSError as e:
    sys.exit(f"failed to open {debug_file_name}: {e.strerror or e}")

step_count = i + 1

# Output values and rates.
prop_names = sorted(properties.keys())
rates = [
    {'property': prop, 'function': fn}
    for prop in prop_names
    for fn in sorted(properties[prop]['rate'].keys())
]

MISSING = "0.000000E+00"

with open("debugValues.csv", "w") as value_log:
    value_log.write("time , " + " , ".join(prop_names) + "\n")
    for j in range(step_count):
        values = [properties[p]['value'].get(j, MISSING) for p in prop_names]
        value_log.write(times[j] + " , " + " , ".join(values) + "\n")

with open("debugRates.csv", "w") as rates_log:
    headers = [r['property'] + ":" + r['function'] for r in rates]
    rates_log.write("time , " + " , ".join(headers) + "\n")
    for j in range(step_count):
        row = [properties[r['property']]['rate'][r['function']].get(j, MISSING) for r in rates]
        rates_log.write(times[j] + " , " + " , ".join(row) + "\n")
