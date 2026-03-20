#!/usr/bin/env python3
import subprocess
import sys

# Test Perl modules.
# Andrew Benson (ported to Python)

# Check individual compilation of .pl and .pm files.
subprocess.run("cd ..; find . -type f -name '*.pl' | xargs -r -n 1 perl -c 1>perl.tmp 2>&1; sed -r 's/failed/FAILED/g' perl.tmp; rm perl.tmp", shell=True)
subprocess.run("cd ..; find perl -type f -name '*.pm' | xargs -r -n 1 perl -c 1>perl.tmp 2>&1; sed -r 's/failed/FAILED/g' perl.tmp; rm perl.tmp", shell=True)
