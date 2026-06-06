#!/usr/bin/env python3
import sys
import re

# Postprocess output from the linker to remove irrelevant warnings.
# Andrew Benson (12-April-2024)

if len(sys.argv) != 1:
    print("Usage: postprocessLinker.py", file=sys.stderr)
    sys.exit(1)

status = 0
for line in sys.stdin:
    if (
        re.search(r"warning: the use of `mktemp' is dangerous, better use `mkstemp'", line)
        or re.search(r"warning: Using 'dlopen' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking", line)
        or re.search(r"warning: [^:]+: requires executable stack \(because the \.note\.GNU\-stack section is executable\)", line)
    ):
        continue
    sys.stdout.write(line)
    if re.match(r'^Error:', line):
        status = 1

sys.exit(status)
