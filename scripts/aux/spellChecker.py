#!/usr/bin/env python3
import os
import sys
import latex_spellcheck

# Perform spell checking of files.
# Andrew Benson (28-February-2023) [Python port]

if len(sys.argv) != 3:
    print("Usage: spellChecker.py <fileName> <warningFile>", file=sys.stderr)
    raise SystemExit(1)

file_name         = sys.argv[1]
warning_file_name = sys.argv[2]

# Check the file.
warnings = latex_spellcheck.spell_check_file(file_name, file_name)

# Append any warnings to the warning file.
if warnings:
    with open(warning_file_name, 'a') as fh:
        fh.write(warnings)
