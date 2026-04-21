#!/usr/bin/env python3
# Format module use sections of Fortran source files.
# Andrew Benson (ported to Python 2026)
#
# Mirrors scripts/aux/formatModuleUses.pl

import sys
import os
import shutil
import argparse

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.SourceTree import parse_file, walk_tree, serialize
from Galacticus.Build.SourceTree.Parse.ModuleUses import update_uses

parser = argparse.ArgumentParser(description='Format module use sections of Fortran source files')
parser.add_argument('infile', help='Input Fortran source file')
parser.add_argument('--suffix', default='~', help='Suffix for backup file (default: ~)')
args = parser.parse_args()

tree = parse_file(args.infile)
for node in walk_tree(tree):
    if node.get('type') == 'moduleUse':
        update_uses(node)

shutil.move(args.infile, args.infile + args.suffix)
with open(args.infile, 'w') as fh:
    fh.write(serialize(tree['firstChild']))
