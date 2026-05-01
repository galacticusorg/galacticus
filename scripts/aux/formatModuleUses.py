#!/usr/bin/env python3
# Format module use sections of Fortran source files.
# Andrew Benson (ported to Python 2026)
#
# Mirrors scripts/aux/formatModuleUses.pl

import sys
import os
import tempfile
import argparse


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

content  = serialize(tree['firstChild'])
dirpath  = os.path.dirname(os.path.abspath(args.infile))
tmp_path = None
try:
    fd, tmp_path = tempfile.mkstemp(dir=dirpath)
    with os.fdopen(fd, 'w') as fh:
        fh.write(content)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(args.infile, args.infile + args.suffix)
    os.replace(tmp_path, args.infile)
    tmp_path = None
except Exception:
    if tmp_path is not None:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
    raise
