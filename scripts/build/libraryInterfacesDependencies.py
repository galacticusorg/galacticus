#!/usr/bin/env python3
import os
import sys
import xml.etree.ElementTree as ET

# Build a Makefile with dependencies for libgalacticus.
# Andrew Benson (ported to Python 2026)

exec_path  = os.environ['GALACTICUS_EXEC_PATH']
build_path = os.environ['BUILDPATH'].rstrip('/') + '/'

# Read the list of functionClasses to compile into the library.
classes_xml = os.path.join(exec_path, "source", "libraryClasses.xml")
tree = ET.parse(classes_xml)
root = tree.getroot()

# Collect function class names, sorted for deterministic output.
# The XML structure is <libraryClasses><classes><name1/><name2/>...</classes></libraryClasses>
# or <libraryClasses><classes name1="" name2="..."/></libraryClasses> — handle both.
# Based on the Perl code: $libraryFunctionClasses->{'classes'} is iterated with hashList(..., keyAs => "name"),
# which adds a 'name' key equal to the hash key. So the XML has <classes> containing child elements
# whose tag names are the class names.
classes_el = root.find('classes')
if classes_el is None:
    print("libraryInterfacesDependencies.py: could not find <classes> in libraryClasses.xml",
          file=sys.stderr)
    sys.exit(1)

# Each child tag is a class name; sort for deterministic output.
class_names = sorted(child.tag for child in classes_el)

rules = []

for name in class_names:
    rules.append(
        f"{build_path}libgalacticus/{name}.F90:\n"
        f"\t./scripts/build/libraryInterfaces.pl\n"
        f"{build_path}libgalacticus/{name}.p.F90.up:"
        f" {build_path}libgalacticus/{name}.F90"
        f" {build_path}hdf5FCInterop.dat"
        f" {build_path}openMPCriticalSections.xml\n"
        f"\t./scripts/build/preprocess.pl"
        f" {build_path}libgalacticus/{name}.F90"
        f" {build_path}libgalacticus/{name}.p.F90\n"
        f"{build_path}libgalacticus/{name}.p.F90 :"
        f" {build_path}libgalacticus/{name}.p.F90.up\n"
        f"\t@true\n"
        f"{build_path}libgalacticus/{name}.o :"
        f" {build_path}libgalacticus/{name}.p.F90"
        f" {build_path}libgalacticus/{name}.d Makefile\n"
        f"\t@mkdir -p {build_path}/moduleBuild\n"
        f"\t$(FCCOMPILER) -c {build_path}libgalacticus/{name}.p.F90"
        f" -o {build_path}libgalacticus/{name}.o $(FCFLAGS) 2>&1"
        f" | ./scripts/build/postprocess.py {build_path}libgalacticus/{name}.p.F90\n"
        f"\n"
    )

# Rule for libgalacticus_classes.d.
dep_d_files = ' '.join(f"{build_path}libgalacticus/{n}.d" for n in class_names)
rules.append(f"{build_path}libgalacticus_classes.d: {dep_d_files}\n")
rules.append(f"\t@rm -f {build_path}libgalacticus_classes.d~\n")
for name in class_names:
    rules.append(f"\t@cat {build_path}libgalacticus/{name}.d >> {build_path}libgalacticus_classes.d~\n")
rules.append(
    f"\t@sort -u {build_path}libgalacticus_classes.d~ -o {build_path}libgalacticus_classes.d~\n"
    f"\t@if cmp -s {build_path}libgalacticus_classes.d {build_path}libgalacticus_classes.d~ ; then \\\n"
    f"\t rm {build_path}libgalacticus_classes.d~ ; \\\n"
    f"\telse \\\n"
    f"\t mv {build_path}libgalacticus_classes.d~ {build_path}libgalacticus_classes.d ; \\\n"
    f"\tfi\n"
)

# Rule for libgalacticus.so.
dep_o_files = ' '.join(f"{build_path}libgalacticus/{n}.o" for n in class_names)
rules.append(f"libgalacticus.so: {dep_o_files}\n")

with open(os.path.join(build_path, "Makefile_Library_Dependencies"), 'w') as fh:
    fh.write(''.join(rules))
