#!/usr/bin/env python3
"""Verify the units (compound type) attribute written by tests.IO.HDF5.

Reads the `unitsAttribute` compound attribute from the group `myGroup` in the
file written by tests.IO.HDF5 and checks that each field round-tripped
correctly. Exits with status 0 on success and 1 on any mismatch, so that the
Fortran test can assert on the exit status.
"""
import sys
import h5py

fileName = "testSuite/outputs/test.IO.HDF5.hdf5"

try:
    with h5py.File(fileName, "r") as fileObject:
        attribute = fileObject["myGroup"].attrs["unitsAttribute"]
except (OSError, KeyError) as error:
    sys.stderr.write(f"verify_units.py: unable to read units attribute: {error}\n")
    sys.exit(1)

# The attribute is a scalar compound; extract each named field.
unitsInSI   = float(attribute["unitsInSI"])
description = attribute["description"]
quantity   = attribute["quantity"]
isComoving  = int(attribute["isComoving"])

# Fixed-length string fields come back as bytes; decode and strip the null terminator and any padding.
if isinstance(description, bytes):
    description = description.decode("ascii")
if isinstance(quantity, bytes):
    quantity = quantity.decode("ascii")
description = description.split("\x00", 1)[0]
quantity   = quantity  .split("\x00", 1)[0]

expectedUnitsInSI = 3.0856775814913673e22
failures = []
if abs(unitsInSI - expectedUnitsInSI) > 1.0e-6 * expectedUnitsInSI:
    failures.append(f"unitsInSI: got {unitsInSI!r}, expected {expectedUnitsInSI!r}")
if description != "megaparsec":
    failures.append(f"description: got {description!r}, expected 'megaparsec'")
if quantity != "length":
    failures.append(f"quantity: got {quantity!r}, expected 'length'")
if isComoving != 1:
    failures.append(f"isComoving: got {isComoving!r}, expected 1")

if failures:
    for failure in failures:
        sys.stderr.write(f"verify_units.py: {failure}\n")
    sys.exit(1)

sys.exit(0)
