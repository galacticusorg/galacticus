#!/usr/bin/env python3
"""Verify write-only objects written by tests.IO.HDF5.

Reads objects written by tests.IO.HDF5 for which no Fortran reader exists and
checks their shapes and values: the 2-D double attribute, the 3-D integer
dataset, and the 4-D long integer dataset. HDF5 stores dimensions in C order,
so h5py sees each Fortran array transposed, with the Fortran first dimension
varying fastest. Exits with status 0 on success and 1 on any mismatch, so
that the Fortran test can assert on the exit status.
"""
import sys
import h5py
import numpy as np

fileName = "testSuite/outputs/test.IO.HDF5.hdf5"

# Expected values matching those written by tests.IO.HDF5.
expected = {
    # reshape([1..6],[3,2]) -> h5py shape (2,3).
    "doubleAttribute2dArray": np.arange(1.0, 7.0).reshape(2, 3),
    # reshape([(2*i+1,i=1,60)],[5,4,3]) -> h5py shape (3,4,5).
    "integerDataset3dArray": 2 * np.arange(1, 61).reshape(3, 4, 5) + 1,
    # reshape([(2**32+7*i,i=1,120)],[3,4,2,5]) -> h5py shape (5,2,4,3).
    "integer8Dataset4dArray": 2**32 + 7 * np.arange(1, 121).reshape(5, 2, 4, 3),
}

failures = []
try:
    with h5py.File(fileName, "r") as fileObject:
        group = fileObject["myGroup"]
        values = {
            "doubleAttribute2dArray": group.attrs["doubleAttribute2dArray"],
            "integerDataset3dArray": group["integerDataset3dArray"][:],
            "integer8Dataset4dArray": group["integer8Dataset4dArray"][:],
        }
except (OSError, KeyError) as error:
    sys.stderr.write(f"verify_attributes.py: unable to read object: {error}\n")
    sys.exit(1)

for name, expectedValue in expected.items():
    value = values[name]
    if value.shape != expectedValue.shape:
        failures.append(f"{name}: shape: got {value.shape!r}, expected {expectedValue.shape!r}")
    elif not np.array_equal(value, expectedValue):
        failures.append(f"{name}: values: got {value!r}, expected {expectedValue!r}")

if failures:
    for failure in failures:
        sys.stderr.write(f"verify_attributes.py: mismatch: {failure}\n")
    sys.exit(1)
sys.exit(0)
