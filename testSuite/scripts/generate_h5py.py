#!/usr/bin/env python3
import h5py

# Generate an HDF5 file using h5py which contains a single string attribute. This is used to test reading of variable length
# strings as written by h5py.
f                          = h5py.File("testSuite/outputs/h5py.hdf5","w")
f.attrs["stringAttribute"] = "this is a variable length string"

