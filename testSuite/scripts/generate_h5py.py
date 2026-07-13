#!/usr/bin/env python3
import h5py
import numpy as np

# Generate an HDF5 file using h5py which contains a string attribute. This is used to test reading of variable length
# strings as written by h5py.
f                          = h5py.File("testSuite/outputs/h5py.hdf5","w")
f.attrs["stringAttribute"] = "this is a variable length string"

# Add a table (a compound-type dataset following the H5TB conventions) used to test reading of table columns. The long integer
# column uses values exceeding the 32-bit range to confirm genuine 64-bit reads.
tableData                    = np.zeros(8,dtype=[("realColumn","<f4"),("integerColumn","<i4"),("integer8Column","<i8"),("characterColumn","S8")])
tableData["realColumn"     ] =         np.arange(1,9)+0.5
tableData["integerColumn"  ] =      10*np.arange(1,9)
tableData["integer8Column" ] = 2**32+  np.arange(1,9)
tableData["characterColumn"] = [f"rec{i}".encode("ascii") for i in range(1,9)]
table                        = f.create_dataset("myTable",data=tableData)
table.attrs.create("CLASS"  ,np.bytes_("TABLE"     ))
table.attrs.create("VERSION",np.bytes_("3.0"       ))
table.attrs.create("TITLE"  ,np.bytes_("Test table"))
for i, name in enumerate(tableData.dtype.names):
    table.attrs.create(f"FIELD_{i}_NAME",np.bytes_(name))
f.close()
