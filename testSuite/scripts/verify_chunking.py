#!/usr/bin/env python3
"""Verify the chunk shapes chosen for large datasets written by tests.IO.HDF5.

Datasets too large to be held in a single HDF5 chunk have their chunk
dimensions reduced until the chunk fits within HDF5's 4GB limit. The reduced
chunk must still tile the extent of the dataset - if it does not, HDF5
allocates an extra chunk along the affected dimension of which only a small
fraction is ever written, wasting up to ~50% of the allocated space per
dimension (see https://github.com/galacticusorg/galacticus/issues/1365). The
limit must also be evaluated using the size of the actual datatype, not an
assumed 8 bytes.

HDF5 stores dimensions in C order, so h5py sees each Fortran array transposed,
with the Fortran first dimension varying fastest. Exits with status 0 on
success and 1 on any mismatch, so that the Fortran test can assert on the exit
status.
"""
import sys
import h5py

fileName = "testSuite/outputs/test.IO.HDF5.hdf5"

# Maximum chunk size (in bytes) permitted by HDF5.
chunkSizeMaximum = 4294967296

# Expected shapes and chunk shapes, in h5py (C) order, i.e. reversed relative to the Fortran extents used to create them.
expected = {
    # Fortran (600,100,100,100) doubles: 4.80GB, halved once along the first dimension to give a chunk of 2.40GB.
    "bigDataset"       : {"shape": (100, 100, 100, 600), "chunks": (100, 100, 100, 300)},
    # Fortran (601,100,100,100) doubles: 4.81GB, halved once along the first dimension. Rounding down would give a chunk of
    # 300, which tiles 601 only in three chunks (900 elements allocated for 601 held).
    "bigDatasetOdd"    : {"shape": (100, 100, 100, 601), "chunks": (100, 100, 100, 301)},
    # Fortran (601,201,100,100) doubles: 9.66GB, requiring a reduction along each of the first two dimensions.
    "bigDatasetOdd2D"  : {"shape": (100, 100, 201, 601), "chunks": (100, 100, 101, 301)},
    # Fortran (600,100,100,100) 4-byte integers: 2.40GB, which fits within a single chunk and so needs no reduction at all.
    "bigDatasetInteger": {"shape": (100, 100, 100, 600), "chunks": (100, 100, 100, 600)},
}

failures = []
try:
    with h5py.File(fileName, "r") as fileObject:
        for name, expectedValue in expected.items():
            dataset = fileObject[name]
            if dataset.shape != expectedValue["shape"]:
                failures.append(f"{name}: shape: got {dataset.shape!r}, expected {expectedValue['shape']!r}")
                continue
            if dataset.chunks is None:
                failures.append(f"{name}: dataset is not chunked")
                continue
            if dataset.chunks != expectedValue["chunks"]:
                failures.append(f"{name}: chunks: got {dataset.chunks!r}, expected {expectedValue['chunks']!r}")
                continue
            # Check that the chunk is within HDF5's maximum chunk size, as measured using the actual datatype.
            chunkSize = dataset.dtype.itemsize
            for chunkDimension in dataset.chunks:
                chunkSize *= chunkDimension
            if chunkSize > chunkSizeMaximum:
                failures.append(f"{name}: chunk of {chunkSize} bytes exceeds the maximum of {chunkSizeMaximum} bytes")
                continue
            # Check that the chunks tile the extent with no more than one chunk's worth of overshoot in total - i.e. that the
            # space HDF5 would allocate exceeds that needed by less than the size of a single chunk.
            elementsAllocated = 1
            elementsHeld = 1
            for extent, chunkDimension in zip(dataset.shape, dataset.chunks):
                elementsAllocated *= -(-extent // chunkDimension) * chunkDimension
                elementsHeld *= extent
            elementsChunk = 1
            for chunkDimension in dataset.chunks:
                elementsChunk *= chunkDimension
            if elementsAllocated - elementsHeld >= elementsChunk:
                failures.append(
                    f"{name}: chunks {dataset.chunks!r} do not tile extent {dataset.shape!r}:"
                    f" {elementsAllocated} elements allocated for {elementsHeld} held"
                )
except (OSError, KeyError) as error:
    sys.stderr.write(f"verify_chunking.py: unable to read dataset: {error}\n")
    sys.exit(1)

if failures:
    for failure in failures:
        sys.stderr.write(f"verify_chunking.py: mismatch: {failure}\n")
    sys.exit(1)
sys.exit(0)
