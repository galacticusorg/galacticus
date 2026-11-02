#!/usr/bin/env python3
import subprocess
import sys
import os
import argparse
import h5py
import numpy as np

# Test that the Farahi excursion set first crossing solvers give the same tables when their work is shared between MPI
# processes as when run on a single process. The models span several epochs, so that the tables are extended in time, carrying
# over solutions already found.
# Andrew Benson, Claude

# Parse command line options.
parser = argparse.ArgumentParser()
parser.add_argument("--processesPerNode", type=int, default=1   )
parser.add_argument("--allowRunAsRoot"  , type=str, default="no")
args, _ = parser.parse_known_args()

# We need at least 3 processes to run this test.
countProcesses = 3
if args.processesPerNode < countProcesses:
    print(f"SKIPPED: at least {countProcesses} processes per node are required for this test")
    sys.exit(0)

allowRunAsRoot = " --allow-run-as-root" if args.allowRunAsRoot == "yes" else ""

# Run each class on one process, and on several. The first run, on a single process, is a warm-up whose results are discarded
# (they are overwritten by the second run): it ensures that the cached tabulations (of sigma(M) and spherical collapse) exist on
# file, so that every compared run reads the same tabulations. A run which builds those tabulations itself can differ from one
# which reads them by ~1e-7, because the result depends on how the tabulations were extended during the run.
classes       = [ "farahi", "farahiMidpoint" ]
probabilities = {}
for className in classes:
    probabilities[className] = {}
    for count in [ 1, 1, countProcesses ]:
        subprocess.run(f"rm -f outputs/excursionSetsFarahiMPI_{className}:MPI*.hdf5", shell=True)
        status = subprocess.run(
            f"export OMP_NUM_THREADS=1; cd ..; mpirun --oversubscribe -np {count}{allowRunAsRoot} Galacticus.exe testSuite/parameters/excursionSetsFarahiMPI/{className}.xml",
            shell=True
        )
        if status.returncode != 0:
            print(f"FAILED: failed to run model for class '{className}' on {count} process(es)")
            sys.exit(0)
        probabilities[className][count] = []
        for rank in range(count):
            fileName = f"outputs/excursionSetsFarahiMPI_{className}:MPI{rank:04d}.hdf5"
            if not os.path.exists(fileName):
                print(f"FAILED: output file '{fileName}' is missing")
                sys.exit(0)
            with h5py.File(fileName, "r") as model:
                probabilities[className][count].append(model["excursionSets/firstCrossingProbability"][:])

failed = False

# Every process must hold the same table as the single process run.
for className in classes:
    reference = probabilities[className][1][0]
    for rank, probability in enumerate(probabilities[className][countProcesses]):
        if not np.allclose(probability, reference, rtol=1.0e-12, atol=0.0):
            print(f"FAILED: class '{className}' first crossing probabilities on process {rank} of {countProcesses} differ from those on a single process")
            print(f"   single process --> {reference}")
            print(f"   process {rank}      --> {probability}")
            failed = True

# The two classes solve the same problem, differing only in their discretization, so must agree closely. This catches errors
# common to every number of processes.
if not np.allclose(probabilities["farahiMidpoint"][1][0], probabilities["farahi"][1][0], rtol=5.0e-2, atol=0.0):
    print("FAILED: 'farahi' and 'farahiMidpoint' first crossing probabilities differ")
    print(f"   farahi         --> {probabilities['farahi'        ][1][0]}")
    print(f"   farahiMidpoint --> {probabilities['farahiMidpoint'][1][0]}")
    failed = True

if not failed:
    print("SUCCESS: Farahi first crossing probabilities are independent of the number of MPI processes")
