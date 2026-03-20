#!/usr/bin/env python3
import subprocess
import sys
import os
import re
import argparse
import h5py
import numpy as np

# Run a set of Galacticus models to test the state store/retrieve functionality under MPI.
# Andrew Benson (ported to Python)

# Parse command line options.
parser = argparse.ArgumentParser()
parser.add_argument("--processesPerNode", type=int, default=1   )
parser.add_argument("--allowRunAsRoot"  , type=str, default="no")
args, _ = parser.parse_known_args()

# We need at least 8 processes to run this test.
if args.processesPerNode < 8:
    print("SKIPPED: at least 8 processes per node are required for this test")
    sys.exit(0)

allowRunAsRoot = " --allow-run-as-root" if args.allowRunAsRoot == "yes" else ""

# Allow two run-throughs.
for runPass in range(2):
    print(f"State store/restore test pass #{runPass}")

    # Run full store model.
    subprocess.run("rm -f outputs/state.state*:MPI* outputs/state.gsl.state*:MPI*", shell=True)
    status = subprocess.run(
        f"export OMP_NUM_THREADS=1; cd ..; mpirun --oversubscribe -np 8{allowRunAsRoot} Galacticus.exe testSuite/parameters/state/store.xml",
        shell=True
    )
    if status.returncode != 0:
        print("FAILED: failed to run store model")
        sys.exit(0)

    # Find which MPI process ran the final tree.
    finalTreeProcessMPI = None
    for fileName in os.listdir("outputs"):
        m = re.match(r'state\.state\.log:MPI(\d+)', fileName)
        if m:
            processMPI = m.group(1)
            with open(f"outputs/{fileName}") as stateLogFile:
                for line in stateLogFile:
                    m2 = re.match(r'\s*Storing state for tree #(\d+)', line)
                    if m2 and int(m2.group(1)) == 15:
                        finalTreeProcessMPI = processMPI

    if finalTreeProcessMPI is not None:
        print(f"Final tree was run on MPI process {finalTreeProcessMPI}")
        if finalTreeProcessMPI != "0000":
            subprocess.run(f"cp -f outputs/state.state:MPI{finalTreeProcessMPI} outputs/state.state:MPI0000", shell=True)
            subprocess.run(f"cp -f outputs/state.gsl.state:MPI{finalTreeProcessMPI} outputs/state.gsl.state:MPI0000", shell=True)
    else:
        print("FAILED: failed to identify which thread/process ran final tree")
        sys.exit(0)

    # Run the restore model.
    status = subprocess.run(
        f"export OMP_NUM_THREADS=1; cd ..; mpirun --oversubscribe -np 1{allowRunAsRoot} Galacticus.exe testSuite/parameters/state/retrieve.xml",
        shell=True
    )
    if status.returncode != 0:
        print("FAILED: failed to run retrieve model")
        sys.exit(0)

    if not os.path.exists(f"outputs/stateStore:MPI{finalTreeProcessMPI}.hdf5"):
        print(f"FAILED: stateStore:MPI{finalTreeProcessMPI}.hdf5 file is missing")
        sys.exit(0)
    if not os.path.exists("outputs/stateRetrieve:MPI0000.hdf5"):
        print("FAILED: stateRetrieve:MPI0000.hdf5 file is missing")
        sys.exit(0)

    with h5py.File(f"outputs/stateStore:MPI{finalTreeProcessMPI}.hdf5", "r") as store, \
         h5py.File("outputs/stateRetrieve:MPI0000.hdf5", "r") as retrieve:

        storeOutput    = store["Outputs/Output1"]
        retrieveOutput = retrieve["Outputs/Output1"]
        storeData      = storeOutput["nodeData"]
        retrieveData   = retrieveOutput["nodeData"]

        storeTreeIndex = storeOutput["mergerTreeIndex"][:]
        treeFinal      = np.where(storeTreeIndex == 15)[0]
        if len(treeFinal) != 1:
            print("FAILED: unable to (uniquely) identify final tree in stored model output")
            sys.exit(0)
        treeFinalIndex = treeFinal[0]

        storeTreeStart   = storeOutput["mergerTreeStartIndex"][treeFinalIndex]
        storeTreeSize    = storeOutput["mergerTreeCount"][treeFinalIndex]
        retrieveTreeSize = retrieveOutput["mergerTreeCount"][-1]

        failed        = False
        statusMessage = ""

        if storeTreeSize != retrieveTreeSize:
            statusMessage += "FAILED: number of nodes in output changed after state retrieve\n"
            failed         = True

        datasets = list(storeData.keys())
        for dataset in datasets:
            storeDataset    = storeData[dataset][storeTreeStart:storeTreeStart + storeTreeSize]
            retrieveDataset = retrieveData[dataset][:]
            equal           = np.all(storeDataset == retrieveDataset)
            if not equal:
                statusMessage += f"FAILED: dataset '{dataset}' changed after state retrieve\n"
                statusMessage += f"   before --> {storeDataset}\n"
                statusMessage += f"   after  --> {retrieveDataset}\n"
                failed         = True

    if failed:
        if runPass == 1:
            print(statusMessage)
    else:
        print("SUCCESS!")
        break
