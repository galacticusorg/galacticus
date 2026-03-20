#!/usr/bin/env python3
import subprocess
import sys
import os
import re
import h5py
import numpy as np

# Run a set of Galacticus models to test the state store/retrieve functionality.
# Andrew Benson (ported to Python)

# Allow two run-throughs. If tabulation files are updated during the "store" run it can result in divergent behavior in the
# "restore" run. This is a limitation of our current state store/restore infrastructure.
for runPass in range(2):
    print(f"State store/restore test pass #{runPass}")

    # Run full store model.
    subprocess.run("rm -f outputs/state.state*:openMP* outputs/state.gsl.state*:openMP*", shell=True)
    status = subprocess.run(
        "export OMP_NUM_THREADS=12; cd ..; ./Galacticus.exe testSuite/parameters/state/store.xml",
        shell=True
    )
    if status.returncode != 0:
        print("FAILED: failed to run store model")
        sys.exit(0)

    # Find which thread ran the final tree.
    finalTreeThread = None
    for fileName in os.listdir("outputs"):
        m = re.match(r'state\.state\.log:openMP(\d+)', fileName)
        if m:
            thread = m.group(1)
            with open(f"outputs/{fileName}") as stateLogFile:
                for line in stateLogFile:
                    m2 = re.match(r'\s*Storing state for tree #(\d+)', line)
                    if m2 and int(m2.group(1)) == 15:
                        finalTreeThread = thread

    if finalTreeThread is not None:
        print(f"Final tree was run by thread {finalTreeThread}")
        subprocess.run(f"cp -f outputs/state.state:openMP{finalTreeThread} outputs/state.state", shell=True)
        subprocess.run(f"cp -f outputs/state.gsl.state:openMP{finalTreeThread} outputs/state.gsl.state", shell=True)
    else:
        print("FAILED: failed to identify which thread ran final tree")
        sys.exit(0)

    # Run the restore model.
    status = subprocess.run(
        "export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/state/retrieve.xml",
        shell=True
    )
    if status.returncode != 0:
        print("FAILED: failed to run retrieve model")
        sys.exit(0)

    if not os.path.exists("outputs/stateRetrieve.hdf5"):
        print("FAILED: stateRetrieve.hdf5 file is missing")
        sys.exit(0)
    if not os.path.exists("outputs/stateStore.hdf5"):
        print("FAILED: stateStore.hdf5 file is missing")
        sys.exit(0)

    with h5py.File("outputs/stateStore.hdf5", "r") as store, \
         h5py.File("outputs/stateRetrieve.hdf5", "r") as retrieve:

        storeOutput    = store["Outputs/Output1"]
        retrieveOutput = retrieve["Outputs/Output1"]
        storeData      = storeOutput["nodeData"]
        retrieveData   = retrieveOutput["nodeData"]

        # Find the tree in the store model.
        storeTreeIndex  = storeOutput["mergerTreeIndex"][:]
        treeFinal       = np.where(storeTreeIndex == 15)[0]
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
