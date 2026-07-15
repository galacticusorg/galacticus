#!/usr/bin/env python3
"""State store/restore round-trip test for the recursive="yes" machinery
(issue #695).

This is the state-store analogue of ``test-stateRestore.py``, but the model
uses a percolation virial density contrast with the default
``virialDensityContrastDefinition`` halo scale, forming the bounded
construction cycle served by the generated shim. The shim overrides
``stateStore``/``stateRestore`` to forward to the real object under
construction, so this test checks that state storage and retrieval of a model
whose object graph contains a recursion shim round-trips exactly: storing the
state after the penultimate tree and restoring it to run the final tree must
reproduce the store run's final-tree output bit-for-bit.

Andrew Benson
"""

import subprocess
import sys
import os
import re
import h5py
import numpy as np

subprocess.run("mkdir -p outputs", shell=True)

# Allow two run-throughs. If tabulation files are updated during the "store" run
# it can result in divergent behavior in the "restore" run. This is a limitation
# of our current state store/restore infrastructure (see test-stateRestore.py).
for runPass in range(2):
    print(f"Recursive state store/restore test pass #{runPass}")

    # Run full store model.
    subprocess.run(
        "rm -f outputs/stateRecursive.state*:openMP* "
        "outputs/stateRecursive.gsl.state*:openMP*",
        shell=True,
    )
    status = subprocess.run(
        "export OMP_NUM_THREADS=12; cd ..; "
        "./Galacticus.exe testSuite/parameters/state/storeRecursive.xml",
        shell=True,
    )
    if status.returncode != 0:
        print("FAILED: failed to run store model")
        sys.exit(1)

    # Find which thread ran the final tree.
    finalTreeThread = None
    for fileName in os.listdir("outputs"):
        m = re.match(r'stateRecursive\.state\.log:openMP(\d+)', fileName)
        if m:
            thread = m.group(1)
            with open(f"outputs/{fileName}") as stateLogFile:
                for line in stateLogFile:
                    m2 = re.match(r'\s*Storing state for tree #(\d+)', line)
                    if m2 and int(m2.group(1)) == 15:
                        finalTreeThread = thread

    if finalTreeThread is not None:
        print(f"Final tree was run by thread {finalTreeThread}")
        subprocess.run(
            f"cp -f outputs/stateRecursive.state:openMP{finalTreeThread} "
            f"outputs/stateRecursive.state",
            shell=True,
        )
        subprocess.run(
            f"cp -f outputs/stateRecursive.gsl.state:openMP{finalTreeThread} "
            f"outputs/stateRecursive.gsl.state",
            shell=True,
        )
    else:
        print("FAILED: failed to identify which thread ran final tree")
        sys.exit(1)

    # Run the restore model.
    status = subprocess.run(
        "export OMP_NUM_THREADS=1; cd ..; "
        "./Galacticus.exe testSuite/parameters/state/retrieveRecursive.xml",
        shell=True,
    )
    if status.returncode != 0:
        print("FAILED: failed to run retrieve model")
        sys.exit(1)

    if not os.path.exists("outputs/stateRetrieveRecursive.hdf5"):
        print("FAILED: stateRetrieveRecursive.hdf5 file is missing")
        sys.exit(1)
    if not os.path.exists("outputs/stateStoreRecursive.hdf5"):
        print("FAILED: stateStoreRecursive.hdf5 file is missing")
        sys.exit(1)

    with h5py.File("outputs/stateStoreRecursive.hdf5", "r") as store, \
         h5py.File("outputs/stateRetrieveRecursive.hdf5", "r") as retrieve:

        storeOutput    = store["Outputs/Output1"]
        retrieveOutput = retrieve["Outputs/Output1"]
        storeData      = storeOutput["nodeData"]
        retrieveData   = retrieveOutput["nodeData"]

        # Find the final tree in the store model.
        storeTreeIndex = storeOutput["mergerTreeIndex"][:]
        treeFinal      = np.where(storeTreeIndex == 15)[0]
        if len(treeFinal) != 1:
            print("FAILED: unable to (uniquely) identify final tree in stored "
                  "model output")
            sys.exit(1)
        treeFinalIndex = treeFinal[0]

        storeTreeStart   = storeOutput["mergerTreeStartIndex"][treeFinalIndex]
        storeTreeSize    = storeOutput["mergerTreeCount"][treeFinalIndex]
        retrieveTreeSize = retrieveOutput["mergerTreeCount"][-1]

        failed        = False
        statusMessage = ""

        if storeTreeSize != retrieveTreeSize:
            statusMessage += ("FAILED: number of nodes in output changed after "
                              "state retrieve\n")
            failed         = True

        if not failed:
            for dataset in list(storeData.keys()):
                storeDataset    = storeData[dataset][
                    storeTreeStart:storeTreeStart + storeTreeSize]
                retrieveDataset = retrieveData[dataset][:]
                if not np.all(storeDataset == retrieveDataset):
                    statusMessage += (f"FAILED: dataset '{dataset}' changed "
                                      f"after state retrieve\n")
                    statusMessage += f"   before --> {storeDataset}\n"
                    statusMessage += f"   after  --> {retrieveDataset}\n"
                    failed         = True

    if failed:
        if runPass == 1:
            print(statusMessage)
            print("FAILED: recursive state store/restore test")
            sys.exit(1)
    else:
        print("SUCCESS: recursive state store/restore test")
        break
