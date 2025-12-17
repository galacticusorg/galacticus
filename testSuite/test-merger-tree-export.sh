#!/usr/bin/env bash

# Import and then export merger trees, verifying consistency between the exported and original tree files.
# Andrew Benson (16-December-2025)

status="SUCCESS"
for tree in testTreeExportBranchJumps                         \
	    testTreeExportInitialSubhalo                      \
	    testTreeExportSubhaloPromotions                   \
	    testTreeExportSubhaloSubhaloMergers               \
	    testTreeExportSubSubhaloToSubhalo                 \
	    testTreeExportSubhaloHostIsClone                  \
	    testTreeExportDoubleBranchJumpAndMergeWithClone   \
	    testTreeExportInitialSubhaloWithInitialBranchJump \
	    testTreeExportDeepHierarchyBranchJump             \
	    testTreeExportMisorderedBranchJumps; do
    echo Exporting tree from file: ${tree}
    mkdir -p outputs
    cp data/mergerTrees/${tree}.hdf5 outputs/testTreeExportOriginal.hdf5
    rm -f outputs/testTreeExportExported.hdf5
    pushd ..
    ./Galacticus.exe testSuite/parameters/testTreeExport.xml
    if [ $? -eq 0 ]; then
	./scripts/aux/mergerTreeExportVerify.py testSuite/outputs/testTreeExportOriginal.hdf5 testSuite/outputs/testTreeExportExported.hdf5
	if [ $? -ne 0 ]; then
	    status="FAILED {verify}"
	fi
    else
	status="FAILED {model run}"
    fi
    popd
done
echo ${status}

exit 0
