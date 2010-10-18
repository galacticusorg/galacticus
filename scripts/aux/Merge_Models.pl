#!/usr/bin/env perl
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::IO::HDF5::Group;
use Data::Dumper;

# Merges Galacticus models.
# Andrew Benson (8-July-2010)

if ($#ARGV < 1) {die("Usage: Merge_Models.pl <model1> ...... <outputModel>")};
@mergeFileNames = @ARGV[0..$#ARGV-1];
$outputFileName = $ARGV[$#ARGV];

# Open all files.
print "-> Opening output file: ".$outputFileName."\n";
$outputFile = new PDL::IO::HDF5(">".$outputFileName);
foreach $mergeFileName ( @mergeFileNames ) {
    print "-> Opening file for merge: ".$mergeFileName."\n";
    $mergeFiles[++$#mergeFiles] = new PDL::IO::HDF5($mergeFileName);
}

# Copy parameters and version groups from merge files to output file.
print "-> Copying Parameters and Version groups\n";
foreach $groupName ( "Parameters", "Version" ) {
    # Open the groups.
    $group    = $mergeFiles[0]->group($groupName);
    $newGroup = $outputFile   ->group($groupName);
    # Read and write the attributes.
    @attributes = $group->attrs();
    foreach $attribute ( @attributes ) {
	@attrValue = $group->attrGet($attribute);
	$newGroup->attrSet($attribute => $attrValue[0]);
    }
}

# Copy output group attributes.
$outputsGroup = $mergeFiles[0]->group("Outputs");
@outputGroups = $outputsGroup->groups;
# Loop over all outputs groups.
foreach $outputGroupName ( @outputGroups ) {
    $outputGroup    = $outputsGroup->group($outputGroupName);
    $newOutputGroup = $outputFile->group("Outputs")->group($outputGroupName);
    @attributes     = $outputGroup->attrs();
    foreach $attribute ( @attributes ) {
	@attrValue = $outputGroup->attrGet($attribute);
	$newOutputGroup->attrSet($attribute => $attrValue[0]);
    }
}

# Modify the parameters to reflect that this is now a complete model.
print "-> Modifying parameter attributes of base file\n";
$group = $outputFile->group("Parameters");
$group->attrSet(
		"treeEvolveWorkerNumber" => 1,
		"treeEvolveWorkerCount" => 1
		);

# Combine the datasets in the nodeData group.
print "-> Combining main datasets\n";
# Get the Outputs group from the output file.
$outputsGroup = $mergeFiles[0]->group("Outputs");
# Get a list of subgroups to search.
@outputGroups    = $outputsGroup->groups;
# Loop over all outputs groups.
foreach $outputGroup ( @outputGroups ) {
    # Get the node data group.
    $nodeDataGroup = $outputsGroup->group($outputGroup."/nodeData");
    $newNodeDataGroup = $outputFile->group("Outputs")->group($outputGroup)->group("nodeData");
    # Get a list of all datasets.
    @datasets = $nodeDataGroup->datasets;
    # Loop over datasets.
    foreach $dataset ( @datasets ) {
	$newDatasetValues = pdl [];
	# Loop over merge files.
	$iFile = -1;
	foreach $mergeFile ( @mergeFiles ) {
	    # Read the dataset from this file.
	    ++$iFile;
	    $outputsGroup  = $mergeFile    ->group("Outputs");
	    $nodeDataGroup = $outputsGroup ->group($outputGroup."/nodeData");
	    $thisDataset   = $nodeDataGroup->dataset($dataset);
	    $datasetValues = $thisDataset  ->get();
	    $offsets[$iFile] = nelem($datasetValues);
	    # Append to the combined dataset.
	    $newDatasetValues = $newDatasetValues->append($datasetValues);
	}
	# Create a new dataset in the output file.
	$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
						  parent  => $newNodeDataGroup,
						  fileObj => $outputFile
						  );
	# Write the new dataset back to the output file.
	$newDataset->set($newDatasetValues);
    }
}

# Combine tree datasets.
print "-> Combining tree datasets\n";
# Get the Outputs group from the output file.
$outputsGroup = $mergeFiles[0]->group("Outputs");
# Get a list of subgroups to search.
@outputGroups    = $outputsGroup->groups;
# Loop over all outputs groups.
foreach $outputGroup ( @outputGroups ) {
    # Get the new output group.
    $newOutputGroup = $outputFile->group("Outputs")->group($outputGroup);
    # Loop over the tree datasets.
    foreach $dataset ( "mergerTreeIndex", "mergerTreeCount", "mergerTreeStartIndex", "mergerTreeWeight" ) {
	# Reset values.
	$newDatasetValues = pdl [];
	# Loop over merge files.
	$iFile = -1;
	$offset = 0;
	foreach $mergeFile ( @mergeFiles ) {
	    # Read the dataset from this file.
	    ++$iFile;
	    $thisOutputGroup= $mergeFile  ->group("Outputs")->group($outputGroup);
	    $thisDataset    = $thisOutputGroup->dataset($dataset);
	    $datasetValues  = $thisDataset->get();
	    # Adjust dataset if necessary.
	    $datasetValues += $offset if ( $dataset eq "mergerTreeStartIndex" );
	    # Append to the combined dataset.
	    $newDatasetValues = $newDatasetValues->append($datasetValues);
            # Adjust the offset.
	    $offset += $offsets[$iFile];
	}
	# Create a new dataset in the output file.
	$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
						  parent  => $newOutputGroup,
						  fileObj => $outputFile
						  );
	# Write the new dataset back to the output file.
	$newDataset->set($newDatasetValues);
    }
}

# Create references for each tree.
# Get the Outputs group.
$outputsGroup = $outputFile->group("Outputs");
# Get a list of subgroups to search.
@outputGroups = $outputsGroup->groups();
# Loop over all outputs groups.
foreach $outputGroup ( @outputGroups ) {
    # Get the nodeData group.
    $nodeDataGroup = $outputFile->group("Outputs")->group($outputGroup)->group("nodeData");
    # Get a list of all datasets in the nodeData group.
    @datasets = $nodeDataGroup->datasets();
    # Read the merger tree indices, starts, offsets and weights.
    $mergerTreeIndex      = $outputFile->group("Outputs")->group($outputGroup)->dataset("mergerTreeIndex"     )->get();
    $mergerTreeStartIndex = $outputFile->group("Outputs")->group($outputGroup)->dataset("mergerTreeStartIndex")->get();
    $mergerTreeCount      = $outputFile->group("Outputs")->group($outputGroup)->dataset("mergerTreeCount"     )->get();
    $mergerTreeWeight     = $outputFile->group("Outputs")->group($outputGroup)->dataset("mergerTreeWeight"    )->get();
    # Loop over each merger tree.
    for($iTree=0;$iTree<nelem($mergerTreeIndex);++$iTree) {
	# Create a group for this tree.
	$treeGroup = $outputFile->group("Outputs")->group($outputGroup)->group("mergerTree".$mergerTreeIndex->index($iTree));
	# Create start and count arrays.
	$start[0] = $mergerTreeStartIndex->index($iTree);
	$count[0] = $mergerTreeCount     ->index($iTree);
	# Loop over all datasets.
	foreach $dataset ( @datasets ) {
	    # Create reference.
	    $nodeData = $nodeDataGroup->dataset($dataset);
	    $treeGroup->reference($nodeData,$dataset,@start,@count);
	}
	# Set the volumeWeight attribute.
	$treeGroup->attrSet("volumeWeight" => $mergerTreeWeight->index($iTree));
    }
}

# Hash of groups that contain trees.
%groupsWithTrees = (
		    "mergerTreeStructures" => "."
		    );

# Loop over all files to merge, scan for groups containing trees and merge into the output file.
foreach $mergeFile ( @mergeFiles ) {
    # Get list of all groups in root group.
    @rootGroups = $mergeFile->groups;
    # Loop over all groups in root group.
    foreach $rootGroup ( @rootGroups ) {
	# Test if this group may contain trees.
	if ( exists($groupsWithTrees{$rootGroup}) ) {
	    print "-> Combining trees in ".$rootGroup." group\n";
	    # Get the group in the merging and output files.
	    $outputsGroup    = $mergeFile->group($rootGroup);
	    $newOutputsGroup = $outputFile->group($rootGroup);
	    # Get a list of subgroups to search.
	    if ( $groupsWithTrees{$rootGroup} eq "." ) {
		undef(@outputGroups);
		@outputGroups[0] = $outputsGroup;
	    } else {
		@outputGroups    = $outputsGroup->groups;
	    }
	    # Loop over subgroups.
	    foreach $outputGroup ( @outputGroups ) {
		# Get groups into which trees will be stored.
		if ( $groupsWithTrees{$rootGroup} eq "." ) {
		    $thisOutputGroup = $outputsGroup;
		    $newOutputGroup  = $newOutputsGroup;
		} else {
		    $thisOutputGroup = $outputsGroup->group($outputGroup);
		    $newOutputGroup  = $newOutputsGroup->group($outputGroup);
		}
		# Get a list of trees to merge.
		@treeGroups = $thisOutputGroup->groups;
		# Loop over trees.
		foreach $treeGroup ( @treeGroups ) {
		    # Get tree from merge file.
		    $thisTree = $thisOutputGroup->group($treeGroup);
		    # Create group for tree in output file.
		    $newTree = new PDL::IO::HDF5::Group( name    => $treeGroup,
							 parent  => $newOutputGroup,
							 fileObj => $outputFile
							 );
		    # Get a list of datasets in this tree.
		    @datasets = $thisTree->datasets;
		    foreach $dataset ( @datasets ) {
			# Get the dataset.
			$thisDataset = $thisTree->dataset($dataset);
			$datasetValues = $thisDataset->get;
			# Create a new dataset in the output file.
			$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
								  parent  => $newTree,
								  fileObj => $outputFile
								  );
			$newDataset->set($datasetValues);
		    }
		}
	    }
	}
    }
}

# Cumulate global histories from all files.
@histories = ( "historyGasDensity", "historyNodeDensity", "historyStellarDensity", "historyStarFormationRate" );
# Loop over histories.
print "-> Combining histories\n";
$newHistoryGroup = $outputFile->group("globalHistory");
foreach $history ( @histories ) {
    # Add datasets from all files to merge.
    undef($newDataValues);
    foreach $mergeFile ( @mergeFiles ) {
	$dataValues     = $mergeFile->dataset("globalHistory/".$history)->get;
	$newDataValues += $dataValues;
    }
    # Store in output file.
    $newDataset = new PDL::IO::HDF5::Dataset( name    => $history,
					      parent  => $newHistoryGroup,
					      fileObj => $outputFile
					      );
    $newDataset->set($newDataValues);
}
# Non-cumulative history datasets.
@histories = ( "historyExpansion", "historyTime" );
$newHistoryGroup = $outputFile->group("globalHistory");
foreach $history ( @histories ) {
    # Copy dataset from first file to merge.
    $newDataValues = $mergeFiles[0]->dataset("globalHistory/".$history)->get;
    # Store in output file.
    $newDataset = new PDL::IO::HDF5::Dataset( name    => $history,
					      parent  => $newHistoryGroup,
					      fileObj => $outputFile
					      );
    $newDataset->set($newDataValues);
}

# Write completion messsage.
print "-> Done\n";

exit;
