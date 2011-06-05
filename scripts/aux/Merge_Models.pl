#!/usr/bin/env perl
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::IO::HDF5::Group;
use File::Copy;

# Merges Galacticus models.
# Andrew Benson (8-July-2010)

if ($#ARGV < 1) {die("Usage: Merge_Models.pl <model1> ...... <outputModel>")};
$baseFileName = $ARGV[0];
@mergeFileNames = @ARGV[1..$#ARGV-1] if ($#ARGV > 1);
$outputFileName = $ARGV[$#ARGV];

# Copy the base file to the output file.
copy($baseFileName,$outputFileName);

# Open all files.
$baseFile   = new PDL::IO::HDF5($baseFileName);
$outputFile = new PDL::IO::HDF5(">".$outputFileName);
foreach $mergeFileName ( @mergeFileNames ) {
    $mergeFiles[++$#mergeFiles] = new PDL::IO::HDF5($mergeFileName);
}

# Modify the parameters to reflect that this is now a complete model.
$dataset = $outputFile->dataset("Parameters/treeEvolveWorkerNumber");
$value = $dataset->get;
$value->index(0) = 1;
$dataset->set($value);
$dataset = $outputFile->dataset("Parameters/treeEvolveWorkerCount");
$value = $dataset->get;
$value->index(0) = 1;
$dataset->set($value);

# Hash of groups that contain trees.
%groupsWithTrees = (
		    "Outputs"              => "Output",
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
foreach $history ( @histories ) {
    # Get dataset from base file.
    $newDataValues = $baseFile->dataset("globalHistory/".$history)->get;
    # Add datasets from all files to merge.
    foreach $mergeFile ( @mergeFiles ) {
	$dataValues = $mergeFile->dataset("globalHistory/".$history)->get;
	$newDataValues += $dataValues;
    }
    # Store in output file.
    $dataset = $outputFile->dataset("globalHistory/".$history);
    $dataset->set($newDataValues);
}

exit;
