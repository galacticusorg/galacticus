#!/usr/bin/env perl
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::IO::HDF5::Group;
use Switch;
use Data::Dumper;

# Merges Galacticus models.
# Andrew Benson (8-July-2010)

if ($#ARGV < 1) {die("Usage: Merge_Models.pl <model1> ...... <outputModel>")};
@mergeFileNames = @ARGV[0..$#ARGV-1];
$outputFileName = $ARGV[$#ARGV];

# Quit if output file already exists.
die ("Merge_Models.pl: output file already exists") if ( -e $outputFileName );

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
    &Copy_Attributes($group,$newGroup);
}

# Copy output group attributes.
$outputsGroup = $mergeFiles[0]->group("Outputs");
@outputGroups = $outputsGroup->groups;
# Loop over all outputs groups.
foreach $outputGroupName ( @outputGroups ) {
    $outputGroup    = $outputsGroup->group($outputGroupName);
    $newOutputGroup = $outputFile->group("Outputs")->group($outputGroupName);
    &Copy_Attributes($outputGroup,$newOutputGroup);
}

# Modify the parameters to reflect that this is now a complete model.
print "-> Modifying parameter attributes of base file\n";
$group = $outputFile->group("Parameters");
$group->attrSet(
		"treeEvolveWorkerNumber" => pdl 1,
		"treeEvolveWorkerCount" => pdl 1
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

	# Copy any attributes to the new dataset.
	$thisDataset = $mergeFiles[0]->group("Outputs")->group($outputGroup)->group("nodeData")->dataset($dataset);
	&Copy_Attributes($thisDataset,$newDataset);


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
print "-> Creating merger tree references\n";
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
	if ( $count[0] > 0 ) {
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
}

# Other datasets.
%outputRules = (
    globalHistory => {
	"^historyExpansion\$"                                 => "singleCopy",
	"^historyTime\$"                                      => "singleCopy",
	"^historyGasDensity\$"                                => "cumulate",
	"^historyNodeDensity\$"                               => "cumulate",
	"^historyStellarDensity\$"                            => "cumulate",
	"^historyStarFormationRate\$"                         => "cumulate"
    },
    haloModel => {
	"^powerSpectrum\$"                                    => "singleCopy",
	"^wavenumber\$"                                       => "singleCopy",
	"^Output\\d+\\/mergerTree\\d+\\/fourierProfile\\d+\$" => "copy"
    },
    massAccretionHistories => {
	"^mergerTree\\d+\\/node*"                             => "copy"
    },
    mergerTreeStructures => {
	"^mergerTree\\d+\\/*"                                 => "copy"
    }
    );
@availableGroups = $mergeFiles[0]->groups();
foreach $availableGroup ( @availableGroups ) {
    $availableGroupNames{$availableGroup} = 1;
}
foreach $outputGroup ( keys(%outputRules) ) {
    # Write message.
    print "-> Checking for group ".$outputGroup."\n";
    if ( exists( $availableGroupNames{$outputGroup}) ) {
	$foundGroup = 0;
	# Loop over all merge files.
	$isFirstFile = 1;
	foreach $mergeFile ( @mergeFiles ) {
	    # Place the base group on a stack.
	    @groupStack = ( $outputGroup );
	    # Process groups until none remain.
	    while ( $#groupStack >= 0 ) {
		# Pop a group from the stack.
		$thisGroup = shift @groupStack;
		# Get a list of groups.
		@groups = $mergeFile->group($thisGroup)->groups;
		# Ensure all groups exist.
		foreach $group ( @groups ) {
		    $parentGroup = $outputFile->group($thisGroup);
		    unless ( $testID = $parentGroup->group($group)->IDget ) {
			$newGroup = new PDL::IO::HDF5::Group( name    => $group,
							      parent  => $parentGroup,
							      fileObj => $outputFile
			    );
			$originalGroup = $mergeFile->group($thisGroup);
			&Copy_Attributes($originalGroup,$newGroup);		    
		    }
		}
		# Push groups to the stack.
		push(@groupStack,map {$thisGroup."/".$_} @groups);
		# Get list of datasets.
		@datasets = $mergeFile->group($thisGroup)->datasets;
		# Loop over all datasets.
		foreach $dataset ( @datasets ) {
		    # Write a message.
		    print "  -> processing\n" if ( $foundGroup == 0 );
		    $foundGroup = 1;
		    # Check for a match to an instruction.
		    %instructions = %{$outputRules{$outputGroup}};
		    $action = "";
		    ($fullPath = $thisGroup) =~ s/^$outputGroup//;
		    ($fullPath .= "/".$dataset) =~ s/^\///;
		    foreach $instruction ( keys(%instructions) ) {
			$action = $instructions{$instruction} if ( $fullPath =~ m/$instruction/ );
		    }
		    switch ($action) {
			case ( "singleCopy" ) {
			    if ( $isFirstFile == 1 ) {
				$thisGroupObject = $outputFile->group($thisGroup);
				$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									  parent  => $thisGroupObject,
									  fileObj => $outputFile
				    );
				$originalDataset = $mergeFile->group($thisGroup)->dataset($dataset);
				$datasetValues = $originalDataset->get;
				$newDataset->set($datasetValues);
				&Copy_Attributes($originalDataset,$newDataset);
			    }
			}
			case ( "copy" ) {
			    $thisGroupObject = $outputFile->group($thisGroup);
			    $newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
								      parent  => $thisGroupObject,
								      fileObj => $outputFile
				);
			    $originalDataset = $mergeFile->group($thisGroup)->dataset($dataset);
			    $datasetValues = $originalDataset->get;
			    $newDataset->set($datasetValues);
			    &Copy_Attributes($originalDataset,$newDataset);
			}
			case ( "cumulate" ) {
			    $thisGroupObject = $outputFile->group($thisGroup);
			    $originalDataset = $mergeFile->group($thisGroup)->dataset($dataset);
			    $datasetValues = $originalDataset->get;
			    if ( $isFirstFile == 1 ) {
				$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									  parent  => $thisGroupObject,
									  fileObj => $outputFile
				    );
			    } else {
				$newDataset = $thisGroupObject->dataset($dataset);
				$datasetValues += $newDataset->get;
			    }
			    $newDataset->set($datasetValues);
			    &Copy_Attributes($originalDataset,$newDataset) if ( $isFirstFile == 1 );
			}
		    }
		}
	    }
	    
	    # Mark that we have processed the first file.
	    $isFirstFile = 0;
	    
	}
    }
}
# Write completion messsage.
print "-> Done\n";

exit;

sub Copy_Attributes {
    # Copy all attributes from an object in the input files to the output file.
    # Get the objects from and to.
    $objectFrom = shift;
    $objectTo   = shift;
  
    # Copy all attributes.
    @attributes = $objectFrom->attrs();
    foreach $attribute ( @attributes ) {
	@attrValue = $objectFrom->attrGet($attribute);
	$objectTo->attrSet($attribute => $attrValue[0]);
    }
}
