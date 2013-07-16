#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::IO::HDF5::Group;
use Switch;
use Data::UUID;
use Data::Dumper;

# Merges Galacticus models.
# Andrew Benson (8-July-2010)

die("Usage: Merge_Models.pl <model1> ...... <outputModel>")
    unless ( scalar(@ARGV) > 1 );
my @mergeFileNames = @ARGV[0..scalar(@ARGV)-2];
my $outputFileName = $ARGV[-1];

# Quit if output file already exists.
die ("Merge_Models.pl: output file already exists") if ( -e $outputFileName );

# Open all files.
print "-> Opening output file: ".$outputFileName."\n";
my $outputFile = new PDL::IO::HDF5(">".$outputFileName);
my @mergeFiles;
foreach my $mergeFileName ( @mergeFileNames ) {
    print "-> Opening file for merge: ".$mergeFileName."\n";
    push(@mergeFiles,new PDL::IO::HDF5($mergeFileName));
}

# Copy UUID attributes.
my @mergeUUIDs;
foreach my $mergeFile ( @mergeFiles ) {
	push(@mergeUUIDs,$mergeFile->attrGet("UUID"));
}
$outputFile->attrSet("UUIDs" => join(":",@mergeUUIDs));
my $ug    = new Data::UUID;
my $uuid1 = $ug->create_str();
$outputFile->attrSet("UUID" => $uuid1);

# Copy parameters and version groups from merge files to output file.
print "-> Copying Parameters and Version groups\n";
foreach my $groupName ( "Parameters", "Version" ) {
    # Open the groups.
    my $group    = $mergeFiles[0]->group($groupName);
    my $newGroup = $outputFile   ->group($groupName);
    &Copy_Attributes($group,$newGroup);
}

# Copy output group attributes.
my $outputsGroup = $mergeFiles[0]->group("Outputs");
my @outputGroups = $outputsGroup->groups;
# Loop over all outputs groups.
foreach my $outputGroupName ( @outputGroups ) {
    my $outputGroup    = $outputsGroup->group($outputGroupName);
    my $newOutputGroup = $outputFile->group("Outputs")->group($outputGroupName);
    &Copy_Attributes($outputGroup,$newOutputGroup);
}

# Modify the parameters to reflect that this is now a complete model.
print "-> Modifying parameter attributes of base file\n";
my $group = $outputFile->group("Parameters");
$group->attrSet(
		"treeEvolveWorkerNumber" => pdl 1,
		"treeEvolveWorkerCount" => pdl 1
		);

# Combine the datasets in the nodeData group.
print "-> Combining main datasets\n";
# Get the Outputs group from the output file.
$outputsGroup = $mergeFiles[0]->group("Outputs");
# Get a list of subgroups to search.
@outputGroups = $outputsGroup->groups;
# Declare a structure to hold merger tree offset data.
my $offsets;
# Loop over all outputs groups.
foreach my $outputGroup ( @outputGroups ) {
    # Get the node data group.
    my $nodeDataGroup    = $outputsGroup->group($outputGroup."/nodeData");
    my $newNodeDataGroup = $outputFile->group("Outputs")->group($outputGroup)->group("nodeData");
    # Get a list of all datasets.
    my @datasets = $nodeDataGroup->datasets;
    # Loop over datasets.
    foreach my $dataset ( @datasets ) {
	my $newDatasetValues = pdl [];
	# Loop over merge files.
	my $iFile = -1;
	foreach my $mergeFile ( @mergeFiles ) {
	    # Read the dataset from this file.
	    ++$iFile;
	    my $outputsGroup  = $mergeFile    ->group("Outputs");
	    my $nodeDataGroup = $outputsGroup ->group($outputGroup."/nodeData");
	    my $thisDataset   = $nodeDataGroup->dataset($dataset);
	    # Test if the dataset exists.
	    if ( my $testID = $thisDataset->IDget ) {
		my $datasetValues = $thisDataset  ->get();
		$offsets->{$outputGroup}->{$iFile} = nelem($datasetValues);
		# Append to the combined dataset.
		$newDatasetValues = $newDatasetValues->append($datasetValues);
	    }
	}
	# Create a new dataset in the output file.
	my $newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
						     parent  => $newNodeDataGroup,
						     fileObj => $outputFile
						     );
	# Write the new dataset back to the output file.
	$newDataset->set($newDatasetValues);

	# Copy any attributes to the new dataset.
	my $thisDataset = $mergeFiles[0]->group("Outputs")->group($outputGroup)->group("nodeData")->dataset($dataset);
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
foreach my $outputGroup ( @outputGroups ) {
    # Get the new output group.
    my $newOutputGroup = $outputFile->group("Outputs")->group($outputGroup);
    # Loop over the tree datasets.
    foreach my $dataset ( "mergerTreeIndex", "mergerTreeCount", "mergerTreeStartIndex", "mergerTreeWeight" ) {
	# Reset values.
	my $newDatasetValues = pdl [];
	# Loop over merge files.
	my $iFile = -1;
	my $offset = 0;
	foreach my $mergeFile ( @mergeFiles ) {
	    # Read the dataset from this file.
	    ++$iFile;
	    my $thisOutputGroup= $mergeFile  ->group("Outputs")->group($outputGroup);
	    my $thisDataset    = $thisOutputGroup->dataset($dataset);
	    my $datasetValues  = $thisDataset->get();
	    # Adjust dataset if necessary.
	    $datasetValues += $offset if ( $dataset eq "mergerTreeStartIndex" );
	    # Append to the combined dataset.
	    $newDatasetValues = $newDatasetValues->append($datasetValues);
            # Adjust the offset.
	    $offset += $offsets->{$outputGroup}->{$iFile} if ( exists($offsets->{$outputGroup}->{$iFile}) );
	}
	# Create a new dataset in the output file.
	my $newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
						     parent  => $newOutputGroup,
						     fileObj => $outputFile
						     );
	# Write the new dataset back to the output file.
	$newDataset->set($newDatasetValues);
    }
}

# Create references for each tree.
my @parameters = $mergeFiles[0]->group('Parameters')->attrGet("mergerTreeOutputReferences");
if ( $parameters[0] eq "true" ) {
    print "-> Creating merger tree references\n";
    # Get the Outputs group.
    $outputsGroup = $outputFile->group("Outputs");
    # Get a list of subgroups to search.
    @outputGroups = $outputsGroup->groups();
    # Loop over all outputs groups.
    foreach my $outputGroup ( @outputGroups ) {
	# Get the nodeData group.
	my $nodeDataGroup = $outputFile->group("Outputs")->group($outputGroup)->group("nodeData");
	# Get a list of all datasets in the nodeData group.
	my @datasets = $nodeDataGroup->datasets();
	# Read the merger tree indices, starts, offsets and weights.
	my $mergerTreeIndex      = $outputFile->group("Outputs")->group($outputGroup)->dataset("mergerTreeIndex"     )->get();
	my $mergerTreeStartIndex = $outputFile->group("Outputs")->group($outputGroup)->dataset("mergerTreeStartIndex")->get();
	my $mergerTreeCount      = $outputFile->group("Outputs")->group($outputGroup)->dataset("mergerTreeCount"     )->get();
	my $mergerTreeWeight     = $outputFile->group("Outputs")->group($outputGroup)->dataset("mergerTreeWeight"    )->get();
	# Loop over each merger tree.
	for(my $iTree=0;$iTree<nelem($mergerTreeIndex);++$iTree) {
	    # Create a group for this tree.
	    my $treeGroup = $outputFile->group("Outputs")->group($outputGroup)->group("mergerTree".$mergerTreeIndex->index($iTree));
	    # Create start and count arrays.
	    my @start;
	    my @count;
	    $start[0] = $mergerTreeStartIndex->index($iTree);
	    $count[0] = $mergerTreeCount     ->index($iTree);
	    if ( $count[0] > 0 ) {
		# Loop over all datasets.
		foreach my $dataset ( @datasets ) {
		    # Create reference.
		    my $nodeData = $nodeDataGroup->dataset($dataset);
		    $treeGroup->reference($nodeData,$dataset,@start,@count);
		}
		# Set the volumeWeight attribute.
		$treeGroup->attrSet("volumeWeight" => $mergerTreeWeight->index($iTree));
	    }
	}
    }
}

# Other datasets.
my %outputRules = (
		   globalHistory          => {
		       "^historyExpansion\$"                                   => "singleCopy",
		       "^historyTime\$"                                        => "singleCopy",
		       "^historyGasDensity\$"                                  => "cumulate",
		       "^historyHotGasDensity\$"                               => "cumulate",
		       "^historyNodeDensity\$"                                 => "cumulate",
		       "^historyStellarDensity\$"                              => "cumulate",
		       "^historyDiskStellarDensity\$"                          => "cumulate",
		       "^historySpheroidStellarDensity\$"                      => "cumulate",
		       "^historyStarFormationRate\$"                           => "cumulate",
		       "^historyDiskStarFormationRate\$"                       => "cumulate",
		       "^historySpheroidStarFormationRate\$"                   => "cumulate"
		       },
		   haloModel              => {
		       "^powerSpectrum\$"                                      => "singleCopy",
		       "^wavenumber\$"                                         => "singleCopy",
		       "^Output\\d+\\/mergerTree\\d+\\/fourierProfile\\d+\$"   => "copy"
		       },
		   massAccretionHistories => {
		       "^mergerTree\\d+\\/node*"                               => "copy"
		       },
		   mergerTreeStructures   => {
		       "^mergerTree\\d+\\/*"                                   => "copy"
		       },
		   starFormationHistories => {
		       "^metallicities\$"                                      => "singleCopy",
		       "^Output\\d+\\/mergerTree\\d+\\/(disk|spheroid).*\$"    => "copy"
                       },
		   blackHole              => {
		       "^Output\\d+\\/.*\$"                                    => "append"
		       },
		   blackHoleMergers       => {
		       "^massBlackHole1\$"                                     => "append",
		       "^massBlackHole2\$"                                     => "append",
		       "^timeOfMerger\$"                                       => "append",
		       "^volumeWeight\$"                                       => "append"
                       },
		   grasilSEDs             => {
		       "^Output\\d+\\/mergerTree\\d+\\/node\\d+\\/SED"         => "copy",
		       "^Output\\d+\\/mergerTree\\d+\\/node\\d+\\/inclination" => "copy",
		       "^Output\\d+\\/mergerTree\\d+\\/node\\d+\\/wavelength"  => "copy"
		       }
		   );
my %availableGroupNames;
foreach my $mergeFile ( @mergeFiles ) {
    my @availableGroups = $mergeFile->groups();
    foreach my $availableGroup ( @availableGroups ) {
	$availableGroupNames{$availableGroup} = 1;
    }
}
foreach my $outputGroup ( keys(%outputRules) ) {
    # Write message.
    print "-> Checking for group ".$outputGroup."\n";
    if ( exists( $availableGroupNames{$outputGroup}) ) {
	my $foundGroup = 0;
	# Loop over all merge files.
	my $isFirstFile = 1;
	foreach my $mergeFile ( @mergeFiles ) {
	    # Check if the group exists in this file.
	    my @availableGroups = $mergeFile->groups();
	    my $groupIsPresent = 0;
	    foreach my $availableGroup ( @availableGroups ) {
		$groupIsPresent = 1 if ( $availableGroup eq $outputGroup );
	    }
	    if ( $groupIsPresent == 1 ) {
		# Place the base group on a stack.
		my @groupStack = ( $outputGroup );
		# Process groups until none remain.
		while ( scalar(@groupStack) > 1 ) {
		    # Pop a group from the stack.
		    my $thisGroup = shift @groupStack;
		    # Get a list of groups.
		    my @groups = $mergeFile->group($thisGroup)->groups;
		    # Ensure all groups exist.
		    foreach my $group ( @groups ) {
			my $parentGroup = $outputFile->group($thisGroup);
			unless ( my $testID = $parentGroup->group($group)->IDget ) {
			    my $newGroup = new PDL::IO::HDF5::Group( name    => $group,
								     parent  => $parentGroup,
								     fileObj => $outputFile
								     );
			    my $originalGroup = $mergeFile->group($thisGroup);
			    &Copy_Attributes($originalGroup,$newGroup);		    
			}
		    }
		    # Push groups to the stack.
		    push(@groupStack,map {$thisGroup."/".$_} @groups);
		    # Get list of datasets.
		    my @datasets = $mergeFile->group($thisGroup)->datasets;
		    # Loop over all datasets.
		    foreach my $dataset ( @datasets ) {
			# Write a message.
			print "  -> processing\n" if ( $foundGroup == 0 );
			$foundGroup = 1;
			# Check for a match to an instruction.
			my %instructions = %{$outputRules{$outputGroup}};
			my $action = "";
			(my $fullPath = $thisGroup) =~ s/^$outputGroup//;
			($fullPath .= "/".$dataset) =~ s/^\///;
			foreach my $instruction ( keys(%instructions) ) {
			    $action = $instructions{$instruction} if ( $fullPath =~ m/$instruction/ );
			}
			switch ($action) {
			    case ( "singleCopy" ) {
				if ( $isFirstFile == 1 ) {
				    my $thisGroupObject = $outputFile->group($thisGroup);
				    my $newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
										 parent  => $thisGroupObject,
										 fileObj => $outputFile
										 );
				    my $originalDataset = $mergeFile->group($thisGroup)->dataset($dataset);
				    my $datasetValues = $originalDataset->get;
				    $newDataset->set($datasetValues);
				    &Copy_Attributes($originalDataset,$newDataset);
				}
			    }
			    case ( "copy" ) {
				my $thisGroupObject = $outputFile->group($thisGroup);
				my $newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									     parent  => $thisGroupObject,
									     fileObj => $outputFile
									     );
				my $originalDataset = $mergeFile->group($thisGroup)->dataset($dataset);
				my $datasetValues = $originalDataset->get;
				$newDataset->set($datasetValues);
				&Copy_Attributes($originalDataset,$newDataset);
			    }
			    case ( "cumulate" ) {
				my $thisGroupObject = $outputFile->group($thisGroup);
				my $originalDataset = $mergeFile->group($thisGroup)->dataset($dataset);
				my $datasetValues = $originalDataset->get;
				my $newDataset;
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
			    case ( "append" ) {
				my $thisGroupObject = $outputFile->group($thisGroup);
				my $originalDataset = $mergeFile ->group($thisGroup)->dataset($dataset);
				my $datasetValues   = $originalDataset->get();
				my $newDataset;
				if ( $isFirstFile == 1 ) {
				    $newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									      parent  => $thisGroupObject,
									      fileObj => $outputFile
									      );
				} else {
				    $newDataset    = $thisGroupObject->dataset($dataset);
				    $datasetValues = $datasetValues->append($newDataset->get());
				}
				$newDataset->set($datasetValues,unlimited => 1);
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
}
# Write completion messsage.
print "-> Done\n";

exit;

sub Copy_Attributes {
    # Copy all attributes from an object in the input files to the output file.
    # Get the objects from and to.
    my $objectFrom = shift;
    my $objectTo   = shift;
  
    # Copy all attributes.
    my @attributes = $objectFrom->attrs();
    foreach my $attribute ( @attributes ) {
	my @attrValue = $objectFrom->attrGet($attribute);
	$objectTo->attrSet($attribute => $attrValue[0]);
    }
}
