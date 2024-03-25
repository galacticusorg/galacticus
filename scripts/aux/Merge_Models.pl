#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::IO::HDF5::Group;
use Data::UUID;
use Data::Dumper;
use Scalar::Util qw(reftype);
use List::Uniq ':all';

# Merges Galacticus models.
# Andrew Benson (8-July-2010)

# Storage for weight datasets.
my $weights;

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
    die("Merge_Models.pl: file '".$mergeFileName."' does not exist")
	unless ( -e $mergeFileName );
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
foreach my $mergeFile ( @mergeFiles ) {
    foreach my $topLevelGroupName ( "Parameters", "Version" ) {
	my @groupStack = ( $topLevelGroupName );
	while ( scalar(@groupStack) > 0 ) {
	    # Pop a group from the stack.
	    my $groupName = pop(@groupStack);
	    # Open the groups.
	    my $group    = $mergeFile ->group($groupName);
	    my $newGroup = $outputFile->group($groupName);
	    # Find any subgroups.
	    push(@groupStack,map {$groupName."/".$_} $group->groups());
	    &Copy_Attributes($group,$newGroup);
	}
    }
}

# Modify the parameters to reflect that this is now a complete model.
print "-> Modifying parameter attributes of base file\n";
my $group = $outputFile->group("Parameters");
$group->attrSet("treeEvolveWorkerNumber" => pdl 1);
$group->attrSet("treeEvolveWorkerCount"  => pdl 1);

# Find all output times and datasets across all outputs.
print "-> Identifying outputs and datasets\n";
my %outputGroupsAvailable;
my $datasetsAvailable;
my $iFile = -1;
foreach my $mergeFile ( @mergeFiles ) {
    # Read the dataset from this file.
    ++$iFile;
    next
	unless ( grep {$_ eq "Outputs"} $mergeFile->groups() );
    print "  -> file: ".$mergeFileNames[$iFile]."\n";
    my $outputsGroup = $mergeFile->group("Outputs");
    my @outputGroups = $outputsGroup->groups();
    push(@{$outputGroupsAvailable{$_}},$mergeFile)
	foreach ( @outputGroups );
    foreach my $outputGroup ( @outputGroups ) {
	print "    -> output: ".$outputGroup."\n";
	my $nodeDataGroup = $outputsGroup->group($outputGroup)->group("nodeData");
	push(@{$datasetsAvailable->{$outputGroup}->{$_}},$mergeFile)
	    foreach ( $nodeDataGroup->datasets() );
    }
}
# Copy output group content.
foreach my $outputGroupName ( sort(keys(%outputGroupsAvailable)) ) {
    print "--> Copy attributes for ".$outputGroupName."\n";
    my $outputGroup    = $outputGroupsAvailable{$outputGroupName}->[0]->group("Outputs")->group($outputGroupName);
    my $newOutputGroup = $outputFile->group("Outputs")->group($outputGroupName);
    &Copy_Attributes($outputGroup,$newOutputGroup);
    # Get the new nodeData group.
    my $newNodeDataGroup = $outputFile->group("Outputs")->group($outputGroupName)->group("nodeData");
    # Combine the datasets in the nodeData group.
    print "---> Combining main datasets\n";
    foreach my $datasetName ( sort(keys(%{$datasetsAvailable->{$outputGroupName}})) ) {
	print "----> Combining dataset ".$datasetName."\n";
	# Create a new dataset in the output file.
	my $newDataset = new PDL::IO::HDF5::Dataset( name    => $datasetName,
						     parent  => $newNodeDataGroup,
						     fileObj => $outputFile
	    );
	my $newDatasetValues = pdl->null();
	# Iterare over all merge models with this dataset.
	my $iFile = -1;
	foreach my $mergeFile ( @{$datasetsAvailable->{$outputGroupName}->{$datasetName}} ) {
	    ++$iFile;
	    last
		if ( $datasetName =~ m/Columns$/ && $iFile > 0 );
	    # Get the nodeData group from the merge file.
	    my $nodeDataGroup = $mergeFile->group("Outputs")->group($outputGroupName)->group("nodeData");
	    my $dataset       = $nodeDataGroup->dataset($datasetName);
	    my $datasetValues = $dataset      ->get();
	    # Append to the combined dataset.
	    $newDatasetValues = $newDatasetValues->glue(-1,$datasetValues);
	}
	# Write the new dataset back to the output file.
	$newDataset->set($newDatasetValues);
	# Copy any attributes to the new dataset.
	my $thisDataset = $datasetsAvailable->{$outputGroupName}->{$datasetName}->[0]->group("Outputs")->group($outputGroupName)->group("nodeData")->dataset($datasetName);
	&Copy_Attributes($thisDataset,$newDataset);
    }
}
# Combine tree datasets.
print "-> Combining tree datasets\n";
# Iterate over output groups.
foreach my $outputGroupName ( sort(keys(%outputGroupsAvailable)) ) {
    # Get the new output group.
    my $newOutputGroup = $outputFile->group("Outputs")->group($outputGroupName);
    # Loop over the tree datasets.
    my @offsets;
    foreach my $datasetName ( "mergerTreeIndex", "mergerTreeCount", "mergerTreeStartIndex", "mergerTreeWeight" ) {
	# Create new dataset.
	my $newDatasetValues = pdl->null();
	# Iterate over merge files.
	my $offset = 0;
	my $iFile = -1;
	foreach my $mergeFile ( @{$outputGroupsAvailable{$outputGroupName}} ) {
	    ++$iFile;
	    # Read the dataset from this file.
	    my $thisOutputGroup= $mergeFile      ->group("Outputs")->group($outputGroupName);
	    my $dataset        = $thisOutputGroup->dataset($datasetName);
	    my $datasetValues  = $dataset->get();
	    # Adjust dataset if necessary.
	    $datasetValues += $offsets[$iFile]
		if ( $datasetName eq "mergerTreeStartIndex" );
	    # Append to the combined dataset.
	    $newDatasetValues = $newDatasetValues->append($datasetValues);
	    # Adjust the offset.
	    if ( $datasetName eq "mergerTreeCount" ) {
		push(@offsets,$offset);
		$offset += sum($datasetValues);
	    }
	}
	# Create a new dataset in the output file.
	my $newDataset = new PDL::IO::HDF5::Dataset( name    => $datasetName,
						     parent  => $newOutputGroup,
						     fileObj => $outputFile
	    );
	# Write the new dataset back to the output file.
	$newDataset->set($newDatasetValues);
	# Create references for each tree.
	if ( scalar(keys(%{$datasetsAvailable->{$outputGroupName}})) > 0 ) {
	    (my $datasetFirstName) = sort(keys(%{$datasetsAvailable->{$outputGroupName}}));
	    my @parameters = $datasetsAvailable->{$outputGroupName}->{$datasetFirstName}->[0]->group('Parameters')->group('mergerTreeOutputter')->attrGet("outputReferences");
	    if ( $parameters[0] eq "true" ) {
		print "-> Creating merger tree references\n";
		# Get the nodeData group.
		my $nodeDataGroup = $outputFile->group("Outputs")->group($outputGroupName)->group("nodeData");
		# Read the merger tree indices, starts, offsets and weights.
		my $mergerTreeIndex      = $outputFile->group("Outputs")->group($outputGroupName)->dataset("mergerTreeIndex"     )->get();
		my $mergerTreeStartIndex = $outputFile->group("Outputs")->group($outputGroupName)->dataset("mergerTreeStartIndex")->get();
		my $mergerTreeCount      = $outputFile->group("Outputs")->group($outputGroupName)->dataset("mergerTreeCount"     )->get();
		my $mergerTreeWeight     = $outputFile->group("Outputs")->group($outputGroupName)->dataset("mergerTreeWeight"    )->get();
		# Loop over each merger tree.
		for(my $iTree=0;$iTree<nelem($mergerTreeIndex);++$iTree) {
		    # Create a group for this tree.
		    my $treeGroup = $outputFile->group("Outputs")->group($outputGroupName)->group("mergerTree".$mergerTreeIndex->index($iTree));
		    # Create start and count arrays.
		    my @start;
		    my @count;
		    $start[0] = $mergerTreeStartIndex->index($iTree);
		    $count[0] = $mergerTreeCount     ->index($iTree);
		    if ( $count[0] > 0 ) {
			foreach my $propertyDatasetName ( sort(keys(%{$datasetsAvailable->{$outputGroupName}})) ) {
			    # Create reference.
			    my $nodeData = $nodeDataGroup->dataset($propertyDatasetName);
			    $treeGroup->reference($nodeData,$propertyDatasetName,@start,@count);
			}
		    }
		    # Set the volumeWeight attribute.
		    $treeGroup->attrSet("volumeWeight" => $mergerTreeWeight->index($iTree));
		}
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
    },
    analysis => {
	"^haloMassFunctionZ[\\d\\.]+\\/mass\$"                           => "singleCopy",
	"^haloMassFunctionZ[\\d\\.]+\\/massFunction\$"                   => "cumulate",
	"^haloMassFunctionZ[\\d\\.]+\\/massFunctionCovariance\$"         => "cumulate",
	"^isolatedHaloMassFunctionZ[\\d\\.]+\\/mass\$"                   => "singleCopy",
	"^isolatedHaloMassFunctionZ[\\d\\.]+\\/massFunction\$"           => "cumulate",
	"^isolatedHaloMassFunctionZ[\\d\\.]+\\/massFunctionCovariance\$" => "cumulate",	
    },
    conditionalMassFunction => {
	"^normalization\$"                         => "cumulate"                         ,
	"^normalizationError\$"                    => "cumulateQuadrature"               ,
	"^normalizationCovariance\$"               => "cumulate"                         ,
	"^normalizationSubhaloMassFunction\$"      => "cumulate"                         ,
	"^normalizationSubhaloMassFunctionError\$" => "cumulateQuadrature"               ,
	"^massParent\$"                            => "singleCopy"                       ,
	"^massRatio\$"                             => "singleCopy"                       ,
	"^redshiftParent\$"                        => "singleCopy"                       ,
	"^redshiftProgenitor\$"                    => "singleCopy"                       ,
	"^conditionalMassFunction\$"               => \&conditionalMassFunction          ,
	"^conditionalMassFunctionError\$"          => \&conditionalMassFunctionError     ,
	"^conditionalMassFunctionCovariance\$"     => \&conditionalMassFunctionCovariance,
	"^primaryProgenitorMassFunction\$"         => \&conditionalMassFunction          ,
	"^primaryProgenitorMassFunctionError\$"    => \&conditionalMassFunctionError     ,
	"^formationRateFunction\$"                 => \&conditionalMassFunction          ,
	"^formationRateFunctionError\$"            => \&conditionalMassFunctionError     ,
	"^subhaloMassFunction\$"                   => \&conditionalMassFunction          ,
	"^subhaloMassFunctionError\$"              => \&conditionalMassFunctionError
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
	my $iFile       = -1;
	foreach my $mergeFile ( @mergeFiles ) {
	    ++$iFile;
	    my $isFirstFile = $iFile == 0            ? 1 : 0;
	    my $isLastFile  = $iFile == $#mergeFiles ? 1 : 0;
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
		while ( scalar(@groupStack) > 0 ) {
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
			print "  -> processing\n"
			    if ( $foundGroup == 0 );
			$foundGroup = 1;
			# Check for a match to an instruction.
			my %instructions = %{$outputRules{$outputGroup}};
			my $action;
			(my $fullPath = $thisGroup) =~ s/^$outputGroup//;
			($fullPath .= "/".$dataset) =~ s/^\///;
			foreach my $instruction ( keys(%instructions) ) {
			    $action = $instructions{$instruction} if ( $fullPath =~ m/$instruction/ );
			}
			next
			    unless ( defined($action) );
			# Process the dataset.
			if ( reftype($action) && reftype($action) eq "CODE" ) {
			    &{$action}($thisGroup,$dataset,$isFirstFile,$isLastFile,$outputFile,$mergeFile);
			} elsif ( $action eq "singleCopy" ) {
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
			} elsif ( $action eq "copy" ) {
			    my $thisGroupObject = $outputFile->group($thisGroup);
			    my $newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									 parent  => $thisGroupObject,
									 fileObj => $outputFile
				);
			    my $originalDataset = $mergeFile->group($thisGroup)->dataset($dataset);
			    my $datasetValues = $originalDataset->get;
			    $newDataset->set($datasetValues);
			    &Copy_Attributes($originalDataset,$newDataset);
			} elsif ( $action eq "cumulate" ) {
			    # Check for existance of dataset in the file to be merged.
			    if ( grep {$_ eq $dataset} $mergeFile->group($thisGroup)->datasets() ) {
				my $originalDataset = $mergeFile->group($thisGroup)->dataset($dataset);
				my $datasetValues   = $originalDataset->get();
				my $thisGroupObject = $outputFile->group($thisGroup);
				# Check for existance of dataset in output file.
				my $newDataset;
				my $datasetExistsInOutput = grep {$_ eq $dataset } $thisGroupObject->datasets();
				unless ( $datasetExistsInOutput ) {
				    $newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									      parent  => $thisGroupObject,
									      fileObj => $outputFile
					);
				} else {
				    $newDataset     = $thisGroupObject->dataset($dataset);
				    $datasetValues += $newDataset->get;
				}
				$newDataset->set($datasetValues);
				&Copy_Attributes($originalDataset,$newDataset) unless ( $datasetExistsInOutput );
			    }
			} elsif ( $action eq "cumulateQuadrature" ) {
			    my $thisGroupObject = $outputFile->group($thisGroup);
			    my $originalDataset = $mergeFile ->group($thisGroup)->dataset($dataset);
			    my $datasetValues   = $originalDataset->get;
			    $datasetValues .= $datasetValues**2;
			    my $newDataset;
			    if ( $isFirstFile == 1 ) {
				$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									  parent  => $thisGroupObject,
									  fileObj => $outputFile
				    );
			    } else {
				$newDataset     = $thisGroupObject->dataset($dataset);
				$datasetValues += $newDataset->get;
			    }
			    $datasetValues .= sqrt($datasetValues)
				if ( $isLastFile == 1 );
			    $newDataset->set($datasetValues);
			    &Copy_Attributes($originalDataset,$newDataset) if ( $isFirstFile == 1 );
			} elsif ( $action eq "average" ) {
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
				$newDataset     = $thisGroupObject->dataset($dataset);
				$datasetValues += $newDataset->get;
			    }
			    $datasetValues /= scalar(@mergeFiles)
				if ( $isLastFile == 1 );
			    $newDataset->set($datasetValues);
			    &Copy_Attributes($originalDataset,$newDataset) if ( $isFirstFile == 1 );
			} elsif ( $action eq "averageQuadrature" ) {
			    my $thisGroupObject = $outputFile->group($thisGroup);
			    my $originalDataset = $mergeFile ->group($thisGroup)->dataset($dataset);
			    my $datasetValues   = $originalDataset->get;
			    $datasetValues .= $datasetValues**2;
			    my $newDataset;
			    if ( $isFirstFile == 1 ) {
				$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									  parent  => $thisGroupObject,
									  fileObj => $outputFile
				    );
			    } else {
				$newDataset     = $thisGroupObject->dataset($dataset);
				$datasetValues += $newDataset->get;
			    }
			    $datasetValues .= sqrt($datasetValues)/scalar(@mergeFiles)
				if ( $isLastFile == 1 );
			    $newDataset->set($datasetValues);
			    &Copy_Attributes($originalDataset,$newDataset) if ( $isFirstFile == 1 );
			} elsif ( $action eq "averageWeighted" ) {
			    my $thisGroupObject = $outputFile->group($thisGroup);
			    my $originalDataset = $mergeFile->group($thisGroup)->dataset($dataset         );
			    my $originalWeights = $mergeFile->group($thisGroup)->dataset($dataset."Weight");
			    my $datasetValues   = $originalDataset->get;
			    my $datasetWeights  = $originalWeights->get;
			    $datasetValues     *= $datasetWeights;
			    my $newDataset;
			    my $newWeights;
			    if ( $isFirstFile == 1 ) {
				$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									  parent  => $thisGroupObject,
									  fileObj => $outputFile
				    );
				$newWeights = new PDL::IO::HDF5::Dataset( name    => $dataset."Weight",
									  parent  => $thisGroupObject,
									  fileObj => $outputFile
				    );
			    } else {
				$newDataset      = $thisGroupObject->dataset($dataset         );
				$newWeights      = $thisGroupObject->dataset($dataset."Weight");
				my $oldValues    = $newDataset->get;
				my $oldWeights   = $newWeights->get;
				$datasetValues  += $oldValues ;
				$datasetWeights += $oldWeights;
			    }
			    if ( $isLastFile == 1 ) {
				my $nonZeroWeights = which($datasetWeights > 0.0);
				$datasetValues->flat()->($nonZeroWeights) /= $datasetWeights->flat()->($nonZeroWeights);
			    }
			    $newDataset->set($datasetValues );
			    $newWeights->set($datasetWeights);
			    &Copy_Attributes($originalDataset,$newDataset) if ( $isFirstFile == 1 );
			} elsif ( $action eq "averageWeightedQuadrature" ) {
			    my $weightDataset;
			    if ( $dataset =~ m/Error$/ ) {
				($weightDataset = $dataset) =~ s/Error/Weight/;
			    } else {
				$weightDataset = $dataset."Weight";
			    }
			    my $thisGroupObject = $outputFile->group($thisGroup);
			    my $originalDataset = $mergeFile ->group($thisGroup)->dataset($dataset      );
			    my $originalWeights = $mergeFile ->group($thisGroup)->dataset($weightDataset);
			    my $datasetValues   = $originalDataset->get;
			    my $datasetWeights  = $originalWeights->get;
			    $datasetValues     .= ($datasetValues*$datasetWeights)**2;
			    my $newDataset;
			    my $newWeights;
			    if ( $isFirstFile == 1 ) {
				$newDataset = new PDL::IO::HDF5::Dataset( name    => $dataset,
									  parent  => $thisGroupObject,
									  fileObj => $outputFile
				    );
			    } else {
				$newDataset      = $thisGroupObject->dataset($dataset      );
				$newWeights      = $thisGroupObject->dataset($weightDataset);
				my $oldValues    = $newDataset->get;
				my $oldWeights   = $newWeights->get;
				$datasetValues  += $oldValues ;
				$datasetWeights += $oldWeights;
			    }
			    if ( $isLastFile == 1 ) {
				my $nonZeroWeights = which($datasetValues > 0.0);
				$datasetValues->flat()->($nonZeroWeights) .= sqrt($datasetValues->flat()->($nonZeroWeights))/$datasetWeights->flat()->($nonZeroWeights);
			    }
			    $newDataset->set($datasetValues);
			    &Copy_Attributes($originalDataset,$newDataset) if ( $isFirstFile == 1 );
			} elsif ( $action eq "append" ) {
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

sub conditionalMassFunction {
    # Accumulate conditional mass function (or similar), weighting by the relevant normalization.
    my $group       = shift();
    my $dataset     = shift();
    my $isFirstFile = shift();
    my $isLastFile  = shift();
    my $outputFile  = shift();
    my $mergeFile   = shift();
    # Determine the normalization dataset.
    my $normalizationName = $dataset =~ m/subhalo/ ? "normalizationSubhaloMassFunction" : "normalization";
    # Read the dataset and normalization.
    my $values        = $mergeFile->group($group)->dataset($dataset          )->get();
    my $normalization = $mergeFile->group($group)->dataset($normalizationName)->get();
    # Scale by the normalization.
    my @order       = 0..$values->ndims()-1;
    my @reorder     = reverse(@order);
    my @postOrder   = 0..$normalization->ndims()-1;
    my @postReorder = reverse(@postOrder);
    my $iterator    = pdl zeroes($values->ndims()-$normalization->ndims());
    my $shape       = $values->reorder(@reorder)->shape();
    do {
	$values->reorder(@reorder)->range($iterator)->reorder(@postReorder) *= $normalization;
    } until ( &iteratorUpdate($iterator,$shape) );
    # Read the current cumulated result if necessary, and accumulate the normalization.
    if ( $isFirstFile == 1 ) {
	$weights->{$group.":".$dataset} = $normalization;
    } else {
	my $valuesPrevious = $outputFile->group($group)->dataset($dataset)->get();
	$values                         += $valuesPrevious;
	$weights->{$group.":".$dataset} += $normalization;
    }
    # If this is the final file, normalize the result.
    if ( $isLastFile == 1 ) {
	$iterator .= 0;
	do {
	    my $nonZero = which($weights->{$group.":".$dataset} > 0.0);
	    $values->reorder(@reorder)->range($iterator)->reorder(@postReorder)->flat()->($nonZero) /= $weights->{$group.":".$dataset}->flat()->($nonZero);
	} until ( &iteratorUpdate($iterator,$shape) );
    }
    # Store the dataset back to the output file.
    $outputFile->group($group)->dataset($dataset)->set($values);
}

sub conditionalMassFunctionError {
    # Accumulate conditional mass function (or similar) error, weighting by the relevant normalization.
    my $group       = shift();
    my $dataset     = shift();
    my $isFirstFile = shift();
    my $isLastFile  = shift();
    my $outputFile  = shift();
    my $mergeFile   = shift();
    # Determine the normalization dataset.
    my $normalizationName = $dataset =~ m/subhalo/ ? "normalizationSubhaloMassFunction" : "normalization";
    # Determine base dataset.
    (my $baseName = $dataset) =~ s/Error$//;
    # Read the dataset and normalization.
    my $values             = $mergeFile->group($group)->dataset($dataset                  )->get();
    my $base               = $mergeFile->group($group)->dataset($baseName                 )->get();
    my $normalization      = $mergeFile->group($group)->dataset($normalizationName        )->get();
    my $normalizationError = $mergeFile->group($group)->dataset($normalizationName."Error")->get();
    # Scale by the normalization.
    my @order       = 0..$values->ndims()-1;
    my @reorder     = reverse(@order);
    my @postOrder   = 0..$normalization->ndims()-1;
    my @postReorder = reverse(@postOrder);
    my $iterator    = pdl zeroes($values->ndims()-$normalization->ndims());
    my $shape       = $values->reorder(@reorder)->shape();
    do {
	my $baseRange   = $base  ->reorder(@reorder)->range($iterator)->reorder(@postReorder);
	my $valuesRange = $values->reorder(@reorder)->range($iterator)->reorder(@postReorder);
	my $nonZero = which($normalization > 0.0);
	$valuesRange->flat()->($nonZero) .= 
	    +(
	      +(
	        +$valuesRange       ->flat()->($nonZero)
	        *$normalization     ->flat()->($nonZero)
	       )**2
	      -(
		+$normalizationError->flat()->($nonZero)
		*$baseRange         ->flat()->($nonZero)
	       )**2
	     );
 	$baseRange->flat()->($nonZero) *= $normalization->flat()->($nonZero);
    } until ( &iteratorUpdate($iterator,$shape) );
    # Read the current cumulated result if necessary, and accumulate the normalization.
    if ( $isFirstFile == 1 ) {
	$weights->{$group.":".$dataset        } =  $normalization;
	$weights->{$group.":".$dataset."Error"} =  $normalizationError**2;
	$weights->{$group.":".$dataset."Base" } =  $base;
    } else {
	my $valuesPrevious = $outputFile->group($group)->dataset($dataset)->get();
	$values                                 +=  $valuesPrevious;
	$weights->{$group.":".$dataset        } +=  $normalization;
	$weights->{$group.":".$dataset."Error"} +=  $normalizationError**2;
	$weights->{$group.":".$dataset."Base" } +=  $base;
    }
    # If this is the final file, normalize the result.
    if ( $isLastFile == 1 ) {
    	$iterator .= 0;
    	do {
    	    my $nonZero = which($weights->{$group.":".$dataset} > 0.0);
    	    my $baseRange   = $weights->{$group.":".$dataset."Base" }->reorder(@reorder)->range($iterator)->reorder(@postReorder);
    	    my $valuesRange = $values                                ->reorder(@reorder)->range($iterator)->reorder(@postReorder);
    	    $valuesRange->flat()->($nonZero) .= 
    		+sqrt(
    		      +$valuesRange                               ->flat()->($nonZero)
    		      +$weights    ->{$group.":".$dataset."Error"}->flat()->($nonZero)
    		      /$weights    ->{$group.":".$dataset        }->flat()->($nonZero)**2
     		      *$baseRange                                 ->flat()->($nonZero)**2
   		     )
		/      $weights    ->{$group.":".$dataset        }->flat()->($nonZero);
	} until ( &iteratorUpdate($iterator,$shape) );
    }
    # Store the dataset back to the output file.
    $outputFile->group($group)->dataset($dataset)->set($values);
}

sub conditionalMassFunctionCovariance {
    # Accumulate conditional mass function (or similar) covariance, weighting by the relevant normalization.
    my $group       = shift();
    my $dataset     = shift();
    my $isFirstFile = shift();
    my $isLastFile  = shift();
    my $outputFile  = shift();
    my $mergeFile   = shift();
    # Determine the normalization dataset.
    my $normalizationName = $dataset =~ m/subhalo/ ? "normalizationSubhaloMassFunction" : "normalization";
    # Determine base dataset.
    (my $baseName = $dataset) =~ s/Covariance$//;
    # Read the dataset and normalization.
    my $values                  = $mergeFile->group($group)->dataset($dataset                       )->get();
    my $base                    = $mergeFile->group($group)->dataset($baseName                      )->get();
    my $normalization           = $mergeFile->group($group)->dataset($normalizationName             )->get();
    my $normalizationCovariance = $mergeFile->group($group)->dataset($normalizationName."Covariance")->get();
    # Scale by the normalization.
    my @order       = 0..$base->ndims()-1;
    my @reorder     = reverse(@order);
    my @postOrder   = 0..$normalization->ndims()-1;
    my @postReorder = reverse(@postOrder);
    my $shape       = $base->reorder(@reorder)->shape();
    my @orderCov    = 0..$values->ndims()-1;
    my @reorderCov  = reverse(@orderCov);
    my $shape1      = $base->shape();
    my $iterator1   = pdl zeroes($base->ndims());
    do {
    	# Extract base and normalization values corresponding to iterator #1.
     	my $base1                        = $base                   ->indexND($iterator1                               );
    	my $normalization1               = $normalization          ->indexND($iterator1->(0:$normalization->ndims()-1));
    	if ( $normalization1 > 0.0 ) {
    	    my $normalizationCovariance1 = $normalizationCovariance->range  ($iterator1->(0:$normalization->ndims()-1));
    	    my $values1                  = $values                 ->range  ($iterator1                               );
    	    my $iterator                 = pdl zeroes($base->ndims()-$normalization->ndims());
    	    do {
    		my $baseRange   = $base   ->reorder(@reorder)->range($iterator)->reorder(@postReorder);
    		my $valuesRange = $values1->reorder(@reorder)->range($iterator)->reorder(@postReorder);
    		my $index       = $iterator1->((0));
    		$valuesRange->(($index),:) .= 
    		    +$valuesRange             ->(($index),:)
    		    *$normalization1
    		    *$normalization           ->(($index),:)   
    		    -$normalizationCovariance1->(($index),:)
    		    *$base1
    		    *$baseRange               ->(($index),:);
    	    } until ( &iteratorUpdate($iterator,$shape) );
    	}
    } until ( &iteratorUpdate($iterator1,$shape1) );
    # Multiply the base values by the normalization.
    my $iterator = pdl zeroes($base->ndims()-$normalization->ndims());
    do {
	my $baseRange = $base->reorder(@reorder)->range($iterator)->reorder(@postReorder);
	my $nonZero   = which($normalization > 0.0);
 	$baseRange->flat()->($nonZero) *= $normalization->flat()->($nonZero);
    } until ( &iteratorUpdate($iterator,$shape) );
    # Read the current cumulated result if necessary, and accumulate the normalization.
    if ( $isFirstFile == 1 ) {
	$weights->{$group.":".$dataset             } =  $normalization;
	$weights->{$group.":".$dataset."Covariance"} =  $normalizationCovariance;
	$weights->{$group.":".$dataset."Base"      } =  $base;
    } else {
	my $valuesPrevious = $outputFile->group($group)->dataset($dataset)->get();
	$values                                      +=  $valuesPrevious;
	$weights->{$group.":".$dataset             } +=  $normalization;
	$weights->{$group.":".$dataset."Covariance"} +=  $normalizationCovariance;
	$weights->{$group.":".$dataset."Base"      } +=  $base;
    }
    # If this is the final file, normalize the result.
    if ( $isLastFile == 1 ) {
	$iterator1 .= 0;
	do {
	    # Extract base and normalization values corresponding to iterator #1.
	    my $base1                        = $weights->{$group.":".$dataset."Base"      }->indexND($iterator1                               );
	    my $normalization1               = $weights->{$group.":".$dataset             }->indexND($iterator1->(0:$normalization->ndims()-1));
	    if ( $normalization1 > 0.0 ) {
		my $normalizationCovariance1 = $weights->{$group.":".$dataset."Covariance"}->range  ($iterator1->(0:$normalization->ndims()-1));
		my $values1                  = $values                                     ->range  ($iterator1                               );
		my $iterator                 = pdl zeroes($base->ndims()-$normalization->ndims());
		do {
		    my $baseRange   = $weights->{$group.":".$dataset."Base"}->reorder(@reorder)->range($iterator)->reorder(@postReorder);
		    my $valuesRange = $values1                              ->reorder(@reorder)->range($iterator)->reorder(@postReorder);
		    my $index       = $iterator1->((0));
		    my $nonZero     = which($weights->{$group.":".$dataset}->(($index),:) > 0.0);
		    $valuesRange->(($index),$nonZero) .=
			+$valuesRange                                    ->(($index),$nonZero)
			/$normalization1
			/$weights                 ->{$group.":".$dataset}->(($index),$nonZero)
			+$normalizationCovariance1                       ->(($index),$nonZero)
			*$base1
			*$baseRange                                      ->(($index),$nonZero)
			/$normalization1                                                      **2
			/$weights                 ->{$group.":".$dataset}->(($index),$nonZero)**2;
		} until ( &iteratorUpdate($iterator,$shape) );
	    }
	} until ( &iteratorUpdate($iterator1,$shape1) );
    }
    # Store the dataset back to the output file.
    $outputFile->group($group)->dataset($dataset)->set($values);
}

sub iteratorUpdate {
    my $iterator = shift();
    my $shape    = shift();
    my $dim      = $iterator->dim(0);
    while ( $dim > -1 ) {
	--$dim;
	$iterator->(($dim)) += 1;
	if ( $dim > 0 && $iterator->(($dim)) >= $shape->(($dim)) ) {
	    $iterator->(($dim)) .= 0;
	} else {
	    last;
	}
    }
    my $result = $iterator->((0)) >= $shape->((0)) ? 1 : 0;
    return $result;
}
