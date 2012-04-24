# Contains a Perl module which implements various useful functionality for extracting data from Galacticus HDF5 files.

package HDF5;
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::NiceSlice;
use Data::Dumper;

our %galacticusFunctions = ();

# The following is a hash which contains knowledge on how to resolve issues relating to missing datasets. The hash keys are
# regular expressions which match to data set names. The hash keys are functions which should return a hash describing the
# parameters, and values to be set or appended, to those parameters.
my %propertyKnowledgeBase = (
    "^nodeVirial(Radius|Velocity)\$" => sub {my ($d,$r)=@_;return {add => {outputVirialData => "true"}} if ($d =~ m/$r/)},
    "^nodeBias\$" => sub {my ($d,$r)=@_;return {add => {outputHaloModelData => "true"}} if ($d =~ m/$r/)},
    "^(disk|spheroid)StellarLuminosity:([^:]+):([^:]+):z([0-9\.]+)\$" => sub {my ($d,$r)=@_;return {append => {luminosityFilter => $2, luminosityType => $3, luminosityRedshift => $4}} if ($d =~ m/$r/)}
    );
our @missingDatasets;
our %parameterList;

sub Open_File {
    my $dataBlock = shift;
    unless ( exists($dataBlock->{'hdf5File'}) ) {
	$dataBlock->{'hdf5File'} = new PDL::IO::HDF5(">".$dataBlock->{'file'});
	$dataBlock->{'hdf5File'};
    }
}

sub Get_Times {
    my $dataBlock       = shift;
    my $outputNumbers   = pdl [];
    my $times           = pdl [];
    my $expansionFactor = pdl [];
    &Open_File($dataBlock);
    my @outputs = sort $dataBlock->{'hdf5File'}->group("Outputs")->groups;
    foreach my $output ( @outputs ) {
	if ( $output =~ m/Output(\d+)/ ) {
	    my $outputNumber = $1;
	    $outputNumbers   = $outputNumbers  ->append($outputNumber);
	    $times           = $times          ->append($dataBlock->{'hdf5File'}->group("Outputs/".$output)->attrGet("outputTime")           );
	    $expansionFactor = $expansionFactor->append($dataBlock->{'hdf5File'}->group("Outputs/".$output)->attrGet("outputExpansionFactor"));
	}
    }
    $dataBlock->{'outputs'}->{'outputNumber'}    = $outputNumbers;
    $dataBlock->{'outputs'}->{'time'}            = $times;
    $dataBlock->{'outputs'}->{'expansionFactor'} = $expansionFactor;
    $dataBlock->{'outputs'}->{'redshift'}        = 1.0/$expansionFactor-1.0;
}

sub Select_Output {
    my $dataBlock = shift;
    my $redshift = $_[0];
    &Get_Times($dataBlock) unless ( exists($dataBlock->{'outputs'}) );
    my $outputs = $dataBlock->{'outputs'};
    my $expansionFactor = pdl 1.0/($redshift+1.0);
    my $tolerance = pdl 1.0e-3;
    my $foundMatch = 0;
    for(my $i=0;$i<nelem($outputs->{'expansionFactor'});++$i) {
	if ( abs($outputs->{'expansionFactor'}->index($i)-$expansionFactor) < $tolerance ) {
	    $dataBlock->{'output'} = $outputs->{'outputNumber'}->index($i);
	    $foundMatch = 1;
	}
    }
    if ( $foundMatch == 0 ) {
	my $redshiftsAvailable = 1.0/$outputs->{'expansionFactor'}-1.0;
	my $message  = "Select_Output(): Unable to find matching redshift.\n";
	$message .= "                 Requested redshift was: ".$redshift."\n";
	$message .= "                 Available redshifts are: ".$redshiftsAvailable."\n";
	$message .= "                 Try adding the requested redshift to the 'outputRedshifts' parameter in Galacticus.\n";
	die($message);
    }
    # Ensure that the data sets available gets reset for this new output.
    delete($dataBlock->{'dataSetsAvailable'});
}

sub Get_History {
    my $dataBlock = shift;
    my @dataNames = @{$_[0]};
    &Open_File($dataBlock);
    foreach my $dataSetName ( @dataNames ) {
	$dataBlock->{'history'}->{$dataSetName} = $dataBlock->{'hdf5File'}->group("globalHistory")->dataset($dataSetName)->get;
    }
}

sub Get_Datasets_Available {
    my $dataBlock = shift;
    unless ( exists($dataBlock->{'dataSetsAvailable'}) ) {
	&Open_File($dataBlock);
	my @dataSets = $dataBlock->{'hdf5File'}->group("Outputs/Output".$dataBlock->{'output'}."/nodeData")->datasets;
	foreach my $dataSet ( @dataSets ) {$dataBlock->{'dataSetsAvailable'}->{$dataSet} = 1};
    }
}

sub Count_Trees {
    my $dataBlock = shift;
    &Open_File($dataBlock);
    my $outputIndex = 1;
    $outputIndex = $dataBlock->{'output'} if ( exists($dataBlock->{'output'}) );
    # Only read the list of available trees if we have yet to do so, or if it was previously read for a different output.
    unless ( exists($dataBlock->{'mergerTreesAvailable'}) && $dataBlock->{'mergerTreesAvailableUID'} == $outputIndex ) {
	my $treesAvailable = $dataBlock->{'hdf5File'}->group("Outputs/Output".$outputIndex)->dataset("mergerTreeIndex")->get;
	@{$dataBlock->{'mergerTreesAvailable'}} = $treesAvailable->list();
	$dataBlock->{'mergerTreesAvailableUID'} = $outputIndex;
    }
}

sub Get_Parameters {
    my $dataBlock = shift;
    unless ( exists($dataBlock->{'parameters'}) ) {
	&Open_File($dataBlock);
	my @parameterNames = $dataBlock->{'hdf5File'}->group("Parameters")->attrs;
	my @parameterValues = $dataBlock->{'hdf5File'}->group("Parameters")->attrGet(@parameterNames);
	for(my $iParameter=0;$iParameter<=$#parameterNames;++$iParameter) {
	    $dataBlock->{'parameters'}->{$parameterNames[$iParameter]} = $parameterValues[$iParameter];
	}
    }
}

sub Get_Dataset {
    my $dataBlock = shift;
    my @dataNames = @{$_[0]};
    my @mergerTrees;
    if ( $dataBlock->{'tree'} eq "all" ) {
	&Count_Trees($dataBlock);
	@mergerTrees = @{$dataBlock->{'mergerTreesAvailable'}};
	my $treeCount = scalar(@mergerTrees);
    } else {
	$mergerTrees[0] = $dataBlock->{'tree'};
    }

    # Open the HDF5 file.
    &Open_File($dataBlock);

    # Extract a list of parameters.
    &Get_Parameters        ($dataBlock);
    # Extract a list of available datasets.
    &Get_Datasets_Available($dataBlock);
    # Determine the range of data to be extracted.
    my ($dataStart, $dataEnd);
    if ( exists($dataBlock->{'dataRange'}) ) {
	$dataStart = ${$dataBlock->{'dataRange'}}[0];
	$dataEnd   = ${$dataBlock->{'dataRange'}}[1];
    } else {
	$dataStart = -1;
	$dataEnd   = -1;
    }

    # Determine if we are to store derived quantities in the HDF5 file.
    my $storeDataSets;
    if ( exists($dataBlock->{'store'}) ) {
    	$storeDataSets = $dataBlock->{'store'};
    	if ( $storeDataSets == 1 ) {
    	    die ("Get_Dataset(): store only allowed if reading full output at present") unless ( $dataStart == -1 && $dataEnd == -1 && $dataBlock->{'tree'} eq "all" );
    	}
    } else {
    	$storeDataSets = 0;
    }

    foreach my $dataSetName ( @dataNames ) {
    	unless ( exists($dataBlock->{'dataSets'}->{$dataSetName}) ) {
     	    if ( exists($dataBlock->{'dataSetsAvailable'}->{$dataSetName}) || $dataSetName eq "volumeWeight" ) {
     		# Dataset exists in the output file, so simply read it.
     		my $data     = pdl [];
     		my $dataTree = pdl [];
                 # Get merger tree indexing information.
     		my $mergerTreeIndex      = $dataBlock->{'hdf5File'}->dataset("Outputs/Output".$dataBlock->{'output'}."/mergerTreeIndex"     )->get;
     		my $mergerTreeStartIndex = $dataBlock->{'hdf5File'}->dataset("Outputs/Output".$dataBlock->{'output'}."/mergerTreeStartIndex")->get;
     		my $mergerTreeCount      = $dataBlock->{'hdf5File'}->dataset("Outputs/Output".$dataBlock->{'output'}."/mergerTreeCount"     )->get;
     		my $mergerTreeWeight     = $dataBlock->{'hdf5File'}->dataset("Outputs/Output".$dataBlock->{'output'}."/mergerTreeWeight"    )->get;
     		foreach my $mergerTree ( @mergerTrees ) {
     		    # Check that this tree contains some nodes at this output. If it does not, skip it.
     		    my $treeIndex      = which($mergerTreeIndex == $mergerTree);
		    if ( nelem($treeIndex) > 1 ) {
			print "Galacticus::HDF5 - Warning: apparent repeated merger tree index - taking the first instance (** could be a PDL bug**)\n";
			$treeIndex = $treeIndex(0:0);
		    }
     		    my $treeStartIndex = $mergerTreeStartIndex->index($treeIndex);
     		    my $treeCount      = $mergerTreeCount     ->index($treeIndex)->squeeze;
     		    my $treeWeight     = $mergerTreeWeight    ->index($treeIndex)->squeeze;
     		    my $treeEndIndex   = $treeStartIndex+$treeCount-1;
     		    if ( $treeCount > 0 ) {
     			if ( $dataSetName eq "volumeWeight" ) {
     			    $data = $data->append($treeWeight*ones($treeCount->list));
     			} else {
     			    # Read the dataset.
     			    my $thisTreeData = $dataBlock->{'hdf5File'}->group("Outputs/Output".$dataBlock->{'output'}."/nodeData")->dataset($dataSetName)->get($treeStartIndex,$treeEndIndex);
     			    # Append the dataset.
     			    $data = $data->append($thisTreeData);
     			}
     			# Append the merger tree index.
     			unless ( exists($dataBlock->{'dataSets'}->{'mergerTreeIndex'}) ) {
     			    $dataTree = $dataTree->append($mergerTree*ones($treeCount->list));	
     			}
     		    }
     		}
     		$dataBlock->{'dataSets'}->{$dataSetName} = $data;
     		undef($data);
     		unless ( exists($dataBlock->{'dataSets'}->{'mergerTreeIndex'}) ) {
     		    $dataBlock->{'dataSets'}->{'mergerTreeIndex'} = $dataTree;
     		    undef($dataTree);
     		}
	     } else {
	     	# Dataset is not present in the output file, search for a match to a derived property.
	     	my $foundMatch = 0;
	    	foreach my $regEx ( keys(%galacticusFunctions) ) {
	    	    if ( $dataSetName =~ m/^$regEx$/ ) {
	    		$foundMatch = 1;
	    		my $getFunc = $galacticusFunctions{$regEx};
	    		&{$getFunc}($dataBlock,$dataSetName);
	    		if ( $storeDataSets == 1 ) {
	    		    my $dataSets = $dataBlock->{'dataSets'};
	    		    my $nodeDataGroup = $dataBlock->{'hdf5File'}->group("Outputs/Output".$dataBlock->{'output'}."/nodeData");
	    		    my $outputDataSet = new PDL::IO::HDF5::Dataset( name    => $dataSetName,
	    								 parent  => $nodeDataGroup,
	    								 fileObj => $dataBlock->{'hdf5File'}
	    			);
	    		    $outputDataSet->set(${$dataSets->{$dataSetName}});

	    		    # Determine if merger tree references need to be written for this model.
	    		    my $createReference;
	    		    if ( exists($dataBlock->{'parameters'}->{'mergerTreeOutputReferences'}) ) {
	    			if ( $dataBlock->{'parameters'}->{'mergerTreeOutputReferences'} eq "true" ) {
	    			    $createReference = 1;
	    			} else {
	    			    $createReference = 0;
	    			}
	    		    } else {
	    			$createReference = 0;
	    		    }

	    		    # Write merger tree references if necessary.
	    		    if ( $createReference == 1 ) {
	    			my $dataIndexStart = 0;
	    			my $mergerTreeIndex = $dataBlock->{'hdf5File'}->dataset("Outputs/Output".$dataBlock->{'output'}."/mergerTreeIndex")->get;
	    			my $mergerTreeCount = $dataBlock->{'hdf5File'}->dataset("Outputs/Output".$dataBlock->{'output'}."/mergerTreeCount")->get;
	    			my $start = pdl [-1];
	    			foreach my $mergerTree ( @mergerTrees ) {
	    			    # Count up the number of entries in this tree.
	    			    my $mergerTreeGroup = $dataBlock->{'hdf5File'}->group("Outputs/Output".$dataBlock->{'output'}."/mergerTree".$mergerTree);
	    			    $start          += 1;
	    			    my $treeIndex       = which($mergerTreeIndex == $mergerTree);
	    			    my $dataCount       = $mergerTreeCount->index($treeIndex)->squeeze;
	    			    $mergerTreeGroup->reference($outputDataSet,$dataSetName,[$dataIndexStart],[$dataCount]);
	    			    $dataIndexStart += $dataCount;
	    			}
	    		    }
	    		}
	    	    }
	    	}		
	    	# Exit if the dataset was not matched.
		unless ( $foundMatch == 1 ) {
		    push(@missingDatasets,$dataSetName);
		    foreach my $regEx ( keys(%propertyKnowledgeBase) ) {
			if ( $dataSetName =~ m/$regEx/ ) {
			    (my $result = $dataSetName) =~ s/$regEx/$propertyKnowledgeBase{$regEx}/;
			    my $solutions = &{$propertyKnowledgeBase{$regEx}}($dataSetName,$regEx);
			    if ( exists($solutions->{'add'}) ) {
				foreach my $parameter ( keys(%{$solutions->{'add'}}) ) {
				    $parameterList{$parameter} = $solutions->{'add'}->{$parameter};
				}
			    }
			    if ( exists($solutions->{'append'}) ) {
				foreach my $parameter ( keys(%{$solutions->{'append'}}) ) {
				    $parameterList{$parameter} .= " " if ( exists($parameterList{$parameter}) );
				    $parameterList{$parameter} .= $solutions->{'append'}->{$parameter};
				}
			    }
			}
		    }
		}		
	     }
	}
    }
    if ( scalar(@missingDatasets) > 0 ) {
	print "Some requested datasets were not found: ".join(" ",@missingDatasets)."\n";
	print "To resolve this issue try adding the following to your input parameter file:\n";
	foreach my $parameter ( keys(%parameterList) ) {
	    print "<parameter>\n";
	    print "  <name>".$parameter."</name>\n";
	    print "  <value>".$parameterList{$parameter}."</value>\n";
	    print "</parameter>\n";
	}
	print "\nThis may not be a complete list of additional parameters required to solve this issue - check carefully what outputs are required and ensure that you have requested these from Galacticus.\n";
	die();
    }
}

sub Reset_Structure {
    # Attempts to remove all data from the data structure, thereby resetting it.
    my $dataBlock = shift;
    
    foreach my $element ( keys(%{$dataBlock}) ) {
	if ( $element eq "hdf5File" ) {
	    undef($dataBlock->{$element});
	} else {
	    if ( UNIVERSAL::isa( $dataBlock->{$element}, "HASH" ) ) {
		foreach my $key ( keys(%{$dataBlock->{$element}}) ) {
		    undef($dataBlock->{$element}->{$key});
		}
	    }
	    undef($dataBlock->{$element});
	}
    }
}

1;
