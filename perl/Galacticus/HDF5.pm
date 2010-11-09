# Contains a Perl module which implements various useful functionality for extracting data from Galacticus HDF5 files.

package HDF5;
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::NiceSlice;
use Data::Dumper;

%galacticusFunctions = ();

my $status = 1;
$status;

sub Get_Times {
    my $dataHash = shift;
    $outputNumbers   = pdl [];
    $times           = pdl [];
    $expansionFactor = pdl [];
    $HDFfile = new PDL::IO::HDF5(">".${$dataHash}{'file'});
    @outputs = sort $HDFfile->group("Outputs")->groups;
    foreach $output ( @outputs ) {
	if ( $output =~ m/Output(\d+)/ ) {$outputNumber = $1};
	$outputNumbers   = $outputNumbers  ->append($outputNumber);
	$times           = $times          ->append($HDFfile->group("Outputs/".$output)->attrGet("outputTime")           );
	$expansionFactor = $expansionFactor->append($HDFfile->group("Outputs/".$output)->attrGet("outputExpansionFactor"));
    }
    ${${${$dataHash}{'outputs'}}{'outputNumber'}}    = $outputNumbers;
    ${${${$dataHash}{'outputs'}}{'time'}}            = $times;
    ${${${$dataHash}{'outputs'}}{'expansionFactor'}} = $expansionFactor;
    ${${${$dataHash}{'outputs'}}{'redshift'}}        = 1.0/$expansionFactor-1.0;
}

sub Select_Output {
    my $dataHash = shift;
    my $redshift = $_[0];
    &Get_Times($dataHash) unless ( exists(${$dataHash}{'outputs'}) );
    $outputs = $dataHash->{'outputs'};
    $expansionFactor = pdl 1.0/($redshift+1.0);
    $tolerance = pdl 1.0e-3;
    $foundMatch = 0;
    for($i=0;$i<nelem(${$outputs->{'expansionFactor'}});++$i) {
	if ( abs(${$outputs->{'expansionFactor'}}->index($i)-$expansionFactor) < $tolerance ) {
	    $dataHash->{'output'} = ${$outputs->{'outputNumber'}}->index($i);
	    $foundMatch = 1;
	}
    }
    die("Select_Output(): Unable to find matching redshift.\n") if ( $foundMatch == 0 );
    # Ensure that the data sets available gets reset for this new output.
    delete(${$dataHash}{'dataSetsAvailable'});
}

sub Get_History {
    my $dataHash = shift;
    my @dataNames = @{$_[0]};
    $HDFfile = new PDL::IO::HDF5(">".${$dataHash}{'file'});
    foreach $dataSetName ( @dataNames ) {
	${${${$dataHash}{'history'}}{$dataSetName}} = $HDFfile->group("globalHistory")->dataset($dataSetName)->get;
    }
}

sub Count_Trees {
    $dataHash = shift;
    unless ( exists(${$dataHash}{'mergerTreesAvailable'}) ) {
	$HDFfile = new PDL::IO::HDF5(">".${$dataHash}{'file'});
	@treesAvailable = $HDFfile->group("Outputs/Output1")->groups;
	undef(@filteredTrees);
	foreach $tree ( @treesAvailable) {
	    if ( $tree =~ m/mergerTree(\d+)/ ) {
		@dataSets = $HDFfile->group("Outputs/Output1")->group($tree)->datasets;
		push(@filteredTrees,$1) if ( $#dataSets >= 0);
	    }
	}
	@{${$dataHash}{'mergerTreesAvailable'}} = sort {$a <=> $b} @filteredTrees;
    }
}

sub Get_Parameters {
    my $dataHash = shift;
    unless ( exists(${$dataHash}{'parameters'}) ) {
	$HDFfile = new PDL::IO::HDF5(">".${$dataHash}{'file'});
	@parameterNames = $HDFfile->group("Parameters")->attrs;
	@parameterValues = $HDFfile->group("Parameters")->attrGet(@parameterNames);
	for($iParameter=0;$iParameter<=$#parameterNames;++$iParameter) {
	    if ( ref($parameterValues[$iParameter]) eq "PDL::Char" ) {
		${${$dataHash}{'parameters'}}{$parameterNames[$iParameter]} = $parameterValues[$iParameter];
	    } else {
		${${$dataHash}{'parameters'}}{$parameterNames[$iParameter]} = $parameterValues[$iParameter];
	    }
	}
    }
}

sub Get_Dataset {
    my $dataHash = shift;
    my @dataNames = @{$_[0]};
    undef(@mergerTrees);
    if ( ${$dataHash}{'tree'} eq "all" ) {
	&Count_Trees($dataHash);
	@mergerTrees = @{${$dataHash}{'mergerTreesAvailable'}};
    } else {
	$mergerTrees[0] = ${$dataHash}{'tree'};
    }

    # Open the HDF5 file.
    $HDFfile = new PDL::IO::HDF5(">".${$dataHash}{'file'});

    # Extract a list of available datasets.
    unless ( exists(${$dataHash}{'dataSetsAvailable'}) ) {
	# Find a merger tree that contains some output.
	foreach $mergerTree ( @mergerTrees ) {
	    if ( $testID = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->dataset("nodeIndex")->IDget ) {	
		@dataSets = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->datasets;
		foreach $dataSet ( @dataSets ) {${${$dataHash}{'dataSetsAvailable'}}{$dataSet} = 1};
		last;
	    }
	}
    }

    # Determine the range of data to be extracted.
    if ( exists(${$dataHash}{'dataRange'}) ) {
	$dataStart = ${${$dataHash}{'dataRange'}}[0];
	$dataEnd = ${${$dataHash}{'dataRange'}}[1];
    } else {
	$dataStart = -1;
	$dataEnd = -1;
    }

    # Determine if we are to store derived quantities in the HDF5 file.
    if ( exists(${$dataHash}{'store'}) ) {
	$storeDataSets = ${$dataHash}{'store'};
	if ( $storeDataSets == 1 ) {
	    die ("Get_Dataset(): store only allowed if reading full output at present") unless ( $dataStart == -1 && $dataEnd == -1 && ${$dataHash}{'tree'} eq "all" );
	}
    } else {
	$storeDataSets = 0;
    }
    
    foreach $dataSetName ( @dataNames ) {
	unless ( exists(${${$dataHash}{'dataSets'}}{$dataSetName}) ) {
	    if ( exists(${${$dataHash}{'dataSetsAvailable'}}{$dataSetName}) || $dataSetName eq "volumeWeight" ) {
		# Dataset exists in the output file, so simply read it.
		$data = pdl [];
		$dataTree = pdl [];
		foreach $mergerTree ( @mergerTrees ) {
		    # Check that this tree contains some nodes at this output. If it does not, skip it.
		    if ( $testID = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->dataset("nodeIndex")->IDget ) {
			if ( $dataSetName eq "volumeWeight" ) {
                            # Get the volume weight.
			    @volumeWeight = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->attrGet($dataSetName);
			    # Append the volumeWeight property once for each galaxy in this tree.
			    $dataSetTemporary = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->dataset("nodeIndex")->get;
			    $nodeCount = nelem($dataSetTemporary);
			    undef($dataSetTemporary);
			    $data = $data->append($volumeWeight[0]*ones($nodeCount));
			} else {
			    # Read the dataset.
			    $thisTreeData = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->dataset($dataSetName)->get;
			    # Append the dataset.
			    $data = $data->append($thisTreeData);
			}
			# Append the merger tree index.
			unless ( exists(${${$dataHash}{'dataSets'}}{'mergerTreeIndex'}) ) {
			    $dataSetTemporary = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->dataset("nodeIndex")->get;
			    $nodeCount = nelem($dataSetTemporary);
			    undef($dataSetTemporary);			   
			    $dataTree = $dataTree->append($mergerTree*ones($nodeCount));	
			}
		    }
		}
		${${${$dataHash}{'dataSets'}}{$dataSetName}} = $data;
		undef($data);
		unless ( exists(${${$dataHash}{'dataSets'}}{'mergerTreeIndex'}) ) {
		    ${${${$dataHash}{'dataSets'}}{'mergerTreeIndex'}} = $dataTree;
		    undef($dataTree);
		}
	    } else {
		# Dataset is not present in the output file, search for a match to a derived property.
		$foundMatch = 0;
		foreach $regEx ( keys(%galacticusFunctions) ) {
		    if ( $dataSetName =~ m/^$regEx$/ ) {
			$foundMatch = 1;
			$getFunc = $galacticusFunctions{$regEx};
			&$getFunc(\%{$dataHash},$dataSetName);
			if ( $storeDataSets == 1 ) {
			    $dataSets = $dataHash->{'dataSets'};
			    $nodeDataGroup = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/nodeData");
			    $outputDataSet = new PDL::IO::HDF5::Dataset( name    => $dataSetName,
									 parent  => $nodeDataGroup,
									 fileObj => $HDFfile
				);
			    $outputDataSet->set(${$dataSets->{$dataSetName}});
			    $dataIndexStart = 0;
			    foreach $mergerTree ( @mergerTrees ) {
				# Count up the number of entries in this tree.
				$mergerTreeGroup = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree);
				$start = pdl [$mergerTree-1];
				$dataCount = $HDFfile->group("Outputs/Output".${$dataHash}{'output'})->dataset("mergerTreeCount")->get($start,$start);
				$mergerTreeGroup->reference($outputDataSet,$dataSetName,[$dataIndexStart],[$dataCount]);
				$dataIndexStart += $dataCount;
			    }
			}
		    }
		}	
		# Exit if the dataset was not matched.
		die("Dataset ".$dataSetName." was not found or matched to any derived property") unless ( $foundMatch == 1 );		
	    }
	}
    }
}
