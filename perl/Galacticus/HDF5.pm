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
	$times           = $times          ->append($HDFfile->group("Outputs/".$output)->dataset("outputTime")           ->get);
	$expansionFactor = $expansionFactor->append($HDFfile->group("Outputs/".$output)->dataset("outputExpansionFactor")->get);
    }
    ${${${$dataHash}{'outputs'}}{'outputNumber'}}    = $outputNumbers;
    ${${${$dataHash}{'outputs'}}{'time'}}            = $times;
    ${${${$dataHash}{'outputs'}}{'expansionFactor'}} = $expansionFactor;
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
	foreach $tree ( @treesAvailable) {
	    $tree =~ s/mergerTree(\d+)/$1/;
	}
	@{${$dataHash}{'mergerTreesAvailable'}} = sort {$a <=> $b} @treesAvailable;
    }
}

sub Get_Parameters {
    my $dataHash = shift;
    unless ( exists(${$dataHash}{'parameters'}) ) {
	$HDFfile = new PDL::IO::HDF5(">".${$dataHash}{'file'});
	@parameterNames = $HDFfile->group("Parameters")->datasets;
	foreach $parameterName ( @parameterNames ) {
	    $parameterValue = $HDFfile->group("Parameters")->dataset($parameterName)->get;
	    if ( ref($parameterValue) eq "PDL::Char" ) {
		${${$dataHash}{'parameters'}}{$parameterName} = $parameterValue->atstr(0);
	    } else {
		${${$dataHash}{'parameters'}}{$parameterName} = $parameterValue->index(0);
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
	@dataSets = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTrees[0])->datasets;
	foreach $dataSet ( @dataSets ) {${${$dataHash}{'dataSetsAvailable'}}{$dataSet} = 1};
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
	    if ( exists(${${$dataHash}{'dataSetsAvailable'}}{$dataSetName}) ) {
		$data = pdl [];
		$dataTree = pdl [];
		foreach $mergerTree ( @mergerTrees ) {
		    $thisTreeData = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->dataset($dataSetName)->get;
		    if ( $dataSetName eq "volumeWeight" ) {
			# Append the volumeWeight property once for each galaxy in this tree.
			@galaxyCount = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree)->dataset("nodeIndex")->dims;
			$data = $data->append($thisTreeData*ones($galaxyCount[0]));
		    } else {
			# Append the dataset.
			$data = $data->append($thisTreeData);
		    }
		    unless ( exists(${${$dataHash}{'dataSets'}}{'mergerTreeIndex'}) ) {
			$dataTree = $dataTree->append($mergerTree*ones($galaxyCount[0]));	
		    }
		}
		${${${$dataHash}{'dataSets'}}{$dataSetName}} = $data;
		undef($data);
		unless ( exists(${${$dataHash}{'dataSets'}}{'mergerTreeIndex'}) ) {
		    ${${${$dataHash}{'dataSets'}}{'mergerTreeIndex'}} = $dataTree;
		    undef($dataTree);
		}
	    } else {
		foreach $regEx ( keys(%galacticusFunctions) ) {
		    if ( $dataSetName =~ m/^$regEx$/ ) {
			$getFunc = $galacticusFunctions{$regEx};
			&$getFunc(\%{$dataHash},$dataSetName);
			if ( $storeDataSets == 1 ) {
			    $dataSets = $dataHash->{'dataSets'};
			    $dataIndexStart = 0;
			    foreach $mergerTree ( @mergerTrees ) {
				# Count up the number of entries in this tree. (We do this by counting the nodeIndex property.)
				$mergerTreeGroup = $HDFfile->group("Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree);
				@galaxyCount = $mergerTreeGroup->dataset("nodeIndex")->dims;
				$dataCount = $galaxyCount[0];
				$dataIndexEnd = $dataIndexStart+$dataCount-1;
				$dataSlice = ${$dataSets->{$dataSetName}}->($dataIndexStart:$dataIndexEnd);
				$outputDataSet = new PDL::IO::HDF5::Dataset( name    => $dataSetName,
									     parent  => $mergerTreeGroup,
									     fileObj => $HDFfile
									     );
				$outputDataSet->set($dataSlice);
				$dataIndexStart += $dataCount;
			    }
			}
		    }
		}
	    }
	}
    }
}
