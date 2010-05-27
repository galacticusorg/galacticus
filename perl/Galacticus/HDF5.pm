# Contains a Perl module which implements various useful functionality for extracting data from Galacticus HDF5 files.

package HDF5;
use PDL;
use Data::Dumper;

%galacticusFunctions = ();

my $status = 1;
$status;

sub Get_Times {
    my $dataHash = shift;
    $hdf5Path = ${$dataHash}{'file'}."/Outputs";
    open(hPipe,"h5ls -S $hdf5Path |");
    $outputNumbers = pdl [];
    $times = pdl [];
    $expansionFactor = pdl [];
    while ( $line = <hPipe> ) {
	if ( $line =~ m/Output(\d+)/ ) {
	    $output = $1;
	    $outputNumbers = $outputNumbers->append($output);
	    $hdf5Path = ${$dataHash}{'file'}."/Outputs/Output".$output."/outputTime";
	    open(hdPipe,"h5ls -S -d $hdf5Path |");
	    $dline = <hdPipe>;
	    $dline = <hdPipe>;
	    $dline = <hdPipe>;
	    close(hdPipe);
	    $dline =~ s/^\s*//;
	    $dline =~ s/\s*$//;
	    $times   = $times->append($dline);
	    $hdf5Path = ${$dataHash}{'file'}."/Outputs/Output".$output."/outputExpansionFactor";
	    open(hdPipe,"h5ls -S -d $hdf5Path |");
	    $dline = <hdPipe>;
	    $dline = <hdPipe>;
	    $dline = <hdPipe>;
	    close(hdPipe);
	    $dline =~ s/^\s*//;
	    $dline =~ s/\s*$//;
	    $expansionFactor = $expansionFactor->append($dline);
	}
    }
    close(hPipe);
    ${${${$dataHash}{'outputs'}}{'outputNumber'}} = $outputNumbers;
    ${${${$dataHash}{'outputs'}}{'time'}} = $times;
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
    foreach $dataSetName ( @dataNames ) {
	$hdf5Path = ${$dataHash}{'file'}."/globalHistory/".$dataSetName;
	open(hPipe,"h5ls -S -d $hdf5Path |");
	$line = <hPipe>;
	$line = <hPipe>;
	$data = pdl [];
	while ( $line = <hPipe> ) {
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    $data = $data->append($line);
	}
	close(hPipe);
	${${${$dataHash}{'history'}}{$dataSetName}} = $data;
	undef(@data);
    }
}

sub Count_Trees {
    $dataHash = shift;
    unless ( exists(${$dataHash}{'mergerTreesAvailable'}) ) {
	$hdf5Path = ${$dataHash}{'file'}."/Outputs/Output1";
	open(hPipe,"h5ls $hdf5Path |");
	while ( $line = <hPipe> ) {
	    if ( $line =~ m/^mergerTree(\d+)/ ) {
		push(@{${$dataHash}{'mergerTreesAvailable'}},$1);
	    }
	}
	close(hPipe);
    }
}

sub Get_Parameters {
    my $dataHash = shift;
    unless ( exists(${$dataHash}{'parameters'}) ) {
	$hdf5Path = ${$dataHash}{'file'}."/Parameters";
	open(hPipe,"h5ls -S -d $hdf5Path |");
	while ( $line = <hPipe> ) {
	    if ( $line =~ m/^\s*(\S+)\s*Dataset/ ) {
		$parameterName = $1;
		$line = <hPipe>;
		$line = <hPipe>;
		$line =~ s/\"//g;
		$line =~ s/^\s*//;
		$line =~ s/\s*$//;
		${${$dataHash}{'parameters'}}{$parameterName} = $line;
	    }
	}
	close(hPipe);
    }
}

sub Get_Dataset {
    my $dataHash = shift;
    my @dataNames = @{$_[0]};
    undef(@mergerTrees);
    if ( ${$dataHash}{'tree'} eq "all" ) {
	unless ( exists(${$dataHash}{'mergerTreesAvailable'}) ) {
	    $hdf5Path = ${$dataHash}{'file'}."/Outputs/Output".${$dataHash}{'output'};
	    open(hPipe,"h5ls $hdf5Path |");
	    while ( $line = <hPipe> ) {
		if ( $line =~ m/^mergerTree(\d+)/ ) {
		    push(@{${$dataHash}{'mergerTreesAvailable'}},$1);
		}
	    }
	    close(hPipe);
	}
	@mergerTrees = @{${$dataHash}{'mergerTreesAvailable'}};
	@mergerTrees = sort {$a <=> $b} @mergerTrees;
    } else {
	$mergerTrees[0] = ${$dataHash}{'tree'};
    }

    unless ( exists(${$dataHash}{'dataSetsAvailable'}) ) {
	$hdf5Path = ${$dataHash}{'file'}."/Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTrees[0];
	open(hPipe,"h5ls $hdf5Path |");
	while ( $line = <hPipe> ) {
	    if ( $line =~ m/^(\S+)/ ) {
		${${$dataHash}{'dataSetsAvailable'}}{$1} = 1;
	    }
	}
	close(hPipe);
    }

    if ( exists(${$dataHash}{'dataRange'}) ) {
	$dataStart = ${${$dataHash}{'dataRange'}}[0];
	$dataEnd = ${${$dataHash}{'dataRange'}}[1];
    } else {
	$dataStart = -1;
	$dataEnd = -1;
    }

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
		$dataCount = 0;
		foreach $mergerTree ( @mergerTrees ) {
		    if ( $dataSetName eq "volumeWeight" ) {
			$hdf5Path = ${$dataHash}{'file'}."/Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree."/nodeIndex";
			open(hPipe,"h5ls $hdf5Path |");
			$line = <hPipe>;
			if ( $line =~ m/\{(\d+)/ ) {
			    $galaxyCount = $1;
			} else {
			    $galaxyCount = 0;
			}
			close(hPipe);
		    }
		    $hdf5Path = ${$dataHash}{'file'}."/Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree."/".$dataSetName;
		    open(hPipe,"h5ls -S -d $hdf5Path |");
		    $line = <hPipe>;
		    $line = <hPipe>;
		    while ( $line = <hPipe> ) {
			++$dataCount;
			if ( $dataStart <= 0 || ( $dataCount >= $dataStart && $dataCount <= $dataEnd ) ) {
			    $line =~ s/^\s*//;
			    $line =~ s/\s*$//;
			    if ( $dataSetName eq "volumeWeight" && $galaxyCount > 0 ) {
				$data = $data->append($line*ones($galaxyCount));
			    } else {
				$data = $data->append($line);
			    }
			}
		    }
		    close(hPipe);
		}
		${${${$dataHash}{'dataSets'}}{$dataSetName}} = $data;
		undef(@data);
	    } else {
		foreach $regEx ( keys(%galacticusFunctions) ) {
		    if ( $dataSetName =~ m/^$regEx$/ ) {
			$getFunc = $galacticusFunctions{$regEx};
			&$getFunc(\%{$dataHash},$dataSetName);
			if ( $storeDataSets == 1 ) {
			    $tmpConfigFile = "config.tmp";
			    $tmpDataFile = "data.tmp";
			    $dataSets = $dataHash->{'dataSets'};
			    $dataIndex = -1;
			    foreach $mergerTree ( @mergerTrees ) {
				# Count up the number of entries in this tree. (We do this by counting the nodeIndex property.)
				open(tmpData,">".$tmpDataFile);
				$hdf5Path = ${$dataHash}{'file'}."/Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree."/nodeIndex";
				open(hPipe,"h5ls -S -d $hdf5Path |");
				$line = <hPipe>;
				$line = <hPipe>;
				$dataCount = 0;
				while ( $line = <hPipe> ) {
				    ++$dataCount;
				    ++$dataIndex;
				    print tmpData ${$dataSets->{$dataSetName}}->index($dataIndex)."\n";
				}
				close(hPipe);
				close(tmpData);
				open(configFile,">".$tmpConfigFile);
				print configFile "PATH /Outputs/Output".${$dataHash}{'output'}."/mergerTree".$mergerTree."/".$dataSetName."\n";
				print configFile "INPUT-CLASS TEXTFP\n";
				print configFile "RANK 1\n";
				print configFile "DIMENSION-SIZES ".$dataCount."\n";
				print configFile "OUTPUT-CLASS FP\n";
				print configFile "OUTPUT-SIZE 64\n";
				print configFile "OUTPUT-ARCHITECTURE NATIVE\n";
				close(configFile);
				system("h5import $tmpDataFile -c $tmpConfigFile -o ".${$dataHash}{'file'});
				unlink($tmpConfigFile);
				unlink($tmpDataFile);
			    }
			}
		    }
		}
	    }
	}
    }
}
