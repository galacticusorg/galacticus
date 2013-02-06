# Contains a Perl module which implements properties derived from star formation histories.

package SFH;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
require Galacticus::HDF5;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
			       "^peakSFR\$" => \&SFH::Get_Peak_SFR
			       );

sub Get_Peak_SFR {
    my $model       = shift;
    my $dataSetName = shift;

    # Define prefixes.
    my $giga = pdl 1.0e9;

    # Set stellar mass threshold for computation.
    my $stellarMassThreshold = pdl 0.0;
    $stellarMassThreshold = $model->{'SFH'}->{'stellarMassThreshold'}
    if ( exists($model->{'SFH'}->{'stellarMassThreshold'}) );

    # Get required datasets.
    &HDF5::Get_Dataset($model,[
			       "nodeIndex"          ,
			       "mergerTreeIndex"    ,
			       "diskStellarMass"    ,
			       "spheroidStellarMass"
			       ]
		       );
    my $dataSets = $model->{'dataSets'};

    # Find the total stellar mass.
    my $stellarMass = $dataSets->{'diskStellarMass'    }
    +$dataSets->{'spheroidStellarMass'};

    # The the output index.
    my $outputIndex = $model->{'output'};

    # Loop through galaxies and process their SFH.
    my $sfrPeak = pdl [];
    for(my $i=0;$i<nelem($dataSets->{'nodeIndex'});++$i) {

	# Filter on stellar mass.
	if ( $stellarMass->index($i) > $stellarMassThreshold ) {

	    # Extract tree and node indices.
	    my $mergerTreeIndex = $dataSets->{'mergerTreeIndex'}->index($i);
	    my $nodeIndex       = $dataSets->{'nodeIndex'      }->index($i);
	    # Find the star formation histories.
	    my $SFH;
	    foreach my $component ( "disk", "spheroid" ) {
		my $treeGroup   =
		    "starFormationHistories/"     .
		    "Output"     .$outputIndex    .
		    "/mergerTree".$mergerTreeIndex;
		my @treeDatasets = $model->{'hdf5File'}->group($treeGroup)->datasets();
		my $timeDataset = $component."Time".$nodeIndex;
		if ( grep {$_ eq $timeDataset} @treeDatasets ) {
		    $SFH->{$component}->{'time'} = $model->{'hdf5File'}->group($treeGroup)->dataset($timeDataset)->get();
		} else {
		    $SFH->{$component}->{'time'} = pdl (0.0001,0.0002);
		}
		my $sfhDataset  = $component."SFH".$nodeIndex;
		if ( grep {$_ eq $sfhDataset} @treeDatasets ) {
		    $SFH->{$component}->{'sfh' } = $model->{'hdf5File'}->group($treeGroup)->dataset($sfhDataset )->get();
		} else {
		    $SFH->{$component}->{'sfh' } = pdl zeroes(2);
		}
		# Convert times to timesteps.
		my $timeEnd   = $SFH->{$component}->{'time'};
		my $timeBegin = pdl [0.0];
		$timeBegin    = $timeBegin->append($timeEnd(0:nelem($timeEnd)-2));
		$SFH->{$component}->{'timeStep'} = $timeEnd-$timeBegin;
		# Adjust time to bin center.
		$SFH->{$component}->{'time'} .= ($timeBegin+$timeEnd)/2.0;
		# Convert SFH to SFR.
		$SFH->{$component}->{'sfr'} = $SFH->{$component}->{'sfh'}/$SFH->{$component}->{'timeStep'}/$giga;
	    }
	    # Add disk and spheroid rates.
	    my $timesCombined = pdl [];
	    $timesCombined  = $timesCombined->append($SFH->{'disk'    }->{'time'});
	    $timesCombined  = $timesCombined->append($SFH->{'spheroid'}->{'time'});
	    $timesCombined .= $timesCombined->qsort();
	    $timesCombined  = $timesCombined->uniq();
	    (my $sfrDisk    , my $sfrDiskError    ) =
		interpolate($timesCombined,$SFH->{'disk'    }->{'time'},$SFH->{'disk'    }->{'sfr'});
	    (my $sfrSpheroid, my $sfrSpheroidError) = 
		interpolate($timesCombined,$SFH->{'spheroid'}->{'time'},$SFH->{'spheroid'}->{'sfr'});
	    my $sfrCombined = $sfrDisk+$sfrSpheroid;
	    # Find the peak SFR.
	    my $sfrIndex = $sfrCombined->qsorti();
	    $sfrPeak     = $sfrPeak->append($sfrCombined->index($sfrIndex((-1))));
	} else {
	    $sfrPeak     = $sfrPeak->append(-1.0                                );
	}

	# Transfer to the output data structure.
	$model->{'dataSets'}->{"peakSFR"} = $sfrPeak;
    }
}

1;
