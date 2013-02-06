# Contains a Perl module which implements total stellar mass calculations for Galacticus.

package StellarMass;
use strict;
use warnings;
use PDL;
require Galacticus::HDF5;
use Data::Dumper;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "stellarMass"       => \&StellarMass::Get_StellarMass      ,
    "starFormationRate" => \&StellarMass::Get_StarFormationRate
    );

sub Get_StellarMass {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Get available datasets.
    &HDF5::Get_Datasets_Available($model);

    # Decide which datasets to get.
    my @dataSetsRequired = ( "nodeIndex" );
    my @stellarMassComponents;
    push(@stellarMassComponents,"diskMassStellar"    )
	if ( exists($model->{'dataSetsAvailable'}->{'diskMassStellar'    }) );
    push(@stellarMassComponents,"spheroidMassStellar")
	if ( exists($model->{'dataSetsAvailable'}->{'spheroidMassStellar'}) );
    push(@dataSetsRequired,@stellarMassComponents);

    # Get the datasets.
    &HDF5::Get_Dataset($model,\@dataSetsRequired);

    # Sum the stellar masses.
    $model->{'dataSets'}->{$dataSetName} = pdl zeroes(nelem($model->{'dataSets'}->{'nodeIndex'}));
    foreach my $component ( @stellarMassComponents ) {
	$model->{'dataSets'}->{$dataSetName} += $model->{'dataSets'}->{$component};
    }

}

sub Get_StarFormationRate {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Get available datasets.
    &HDF5::Get_Datasets_Available($model);

    # Decide which datasets to get.
    my @dataSetsRequired = ( "nodeIndex" );
    my @starFormationRateComponents;
    push(@starFormationRateComponents,"diskStarFormationRate"    )
	if ( exists($model->{'dataSetsAvailable'}->{'diskStarFormationRate'    }) );
    push(@starFormationRateComponents,"spheroidStarFormationRate")
	if ( exists($model->{'dataSetsAvailable'}->{'spheroidStarFormationRate'}) );
    push(@dataSetsRequired,@starFormationRateComponents);

    # Get the datasets.
    &HDF5::Get_Dataset($model,\@dataSetsRequired);

    # Sum the stellar masses.
    $model->{'dataSets'}->{$dataSetName} = pdl zeroes(nelem($model->{'dataSets'}->{'nodeIndex'}));
    foreach my $component ( @starFormationRateComponents ) {
	$model->{'dataSets'}->{$dataSetName} += $model->{'dataSets'}->{$component};
    }

}

1;
