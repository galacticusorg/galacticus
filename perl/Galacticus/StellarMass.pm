# Contains a Perl module which implements total stellar mass calculations for Galacticus.

package Galacticus::StellarMass;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use Galacticus::HDF5;
use Data::Dumper;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
    "massStellar"       => \&Galacticus::StellarMass::Get_StellarMass      ,
    "starFormationRate" => \&Galacticus::StellarMass::Get_StarFormationRate
    );

sub Get_StellarMass {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Get available datasets.
    &Galacticus::HDF5::Get_Datasets_Available($model);

    # Decide which datasets to get.
    my @dataSetsRequired = ( "mergerTreeWeight" );
    my @stellarMassComponents;
    push(@stellarMassComponents,"diskMassStellar"    )
	if ( exists($model->{'dataSetsAvailable'}->{'diskMassStellar'    }) );
    push(@stellarMassComponents,"spheroidMassStellar")
	if ( exists($model->{'dataSetsAvailable'}->{'spheroidMassStellar'}) );
    push(@dataSetsRequired,@stellarMassComponents);

    # Get the datasets.
    &Galacticus::HDF5::Get_Dataset($model,\@dataSetsRequired);

    # Sum the stellar masses.
    $model->{'dataSets'}->{$dataSetName} = pdl zeroes(nelem($model->{'dataSets'}->{'mergerTreeWeight'}));
    foreach my $component ( @stellarMassComponents ) {
	$model->{'dataSets'}->{$dataSetName} += $model->{'dataSets'}->{$component};
    }

}

sub Get_StarFormationRate {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Get available datasets.
    &Galacticus::HDF5::Get_Datasets_Available($model);

    # Decide which datasets to get.
    my @dataSetsRequired = ( "mergerTreeWeight" );
    my @starFormationRateComponents;
    push(@starFormationRateComponents,"diskStarFormationRate"    )
	if ( exists($model->{'dataSetsAvailable'}->{'diskStarFormationRate'    }) );
    push(@starFormationRateComponents,"spheroidStarFormationRate")
	if ( exists($model->{'dataSetsAvailable'}->{'spheroidStarFormationRate'}) );
    push(@dataSetsRequired,@starFormationRateComponents);

    # Get the datasets.
    &Galacticus::HDF5::Get_Dataset($model,\@dataSetsRequired);

    # Sum the stellar masses.
    $model->{'dataSets'}->{$dataSetName} = pdl zeroes(nelem($model->{'dataSets'}->{'mergerTreeWeight'}));
    foreach my $component ( @starFormationRateComponents ) {
	$model->{'dataSets'}->{$dataSetName} += $model->{'dataSets'}->{$component};
    }

}

1;
