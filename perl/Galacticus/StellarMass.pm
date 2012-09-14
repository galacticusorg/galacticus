# Contains a Perl module which implements total stellar mass calculations for Galacticus.

package StellarMass;
use strict;
use warnings;
use PDL;
require Galacticus::HDF5;
use Data::Dumper;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "stellarMass" => \&StellarMass::Get_StellarMass
    );

sub Get_StellarMass {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Get available datasets.
    &HDF5::Get_Datasets_Available($model);

    # Decide which datasets to get.
    my @dataSetsRequired = ( "nodeIndex" );
    my @stellarMassComponents;
    push(@stellarMassComponents,"diskStellarMass"    )
	if ( exists($model->{'dataSetsAvailable'}->{'diskStellarMass'    }) );
    push(@stellarMassComponents,"spheroidStellarMass")
	if ( exists($model->{'dataSetsAvailable'}->{'spheroidStellarMass'}) );
    push(@dataSetsRequired,@stellarMassComponents);

    # Get the datasets.
    &HDF5::Get_Dataset($model,\@dataSetsRequired);

    # Sum the stellar masses.
    $model->{'dataSets'}->{$dataSetName} = pdl zeroes(nelem($model->{'dataSets'}->{'nodeIndex'}));
    foreach my $component ( @stellarMassComponents ) {
	$model->{'dataSets'}->{$dataSetName} += $model->{'dataSets'}->{$component};
    }

}

1;
