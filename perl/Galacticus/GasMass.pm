# Contains a Perl module which implements total gas mass calculations for Galacticus.

package GasMass;
use strict;
use warnings;
use PDL;
use Data::Dumper;
require Galacticus::HDF5;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "massColdGas" => \&GasMass::Get_Cold_Gas_Mass,
    );

sub Get_Cold_Gas_Mass {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Get available datasets.
    &HDF5::Get_Datasets_Available($model);

    # Decide which datasets to get.
    my @dataSetsRequired = ( "nodeIndex" );
    my @gasMassComponents;
    push(@gasMassComponents,"diskMassGas"    )
	if ( exists($model->{'dataSetsAvailable'}->{'diskMassGas'    }) );
    push(@gasMassComponents,"spheroidMassGas")
	if ( exists($model->{'dataSetsAvailable'}->{'spheroidMassGas'}) );
    push(@dataSetsRequired,@gasMassComponents);

    # Get the datasets.
    &HDF5::Get_Dataset($model,\@dataSetsRequired);

    # Sum the gas masses.
    $model->{'dataSets'}->{$dataSetName} = pdl zeroes(nelem($model->{'dataSets'}->{'nodeIndex'}));
    foreach my $component ( @gasMassComponents ) {
	$model->{'dataSets'}->{$dataSetName} += $model->{'dataSets'}->{$component};
    }

}

1;
