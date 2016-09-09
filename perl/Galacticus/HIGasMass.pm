# Contains a Perl module which implements total HI gas mass calculations for Galacticus.

package Galacticus::HIGasMass;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use Data::Dumper;
use Galacticus::HDF5;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
    "hiGasMass" => \&Galacticus::HIGasMass::Get_HIGasMass
    );

sub Get_HIGasMass {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Define a constant HI mass fraction (Power, Baugh & Lacey; 2009; http://adsabs.harvard.edu/abs/2009arXiv0908.1396P).
    my $hiMassFraction = pdl 0.54;

    # Get available datasets.
    &Galacticus::HDF5::Get_Datasets_Available($model);

    # Decide which datasets to get.
    my @dataSetsRequired = ( "nodeIndex" );
    my @gasMassComponents;
    push(@gasMassComponents,"diskMassGas"    )
	if ( exists($model->{'dataSetsAvailable'}->{'diskMassGas'    }) );
    push(@gasMassComponents,"spheroidMassGas")
	if ( exists($model->{'dataSetsAvailable'}->{'spheroidMassGas'}) );
    push(@dataSetsRequired,@gasMassComponents);

    # Get the datasets.
    &Galacticus::HDF5::Get_Dataset($model,\@dataSetsRequired);

    # Sum the stellar masses.
    $model->{'dataSets'}->{$dataSetName} = pdl zeroes(nelem($model->{'dataSets'}->{'nodeIndex'}));
    foreach my $component ( @gasMassComponents ) {
	$model->{'dataSets'}->{$dataSetName} += $hiMassFraction*$model->{'dataSets'}->{$component};
    }

}

1;
