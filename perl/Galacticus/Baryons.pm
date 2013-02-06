# Contains a Perl module which implements baryon fraction calculations for Galacticus.

package Baryons;
use PDL;
require Galacticus::HDF5;
use Data::Dumper;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "hotHalo(Fraction|Frac)" => \&Baryons::Get_hotHaloFraction
    );

sub Get_hotHaloFraction {
    $dataSet = shift;
    $dataSetName = $_[0];
    &HDF5::Get_Dataset($dataSet,['hotHaloMass','nodeMass']);
    $dataSets = $dataSet->{'dataSets'};
    $dataSets->{$dataSetName} = $dataSets->{'hotHaloMass'}/$dataSets->{'nodeMass'};
}

1;
