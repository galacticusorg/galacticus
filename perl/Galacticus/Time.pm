# Contains a Perl module which implements cosmic time calculations for Galacticus.

package Time;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use Galacticus::HDF5;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "redshift"        => \&Time::Get_Time,
    "expansionFactor" => \&Time::Get_Time,
    "time"            => \&Time::Get_Time
    );

sub Get_Time {
    # Returns redshift, time or expansion factor for all galaxies.
    my $model       = shift;
    my $dataSetName = $_[0];
    &HDF5::Get_Times  ($model              );
    &HDF5::Get_Dataset($model,["nodeIndex"]);
    my $dataSets = $model->{'dataSets'};
    $dataSets->{$dataSetName}  = pdl ones(nelem($dataSets->{"nodeIndex"}));
    $dataSets->{$dataSetName} *= $model->{'outputs'}->{$dataSetName}->(($model->{'outputIndex'}));
}

1;
