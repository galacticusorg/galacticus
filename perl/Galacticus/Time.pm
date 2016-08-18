# Contains a Perl module which implements cosmic time calculations for Galacticus.

package Galacticus::Time;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use Galacticus::HDF5;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
    "redshift"        => \&Galacticus::Time::Get_Time,
    "expansionFactor" => \&Galacticus::Time::Get_Time,
    "time"            => \&Galacticus::Time::Get_Time
    );

sub Get_Time {
    # Returns redshift, time or expansion factor for all galaxies.
    my $model       = shift;
    my $dataSetName = $_[0];
    &Galacticus::HDF5::Get_Times  ($model              );
    &Galacticus::HDF5::Get_Dataset($model,["nodeIndex"]);
    my $dataSets = $model->{'dataSets'};
    $dataSets->{$dataSetName}  = pdl ones(nelem($dataSets->{"nodeIndex"}));
    $dataSets->{$dataSetName} *= $model->{'outputs'}->{$dataSetName}->(($model->{'outputIndex'}));
}

1;
