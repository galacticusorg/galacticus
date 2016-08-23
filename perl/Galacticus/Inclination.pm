# Contains a Perl module which implements disk inclination calculations for Galacticus.

package Galacticus::Inclination;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use Galacticus::HDF5;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
    "inclination" => \&Galacticus::Inclination::Get_Inclination
    );

sub Get_Inclination {
    my $dataSet = shift;
    my $dataSetName = $_[0];
    &Galacticus::HDF5::Get_Dataset($dataSet,["nodeIndex"]);
    my $dataSets = $dataSet->{'dataSets'};
    $dataSets->{"inclination"} = 180.0*acos(random(nelem($dataSets->{"nodeIndex"})))/3.1415927;
}

1;
