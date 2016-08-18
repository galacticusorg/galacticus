# Contains a Perl module which implements calculations luminosities of various sub-mm lines.

package Galacticus::SubMmLines;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use Galacticus::HDF5;
use Galacticus::Grasil;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
			       "luminosityCII"   => \&Galacticus::SubMmLines::Get_CIIline,
    );

sub Get_CIIline {
    # The [CII]157.7um line is computed using a fit to the mean relation and scatter in the results of Graciá-Carpio et al. (2011,
    # http://adsabs.harvard.edu/abs/2011ApJ...728L...7G) as shown in their Figure 1. The fit to the trend and scatter was provided
    # by Peter Capak <capak@astro.caltech.edu>. The luminosity is returned in Solar luminosities, and the input infrared
    # luminosity (as computed from the Grasil SED) is expected to be in Solar luminosities.

    # Get the data structure and the dataset name.
    my $dataBlock   = shift;
    my $dataSetName = $_[0];

    # Ensure that we have the total infrared luminosity.
    &Galacticus::HDF5::Get_Dataset($dataBlock,["grasilInfraredLuminosity"]);
    my $dataSets = $dataBlock->{'dataSets'};

    my $logLuminosity = log10($dataSets->{'grasilInfraredLuminosity'});
    my $dispersion = 0.0467*$logLuminosity-0.1712;
    my $scatter    = $dispersion*grandom(nelem($dispersion));

    # Compute the CII line luminosity.
    $dataSets->{$dataSetName} = 10.0**($logLuminosity-0.1738*$logLuminosity-0.7882+$scatter);

}

1;
