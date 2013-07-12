# Contains a Perl module which implements calculations luminosities of various sub-mm lines.

package SubMmLines;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
 $ENV{"GALACTICUS_ROOT_V092"} = getcwd()."/";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
require Galacticus::HDF5;
require Galacticus::Grasil;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
			       "luminosityCII"   => \&SubMmLines::Get_CIIline,
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
    &HDF5::Get_Dataset($dataBlock,["grasilInfraredLuminosity"]);
    my $dataSets = $dataBlock->{'dataSets'};

    my $logLuminosity = log10($dataSets->{'grasilInfraredLuminosity'});
    my $dispersion = 0.0467*$logLuminosity-0.1712;
    my $scatter    = $dispersion*grandom(nelem($dispersion));

    # Compute the CII line luminosity.
    $dataSets->{$dataSetName} = 10.0**($logLuminosity-0.1738*$logLuminosity-0.7882+$scatter);

}

1;
