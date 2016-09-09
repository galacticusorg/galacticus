# Contains a Perl module which implements calculations of 850um sub-mm fluxes using
# the simple fitting formula from Hayward et al. (2011; http://adsabs.harvard.edu/abs/2011arXiv1101.0002H).

package Galacticus::SubMmFluxesHayward;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use Galacticus::HDF5;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
			       "flux850micronHayward"   => \&Galacticus::SubMmFluxesHayward::Get_850micron,
    );

sub Get_850micron {
    # Get the data structure and the dataset name.
    my $dataBlock   = shift;
    my $dataSetName = $_[0];

    # Specify parameters of the fitting function.
    my $dustToMetalsRatio         = pdl 0.61;
    my $fitNormalization          = pdl 0.65e-3;
    my $starFormationRateExponent = pdl 0.42;
    my $dustMassExponent          = pdl 0.58;
    if ( exists($dataBlock->{'haywardSubMmFit'}) ) {
	$dustToMetalsRatio         = $dataBlock->{'haywardSubMmFit'}->{'dustToMetalsRatio'        }
	if ( exists($dataBlock->{'haywardSubMmFit'}->{'dustToMetalsRatio'         }) );
	$fitNormalization          = $dataBlock->{'haywardSubMmFit'}->{'fitNormalization'         }
	if ( exists($dataBlock->{'haywardSubMmFit'}->{'fitNormalization'          }) );
	$starFormationRateExponent = $dataBlock->{'haywardSubMmFit'}->{'starFormationRateExponent'}
	if ( exists($dataBlock->{'haywardSubMmFit'}->{'starFormationRateExponent' }) );
	$dustMassExponent          = $dataBlock->{'haywardSubMmFit'}->{'dustmassExponent'         }
	if ( exists($dataBlock->{'haywardSubMmFit'}->{'dustmassExponent'          }) );
    }

    # Ensure that we have required datasets.
    &Galacticus::HDF5::Get_Dataset($dataBlock,[
			   "diskStarFormationRate"  ,"spheroidStarFormationRate"  ,
			   "diskAbundancesGasMetals","spheroidAbundancesGasMetals"
		       ]);
    my $dataSets = $dataBlock->{'dataSets'};
    # Compute dust mass and total star formation rate.
    my $starFormationRate = $dataSets->{"diskStarFormationRate"}+$dataSets->{"spheroidStarFormationRate"};
    my $dustMass          = ($dataSets->{"diskAbundancesGasMetals"}+$dataSets->{"spheroidAbundancesGasMetals"})*$dustToMetalsRatio;
    # Select comoving distances at random for the galaxies.
    $dataSets->{$dataSetName} = 
	$fitNormalization
	*($starFormationRate/100.0e9)**$starFormationRateExponent
	*($dustMass         /  1.0e8)**$dustMassExponent;

}

1;
