# Contains a Perl module which implements calculations of redshift survey quantities for Galacticus.

package Survey;
use strict;
use warnings;
use PDL;
require Galacticus::HDF5;
use Astro::Cosmology;
use Data::Dumper;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
			       "comovingDistance"   => \&Survey::Get_SurveyProperties,
			       "luminosityDistance" => \&Survey::Get_SurveyProperties,
			       "redshift"           => \&Survey::Get_SurveyProperties,
			       "angularWeight"      => \&Survey::Get_SurveyProperties
    );

sub Get_SurveyProperties {
    # Get the data structure and the dataset name.
    my $dataBlock     = shift;
    my $dataSetName = $_[0];

    # Check if we have lightcone data available.
    my $useLightcone = 0;
    $useLightcone = 1 if ( exists($dataBlock->{'dataSetsAvailable'}->{'lightconeRedshift'}) );

    # If no lightcone data is available, generate cosmology lookup tables.
    my ($comovingDistanceMinimum, $comovingDistanceMaximum, $comovingDistances, $redshifts);
    if ( $useLightcone == 0 ) {
	
	# Ensure that times and parameters have been read.
	&HDF5::Get_Times     ($dataBlock) unless ( exists($dataBlock->{'outputs'   }) );
	&HDF5::Get_Parameters($dataBlock) unless ( exists($dataBlock->{'parameters'}) );
	
	# Initialize a cosmology.
	my $cosmology = Astro::Cosmology->new(
					      omega_matter => $dataBlock->{'parameters'}->{'Omega_Matter' },
					      omega_lambda => $dataBlock->{'parameters'}->{'Omega_DE'},
					      H0           => $dataBlock->{'parameters'}->{'H_0'     }
					      );
	
	# Get list of redshifts.
	my $redshiftsAvailable     = $dataBlock->{'outputs'}->{'redshift'}->qsort();
	my $redshiftIndex = nelem($redshiftsAvailable)-$dataBlock->{'output'};
	die("Galacticus::Survey requires more than one output time") if ( nelem($redshiftsAvailable) <= 1 );
	
	# Find the range of comoving distances available for this output.
	my $redshiftMinimum;
	my $redshiftMaximum;
	if ( $redshiftIndex == 0                 ) {
	    $redshiftMinimum = pdl      $redshiftsAvailable->index($redshiftIndex);
	} else {
	    $redshiftMinimum = pdl 0.5*($redshiftsAvailable->index($redshiftIndex)+$redshiftsAvailable->index($redshiftIndex-1));
	}
	if ( $redshiftIndex == nelem($redshiftsAvailable)-1 ) {
	    $redshiftMaximum = pdl      $redshiftsAvailable->index($redshiftIndex);
	} else {
	    $redshiftMaximum = pdl 0.5*($redshiftsAvailable->index($redshiftIndex)+$redshiftsAvailable->index($redshiftIndex+1));
	}
	$comovingDistanceMinimum = $cosmology->comov_dist($redshiftMinimum);
	$comovingDistanceMaximum = $cosmology->comov_dist($redshiftMaximum);
	
	# Create a look-up table of comoving distance vs. redshift.
	my $pointsPerRedshift = pdl 1000;
	my $points            = pdl int(($redshiftMaximum-$redshiftMinimum)*$pointsPerRedshift)+1;
	$redshifts            = pdl (0..$points)*($redshiftMaximum-$redshiftMinimum)/$points+$redshiftMinimum;
	$comovingDistances    = pdl [];
	for(my $i=0;$i<nelem($redshifts);++$i) {
	    $comovingDistances = $comovingDistances->append($cosmology->comov_dist($redshifts->index($i)));
	}
    }

    # Compute the requested dataset.
    if ( $dataSetName eq "comovingDistance"   ) {
	if ( $useLightcone == 0 ) {
	    # Ensure that we have the "nodeIndex" property.
	    &HDF5::Get_Dataset($dataBlock,["nodeIndex"]);
	    my $dataSets = $dataBlock->{'dataSets'};
	    # Select comoving distances at random for the galaxies.
	    $dataSets->{$dataSetName} = (($comovingDistanceMaximum**3-$comovingDistanceMinimum**3)*random(nelem($dataSets->{"nodeIndex"}))+$comovingDistanceMinimum**3)**(1.0/3.0);
	} else {
	    # Ensure that we have the lightcone position properties.
	    &HDF5::Get_Dataset($dataBlock,["lightconePositionX","lightconePositionY","lightconePositionZ"]);
	    my $dataSets = $dataBlock->{'dataSets'};
	    # Compute comoving distances for galaxies.
	    $dataSets->{$dataSetName} = sqrt($dataSets->{"lightconePositionX"}**2+$dataSets->{"lightconePositionY"}**2+$dataSets->{"lightconePositionZ"}**2);
	}
    }
    if ( $dataSetName eq "redshift"           ) {
	if ( $useLightcone == 0 ) {
	    # Ensure that we have the "comovingDistance" property.
	    &HDF5::Get_Dataset($dataBlock,["comovingDistance"]);
	    my $dataSets = $dataBlock->{'dataSets'};
	    # Compute the redshift for each galaxy.
	    ($dataSets->{$dataSetName},my $error) = interpolate($dataSets->{"comovingDistance"},$comovingDistances,$redshifts);
	} else {
	    # Ensure that we have the "lightconeRedshift" property.
	    &HDF5::Get_Dataset($dataBlock,["lightconeRedshift"]);
	    my $dataSets = $dataBlock->{'dataSets'};
	    # Compute the redshift for each galaxy.
	    $dataSets->{$dataSetName} = $dataSets->{"lightconeRedshift"};
	}
    }
    if ( $dataSetName eq "luminosityDistance" ) {
	# Ensure that we have the "comovingDistance" property.
	&HDF5::Get_Dataset($dataBlock,["comovingDistance","redshift"]);
	my $dataSets = $dataBlock->{'dataSets'};
	# Compute the luminosity distance for each galaxy.
	$dataSets->{$dataSetName} = $dataSets->{"comovingDistance"}*(1.0+$dataSets->{"redshift"});
    }
    if ( $dataSetName eq "angularWeight"      ) {
	if ( $useLightcone == 0 ) {
	    # Ensure that we have the "mergerTreeWeight" property.
	    &HDF5::Get_Dataset($dataBlock,["mergerTreeWeight"]);
	    my $dataSets = $dataBlock->{'dataSets'};
	    # Compute the angular weight for each galaxy.
	    $dataSets->{$dataSetName} = $dataSets->{"mergerTreeWeight"}*($comovingDistanceMaximum**3-$comovingDistanceMinimum**3)/3.0;
	} else {
	    die("Galacticus::Survey: lightcone data was found, but no 'angularWeight' dataset is present");
	}
    }
}

1;
