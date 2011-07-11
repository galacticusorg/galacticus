# Contains a Perl module which implements calculations of redshift survey quantities for Galacticus.

package Survey;
use PDL;
use Galacticus::HDF5;
use Astro::Cosmology;
use Data::Dumper;
use Switch (Perl6);

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
			       "comovingDistance"   => "Survey::Get_SurveyProperties",
			       "luminosityDistance" => "Survey::Get_SurveyProperties",
			       "redshift"           => "Survey::Get_SurveyProperties",
			       "angularWeight"      => "Survey::Get_SurveyProperties"
    );

my $status = 1;
$status;

sub Get_SurveyProperties {
    # Get the data structure and the dataset name.
    my $dataSet     = shift;
    my $dataSetName = $_[0];

    # Ensure that times and parameters have been read.
    &HDF5::Get_Times     ($dataSet) unless ( exists($dataSet->{'outputs'   }) );
    &HDF5::Get_Parameters($dataSet) unless ( exists($dataSet->{'parameters'}) );

    # Initialize a cosmology.
    my $cosmology = Astro::Cosmology->new(
	omega_matter => $dataSet->{'parameters'}->{'Omega_0' },
	omega_lambda => $dataSet->{'parameters'}->{'Lambda_0'},
	H0           => $dataSet->{'parameters'}->{'H_0'     }
	);

    # Get list of redshifts.
    my $redshifts     = ${$dataSet->{'outputs'}->{'redshift'}}->qsort();
    my $redshiftIndex = nelem($redshifts)-$dataSet->{'output'};
    die("Galacticus::Survey requires more than one output time") if ( nelem($redshifts) <= 1 );

    # Find the range of comoving distances available for this output.
    my $redshiftMinimum;
    my $redshiftMaximum;
    if ( $redshiftIndex == 0                 ) {
	$redshiftMinimum = pdl      $redshifts->index($redshiftIndex);
    } else {
	$redshiftMinimum = pdl 0.5*($redshifts->index($redshiftIndex)+$redshifts->index($redshiftIndex-1));
    }
    if ( $redshiftIndex == nelem($redshifts)-1 ) {
	$redshiftMaximum = pdl      $redshifts->index($redshiftIndex);
    } else {
	$redshiftMaximum = pdl 0.5*($redshifts->index($redshiftIndex)+$redshifts->index($redshiftIndex+1));
    }
    my $comovingDistanceMinimum = $cosmology->comov_dist($redshiftMinimum);
    my $comovingDistanceMaximum = $cosmology->comov_dist($redshiftMaximum);

    # Create a look-up table of comoving distance vs. redshift.
    my $pointsPerRedshift = pdl 1000;
    my $points            = pdl int(($redshiftMaximum-$redshiftMinimum)*$pointsPerRedshift)+1;
    my $redshifts = pdl (0..$points)*($redshiftMaximum-$redshiftMinimum)/$points+$redshiftMinimum;
    my $comovingDistances = pdl [];
    for(my $i=0;$i<nelem($redshifts);++$i) {
	$comovingDistances = $comovingDistances->append($cosmology->comov_dist($redshifts->index($i)));
    }

    # Compute the requested dataset.
    given ($dataSetName) {
	when ( "comovingDistance"   ) {
	    # Ensure that we have the "nodeIndex" property.
	    &HDF5::Get_Dataset($dataSet,["nodeIndex"]);
	    my $dataSets = \%{${$dataSet}{'dataSets'}};
	    # Select comoving distances at random for the galaxies.
	    ${$dataSets->{$dataSetName}} = (($comovingDistanceMaximum**3-$comovingDistanceMinimum**3)*random(nelem(${$dataSets->{"nodeIndex"}}))+$comovingDistanceMinimum**3)**(1.0/3.0);
	}
	when ( "redshift"           ) {
	    # Ensure that we have the "comovingDistance" property.
	    &HDF5::Get_Dataset($dataSet,["comovingDistance"]);
	    my $dataSets = \%{${$dataSet}{'dataSets'}};
	    # Compute the redshift for each galaxy.
	    (${$dataSets->{$dataSetName}},my $error) = interpolate(${$dataSets->{"comovingDistance"}},$comovingDistances,$redshifts);

	}
	when ( "luminosityDistance" ) {
	    # Ensure that we have the "comovingDistance" property.
	    &HDF5::Get_Dataset($dataSet,["comovingDistance","redshift"]);
	    my $dataSets = \%{${$dataSet}{'dataSets'}};
	    # Compute the luminosity distance for each galaxy.
	    ${$dataSets->{$dataSetName}} = ${$dataSets->{"comovingDistance"}}*(1.0+${$dataSets->{"redshift"}});

	}
	when ( "angularWeight"      ) {
	    # Ensure that we have the "volumeWeight" property.
	    &HDF5::Get_Dataset($dataSet,["volumeWeight"]);
	    my $dataSets = \%{${$dataSet}{'dataSets'}};
	    # Compute the angular weight for each galaxy.
	    ${$dataSets->{$dataSetName}} = ${$dataSets->{"volumeWeight"}}*($comovingDistanceMaximum**3-$comovingDistanceMinimum**3)/3.0;
	}
    }
}
