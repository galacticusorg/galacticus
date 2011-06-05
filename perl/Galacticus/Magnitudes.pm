# Contains a Perl module which implements magnitude calculations for Galacticus.

package Magnitudes;
use PDL;
use Galacticus::HDF5;
use Galacticus::DustAttenuation;
use Galacticus::Luminosities;
use Data::Dumper;
use XML::Simple;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^magnitude([^:]+):([^:]+):([^:]+):z([\\d\\.]+)(:dust[^:]+)?(:vega|:AB)?" => "Magnitudes::Get_Magnitude"
    );

my $status = 1;
$status;

sub Get_Magnitude {
    $dataSet = shift;
    $dataSetName = $_[0];
    # Check that the dataset name matches the expected regular expression.
    if ( $dataSetName =~ m/^magnitude([^:]+):([^:]+):([^:]+):z([\d\.]+)(:dust[^:]+)?(:vega|:AB)?/ ) {
	# Extract the dataset name information.
	$component     = $1;
	$filter        = $2;
	$frame         = $3;
	$redshift      = $4;
	$dustExtension = $5;
	if ( $6 eq ":vega" ) {
	    $vegaMagnitude = 1;
	} else {
	    $vegaMagnitude = 0;
	}
	# Construct the name of the corresponding luminosity property.
	$luminosityDataset = lc($component)."StellarLuminosity:".$filter.":".$frame.":z".$redshift.$dustExtension;
	&HDF5::Get_Dataset($dataSet,[$luminosityDataset]);
	$dataSets = \%{${$dataSet}{'dataSets'}};
	${$dataSets->{$dataSetName}} = -2.5*log10(${$dataSets->{$luminosityDataset}}+1.0e-40);
	# If a Vega magnitude was requested, add the appropriate offset.
	if ( $vegaMagnitude == 1 ) {
	    unless ( exists($vegaOffsets{$filter}) ) {
		$filterPath = "./data/filters/".$filter.".xml";
		die("Get_Magnitudes(): can not find filter file for: ".$filter) unless ( -e $filterPath );
		$xml = new XML::Simple;
		$filterData = $xml->XMLin($filterPath);
		unless ( exists($filterData->{'vegaOffset'}) ) {
		    # No Vega offset data available for filter - run the script that computes it.
		    system("scripts/filters/vega_offset_effective_lambda.pl");
		    $filterData = $xml->XMLin($filterPath);
		    die ("Get_Magnitudes(): failed to compute Vega offsets for filters") unless ( exists($filterData->{'vegaOffset'}) );
		}
		$vegaOffsets{$filter} = pdl $filterData->{'vegaOffset'};
	    }
	    ${$dataSets->{$dataSetName}} += $vegaOffsets{$filter};
	}
    } else {
	die("Get_Magnitude(): unable to parse data set: ".$dataSetName);
    }
}
