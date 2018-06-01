# Contains a Perl module which implements calculation of AGN luminosities.

package Galacticus::AGNLuminosities;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::IO::Misc;
use XML::Simple;
use DateTime;
use Galacticus::HDF5;
use Galacticus::ColumnDensity;
use Galacticus::ISMCrossSections;
use Galacticus::Path;
use Galacticus::Filters;

%Galacticus::HDF5::galacticusFunctions = 
    (
     %Galacticus::HDF5::galacticusFunctions,
     "^agnLuminosity:[^:]+:[^:]+:z[\\d\\.]+(:noAbsorption)?(:alpha[0-9\\-\\+\\.]+)??\$" => \&Galacticus::AGNLuminosities::Get_AGN_Luminosity
    );

sub Get_AGN_Luminosity {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Define constants.
    my $speedOfLight      = pdl 2.99792458e+08; # m/s
    my $angstroms         = pdl 1.00000000e-10; # m
    my $luminositySolar   = pdl 3.84500000e+26; # W
    my $luminosityAB      = pdl 4.46592015e+13; # W/Hz
    my $massSolar         = pdl 1.98892000e+30; # kg
    my $gigaYear          = pdl 3.15569260e+16; # s
    my $luminosityABSolar = $luminosityAB/$luminositySolar;
    my $plancksConstant   = pdl 6.62606800000e-34; # J s.
    my $kilo              = pdl 1.00000000000e+03;
    my $electronVolt      = pdl 1.60217646000e-19; # J.
    my $fileFormatCurrent = pdl 2;
    
    # Ensure spectra file exists.
    system(&galacticusPath()."Galacticus.exe ".&galacticusPath()."parameters/accretion_disks.spectra.Hopkins2007.build_file.xml")
	if ( -e &galacticusPath()."Galacticus.exe" );
    
    # Read the AGN SEDs from file.
    my $hdfFile                        = new PDL::IO::HDF5(&galacticusPath()."data/blackHoles/AGN_SEDs_Hopkins2007.hdf5");
    (my $fileFormat                 )  = $hdfFile->attrGet('fileFormat'          )       ;
    die("Galacticus::AGNLuminosities::Get_AGN_Luminosity: file format is version ".$fileFormat." but was expecting version ".$fileFormatCurrent)
	unless ( $fileFormat == $fileFormatCurrent );
    my $wavelengths                    = $hdfFile->dataset('wavelength'          )->get();
    my $luminositiesBolometricLinear   = $hdfFile->dataset('bolometricLuminosity')->get();
    my $SEDs                           = $hdfFile->dataset('SED'                 )->get();
    my $luminositiesBolometric         = log10($luminositiesBolometricLinear);

    # Determine the filter, frame and redshift for which the luminosity is required.
    if ( $dataSetName =~ m/^agnLuminosity:([^:]+):([^:]+):z([\d\.]+)(:noAbsorption)?(:alpha[0-9\-\+\.]+)??$/ ) {
	# Extract the name of the line and redshift.
	my $filterName   = $1;
	my $frame        = $2;
	my $redshift     = $3;
	my $noAbsorption = $4;
	my $alpha        = $5;
	$alpha =~ s/:alpha//
	    if ( defined($alpha) );
	$noAbsorption = ""
	    unless ( defined($noAbsorption) );

	# Get the AGN bolometric luminosities (in units of Solar luminosities).
	&Galacticus::HDF5::Get_Dataset($model,['blackHoleAccretionRate','blackHoleRadiativeEfficiency', 'columnDensityDisk', 'columnDensitySpheroid', 'diskMassGas', 'diskAbundancesGasMetals', 'spheroidMassGas', 'spheroidAbundancesGasMetals' ]);
	my $dataSets = $model->{'dataSets'};
	my $columnDensityDisk     = $dataSets->{'columnDensityDisk'    };
	my $columnDensitySpheroid = $dataSets->{'columnDensitySpheroid'};
	my $metallicityDisk       = $dataSets->{'diskAbundancesGasMetals'    }/$dataSets->{'diskMassGas'    };
	my $metallicitySpheroid   = $dataSets->{'spheroidAbundancesGasMetals'}/$dataSets->{'spheroidMassGas'};
	my $bolometricLuminosity = 
	    $dataSets->{'blackHoleRadiativeEfficiency'}
	*$dataSets->{'blackHoleAccretionRate'      }
	*$massSolar                                 
	    /$gigaYear                                  
	    *$speedOfLight**2                           
	    /$luminositySolar;
	# Convert to base-10 logarithm of bolometric luminosity.
	my $zeroLuminosities = which($bolometricLuminosity <= 0.0);
	$bolometricLuminosity .= log10($bolometricLuminosity);
	$bolometricLuminosity->index($zeroLuminosities) .= -10.0;

	# Load the filter.
	(my $filterWavelengths, my $filterResponse) = &Galacticus::Filters::Load($filterName);
	
	# Make a joint set of filter and SED wavelengths.
	my $jointWavelengths = $wavelengths->copy();
	$jointWavelengths *= (1.0+$redshift)
	    if ( $frame eq "observed" );
	$jointWavelengths  = $jointWavelengths->append($filterWavelengths);
	my $nonZero        = which(($jointWavelengths >= $filterWavelengths((0))) & ($jointWavelengths <= $filterWavelengths((-1))));
	$jointWavelengths  = $jointWavelengths->index($nonZero)->qsort();

	# Interpolate the filter response onto the joint wavelengths.
	(my $jointResponse, my $interpolateError) = interpolate($jointWavelengths,$filterWavelengths,$filterResponse);

	# Generate a set of delta wavelengths for use in integrations.
	my $deltaWavelengths = pdl [];
	for(my $i=0;$i<nelem($jointWavelengths);++$i) {
	    my $deltaWavelength;
	    if ( $i == 0 ) {
		$deltaWavelength = $jointWavelengths((1))-$jointWavelengths((0));
	    } elsif ( $i == nelem($jointWavelengths)-1 ) {
		$deltaWavelength = $jointWavelengths((-1))-$jointWavelengths((-2));
	    } else {
		$deltaWavelength = $jointWavelengths(($i+1))-$jointWavelengths(($i-1));
	    }
	    $deltaWavelengths = $deltaWavelengths->append($deltaWavelength/2.0);
	}
	
	# If an alpha parameter has been specified, then we're being asked for a broad-band luminosity, defined assuming a
	# particular spectral shape for the AGN spectrum (f_ν∝ν^α) to convert from photon counts to a luminosity. In this case,
	# compute the spectral shape correction factor using this spectrum and include factors to convert from the "AB-magnitude
	# system luminosity zero point" luminosity computed below, to a broad band luminosity in Watts.
	my $spectralShapeCorrection = pdl 1.0;
	if ( defined($alpha) ) {
	    my $photonIntegral = sum($jointResponse*$deltaWavelengths/$jointWavelengths**(1.0+$alpha));
	    my $energyIntegral = sum($jointResponse*$deltaWavelengths/$jointWavelengths**(2.0+$alpha));
	    $spectralShapeCorrection = ($energyIntegral/$photonIntegral)*$luminosityAB*$speedOfLight/$angstroms;
	}

	# Get energies (in keV) corresponding to wavelengths.	
	my $jointEnergies =
	    $speedOfLight
	    *$plancksConstant
	    /$angstroms
	    /$electronVolt
	    /$kilo
	    /$jointWavelengths;
	
	# Integrate SEDs under the filter.	
	my $sedWavelengths = $wavelengths->copy();
	$sedWavelengths *= (1.0+$redshift)
	    if ( $frame eq "observed" );

	# Get SEDs for the joint wavelengths at the model luminosities.
	my $wavelengthSize = $SEDs->getdim(1);
	(my $jointSED, $interpolateError) = interpolate($jointWavelengths(*1),$sedWavelengths(*$wavelengthSize)->xchg(0,1),$SEDs);
	# Determine cross-section model to use.
	my $crossSectionModel = "Wilms2000";
	$crossSectionModel = $model->{'agnLuminosities'}->{'crossSectionsModel'}
	if ( exists($model->{'agnLuminosities'}->{'crossSectionsModel'}) );

	# Construct the AGN luminosities.
	my $luminosityAGN = pdl zeroes(nelem($bolometricLuminosity));
	my $inRange       = 
	    intersect
	    (
	     which($bolometricLuminosity > $luminositiesBolometric(( 0))),
	     which($bolometricLuminosity < $luminositiesBolometric((-1)))
	    );
	$inRange = 
	    intersect
	    (
	     $inRange,
	     $model->{'selection'}
	    )
	    if ( exists($model->{'selection'}) );
	my $nonZeroLuminosity  = $luminosityAGN        ->index($inRange);
	$bolometricLuminosity  = $bolometricLuminosity ->index($inRange);
	$columnDensityDisk     = $columnDensityDisk    ->index($inRange);
	$columnDensitySpheroid = $columnDensitySpheroid->index($inRange);
	$metallicityDisk       = $metallicityDisk      ->index($inRange);
	$metallicitySpheroid   = $metallicitySpheroid  ->index($inRange);
	my $jointSize          = $jointSED->getdim(1);
	my $x                  = $luminositiesBolometric(*$jointSize)->xchg(0,1);
	for (my $i=0;$i<nelem($nonZeroLuminosity);++$i) {
	    (my $interpolatedSED, $interpolateError) = interpolate($bolometricLuminosity($i),$x,$jointSED);
	    my $crossSectionsDisk                    = zeroes(nelem($jointWavelengths));
	    my $crossSectionsSpheroid                = zeroes(nelem($jointWavelengths));
	    unless ( $noAbsorption eq ":noAbsorption" ) {
		my %options = 
		    (
		     model => $crossSectionModel
		    );
		$options{'metallicity'} = $metallicityDisk(($i))
		    if ( 
			exists($model->{'agnLuminosities'}->{'useMetallicity'}) &&
			$model->{'agnLuminosities'}->{'useMetallicity'} == 1
		    );
		$crossSectionsDisk                  .=
		    &Galacticus::ISMCrossSections::Cross_Sections
		    (
		     $jointEnergies,
		     %options
		    );
		delete($options{'metallicity'});
		$options{'metallicity'} = $metallicitySpheroid(($i))
		    if ( 
			exists($model->{'agnLuminosities'}->{'useMetallicity'}) &&
			$model->{'agnLuminosities'}->{'useMetallicity'} == 1
		    );
		$crossSectionsSpheroid              .=
		    &Galacticus::ISMCrossSections::Cross_Sections
		    (
		     $jointEnergies,
		     %options
		    );
	    }
	    my $absorption            = 
		exp
		(
		 -$crossSectionsDisk    *$columnDensityDisk    (($i))
		 -$crossSectionsSpheroid*$columnDensitySpheroid(($i))
		);
	    $nonZeroLuminosity(($i)) .= 
		sum($interpolatedSED  *$jointResponse*$deltaWavelengths*$absorption/$jointWavelengths)/
		sum($luminosityABSolar*$jointResponse*$deltaWavelengths            /$jointWavelengths);
	}
	$luminosityAGN *= $spectralShapeCorrection;
	$dataSets->{$dataSetName} = $luminosityAGN;
    } else {
	die("Get_AGN_Luminosity(): unable to parse data set: ".$dataSetName);
    }
}

1;
