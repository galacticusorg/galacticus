# Contains a Perl module which implements calculation of AGN luminosities.

package AGNLuminosities;
use strict;
use warnings;
use utf8;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::IO::Misc;
use XML::Simple;
use DateTime;
require Galacticus::HDF5;
require Galacticus::ColumnDensity;
require Galacticus::ISMCrossSections;
require Galacticus::Filters;

%HDF5::galacticusFunctions = 
    (
     %HDF5::galacticusFunctions,
     "^agnLuminosity:[^:]+:[^:]+:z[\\d\\.]+(:noAbsorption)?(:alpha[0-9\\-\\+\\.]+)??\$" => \&AGNLuminosities::Get_AGN_Luminosity
    );

# AGN SED data.
our $SEDs;
our $wavelengths;
our $luminositiesBolometric;
our $fileFormatCurrent = pdl long(1);

sub Build_AGN_Spectra {
    # Get file format.
    my $fileFormat = shift;
    die('Build_AGN_Spectra: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
	unless ( $fileFormat == $fileFormatCurrent );

    # Define constants.
    my $speedOfLight        = pdl 2.99792458000e+08; # m/s
    my $angstroms           = pdl 1.00000000000e-10; # m
    my $plancksConstant     = pdl 6.62606800000e-34; # J s.
    my $electronVolt        = pdl 1.60217646000e-19; # J.
    my $luminositySolar     = pdl 3.83900000000e+26; # W
    my $luminosityAB        = pdl 4.46592015000e+13; # W/Hz
    my $massSolar           = pdl 1.98892000000e+30; # kg
    my $gigaYear            = pdl 3.15569260000e+16; # s
    my $kilo                = pdl 1.00000000000e+03;
    my $luminosityABSolar   = $luminosityAB/$luminositySolar;

    # Ensure AGN SED data is loaded.
    unless ( defined($SEDs) ) {

	# Determine whether or not we should (re)make the file.
	my $makeFile = 0;
	if ( -e "data/blackHoles/AGN_SEDs.hdf5" ) {
	    my $hdfFile = new PDL::IO::HDF5("data/blackHoles/AGN_SEDs.hdf5");
	    my @attributes = $hdfFile->attrs();
	    if ( grep {$_ eq "fileFormat"} @attributes  ) {
		my @fileFormatCurrentFile = $hdfFile->attrGet('fileFormat');
		$makeFile = 1 unless ( $fileFormatCurrentFile[0] == $fileFormatCurrent );
	    } else {
		$makeFile = 1;
	    }
	} else {
	    $makeFile = 1;
	}

	# Make the file if necessary.
	if ( $makeFile == 1 ) {

	    # Download the AGN SED code.
	    unless ( -e "aux/AGN_Spectrum/agn_spectrum.c" ) {
		system("mkdir -p aux/AGN_Spectrum; wget --no-check-certificate http://www.cfa.harvard.edu/~phopkins/Site/qlf_files/agn_spectrum.c -O aux/AGN_Spectrum/agn_spectrum.c");
	    }
	    die("Get_AGN_Luminosity(): failed to download agn_spectrum.c")
		unless ( -e "aux/AGN_Spectrum/agn_spectrum.c" );
	    
	    # Compile the AGN SED code.
	    unless ( -e "aux/AGN_Spectrum/agn_spectrum.x" ) {
		system("cd aux/AGN_Spectrum; gcc agn_spectrum.c -o agn_spectrum.x -lm");
	    }
	    die("Get_AGN_Luminosity(): failed to compile agn_spectrum.c")
		unless ( -e "aux/AGN_Spectrum/agn_spectrum.x" );
	    
	    # Generate a tabulation of AGN spectra over a sufficiently large range of AGN luminosity.
	    system("mkdir -p data/blackHoles");
	    my $luminosityBolometricMinimum = pdl  6.0;
	    my $luminosityBolometricMaximum = pdl 28.0;
	    my $luminosityBolometricCount   = 200;
	    $luminositiesBolometric = pdl [];
	    for(my $i=0;$i<$luminosityBolometricCount;++$i) {
		my $luminosityBolometric = ($luminosityBolometricMaximum-$luminosityBolometricMinimum)*$i/($luminosityBolometricCount-1)+$luminosityBolometricMinimum;
		$luminositiesBolometric = $luminositiesBolometric->append($luminosityBolometric);
		my $wavelength = pdl [];
		my $SED        = pdl [];
		open(pHndl,"aux/AGN_Spectrum/agn_spectrum.x ".$luminosityBolometric."|");
		while ( my $line = <pHndl> ) {
		    unless ( $line =~ m/^\s*\;/ ) {
			$line =~ s/^\s*//;
			$line =~ s/\s*$//;
			my @columns = split(/\s+/,$line);
			my $frequency  = $columns[0];
			my $nuLnuSolar = $columns[1];
			if ( $frequency > 0.0 ) {
			    $wavelength = $wavelength->append(
				$speedOfLight/
				(10.0**$frequency)/
				$angstroms
				);
			    $SED        = $SED       ->append(
				(10.0**$nuLnuSolar)/
				(10.0**$frequency)
				);
			}
		    }
		}
		close(pHndl);
		# Construct a PDL to hold the SEDs.
		unless ( defined($SEDs) ) {
		    $SEDs        = pdl zeroes(nelem($SED),$luminosityBolometricCount);
		    $wavelengths = $wavelength(-1:0);
		}
		# Store the SED.
		$SEDs(:,($i)) .= $SED(-1:0);
	    }
	    
	    # Store the data to file.
	    my $hdfFile = new PDL::IO::HDF5(">data/blackHoles/AGN_SEDs.hdf5");
	    my $wavelengthDataSet = new PDL::IO::HDF5::Dataset(
		name    => "wavelength",
		parent  => $hdfFile,
		fileObj => $hdfFile
		);
	    $wavelengthDataSet->set($wavelengths);
	    $wavelengthDataSet->attrSet(
		units     => "Angstroms (Å)",
		unitsInSI => 1.0e-10
		);
	    my $luminosityDataSet = new PDL::IO::HDF5::Dataset(
		name    => "bolometricLuminosity",
		parent  => $hdfFile,
		fileObj => $hdfFile
		);
	    $luminosityDataSet->set($luminositiesBolometric);
	    $luminosityDataSet->attrSet(
		units     => "L☉",
		unitsInSI => 3.827e33
		);
	    my $sedDataSet = new PDL::IO::HDF5::Dataset(
		name    => "SED",
		parent  => $hdfFile,
		fileObj => $hdfFile
		);
	    $sedDataSet->set($SEDs);
	    $sedDataSet->attrSet(
		units     => "L☉/Hz",
		unitsInSI => 3.827e33
		);
	    
	    # Add some metadata.
	    my $dt = DateTime->now->set_time_zone('local');
	    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
	    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
	    $hdfFile->attrSet(
		source       => "Computed using agn_spectrum.c downloaded from https://www.cfa.harvard.edu/~phopkins/Site/qlf.html",
		URL          => "http://adsabs.harvard.edu/abs/2007ApJ...654..731H",
		reference    => "Hopkins et al. (2007)",
		creationTime => $now,
		fileFormat   => $fileFormatCurrent
		);
	} else {
	    # Read the AGN SEDs from file.
	    my $hdfFile             = new PDL::IO::HDF5("data/blackHoles/AGN_SEDs.hdf5");
	    $wavelengths            = $hdfFile->dataset('wavelength'          )->get();
	    $luminositiesBolometric = $hdfFile->dataset('bolometricLuminosity')->get();
	    $SEDs                   = $hdfFile->dataset('SED'                 )->get();
	}
	die("Build_AGN_Spectra(): failed to created data/blackHoles/AGN_SEDs.hdf5")
	    unless ( -e "data/blackHoles/AGN_SEDs.hdf5" );
    }
}

sub Get_AGN_Luminosity {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Define constants.
    my $speedOfLight      = pdl 2.99792458e+08; # m/s
    my $angstroms         = pdl 1.00000000e-10; # m
    my $luminositySolar   = pdl 3.83900000e+26; # W
    my $luminosityAB      = pdl 4.46592015e+13; # W/Hz
    my $massSolar         = pdl 1.98892000e+30; # kg
    my $gigaYear          = pdl 3.15569260e+16; # s
    my $luminosityABSolar = $luminosityAB/$luminositySolar;
    my $plancksConstant   = pdl 6.62606800000e-34; # J s.
    my $kilo              = pdl 1.00000000000e+03;
    my $electronVolt      = pdl 1.60217646000e-19; # J.

    # Ensure spectra file exists.
    &Build_AGN_Spectra($fileFormatCurrent);

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
	&HDF5::Get_Dataset($model,['blackHoleAccretionRate','blackHoleRadiativeEfficiency', 'columnDensityDisk', 'columnDensitySpheroid', 'diskMassGas', 'diskAbundancesGasMetals', 'spheroidMassGas', 'spheroidAbundancesGasMetals' ]);
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
	(my $filterWavelengths, my $filterResponse) = &Filters::Load($filterName);
	
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
	     which($bolometricLuminosity > $luminositiesBolometric(0)),
	     which($bolometricLuminosity < $luminositiesBolometric(-1))
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
		    &ISMCrossSections::Cross_Sections
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
		    &ISMCrossSections::Cross_Sections
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
