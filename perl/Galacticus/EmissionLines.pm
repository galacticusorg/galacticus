# Contains a Perl module which implements calculation of emission line luminosities, based on
# the methodology of Panuzzo et al. (2003; http://adsabs.harvard.edu/abs/2003A%26A...409...99P).

package EmissionLines;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::GSL::INTEG;
use Data::Dumper;
require Cloudy;
require Galacticus::HDF5;
require Galacticus::IonizingContinuua;
require Galacticus::Launch::PBS;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^(disk|spheroid)LineLuminosity:[^:]+(:[^:]+){0,2}:z[\\d\\.]+\$" => \&EmissionLines::Get_Line_Luminosity      ,
    "^totalLineLuminosity:[^:]+(:[^:]+){0,2}:z[\\d\\.]+\$"           => \&EmissionLines::Get_Total_Line_Luminosity
    );

# Module data.
our $emissionLines;
our $energyThermal;
our $lineData;

# Line list.
our %lineList =
    (
     "H  1  6563A" => "balmerAlpha6563",
     "H  1  4861A" => "balmerBeta4861" ,
     "O II  3726A" => "oxygenII3726"   ,
     "O II  3729A" => "oxygenII3729"   ,
     "O  3  4959A" => "oxygenIII4959"  ,
     "O  3  5007A" => "oxygenIII5007"  ,
     "N  2  6584A" => "nitrogenII6584" ,
     "S II  6731A" => "sulfurII6731"   ,
     "S II  6716A" => "sulfurII6716"
    );

# Constants used in these calculations.
our $Pi                          = pdl  3.1415927000000000e+00;
our $centi                       = pdl  1.0000000000000000e-02;
our $erg                         = pdl  1.0000000000000000e-07; # J
our $plancksConstant             = pdl  6.6260680000000000e-34; # J s
our $speedOfLight                = pdl  2.9979245800000000e+08; # m/s
our $angstroms                   = pdl  1.0000000000000000e-10; # m
our $luminositySolar             = pdl  3.8390000000000000e+26; # W
our $luminosityAB                = pdl  4.4659201576470211e+13; # W/Hz
our $massSolar                   = pdl  1.9891000000000000e+30; # kg
our $megaParsec                  = pdl  3.0856775800000000e+22; # m
our $massAtomic                  = pdl  1.6605389200000000e-27; # kg
our $atomicMassHydrogen          = pdl  1.0079400000000000e+00; # amu
our $massFractionHydrogen        = pdl  0.7070000000000000e+00; # Solar composition
our $boltzmannsConstant          = pdl  1.3810000000000000e-23; # J/K
our $electronVolt                = pdl  1.6020000000000000e-19; # J
our $hydrogenOneIonizationEnergy = pdl 13.5990000000000000e+00; # eV
our $heliumOneIonizationEnergy   = pdl 24.5880000000000000e+00; # eV
our $heliumTwoIonizationEnergy   = pdl 54.4180000000000000e+00; # eV
our $oxygenTwoIonizationEnergy   = pdl 35.1180000000000000e+00; # eV

sub Get_Total_Line_Luminosity {
    my $model       = shift;
    my $dataSetName = $_[0];
    (my $diskDataset     = $dataSetName) =~ s/^total/disk/;
    (my $spheroidDataset = $dataSetName) =~ s/^total/spheroid/;
    &HDF5::Get_Dataset($model,[$diskDataset,$spheroidDataset]);
    my $dataSets = $model->{'dataSets'};
    $dataSets->{$dataSetName} =
	+$dataSets->{$diskDataset}
        +$dataSets->{$spheroidDataset};
}

sub Get_Line_Luminosity {
    my $model       = shift;
    my $dataSetName = $_[0];
    # Define interpolants.
    my @interpolant = 
	(
	 'metallicity'                 ,
	 'densityHydrogen'             ,
	 'ionizingFluxHydrogen'        ,
	 'ionizingFluxHeliumToHydrogen',
	 'ionizingFluxOxygenToHydrogen'
	);
    # Define HII region properties.
    my $massHIIRegion     = pdl 7.5e+3; # Solar masses.
    my $durationHIIRegion = pdl 1.0e-3; # Gyr.
    $massHIIRegion     = $model->{'emissionLines'}->{'hiiRegion'}->{'mass'    }
        if ( exists($model->{'emissionLines'}->{'hiiRegion'}->{'mass'    }) );
    $durationHIIRegion = $model->{'emissionLines'}->{'hiiRegion'}->{'lifetime'}
        if ( exists($model->{'emissionLines'}->{'hiiRegion'}->{'lifetime'}) );
    # Read the emission lines file if necessary.
    unless ( defined($emissionLines) ) {
	my $tableFileName = $galacticusPath."data/hiiRegions/emissionLines.hdf5";
	&Generate_Tables($tableFileName,$model)
	    unless ( -e $tableFileName );
	$emissionLines->{'file'} = new PDL::IO::HDF5($tableFileName);
	$emissionLines->{$_} = log10($emissionLines->{'file'}->dataset($_)->get())
	    foreach ( @interpolant );
	@{$emissionLines->{'lineNames'}} = $emissionLines->{'file'}->group('lines')->datasets();
    }
    # Check that the dataset name matches the expected regular expression.
    if ( $dataSetName =~ m/^(disk|spheroid)LineLuminosity:([^:]+)(:[^:]+)??(:[^:]+)??:z([\d\.]+)$/ ) {
	# Extract the name of the line and redshift.
	my $component  = $1;
	my $lineLabel  = $2;
	my $filterName = $3;
	my $frame      = $4;
	my $redshift   = $5;
	# Load the line if necessary.
	unless ( exists($emissionLines->{'lines'}->{$lineLabel}) ) {
	    die("Galacticus::EmissionLines: unable to find line '".$lineLabel."' in database")
		unless ( grep {$_ eq $lineLabel} @{$emissionLines->{'lineNames'}} );
	    $emissionLines->{'lines'}->{$lineLabel}->{'luminosity'} = $emissionLines->{'file'}->group('lines')->dataset($lineLabel)->get();
	    $emissionLines->{'lines'}->{$lineLabel}->{'wavelength'} = $emissionLines->{'file'}->group('lines')->dataset($lineLabel)->attrGet('wavelength');
	}
	# Construct the name of the corresponding luminosity property.
	my @properties = (
	    $component."LymanContinuumLuminosity:z" .$redshift,
	    $component."HeliumContinuumLuminosity:z".$redshift,
	    $component."OxygenContinuumLuminosity:z".$redshift,
	    $component."MassGas"                              ,
	    $component."AbundancesGasMetals"                  ,
	    $component."Radius"                               ,
	    $component."StarFormationRate"
	    );
	&HDF5::Get_Dataset($model,\@properties);
	my $dataSets = $model->{'dataSets'};
	my $properties;
	# Compute the metallicity in Solar units.
	my $hasGas = which($dataSets->{$component."MassGas"} > 0.0);
	$properties->{'metallicity'} = pdl zeroes(nelem($dataSets->{$component."MassGas"}));
	$properties->{'metallicity'}->($hasGas) .=
	    log10(
		+$dataSets->{$component."AbundancesGasMetals"}->($hasGas)
		/$dataSets->{$component."MassGas"            }->($hasGas)
	    );
	# Compute hydrogen density.
	$properties->{'densityHydrogen'} = pdl zeroes(nelem($dataSets->{$component."MassGas"}));
	$properties->{'densityHydrogen'}->($hasGas) .=
	    log10(
		$massSolar
		*$dataSets->{$component."MassGas"}->($hasGas)
		/$dataSets->{$component."Radius" }->($hasGas)**3
		/$megaParsec                                 **3
		*$centi                                      **3
		/4.0
		/$Pi
		/$massAtomic
		/$atomicMassHydrogen
		/$massFractionHydrogen
	    );
	# Compute Lyman continuum luminosity.
	my $hasFlux = which($dataSets->{$component."LymanContinuumLuminosity:z".$redshift} > 0.0);
	$properties->{'ionizingFluxHydrogen'} = pdl zeroes(nelem($dataSets->{$component."MassGas"}));
	$properties->{'ionizingFluxHydrogen'}->($hasFlux) .=
	    log10($dataSets->{$component."LymanContinuumLuminosity:z".$redshift}->($hasFlux))+50.0;
	# Compute helium to Lyman continuum luminosity ratio.
	$properties->{'ionizingFluxHeliumToHydrogen'} = pdl zeroes(nelem($dataSets->{$component."MassGas"}));
	$properties->{'ionizingFluxHeliumToHydrogen'}->($hasFlux) .=
	    log10(
		+$dataSets->{$component."HeliumContinuumLuminosity:z".$redshift}->($hasFlux)
		/$dataSets->{$component."LymanContinuumLuminosity:z" .$redshift}->($hasFlux)
	    );
	# Compute oxygen to Lyman continuum luminosity ratio.
	$properties->{'ionizingFluxOxygenToHydrogen'} = pdl zeroes(nelem($dataSets->{$component."MassGas"}));
	$properties->{'ionizingFluxOxygenToHydrogen'}->($hasFlux) .=
	    log10(
		+$dataSets->{$component."OxygenContinuumLuminosity:z".$redshift}->($hasFlux)
		/$dataSets->{$component."LymanContinuumLuminosity:z" .$redshift}->($hasFlux)
	    );
	# Check whether a raw line luminosity, or the luminosity under a filter is required.
	my $luminosityMultiplier = pdl 1.0;
	if ( defined($filterName) ) {
	    # A filter was specified.
	    $filterName =~ s/^://;
	    # A frame must also be specified in this case.
	    if ( defined($frame) ) {
		$frame =~ s/^://;
		die("Get_Line_Luminosity(): frame must be either 'rest' or 'observed'")
		    unless ( $frame eq "rest" || $frame eq "observed" );
	    } else {
		die("Get_Line_Luminosity(): a frame ('rest' or 'observed') must be specified");
	    }
	    # Load the filter transmission curve.
	    my $filterFile = "./data/filters/".$filterName.".xml";
	    my $xml        = new XML::Simple;
	    my $filter     = $xml->XMLin($filterFile);
	    my $wavelengths = pdl [];
	    my $response    = pdl [];
	    foreach my $datum ( @{$filter->{'response'}->{'datum'}} ) {
		$datum =~ s/^\s*//;
		$datum =~ s/\s*$//;
		my @columns = split(/\s+/,$datum);
		$wavelengths = append($wavelengths,$columns[0]);
		$response    = append($response   ,$columns[1]);
	    }
	    # Check if the line lies within the extent of the filter.
	    my $observedWavelength = $emissionLines->{'lines'}->{$lineLabel}->{'wavelength'};
	    $observedWavelength *= (1.0+$redshift) if ( $frame eq "observed" );
	    if ( $observedWavelength >= $wavelengths((0)) && $observedWavelength <= $wavelengths((-1)) ) {
		# Interpolate the transmission to the line wavelength.
		(my $transmission, my $interpolateError)
		    = interpolate($observedWavelength,$wavelengths,$response);
		# Integrate a zero-magnitude AB source under the filter.
		my $filterLuminosityAB = pdl 0.0;
		for(my $i=0;$i<nelem($wavelengths);++$i) {
		    my $deltaWavelength;
		    if ( $i == 0 ) {
			$deltaWavelength = $wavelengths->index($i+1)-$wavelengths->index($i  );
		    } elsif ( $i == nelem($wavelengths)-1 ) {
			$deltaWavelength = $wavelengths->index($i  )-$wavelengths->index($i-1);
		    } else {
			$deltaWavelength = $wavelengths->index($i+1)-$wavelengths->index($i-1);
		    }
		    $filterLuminosityAB 
			+=
			+$speedOfLight
			/$angstroms
			*$luminosityAB
			/$luminositySolar
			*$response       ->index($i)
			*$deltaWavelength
			/$wavelengths    ->index($i)**2
			/2.0;
		}
		# Compute the multiplicative factor to convert line luminosity to luminosity in
		# AB units in the filter.
		$luminosityMultiplier .= $transmission/$filterLuminosityAB;
		# Galacticus defines observed-frame luminosities by simply redshifting the
		# galaxy spectrum without changing the amplitude of F_ν (i.e. the compression of
		# the spectrum into a smaller range of frequencies is not accounted for). For a
		# line, we can understand how this should affect the luminosity by considering
		# the line as a Gaussian with very narrow width (such that the full extent of
		# the line always lies in the filter). In this case, when the line is redshifted
		# the width of the Gaussian (in frequency space) is reduced, while the amplitude
		# is unchanged (as, once again, we are not taking into account the compression
		# of the spectrum into the smaller range of frequencies). The integral over the
		# line will therefore be reduced by a factor of (1+z) - this factor is included
		# in the following line. Note that, when converting this observed luminosity
		# into an observed flux a factor of (1+z) must be included to account for
		# compression of photon frequencies (just as with continuum luminosities in
		# Galacticus) which will counteract the effects of the 1/(1+z) included below.
		$luminosityMultiplier /= (1.0+$redshift)
		    if ( $frame eq "observed" );
	    } else {
		# Line lies outside of the filter, so it contributes zero luminosity.
		$luminosityMultiplier .= 0.0;
	    }
	}
	# Find the number of HII regions.
	my $numberHIIRegion   = 
	    +$dataSets->{$component."StarFormationRate"}
	    *$durationHIIRegion
	    /$massHIIRegion;
	# Convert the hydrogen ionizing luminosity to be per HII region.
	$properties->{'ionizingFluxHydrogen'} -= log10($numberHIIRegion);
	# Find interpolation indices for all five interpolants.
	my $index;
	my $error;
	($index->{$_}, $error->{$_}) = 
	    interpolate(
		               $properties   ->{$_}  ,
		               $emissionLines->{$_}  ,
		sequence(nelem($emissionLines->{$_}))
	    )
	    foreach ( @interpolant );
	my $indices = 
	    cat(
		$index->{'metallicity'                 },
		$index->{'densityHydrogen'             },
		$index->{'ionizingFluxHydrogen'        },
		$index->{'ionizingFluxHeliumToHydrogen'},
		$index->{'ionizingFluxOxygenToHydrogen'}
	    )->xchg(0,1);
	# Interpolate in all five parameters to get the line luminosity per HII region.
	my $lineLuminosityPerHIIRegion = 
	    $emissionLines->{'lines'}->{$lineLabel}->{'luminosity'}->interpND($indices);
	# Convert to line luminosity in Solar luminosities (or AB maggies if a filter was specified).
	$dataSets->{$dataSetName} = 
	    +$luminosityMultiplier
	    *$lineLuminosityPerHIIRegion
	    *$erg
	    *$numberHIIRegion
	    /$luminositySolar;
    } else {
	die("Get_Line_Luminosity(): unable to parse data set: ".$dataSetName);
    }
}

sub Generate_Tables {
    # Generate tables of line luminosities using Cloudy.
    my $tableFileName = shift();
    my $model         = shift();
    # Get path to Cloudy.
    (my $cloudyPath, my $cloudyVersion) = &Cloudy::Initialize();
    # Generate look-up table of oxygen-to-helium ionizing photons emission rate ratio as a function of temperature.
    my $integrationRuleOrder         = 6;
    my $integrationToleranceRelative = 1.0e-3;
    my $integrationToleranceAbsolute = 0.0;
    my $intervalsMaximum             = 1000;
    # Iterate over temperature.
    my $temperatureMinimum = pdl   2000.0;
    my $temperatureMaximum = pdl 250000.0;
    my $temperatureStep    = pdl     50.0;
    my $temperatureCount   = int(($temperatureMaximum-$temperatureMinimum)/$temperatureStep)+1;
    my $temperatures       = pdl ($temperatureMaximum-$temperatureMinimum)*sequence($temperatureCount)/($temperatureCount-1)+$temperatureMinimum;
    my $ratiosOxygenHelium = pdl zeroes($temperatureCount);
    for(my $i=0;$i<$temperatureCount;++$i) {
    	$energyThermal = pdl $boltzmannsConstant*$temperatures->(($i))/$electronVolt;
    	(my $oxygenRate, my $oxygenError, my $oxygenStatus) 
    	    = gslinteg_qag
    	    (
    	     \&planckFunction              ,
    	     $oxygenTwoIonizationEnergy   ,
    	     $heliumTwoIonizationEnergy   ,
    	     $integrationToleranceRelative,
    	     $integrationToleranceAbsolute,
    	     $intervalsMaximum            ,
    	     $integrationRuleOrder
    	    );
    	(my $heliumRate, my $heliumError, my $heliumStatus) 
    	    = gslinteg_qag
    	    (
    	     \&planckFunction              ,
    	     $heliumOneIonizationEnergy   ,
    	     $heliumTwoIonizationEnergy   ,
    	     $integrationToleranceRelative,
    	     $integrationToleranceAbsolute,
    	     $intervalsMaximum            ,
    	     $integrationRuleOrder
    	    );
    	$ratiosOxygenHelium->(($i)) .= $oxygenRate/$heliumRate;
    }
    # Find temporary workspace.
    my $workspace = "./emissionLinesWork/";
    $workspace = $model->{'emissionLines'}->{'workspace'}
        if ( exists($model->{'emissionLines'}->{'workspace'}) );
    $workspace .= "/"
	unless ( $workspace =~ m/\/$/ );
    system("mkdir -p ".$workspace);
    # Initialize PBS job stack.
    my @pbsStack;
    # Specify ranges of metallicity, density, and ionizing fluxes.
    my $metallicities           = pdl [ 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 1.0, 2.0 ];
    my $logHydrogenDensities    = pdl sequence(6)        + 1.00;
    my $logHydrogenLuminosities = pdl sequence(9)*4.0/8.0+48.00;
    my $heliumToHydrogenRatios  = pdl sequence(5)*0.05   + 0.05;
    my $oxygenToHydrogenRatios  = pdl sequence(5)*0.06   + 0.06;
    # Initialize line arrays.
    $lineData->{$lineList{$_}}->{'luminosity'} =
	pdl zeroes
	(
	 nelem($oxygenToHydrogenRatios ),
	 nelem($heliumToHydrogenRatios ),
	 nelem($logHydrogenLuminosities),
	 nelem($logHydrogenDensities   ),
	 nelem($metallicities          )
	)
	foreach ( keys(%lineList) );
    # Iterate over ranges.
    my $jobCount = -1;
    my $iLogHydrogenLuminosity = -1;
    foreach my $logHydrogenLuminosity ( $logHydrogenLuminosities->list() ) {
	++$iLogHydrogenLuminosity;
	my $iHeliumToHydrogenRatio = -1;
	foreach my $heliumToHydrogenRatio ( $heliumToHydrogenRatios->list() ) {
	    ++$iHeliumToHydrogenRatio;
	    my $iOxygenToHydrogenRatio = -1;
	    foreach my $oxygenToHydrogenRatio ( $oxygenToHydrogenRatios->list() ) {
		++$iOxygenToHydrogenRatio;
		# Find helium-region temperature.
		my $ratioOxygenHelium = $oxygenToHydrogenRatio/$heliumToHydrogenRatio;
		(my $temperatureHelium, my $temperatureError) = interpolate($ratioOxygenHelium,$ratiosOxygenHelium,$temperatures);
		$temperatureHelium = $temperatureMaximum
		    if ( $temperatureHelium > $temperatureMaximum );
		# Find helium-region normalization.
		$energyThermal = $boltzmannsConstant*$temperatureHelium/$electronVolt;
		(my $normalizationHelium, my $normalizationHeliumError, my $normalizationHeliumStatus) 
		    = gslinteg_qag
		    (
		     \&planckFunction             ,
		     $heliumOneIonizationEnergy   ,
		     $heliumTwoIonizationEnergy   ,
		     $integrationToleranceRelative,
		     $integrationToleranceAbsolute,
		     $intervalsMaximum            ,
		     $integrationRuleOrder
		    );
		$normalizationHelium = $heliumToHydrogenRatio/$normalizationHelium;
		# Find hydrogen-region temperature.
		my $temperatureHydrogen;
		if ( $heliumToHydrogenRatio >= 0.005 ) {
		    $temperatureHydrogen = 3.0e4+4.0e4*$heliumToHydrogenRatio;
		} else {
		    $temperatureHydrogen = 4.0e4+0.5e4*log10($heliumToHydrogenRatio);
		}
		# Find hydrogen-region normalization.
		$energyThermal = $boltzmannsConstant*$temperatureHydrogen/$electronVolt;
		(my $normalizationHydrogen, my $normalizationHydrogenError, my $normalizationHydrogenStatus) 
		    = gslinteg_qag
		    (
		     \&planckFunction             ,
		     $hydrogenOneIonizationEnergy ,
		     $heliumOneIonizationEnergy   ,
		     $integrationToleranceRelative,
		     $integrationToleranceAbsolute,
		     $intervalsMaximum            ,
		     $integrationRuleOrder
		    );
		$normalizationHydrogen = (1.0-$heliumToHydrogenRatio)/$normalizationHydrogen;
		# Get non-ionizing temperature and normalization.
		my $temperatureNonIonizing;
		my $normalizationNonIonizing;
		if ( $heliumToHydrogenRatio >= 0.005 ) {
		    if ( $temperatureHydrogen > 40000.0 ) {
			$temperatureNonIonizing   = 40000.0;
		    } else {
			$temperatureNonIonizing   = $temperatureHydrogen;
		    }
		    $normalizationNonIonizing = 2.5*$normalizationHydrogen;
		} else {
		    $temperatureNonIonizing = 30000.0;
		    $normalizationNonIonizing = 3.5*$normalizationHydrogen;
		}
		my $iMetallicity = -1;
		foreach my $metallicity ( $metallicities->list() ) {
		    ++$iMetallicity;
		    my $iLogHydrogenDensity = -1;
		    foreach my $logHydrogenDensity ( $logHydrogenDensities->list() ) {
			++$iLogHydrogenDensity;
			++$jobCount;
			# Generate Cloudy script.
			my $cloudyScript;
			$cloudyScript .= "title emission line job number ".$jobCount."\n";
			# Generate spectrum table.
			my $counter = -1;
			for(my $logWavelength=11.0;$logWavelength>=2.0;$logWavelength-=0.005) {
			    my $wavelength = 10.0**$logWavelength;
			    my $normalization;
			    my $temperature;
			    if      (                         $wavelength < 227.80 ) {
				$normalization = 0.0;
			    } elsif ( $wavelength >= 227.8 && $wavelength < 504.10 ) {
				$normalization = $normalizationHelium;
				$temperature   = $temperatureHelium;
			    } elsif ( $wavelength >= 504.1 && $wavelength < 911.76 ) {
				$normalization = $normalizationHydrogen;
				$temperature   = $temperatureHydrogen;
			    } else                                           {
				$normalization = $normalizationNonIonizing;
				$temperature    =$temperatureNonIonizing;
			    }
			    if ( $normalization > 0.0 ) {
				++$counter;
				if ( $counter % 3 == 0 ) {
				    if ( $counter == 0 ) {
					$cloudyScript .= "interpolate";
				    } else {
					$cloudyScript .= "\ncontinue";
				    }
				}
				my $energy          = $plancksConstant*$speedOfLight/$wavelength/$angstroms/$electronVolt;
				my $energyRydbergs  = $energy/$hydrogenOneIonizationEnergy;
				$energyThermal      = $boltzmannsConstant*$temperature/$electronVolt;
				my $luminosity      = $normalization*$energy**3/(exp($energy/$energyThermal)-1.0);
				my $logLuminosity   = log10($luminosity);
				$cloudyScript      .= " (".$energyRydbergs." ".$logLuminosity.")";
			
			    }
			}
			$cloudyScript .= "\n";
			$cloudyScript .= "q(h) = ".$logHydrogenLuminosity."\n";
			$cloudyScript .= "metals ".$metallicity          ."\n";
			$cloudyScript .= "hden   ".$logHydrogenDensity   ."\n";
			$cloudyScript .= "sphere expanding\n";
			$cloudyScript .= "radius 16.0\n";
			$cloudyScript .= "stop temperature 1000 k\n";
			$cloudyScript .= "iterate to convergence\n";
			# Write save location.
			$cloudyScript .= "print lines faint _off\n";
			$cloudyScript .= "save lines, array \"".$workspace."lines".$jobCount.".out\"";
			# Write the Cloudy script to file.
			my $cloudyScriptFileName = $workspace."cloudyInput".$jobCount.".txt";
			open(my $cloudyScriptFile,">".$cloudyScriptFileName);
			print $cloudyScriptFile $cloudyScript;
			close($cloudyScriptFile);
			# Generate PBS job.
			my %pbsJob = 
			    (
			     launchFile   => $workspace."emissionLines".$jobCount.".pbs",
			     label        =>            "emissionLines".$jobCount       ,
			     logFile      => $workspace."emissionLines".$jobCount.".log",
			     command      => 
			                     "ulimit -c 0\n"                                               .
       			                     $cloudyPath."/source/cloudy.exe < ".$cloudyScriptFileName."\n".
			                     "if [ $? != 0 ]; then\necho CLOUDY FAILED\nfi\n"              ,
			     ppn          => 1,
			     walltime     => "1:00:00",
			     onCompletion => 
			     {
				 function  => \&linesParse,
				 arguments => 
				     [
				      $workspace."lines"        .$jobCount.".out",
				      $workspace."emissionLines".$jobCount.".pbs",
				      $workspace."emissionLines".$jobCount.".log",
				      $cloudyScriptFileName                      ,
				      $iLogHydrogenLuminosity                    ,
				      $iHeliumToHydrogenRatio                    ,
				      $iOxygenToHydrogenRatio                    ,
				      $iMetallicity                              ,
				      $iLogHydrogenDensity
				     ]
			     }
			    );
			# Push job to PBS queue stack.
			push(@pbsStack,\%pbsJob);
		    }
		}
	    }
	}
    }
    # Submit jobs.
    my %arguments =
	(
	 pbsJobMaximum => 250
	);
    &PBS::SubmitJobs(\%arguments,@pbsStack);
    # Write the line data to file.
    my $tableFile = new PDL::IO::HDF5(">".$tableFileName);
    # Write parameter grid points.
    $tableFile->dataset('metallicity'                 )->set(      $metallicities          );
    $tableFile->dataset('densityHydrogen'             )->set(10.0**$logHydrogenDensities   );
    $tableFile->dataset('ionizingFluxHydrogen'        )->set(10.0**$logHydrogenLuminosities);
    $tableFile->dataset('ionizingFluxHeliumToHydrogen')->set(      $heliumToHydrogenRatios );
    $tableFile->dataset('ionizingFluxOxygenToHydrogen')->set(      $oxygenToHydrogenRatios );
    # Write parameter grid point attributes.
    $tableFile->dataset('metallicity'                 )->attrSet(description => "Metallicity relative to Solar."                                    );
    $tableFile->dataset('densityHydrogen'             )->attrSet(description => "Total hydrogen number density."                                    );
    $tableFile->dataset('densityHydrogen'             )->attrSet(units       => "cm⁻³"                                                              );
    $tableFile->dataset('densityHydrogen'             )->attrSet(unitsInSI   => 1.0e6                                                               );
    $tableFile->dataset('ionizingFluxHydrogen'        )->attrSet(description => "Hydrogen ionizing photon emission rate."                           );
    $tableFile->dataset('ionizingFluxHydrogen'        )->attrSet(units       => "photons/s"                                                         );
    $tableFile->dataset('ionizingFluxHydrogen'        )->attrSet(unitsInSI   => 1.0                                                                 );
    $tableFile->dataset('ionizingFluxHeliumToHydrogen')->attrSet(description => "Helium ionizing photon emission rate relative to that of hydrogen.");
    $tableFile->dataset('ionizingFluxOxygenToHydrogen')->attrSet(description => "Oxygen ionizing photon emission rate relative to that of hydrogen.");
    # Write line data.
    my $lineGroup = $tableFile->group('lines');
    foreach ( keys(%lineList) ) {
	my $lineName = $lineList{$_};
	$lineGroup->dataset($lineName)->set($lineData->{$lineName}->{'luminosity'});
	$lineGroup->dataset($lineName)->attrSet(description => "Luminosity of the line."             );
	$lineGroup->dataset($lineName)->attrSet(units       => "erg/s"                               );
	$lineGroup->dataset($lineName)->attrSet(unitsInSI   => 1.0e-7                                );
	$lineGroup->dataset($lineName)->attrSet(wavelength  => $lineData->{$lineName}->{'wavelength'});

    }
}

sub planckFunction {
    # Blackbody spectrum for integration.
    my ($energy) = @_;
    return $energy**2/(exp($energy/$energyThermal)-1.0);
}

sub linesParse {
    # Parse output from a Cloudy job to extract line data.
    my $linesFileName          = shift();
    my $pbsLaunchFileName      = shift();
    my $pbsLogFileName         = shift();
    my $cloudyScriptFileName   = shift();
    my $iLogHydrogenLuminosity = shift();
    my $iHeliumToHydrogenRatio = shift();
    my $iOxygenToHydrogenRatio = shift();
    my $iMetallicity           = shift();
    my $iLogHydrogenDensity    = shift();
    # Check for successful completion.
    system("grep -q FAILED ".$pbsLogFileName);
    die("EmissionLines::linesParse(): Cloudy failed - see ".$pbsLogFileName)
	if ( $? == 0 );
    # Allow multiple attempts to read the lines file in case the file is written over NFS and we must wait for it to catch up.
    my $attemptsMaximum = 10;
    for(my $attempt=0;$attempt<$attemptsMaximum;++$attempt) {
	# Read the lines file.
	my $badFile = 0;
	open(my $linesFile,$linesFileName);
	while ( my $line = <$linesFile> ) {
	    unless ( $line =~ m/^#/ ) {
		my @columns        = split(/\t/,$line);
		$badFile = 1
		    unless ( scalar(@columns) == 5 );
		my $lineEnergy     = $columns[0];
		my $lineLabel      = $columns[1];
		my $lineLuminosity = $columns[3];
		if ( exists($lineList{$lineLabel}) ) {
		    my $lineName = $lineList{$lineLabel};
		    $lineLuminosity =~ s/^\s+//;
		    $lineLuminosity =~ s/\s+$//;
		    $lineEnergy     =~ s/^\s+//;
		    $lineEnergy     =~ s/\s+$//;
		    my $lineWavelength = $plancksConstant*$speedOfLight/$hydrogenOneIonizationEnergy/$electronVolt/$lineEnergy/$angstroms;
		    $lineData
			->{$lineName}
		    ->{'luminosity'}
		    ->(
			($iOxygenToHydrogenRatio),
			($iHeliumToHydrogenRatio),
			($iLogHydrogenLuminosity),
			($iLogHydrogenDensity   ),
			($iMetallicity          )
			)
			.= 10.0**$lineLuminosity;
		    $lineData
			->{$lineName}
		    ->{'wavelength'}
		    = $lineWavelength;
		}
	    }
	}
	close($linesFile);
	if ( $badFile == 0 ) {
	    last;
	} else {
	    if ( $attempt == $attemptsMaximum-1 ) {
		die("EmissionLines::linesParse(): file '".$linesFileName."' seems to be corrupted");
	    } else {
		sleep(10);
	    }
	}
    }
    # Clean up.
    unlink($linesFileName,$pbsLaunchFileName,$pbsLogFileName,$cloudyScriptFileName);
}

1;
