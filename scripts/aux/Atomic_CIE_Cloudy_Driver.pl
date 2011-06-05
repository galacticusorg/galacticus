#!/usr/bin/env perl
use XML::Simple;
use Data::Dumper;
use Fcntl qw (:flock);

# Driver script for Cloudy.
# Andrew Benson (26-Jan-2010)

# Get arguments.
if ( $#ARGV != 2 ) {die "Usage: Atomic_CIE_Cloudy_Driver.pl <logMetallicityMaximum> <coolingFunctionFile> <ionizationStateFile>"};
$logMetallicityMaximum = $ARGV[0];
$coolingFunctionFile   = $ARGV[1];
$ionizationStateFile   = $ARGV[2];

# Determine if we need to compute cooling functions.
if ( -e $coolingFunctionFile ) {
    $computeCoolingFunctions = 0;
} else {
    $computeCoolingFunctions = 1;
}

# Determine if we need to compute cooling functions.
if ( -e $ionizationStateFile ) {
    $computeIonizationStates = 0;
} else {
    $computeIonizationStates = 1;
}

# Check if we need to do calculations.
if ( $computeCoolingFunctions == 1 || $computeIonizationStates == 1 ) {

    # Open and lock the cooling function and ionization state files.
    if ( $computeCoolingFunctions == 1 ) {
	open(coolingFunctionOutHndl,">".$coolingFunctionFile);
	flock(coolingFunctionOutHndl,LOCK_EX);
    }
    if ( $computeIonizationStates == 1 ) {
	open(ionizationStateOutHndl,">".$ionizationStateFile);
	flock(ionizationStateOutHndl,LOCK_EX);
    }

    # (Logarithmic) temperature range.
    $logTemperatureMinimum = 2.500;
    $logTemperatureMaximum = 9.000;
    $logTemperatureDelta   = 0.025;
    
    # Logarithmic metallicity range (a zero metallicity case is always included too).
    $logMetallicityMinimum = -4.00;
    $logMetallicityDelta   =  0.25;
    
    # Generate metallicities array.
    @logMetallicities = ( -999.0 ); # Proxy for zero metallicity.
    $logMetallicity = $logMetallicityMinimum-$logMetallicityDelta;
    while ( $logMetallicity < $logMetallicityMaximum ) {
	$logMetallicity += $logMetallicityDelta;	
	$logMetallicities[++$#logMetallicities] = $logMetallicity;
    }

    # Specify Solar and primodial helium abundances (as used in Cloudy).
    $heliumAbundancePrimordial = 0.072;
    $heliumAbundanceSolar      = 0.100;
    
    # Specify Cloudy version.
    $cloudyVersion = "c10.00_rc1";

    # Download the code.
    unless ( -e "aux/".$cloudyVersion.".tar.gz" ) {
	print "Cloudy_Driver.pl: downloading Cloudy code.\n";
	system("wget \"http://viewvc.nublado.org/index.cgi/tags/develop/".$cloudyVersion.".tar.gz?root=cloudy&view=tar\" -O aux/".$cloudyVersion.".tar.gz");
	die("Cloudy_Driver.pl: FATAL - failed to download Cloudy code.") unless ( -e "aux/".$cloudyVersion.".tar.gz" );
    }
    
    # Unpack the code.
    unless ( -e "aux/".$cloudyVersion ) {
	print "Cloudy_Driver.pl: unpacking Cloudy code.\n";
	system("tar -x -v -z -C aux -f aux/".$cloudyVersion.".tar.gz");
	die("Cloudy_Driver.pl: FATAL - failed to unpack Cloudy code.") unless ( -e "aux/".$cloudyVersion );
    }
    
    # Build the code.
    unless ( -e "aux/".$cloudyVersion."/source/cloudy.exe" ) {
	print "Cloudy_Driver.pl: compiling Cloudy code.\n";
	system("cd aux/".$cloudyVersion."/source; chmod u=wrx configure.sh capabilities.pl; make");
	die("Cloudy_Driver.pl: FATAL - failed to build Cloudy code.") unless ( -e "aux/".$cloudyVersion."/source/cloudy.exe" );
    }
   
    # Temporary files for Cloudy data.
    $coolingTempFile    = "./cloudy_cooling.tmp";
    $overviewTempFile   = "./cloudy_overview.tmp";
    
    # Counter for number of cooling functions and ionization states tabulated.
    $iCoolingFunction = -1;
    $iIonizationState = -1;

    # Write message.
    print "Computing cooling functions and ionization states using Cloudy (this may take a long time)...\n";
    
    # Loop over metallicities.
    foreach $logMetallicity ( @logMetallicities ) {
	
	# Increment cooling function and ionization state counter.
	++$iCoolingFunction;
	++$iIonizationState;
	
	# Destroy the previous cooling function data.
	undef(@temperatures     );
	undef(@coolingRates     );
	undef(@electronDensities);
	undef(@hiDensities      );
	undef(@hiiDensities     );
	undef(@heiDensities     );
	undef(@heiiDensities    );
	undef(@heiiDensities    );
	
	# Store the metallicity for this cooling function and ionization state.
	${${$coolingFunctions{'coolingFunction'}}[$iCoolingFunction]}{'metallicity'} = $logMetallicity;
	${${$ionizationStates{'ionizationState'}}[$iIonizationState]}{'metallicity'} = $logMetallicity;
	
	# Run Cloudy.
	open(cloudyPipe,"|aux/".$cloudyVersion."/source/cloudy.exe 1> /dev/null");
	print cloudyPipe "print off\n";
	print cloudyPipe "background, z=0\n";            # Use a very low level incident continuum.
	print cloudyPipe "cosmic rays background\n";     # Include cosmic ray background ionization rate.
	print cloudyPipe "stop zone 1\n";                # Stop after a single zone.
	print cloudyPipe "no photoionization\n";         # Do threee iterations to ensure convergence is reached.
	print cloudyPipe "hden 0.0\n";
	if ( $logMetallicity <= -999.0 ) {
	    print cloudyPipe "abundances primordial\n";
	} else {
	    print cloudyPipe "metals _log ".$logMetallicity."\n";
	    # Assume a linear growth of helium abundance with metallicity.
	    $heliumAbundance = $heliumAbundancePrimordial+($heliumAbundanceSolar-$heliumAbundancePrimordial)*(10.0**$logMetallicity);
	    print cloudyPipe "element abundance linear helium ".$heliumAbundance."\n";
	}
	print cloudyPipe "constant temper ".$logTemperatureMinimum." vary\n";
	print cloudyPipe "grid ".$logTemperatureMinimum." to ".$logTemperatureMaximum." step ".$logTemperatureDelta."\n";
	print cloudyPipe "no molecules\n";
	print cloudyPipe "set trim -20\n";
	print cloudyPipe "punch cooling \"".$coolingTempFile."\"\n";
	print cloudyPipe "punch overview \"".$overviewTempFile."\"\n";
	close(cloudyPipe);

	# Extract the cooling rate.
	open(coolHandle,$coolingTempFile);
	$headerLine = <coolHandle>;
	while ( $dataLine = <coolHandle> ) {
	    $separator   = <coolHandle>;
	    @dataColumns = split(/\s+/,$dataLine);
	    $temperature = $dataColumns[1];
	    $heatingRate = $dataColumns[2];
	    $coolingRate = $dataColumns[3];
	    $temperatures[++$#temperatures] = $temperature;
	    $coolingRates[++$#coolingRates] = $coolingRate;
	}
	close(coolHandle);
	unlink($coolingTempFile);

	# Extract the electron and hydrogen density.
	open(overviewHandle,$overviewTempFile);
	$headerLine = <overviewHandle>;
	while ( $dataLine = <overviewHandle> ) {
	    $separator   = <overviewHandle>;
	    @dataColumns = split(/\s+/,$dataLine);
	    $electronDensity = 10.0**$dataColumns[4];
	    $electronDensities[++$#electronDensities] = $electronDensity;
	    $hiDensity       = 10.0**$dataColumns[6];
	    $hiDensities      [++$#hiDensities      ] = $hiDensity;
	    $hiiDensity      = 10.0**$dataColumns[7];
	    $hiiDensities     [++$#hiiDensities     ] = $hiiDensity;
	}
	close(overviewHandle);
	unlink($overviewTempFile);

	# Store cooling function data.
	@{${${$coolingFunctions{'coolingFunction'}}[$iCoolingFunction]}{'coolingRate'}->{'datum'}}     = @coolingRates;
	@{${${$coolingFunctions{'coolingFunction'}}[$iCoolingFunction]}{'temperature'}->{'datum'}}     = @temperatures;

	# Store ionization state data.
	@{${${$ionizationStates{'ionizationState'}}[$iIonizationState]}{'electronDensity'}->{'datum'}} = @electronDensities;
	@{${${$ionizationStates{'ionizationState'}}[$iIonizationState]}{'hiDensity'      }->{'datum'}} = @hiDensities;
	@{${${$ionizationStates{'ionizationState'}}[$iIonizationState]}{'hiiDensity'     }->{'datum'}} = @hiiDensities;
	@{${${$ionizationStates{'ionizationState'}}[$iIonizationState]}{'temperature'    }->{'datum'}} = @temperatures;
	
    }
    
    # Cooling functions:
    # Specify extrapolation methods in temperature.
    ${${${$coolingFunctions{'extrapolation'}}{'temperature'}}[0]}{'limit'}  = "low";
    ${${${$coolingFunctions{'extrapolation'}}{'temperature'}}[0]}{'method'} = "power law";
    ${${${$coolingFunctions{'extrapolation'}}{'temperature'}}[1]}{'limit'}  = "high";
    ${${${$coolingFunctions{'extrapolation'}}{'temperature'}}[1]}{'method'} = "power law";
    # Specify extrapolation methods in metallicity.
    ${${${$coolingFunctions{'extrapolation'}}{'metallicity'}}[0]}{'limit'}  = "low";
    ${${${$coolingFunctions{'extrapolation'}}{'metallicity'}}[0]}{'method'} = "power law";
    ${${${$coolingFunctions{'extrapolation'}}{'metallicity'}}[1]}{'limit'}  = "high";
    ${${${$coolingFunctions{'extrapolation'}}{'metallicity'}}[1]}{'method'} = "power law";
    # Add some description.
    $coolingFunctions{'description'} = "CIE cooling functions computed by Cloudy ".$cloudyVersion;
    ${$coolingFunctions{'units'}}[0] = "Temperature: Kelvin";
    ${$coolingFunctions{'units'}}[1] = "Cooling rate: Lambda(T)/ergs cm^3 s^-1";
  
    # Ionization states:
    # Specify extrapolation methods in temperature.
    ${${${$ionizationStates{'extrapolation'}}{'temperature'}}[0]}{'limit'}  = "low";
    ${${${$ionizationStates{'extrapolation'}}{'temperature'}}[0]}{'method'} = "fixed";
    ${${${$ionizationStates{'extrapolation'}}{'temperature'}}[1]}{'limit'}  = "high";
    ${${${$ionizationStates{'extrapolation'}}{'temperature'}}[1]}{'method'} = "fixed";
    # Specify extrapolation methods in metallicity.
    ${${${$ionizationStates{'extrapolation'}}{'metallicity'}}[0]}{'limit'}  = "low";
    ${${${$ionizationStates{'extrapolation'}}{'metallicity'}}[0]}{'method'} = "fixed";
    ${${${$ionizationStates{'extrapolation'}}{'metallicity'}}[1]}{'limit'}  = "high";
    ${${${$ionizationStates{'extrapolation'}}{'metallicity'}}[1]}{'method'} = "fixed";
    # Add some description.
    $ionizationStates{'description'} = "CIE ionization states computed by Cloudy 08.00";
    ${$ionizationStates{'units'}}[0] = "Temperature: Kelvin";
    ${$ionizationStates{'units'}}[1] = "Densities: by number relative to hydrogen";

    # Output cooling functions to an XML file.
    if ( $computeCoolingFunctions == 1 ) {
	$coolingFunctions = \%coolingFunctions;
	$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"coolingFunctions");
	print coolingFunctionOutHndl $xmlOutput->XMLout($coolingFunctions);
	close(coolingFunctionOutHndl);
    }
    
    # Output ionization states to an XML file.
    if ( $computeIonizationStates == 1 ) {
	$ionizationStates = \%ionizationStates;
	$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"ionizationStates");
	print ionizationStateOutHndl $xmlOutput->XMLout($ionizationStates);
	close(ionizationStateOutHndl);
    }

    # Write message.
    print "...done\n";
    
}
    
exit;
