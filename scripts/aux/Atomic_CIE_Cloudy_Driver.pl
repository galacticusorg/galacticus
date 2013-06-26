#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V092'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V092'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");
use XML::Simple;
use Data::Dumper;
use Fcntl qw (:flock);

# Driver script for Cloudy.
# Andrew Benson (26-Jan-2010)

# Get arguments.
if ( $#ARGV != 3 ) {die "Usage: Atomic_CIE_Cloudy_Driver.pl <logMetallicityMaximum> <coolingFunctionFile> <chemicalStateFile> <fileFormatVersion>"};
my $logMetallicityMaximum = $ARGV[0];
my $coolingFunctionFile   = $ARGV[1];
my $chemicalStateFile     = $ARGV[2];
my $fileFormat            = $ARGV[3];

# Ensure the requested file format version is compatible.
my $fileFormatCurrent = 1;
die('Atomic_CIE_Cloudy_Driver.pl: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
    unless ( $fileFormat == $fileFormatCurrent );

# Determine if we need to compute cooling functions.
my $computeCoolingFunctions;
if ( -e $coolingFunctionFile ) {
    my $xmlDoc = new XML::Simple;
    my $coolingFunction = $xmlDoc->XMLin($coolingFunctionFile);
    if ( exists($coolingFunction->{'fileFormat'}) ) { 
	$computeCoolingFunctions = 1 unless ( $coolingFunction->{'fileFormat'} == $fileFormatCurrent );
    } else {
	$computeCoolingFunctions = 0;
    }
} else {
    $computeCoolingFunctions = 1;
}

# Determine if we need to compute cooling functions.
my $computeChemicalStates;
if ( -e $chemicalStateFile ) {
    my $xmlDoc = new XML::Simple;
    my $chemicalState = $xmlDoc->XMLin($chemicalStateFile);
    if ( exists($chemicalState->{'fileFormat'}) ) { 
	$computeChemicalStates = 1 unless ( $chemicalState->{'fileFormat'} == $fileFormatCurrent );
    } else {
	$computeChemicalStates = 0;
    }
} else {
    $computeChemicalStates = 1;
}

# Check if we need to do calculations.
if ( $computeCoolingFunctions == 1 || $computeChemicalStates == 1 ) {

    # Open and lock the cooling function and chemical state files.
    if ( $computeCoolingFunctions == 1 ) {
	open(coolingFunctionOutHndl,">".$coolingFunctionFile);
	flock(coolingFunctionOutHndl,LOCK_EX);
    }
    if ( $computeChemicalStates == 1 ) {
	open(chemicalStateOutHndl,">".$chemicalStateFile);
	flock(chemicalStateOutHndl,LOCK_EX);
    }

    # (Logarithmic) temperature range.
    my $logTemperatureMinimum = 2.500;
    my $logTemperatureMaximum = 9.000;
    my $logTemperatureDelta   = 0.025;
    
    # Logarithmic metallicity range (a zero metallicity case is always included too).
    my $logMetallicityMinimum = -4.00;
    my $logMetallicityDelta   =  0.25;
    
    # Generate metallicities array.
    my @logMetallicities = ( -999.0 ); # Proxy for zero metallicity.
    my $logMetallicity = $logMetallicityMinimum-$logMetallicityDelta;
    while ( $logMetallicity < $logMetallicityMaximum ) {
	$logMetallicity += $logMetallicityDelta;	
	$logMetallicities[++$#logMetallicities] = $logMetallicity;
    }

    # Specify Solar and primodial helium abundances (as used in Cloudy).
    my $heliumAbundancePrimordial = 0.072;
    my $heliumAbundanceSolar      = 0.100;
    
    # Specify Cloudy version.
    my $cloudyVersion = "c13.02";

    # Download the code.
    unless ( -e "aux/".$cloudyVersion.".tar.gz" ) {
	print "Atomic_CIE_Cloudy_Driver.pl: downloading Cloudy code.\n";
	system("wget \"http://data.nublado.org/cloudy_releases/c13/".$cloudyVersion.".tar.gz\" -O ".$galacticusPath."aux/".$cloudyVersion.".tar.gz");
	die("Atomic_CIE_Cloudy_Driver.pl: FATAL - failed to download Cloudy code.") unless ( -e $galacticusPath."aux/".$cloudyVersion.".tar.gz" );
    }
    
    # Unpack the code.
    unless ( -e $galacticusPath."aux/".$cloudyVersion ) {
	print "Atomic_CIE_Cloudy_Driver.pl: unpacking Cloudy code.\n";
	system("tar -x -v -z -C aux -f ".$galacticusPath."aux/".$cloudyVersion.".tar.gz");
	die("Atomic_CIE_Cloudy_Driver.pl: FATAL - failed to unpack Cloudy code.") unless ( -e $galacticusPath."aux/".$cloudyVersion );
    }
    
    # Build the code.
    unless ( -e $galacticusPath."aux/".$cloudyVersion."/source/cloudy.exe" ) {
	print "Atomic_CIE_Cloudy_Driver.pl: compiling Cloudy code.\n";
	system("cd ".$galacticusPath."aux/".$cloudyVersion."/source; chmod u=wrx configure.sh capabilities.pl; make");
	die("Atomic_CIE_Cloudy_Driver.pl: FATAL - failed to build Cloudy code.") unless ( -e $galacticusPath."aux/".$cloudyVersion."/source/cloudy.exe" );
    }
   
    # Temporary files for Cloudy data.
    my $coolingTempFile    = "cloudy_cooling.tmp";
    my $overviewTempFile   = "cloudy_overview.tmp";
    
    # Counter for number of cooling functions and chemical states tabulated.
    my $iCoolingFunction = -1;
    my $iChemicalState = -1;

    # Write message.
    print "Computing cooling functions and chemical states using Cloudy (this may take a long time)...\n";
    
    # Initialize data structures.
    my %coolingFunctions;
    my %chemicalStates;
    
    # Loop over metallicities.
    foreach my $logMetallicity ( @logMetallicities ) {
	
	# Increment cooling function and chemical state counter.
	++$iCoolingFunction;
	++$iChemicalState;
	
	# Destroy the previous cooling function data.
	my @temperatures     ;
	my @coolingRates     ;
	my @electronDensities;
	my @hiDensities      ;
	my @hiiDensities     ;
	my @heiDensities     ;
	my @heiiDensities    ;
	
	# Store the metallicity for this cooling function and chemical state.
	${${$coolingFunctions{'coolingFunction'}}[$iCoolingFunction]}{'metallicity'} = $logMetallicity;
	${${$chemicalStates{'chemicalState'}}[$iChemicalState]}{'metallicity'} = $logMetallicity;
	
	# Run Cloudy.
	open(cloudyScript,">".$galacticusPath."aux/".$cloudyVersion."/source/input.in");
	print cloudyScript "print off\n";
	print cloudyScript "background, z=0\n";            # Use a very low level incident continuum.
	print cloudyScript "cosmic rays background\n";     # Include cosmic ray background ionization rate.
	print cloudyScript "stop zone 1\n";                # Stop after a single zone.
	print cloudyScript "no photoionization\n";         # Do three iterations to ensure convergence is reached.
	print cloudyScript "hden 0.0\n";
	if ( $logMetallicity <= -999.0 ) {
	    print cloudyScript "abundances primordial\n";
	} else {
	    print cloudyScript "metals _log ".$logMetallicity."\n";
	    # Assume a linear growth of helium abundance with metallicity.
	    my $heliumAbundance = $heliumAbundancePrimordial+($heliumAbundanceSolar-$heliumAbundancePrimordial)*(10.0**$logMetallicity);
	    print cloudyScript "element abundance linear helium ".$heliumAbundance."\n";
	}
	print cloudyScript "constant temper ".$logTemperatureMinimum." vary\n";
	print cloudyScript "grid ".$logTemperatureMinimum." to ".$logTemperatureMaximum." step ".$logTemperatureDelta."\n";
	print cloudyScript "no molecules\n";
	print cloudyScript "set trim -20\n";
	print cloudyScript "punch cooling \"".$coolingTempFile."\"\n";
	print cloudyScript "punch overview \"".$overviewTempFile."\"\n";
	close(cloudyScript);
	system("cd ".$galacticusPath."aux/".$cloudyVersion."/source; cloudy.exe -r input");

	# Extract the cooling rate.
	open(coolHandle,$galacticusPath."aux/".$cloudyVersion."/source/".$coolingTempFile);
	my $headerLine = <coolHandle>;
	while ( my $dataLine = <coolHandle> ) {
	    my $separator   = <coolHandle>;
	    my @dataColumns = split(/\s+/,$dataLine);
	    my $temperature = $dataColumns[1];
	    my $heatingRate = $dataColumns[2];
	    my $coolingRate = $dataColumns[3];
	    $temperatures[++$#temperatures] = $temperature;
	    $coolingRates[++$#coolingRates] = $coolingRate;
	}
	close(coolHandle);
	unlink($coolingTempFile);

	# Extract the electron and hydrogen density.
	open(overviewHandle,$galacticusPath."aux/".$cloudyVersion."/source/".$overviewTempFile);
	$headerLine = <overviewHandle>;
	while ( my $dataLine = <overviewHandle> ) {
	    my $separator   = <overviewHandle>;
	    my @dataColumns = split(/\s+/,$dataLine);
	    my $electronDensity = 10.0**$dataColumns[4];
	    $electronDensities[++$#electronDensities] = $electronDensity;
	    my $hiDensity       = 10.0**$dataColumns[6];
	    $hiDensities      [++$#hiDensities      ] = $hiDensity;
	    my $hiiDensity      = 10.0**$dataColumns[7];
	    $hiiDensities     [++$#hiiDensities     ] = $hiiDensity;
	}
	close(overviewHandle);
	unlink($overviewTempFile);

	# Store cooling function data.
	@{${${$coolingFunctions{'coolingFunction'}}[$iCoolingFunction]}{'coolingRate'}->{'datum'}}     = @coolingRates;
	@{${${$coolingFunctions{'coolingFunction'}}[$iCoolingFunction]}{'temperature'}->{'datum'}}     = @temperatures;

	# Store chemical state data.
	@{${${$chemicalStates{'ionizationState'}}[$iChemicalState]}{'electronDensity'}->{'datum'}} = @electronDensities;
	@{${${$chemicalStates{'ionizationState'}}[$iChemicalState]}{'hiDensity'      }->{'datum'}} = @hiDensities;
	@{${${$chemicalStates{'ionizationState'}}[$iChemicalState]}{'hiiDensity'     }->{'datum'}} = @hiiDensities;
	@{${${$chemicalStates{'ionizationState'}}[$iChemicalState]}{'temperature'    }->{'datum'}} = @temperatures;
	
    }
    
    # Cooling functions:
    # Specify extrapolation methods in temperature.
    ${${${$coolingFunctions{'extrapolation'}}{'temperature'}}[0]}{'limit'}  = "low";
    ${${${$coolingFunctions{'extrapolation'}}{'temperature'}}[0]}{'method'} = "power law";
    ${${${$coolingFunctions{'extrapolation'}}{'temperature'}}[1]}{'limit'}  = "high";
    ${${${$coolingFunctions{'extrapolation'}}{'temperature'}}[1]}{'method'} = "power law";
    # Specify extrapolation methods in metallicity.
    ${${${$coolingFunctions{'extrapolation'}}{'metallicity'}}[0]}{'limit'}  = "low";
    ${${${$coolingFunctions{'extrapolation'}}{'metallicity'}}[0]}{'method'} = "fixed";
    ${${${$coolingFunctions{'extrapolation'}}{'metallicity'}}[1]}{'limit'}  = "high";
    ${${${$coolingFunctions{'extrapolation'}}{'metallicity'}}[1]}{'method'} = "fixed";
    # Add some description.
    $coolingFunctions{'description'} = "CIE cooling functions computed by Cloudy ".$cloudyVersion;
    ${$coolingFunctions{'units'}}[0] = "Temperature: Kelvin";
    ${$coolingFunctions{'units'}}[1] = "Cooling rate: Lambda(T)/ergs cm^3 s^-1";
  
    # Chemical states:
    # Specify extrapolation methods in temperature.
    ${${${$chemicalStates{'extrapolation'}}{'temperature'}}[0]}{'limit'}  = "low";
    ${${${$chemicalStates{'extrapolation'}}{'temperature'}}[0]}{'method'} = "fixed";
    ${${${$chemicalStates{'extrapolation'}}{'temperature'}}[1]}{'limit'}  = "high";
    ${${${$chemicalStates{'extrapolation'}}{'temperature'}}[1]}{'method'} = "fixed";
    # Specify extrapolation methods in metallicity.
    ${${${$chemicalStates{'extrapolation'}}{'metallicity'}}[0]}{'limit'}  = "low";
    ${${${$chemicalStates{'extrapolation'}}{'metallicity'}}[0]}{'method'} = "fixed";
    ${${${$chemicalStates{'extrapolation'}}{'metallicity'}}[1]}{'limit'}  = "high";
    ${${${$chemicalStates{'extrapolation'}}{'metallicity'}}[1]}{'method'} = "fixed";
    # Add some description.
    $chemicalStates{'description'} = "CIE ionization states computed by Cloudy 08.00";
    ${$chemicalStates{'units'}}[0] = "Temperature: Kelvin";
    ${$chemicalStates{'units'}}[1] = "Densities: by number relative to hydrogen";

    # Output cooling functions to an XML file.
    if ( $computeCoolingFunctions == 1 ) {
	my $coolingFunctions = \%coolingFunctions;
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"coolingFunctions");
	print coolingFunctionOutHndl $xmlOutput->XMLout($coolingFunctions);
	close(coolingFunctionOutHndl);
    }
    
    # Output chemical states to an XML file.
    if ( $computeChemicalStates == 1 ) {
	my $chemicalStateStructure = \%chemicalStates;
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"chemicalStates");
	print chemicalStateOutHndl $xmlOutput->XMLout($chemicalStateStructure);
	close(chemicalStateOutHndl);
    }

    # Write message.
    print "...done\n";
    
}
    
exit;
