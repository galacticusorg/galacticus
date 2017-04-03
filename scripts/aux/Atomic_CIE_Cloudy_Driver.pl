#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
use Fcntl qw (:flock);
use Cloudy;

# Driver script for Cloudy.
# Andrew Benson (26-Jan-2010)

# Get arguments.
die "Usage: Atomic_CIE_Cloudy_Driver.pl <logMetallicityMaximum> <coolingFunctionFile> <chemicalStateFile> <fileFormatVersion>"
    unless ( scalar(@ARGV) == 4 );
my $logMetallicityMaximum = $ARGV[0];
my $coolingFunctionFile   = $ARGV[1];
my $chemicalStateFile     = $ARGV[2];
my $fileFormat            = $ARGV[3];

# Ensure the requested file format version is compatible.
my $fileFormatCurrent = 1;
die('Atomic_CIE_Cloudy_Driver.pl: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
    unless ( $fileFormat == $fileFormatCurrent );

# Determine if we need to compute cooling functions.
my $computeCoolingFunctions = 0;
if ( -e $coolingFunctionFile ) {
    my $coolingFunction = new PDL::IO::HDF5($coolingFunctionFile);
    if ( grep {$_ eq "fileFormat"} $coolingFunction->attrs() ) {
	(my $fileFormatFile) = $coolingFunction->attrGet('fileFormat');
	$computeCoolingFunctions = 1
	    unless ( $fileFormatFile == $fileFormatCurrent );
    }
} else {
    $computeCoolingFunctions = 1;
}

# Determine if we need to compute chemical states.
my $computeChemicalStates = 0;
if ( -e $chemicalStateFile ) {
    my $chemicalState = new PDL::IO::HDF5($chemicalStateFile);
    if ( grep {$_ eq "fileFormat"} $chemicalState->attrs() ) {
	(my $fileFormatFile) = $chemicalState->attrGet('fileFormat');
	$computeChemicalStates = 1
	    unless ( $fileFormatFile == $fileFormatCurrent );
    }
} else {
    $computeChemicalStates = 1;
}

# Check if we need to do calculations.
my $coolingFunctionLock;
my $chemicalStateLock  ;
if ( $computeCoolingFunctions == 1 || $computeChemicalStates == 1 ) {

    # Open and lock the cooling function and chemical state files.
    if ( $computeCoolingFunctions == 1 ) {
	open($coolingFunctionLock,">".$coolingFunctionFile.".lock");
	flock($coolingFunctionLock,LOCK_EX);
    }
    if ( $computeChemicalStates == 1 ) {
	open($chemicalStateLock,">".$chemicalStateFile.".lock");
	flock($chemicalStateLock,LOCK_EX);
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
	push(@logMetallicities,$logMetallicity);
    }

    # Specify Solar and primodial helium abundances (as used in Cloudy).
    my $heliumAbundancePrimordial = 0.072;
    my $heliumAbundanceSolar      = 0.100;

    # Ensure Cloudy is built.
    (my $cloudyPath, my $cloudyVersion) = &Cloudy::Initialize();
    
    # Temporary files for Cloudy data.
    my $coolingTempFile    = "cloudy_cooling.tmp";
    my $overviewTempFile   = "cloudy_overview.tmp";
    
    # Counter for number of metallicities tabulated.
    my $iMetallicity = -1;

    # Write message.
    print "Computing cooling functions and chemical states using Cloudy (this may take a long time)...\n";
    
    # Initialize data structure1.
    my $cloudyData;

    # Store metallicities.
    $cloudyData->{'metallicity'} = pdl @logMetallicities;

    # Initialize table names.
    my @coolingFunctionTables = ( "coolingRate" );
    my @chemicalStateTables   = ( "electronDensity", "hiDensity", "hiiDensity" );
    
    # Loop over metallicities.
    foreach my $logMetallicity ( @logMetallicities ) {
	
	# Increment metallicity counter.
	++$iMetallicity;
	print " ⮡ Computing for log(Z/Z☉)=".$logMetallicity."\n";
	
	# Destroy the previous cooling function data.
	my @temperatures     ;
	my @coolingRates     ;
	my @electronDensities;
	my @hiDensities      ;
	my @hiiDensities     ;
	
	# Run Cloudy.
	open(cloudyScript,">".$cloudyPath."/source/input.in");
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
	system("cd ".$cloudyPath."/source; cloudy.exe -r input");

	# Extract the cooling rate.
	open(coolHandle,$cloudyPath."/source/".$coolingTempFile);
	my $headerLine = <coolHandle>;
	while ( my $dataLine = <coolHandle> ) {
	    my $separator   = <coolHandle>;
	    my @dataColumns = split(/\s+/,$dataLine);
	    my $temperature = $dataColumns[1];
	    my $heatingRate = $dataColumns[2];
	    my $coolingRate = $dataColumns[3];
	    push(@temperatures,$temperature);
	    push(@coolingRates,$coolingRate);
	}
	close(coolHandle);
	unlink($coolingTempFile);

	# Extract the electron and hydrogen density.
	open(overviewHandle,$cloudyPath."/source/".$overviewTempFile);
	$headerLine = <overviewHandle>;
	while ( my $dataLine = <overviewHandle> ) {
	    my $separator   = <overviewHandle>;
	    my @dataColumns = split(/\s+/,$dataLine);
	    my $electronDensity = 10.0**$dataColumns[4];
	    push(@electronDensities,$electronDensity);
	    my $hiDensity       = 10.0**$dataColumns[6];
	    push(@hiDensities      ,      $hiDensity);
	    my $hiiDensity      = 10.0**$dataColumns[7];
	    push(@hiiDensities     ,     $hiiDensity);
	}
	close(overviewHandle);
	unlink($overviewTempFile);

	# Initialize data.
	foreach ( @coolingFunctionTables, @chemicalStateTables ) {
	    $cloudyData->{$_} = pdl zeroes(scalar(@logMetallicities),scalar(@temperatures))
		unless ( defined($cloudyData->{$_}) );
	}

	# Store temperatures.
	$cloudyData->{'temperature'} = pdl @temperatures
	    unless ( defined($cloudyData->{'temperature'}) );
	
	# Store cooling function data.
	$cloudyData->{'coolingRate'    }->(($iMetallicity),:) .= pdl @coolingRates     ;

	# Store chemical state data.
	$cloudyData->{'electronDensity'}->(($iMetallicity),:) .= pdl @electronDensities;
	$cloudyData->{'hiDensity'      }->(($iMetallicity),:) .= pdl @hiDensities      ;
	$cloudyData->{'hiiDensity'     }->(($iMetallicity),:) .= pdl @hiiDensities     ;
    }
    
    # Output cooling functions to an HDF5 file.
    if ( $computeCoolingFunctions == 1 ) {
	my $coolingFunction = new PDL::IO::HDF5(">".$coolingFunctionFile);
	# Store datasets.
	$coolingFunction->dataset($_)->set($cloudyData->{$_})
	    foreach ( 'metallicity', 'temperature', @coolingFunctionTables );
	# Add extrapolation attributes.
	my $temperature = $coolingFunction->dataset('temperature');
	my $metallicity = $coolingFunction->dataset('metallicity');
	$temperature->attrSet(extrapolateLow  => "powerLaw");
	$temperature->attrSet(extrapolateHigh => "powerLaw");
	$metallicity->attrSet(extrapolateLow  => "fix"     );
	$metallicity->attrSet(extrapolateHigh => "fix"     );
	# Add units attributes.
	$temperature->attrSet(units     => "K"    );
	$temperature->attrSet(unitsInSI => pdl 1.0);
	# Add provenance.
	$coolingFunction->attrSet(description => "CIE cooling functions computed by Cloudy ".$cloudyVersion);
	# Add file format.
	$coolingFunction->attrSet(fileFormat => pdl long($fileFormatCurrent));
	# Destroy the lock file.
	close($coolingFunctionLock);
    }
    
    # Output chemical states to an HDF5 file.
    if ( $computeChemicalStates == 1 ) {
	my $chemicalState = new PDL::IO::HDF5(">".$chemicalStateFile);
	$chemicalState->dataset($_)->set($cloudyData->{$_})
	    foreach ( 'metallicity', 'temperature', @chemicalStateTables );
	# Add extrapolation attributes.
	my $temperature = $chemicalState->dataset('temperature');
	my $metallicity = $chemicalState->dataset('metallicity');
	$temperature->attrSet(extrapolateLow  => "fix");
	$temperature->attrSet(extrapolateHigh => "fix");
	$metallicity->attrSet(extrapolateLow  => "fix");
	$metallicity->attrSet(extrapolateHigh => "fix");
	# Add units attributes.
	$temperature->attrSet(units     => "K"    );
	$temperature->attrSet(unitsInSI => pdl 1.0);
	# Add provenance.
	$chemicalState->attrSet(description => "CIE ionization states computed by Cloudy ".$cloudyVersion);
	# Add file format.
	$chemicalState->attrSet(fileFormat => pdl long($fileFormatCurrent));
	close($chemicalStateLock);
    }

    # Write message.
    print "...done\n";
    
}
    
exit;
