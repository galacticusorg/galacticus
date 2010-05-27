#!/usr/bin/env perl
use XML::Simple;
use Data::Dumper;
use Fcntl qw (:flock);

# Driver script for Cloudy.
# Andrew Benson (26-Jan-2010)

# Get arguments.
if ( $#ARGV != 1 ) {die "Usage: Cloudy_Driver.pl <logMetallicityMaximum> <coolingFunctionFile>"};
$logMetallicityMaximum = $ARGV[0];
$coolingFunctionFile   = $ARGV[1];

# Check if the file exists.
unless ( -e $coolingFunctionFile ) {

    # Open and lock the cooling function file.
    open(outHndl,">".$coolingFunctionFile);
    flock(outHndl,LOCK_EX);

    # (Logarithmic) temperature range.
    $logTemperatureMinimum = 2.500;
    $logTemperatureMaximum = 9.000;
    $logTemperatureDelta   = 0.025;
    
    # Logarithmic metallicity range (a zero metallicity case is always included too).
    $logMetallicityMinimum = -4.00;
    $logMetallicityDelta   =  0.25;
    
    # Generate metallicities array.
    @logMetallicities = ( -999.0 ); # Proxy for zero metallicity.
    for($logMetallicity=$logMetallicityMinimum;$logMetallicity<=$logMetallicityMaximum;$logMetallicity+=$logMetallicityDelta) {
	$logMetallicities[++$#logMetallicities] = $logMetallicity;
    }
    
    # Specify Solar and primodial helium abundances (as used in Cloudy).
    $heliumAbundancePrimordial = 0.072;
    $heliumAbundanceSolar      = 0.100;
    
    # Download the code.
    unless ( -e "aux/c08.00.tar.gz" ) {
	print "Cloudy_Driver.pl: downloading Cloudy code.\n";
	system("wget ftp://gradj.pa.uky.edu/gary/cloudy_gold/c08.00.tar.gz -O aux/c08.00.tar.gz");
	die("Cloudy_Driver.pl: FATAL - failed to download Cloudy code.") unless ( -e "aux/c08.00.tar.gz" );
    }
    
    # Unpack the code.
    unless ( -e "aux/c08.00" ) {
	print "Cloudy_Driver.pl: unpacking Cloudy code.\n";
	system("tar -x -v -z -C aux -f aux/c08.00.tar.gz");
	die("Cloudy_Driver.pl: FATAL - failed to unpack Cloudy code.") unless ( -e "aux/c08.00" );
    }
    
    # Build the code.
    unless ( -e "aux/c08.00/source/cloudy.exe" ) {
	print "Cloudy_Driver.pl: compiling Cloudy code.\n";
	system("cd aux/c08.00/source; chmod u=wrx precompile.pl; make");
	die("Cloudy_Driver.pl: FATAL - failed to build Cloudy code.") unless ( -e "aux/c08.00/source/cloudy.exe" );
    }
   
    # Temporary file for Cloudy data.
    $cloudyTempFile = "./cloudy_cooling.tmp";
    
    # Counter for number of cooling functions tabulated.
    $iCoolingFunction = -1;
    
    # Loop over metallicities.
    foreach $logMetallicity ( @logMetallicities ) {
	
	# Increment cooling function counter.
	++$iCoolingFunction;
	
	# Destroy the previous cooling function data.
	undef(@data);
	
	# Store the metallicity for this cooling function.
	${${$coolingFunctions{'coolingFunction'}}[$iCoolingFunction]}{'metallicity'} = $logMetallicity;
	
	# Loop over temperatures.
	for($logTemperature=$logTemperatureMinimum;$logTemperature<=$logTemperatureMaximum;$logTemperature+=$logTemperatureDelta) {
	    $temperature = 10.0**$logTemperature;

	    # Run Cloudy.
	    open(cloudyPipe,"|aux/c08.00/source/cloudy.exe");
	    print cloudyPipe "print off\n";
	    print cloudyPipe "background 0\n";
	    print cloudyPipe "cosmic ray background\n";
	    print cloudyPipe "stop zone 1\n";
	    print cloudyPipe "hden 0.0\n";
	    if ( $logMetallicity <= -999.0 ) {
		print cloudyPipe "abundances primordial\n";
	    } else {
		print cloudyPipe "metals _log ".$logMetallicity."\n";
		# Assume a linear growth of helium abundance with metallicity.
		$heliumAbundance = $heliumAbundancePrimordial+($heliumAbundanceSolar-$heliumAbundancePrimordial)*(10.0**$logMetallicity);
		print cloudyPipe "element abundance linear helium ".$heliumAbundance."\n";
	    }
	    print cloudyPipe "constant temper ".$logTemperature."\n";
	    print cloudyPipe "no molecules\n";
	    print cloudyPipe "set trim -20\n";
	    print cloudyPipe "punch cooling \"".$cloudyTempFile."\"\n";
	    close(cloudyPipe);
	    
	    # Extract the cooling rate.
	    open(coolHandle,$cloudyTempFile);
	    $headerLine = <coolHandle>;
	    $dataLine = <coolHandle>;
	    close(coolHandle);
	    unlink($cloudyTempFile);
	    @dataColumns = split(/\s+/,$dataLine);
	    $coolingRate = $dataColumns[3];
	    
	    $data[++$#data] = $temperature." ".$coolingRate;
	}
	
	@{${${$coolingFunctions{'coolingFunction'}}[$iCoolingFunction]}{'datum'}} = @data;
	
    }
    
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
    ${$coolingFunctions{'description'}}[0] = "CIE cooling functions computed by Cloudy 08.00";
    ${$coolingFunctions{'description'}}[1] = "Data are: log10(T/K) Lambda(T)/ergs cm^3 s^-1";
    
    # Output the data to an XML file.
    $coolingFunctions = \%coolingFunctions;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"coolingFunctions");
    print outHndl $xmlOutput->XMLout($coolingFunctions);
    close(outHndl);
}
    
exit;
