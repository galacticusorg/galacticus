#!/usr/bin/env perl
use XML::Simple;
use Data::Dumper;
use PDL;
use File::Copy;

# Compute Vega-AB offsets and effective wavelengths for any filters which do not already have them.
# Andrew Benson (16-Feb-2010)

# Check arguments make sense.
if ( $#ARGV > 1 ) {die("Usage: vega_offset_effective_lambda.pl [<filtersDirectory> [<vegaSpectrumFile>]]")};
# Get the filters directory if one is specified.
if ( $#ARGV >= 0 ) {
    $filtersDirectory = $ARGV[0];
} else {
    $filtersDirectory = "./data/filters";
}
# Get the Vega spectrum file.
if ( $#ARGV == 1 ) {
    $vegaSpectrumFile = $ARGV[1];
} else {
    # None given, so use the default.
    $vegaSpectrumFile = "./data/vega/A0V_Castelli.xml";
    # Check that the file exists - if not, attempt to download data and create it.
    unless ( -e $vegaSpectrumFile ) {
	print "Cannot find A0V_Castelli.xml file - will attempt to download data and create it....\n";
	$castelliURL = "http://wwwuser.oat.ts.astro.it/castelli/grids/gridp00k2odfnew/fp00t9500g40k2odfnew.dat";
	open(pipeHndl,"wget ".$castelliURL." -O - |");
	while ( $line = <pipeHndl> ) {
	    if ( $line =~ m/^\s*FLUX/ ) {
		$line =~ s/^\s*//;
		$line =~ s/\s*$//;
		@columns = split(/\s+/,$line);
		if ( $#columns > 0) {
		    $columns[2] = $columns[2]*10.0;             # Convert to Angstroms.
		    $columns[4] = $columns[4]/($columns[2]**2); # Convert to F_lambda.
		    ${$vegaSpectrum{'datum'}}[++$#{$vegaSpectrum{'datum'}}] = $columns[2]."\t".$columns[4];
		}
	    }
	}
	close(pipeHndl);
	$vegaSpectrum{'description'} = "Model spectrum of A0V stars from Castelli & Kurucz (2004).";
	$vegaSpectrum{'origin'} = $castelliURL;
	${$vegaSpectrum{'units'}}[0] = "wavelengths: Angstroms";
	${$vegaSpectrum{'units'}}[1] = "fluxes, F_lambda: arbitrary";
	$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"spectrum");
	system("mkdir -p `dirname $vegaSpectrumFile`");
	open(outHndl,">".$vegaSpectrumFile);
	print outHndl $xmlOutput->XMLout(\%vegaSpectrum);
	close(outHndl);
    }
}

# Create XML object to process XML files.
$xml = new XML::Simple;

# Read the spectrum of Vega.
$vegaSpectrum = $xml->XMLin($vegaSpectrumFile);
# Extract the SED to PDL variables.
$first = 1;
foreach $datum ( @{$vegaSpectrum->{'datum'}} ) {
    $datum =~ s/^\s*//;
    $datum =~ s/\s*$//;
    @columns = split(/\s+/,$datum);
    if ( $first == 1 ) {
	$wavelengthsVega = pdl $columns[0];
	$spectrumVega    = pdl $columns[1];
    } else {
	$wavelengthsVega = append($wavelengthsVega,$columns[0]);
	$spectrumVega    = append($spectrumVega   ,$columns[1]);
    }
    $first = 0;
}

# Include the V-band filter used for defining the Vega magnitude system as the first on the list.
$fileNames[0] = "Buser_V.xml";

# Open the filters directory and search for filter files.
opendir(dirHndl,$filtersDirectory);
while ( $fileName = readdir(dirHndl) ) {
    # Select XML files.
    if ( $fileName =~ m/\.xml$/ && $fileName !~ m/^Buser_V\.xml$/ ) {	
	$fileNames[++$#fileNames] = $fileName;
    }
}
closedir(dirHndl);
# Process all filters.
foreach $fileName ( @fileNames ) {
    print "Processing file: $fileName\n";
    # Process the XML data.
    $filePath = $filtersDirectory."/".$fileName;
    $filter = $xml->XMLin($filePath);
    # Extract the filter curve to PDL variables.
    $first = 1;
    foreach $datum ( @{$filter->{'response'}->{'datum'}} ) {
	$datum =~ s/^\s*//;
	$datum =~ s/\s*$//;
	@columns = split(/\s+/,$datum);
	if ( $first == 1 ) {
	    $wavelengths = pdl $columns[0];
	    $response    = pdl $columns[1];
	} else {
	    $wavelengths = append($wavelengths,$columns[0]);
	    $response    = append($response   ,$columns[1]);
	}
	$first = 0;
    }
    $inserts = "";

    # Check if an effective wavelength is listed - if not, compute it.
    unless ( exists($filter->{'effectiveWavelength'}) ) {
	$effectiveWavelength = sum($response*$wavelengths)/sum($response);
	$inserts .= "  <effectiveWavelength>".$effectiveWavelength."</effectiveWavelength>\n";
	print "  -> Computed effective wavelength: ".$effectiveWavelength." Angstroms\n";
    }
    # Check if an Vega-AB offset is listed - if not, compute it.
    unless ( exists($filter->{'vegaOffset'}) && $fileName !~ m/^Buser_V\.xml$/ ) {
	# Get range of wavelengths in the filter.
	$wavelengthMinimum = pdl $wavelengths->index(0);
	$wavelengthMaximum = pdl $wavelengths->index(nelem($wavelengths)-1);
	$wavelengthsVegaInFilter = where($wavelengthsVega,$wavelengthsVega >= $wavelengthMinimum & $wavelengthsVega <= $wavelengthMaximum );
	# Make and sort joint wavelength array.
	$wavelengthsJoint = append($wavelengths,$wavelengthsVegaInFilter);
	$wavelengthsJoint = uniq $wavelengthsJoint;
	# Interpolate filter and Vega spectrum to these wavelengths.
	($responseJoint,$responseError) = interpolate($wavelengthsJoint, $wavelengths    , $response    );
	($spectrumJoint,$spectrumError) = interpolate($wavelengthsJoint, $wavelengthsVega, $spectrumVega);
	$spectrumAB = 1.0/$wavelengthsJoint**2;
	# Get the filtered spectrum.
	$filteredSpectrum   = $responseJoint*$spectrumJoint;
	$filteredSpectrumAB = $responseJoint*$spectrumAB;
	# Compute the integrated flux.
	$fluxVega = pdl 0.0;
	$fluxAB   = pdl 0.0;
	for($i=0;$i<nelem($wavelengthsJoint);++$i) {
	    if ( $i == 0 ) {
		$deltaWavelength = $wavelengthsJoint->index($i+1)-$wavelengthsJoint->index($i  );
	    } elsif ( $i == nelem($wavelengthsJoint)-1 ) {
		$deltaWavelength = $wavelengthsJoint->index($i  )-$wavelengthsJoint->index($i-1);
	    } else {
		$deltaWavelength = $wavelengthsJoint->index($i+1)-$wavelengthsJoint->index($i-1);
	    }
	    $fluxVega += ($filteredSpectrum->index($i))*$deltaWavelength;
	    $fluxAB   += ($filteredSpectrumAB->index($i))*$deltaWavelength;
	}
	if ( $fileName =~ m/^Buser_V\.xml$/ ) {
	    # Keep a copy of the fluxes in the V-band for later use.
	    $fluxVegaV = $fluxVega;
	    $fluxABV   = $fluxAB;
	    $vegaOffset = 0.0; # By definition, no offset in V-band.
	} else {
	    $vegaOffset = 2.5*log10($fluxVega*$fluxABV/$fluxVegaV/$fluxAB);
	}
	unless ( exists($filter->{'vegaOffset'}) ) {
	    $inserts .= "  <vegaOffset>".$vegaOffset."</vegaOffset>\n";
	    print "  -> Computed Vega-AB offset: ".$vegaOffset."\n";
	}
    }

    # Insert new data into file if necessary.
    unless ( $inserts eq "" ) {
	open(tmpHndl,">filter.tmp");
	open(inHndl,$filePath);
	while ( $line = <inHndl> ) {
	    $line =~ s/<\/filter>/$inserts<\/filter>/;
	    print tmpHndl $line;
	}
	close(inHndl);
	close(tmpHndl);
	unlink($filePath);
	move("filter.tmp",$filePath);
    }
}

exit;
