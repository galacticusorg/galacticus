#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use Data::Dumper;
use PDL;
use File::Copy;

# Compute Vega-AB offsets and effective wavelengths for any filters which do not already have them.
# Andrew Benson (16-Feb-2010)

# Check arguments make sense.
die("Usage: vega_offset_effective_lambda.pl [<filtersDirectory> [<vegaSpectrumFile>]]")
    unless ( scalar(@ARGV) <= 3 );
# Get the filters directory if one is specified.
my $filtersDirectory;
if ( scalar(@ARGV) > 0 ) {
    $filtersDirectory = $ARGV[0];
} else {
    $filtersDirectory = "./data/filters";
}
# Get the Vega spectrum file.
my $vegaSpectrumFile;
if ( scalar(@ARGV) == 2 ) {
    $vegaSpectrumFile = $ARGV[1];
} else {
    # None given, so use the default.
    $vegaSpectrumFile = "./data/stellarAstrophysics/vega/A0V_Castelli.xml";
    # Check that the file exists - if not, attempt to download data and create it.
    unless ( -e $vegaSpectrumFile ) {
	print "Cannot find A0V_Castelli.xml file - will attempt to download data and create it....\n";
	my %vegaSpectrum;
	my $castelliURL = "http://wwwuser.oats.inaf.it/castelli/grids/gridp00k2odfnew/fp00t9500g40k2odfnew.dat";
	open(pipeHndl,"wget ".$castelliURL." -O - |");
	while ( my $line = <pipeHndl> ) {
	    if ( $line =~ m/^\s*FLUX/ ) {
		$line =~ s/^\s*//;
		$line =~ s/\s*$//;
		my @columns = split(/\s+/,$line);
		if ( scalar(@columns) > 1 ) {
		    $columns[2] = $columns[2]*10.0;             # Convert to Angstroms.
		    $columns[4] = $columns[4]/($columns[2]**2); # Convert to F_lambda.
		    push(@{$vegaSpectrum{'datum'}},$columns[2]."\t".$columns[4]);
		}
	    }
	}
	close(pipeHndl);
	$vegaSpectrum{'description'} = "Model spectrum of A0V stars from Castelli & Kurucz (2004).";
	$vegaSpectrum{'origin'} = $castelliURL;
	${$vegaSpectrum{'units'}}[0] = "wavelengths: Angstroms";
	${$vegaSpectrum{'units'}}[1] = "fluxes, F_lambda: arbitrary";
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"spectrum");
	system("mkdir -p `dirname $vegaSpectrumFile`");
	open(outHndl,">".$vegaSpectrumFile);
	print outHndl $xmlOutput->XMLout(\%vegaSpectrum);
	close(outHndl);
    }
}

# Create XML object to process XML files.
my $xml = new XML::Simple;

# Read the spectrum of Vega.
my $vegaSpectrum = $xml->XMLin($vegaSpectrumFile);
# Extract the SED to PDL variables.
my $wavelengthsVega = pdl [];
my $spectrumVega    = pdl [];
foreach my $datum ( @{$vegaSpectrum->{'datum'}} ) {
    $datum =~ s/^\s*//;
    $datum =~ s/\s*$//;
    my @columns = split(/\s+/,$datum);
    $wavelengthsVega = append($wavelengthsVega,$columns[0]);
    $spectrumVega    = append($spectrumVega   ,$columns[1]);
}

# Include the V-band filter used for defining the Vega magnitude system as the first on the list.
my @fileNames = ( "Buser_V.xml" );

# Open the filters directory and search for filter files.
opendir(my $dirHndl,$filtersDirectory);
while ( my $fileName = readdir($dirHndl) ) {
    # Select XML files.
    push(@fileNames,$fileName)
	if ( $fileName =~ m/\.xml$/ && $fileName !~ m/^Buser_V\.xml$/ );
}
closedir($dirHndl);
# Process all filters.
my $fluxVegaV;
my $fluxABV;
foreach my $fileName ( @fileNames ) {
    print "Processing file: $fileName\n";
    # Process the XML data.
    my $filePath = $filtersDirectory."/".$fileName;
    my $filter = $xml->XMLin($filePath);
    # Extract the filter curve to PDL variables.
    my $wavelengths = pdl [];
    my $response    = pdl [];
    foreach my $datum ( @{$filter->{'response'}->{'datum'}} ) {
	$datum =~ s/^\s*//;
	$datum =~ s/\s*$//;
	my @columns = split(/\s+/,$datum);
	$wavelengths = append($wavelengths,$columns[0]);
	$response    = append($response   ,$columns[1]);
    }
    my $inserts = "";

    # Check if an effective wavelength is listed - if not, compute it.
    unless ( exists($filter->{'effectiveWavelength'}) ) {
	my $effectiveWavelength = sum($response*$wavelengths)/sum($response);
	$inserts .= "  <effectiveWavelength>".$effectiveWavelength."</effectiveWavelength>\n";
	print "  -> Computed effective wavelength: ".$effectiveWavelength." Angstroms\n";
    }
    # Check if an Vega-AB offset is listed - if not, compute it.
    unless ( exists($filter->{'vegaOffset'}) && $fileName !~ m/^Buser_V\.xml$/ ) {
	# Get range of wavelengths in the filter.
	my $wavelengthMinimum = pdl $wavelengths->index(0);
	my $wavelengthMaximum = pdl $wavelengths->index(nelem($wavelengths)-1);
	my $wavelengthsVegaInFilter =
	    where(
		$wavelengthsVega,
		($wavelengthsVega >= $wavelengthMinimum)
		& 
		($wavelengthsVega <= $wavelengthMaximum)
	    );
	# Make and sort joint wavelength array.
	my $wavelengthsJoint = append($wavelengths,$wavelengthsVegaInFilter);
	$wavelengthsJoint = uniq $wavelengthsJoint;
	# Interpolate filter and Vega spectrum to these wavelengths.
	(my $responseJoint, my $responseError) = interpolate($wavelengthsJoint, $wavelengths    , $response    );
	(my $spectrumJoint, my $spectrumError) = interpolate($wavelengthsJoint, $wavelengthsVega, $spectrumVega);
	my $spectrumAB = 1.0/$wavelengthsJoint**2;
	# Get the filtered spectrum.
	my $filteredSpectrum   = $responseJoint*$spectrumJoint;
	my $filteredSpectrumAB = $responseJoint*$spectrumAB;
	# Compute the integrated flux.
	my $fluxVega = pdl 0.0;
	my $fluxAB   = pdl 0.0;
	for(my $i=0;$i<nelem($wavelengthsJoint);++$i) {
	    my $deltaWavelength;
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
	my $vegaOffset;
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
	open(my $tmpHndl,">filter.tmp");
	open(my $inHndl,$filePath);
	while ( my $line = <$inHndl> ) {
	    $line =~ s/<\/filter>/$inserts<\/filter>/;
	    print $tmpHndl $line;
	}
	close($inHndl);
	close($tmpHndl);
	unlink($filePath);
	move("filter.tmp",$filePath);
    }
}

exit;
