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
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::NiceSlice;
use Fcntl qw(SEEK_SET);

# Convert Bruzual & Charlot stellar population data to Galacticus' format.
# Andrew Benson (19-October-2010).

# Specify the data directory.
my $dataDirectory = "SSP_BC2003";

# Specify models to convert.
my @modelsToConvert = ( "padova_1994_chabrier_imf", "padova_1994_salpeter_imf" );

# Ensure data files have been downloaded.
foreach my $model ( @modelsToConvert ) {
    die("Please download the model archives from the BC2003 website before running this script.") unless ( -e $dataDirectory."/bc03.models.".$model.".tar.gz" );
}

# Lookup table for metallicity labels.
my %metallicityLookup = (
    22 => 0.0001,
    32 => 0.0004,
    42 => 0.0040,
    52 => 0.0080,
    62 => 0.0200,
    72 => 0.0500
    );

# Lookup table for resolutions.
my %resolutionLookup = (
    lr => "lowResolution",
    hr => "highResolution"
    );

# Loop over models.
foreach my $model ( @modelsToConvert ) {

    # Extract the IMF name and other data from the model name.
    my $tracks;
    my $year;
    my $IMF;
    if ( $model =~ m/([a-z]+)_(\d+)_([a-z]+)_imf/ ) {
	$tracks = ucfirst($1);
	$year   = $2;
	$IMF    = $3;
    } else {
	die("Model name is in unrecognized format.");
    }

    # Unpack the data.
    system("cd ".$dataDirectory."; tar xvfz bc03.models.".$model.".tar.gz" );

    # Build the base directory for this model.
    my $baseDirectory = $dataDirectory."/bc03/models/".$tracks.$year."/".$IMF;

    # Unpack the ASCII data files.
    system("find ".$baseDirectory." -name \"*_ASCII.gz\" | xargs gunzip");

    # Loop over low and high-res spectra.
    foreach my $resolution ( "lr", "hr" ) {

	# Build a list of files to process.
	my %files;
	opendir(my $inDir,$baseDirectory);
	while ( my $file = readdir($inDir) ) {
	    if ( $file =~ m/bc2003_$resolution\_m(\d+)_[a-z]+_ssp.ised_ASCII/ ) {	    
		my $metallicity = $1;
		$files{$metallicity} = $file;
	    }
	}
	closedir($inDir);

	# Reset datasets.
	my $metallicities     = pdl [];
	my $ageDataset        = pdl [];
	my $wavelengthDataset = pdl [];
	my $spectra;

	# Process files.
	my $iMetallicity = -1;
	foreach my $file ( sort(keys(%files)) ) {
	    ++$iMetallicity;	

	    # Reset datasets.
	    my $ages        = pdl [];
	    my $wavelengths = pdl [];

	    # Get the file name.
	    my $fileName = $baseDirectory."/".$files{$file};
	    print "Processing file: ".$fileName."\n";
	    # Extract the metallicity.
	    $metallicities = $metallicities->append($metallicityLookup{$file});

	    # Open and process the file.
	    open(my $iHndl,$fileName);
	    my $iLine    =  0;
	    my $ageCount = -1;
	    while ( $ageCount == -1 || nelem($ages) < $ageCount ) {
		++$iLine;
		my $line = <$iHndl>;
		$line =~ s/^\s*//;
		$line =~ s/\s*$//;
		my @columns = split(/\s+/,$line);
		$ageCount = shift(@columns) if ( $ageCount == -1 );
		$ages = $ages->append(pdl @columns);
	    }
	    # Store ages if not already done. (Convert to Gyr as we do so.)
	    my $ageDataset = $ages/1.0e9 unless ( nelem($ageDataset) > 0 );

	    # Skip over lines prior to those beginning with "Padova".
	    my $gotPadova = 0;
	    my $filePos;
	    while ( my $line = <$iHndl> ) {
		++$iLine;
		if ( $line =~ m/^Padova/ ) {
		    $gotPadova = 1;
		} else {
		    if ( $gotPadova == 1 ) {
			# This line does not start with "Padova", but previous lines have. We've skipped enough.
			--$iLine;
			seek($iHndl,$filePos,SEEK_SET);
			last;
		    }
		}
		$filePos = tell($iHndl);
	    }

	    # Get wavelengths from the file.
	    my $wavelengthCount = -1;
	    while ( $wavelengthCount == -1 || nelem($wavelengths) < $wavelengthCount ) {
		++$iLine;
		my $line = <$iHndl>;
		$line =~ s/^\s*//;
		$line =~ s/\s*$//;
		my @columns = split(/\s+/,$line);
		$wavelengthCount = shift(@columns) if ( $wavelengthCount == -1 );
		$wavelengths = $wavelengths->append(pdl @columns);
	    }
	    # Store wavelengths if not already done.
	    my $wavelengthDataset = $wavelengths unless ( nelem($wavelengthDataset) > 0 );

	    # Create a PDL to hold the spectra.
	    my $metallicityCount = keys %files;
	    $spectra = pdl zeroes($metallicityCount,$ageCount,$wavelengthCount) unless (defined($spectra));

	    # Loop over ages.
	    for (my $iAge=0;$iAge<$ageCount;++$iAge) {
		# Grab spectrum for all wavelengths.
		my $spectrum    = pdl [];
		$wavelengthCount = -1;
		while ( $wavelengthCount == -1 || nelem($spectrum) < $wavelengthCount ) {
		    ++$iLine;
		    my $line = <$iHndl>;
		    $line =~ s/^\s*//;
		    $line =~ s/\s*$//;
		    my @columns = split(/\s+/,$line);
		    $wavelengthCount = shift(@columns) if ( $wavelengthCount == -1 );
		    if ( nelem($spectrum) + scalar(@columns) <= $wavelengthCount ) {
			$spectrum = $spectrum->append(pdl @columns);
			$line = "";
		    } else {
			while ( nelem($spectrum) < $wavelengthCount ) {
			    $spectrum = $spectrum->append(shift @columns);
			    $line = join(" ",@columns);
			}
		    }
		}
		# Store the spectrum. (Multiply by wavelength squared in preparation for converting to LSolar/Hz.)
		$spectra(($iMetallicity),($iAge),:) .= $spectrum*$wavelengthDataset**2;

		# Skip the fitting function records.
		my $fitCount = -1;
		my $fit      = pdl [];
		while ( $fitCount == -1 || nelem($fit) < $fitCount ) {
		    ++$iLine;
		    my $line .= <$iHndl>;
		    $line =~ s/^\s*//;
		    $line =~ s/\s*$//;
		    my @columns = split(/\s+/,$line);
		    $fitCount = shift(@columns) if ( $fitCount == -1 );
		    $fit = $fit->append(pdl @columns);
		    $line = "";
		}
	    }
	    close($iHndl);

	}

	# Convert metallicities to logarithmic relative to Solar.
	my $metallicitySolar = pdl 0.0188;
	$metallicities .= log10($metallicities/$metallicitySolar);

	# Convert spectra to correct units. 
	my $angstromsToMeters = pdl 1.0e-10;
	my $speedOfLight      = pdl 2.998e8;
	$spectra             *= $angstromsToMeters/$speedOfLight;

	# Create the HDF5 output file.
	$IMF = ucfirst($IMF);
	my $HDFfile = new PDL::IO::HDF5(">".$galacticusPath."data/stellarPopulations/SSP_Spectra_BC2003_".$resolutionLookup{$resolution}."_imf".$IMF.".hdf5");
	$HDFfile->dataset("ages"         )->set($ageDataset       );
	$HDFfile->dataset("wavelengths"  )->set($wavelengthDataset);
	$HDFfile->dataset("metallicities")->set($metallicities    );
	$HDFfile->dataset("spectra"      )->set($spectra          );
	
	# Write attribute values.
	$HDFfile->group("source")->attrSet(
	    source    => "Bruzual and Charlot, 2003, MNRAS, 344, 1000",
	    sourceURL => "http://adsabs.harvard.edu/abs/2003MNRAS.344.1000B",
	    URL       => "http://www2.iap.fr/users/charlot/bc2003/"
	    );

    }
}

# Write a completion message.
print "Conversion finished. You may delete the ".$dataDirectory." directory now.\n";

exit;
