#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use DateTime;

# Convert BPASS stellar population data to Galacticus' format.
# Andrew Benson (03-June-2014).

# Specify the data directory.
die("Usage: convertBPASSSSPsToGalacticus.pl <dataDirectory>")
    unless ( scalar(@ARGV) == 1 );
my $dataDirectory = $ARGV[0];

# Download the archived data.
my $fileURL = "http://flexiblelearning.auckland.ac.nz/bpass/2/files/";
my @files = 
    (
     {
	 metallicity => 0.001,
	 fileName => "sed_bpass_z001_tar.gz"
     },
     {
     	 metallicity => 0.004,
     	 fileName => "sed_bpass_z004_tar.gz"
     },
     {
     	 metallicity => 0.008,
     	 fileName => "sed_bpass_z008_tar.gz"
     },
     {
     	 metallicity => 0.020,
     	 fileName => "sed_bpass_z020_tar.gz"
     },
     {
     	 metallicity => 0.040,
     	 fileName => "sed_bpass_z040_tar.gz"
     }
);
foreach ( @files ) {
    system("wget ".$fileURL.$_->{'fileName'}." -O ".$dataDirectory."/".$_->{'fileName'})
	unless ( -e $dataDirectory."/".$_->{'fileName'});
    die('Convert_BPASS_SSPs_to_Galacticus.pl: failed to download file '.$_->{'fileName'})
	unless ( -e $dataDirectory."/".$_->{'fileName'});
}

# Extract the archives.
foreach ( @files ) {
    if ( $_->{'fileName'} =~ m/(z\d+)/ ) {
	$_->{'zLabel'} = $1;
    }
    system("cd ".$dataDirectory."; tar xvfz ".$_->{'fileName'})
	unless ( -e $dataDirectory."/SEDS/sed.bpass.constant.cloudy.bin.".$_->{'zLabel'} );
}

# Extract population ages.
my $ages = pdl [];
open(my $readme,$dataDirectory."/SEDS/sed.bpass.readme.txt");
while ( my $line = <$readme> ) {
    if ( $line =~ m/^\s*\d+\)\s*Flux\( log\(Age\/yrs\)=([\d\.]+)\s*\) \/  L\(Sun\)\/A\./ ) {
	$ages = $ages->append($1);
    }
}
close($readme);
my $ageCount = nelem($ages);
# Ages are giving in log10 years.
$ages = 10.0**($ages-9.0);

# Extract metallicities and convert to log10 relative to Solar.
my $metallicities    = pdl map {$_->{'metallicity'}} @files;
my $metallicitySolar = pdl 0.0188;
$metallicities      .= log10($metallicities/$metallicitySolar);

# Iterate over single and binary populations.
foreach my $population ( "single", "binary" ) {

    # Initialize data sets.
    my $wavelengths;
    my $spectra;

    # Iterate over metallicities.
    my $iFile = -1;
    foreach my $file ( @files ) {
	++$iFile;

	# Read the required data.
	my $fileName = $dataDirectory."/SEDS/sed.bpass.instant.nocont.".substr($population,0,3).".".$file->{'zLabel'};
	$wavelengths = pdl [];
	for(my $i=0;$i<2;++$i) {
	    my $j = -1;
	    open(my $ssp,$fileName);
	    while ( my $line = <$ssp> ) {
		++$j;
		$line =~ s/^\s*//;
		$line =~ s/\s*$//;
		my @columns = split(/\s+/,$line);
		if ( $i == 0 ) {
		    $wavelengths = $wavelengths->append($columns[0]);
		} else {
		    $spectra->(($j),:,($iFile)) .= pdl @columns[1..$ageCount];
		}
	    }
	    close($ssp);
	    $spectra = pdl zeroes(nelem($wavelengths),nelem($ages),scalar(@files))
		if ( $i == 0 && $iFile == 0 );
	}

	# Convert units.
	## Spectra are for 10^6Msun burst.
	my $angstromsToMeters = pdl 1.0e-10;
	my $speedOfLight      = pdl 2.998e8;
	for(my $i=0;$i<$ageCount;++$i) {
	    $spectra->(:,($i),($iFile)) *= $wavelengths**2*$angstromsToMeters/$speedOfLight/1.0e6;
	}

    }

    # Write data to Galacticus' format file.
    my $outputFile = new PDL::IO::HDF5(">".$ENV{'GALACTICUS_DATA_PATH'}."/dynamic/stellarPopulations/SSP_Spectra_BPASS_".$population.".hdf5");
    $outputFile->dataset("ages"         )->set($ages);
    $outputFile->dataset("ages")->attrSet(
	units     => "Gigayears",
	unitsInSI => 3.15576e+16
	);
    $outputFile->dataset("wavelengths")->set($wavelengths);
    $outputFile->dataset("wavelengths")->attrSet(
	units     => "Angstroms",
	unitsInSI => 1.0e-10
	);
    $outputFile->dataset("metallicities")->set($metallicities);
    $outputFile->dataset("metallicities")->attrSet(
	units     => "dex",
	type      => "logarithmic, relative to Solar"
	);
    $outputFile->dataset("spectra")->set($spectra);
    $outputFile->dataset("spectra")->attrSet(
	units     => "Lsolar/Hz",
	unitsInSI => 3.827e+26
	);
    
    # Add metadata.
    my @references =
	(
	 "Eldridge & Stanway, 2012, MNRAS, 419, 479",
	 "Eldridge, Langer & Tout, 2011, MNRAS, 414, 3501",
	 "Eldridge & Stanway, 2009, MNRAS, 400, 1019",
	 "Eldridge, Izzard & Tout, 2008, MNRAS, 384, 1109"
	);
    my @referenceURLs =
	(
	 "http://adsabs.harvard.edu/abs/2012MNRAS.419..479E",
	 "http://adsabs.harvard.edu/abs/2011MNRAS.414.3501E",
	 "http://adsabs.harvard.edu/abs/2009MNRAS.400.1019E",
	 "http://adsabs.harvard.edu/abs/2008MNRAS.384.1109E"
	);
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
    my $fileFormat = long(1);
    $outputFile->attrSet(
	createdBy     => "Andrew Benson <abenson\@obs.carnegiescience.edu>",
	description   => "Simple stellar population spectra from the BPASS library for the canonical BPASS IMF.",
	url           => "http://www.bpass.org.uk/",
	fileFormat    => $fileFormat,
	timeStamp     => $now,
	references    => join("; ",@references   ),
	referenceURLs => join("; ",@referenceURLs),
	initialMassFunction => "Slope of -1.3 between 0.1 and 0.5M☉ and a typical Salpeter slope of -2.35 above 0.5M☉ and a maximum mass of 120M☉."
	);

}

exit;
