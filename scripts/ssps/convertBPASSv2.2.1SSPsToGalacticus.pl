#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use DateTime;
$PDL::BIGPDL = 1;

# Convert BPASS v2.2.1 stellar population data to Galacticus' format.
# Andrew Benson (16-February-2022).

# Specify the data directory.
die("Usage: convertBPASSv2.2.1SSPsToGalacticus.pl <dataDirectory>")
    unless ( scalar(@ARGV) == 1 );
my $dataDirectory = $ARGV[0];

# Construct lists of ages and metallicities.
my $metallicities      = pdl [ 1.0e-5, 1.0e-4, 1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2, 1.4e-2, 2.0e-2, 3.0e-2, 4.0e-2 ];
my $metallicitiesSolar = log10($metallicities/0.02);
my $ages               = pdl 10.0**(sequence(51)*0.1-3.0);
my $countMetallicities = nelem($metallicities);
my $countAges          = nelem($ages         );
my $countWavelengths   = 100000;

# Iterate over all IMFs.
opendir(my $models,$dataDirectory);
while ( my $modelName = readdir($models) ) {
    next
	unless ( $modelName =~ m/^BPASSv2.2.1_(bin|sin)\-imf/ );
    my $binarity;
    my $imfLabel;
    if ( $modelName =~ m/^BPASSv2.2.1_(bin|sin)\-imf(\d+(all)??_\d+)/ ) {
	$binarity = $1;
	$imfLabel = $2;
    } elsif ( $modelName =~ m/^BPASSv2.2.1_(bin|sin)\-imf(_chab\d+)/ ) {
	$binarity = $1;
	$imfLabel = $2;
    } else {
	die("can not parse model name");
    }
    next
	if ( -e $ENV{'GALACTICUS_DATA_PATH'}."/dynamic/stellarPopulations/SSP_Spectra_BPASSv2.2.1_".$binarity."-imf".$imfLabel.".hdf5" );
    # Iterate over SSPs.
    my $spectra = pdl zeros($countWavelengths,$countAges,$countMetallicities);
    my $wavelengths = pdl zeros($countWavelengths);
    for(my $iMetallicity=0;$iMetallicity<nelem($metallicities);++$iMetallicity) {
	my $modelFileName = $dataDirectory."/".$modelName."/spectra-".$binarity."-imf".$imfLabel.".z".($metallicities->(($iMetallicity)) < 0.99e-3 ? "em".sprintf("%1.1d",-$metallicities->(($iMetallicity))->log10()->long()->sclr()) : sprintf("%3.3d",long(1000.0*$metallicities->(($iMetallicity))->sclr()))).".dat";
	print $modelFileName."\n"; 
	# Unpack file if necessary.
	system("gunzip ".$modelFileName.".gz")
	    if ( -e $modelFileName.".gz" );
	# Read data.
	my $iWavelength = -1;
	open(my $ssp,$modelFileName);
	for(my $iWavelength=0;$iWavelength<$countWavelengths;++$iWavelength) {
	    my $line = <$ssp>;
	    my @columns = split(" ",$line);
	    $wavelengths->(($iWavelength)) .= $columns[0];
	    $spectra->(($iWavelength),:,($iMetallicity)) .= pdl @columns[1..$countAges];
	}
	close($ssp);
	# Convert units.
	## Spectra are for 10^6Msun burst.
	my $angstromsToMeters = pdl 1.0e-10;
	my $speedOfLight      = pdl 2.998e8;
	for(my $i=0;$i<$countAges;++$i) {
	    $spectra->(:,($i),($iMetallicity)) *= $wavelengths**2*$angstromsToMeters/$speedOfLight/1.0e6;
	}
    }
    # Write data to Galacticus' format file.
    my $outputFile = new PDL::IO::HDF5(">".$ENV{'GALACTICUS_DATA_PATH'}."/dynamic/stellarPopulations/SSP_Spectra_BPASSv2.2.1_".$binarity."-imf".$imfLabel.".hdf5");
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
     $outputFile->dataset("metallicities")->set($metallicitiesSolar);
     $outputFile->dataset("metallicities")->attrSet(
 	units     => "dex",
 	type      => "logarithmic, relative to Solar"
 	);
     $outputFile->dataset("spectra")->set($spectra);
     $outputFile->dataset("spectra")->attrSet(
 	units     => "L☉/Hz",
 	unitsInSI => 3.827e+26
 	);
    
     # Add metadata.
     my @references =
 	(
 	 "Eldridge, Stanway et al. (2017; PASA, 34, 58)",
 	 "Stanway & Eldridge et al. (2018; MNRAS, 479, 75)"
 	);
     my @referenceURLs =
 	(
 	 "https://ui.adsabs.harvard.edu/abs/2017PASA...34...58E",
 	 "https://ui.adsabs.harvard.edu/abs/2018MNRAS.479...75S",
 	);
    my $description = "Simple stellar population spectra from the BPASS v2.2.1 library";
    $description .= $binarity eq "bin" ? "including binary stars" : "not including binary stars";
    my $imf;
    if ( $imfLabel =~ m/chab(\d+)/ ) {
	$description .= " for a Chabrier (2003) IMF with upper mass of ".$1."M☉.";
	$imf          = "Chabrier (2003) IMF with upper mass of ".$1."M☉.";
    } elsif ( $imfLabel =~ m/(\d+)all_(\d+)/ ) {
	my $slope = -$1/100.0-1.0;
	$description .= " for a power law IMF with slope ".sprintf("%5.2f",$slope)." from 0.1M☉ to ".$2."M☉.";
	$imf          = "Power law IMF with slope ".sprintf("%5.2f",$slope)." from 0.1M☉ to ".$2."M☉.";
    } elsif ( $imfLabel =~ m/(\d+)_(\d+)/ ) {
	my $slope1 = -1.30;
	my $slope2 = -$1/100.0-1.0;
	$description .= " for a broken power law IMF with slope ".sprintf("%5.2f",$slope1)." from 0.1M☉ to 0.5M☉ and a slope of ".sprintf("%5.2f",$slope2)." from 0.5M☉ to ".$2."M☉.";
	$imf          = "Broken power law IMF with slope ".sprintf("%5.2f",$slope1)." from 0.1M☉ to 0.5M☉ and a slope of ".sprintf("%5.2f",$slope2)." from 0.5M☉ to ".$2."M☉.";
    } else {
	die('failed to parse IMF label');
    }
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
    my $fileFormat = long(1);
    $outputFile->attrSet(
 	createdBy           => "Andrew Benson <abenson\@obs.carnegiescience.edu>",
 	description         => $description,
 	url                 => "http://bpass.auckland.ac.nz",
 	fileFormat          => $fileFormat,
 	timeStamp           => $now,
 	references          => join("; ",@references   ),
 	referenceURLs       => join("; ",@referenceURLs),
 	initialMassFunction => $imf
 	);
}
closedir($models);

exit;
