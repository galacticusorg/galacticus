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

# Convert Maraston stellar population data to Galacticus' format.
# Andrew Benson (18-October-2010).

# Define the base URL for the data.
my $baseURL = "http://www-astro.physics.ox.ac.uk/~maraston/SSPn/SED/";

# Specify the files to be downloaded.
my @downloadFiles = (
    "AgegridSSP_Mar05",
    "Sed_Mar05_SSP_Kroupa.tar.gz",
    "Sed_Mar05_SSP_Salpeter.tar.gz"
    );

# Create a data directory.
my $dataDirectory = "SSP_Maraston";
system("mkdir -p ".$dataDirectory);

# Download data.
foreach my $file ( @downloadFiles ) {
    system("wget ".$baseURL.$file." -O ".$dataDirectory."/".$file);
    system("cd ".$dataDirectory."; tar xvfz ".$file) if ( $file =~ m/\.tar\.gz$/ );
}

# Specify list of IMFs to convert.
my %IMFs = (
    "Kroupa"   => { label => "kr" },
    "Salpeter" => { label => "ss" }
    );

# Specify horizontal branch morphologies to convert.
my %hbMorphologies = (
    #"Blue" => { label => "bhb" }, # Ignore blue horizontal branch files as they contain only two ages.
    "Red"  => { label => "rhb" }
    );

# Specify list of metallicities. (Exclude those for which only crude time grids are available.)
my @metallicities     = ( "z0001", "z001", "z002", "z004" );
my %metallicityValues = (
    "z007"  => +0.67,
    "z004"  => +0.35,
    "z002"  => +0.00,
    "z001"  => -0.33,
    "z0001" => -1.35,
    "z10m4" => -2.25
    );

# Count lines in the age grid file.
my $ageCount = 0;
open(my $iHndl,$dataDirectory."/AgegridSSP_Mar05");
while ( my $line = <$iHndl> ) {
    ++$ageCount;
}
close($iHndl);

# Loop over all IMFs.
foreach my $IMF ( keys(%IMFs) ) {

    # Loop over all horizontal branch morphologies.
    foreach my $hbMorphology ( keys(%hbMorphologies) ) {

	# Initialize data.
	my $fluxData;

        # Loop over all metallicities.
	my $metallicityData = pdl [];
	my $iMetal = -1;
	my $ages    = pdl [];
	my $lambdas = pdl [];
	foreach my $metallicity ( @metallicities ) {
	    ++$iMetal;
	    $metallicityData = $metallicityData->append($metallicityValues{$metallicity});

	    # Construct file name.
	    my $fileName = $dataDirectory."/sed.".$IMFs{$IMF}->{'label'}.$metallicity.".".$hbMorphologies{$hbMorphology}->{'label'};

	    # Open the file and read data.
	    my $lastAge = -1.0;
	    my $iAge    = -1;
	    my $iLine = 0;
	    my $fluxes;
	    open(my $iHndl,$fileName);
	    while ( my $line = <$iHndl> ) {
		next if ( $line =~ m/^\s*$/ );
		++$iLine;
		$line =~ s/^\s*//;
		$line =~ s/\s*$//;
		my @columns = split(/\s+/,$line);
		my $age    = $columns[0];
		my $lambda = $columns[2];
		my $flux   = $columns[3];
		unless ( $age == $lastAge ) {
		    if ( nelem($lambdas) > 0 ) {
			++$iAge;
			$fluxData = pdl zeroes(scalar(@metallicities),$ageCount,nelem($fluxes)) unless (defined($fluxData));
			$fluxData(($iMetal),($iAge),:) .= $fluxes;
		    }
		    $lambdas = pdl [];
		    $fluxes  = pdl [];
		    $ages = $ages->append($age);
		    $lastAge = $age;
		}
		$lambdas = $lambdas->append($lambda);
		$fluxes  = $fluxes ->append($flux*$lambda**2);
	    }
	    close($iHndl);
	    if (defined($lambdas)) {
		++$iAge;
		$fluxData(($iMetal),($iAge),:) .= $fluxes;
	    }

	}

	# Convert fluxes to Lsun/Hz.
	my $solarLuminosity   = pdl 3.826e33;
	my $angstromsToMeters = pdl 1.0e-10;
	my $speedOfLight      = pdl 2.998e8;
	$fluxData *= $angstromsToMeters/$solarLuminosity/$speedOfLight;

	# Create the HDF5 output file.
	my $HDFfile = new PDL::IO::HDF5(">".$galacticusPath."data/stellarPopulations/SSP_Spectra_Maraston_hbMorphology".$hbMorphology."_imf".$IMF.".hdf5");
	$HDFfile->dataset("ages"         )->set($ages           );
	$HDFfile->dataset("wavelengths"  )->set($lambdas        );
	$HDFfile->dataset("metallicities")->set($metallicityData);
	$HDFfile->dataset("spectra"      )->set($fluxData       );

	# Write attribute values.
	$HDFfile->group("source")->attrSet(
	    source    => "Maraston, C. 2005, MNRAS, 362, 799; Maraston, C. 1998, MNRAS, 300, 872",
	    sourceURL => "http://adsabs.harvard.edu/abs/2005MNRAS.362..799M http://adsabs.harvard.edu/abs/1998MNRAS.300..872M",
	    URL       => $baseURL."Claudia%27s_Stellar_Population_Models.html"
	);

    }

}

# Remove the temporary data directory.
system("rm -rf ".$dataDirectory);

exit;
