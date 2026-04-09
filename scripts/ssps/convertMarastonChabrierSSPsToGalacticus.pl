#!/usr/bin/env perl
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::NiceSlice;

# Convert Maraston stellar population data to Galacticus' format.
# Andrew Benson (18-October-2010).

# Create a data directory.
$dataDirectory = "SSP_Maraston";
system("mkdir -p ".$dataDirectory);

# Specify list of metallicities. (Exclude those for which only crude time grids are available.)
@metallicities     = ( "m0.0010", "m0.0100", "m0.0200", "m0.0400" );
%metallicityValues = (
    "m0.0400" => +0.35,
    "m0.0200" => +0.00,
    "m0.0100" => -0.33,
    "m0.0010" => -1.35
    );

# Specify number of ages in the files.
$ageCount = 220;

# Loop over all metallicities.
$metallicityData = pdl [];
$iMetal = -1;
foreach $metallicity ( @metallicities ) {
    ++$iMetal;
    $metallicityData = $metallicityData->append($metallicityValues{$metallicity});
    
    # Clear data.
    $ages    = pdl [];
    $lambdas = pdl [];
    
    # Construct file name.
    $fileName = $dataDirectory."/M05_Chabrier_FullSED_".$metallicity.".dat";
    print $fileName."\n";
    
    # Open the file and read data.
    $lastAge = -1.0;
    $iAge    = -1;
    $iLine = 0;
    open(iHndl,$fileName);
    while ( $line = <iHndl> ) {
	next if ( $line =~ m/^\s*$/ );
	++$iLine;
	$line =~ s/^\s*//;
	$line =~ s/\s*$//;
	@columns = split(/\s+/,$line);
	$age    = $columns[0];
	$lambda = $columns[2];
	$flux   = $columns[3];
	unless ( $age == $lastAge ) {
	    if ( nelem($lambdas) > 0 ) {
		++$iAge;
		$fluxData = pdl zeroes(nelem($fluxes),$ageCount,$#metallicities+1) unless (defined($fluxData));
		$fluxData(:,($iAge),($iMetal)) .= $fluxes;
	    }
	    $lambdas = pdl [];
	    $fluxes  = pdl [];
	    $ages    = $ages->append($age);
	    $lastAge = $age;
	}
	$lambdas = $lambdas->append($lambda           );
	$fluxes  = $fluxes ->append($flux  *$lambda**2);
    }
    close(iHndl);
    if (defined($lambdas)) {
	++$iAge;
	$fluxData(:,($iAge),($iMetal)) .= $fluxes;
    }
    
}

# Convert fluxes to Lsolar/Hz.
my $solarLuminosity = pdl 3.826e33;
my $angstroms       = pdl 1.0e-10;
my $speedLight      = pdl 2.998e8;
$fluxData *= $angstroms/$speedLight/$solarLuminosity;

# Convert ages to Gyr.
$ages /= 1.0e9;

# Create the HDF5 output file.
$HDFfile = new PDL::IO::HDF5(">data/SSP_Spectra_Maraston_imfChabrier.hdf5");
$HDFfile->dataset("ages"         )->set($ages           );
$HDFfile->dataset("ages")->attrSet(
    units     => "Gigayears",
    unitsInSI => 3.15576e+16
    );
$HDFfile->dataset("wavelengths"  )->set($lambdas        );
$HDFfile->dataset("wavelengths")->attrSet(
    units     => "Angstroms",
    unitsInSI => 1.0e-10
    );
$HDFfile->dataset("metallicities")->set($metallicityData);
$HDFfile->dataset("metallicities")->attrSet(
    units     => "dex",
    type      => "logarithmic, relative to Solar"
    );
$HDFfile->dataset("spectra"      )->set($fluxData       );
$HDFfile->dataset("spectra")->attrSet(
    units     => "Lsolar/Hz",
    unitsInSI => 3.827e+26
    );
my $format = pdl long(1);

# Write attribute values.
$HDFfile->group("source")->attrSet(
    source    => "Maraston, C. 2005, MNRAS, 362, 799; Maraston, C. 1998, MNRAS, 300, 872",
    sourceURL => "http://adsabs.harvard.edu/abs/2005MNRAS.362..799M http://adsabs.harvard.edu/abs/1998MNRAS.300..872M",
    URL       => $baseURL."Claudia%27s_Stellar_Population_Models.html",
    provenance => "Provided by Bruno Henriques <bhenriques@mpa-garching.mpg.de>; April 8, 2015"
    );
$HDFfile->attrSet(
    fileFormat => $format
    );

exit;
