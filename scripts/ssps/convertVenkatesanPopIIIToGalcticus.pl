#!/usr/bin/env perl
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::NiceSlice;
use DateTime;

# Convert Aparna Venkatesan's PopIII SEDs to Galacticus' format.
# Andrew Benson (25-June-2012).

# Specify ages.
my $ages          = pdl [ 0.0, 1.0e-2 ];
my $metallicities = pdl [ -4.0 ];

# Specify files to convert.
my @files = (
    {
	file => "data/pop3-inst-10-140.txt",
	out => "SSP_Spectra_Venkatesan_PopIII_m10:140.hdf5",
	massMinimum => "10",
	massMaximum => "140"
    },
    {
	file => "data/pop3-inst.txt",
	out => "SSP_Spectra_Venkatesan_PopIII_m1:100.hdf5",
	massMinimum => "1",
	massMaximum => "100"
    },
    );

# Loop over files.
foreach my $file ( @files ) {
    
    # Read the text file.
    open(iHndl,$file->{'file'});
    my $line = <iHndl>;
    my $wavelength = pdl [];
    my $youngSED   = pdl [];
    my $oldSED     = pdl [];
    while ( my $line = <iHndl> ) {
	$line =~ s/^\s*//;
	$line =~ s/\s*$//;
	my @columns = split(/\s+/,$line);
	$wavelength = $wavelength->append($columns[0]);
	$youngSED   = $youngSED  ->append($columns[1]);
	$oldSED     = $oldSED    ->append($columns[2]);
    }
    close(iHndl);
    
    # Combine spectra into a single PDL.
    my $spectra = pdl zeroes(nelem($wavelength),2,1);
    $spectra->(:,(0),(0)) .= $youngSED*$wavelength**2;
    $spectra->(:,(1),(0)) .= $oldSED  *$wavelength**2;
    
    # Convert spectra to preferred units.
    my $solarLuminosity = pdl 3.845e33; # ergs/s.
    my $speedLight      = pdl 2.998e18; # Angstroms Hz
    $spectra /= 1.0e6; # Convert to 1 Msun population.
    $spectra /= $solarLuminosity*$speedLight;

    # Write the data to file.
    unlink($file->{'out'});
    my $outputFile = new PDL::IO::HDF5(">".$file->{'out'});
    $outputFile->dataset("ages"         )->set($ages         );
    $outputFile->dataset("metallicities")->set($metallicities);
    $outputFile->dataset("wavelengths"  )->set($wavelength   );
    $outputFile->dataset("spectra"      )->set($spectra      );
    
    # Write attribute values.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
    my $fileFormat = long(1);
    $outputFile->attrSet(
	fileFormat  => $fileFormat,
	createdBy   => "Galacticus",
	description => "Simple stellar population spectra for Population III stars from Venkatesan et al. (2003; ApJ; 584; 621; http://adsabs.harvard.edu/abs/2003ApJ...584..621V) using Salpeter IMF from M=".$file->{'massMinimum'}." to ".$file->{'massMaximum'}." M☉",
	timestamp   => $now,
	);
    $outputFile->dataset("ages"         )->attrSet(
	unitsInSI => 3.15576e+16,
	units     => "Gigayears"
	);
    $outputFile->dataset("metallicities")->attrSet(
	units     => "Solar"
	);
    $outputFile->dataset("wavelengths"  )->attrSet(
	unitsInSI => 1.0e-10,
	units     => "Angstroms"
	);
    $outputFile->dataset("spectra"      )->attrSet(
	unitsInSI => 3.827e+33,
	units     => "Lsolar/Hz"
	);
    $outputFile->group("source")->attrSet(
	source    => "Venkatesan et al. (2003; ApJ; 584; 621)",
	sourceURL => "http://adsabs.harvard.edu/abs/2003ApJ...584..621V"
	);
    $outputFile->group("initialMassFunction")->attrSet(
	shape     => "Salpeter",
	massMinimum => $file->{'massMinimum'}." M☉",
	massMaximum => $file->{'massMaximum'}." M☉"
	);

}

exit;
