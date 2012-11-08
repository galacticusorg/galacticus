#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use DateTime;
use File::Copy;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V092'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V092'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");

# Driver script for FSPS_v2.3 (Conroy, White & Gunn stellar population code).
# Andrew Benson (15-Apr-2011)

# Get arguments.
die "Usage: Conroy_SPS_Driver.pl <imfName> <stellarPopulationFile> <fileFormatVersion>" unless ( scalar(@ARGV) == 3 );
my $imfName               = $ARGV[0];
my $stellarPopulationFile = $ARGV[1];
my $fileFormat            = $ARGV[2];

# Ensure the requested file format version is compatible.
my $fileFormatCurrent = pdl long(1);
die('Conroy_SPS_Driver.pl: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
    unless ( $fileFormat == $fileFormatCurrent );

# Determine if we need to make the file.
my $makeFile = 0;
if ( -e $stellarPopulationFile ) {
    my $hdfFile = new PDL::IO::HDF5($stellarPopulationFile);
    my @attributes = $hdfFile->attrs();
    if ( grep {$_ eq "fileFormat"} @attributes  ) {
	my @fileFormatCurrentFile = $hdfFile->attrGet('fileFormat');
	$makeFile = 1 unless ( $fileFormatCurrentFile[0] == $fileFormatCurrent );
    } else {
	$makeFile = 1;
    }
} else {
    $makeFile = 1;
}

# Check if the file exists.
if ( $makeFile == 1 ) {
    
    # Check out the code.
    unless ( -e $galacticusPath."aux/FSPS_v2.3" ) {
 	print "Conroy_SPS_Driver.pl: downloading source code.\n";
 	system("svn checkout http://fsps.googlecode.com/svn/trunk/ ".$galacticusPath."aux/FSPS_v2.3");
 	die("Conroy_SPS_Driver.pl: FATAL - failed to check out svn repository.") unless ( -e $galacticusPath."aux/FSPS_v2.3" );
    }
    
    # Check for updates to the code.
    my $availableRevision;
    my $currentRevision;
    open(pHndl,"svn info -r HEAD ".$galacticusPath."aux/FSPS_v2.3 |");
    while ( my $line = <pHndl> ) {
 	if ( $line =~ m/Last Changed Rev:\s*(\d+)/ ) {$availableRevision = $1};
    }
    close(pHndl);
    open(pHndl,"svn info ".$galacticusPath."aux/FSPS_v2.3 |");
    while ( my $line = <pHndl> ) {
 	if ( $line =~ m/Last Changed Rev:\s*(\d+)/ ) {$currentRevision = $1};
    }
    close(pHndl);
    if ( $currentRevision < $availableRevision ) {
 	print "Conroy_SPS_Driver.pl: updating source code.\n";
 	system("svn revert -R ".$galacticusPath."aux/FSPS_v2.3"); # Revert the code.
 	system("svn update ".$galacticusPath."aux/FSPS_v2.3"); # Grab updates
 	unlink($galacticusPath."aux/FSPS_v2.3/src/galacticus_IMF.f90") # Remove this file to trigger re-patching of the code.
    }
    
    # Patch the code.
    unless ( -e $galacticusPath."aux/FSPS_v2.3/src/galacticus_IMF.f90" ) {
 	foreach my $file ( "galacticus_IMF.f90", "imf.f90.patch", "Makefile.patch", "ssp_gen.f90.patch", "autosps.f90.patch", "qromb.f90.patch" ) {
	    my $todir = $galacticusPath."aux/FSPS_v2.3/src/";
	    $todir .= "nr/"
		if ( $file eq "qromb.f90.patch" );
 	    copy($galacticusPath."aux/FSPS_v2.3_Galacticus_Modifications/".$file,$todir.$file);
 	    if ( $file =~ m/\.patch$/ ) {
		system("cd ".$todir."; patch < $file");
		die("Conroy_SPS_Driver.pl: unable to patch file: ".$file) unless ( $? == 0 );
	    }
 	    print "$file\n";
 	}
 	unlink($galacticusPath."aux/FSPS_v2.3/src/autosps.exe");
    }
    
    # Build the code.
    unless ( -e $galacticusPath."aux/FSPS_v2.3/src/autosps.exe" ) {
 	print "Conroy_SPS_Driver.pl: compiling autosps.exe code.\n";
 	system("cd ".$galacticusPath."aux/FSPS_v2.3/src; export SPS_HOME=`pwd`; make clean; make -j 1");
 	die("Conroy_SPS_Driver.pl: FATAL - failed to build autosps.exe code.") unless ( -e $galacticusPath."aux/FSPS_v2.3/src/autosps.exe" );
    }
    
    # Initialize the data structure.
    my $spectra;
    my $wavelengths;
    my $ages          = pdl [];
    my $metallicities = pdl [];
    
    # Open the output file.
    my $hdfFile = new PDL::IO::HDF5(">".$stellarPopulationFile);

    # Run the code.
    my $pwd = `pwd`;
    chomp($pwd);
    $ENV{'SPS_HOME'} = $galacticusPath."aux/FSPS_v2.3";
    my $iMetallicity = -1;

    # Add a description and other metadata to the file.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
    $hdfFile->attrSet(
	description   => "Simple stellar population spectra from Conroy, White & Gunn for a ".$imfName." initial mass function",
	timestamp     => $now,
	fspsVersion   => "2.3_r".$availableRevision,
	createdBy     => "Galacticus",
	fileFormat    => $fileFormatCurrent
	);

    # Read the IMF file and store it to our output file.
    if ( -e "galacticus.imf" ) {
	my $imfMass = pdl [];
	my $imfPhi  = pdl [];
	open(imfHndl,"galacticus.imf");
	while ( my $line = <imfHndl> ) {
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    my @columns = split(/\s+/,$line);
	    $imfMass = $imfMass->append($columns[0]);
	    $imfPhi  = $imfPhi ->append($columns[1]);
	}
	close(imfHndl);
	my $imfGroup = new PDL::IO::HDF5::Group(
	    name    => "initialMassFunction",
	    parent  => $hdfFile,
	    fileObj => $hdfFile
	    );
	my $massDataSet = new PDL::IO::HDF5::Dataset(
	    name    => "mass",
	    parent  => $imfGroup,
	    fileObj => $hdfFile
	    );
	$massDataSet->set($imfMass);
	$massDataSet->attrSet(
	    units     => "Msolar",
	    unitsInSI => 1.98892e30
	    );
 	my $imfDataSet = new PDL::IO::HDF5::Dataset(
	    name    => "initialMassFunction",
	    parent  => $imfGroup,
	    fileObj => $hdfFile
	    );
	$imfDataSet->set($imfPhi);
	$imfDataSet->attrSet(
	    units     => "Msolar^-1",
	    unitsInSI => 0.50278543e-30
	    );
    } else {
	die("Conroy_SPS_Driver.pl: 'galacticus.imf' file is missing");
    }

    # Loop over metallicities.
    for(my $iZ=1;$iZ<=22;++$iZ) {
	my $outFile = "imf".$imfName.".iZ".$iZ;
	unless ( -e $galacticusPath."aux/FSPS_v2.3/OUTPUTS/".$outFile.".spec" ) {
	    open(spsPipe,"|aux/FSPS_v2.3/src/autosps.exe");
	    print spsPipe "6\n";         # IMF.
	    print spsPipe "0\n";         # Generate SSP.
	    print spsPipe $iZ."\n";      # Specify metallicity.
	    print spsPipe "no\n";        # Do not include dust.
	    print spsPipe $outFile."\n"; # Specify filename.
	    close(spsPipe);
	}
	
	# Read the file.
	my $ageCount =  0;
	my $iAge     = -1;
	my $gotAge   =  0;
	open(specFile,$galacticusPath."aux/FSPS_v2.3/OUTPUTS/".$outFile.".spec");
	while ( my $line = <specFile> ) {
	    chomp($line);
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    if ( $line =~ m/^\#/ ) {
		# Comment line.
		if ( $line =~ m/Log\(Z\/Zsol\):\s*([\+\-\.0-9]+)/ ) {
		    my $metallicity = $1;
		    print "  -> Reading data for metallicity log10(Z/Z_Solar) = ".$metallicity."\n";
		    ++$iMetallicity;
		    $metallicities = $metallicities->append($metallicity);
		}
	    } else {
		if ( $ageCount == 0 ) {
		    # Have not yet got a count of the number of ages in the file - get it now.
		    if ( $line =~ m/^\s*(\d+)\s+(\d+)/ ) {
			$ageCount           = $1;
			my $wavelengthCount = $2;
			$spectra            = pdl zeroes($wavelengthCount,$ageCount,22) unless ( defined($spectra) );
			print "     -> Found ".$ageCount." ages in the file\n";
			print "     -> Found ".$wavelengthCount." wavelengths in the file\n";
			# Read the wavelength array.
			$line        = <specFile>;
			chomp($line);
			$line        =~ s/^\s*//;
			$line        =~ s/\s*$//;
			my @columns  = split(/\s+/,$line);
			$wavelengths = pdl @columns unless ( defined($wavelengths) );
		    } else {
			die("Conroy_SPS_Driver.pl: format of output file is not recognized");
		    }
		} elsif ( $gotAge == 0 ) {
		    # Not currently processing an SPS age - get the age.
		    ++$iAge;
		    my @columns = split(/\s+/,$line);
		    $ages       = $ages->append(10.0**($columns[0]-9.0)) if ( $iMetallicity == 0 );
		    $gotAge     = 1;
		} else {
		    # Are processing an SPS age - grab the wavelengths and append to the array.
		    my @columns = split(/\s+/,$line);
		    $spectra(:,($iAge),($iMetallicity)) .= pdl @columns;
		    $gotAge                              = 0;
		}
	    }
	}
	close(specFile);
    }
    
    # Wavelengths.
    my $wavelengthsDataSet = new PDL::IO::HDF5::Dataset(
	name    => "wavelengths",
	parent  => $hdfFile,
	fileObj => $hdfFile
	);
    $wavelengthsDataSet->set($wavelengths);
    $wavelengthsDataSet->attrSet(
	units     => "Angstroms",
	unitsInSI => 1.0e-10
	);

    # Ages
    my $agesDataSet = new PDL::IO::HDF5::Dataset(
	name    => "ages",
	parent  => $hdfFile,
	fileObj => $hdfFile
	);
    $agesDataSet->set($ages);
    $agesDataSet->attrSet(
	units     => "Gigayears",
	unitsInSI => 3.15576e16
	);

    # Metallicities.
    my $metallicitiesDataSet = new PDL::IO::HDF5::Dataset(
	name    => "metallicities",
	parent  => $hdfFile,
	fileObj => $hdfFile
	);
    $metallicitiesDataSet->set($metallicities);
    
    # Spectra.
    my $spectraDataSet = new PDL::IO::HDF5::Dataset(
	name    => "spectra",
	parent  => $hdfFile,
	fileObj => $hdfFile
	);
    $spectraDataSet->set($spectra);
    $spectraDataSet->attrSet(
	units     => "Lsolar/Hz",
	unitsInSI => 3.827e33
	);

}
exit;

unlink("galacticus.imf");

exit;
