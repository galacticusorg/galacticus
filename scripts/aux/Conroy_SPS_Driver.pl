#!/usr/bin/env perl
use XML::Simple;
use Data::Dumper;
use File::Copy;

# Driver script for FSPS_v2.2 (Conroy, White & Gunn stellar population code).
# Andrew Benson (15-Apr-2011)

# Get arguments.
if ( $#ARGV != 1 ) {die "Usage: Conroy_SPS_Driver.pl <imfName> <stellarPopulationFile>"};
$imfName               = $ARGV[0];
$stellarPopulationFile = $ARGV[1];

# Check if the file exists.
unless ( -e $stellarPopulationFile ) {
    
    # Check out the code.
    unless ( -e "aux/FSPS_v2.2" ) {
	print "Conroy_SPS_Driver.pl: downloading source code.\n";
	system("svn checkout http://fsps.googlecode.com/svn/trunk/ aux/FSPS_v2.2");
	die("Conroy_SPS_Driver.pl: FATAL - failed to check out svn repository.") unless ( -e "aux/FSPS_v2.2" );
    }

    # Check for updates to the code.
    open(pHndl,"svn info -r HEAD aux/FSPS_v2.2 |");
    while ( $line = <pHndl> ) {
	if ( $line =~ m/Last Changed Rev:\s*(\d+)/ ) {$availableRevision = $1};
    }
    close(pHndl);
    open(pHndl,"svn info aux/FSPS_v2.2 |");
    while ( $line = <pHndl> ) {
	if ( $line =~ m/Last Changed Rev:\s*(\d+)/ ) {$currentRevision = $1};
    }
    close(pHndl);
    if ( $currentRevision < $availableRevision ) {
	print "Conroy_SPS_Driver.pl: updating source code.\n";
	system("svn revert aux/FSPS_v2.2"); # Revert the code.
	system("svn update aux/FSPS_v2.2"); # Grab updates
	unlink("aux/FSPS_v2.2/src/galacticus_IMF.f90") # Remove this file to trigger re-patching of the code.
    }

    # Patch the code.
    unless ( -e "aux/FSPS_v2.2/src/galacticus_IMF.f90" ) {
	foreach $file ( "galacticus_IMF.f90", "imf.f90.patch", "Makefile.patch", "ssp_gen.f90.patch", "autosps.f90.patch" ) {
	    copy("aux/FSPS_v2.2_Galacticus_Modifications/".$file,"aux/FSPS_v2.2/src/".$file);
	    if ( $file =~ m/\.patch$/ ) {system("cd aux/FSPS_v2.2/src; patch < $file")};
	    print "$file\n";
	}
	unlink("aux/FSPS_v2.2/src/autosps.exe");
    }

    # Build the code.
    unless ( -e "aux/FSPS_v2.2/src/autosps.exe" ) {
	print "Conroy_SPS_Driver.pl: compiling autosps.exe code.\n";
	system("cd aux/FSPS_v2.2/src; export SPS_HOME=`pwd`; make");
	die("Conroy_SPS_Driver.pl: FATAL - failed to build autosps.exe code.") unless ( -e "aux/FSPS_v2.2/src/autosps.exe" );
    }
    exit;
    # Read the wavelength array.
    undef(%data);
    open(lambdaFile,"aux/FSPS_v2.2/SPECTRA/BaSeL3.1/basel.lambda");
    while ( $line = <lambdaFile> ) {
	chomp($line);
	$line =~ s/^\s*//;
	$line =~ s/\s*$//;
	push(@{$data{'wavelengths'}->{'wavelength'}},$line);
    }
    close(lambdaFile);

    # Run the code.
    $pwd = `pwd`;
    chomp($pwd);
    $ENV{'SPS_HOME'} = $pwd."/aux/FSPS_v2.2";
    $iMetallicity = -1;

    # Add a description of the file.
    $data{'description'} = "Simple stellar population spectra from Conroy, White & Gunn for an ".$imfName." initial mass function";

    # Add a description of the units used to the file.
    $data{'units'}->{'value'} = 3.827e33; # Values are in Solar luminosities per Hz - this is the Solar luminosity (in erg/s) used.
    $data{'units'}->{'system'} = "ergs/s/Hz";

    # Loop over metallicities.
    for($iZ=1;$iZ<=22;++$iZ) {
	$outFile = "imf".$imfName.".iZ".$iZ;
	unless ( -e "aux/FSPS_v2.2/OUTPUTS/".$outFile.".spec" ) {
	    open(spsPipe,"|aux/FSPS_v2.2/src/autosps.exe");
	    print spsPipe "6\n";        # IMF.
	    print spsPipe "0\n";        # Generate SSP.
	    print spsPipe "$iZ\n";      # Specify metallicity.
	    print spsPipe "no\n";       # Do not include dust.
	    print spsPipe "$outFile\n"; # Specify filename.
	    close(spsPipe);
	}

	# Read the file and convert it to XML format.
	$ageCount = 0;
	$iAge     = -1;
	$gotAge   = 0;
	open(specFile,"aux/FSPS_v2.2/OUTPUTS/".$outFile.".spec");
	while ( $line = <specFile> ) {
	    chomp($line);
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    if ( $line =~ m/^\#/ ) {
		# Comment line.
		if ( $line =~ m/Log\(Z\/Zsol\):\s*([\+\-\.0-9]+)/ ) {
		    $metallicity = $1;
		    print "  -> Reading data for metallicity log10(Z/Z_Solar) = $metallicity\n";
		    ++$iMetallicity;
		    ${$data{'metallicity'}}[$iMetallicity]->{'value'} = $metallicity;
		}
	    } else {
		if ( $ageCount == 0 ) {
		    # Have not yet got a count of the number of ages in the file - get it now.
		    $ageCount = $line;
		    print "     -> Found $ageCount ages in the file\n";
		} elsif ( $gotAge == 0 ) {
		    # Not currently processing an SPS age - get the age.
		    ++$iAge;
		    @columns = split(/\s+/,$line);
		    ${$data{'metallicity'}}[$iMetallicity]->{'age'}[$iAge]->{'value'} = 10.0**($columns[0]-9.0);
		    $gotAge = 1;
		} else {
		    # Are processing an SPS age - grab the wavelengths and append to the array.
		    @columns = split(/\s+/,$line);
		    push(@{${$data{'metallicity'}}[$iMetallicity]->{'age'}[$iAge]->{'flux'}},@columns);
		    # Check if we've read enough entries (to match the wavelength grid).
		    if ( $#{${$data{'metallicity'}}[$iMetallicity]->{'age'}[$iAge]->{'flux'}} == $#{$data{'wavelengths'}->{'wavelength'}} ) {$gotAge = 0};
		}
	    }
	}
	close(specFile);
    }

    # Construct an HDF5 file of the output.
    $tmpDataFile = "data.tmp";
    $tmpConfigFile = "config.tmp";

    # Wavelengths.
    open(tmpFile,">".$tmpDataFile);
    print tmpFile join(" ",@{$data{'wavelengths'}->{'wavelength'}});
    close(tmpFile);
    open(configFile,">".$tmpConfigFile);
    print configFile "PATH wavelengths\n";
    print configFile "INPUT-CLASS TEXTFP\n";
    print configFile "RANK 1\n";
    print configFile "DIMENSION-SIZES ".$#{$data{'wavelengths'}->{'wavelength'}}."\n";
    print configFile "OUTPUT-CLASS FP\n";
    print configFile "OUTPUT-SIZE 64\n";
    print configFile "OUTPUT-ARCHITECTURE NATIVE\n";
    close(configFile);
    system("h5import $tmpDataFile -c $tmpConfigFile -o $stellarPopulationFile");

    # Ages
    open(tmpFile,">".$tmpDataFile);
    for($iAge=0;$iAge<$ageCount;++$iAge) {
	print tmpFile ${$data{'metallicity'}}[0]->{'age'}[$iAge]->{'value'}."\n";
    }
    close(tmpFile);
    open(configFile,">".$tmpConfigFile);
    print configFile "PATH ages\n";
    print configFile "INPUT-CLASS TEXTFP\n";
    print configFile "RANK 1\n";
    print configFile "DIMENSION-SIZES ".$ageCount."\n";
    print configFile "OUTPUT-CLASS FP\n";
    print configFile "OUTPUT-SIZE 64\n";
    print configFile "OUTPUT-ARCHITECTURE NATIVE\n";
    close(configFile);
    system("h5import $tmpDataFile -c $tmpConfigFile -o $stellarPopulationFile");

    # Metallicities.
    open(tmpFile,">".$tmpDataFile);
    for($iMetallicity=0;$iMetallicity<=21;++$iMetallicity) {
	print tmpFile ${$data{'metallicity'}}[$iMetallicity]->{'value'}."\n";
    }
    close(tmpFile);
    open(configFile,">".$tmpConfigFile);
    print configFile "PATH metallicities\n";
    print configFile "INPUT-CLASS TEXTFP\n";
    print configFile "RANK 1\n";
    print configFile "DIMENSION-SIZES 22\n";
    print configFile "OUTPUT-CLASS FP\n";
    print configFile "OUTPUT-SIZE 64\n";
    print configFile "OUTPUT-ARCHITECTURE NATIVE\n";
    close(configFile);
    system("h5import $tmpDataFile -c $tmpConfigFile -o $stellarPopulationFile");

    # Spectra.
    open(tmpFile,">".$tmpDataFile);
    for($iWavelength=0;$iWavelength<$#{$data{'wavelengths'}->{'wavelength'}};++$iWavelength) {
	for($iAge=0;$iAge<$ageCount;++$iAge) {
	    for($iMetallicity=0;$iMetallicity<=21;++$iMetallicity) {
		print tmpFile ${${$data{'metallicity'}}[$iMetallicity]->{'age'}[$iAge]->{'flux'}}[$iWavelength]." ";
	    }
	    print tmpFile "\n";
	}
    }
    close(tmpFile);
    open(configFile,">".$tmpConfigFile);
    print configFile "PATH spectra\n";
    print configFile "INPUT-CLASS TEXTFP\n";
    print configFile "RANK 3\n";
    print configFile "DIMENSION-SIZES ".$#{$data{'wavelengths'}->{'wavelength'}}." ".$ageCount." 22\n";
    print configFile "OUTPUT-CLASS FP\n";
    print configFile "OUTPUT-SIZE 64\n";
    print configFile "OUTPUT-ARCHITECTURE NATIVE\n";
    close(configFile);
    system("h5import $tmpDataFile -c $tmpConfigFile -o $stellarPopulationFile");

    unlink($tmpDataFile);
    unlink($tmpConfigFile);
}

unlink("galacticus.imf");

exit;
