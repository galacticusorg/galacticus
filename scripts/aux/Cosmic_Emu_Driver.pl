#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";

# Driver script for Cosmic Emu power spectrum emulator.
# Andrew Benson (10-July-2012)

# Get arguments.
die "Usage: Cosmic_Emu_Driver.pl <parameterFile> <powerSpectrumFile>"
    unless ( scalar(@ARGV) == 2 );
my $parameterFile     = $ARGV[0];
my $powerSpectrumFile = $ARGV[1];

# Create a directory for the emulator.
system("mkdir -p ". $ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu");

# Check for presence of the executable.
unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/emu.exe" ) {
    
    # Check for presence of the source code.
    unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/emu.c" ) {
	
	# Download the code.
	unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/CosmicEmu_v1.1.tar.gz" ) {
	    print "Cosmic_Emu_Driver.pl: downloading Cosmic_Emu code.\n";
	    system("wget http://www.hep.anl.gov/cosmology/CosmicEmu/CosmicEmu_v1.0.tar.gz -O ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/CosmicEmu_v1.1.tar.gz");
	    die("Cosmic_Emu_Driver.pl: FATAL - failed to download Cosmic_Emu code.")
		unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/CosmicEmu_v1.1.tar.gz" );
	}

	# Unpack the code.
	print "Cosmic_Emu_Driver.pl: unpacking Cosmic_Emu code.\n";
	system("tar -x -v -z -C ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu -f ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/CosmicEmu_v1.1.tar.gz");
	die("Cosmic_Emu_Driver.pl: FATAL - failed to unpack Cosmic_Emu code.")
	    unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/CosmicEmu_v1.0/emu.c" );
    }

    # Build the code.
    print "Cosmic_Emu_Driver.pl: compiling Cosmic_Emu code.\n";
    system("cd ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/CosmicEmu_v1.0/; make");
    die("Cosmic_Emu_Driver.pl: FATAL - failed to build Cosmic_Emu code.")
	unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/CosmicEmu_v1.0/emu.exe" );
}

# Parse the parameter file.
my $xml           = new XML::Simple;
my $parameterData = $xml->XMLin($parameterFile);

# Check that required parameters exist.
my @parameters = ( "OmegaBaryon", "OmegaMatter", "HubbleConstant", "sigma_8", "powerSpectrumIndex", "darkEnergyEquationOfState", "redshift" );
foreach my $parameter ( @parameters ) {
    die("Cosmic_Emu_Driver.pl: FATAL - parameter ".$parameter." can not be found.") 
	unless ( exists($parameterData->{$parameter}) );
}

# Calculate derived parameters.
my $omegaMatter = $parameterData->{'OmegaMatter'}->{'value'}*($parameterData->{'HubbleConstant'}->{'value'}/100.0)**2;
my $omegaBaryon = $parameterData->{'OmegaBaryon'}->{'value'}*($parameterData->{'HubbleConstant'}->{'value'}/100.0)**2;

# Run Cosmic_Emu.
open(emuPipe,"|".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/CosmicEmu/CosmicEmu_v1.0/emu.exe");
print emuPipe $powerSpectrumFile."\n";
print emuPipe $omegaMatter."\n";
print emuPipe $omegaBaryon."\n";
print emuPipe $parameterData->{'powerSpectrumIndex'       }->{'value'}."\n";
print emuPipe $parameterData->{'sigma_8'                  }->{'value'}."\n";
print emuPipe $parameterData->{'darkEnergyEquationOfState'}->{'value'}."\n";
print emuPipe $parameterData->{'redshift'                 }->{'value'}."\n";
print emuPipe "2\n";
close(emuPipe);

exit;
