#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V091'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V091'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");

# Driver script for Cosmic Emu power spectrum emulator.
# Andrew Benson (10-July-2012)

# Get arguments.
die "Usage: Cosmic_Emu_Driver.pl <parameterFile> <powerSpectrumFile>"
    unless ( scalar(@ARGV) == 2 );
my $parameterFile     = $ARGV[0];
my $powerSpectrumFile = $ARGV[1];

# Create a directory for the emulator.
system("mkdir -p ". $galacticusPath."aux/CosmicEmu");

# Check for presence of the executable.
unless ( -e $galacticusPath."aux/CosmicEmu/emu.exe" ) {
    
    # Check for presence of the source code.
    unless ( -e $galacticusPath."aux/CosmicEmu/emu.c" ) {
	
	# Download the code.
	unless ( -e $galacticusPath."aux/CosmicEmu/CosmicEmu_v1.1.tar.gz" ) {
	    print "Cosmic_Emu_Driver.pl: downloading Cosmic_Emu code.\n";
	    system("wget http://www.lanl.gov/projects/cosmology/CosmicEmu/CosmicEmu_v1.1.tar.gz -O ".$galacticusPath."aux/CosmicEmu/CosmicEmu_v1.1.tar.gz");
	    die("Cosmic_Emu_Driver.pl: FATAL - failed to download Cosmic_Emu code.")
		unless ( -e $galacticusPath."aux/CosmicEmu/CosmicEmu_v1.1.tar.gz" );
	}

	# Unpack the code.
	print "Cosmic_Emu_Driver.pl: unpacking Cosmic_Emu code.\n";
	system("tar -x -v -z -C ".$galacticusPath."aux/CosmicEmu -f ".$galacticusPath."aux/CosmicEmu/CosmicEmu_v1.1.tar.gz");
	die("Cosmic_Emu_Driver.pl: FATAL - failed to unpack Cosmic_Emu code.")
	    unless ( -e $galacticusPath."aux/CosmicEmu/emu.c" );
    }

    # Build the code.
    print "Cosmic_Emu_Driver.pl: compiling Cosmic_Emu code.\n";
    system("cd ".$galacticusPath."aux/CosmicEmu/; make");
    die("Cosmic_Emu_Driver.pl: FATAL - failed to build Cosmic_Emu code.")
	unless ( -e $galacticusPath."aux/CosmicEmu/emu.exe" );
}

# Parse the parameter file.
my $xml           = new XML::Simple;
my $parameterData = $xml->XMLin($parameterFile);
my $parameters    = $parameterData->{'parameter'};

# Check that required parameters exist.
my @parameters = ( "Omega_b", "Omega_Matter", "H_0", "sigma_8", "powerSpectrumIndex", "darkEnergyEquationOfState", "redshift" );
foreach my $parameter ( @parameters ) {
    die("Cosmic_Emu_Driver.pl: FATAL - parameter ".$parameter." can not be found.") 
	unless ( exists($parameters->{$parameter}) );
}

# Calculate derived parameters.
my $omegaMatter = $parameters->{'Omega_Matter'}->{'value'}*($parameters->{'H_0'}->{'value'}/100.0)**2;
my $omegaBaryon = $parameters->{'Omega_b'     }->{'value'}*($parameters->{'H_0'}->{'value'}/100.0)**2;

# Run Cosmic_Emu.
open(emuPipe,"|".$galacticusPath."aux/CosmicEmu/emu.exe");
print emuPipe $powerSpectrumFile."\n";
print emuPipe $omegaMatter."\n";
print emuPipe $omegaBaryon."\n";
print emuPipe $parameters->{'powerSpectrumIndex'       }->{'value'}."\n";
print emuPipe $parameters->{'sigma_8'                  }->{'value'}."\n";
print emuPipe $parameters->{'darkEnergyEquationOfState'}->{'value'}."\n";
print emuPipe $parameters->{'redshift'                 }->{'value'}."\n";
print emuPipe "2\n";
close(emuPipe);

exit;
