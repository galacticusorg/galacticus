#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use XML::Simple;
use Data::Dumper;
require Galacticus::Options;
require Galacticus::HDF5;

# Compute likelihood based on the timing argument constraint on Local Group mass (http://adsabs.harvard.edu/abs/2012ApJ...758...24F).
# Andrew Benson (03-February-2015)

# Get name of input and output files.
die("timingArgument.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Evaluate the model likelihood.
if ( exists($arguments{'outputFile'}) ) {
    # Define the constraint.
    my $constraintValue            = 12.72;
    my $constraintStatisticalError = 0.187;
    # Read model parameters.
    my $galacticus;
    $galacticus->{'file' } = $galacticusFile;
    $galacticus->{'store'} = 0;
    &HDF5::Get_Parameters($galacticus);
    # Find combined mass of M31 and Milky Way.
    my $massLocalGroup =
	log10(
	    +$galacticus->{'parameters'}->{'analysisLGSatelliteMFHaloMassMilkyWay'}->{'value'}
	    +$galacticus->{'parameters'}->{'analysisLGSatelliteMFHaloMassM31'     }->{'value'}
	);
    # Compute the likelihood.
    my $constraint;
    $constraint->{'logLikelihood'} = sclr(-0.5*($massLocalGroup-$constraintValue)**2/+$constraintStatisticalError**2);
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

exit;
