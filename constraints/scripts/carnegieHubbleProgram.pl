#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use Astro::Cosmology;
use XML::Simple;
use Data::Dumper;
require Galacticus::Options;
require Galacticus::HDF5;
require Galacticus::Constraints::Covariances;
require Stats::Histograms;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Compute likelihood based on the Carnegie Hubble Program constraint on H_0 (http://adsabs.harvard.edu/abs/2012ApJ...758...24F).
# Andrew Benson (16-October-2014)

# Get name of input and output files.
die("carnegieHubbleProgram.pl <galacticusFile> [options]")
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
    my $constraintValue            = 74.3;
    my $constraintStatisticalError =  1.5;
    my $constraintSystematicError  =  2.1;
    # Read model parameters.
    my $galacticus;
    $galacticus->{'file' } = $galacticusFile;
    $galacticus->{'store'} = 0;
    &HDF5::Get_Parameters($galacticus);
    # Compute the likelihood.
    my $constraint;
    $constraint->{'logLikelihood'} = -0.5*($galacticus->{'parameters'}->{'H_0'}-$constraintValue)**2/($constraintSystematicError**2+$constraintStatisticalError**2);
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

exit;
