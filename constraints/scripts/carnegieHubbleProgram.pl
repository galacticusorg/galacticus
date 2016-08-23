#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use XML::Simple;
use Data::Dumper;
use Galacticus::Options;
use Galacticus::HDF5;

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
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);
  
# Define the constraint.
my $constraintValue            = 74.3;
my $constraintStatisticalError =  1.5;
my $constraintSystematicError  =  2.1;
 
# Read model parameters.
my $galacticus;
$galacticus->{'file' } = $galacticusFile;
$galacticus->{'store'} = 0;
&Galacticus::HDF5::Get_Parameters($galacticus);

# Evaluate the model likelihood.
if ( exists($arguments{'outputFile'}) ) {
    # Compute the likelihood.
    my $constraint;
    $constraint->{'label'                } = "carnegieHubbleProgram";
    $constraint->{'logLikelihood'        } = sclr(-0.5*($galacticus->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant'}->{'value'}-$constraintValue)**2/($constraintSystematicError**2+$constraintStatisticalError**2));
    $constraint->{'logLikelihoodVariance'} = 0.0;
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

# Output the results to file if requested.
if ( exists($arguments{'resultFile'}) ) {
    my $resultsFile     = new PDL::IO::HDF5(">".$arguments{'resultFile'});
    my $x               = pdl ones  (1);
    my $y               = pdl [ $galacticus->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant'}->{'value'}->sclr() ];
    my $yData           = pdl [ $constraintValue ];
    my $covariance      = pdl zeroes(1,1);
    my $covarianceData  = pdl zeroes(1,1);
    $covarianceData    .= $constraintSystematicError**2+$constraintStatisticalError**2;
    $resultsFile->dataset('x'             )->set($x                 );
    $resultsFile->dataset('y'             )->set($y);
    $resultsFile->dataset('covariance'    )->set($covariance);
    $resultsFile->dataset('yData'         )->set($yData);
    $resultsFile->dataset('covarianceData')->set($covarianceData);
}

exit;
