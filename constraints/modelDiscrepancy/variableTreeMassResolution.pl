#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
require Galacticus::Options;
require Galacticus::Constraints::DiscrepancyModels;

# Run calculations to determine the model discrepancy arising from the use of variable mass resolution when building trees.
# Andrew Benson (28-January-2013)

# Get arguments.
die("Usage: variableTreeMassResolution.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = 
    (
     make => "yes",
     plot => "no"
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Parameters::Parse_Config($configFile);

# Determine base parameters to use.
my $baseParameters = exists($arguments{'baseParameters'}) ? $arguments{'baseParameters'} : $config->{'likelihood'}->{'baseParameters'};

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'likelihood'}->{'compilation'},$baseParameters);

# Specify models to run.
my $models = 
{
    default =>
    {
	label      => "defaultResolution",
	parameters => 
	    [
	     # Use variable mass resolution.
	     {
		 name  => "mergerTreeMassResolutionMethod",
		 value => "scaled"
	     },
	     {
		 name  => "mergerTreeMassResolutionMethod->massResolutionFractional",
		 value => $parameters->{'mergerTreeMassResolutionMethod'}->{'massResolutionFractional'}->{'value'}
	     }
	    ]
    },
    alternate =>
    {
	label      => "higherResolution",
	parameters =>
	    [
	     # Continue to used scale resolution, but increase the fractional resolution.
	     {
		 name  => "mergerTreeMassResolutionMethod",
		 value => "scaled"
	     },
	     {
		 name  => "mergerTreeMassResolutionMethod->massResolutionFractional",
		 value => 1.0e-2*$parameters->{'mergerTreeMassResolutionMethod'}->{'massResolutionFractional'}->{'value'}
	     }
	    ]
    },
};

# Run the models.
&DiscrepancyModels::RunModels(
	"variableTreeMassResolution"          ,
	"use of variable tree mass resolution",
	$configFile                           ,
	\%arguments                           ,
	$models
    );

exit;
