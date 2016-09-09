#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use Galacticus::Options;
use Galacticus::Constraints::DiscrepancyModels;

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
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFile);

# Determine base parameters to use.
my $baseParameters = exists($arguments{'baseParameters'}) ? $arguments{'baseParameters'} : $config->{'likelihood'}->{'baseParameters'};

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($config->{'likelihood'}->{'compilation'},$baseParameters);

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
&Galacticus::Constraints::DiscrepancyModels::RunModels(
	"variableTreeMassResolution"          ,
	"use of variable tree mass resolution",
	$configFile                           ,
	\%arguments                           ,
	$models
    );

exit;
