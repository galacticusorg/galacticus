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
require Galacticus::Constraints::Parameters;
require Galacticus::Constraints::DiscrepancyModels;

# Run calculations to determine the model discrepancy arising from the use of randomized halo spins.
# Andrew Benson (23-March-2014)

# Get arguments.
die("Usage: randomSpins.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = 
    (
     make         => "yes" ,
     plot         => "no"
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Parameters::Parse_Config($configFile);

# Determine base parameters to use.
my $baseParameters = exists($arguments{'baseParameters'}) ? $arguments{'baseParameters'} : $config->{'likelihood'}->{'baseParameters'};

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'likelihood'}->{'compilation'},$baseParameters);

# Extract the standard reset factor for spins.
my $resetFactor = $parameters->{'randomSpinResetMassFactor'}->{'value'};

# Specify models to run.
my $models = 
{
    default =>
    {
	label      => "standard",
	parameters =>
	    [
	     # Use the standard factor for spin reset.
	     {
		 name  => "randomSpinResetMassFactor",
		 value => $resetFactor
	     }
	    ]
    },
    alternate =>
    {
	label      => "short",
	parameters => 
	    [
	     # Use a short factor for spin reset.
	     {
		 name  => "randomSpinResetMassFactor",
		 value => $resetFactor/1.5
	     },
	    ]
    }
};

# Run the models.
&DiscrepancyModels::RunModels(
	"randomSpins"                         ,
	"choice of random spin reset interval",
	$configFile                           ,
	\%arguments                           ,
	$models
    );

exit;
