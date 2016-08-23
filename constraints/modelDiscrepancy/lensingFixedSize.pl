#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use Galacticus::Options;
use Galacticus::Constraints::DiscrepancyModels;

# Run calculations to determine the model discrepancy arising from the use of the
# fixed size in gravitational lensing calculations.
# Andrew Benson (18-September-2014)

# Get arguments.
die("Usage: lensingFixedSize.pl <configFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = 
    (
     make => "yes",
     plot => "no"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Specify models to run.
my $models = 
{
    default =>
    {
	label      => "defaultSizes",
	parameters =>
	    [
	     # Use default size for lensing.
	     {
		 name  => "analysisMassFunctionsGravitationalLensingSize",
		 value => 1.0e-3
	     }
	    ]
    },
    alternate =>
    {
	label      => "largeSizes",
	parameters => 
	    [
	     # Use large size for lensing.
	     {
		 name  => "analysisMassFunctionsGravitationalLensingSize",
		 value => 2.0e-3
	     }
	    ]
    }
};

# Run the models.
&Galacticus::Constraints::DiscrepancyModels::RunModels(
	"lensingFixedSize"                                  ,
	"use of fixed source size in lensing magnifications",
	$configFile                                         ,
	\%arguments                                         ,
	$models
    );

exit;
