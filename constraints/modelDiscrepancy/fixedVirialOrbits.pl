#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use Galacticus::Options;
use Galacticus::Constraints::DiscrepancyModels;

# Run calculations to determine the model discrepancy arising from the use of fixed virial orbit
# parameters for satellites.
# Andrew Benson (17-November-2012)

# Get arguments.
die("Usage: fixedVirialOrbits.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
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
	label      => "fixedOrbits",
	parameters =>
	    [
	     # Switch to using fixed virial orbits.
	     {
		 name  => "virialOrbitsMethod",
		 value => "fixed"
	     }
	    ]
    },
    alternate =>
    {
	label      => "variableOrbits",
	parameters => 
	    [
	     # Switch to using variable orbits using the distribution of Jiang et al. (2014).
	     {
		 name  => "virialOrbitMethod",
		 value => "jiang2014"
	     },
	    ]
    }
};

# Run the models.
&Galacticus::Constraints::DiscrepancyModels::RunModels(
	"fixedVirialOrbits"                                     ,
	"use of scatterless virial orbit parameter distribution",
	$configFile                                             ,
	\%arguments                                             ,
	$models
    );

exit;
