#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use Galacticus::Options;
use Galacticus::Constraints::DiscrepancyModels;

# Run calculations to determine the model discrepancy arising from the use of the
# instantaneous recycling approximation.
# Andrew Benson (16-September-2014)

# Get arguments.
die("Usage: instantaneousRecycling.pl <configFile> [options]")
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
	label      => "instantaneousRecycling",
	parameters =>
	    [
	     # Switch to using instantaneous recycling.
	     {
		 name  => "stellarPopulationPropertiesMethod",
		 value => "instantaneous"
	     }
	    ]
    },
    alternate =>
    {
	label      => "noninstantaneousRecycling",
	parameters => 
	    [
	     # Switch to using noninstantaneous recycling.
	     {
		 name  => "stellarPopulationPropertiesMethod",
		 value => "noninstantaneous"
	     },
	     {
		 name  => "noninstantHistoryTimesCount",
		 value => 30
	     },
	     {
		 # Cannot use the very simple disks analytic solver in this case as noninstantaneous recycling violates its assumptions.
		 name  => "diskVerySimpleUseAnalyticSolver",
		 value => "false"
	     },
	    ]
    }
};

# Run the models.
&Galacticus::Constraints::DiscrepancyModels::RunModels(
	"instantaneousRecycling"                      ,
	"use of instantaneous recycling approximation",
	$configFile                                   ,
	\%arguments                                   ,
	$models
    );

exit;
