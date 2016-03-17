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

# Run calculations to determine the model discrepancy arising from the use of less scatter in
# the merging timescale for satellite galaxies in the Jiang2008 model than is recommended.
# Andrew Benson (20-November-2012)

# Get arguments.
die("Usage: jiang2008MergingTimeScatter.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = 
    (
     make => "yes",
     plot => "no"
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Specify models to run.
my $models = 
{
    default =>
    {
	label      => "defaultScatter",
	parameters =>
	    [
	     # Switch to using Jiang2008 implementation.
	     {
		 name  => "satelliteMergingTimescalesMethod",
		 value => "jiang2008"
	     }
	    ]
    },
    alternate =>
    {
	label      => "recommendedScatter",
	parameters => 
	    [
	     # Switch to using Jiang2008 implementation.
	     {
		 name  => "satelliteMergingTimescalesMethod",
		 value => "jiang2008"
	     },
	     # Switch to using recommended scatter.
	     {
		 name  => "satelliteMergingJiang2008Scatter",
		 value => "0.4"
	     },
	    ]
    }
};

# Run the models.
&DiscrepancyModels::RunModels(
	"jiang2008MergingTimeScatter"    ,
	"use of scatterless merging time",
	$configFile                      ,
	\%arguments                      ,
	$models
    );

exit;
