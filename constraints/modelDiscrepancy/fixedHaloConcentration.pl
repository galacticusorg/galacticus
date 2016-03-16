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

# Run calculations to determine the model discrepancy arising from the lack of scatter
# in the halo concentration-mass relation.
# Andrew Benson (17-September-2014)

# Get arguments.
die("Usage: fixedHaloConcentration.pl <configFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
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
	label      => "fixedHaloConcentration",
	parameters =>
	    [
	     # Switch to using a fixed concentration relation.
	     {
		 name  => "darkMatterProfileConcentrationMethod",
		 value => "diemerKravtsov2014"
	     },
	     {
		 name  => "darkMatterProfileConcentrationMethod->scatter",
		 value => 0.0
	     }
	    ]
    },
    alternate =>
    {
	label      => "variableHaloConcentration",
	parameters => 
	    [
	     # Switch to using scatter.
	     {
		 name  => "darkMatterProfileConcentrationMethod",
		 value => "diemerKravtsov2014"
	     },
	     {
		 name  => "darkMatterProfileConcentrationMethod->scatter",
		 value => 0.16
	     }
	    ]
    }
};

# Run the models.
&DiscrepancyModels::RunModels(
	"fixedHaloConcentration"                             ,
	"use of scatterless halo concentration-mass relation",
	$configFile                                          ,
	\%arguments                                          ,
	$models
    );

exit;
