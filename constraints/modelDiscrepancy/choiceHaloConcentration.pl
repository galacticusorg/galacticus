#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use Galacticus::Options;
use Galacticus::Constraints::DiscrepancyModels;

# Run calculations to determine the model discrepancy arising from the choice
# of halo concentration-mass relation.
# Andrew Benson (17-September-2014)

# Get arguments.
die("Usage: choiceHaloConcentration.pl <configFile> [options]")
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
	label      => "defaultHaloConcentration",
	parameters =>
	    [
	     # Use the default concentration method.
	     {
		 name  => "darkMatterProfileConcentrationMethod",
		 value => "diemerKravtsov2014"
	     }
	    ]
    },
    alternate => 
    {
	label      => "alternativeHaloConcentration",
	parameters => 
	    [
	     # Use an alternative concentration method.
	     {
		 name  => "darkMatterProfileConcentrationMethod",
		 value => "duttonMaccio2014"
	     },
	     {
		 name  => "darkMatterProfileConcentrationMethod->fitType",
		 value => "nfwMean200"
	     }
	    ]
    }
};

# Run the models.
&Galacticus::Constraints::DiscrepancyModels::RunModels(
	"choiceHaloConcentration"                   ,
	"choice of halo concentration-mass relation",
	$configFile                                 ,
	\%arguments                                 ,
	$models
    );

exit;
