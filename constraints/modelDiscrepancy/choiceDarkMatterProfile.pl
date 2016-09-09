#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use Galacticus::Options;
use Galacticus::Constraints::DiscrepancyModels;

# Run calculations to determine the model discrepancy arising from the choice
# of dark matter density profile.
# Andrew Benson (17-September-2014)

# Get arguments.
die("Usage: choiceDarkMatterProfile.pl <configFile> [options]")
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
	label      => "defaultDarkMatterProfile",
	parameters =>
	    [
	     # Use the default dark matter profile method.
	     {
		 name  => "darkMatterProfileMethod",
		 value => "NFW"
	     }
	    ]
    },
    alternate =>
    {
	label      => "alternativeDarkMatterProfile",
	parameters => 
	    [
	     # Use an alternative dark matter profile method.
	     {
		 name  => "darkMatterProfileMethod",
		 value => "einasto"
	     },
	     {
		 name  => "treeNodeMethodDarkMatterProfile",
		 value => "scaleShape"
	     }
	    ]
    }
};

# Run the models.
&Galacticus::Constraints::DiscrepancyModels::RunModels(
	"choiceDarkMatterProfile"      ,
	"choice of dark matter profile",
	$configFile                    ,
	\%arguments                    ,
	$models
    );

exit;
