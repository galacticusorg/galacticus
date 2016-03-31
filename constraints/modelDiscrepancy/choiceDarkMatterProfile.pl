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
&Options::Parse_Options(\@ARGV,\%arguments);

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
&DiscrepancyModels::RunModels(
	"choiceDarkMatterProfile"      ,
	"choice of dark matter profile",
	$configFile                    ,
	\%arguments                    ,
	$models
    );

exit;
