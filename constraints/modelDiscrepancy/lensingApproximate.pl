#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use Galacticus::Options;
use Galacticus::Constraints::DiscrepancyModels;

# Run calculations to determine the model discrepancy arising from the approximate nature of the lensing magnification
# distribution fitting function.
# Andrew Benson (19-September-2014)

# Get arguments.
die("Usage: lensingApproximate.pl <configFile> [options]")
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
	label      => "default",
	parameters =>
	    [
	     # Exlcuding baryons.
	     {
		 name  => "gravitationalLensingMethod",
		 value => "takahashi2011"
	     }
	    ]
    },
    alternate =>
    {
	label      => "alternate",
	parameters => 
	    [
	     # Including baryons.
	     {
		 name  => "gravitationalLensingMethod",
		 value => "baryonicModifier"
	     },
	     {
		 name  => "gravitationalLensingBaryonicModifierOriginalDistribution",
		 value => "takahashi2011"
	     },
	     {
		 name  => "gravitationalLensingBaryonicModifierAlpha",
		 value => 1.0e-3
	     },
	     {
		 name  => "gravitationalLensingBaryonicModifierBeta",
		  value => 3.9
	     }
	    ]
    }
};

# Run the models.
&Galacticus::Constraints::DiscrepancyModels::RunModels(
	"lensingApproximate"                                                               ,
	"the approximate nature of the lensing magnification distribution fitting function",
	$configFile                                                                        ,
	\%arguments                                                                        ,
	$models
    );

exit;
