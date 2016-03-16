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

# Run calculations to determine the model discrepancy arising from the lack of
# baryons in lensing calculations
# Andrew Benson (19-September-2014)

# Get arguments.
die("Usage: lensingBaryonicContribution.pl <configFile> [options]")
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
	label      => "excludingBaryons",
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
	label      => "includingBaryons",
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
		 value => 2.05e-4
	     },
	     {
		 name  => "gravitationalLensingBaryonicModifierBeta",
		  value => 0.62
	     }
	    ]
    }
};

# Run the models.
&DiscrepancyModels::RunModels(
	"lensingBaryonicContribution"                            ,
	"lack of baryonic contribution to lensing cross-sections",
	$configFile                                              ,
	\%arguments                                              ,
	$models
    );

exit;
