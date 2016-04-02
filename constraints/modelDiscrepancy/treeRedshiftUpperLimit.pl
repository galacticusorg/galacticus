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

# Run calculations to determine the model discrepancy arising from the use of an upper limit to
# redshifts in merger trees.
# Andrew Benson (12-May-2015)

# Get arguments.
die("Usage: treeRedshiftUpperLimit.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = 
    (
     make         => "yes",
     plot         => "no" 
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Specify models to run.
my $models = 
{
    default =>
    {
	label      => "upperLimit",
	parameters =>
	    [
	     # Use whatever upper limit to tree redshifts is set by default.
	    ]
    },
    alternate =>
    {
	label      => "noUpperLimit",
	parameters => 
	    [
	     # Switch to using no upper limit to tree redshifts.
	     {
		 name  => "mergerTreeBuilderMethod->redshiftMaximum",
		 value => "1.0e5"
	     },
	    ]
    }
};

# Run the models.
&DiscrepancyModels::RunModels(
	"treeRedshiftUpperLimit"                           ,
	"use of an upper limit to redshift in merger trees",
	$configFile                                        ,
	\%arguments                                        ,
	$models
    );

exit;
