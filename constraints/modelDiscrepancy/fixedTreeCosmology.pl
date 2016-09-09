#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use Clone qw(clone);
use Galacticus::Options;
use Galacticus::Constraints::DiscrepancyModels;
use Galacticus::Launch::PBS;
use Galacticus::Constraints::Parameters;

# Run calculations to determine the model discrepancy arising from the use of a cosmology-independent set of merger trees.
# Andrew Benson (09-January-2014)

# Get arguments.
die("Usage: fixedTreeCosmology.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = 
    (
     make => "yes",
     plot => "no"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $xml    = new XML::Simple;
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFile);

# Determine base parameters to use.
my $baseParameters = exists($arguments{'baseParameters'}) 
	             ?
	                    $arguments{'baseParameters'}
                     :
	                    $config->{'likelihood'}->{'baseParameters'};

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($config->{'likelihood'}->{'compilation'},$baseParameters);
my @constraints = @{$constraintsRef};

# Determine the work directory.
my $workDirectory = $config->{'likelihood'}->{'workDirectory'};

# Specify models to run.
my $modelDirectory = $workDirectory."/modelDiscrepancy/fixedTreeCosmology/trees/";
my $treeFile       = $modelDirectory."trees.hdf5";
my @treeConstructParameters =
    (
     # Define a model that will build fixed cosmology trees.
     {
	 name  => "mergerTreeConstructMethod",
	 value => "build"
     },
     {
	 name  => "mergerTreesWrite",
	 value => "true"
     },
     {
	 name  => "mergerTreeOperatorMethod",
	 value => "export"
     },
     {
	 name  => "mergerTreeOperatorMethod->outputFileName",
	 value => $treeFile
     },
     {
	 name  => "mergerTreeOperatorMethod->exportFormat",
	 value => "galacticus"
     },
     {
	 name  => "cosmologyParametersMethod->HubbleConstant",
	 value => 67.8100148730579
     },
     {
	 name  => "powerSpectrumPrimordialMethod->index",
	 value => 0.96763953040293
     },
     {
	 name  => "cosmologyParametersMethod->OmegaBaryon",
	 value => 0.048404634
     },
     {
	 name  => "cosmologyParametersMethod->OmegaMatter",
	 value => 0.30766318
     },
     {
	 name  => "cosmologicalMassVarianceMethod->sigma_8",
	 value => 0.814932725148418
     },
     {
	 name  => "reionizationSuppressionOpticalDepth",
	 value => 0.0658422364961477
     },
     {
	 name  => "cosmologyParametersMethod->OmegaDarkEnergy",
	 value => 0.69233682
     },
     {
	 name  => "accretionHaloMethod",
	 value => "null"
     },
     {
	 name  => "mergerTreePruneBaryons",
	 value => "false"
     },
     {
	 name  => "mergerTreeAnalyses",
	 value => undef()
     },
     {
	 name  => "outputRedshifts",
	 value => "9.0"
     }
    );
# Construct the model.
my $newParameters = clone($parameters);
foreach my $parameter ( @treeConstructParameters ) {
    my $thisParameter = $newParameters;
    foreach ( split(/\-\>/,$parameter->{'name'}) ) {
	$thisParameter->{$_}->{'value'} = undef()
	    unless ( exists($thisParameter->{$_}) );
	$thisParameter = $thisParameter->{$_};
    }
    $thisParameter->{'value'} = $parameter->{'value'}
}
# Remove any analyses.
delete($newParameters->{'mergerTreeAnalyses'})
    if ( exists($newParameters->{'mergerTreeAnalyses'}) );
# Set output file name.
$newParameters->{'galacticusOutputFileName'}->{'value'} = $modelDirectory."galacticus.hdf5";
# Adjust the number of trees to run if specified.
$newParameters->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $arguments{'treesPerDecade'}
   if ( exists($arguments{'treesPerDecade'}) );
# Generate the parameter file.
unless ( -e $treeFile ) {
    system("mkdir -p ".$modelDirectory);
    my $parameterFileName = $modelDirectory."/parameters.xml";
    &Galacticus::Constraints::Parameters::Output($newParameters,$parameterFileName);
    # Create PBS job.
    my $command = "mpirun --bynode -np 1 Galacticus.exe ".$parameterFileName."\n";
    my %job =
	(
	 launchFile => $modelDirectory."/launch.pbs",
	 label      => "fixedTreeCosmologyConstruct",
	 logFile    => $modelDirectory."/launch.log",
	 command    => $command
	);
    foreach ( 'ppn', 'walltime', 'memory' ) {
	$job{$_} = $arguments{$_}
	   if ( exists($arguments{$_}) );
    }
    # Queue the calculation.
    my @pbsStack = ( \%job );
    # Send jobs to PBS.
    &Galacticus::Launch::PBS::SubmitJobs(\%arguments,@pbsStack);
}

# Specify models to run.
my $models = 
{
    default =>
    {
	label      => "fixedCosmology",
	parameters => 
	    [
	     # Read merger trees from file.
	     {
		 name  => "mergerTreeConstructMethod",
		 value => "read"
	     },
	     {
		 name  => "mergerTreeReadFileName",
		 value => $treeFile
	     },
	    ]
    },
    alternate =>
    {
	label      => "consistentCosmology",
	parameters =>
	    [
	     # Trees will be generated automatically in this case.
	    ]
    },
};

# Run the models.
&Galacticus::Constraints::DiscrepancyModels::RunModels(
	"fixedTreeCosmology"              ,
	"use of fixed cosmology for trees",
	$configFile                       ,
	\%arguments                       ,
	$models
    );

exit;
