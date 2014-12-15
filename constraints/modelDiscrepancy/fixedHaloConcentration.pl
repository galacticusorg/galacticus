#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use Clone qw(clone);
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
require List::ExtraUtils;
require Galacticus::Constraints::Parameters;
require Galacticus::Constraints::DiscrepancySystematics;
require Galacticus::Launch::PBS;

# Run calculations to determine the model discrepancy arising from the lack of scatter
# in the halo concentration-mass relation.
# Andrew Benson (17-September-2014)

# Get arguments.
die("Usage: fixedHaloConcentration.pl <configFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make => "yes"
    );
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Parse the Galacticus config file.
my $xml    = new XML::Simple;
my $config = $xml->XMLin($configFile, KeyAttr => 0);

# Validate the config file.
die("fixedHaloConcentration.pl: workDirectory must be specified in config file" )
    unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("fixedHaloConcentration.pl: compilation must be specified in config file"   )
    unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("fixedHaloConcentration.pl: baseParameters must be specified in config file")
    unless ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Determine the scratch and work directories.
my $workDirectory    = $config->{'likelihood'}->{'workDirectory'   };
my $scratchDirectory = $config->{'likelihood'}->{'workDirectory'   };
$scratchDirectory    = $config->{'likelihood'}->{'scratchDirectory'}
    if ( exists($config->{'likelihood'}->{'scratchDirectory'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'});

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("fixedHaloConcentration.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = 
    &Parameters::Compilation
    (
     $config->{'likelihood'}->{'compilation'},
     $config->{'likelihood'}->{'baseParameters'}
    );
my @constraints = @{$constraintsRef};

# Switch off thread locking.
$parameters->{'parameter'}->{'treeEvolveThreadLock'}->{'value'} = "false";

# Initialize a stack for PBS models.
my @pbsStack;

# Specify models to run.
my @models = 
    (
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
		  name  => "darkMatterProfileConcentrationDiemerKravtsov2014Scatter",
		  value => 0
	      }
	     ]
     },
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
		  name  => "darkMatterProfileConcentrationDiemerKravtsov2014Scatter",
		  value => 0.16
	      }
	     ]
     }
    );

# Iterate over models.
foreach my $model ( @models ) {
# Specify the output name.
    my $modelDirectory = $workDirectory."/modelDiscrepancy/fixedHaloConcentration/".$model->{'label'};
    system("mkdir -p ".$modelDirectory);
    my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
    push(
	@{$model->{'parameters'}},
	{
	    name  => "galacticusOutputFileName",
	    value => $galacticusFileName
	}
	);
    my $newParameters = clone($parameters);
    $newParameters->{'parameter'}->{$_->{'name'}}->{'value'} = $_->{'value'}
        foreach ( @{$model->{'parameters'}} );
    # Adjust the number of trees to run if specified.
    $newParameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $arguments{'treesPerDecade'}
        if ( exists($arguments{'treesPerDecade'}) );
    # Run the model.
    unless ( -e $galacticusFileName ) {
	# Generate the parameter file.
	my $parameterFileName = $modelDirectory."/parameters.xml";
	&Parameters::Output($newParameters,$parameterFileName);
	# Create PBS job.
	my $command = "mpirun --bynode -np 1 Galacticus.exe ".$parameterFileName."\n";
	foreach my $constraint ( @constraints ) {
	    # Parse the definition file.
	    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
	    # Insert code to run the analysis code.
	    my $analysisCode = $constraintDefinition->{'analysis'};
	    $command .= $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5\n";
	}
	my %job =
	    (
	     launchFile => $modelDirectory."/launch.pbs",
	     label      => "fixedHaloConcentration".$model->{'label'},
	     logFile    => $modelDirectory."/launch.log",
	     command    => $command
	    );
	foreach ( 'ppn', 'walltime', 'memory' ) {
	    $job{$_} = $arguments{$_}
	    if ( exists($arguments{$_}) );
	}
	# Queue the calculation.
	push(
	    @pbsStack,
	    \%job
	    );   
    }
}
# Send jobs to PBS.
&PBS::SubmitJobs(\%arguments,@pbsStack)
    if ( scalar(@pbsStack) > 0 );
# Iterate over constraints.
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Locate the model results.
    my $variableHaloConcentrationResultFileName = 
	$workDirectory.
	"/modelDiscrepancy/fixedHaloConcentration/variableHaloConcentration/".
	$constraintDefinition->{'label'}.".hdf5";
    my $fixedHaloConcentrationResultFileName    = 
	$workDirectory.
	"/modelDiscrepancy/fixedHaloConcentration/fixedHaloConcentration/"   .
	$constraintDefinition->{'label'}.".hdf5";
    # Read the results.
    my $variableHaloConcentrationResult = new PDL::IO::HDF5($variableHaloConcentrationResultFileName);
    my $fixedHaloConcentrationResult    = new PDL::IO::HDF5($fixedHaloConcentrationResultFileName   );
    # Extract the results.
    my $fixedX             = $fixedHaloConcentrationResult   ->dataset('x'         )->get();
    my $fixedY             = $fixedHaloConcentrationResult   ->dataset('y'         )->get();
    my $fixedCovariance    = $variableHaloConcentrationResult->dataset('covariance')->get();
    my $variableY          = $variableHaloConcentrationResult->dataset('y'         )->get();
    my $variableCovariance = $variableHaloConcentrationResult->dataset('covariance')->get();
    # Apply any systematics models.
    my %systematicResults;
    foreach my $argument ( keys(%arguments) ) {
	if ( $argument =~ m/^systematic(.*)/ ) {
	    my $model = $1;
	    if ( exists($DiscrepancySystematics::models{$model}) ) {
		%{$systematicResults{$model}} =
		    &{$DiscrepancySystematics::models{$model}}(
		    \%arguments          ,
		    $constraintDefinition,
		    $fixedX              ,
		    $fixedY              ,
		    $fixedCovariance     ,
		    $variableY           ,
		    $variableCovariance
		);
	    }
	}
    }
    # Find the multiplicative discrepancy between these two models.
    (my $nonZero, my $zero)                      = which_both($fixedY > 0.0);
    my $modelDiscrepancyMultiplicative           = $variableY->copy();
    $modelDiscrepancyMultiplicative->($nonZero) /= $fixedY->($nonZero);
    $modelDiscrepancyMultiplicative->(   $zero) .= 1.0;
    # Compute the covariance.
    my $modelDiscrepancyCovarianceMultiplicative = 
	 $variableCovariance*outer(       1.0/$fixedY   ,       1.0/$fixedY   )
	+$fixedCovariance   *outer($variableY/$fixedY**2,$variableY/$fixedY**2);
    # Output the model discrepancy to file.
    my $outputFile = 
	new PDL::IO::HDF5(
	    ">".
	    $workDirectory.
	    "/modelDiscrepancy/fixedHaloConcentration/discrepancy".
	    ucfirst($constraintDefinition->{'label'}).
	    ".hdf5"
	);
    $outputFile->dataset('multiplicative'          )->set($modelDiscrepancyMultiplicative          );
    $outputFile->dataset('multiplicativeCovariance')->set($modelDiscrepancyCovarianceMultiplicative);
    $outputFile->attrSet(
	description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to use of fixed concentration-mass relation."
	);
    # Add results of systematics models.
    my $systematicGroup = $outputFile->group("systematicModels");
    foreach my $model ( keys(%systematicResults) ) {
	my $modelGroup = $systematicGroup->group($model);
	my %modelResults = %{$systematicResults{$model}};
	foreach my $parameter ( keys(%modelResults) ) {
	    $modelGroup->attrSet($parameter => $modelResults{$parameter});
	}
    }
}

exit;
