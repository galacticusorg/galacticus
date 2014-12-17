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

# Run calculations to determine the model discrepancy arising from the use of the
# fixed size in gravitational lensing calculations.
# Andrew Benson (18-September-2014)

# Get arguments.
die("Usage: lensingFixedSize.pl <configFile> [options]")
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
die("lensingFixedSize.pl: workDirectory must be specified in config file" )
    unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("lensingFixedSize.pl: compilation must be specified in config file"   )
    unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("lensingFixedSize.pl: baseParameters must be specified in config file")
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
    die("lensingFixedSize.pl: failed to build Galacticus.exe")
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
	 label      => "defaultSizes",
	 parameters =>
	     [
	      # Use default size for lensing.
	      {
		  name  => "analysisMassFunctionsGravitationalLensingSize",
		  value => 1.0e-3
	      }
	     ]
     },
     {
	 label      => "largeSizes",
	 parameters => 
	     [
	      # Use large size for lensing.
	      {
		  name  => "analysisMassFunctionsGravitationalLensingSize",
		  value => 2.0e-3
	      }
	     ]
     }
    );

# Iterate over models.
foreach my $model ( @models ) {
# Specify the output name.
    my $modelDirectory = $workDirectory."/modelDiscrepancy/lensingFixedSize/".$model->{'label'};
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
	     label      => "defaultSizes".$model->{'label'},
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
    my $largeSizesResultFileName = 
	$workDirectory.
	"/modelDiscrepancy/lensingFixedSize/largeSizes/"  .
	$constraintDefinition->{'label'}.".hdf5";
    my $defaultSizesResultFileName    = 
	$workDirectory.
	"/modelDiscrepancy/lensingFixedSize/defaultSizes/".
	$constraintDefinition->{'label'}.".hdf5";
    # Read the results.
    my $largeSizesResult   = new PDL::IO::HDF5($largeSizesResultFileName  );
    my $defaultSizesResult = new PDL::IO::HDF5($defaultSizesResultFileName);
    # Extract the results.
    my $defaultX          = $defaultSizesResult->dataset('x'         )->get();
    my $defaultY          = $defaultSizesResult->dataset('y'         )->get();
    my $defaultCovariance = $largeSizesResult  ->dataset('covariance')->get();
    my $largeY            = $largeSizesResult  ->dataset('y'         )->get();
    my $largeCovariance   = $largeSizesResult  ->dataset('covariance')->get();
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
		    $defaultX            ,
		    $defaultY            ,
		    $defaultCovariance   ,
		    $largeY              ,
		    $largeCovariance
		);
	    }
	}
    }
    # Find the multiplicative discrepancy between these two models.
    (my $nonZero, my $zero)                      = which_both($defaultY > 0.0);
    my $modelDiscrepancyMultiplicative           = $largeY->copy();
    $modelDiscrepancyMultiplicative->($nonZero) /= $defaultY->($nonZero);
    $modelDiscrepancyMultiplicative->(   $zero) .= 1.0;
    # Compute the covariance.
    my $modelDiscrepancyCovarianceMultiplicative = 
	 $largeCovariance  *outer(    1.0/$defaultY   ,    1.0/$defaultY   )
	+$defaultCovariance*outer($largeY/$defaultY**2,$largeY/$defaultY**2);
    # Output the model discrepancy to file.
    my $outputFile = 
	new PDL::IO::HDF5(
	    ">".
	    $workDirectory.
	    "/modelDiscrepancy/lensingFixedSize/discrepancy".
	    ucfirst($constraintDefinition->{'label'}).
	    ".hdf5"
	);
    $outputFile->dataset('multiplicative'          )->set($modelDiscrepancyMultiplicative          );
    $outputFile->dataset('multiplicativeCovariance')->set($modelDiscrepancyCovarianceMultiplicative);
    $outputFile->attrSet(
	description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to use of fixed source sizes in gravitational lensing."
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
