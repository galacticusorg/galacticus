#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
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
require Galacticus::Constraints::Parameters;

# Run calculations to determine the model discrepancy arising from the use of variable mass resolution when building trees.
# Andrew Benson (28-January-2013)

# Get arguments.
die("Usage: variableTreeMassResolution.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
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

# Parse the constraint config file.
my $xml    = new XML::Simple;
my $config = $xml->XMLin($configFile, KeyAttr => 0);

# Validate the config file.
die("variableTreeMassResolution.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'workDirectory' }) );
die("variableTreeMassResolution.pl: compilation must be specified in config file"   ) unless ( exists($config->{'compilation'   }) );
die("variableTreeMassResolution.pl: baseParameters must be specified in config file") unless ( exists($config->{'baseParameters'}) );

# Determine the scratch and work directories.
my $workDirectory    = $config->{'workDirectory'};
my $scratchDirectory = $config->{'workDirectory'};
$scratchDirectory    = $config->{'scratchDirectory'} if ( exists($config->{'scratchDirectory'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'workDirectory'});

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("variableTreeMassResolution.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'compilation'},$config->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Initialize a stack for PBS models.
my @pbsStack;

# Specify models to run.
my @models = 
    (
     {
	 label      => "fixedResolution",
	 parameters =>
	     [
	      # Switch to using fixed mass resolution.
	      {
		  name  => "mergerTreesBuildMassResolutionMethod",
		  value => "fixed"
	      }
	     ]
     },
     {
	 label      => "variableResolution",
	 parameters => 
	     [
	      # Switch to using variable mass resolution.
	      {
		  name  => "mergerTreesBuildMassResolutionMethod",
		  value => "scaled"
	      },
	     ]
     }
    );

# Iterate over models.
foreach my $model ( @models ) {
# Specify the output name.
    my $modelDirectory = $workDirectory."/modelDiscrepancy/variableTreeMassResolution/".$model->{'label'};
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
    # Set the fixed mass resolution.
    $newParameters->{'parameter'}->{'mergerTreeBuildMassResolutionScaledMinimum'}->{'value'} =
	$newParameters->{'parameter'}->{'mergerTreeBuildMassResolutionScaledMinimum'}->{'value'};
    $newParameters->{'parameter'}->{'mergerTreeBuildMassResolutionScaledMinimum'}->{'value'} =
	$arguments{'massResolutionFixed'}
           if ( exists($arguments{'massResolutionFixed'}) );
    # Run the model.
    unless ( -e $galacticusFileName ) {
	# Generate the parameter file.
	my $parameterFileName = $modelDirectory."/parameters.xml";
	&Parameters::Output($newParameters,$parameterFileName);
	# Create a batch script for PBS.
	my $batchScriptFileName = $modelDirectory."/launch.pbs";
	open(oHndl,">".$batchScriptFileName);
	print oHndl "#!/bin/bash\n";
	print oHndl "#PBS -N variableTreeMassResolution".$model->{'label'}."\n";
	print oHndl "#PBS -l walltime=3:00:00\n";
	print oHndl "#PBS -l mem=4gb\n";
	print oHndl "#PBS -l nodes=1:ppn=12\n";
	print oHndl "#PBS -j oe\n";
	print oHndl "#PBS -o ".$modelDirectory."/launch.log\n";
	print oHndl "#PBS -V\n";
	print oHndl "cd \$PBS_O_WORKDIR\n";
	print oHndl "export LD_LIBRARY_PATH=/home/abenson/Galacticus/Tools/lib:/home/abenson/Galacticus/Tools/lib64:\$LD_LIBRARY_PATH\n";
	print oHndl "export PATH=/home/abenson/Galacticus/Tools/bin:\$PATH\n";
	print oHndl "export GFORTRAN_ERROR_DUMPCORE=YES\n";
	print oHndl "ulimit -t unlimited\n";
	print oHndl "ulimit -c unlimited\n";
	print oHndl "export OMP_NUM_THREADS=12\n";
	print oHndl "mpirun --bynode -np 1 Galacticus.exe ".$parameterFileName."\n";
	foreach my $constraint ( @constraints ) {
	    # Parse the definition file.
	    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
	    # Insert code to run the analysis code.
	    my $analysisCode = $constraintDefinition->{'analysis'};
	    print oHndl $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".xml\n";
	}
	close(oHndl);
	# Queue the calculation.
	push(
	    @pbsStack,
	    $batchScriptFileName
	    );   
    }
}
# Send jobs to PBS.
&PBS_Submit(@pbsStack)
    if ( scalar(@pbsStack) > 0 );
# Iterate over constraints.
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Locate the model results.
    my $variableResolutionResultFileName = $workDirectory."/modelDiscrepancy/variableTreeMassResolution/variableResolution/".$constraintDefinition->{'label'}.".xml";
    my $fixedResolutionResultFileName    = $workDirectory."/modelDiscrepancy/variableTreeMassResolution/fixedResolution/"   .$constraintDefinition->{'label'}.".xml";
    # Read the results.
    my $variableResolutionResult = $xml->XMLin($variableResolutionResultFileName);
    my $fixedResolutionResult    = $xml->XMLin($fixedResolutionResultFileName   );
    # Extract the results.
    my $fixedX             = pdl @{$fixedResolutionResult   ->{'x'         }};
    my $fixedY             = pdl @{$fixedResolutionResult   ->{'y'         }};
    my $fixedCovariance    = pdl @{$variableResolutionResult->{'covariance'}};
    my $variableY          = pdl @{$variableResolutionResult->{'y'         }};
    my $variableCovariance = pdl @{$variableResolutionResult->{'covariance'}};
    my $ySize              = nelem($fixedY);
    $fixedCovariance       = reshape($fixedCovariance   ,$ySize,$ySize);
    $variableCovariance    = reshape($variableCovariance,$ySize,$ySize);
    # Find the multiplicative discrepancy between these two models.
    (my $nonZero, my $zero)            = which_both($fixedY > 0.0);
    my $modelDiscrepancyMultiplicative = $variableY;
    $modelDiscrepancyMultiplicative->($nonZero) /= $fixedY->($nonZero);
    $modelDiscrepancyMultiplicative->($zero   ) .= 1.0;
    # Compute the covariance arising from finite number of trees in the variable resolution model.
    my $modelDiscrepancyCovariance                      = $variableCovariance+$fixedCovariance;
    my $fixedYOuter                                     = outer($fixedY,$fixedY);
    $fixedYOuter($modelDiscrepancyCovariance == 0.0;?) .= 1.0;
    my $modelDiscrepancyCovarianceMultiplicative        = $modelDiscrepancyCovariance->copy();
    $modelDiscrepancyCovarianceMultiplicative          .=
	log
	(
	 1.0
	 +$modelDiscrepancyCovariance
	 /$fixedYOuter
	);
    # Output the model discrepancy to file.
    my $outputFile = new PDL::IO::HDF5(">".$workDirectory."/modelDiscrepancy/variableTreeMassResolution/discrepancy".ucfirst($constraintDefinition->{'label'}).".hdf5");
    $outputFile->dataset('multiplicative'          )->set($modelDiscrepancyMultiplicative          );
    $outputFile->dataset('multiplicativeCovariance')->set($modelDiscrepancyCovarianceMultiplicative);
    $outputFile->attrSet(
	description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to use of variable mass resolution in merger trees."
	);
}

exit;

sub PBS_Submit {
    # Submit jobs to PBS and wait for them to finish.
    my @pbsStack = @_;
    my %pbsJobs;
    # Determine maximum number allowed in queue at once.
    my $jobMaximum = 10;
    $jobMaximum = $arguments{'pbsJobMaximum'}
    if ( exists($arguments{'pbsJobMaximum'}) );
    # Submit jobs and wait.
    print "Waiting for PBS jobs to finish...\n";
    while ( scalar(keys %pbsJobs) > 0 || scalar(@pbsStack) > 0 ) {
	# Find all PBS jobs that are running.
	my %runningPBSJobs;
	undef(%runningPBSJobs);
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^Job\sId:\s+(\S+)/ ) {$runningPBSJobs{$1} = 1};
	}
	close(pHndl);
	foreach my $jobID ( keys(%pbsJobs) ) {
	    unless ( exists($runningPBSJobs{$jobID}) ) {
		print "PBS job ".$jobID." has finished.\n";
		# Remove the job ID from the list of active PBS jobs.
		delete($pbsJobs{$jobID});
	    }
	}
	# If fewer than ten jobs are in the queue, pop one off the stack.
	if ( scalar(@pbsStack) > 0 && scalar(keys %pbsJobs) < 20 ) {
	    my $batchScript = pop(@pbsStack);
	    # Submit the PBS job.
	    open(pHndl,"qsub ".$batchScript."|");
	    my $jobID = "";
	    while ( my $line = <pHndl> ) {
	    	if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
	    }
	    close(pHndl);	    
	    # Add the job number to the active job hash.
	    unless ( $jobID eq "" ) {
	    	$pbsJobs{$jobID} = 1;
	    }
	    sleep 5;
	} else {
	    # Wait.
	    sleep 60;
	}
    }
}
