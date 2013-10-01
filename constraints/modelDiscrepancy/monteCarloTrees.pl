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

# Run calculations to determine the model discrepancy arising from the use of Monte Carlo merger trees.
# Andrew Benson (16-November-2012)

# Get arguments.
die("Usage: monteCarloTrees.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
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
die("monteCarloTrees.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'workDirectory' }) );
die("monteCarloTrees.pl: compilation must be specified in config file"   ) unless ( exists($config->{'compilation'   }) );
die("monteCarloTrees.pl: baseParameters must be specified in config file") unless ( exists($config->{'baseParameters'}) );

# Determine the scratch and work directories.
my $workDirectory    = $config->{'workDirectory'};
my $scratchDirectory = $config->{'workDirectory'};
$scratchDirectory    = $config->{'scratchDirectory'} if ( exists($config->{'scratchDirectory'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'workDirectory'});

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("monteCarloTrees.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get all required merger trees from the Millennium Simulation database.
my $millenniumDbCommand = $galacticusPath."/scripts/aux/Millennium_Trees_Grab.pl --table MPAHaloTrees..MHalo --traceParticles no";
# Get Millennium database username and password if available.
$millenniumDbCommand .= " --user "    .$arguments{'user'    } 
    if ( exists($arguments{'user'    }) );
$millenniumDbCommand .= " --password ".$arguments{'password'} 
    if ( exists($arguments{'password'}) );
# Determine the number of subvolumes to obtain.
my $subVolumeCount    = 32;
$subVolumeCount = $arguments{'subVolumeCount'}
    if ( exists($arguments{'subVolumeCount'}) );
die('monteCarloTrees.pl: subVolumeCount must be between 4 and 512 inclusive')
    if ($subVolumeCount < 4 || $subVolumeCount > 512);
die('monteCarloTrees.pl: subVolumeCount must be a power of 2')
    unless(($subVolumeCount & ($subVolumeCount-1)) == 0);

# Determine the path to which merger trees should be stored.
my $treeDirectory     = $workDirectory;
if ( -e $galacticusPath."/galacticusConfig.xml" ) {
    my $xml    = new XML::Simple;
    my $config = $xml->XMLin($galacticusPath."/galacticusConfig.xml");
    if ( exists($config->{'millenniumDB'}->{'host'}) ) {
	foreach ( keys(%{$config->{'millenniumDB'}->{'host'}}) ) {
	    $treeDirectory = $config->{'millenniumDB'}->{'host'}->{$_}->{'treePath'}
	        if ( ( $ENV{'HOSTNAME'} =~ m/$_/ || $_ eq "default" ) && exists($config->{'millenniumDB'}->{'host'}->{$_}->{'treePath'}) );
	}
    }
}

# Retrieve merger trees from the Millennium database.
system("cd ".$galacticusPath."; make Millennium_Merger_Tree_File_Maker.exe");
my $failuresEncountered     = 0;
my $subVolumeRetrievedCount = -1;
while ( $subVolumeRetrievedCount < $subVolumeCount-1 ) {
    ++$subVolumeRetrievedCount;
    # Check if we already have this file.
    unless ( -e $treeDirectory."/subvolume".$subVolumeRetrievedCount.".hdf5" ) {
	# Retrieve the raw data.
	my $thisMillenniumDbCommand = $millenniumDbCommand." --select \"root.fileNr = ".$subVolumeRetrievedCount." and root.snapNum = 63\" --output ".$treeDirectory."/subvolume".$subVolumeRetrievedCount.".csv";
	system($thisMillenniumDbCommand);
	unless ( $? == 0 ) {
	    print "monteCarloTrees.pl: failed to retrieve subvolume ".$subVolumeRetrievedCount." from Millennium database\n";
	    ++$failuresEncountered;
	}
	# Convert to Galacticus format.
	system($galacticusPath."/Millennium_Merger_Tree_File_Maker.exe ".$treeDirectory."/subvolume".$subVolumeRetrievedCount.".csv none ".$treeDirectory."/subvolume".$subVolumeRetrievedCount.".hdf5 galacticus 1");
	unless ( $? == 0 ) {
	    print "monteCarloTrees.pl: failed to convert subvolume ".$subVolumeRetrievedCount." to Galacticus format\n";
	    ++$failuresEncountered;
	} else {
	    unlink($treeDirectory."/subvolume".$subVolumeRetrievedCount.".csv");
	}
    }
}
die("monteCarloTrees.pl: failures encountered when downloading from Millennium database")
    unless ( $failuresEncountered == 0 );

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'compilation'},$config->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Initialize a stack for PBS models.
my @pbsStack;

# Specify models to run.
my @models;
my @monteCarloParameters =
    (
     # Switch to using Monte-Carlo trees.
     {
	 name  => "mergerTreeConstructMethod",
	 value => "build"
     },
     # Match the mass resolution and timestepping of the N-body merger trees.
     {
	 name  => "mergerTreeBuildMassResolutionMethod",
	 value => "fixed"
     },
     {
	 name  => "mergerTreeBuildMassResolutionFixed",
	 value => 2.3578789e10
     },
     {
	 name  => "mergerTreeRegridTimes",
	 value => "true"
     },
     {
	 name  => "mergerTreeRegridCount",
	 value => 60
     },
     {
	 name  => "mergerTreeRegridSpacing",
	 value => "millennium"
     },
     # Set cosmological parameters to match the Millennium Simulation.
     {
	 name  => "H_0",
	 value => 73.0
     },
     {
	 name  => "Omega_Matter",
	 value => 0.25
     },
     {
	 name  => "Omega_DE",
	 value => 0.75
     },
     {
	 name  => "Omega_b",
	 value => 0.0455
     },
     {
	 name  => "sigma_8",
	 value => 0.9
     },
     {
	 name  => "powerSpectrumIndex",
	 value => 1.0
     },
     {
	 name  => "powerSpectrumReferenceWavenumber",
	 value => 1.0
     },
     {
	 name  => "powerSpectrumRunning",
	 value => 0.0
     },
     {
	 name  => "transferFunctionMethod",
	 value => "file"
     },
     {
	 name  => "transferFunctionFile",
	 value => "data/largeScaleStructure/transferFunctionMillenniumSimulation.xml"
     },
     # Set an initial random number seed.
     {
	 name  => "randomSeed",
	 value => 826
     }
    );
# Adjust the number of trees to run if specified.
push(
    @monteCarloParameters,
    {
	name  => "mergerTreeBuildTreesPerDecade",
	value => $arguments{'treesPerDecade'}
    }
    )
    if ( exists($arguments{'treesPerDecade'}) );
my %monteCarloModel =
    (
     label      => "monteCarlo",
     parameters => \@monteCarloParameters
    );
push(@models,\%monteCarloModel);
my $randomSeed = 826;
for(my $iSubvolume=0;$iSubvolume<$subVolumeCount;++$iSubvolume) {
    ++$randomSeed;
    push(
	@models,
	{
	    label      => "nBody".$iSubvolume,
	    parameters => 
		[
		 # Switch to using N-body trees.
		 {
		  name  => "mergerTreeConstructMethod",
		  value => "read"
		 },
		 {
		     name  => "mergerTreeReadFileName",
		     value => $treeDirectory."/subvolume".$iSubvolume.".hdf5"
		 },
		 # Prune branches of the tree to keep only those with at least 20 particles.
		 {
		     name  => "mergerTreePruneBranches",
		     value => "true"
		 },
		 {
		     name  => "mergerTreePruningMassThreshold",
		     value => 2.3578789e10
		 },
		 # Set merger tree reading options.
		 {
		     name  => "mergerTreeReadTreeIndexToRootNodeIndex",
		     value => "true"
		 },
		 {
		     name  => "allTreesExistAtFinalTime",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadOutputTimeSnapTolerance",
		     value => 1.0e-3
		 },
		 {
		     name  => "mergerTreeReadPresetPositions",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetOrbits",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetMergerTimes",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetMergerNodes",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetSubhaloMasses",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetSpins",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadPresetScaleRadii",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadAllowBranchJumps",
		     value => "false"
		 },
		 {
		     name  => "mergerTreeReadAllowSubhaloPromotions",
		     value => "false"
		 },
		 # Set cosmological parameters.
		 {
		     name  => "H_0",
		     value => 73.0
		 },
		 {
		     name  => "Omega_Matter",
		     value => 0.25
		 },
		 {
		     name  => "Omega_DE",
		     value => 0.75
		 },
		 {
		     name  => "Omega_b",
		     value => 0.0455
		 },
		 {
		     name  => "sigma_8",
		     value => 0.9
		 },
		 {
		     name  => "powerSpectrumIndex",
		     value => 1.0
		 },
		 {
		     name  => "powerSpectrumReferenceWavenumber",
		     value => 1.0
		 },
		 {
		     name  => "powerSpectrumRunning",
		     value => 0.0
		 },
		 # Set an initial random number seed.
		 {
		     name  => "randomSeed",
		     value => $randomSeed
		 }
		]
	}
	);
}

# Iterate over models.
foreach my $model ( @models ) {
# Specify the output name.
    my $modelDirectory = $workDirectory."/modelDiscrepancy/monteCarloTrees/".$model->{'label'};
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
    # Run the model.
    unless ( -e $galacticusFileName ) {
	# Generate the parameter file.
	my $parameterFileName = $modelDirectory."/parameters.xml";
	&Parameters::Output($newParameters,$parameterFileName);
	# Create a batch script for PBS.
	my $batchScriptFileName = $modelDirectory."/launch.pbs";
	open(oHndl,">".$batchScriptFileName);
	print oHndl "#!/bin/bash\n";
	print oHndl "#PBS -N monteCarloTrees".$model->{'label'}."\n";
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
    # Locate the Monte Carlo model results.
    my $monteCarloResultFileName = $workDirectory."/modelDiscrepancy/monteCarloTrees/monteCarlo/".$constraintDefinition->{'label'}.".xml";
    # Read the results.
    my $monteCarloResult         = $xml->XMLin($monteCarloResultFileName);
    # Extract the results.
    my $monteCarloMass           = pdl @{$monteCarloResult->{'x'         }};
    my $monteCarloMassFunction   = pdl @{$monteCarloResult->{'y'         }};
    my $monteCarloCovariance     = pdl @{$monteCarloResult->{'covariance'}};
    my $ySize                    = nelem($monteCarloMassFunction);
    $monteCarloCovariance        = reshape($monteCarloCovariance,$ySize,$ySize);
    # Iterate over N-body models.
    my $nBodyMass;
    my @nBodyMassFunctions;
    for(my $iSubvolume=0;$iSubvolume<$subVolumeCount;++$iSubvolume) {
	# Locate the N-body model results.
	my $nbodyResultFileName = $workDirectory."/modelDiscrepancy/monteCarloTrees/nBody".$iSubvolume."/".$constraintDefinition->{'label'}.".xml";
	# Read the results.
	my $nbodyResult         = $xml->XMLin($nbodyResultFileName);
	# Extract the results.
	$nBodyMass              = pdl @{$nbodyResult->{'x'}} ;
	push(@nBodyMassFunctions, pdl @{$nbodyResult->{'y'}});
    }
    # Construct samples
    my @covariances;
    my $combineCount = 1;
    my $nBodyMassFunction;
    while ( $combineCount <= $subVolumeCount ) {
	# Combine mass functions.
	my @combinedMassFunctions;
	my $subVolumeBegin = 0;
	while ( $subVolumeBegin+$combineCount-1 < $subVolumeCount ) {
	    my $massFunction = pdl zeroes(nelem($nBodyMass));
	    for(my $subVolume=$subVolumeBegin;$subVolume<$subVolumeBegin+$combineCount;++$subVolume) {
		$massFunction += $nBodyMassFunctions[$subVolume];
	    }
	    $massFunction      *= 512.0/$combineCount;
	    push(@combinedMassFunctions,$massFunction);
	    $nBodyMassFunction  = $massFunction
		if ( $combineCount == $subVolumeCount );
	    $subVolumeBegin    += $combineCount
	}
	# Measure the covariance.
	if ( $combineCount < $subVolumeCount ) {
	    my $covariance = pdl zeroes(nelem($nBodyMass),nelem($nBodyMass));
	    for(my $i=0;$i<nelem($nBodyMass);++$i) {
		for(my $j=0;$j<nelem($nBodyMass);++$j) {
		    my $iMean = pdl 0.0;
		    my $jMean = pdl 0.0;
		    foreach my $combinedMassFunction ( @combinedMassFunctions ) {
			$iMean                   += $combinedMassFunction->($i)                            ;
			$jMean                   +=                             $combinedMassFunction->($j);
			$covariance->(($i),($j)) += $combinedMassFunction->($i)*$combinedMassFunction->($j);
		    }
		    $covariance->(($i),($j)) .= ($covariance->(($i),($j))-$iMean*$jMean/scalar(@combinedMassFunctions))/(scalar(@combinedMassFunctions)-1);
		}
	    }
	    push(@covariances,$covariance);
	}
	# Double the number of subvolumes to combined.
	$combineCount *= 2;
    }
    # Find the multiplicative between these two models.
    (my $nonZero, my $zero)                      = which_both($monteCarloMassFunction > 0.0);
    my $modelDiscrepancyMultiplicative           = $nBodyMassFunction->copy();
    $modelDiscrepancyMultiplicative->($nonZero) /= $monteCarloMassFunction->($nonZero);
    $modelDiscrepancyMultiplicative->($zero   ) .= 1.0;
    # Limit the multiplicative correction to 1.5 - this is to avoid massive corrections that we think arise from differences in
    # merging, but which need to be understood more carefully.
    my $tooHigh = which($modelDiscrepancyMultiplicative > 1.5);
    $modelDiscrepancyMultiplicative->($tooHigh) .= 1.5;
    # Find the covariance remaining in the complete set of subvolumes.
    my $extrapolateFrom                          = scalar(@covariances)-3;
    $extrapolateFrom                             = 0
	if ( $extrapolateFrom < 0 );
    my $modelDiscrepancyCovariance               = $covariances[$extrapolateFrom]/(2**(scalar(@covariances)-$extrapolateFrom));
    $modelDiscrepancyCovariance                 += $monteCarloCovariance;
    my $monteCarloMassFunctionOuter              = outer($monteCarloMassFunction,$monteCarloMassFunction);
    $monteCarloMassFunctionOuter($modelDiscrepancyCovariance == 0.0;?) .= 1.0;
    my $modelDiscrepancyCovarianceMultiplicative = $modelDiscrepancyCovariance->copy();
    $modelDiscrepancyCovarianceMultiplicative .=
	log
	(
	 1.0
	 +$modelDiscrepancyCovariance
	 /$monteCarloMassFunctionOuter
	);
    # Output the model discrepancy to file.
    my $outputFile = new PDL::IO::HDF5(">".$workDirectory."/modelDiscrepancy/monteCarloTrees/discrepancy".ucfirst($constraintDefinition->{'label'}).".hdf5");
    $outputFile->dataset('multiplicative'          )->set($modelDiscrepancyMultiplicative          );
    $outputFile->dataset('multiplicativeCovariance')->set($modelDiscrepancyCovarianceMultiplicative);
    $outputFile->attrSet(
     	description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to use of Monte Carlo merger trees."
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
