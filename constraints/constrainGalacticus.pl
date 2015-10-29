#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use Config;
use if $Config{'useithreads'}, 'threads';
use XML::Simple;
use Data::Dumper;
use Fcntl qw(:DEFAULT :flock);
use MIME::Lite;
use List::Util;
use Storable;
use Fcntl;
use POSIX::RT::Semaphore;
use PDL;
use Storable qw(dclone);require Fortran::Utils;
require File::Which;
require File::NFSLock;
require System::Redirect;
require Galacticus::Constraints::Parameters;
require List::ExtraUtils;
my $useThreads = 0;
if ( $^V gt v5.10.0 && $Config{useithreads} ) {
    $useThreads = 1;
    &PDL::no_clone_skip_warning;
}
# Script timing
# use Time::HiRes qw(gettimeofday tv_interval);
# my $self      = $0;
# my $timeStart = [gettimeofday];

# Finds constraints on Galacticus parameters.
# Andrew Benson (02-September-2011)

# Get command-line parameters.
die("Usage: constrainGalacticus.pl <configFile> <mpiRank> <likelihoodFile> <temperature> <store> <param1> [<param2>......]") 
    unless ( scalar(@ARGV) > 5 );
my $configFile = $ARGV[0];
my $config     = &Parameters::Parse_Config($configFile,useStored => 1);

# Error codes.
my $errorStatusXCPU = 1025; # CPU time limit exceeded.

# Obtain a lock if we are to run models sequentially.
my $runSequential = "no";
$runSequential = $config->{'likelihood'}->{'sequentialModels'}
    if ( exists($config->{'likelihood'}->{'sequentialModels'}) );
my $modelLockFile = "/dev/shm/glcLock";
my $modelLock = new File::NFSLock {
    file               => $modelLockFile,
    lock_type          => LOCK_EX
}
if ( $runSequential eq "yes" );

# Get the MPI rank.
my $mpiRank        = $ARGV[1];

# Get the name for the likelihood file.
my $likelihoodFile = $ARGV[2];

# Get the temperature.
my $temperature    = $ARGV[3];

# Get the store status.
my $store          = $ARGV[4];
# Get a hash of the new parameters.
my $newParameters = &Parameters::Convert_Parameters_To_Galacticus($config,@ARGV[5..$#ARGV]);

# Find the scratch directory.
my $scratchDirectory = $config->{'likelihood'}->{'workDirectory'}."/mcmc";
$scratchDirectory = $config->{'likelihood'}->{'scratchDirectory'}
    if ( exists($config->{'likelihood'}->{'scratchDirectory'}) );

# Expand any environment variable names in the scratch directory.
while ( $scratchDirectory =~ m/\$([_A-Z]+)/ ) {
    my $environmentVariableName  = $1;
    my $environmentVariableValue = $ENV{$environmentVariableName};
    $scratchDirectory =~ s/\$$environmentVariableName/$environmentVariableValue/g;
}

# Ensure scratch and work directories exist.
system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'}."/mcmc")
    unless ( -e $config->{'likelihood'}->{'workDirectory'}."/mcmc" );
system("mkdir -p ".$scratchDirectory)
    unless ( -e $scratchDirectory );

# Report.
if ( exists($config->{'likelihood'}->{'report'}) ) {
    if ( $config->{'likelihood'}->{'report'} eq "yes" ) {
	print "Report from constrainGalacticus.pl:\n";
	print "  MPI rank is : ".$mpiRank."\n";
	print "  Output in   : ".$scratchDirectory."/newParameters_".$mpiRank.".xml\n";
	print "  Parameters  : \n";
	print Dumper($newParameters);
    }
}

# Find the set of base parameters to use.
my $baseParameters = "parameters.xml";
$baseParameters = $config->{'likelihood'}->{'baseParameters'}
   if ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Find the CPU limit.
my $cpuLimit;
if ( exists($config->{'likelihood'}->{'cpulimit'}) ) {
    $cpuLimit = $config->{'likelihood'}->{'cpulimit'};
    if ( exists($config->{'likelihood'}->{'threads'}) ) { 
	$cpuLimit *= $config->{'likelihood'}->{'threads'};
    } else {
	$cpuLimit *= Sys::CPU::cpu_count();
    }
}

# Find the memory limit.
my $memoryLimit;
$memoryLimit = $config->{'likelihood'}->{'memoryLimit'}
   if ( exists($config->{'likelihood'}->{'memoryLimit'}) );

# Extract compilation file and project directories.
my $compilationFile  = $config->{'likelihood'}->{'compilation'  };
my $projectDirectory = $config->{'likelihood'}->{'workDirectory'};
my $galacticusFile   = "constrainGalacticus_".$mpiRank.".hdf5";

# Bad log likelihood (highly improbable) which we will return in failure conditions.
my $badLogLikelihood         = -1.0e30;
my $badLogLikelihoodVariance =  0.0;

# Initialize a list of temporary files to remove after we're finished.
my @temporaryFiles;

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($compilationFile,$baseParameters);
my @constraints = @{$constraintsRef};

# Remove any old semaphore file.
my $semaphoreName = "galacticus";
$semaphoreName = $parameters->{'treeEvolveThreadLockName'}->{'value'}
    if ( exists($parameters->{'treeEvolveThreadLockName'}) );
unlink("/dev/shm/sem.".$semaphoreName)
    if ( -e "/dev/shm/sem.".$semaphoreName );

# Set an output file name.
$parameters->{'galacticusOutputFileName'}->{'value'} = $scratchDirectory."/".$galacticusFile;
push(@temporaryFiles,$parameters->{'galacticusOutputFileName'}->{'value'});

# Set state file names.
my $stateFileRoot;
if ( exists($config->{'likelihood'}->{'saveState'}) && $config->{'likelihood'}->{'saveState'} eq "yes" ) {
    $stateFileRoot = $scratchDirectory."/".$galacticusFile;
    $stateFileRoot =~ s/\.hdf5//;
    $parameters->{'stateFileRoot'}->{'value'} = $stateFileRoot;
    push(@temporaryFiles,$stateFileRoot."*");
}

# Set a random number seed.
$parameters->{'randomSeed'}->{'value'} = int(rand(10000))+1
    unless ( exists($config->{'likelihood'}->{'randomize'}) && $config->{'likelihood'}->{'randomize'} eq "no" );

# If running at a high temperature, modify the number of merger trees per decade.
my $temperatureEffective = 1.0;
if ( 
    exists($parameters->{'mergerTreeBuildTreesPerDecade'}) 
    &&
    $parameters->{'mergerTreeConstructMethod'}->{'value'} eq "build" 
    ) {
    my $treesPerDecadeEffective =
	&List::Util::max(
	    int(
		$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}
		/$temperature
	    )
	    ,$config->{'likelihood'}->{'treesPerDecadeMinimum'}
	);
    $temperatureEffective = $parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}/$treesPerDecadeEffective;
    $parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $treesPerDecadeEffective;
}

# If fixed sets of trees are to be used, create them as necessary, and store to a file.
if ( exists($config->{'likelihood'}->{'useFixedTrees'}) && $config->{'likelihood'}->{'useFixedTrees'} eq "yes" ) {
    # Record the required set of output redshifts.
    my $outputRedshifts = $parameters->{'outputRedshifts'   }->{'value'};
    # Record and remove any analyses.
    my $savedAnalyses   = $parameters->{'mergerTreeAnalyses'}->{'value'};
    delete($parameters->{'mergerTreeAnalyses'});
    # Record merger tree operators.
    my $savedMergerTreeOperator;
    $savedMergerTreeOperator = dclone($parameters->{'mergerTreeOperatorMethod'})
	if ( exists($parameters->{'mergerTreeOperatorMethod'}) );
    # Reset parameters.
    my %parametersToSave =
	(
	 mergerTreePruneBaryons => "false"
	);
    my @savedParameters;
    foreach my $parameterName ( keys(%parametersToSave) ) {
	if ( exists($parameters->{$parameterName}) ) {
	    push
		(
		 @savedParameters,
		 {
		     name  =>               $parameterName            ,
		     value => $parameters->{$parameterName}->{'value'}
		 }
		);
	}
    }
    # Get a lock on the tree file.
    my $fixedTreeDirectory;
    if ( exists($config->{'likelihood'}->{'fixedTreesInScratch'}) && $config->{'likelihood'}->{'fixedTreesInScratch'} eq "yes" ) {
	$fixedTreeDirectory = $config->{'likelihood'}->{'scratchDirectory'}."/";
    } else {
	$fixedTreeDirectory = $config->{'likelihood'}->{'workDirectory'   }."/";
    }
    my $fixedTreeFile      = $fixedTreeDirectory                       .       "fixedTrees".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".hdf5";
    my $buildFixedTreeFile = $config->{'likelihood'}->{'workDirectory'}."/trees/fixedTrees".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".hdf5";
    system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'}."/trees");
    if ( 
	my $lock = new File::NFSLock {
	    file               => $fixedTreeFile,
	    lock_type          => LOCK_EX
	}
	)
    {
	unless ( -e $fixedTreeFile ) {
	    if ( 
		my $buildLock = new File::NFSLock {
		    file               => $buildFixedTreeFile,
		    lock_type          => LOCK_EX
		}
		)
	    {
		unless ( -e $buildFixedTreeFile ) {
		    # Create the tree file if necessary. (Set output redshift to a very large value to avoid any galaxy formation
		    # calculation being carried out - we only want to build the trees.)
		    $parameters->{'outputRedshifts'             }                                  ->{'value'} = "10000.0";
		    $parameters->{'mergerTreePruneBaryons'      }                                  ->{'value'} = "false";
		    my $newOperator;
		    $newOperator                                  ->{'value'} = "export";
		    $newOperator->{'outputFileName'              }->{'value'} = $buildFixedTreeFile;
		    $newOperator->{'mergerTreeExportOutputFormat'}->{'value'} = "galacticus";
		    if ( exists($parameters->{'mergerTreeOperatorMethod'}) ) {
			my @operators;
			if ( $parameters->{'mergerTreeOperatorMethod'}->{'value'} eq "sequence" ) {
			    @operators = &ExtraUtils::as_array($parameters->{'mergerTreeOperatorMethod'}->{'mergerTreeOperatorMethod'});
			    push(@operators,$newOperator);
			} else {
			    @operators = ( dclone($parameters->{'mergerTreeOperatorMethod'}), $newOperator );
			    $parameters->{'mergerTreeOperatorMethod'}->{'value'} = "sequence";
			}
			@{$parameters->{'mergerTreeOperatorMethod'}->{'mergerTreeOperatorMethod'}} = @operators;
		    } else {
			$parameters->{'mergerTreeOperatorMethod'} = $newOperator;
		    }
		    my $treeXML = new XML::Simple (RootName=>"parameters", NoAttr => 1);
		    open(pHndl,">".$config->{'likelihood'}->{'workDirectory'}."/trees/treeBuildParameters".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".xml");
		    print pHndl $treeXML->XMLout($parameters);
		    close pHndl;		  
		    my $treeCommand;
		    $treeCommand .= "ulimit -t ".$cpuLimit."; "
			if ( defined($cpuLimit) );
		    $treeCommand .= "ulimit -v ".$memoryLimit."; "
			if ( defined($memoryLimit) );
		    my $coredump = "YES";
		    $coredump = $config->{'likelihood'}->{'coredump'}
		        if ( exists($config->{'likelihood'}->{'coredump'}) );
		    my $coreDumpSize = "unlimited";
		    $coreDumpSize = 0
			if ( $coredump eq "NO" );
		    $coreDumpSize = $config->{'likelihood'}->{'coredumpsize'}
		        if ( exists($config->{'likelihood'}->{'coredumpsize'}) );
		    $treeCommand .= "ulimit -c ".$coreDumpSize."; export GFORTRAN_ERROR_DUMPCORE=".$coredump."; ulimit -a; date; ./Galacticus.exe ".$config->{'likelihood'}->{'workDirectory'}."/trees/treeBuildParameters".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".xml";
		    my $treeLog = $config->{'likelihood'}->{'workDirectory'}."/trees/treeBuildParameters".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".log";
		    SystemRedirect::tofile($treeCommand,$treeLog);
		    unless ( $? == 0 ) {
			system("mv ".$treeLog." ".$treeLog.".failed.".$$);
			die("constrainGalacticus.pl: Galacticus model failed");		
		    }
		    sleep(1)
			until ( -e $buildFixedTreeFile );
		}
		$buildLock->unlock();
	    }
	    system("cp -f ".$buildFixedTreeFile." ".$fixedTreeFile);
	}
	$lock->unlock();
    }
    # Modify parameters to use the tree file.
    $parameters->{'mergerTreeConstructMethod'}->{'value'} = "read";
    $parameters->{'mergerTreeReadFileName'   }->{'value'} = $fixedTreeFile;
    $parameters->{'outputRedshifts'          }->{'value'} = $outputRedshifts;
    $parameters->{'mergerTreeAnalyses'       }->{'value'} = $savedAnalyses;
    $parameters->{'mergerTreeOperatorMethd'  }            = $savedMergerTreeOperator
	if ( $savedMergerTreeOperator );
    # Restore parameters.
    foreach my $parameter ( @savedParameters ) {
	$parameters->{$parameter->{'name'}}->{'value'} = $parameter->{'value'};
    }
}

# Extract new parameters.
if ( defined($newParameters) ) {
    for my $newParameterName ( keys(%{$newParameters}) ) {
	my $parameter = $parameters;
	foreach ( split(/\-\>/,$newParameterName) ) {
	     $parameter->{$_}->{'value'} = undef()
		 unless ( exists($parameter->{$_}) );
	     $parameter = $parameter->{$_};
    }
	$parameter->{'value'} = $newParameters->{$newParameterName};
    }
}

# Write the modified parameters to file.
&Parameters::Output($parameters,$scratchDirectory."/constrainGalacticusParameters".$mpiRank.".xml");
push(@temporaryFiles,$scratchDirectory."/constrainGalacticusParameters".$mpiRank.".xml");

# Run the Galacticus model.
my $glcCommand;
$glcCommand .= "ulimit -t ".$cpuLimit."; "
    if ( defined($cpuLimit) );
$glcCommand .= "ulimit -v ".$memoryLimit."; "
    if ( defined($memoryLimit) );
$glcCommand .= "export OMP_NUM_THREADS=".$config->{'likelihood'}->{'threads'}."; "
    if ( exists($config->{'likelihood'}->{'threads'}) );
$glcCommand .= "export ".$_."; "
    foreach ( &ExtraUtils::as_array($config->{'likelihood'}->{'environment'}) );
my $coredump = "YES";
$coredump = $config->{'likelihood'}->{'coredump'}
    if ( exists($config->{'likelihood'}->{'coredump'}) );
my $coreDumpSize = "unlimited";
$coreDumpSize = 0
    if ( $coredump eq "NO" );
$coreDumpSize = $config->{'likelihood'}->{'coredumpsize'}
    if ( exists($config->{'likelihood'}->{'coredumpsize'}) );
$glcCommand .= "ulimit -c ".$coreDumpSize."; export GFORTRAN_ERROR_DUMPCORE=".$coredump."; ulimit -a; date;";
# Determine if walltime should be limited using the "timelimit" tool (http://devel.ringlet.net/sysutils/timelimit/).
if (
    exists($config->{'likelihood'}->{'wallTimeLimit'}) 
    && 
    $config->{'likelihood'}->{'wallTimeLimit'} eq "yes"
    &&
    &File::Which::which('timelimit')
    &&
    defined($cpuLimit)
    ) {
    my $wallTimeLimit = int(1.1*$cpuLimit);
    $glcCommand .= " timelimit -t ".$wallTimeLimit." -T 30";
}
$glcCommand .= " ./Galacticus.exe ".$scratchDirectory."/constrainGalacticusParameters".$mpiRank.".xml";
my $logFile = $scratchDirectory."/constrainGalacticusParameters".$mpiRank.".log";
push(@temporaryFiles,$logFile);
# my $timeGalacticusStart = [gettimeofday];
SystemRedirect::tofile($glcCommand,$logFile);
unless ( $? == 0 ) {
    # Since a failure occurred, post to the semaphore to avoid blocking other jobs.
    &semaphorePost($config->{'likelihood'}->{'threads'},$semaphoreName);
    # Check for CPU time being exceeded.
    my $cpuTimeExceeded = 0;
    open(my $log,$logFile);
    while ( my $line = <$log> ) {
	$cpuTimeExceeded = 1
	    if ( $line =~ m/Galacticus exceeded available CPU time/ || $line =~ m/Killed/ );
    }
    close($log);    
    # Check for CPU time limit exceeded.
    if ( $cpuTimeExceeded == 1 ) {
	# Job died due to CPU time being exceeded. No reason to try running it again.
	&reportFailure($config,$scratchDirectory,$logFile,$stateFileRoot,$cpuTimeExceeded);
	# Display the final likelihood.
	&outputLikelihood($config,$badLogLikelihood,$badLogLikelihoodVariance);
	print "constrainGalacticus.pl: Galacticus model failed to complete - second attempt";
	unlink(@temporaryFiles)
	    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
	exit;
    } else {
	# Job died for some other reason, try it again. Issue a failure.
	&reportFailure($config,$scratchDirectory,$logFile,$stateFileRoot,$cpuTimeExceeded);
	# Try running the model again - in case this was a random error.
	if ( exists($config->{'likelihood'}->{'rerunOnError'}) && $config->{'likelihood'}->{'rerunOnError'} eq "yes" ) {
	    print "ERROR: Galacticus model failed to complete - retrying\n";
	    SystemRedirect::tofile($glcCommand,$logFile);
	    unless ( $? == 0 ) {
		# Since a failure occurred, post to the semaphore to avoid blocking other jobs.
		&semaphorePost($config->{'likelihood'}->{'threads'},$semaphoreName);
		# Display the final likelihood.
		&outputLikelihood($config,$badLogLikelihood,$badLogLikelihoodVariance);
		print "constrainGalacticus.pl: Galacticus model failed to complete - second attempt";
		unlink(@temporaryFiles)
		    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
		exit;
	    }
	} else {
		# Since a failure occurred, post to the semaphore to avoid blocking other jobs.
		&semaphorePost($config->{'likelihood'}->{'threads'},$semaphoreName);
		# Display the final likelihood.
		&outputLikelihood($config,$badLogLikelihood,$badLogLikelihoodVariance);
		print "constrainGalacticus.pl: Galacticus model failed to complete";
		unlink(@temporaryFiles)
		    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
		exit;
	}
    }
}
# my $timeGalacticusElapsed = tv_interval($timeGalacticusStart,[gettimeofday]);
# print "%% timing : Galacticus : ".$timeGalacticusElapsed."\n";

# Free model lock if necessary.
$modelLock->unlock()
    if ( $runSequential eq "yes" );

# Perform processing of the model, accumulating likelihood as we go.
my $logLikelihood         =  0.0;
my $logLikelihoodVariance =  0.0;
my $i                     = -1;
my $xml                   = new XML::Simple;
my @threads;
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $constraintDefinition;
    if ( -e $constraint->{'definition'}.".store" ) {
	$constraintDefinition = retrieve($constraint->{'definition'}.".store");
    } else {
	$constraintDefinition = $xml->XMLin($constraint->{'definition'});
    }
    # Run the analysis code.
    ++$i;
    my $analysisCommand = $constraintDefinition->{'analysis'};
    $analysisCommand   .= " ".$scratchDirectory."/".$galacticusFile." --outputFile ".$scratchDirectory."/likelihood".$mpiRank.":".$i.".xml --quiet 1";
    $analysisCommand .= " --temperature ".$temperatureEffective;
    $analysisCommand .= " --modelDiscrepancies ".$projectDirectory."/modelDiscrepancy"
	if ( -e $projectDirectory."/modelDiscrepancy" );
    unless ( $store eq "none" ) {
	my $resultFile = $scratchDirectory."/results".$mpiRank.".xml";
	$analysisCommand .= " --resultFile ".$resultFile;
	push(@temporaryFiles,$resultFile);
    }
    $analysisCommand .= " ".$constraintDefinition->{'analysisArguments'}
        if ( exists($constraintDefinition->{'analysisArguments'}) );
    if ( $useThreads ) {
	my ($thread) = threads->create(\&systemThread,$analysisCommand);
	push(
	    @threads,
	    {
		thread     => $thread,
		constraint => $constraintDefinition
	    }
	    );
    } else {
	system($analysisCommand);
	push(
	    @threads,
	    {
		result     => $?,
		constraint => $constraintDefinition
	    }
	    );
    }
}
$i = -1;
foreach my $constraint ( @constraints ) {
    my $descriptor           = pop(@threads);
    my $constraintDefinition = $descriptor->{'constraint'};
    my $result;
    if ( $useThreads ) {
	my $thread  = $descriptor->{'thread'};
	my @results = $thread->join();
	$result     = $results[0];
    } else {
	$result = $descriptor->{'result'};
    }
    unless ( $result == 0 ) {
	# Issue a failure.
	print "ERROR: Analysis script failed to complete [".$constraint->{'definition'}."]\n";
	&reportFailure($config,$scratchDirectory,$logFile,$stateFileRoot,0);
	# Display the final likelihood.
	&outputLikelihood($config,$badLogLikelihood,$badLogLikelihoodVariance);
	print "constrainGalacticus.pl: analysis code failed";
	unlink(@temporaryFiles)
	    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
	exit;
    }
    # Read the likelihood.
    ++$i;
    my $likelihood = $xml->XMLin($scratchDirectory."/likelihood".$mpiRank.":".$i.".xml");
    if ( $likelihood->{'logLikelihood'} =~ m/nan/i ) {
	# Issue a failure.
        print "ERROR: Likelihood is NaN\n";
	&reportFailure($config,$scratchDirectory,$logFile,$stateFileRoot,0);
	# Display the final likelihood.
	&outputLikelihood($config,$badLogLikelihood,$badLogLikelihoodVariance);
	print "constrainGalacticus.pl: likelihood calculation failed";
	system("rm ".join(" ",@temporaryFiles))
	    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
	exit;
    }
    if ( $likelihood->{'logLikelihood'} =~ m/inf/i ) {
	# Issue a failure.
        print "ERROR: Likelihood is infinity\n";
	&reportFailure($config,$scratchDirectory,$logFile,$stateFileRoot,0);
	# Display the final likelihood.
	&outputLikelihood($config,$badLogLikelihood,$badLogLikelihoodVariance);
	print "constrainGalacticus.pl: likelihood calculation failed";
	system("rm ".join(" ",@temporaryFiles))
	    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
	exit;
    }
    # Extract the likelihood (and variance) and weight it.
    my $thisLogLikelihood         = $likelihood->{'logLikelihood'        };
    my $thisLogLikelihoodVariance = $likelihood->{'logLikelihoodVariance'};
    if ( defined($constraintDefinition->{'weight'}) ) {
	$thisLogLikelihood         *= $constraintDefinition->{'weight'};
	$thisLogLikelihoodVariance *= $constraintDefinition->{'weight'};
    }
    # Accumulate the likelihood and variance.
    $logLikelihood         += $thisLogLikelihood;
    $logLikelihoodVariance += $thisLogLikelihoodVariance;
    # Clean up.
    unlink($scratchDirectory."/likelihood".$mpiRank.":".$i.".xml")
	if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" );
}

# Remove or store the model.
if ( $store eq "none" ) {
    system("rm ".join(" ",@temporaryFiles))
	if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} && scalar(@temporaryFiles) > 0 );
} else {
    my $storeDirectory = $config->{'likelihood'}->{'workDirectory'}."/mcmc/store/model_".$mpiRank."_".$store;
    system("mkdir -p ".$storeDirectory);
    foreach my $file ( @temporaryFiles ) {
	system("mv ".$file." ".$storeDirectory."/");
    }
}

# Display the final likelihood.
&outputLikelihood($config,$logLikelihood,$logLikelihoodVariance);

# Script timing.
# my $timeElapsed = tv_interval($timeStart,[gettimeofday]);
# print "%% timing : ".$self." : ".$timeElapsed."\n";
exit;

sub outputLikelihood {
    # Output the log likelihood.
    my $config                = shift;
    my $logLikelihood         = shift;
    my $logLikelihoodVariance = shift;
    open(oHndl,">".$likelihoodFile);
    print oHndl $logLikelihood        ."\n";
    print oHndl $logLikelihoodVariance."\n";
    close(oHndl);
}

sub reportFailure {
    # Handle failures of model or analysis.
    my $config           = shift;
    my $scratchDirectory = shift;
    my $logFile          = shift;
    my $stateFileRoot    = shift;
    my $cpuTimeExceeded  = shift;
    if ( exists($config->{'likelihood'}->{'failDirectory'}) ) {
	system("mkdir -p ".$config->{'likelihood'}->{'failDirectory'});
	my $failArchiveName = $config->{'likelihood'}->{'failDirectory'}."/archive.tar.bz2";
	$failArchiveName = $config->{'likelihood'}->{'failDirectory'}."/archiveXCPU.tar.bz2"
	    if ( $cpuTimeExceeded == 1 );
	if ( ! -e $failArchiveName && -e "galacticusConfig.xml" ) {
	    # Send an email if possible.
	    my $xml     = new XML::Simple;
	    my $galacticusConfig  = $xml->XMLin("galacticusConfig.xml");
	    if ( exists($galacticusConfig->{'contact'}->{'email'}) && &File::Which::which('sendmail') ) {
		my $message  = "A Galacticus model failed while seeking constraints.\n";
		$message    .= "Failed model is in: ".$failArchiveName."\n";
		my $msg = MIME::Lite->new(
		    From    => '',
		    To      => $galacticusConfig->{'contact'}->{'email'},
		    Subject => 'Galacticus model failed while seeking constraints',
		    Type    => 'TEXT',
		    Data    => $message
		    );
		$msg->send();
	    } else {
		open(my $eHndl,">".$config->{'likelihood'}->{'failDirectory'}."/constrainGalacticusError.txt");
		print $eHndl "A Galacticus model failed while seeking constraints.\n";
		print $eHndl "Failed model is in: ".$failArchiveName."\n";
		close($eHndl);
	    }
	}
	if ( ! -e $failArchiveName ) {
	    system("touch ".$failArchiveName);
	    my $tarCommand = "tar cvfj ".$failArchiveName." ".$scratchDirectory."/constrainGalacticusParameters".$mpiRank.".xml ".$logFile." ".$scratchDirectory."/".$galacticusFile;
	    opendir(dHndl,".");
	    while ( my $fileName = readdir(dHndl) ) {
		$tarCommand .= " ".$fileName
		    if ( $fileName eq "core.".$$ );
	    }
	    closedir(dHndl);
	    $tarCommand .= " ".$stateFileRoot.".*state*"
		if ( $config->{'likelihood'}->{'saveState'} eq "yes" );
	    system($tarCommand);
	}
    }
    if ( exists($config->{'likelihood'}->{'failDirectory'}) ) {
	system("mkdir -p ".$config->{'likelihood'}->{'failDirectory'});
	my $count = 0;
	my $failCountFileName = $config->{'likelihood'}->{'failDirectory'}."/failCount.txt";
	$failCountFileName    = $config->{'likelihood'}->{'failDirectory'}."/failCountXCPU.txt"
	    if ( $cpuTimeExceeded == 1 );
	if ( -e $failCountFileName ) {
	    open(iHndl,$failCountFileName);
	    $count = <iHndl>;
	    close(iHndl);
	    if (defined($count) ) {
		chomp($count);
	    } else {
		$count = 0;
	    }
	}
	++$count;
	open(oHndl,">".$failCountFileName);
	print oHndl $count."\n";
	close(oHndl);
    }
}

sub semaphorePost {
    # When a model fails, it may not be able to release semaphores that it was holding. Therefore, we post to the semaphore to
    # free it.
    my $threadCount   = shift();
    my $semaphoreName = shift();

    # Open the semaphore.
    my $galacticusSemaphore = POSIX::RT::Semaphore->open("/".$semaphoreName, O_CREAT, 0660, $threadCount);
    # Repeatedly post.
    print "Posting to semaphore....\n";
    for(my $i=0;$i<$threadCount;++$i) {
	my $count = $galacticusSemaphore->getvalue();
	print "Current semaphore value is: ".$count."\n";
	$galacticusSemaphore->post();
    }
    my $count = $galacticusSemaphore->getvalue();
    print "Final semaphore value is: ".$count."\n";    
}

sub systemThread {
    # Launch a system command.
    my $command = shift();    
    system($command);
    my $status = $?;
    return $status;
}
