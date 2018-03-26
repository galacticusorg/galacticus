#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
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
use PDL::IO::HDF5;
use File::Slurp;
use Storable qw(dclone);
use Fortran::Utils;
use File::Which;
use File::NFSLock;
use System::Redirect;
use Galacticus::Constraints::Parameters;
use Galacticus::Options;
use List::ExtraUtils;
use Sys::CPU;
use Scalar::Util 'reftype';
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
die("Usage: galacticusLikelihood.pl [options....]") 
    unless ( scalar(@ARGV) > 0 );
my %options =
    (
     likelihoodFile        => "./likelihood.dat"    ,
     name                  => "galacticus"          ,
     executable            => "Galacticus.exe"      ,
     workPath              => "./"                  ,
     scratchPath           => "./"                  ,
     threads               => &Sys::CPU::cpu_count(),
     treesPerDecadeMinimum => 10                    ,
     temperature           => 1.0                   ,
     report                => "F"                   ,
     randomize             => "F"                   ,
     saveState             => "F"                   ,
     cleanUp               => "F"                   ,
     coreDump              => "F"                   ,
     coreDumpSize          => "unlimited"           ,
     useFixedTrees         => "F"                   ,
     fixedTreesInScratch   => "F"                   ,
     adjustMasses          => "F"                   ,
     adjustOutputs         => "F"                   ,
     runSequential         => "F"                   ,
     store                 => "F"                   ,
     wallTimeLimit         => "F"                   ,
     rerunOnError          => "F"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Validate options.
foreach ( 'mpiRank', 'parameter', 'compilation' ) {
    die('galacticusLikelihood.xml: option "'.$_.'" must be specified')
	unless ( exists($options{$_}) );
}

# Error codes.
my $errorStatusXCPU = 1025; # CPU time limit exceeded.

# Obtain a lock if we are to run models sequentially.
my $modelLockFile = "/dev/shm/glcLock";
my $modelLock = new File::NFSLock {
    file               => $modelLockFile,
    lock_type          => LOCK_EX
}
if ( $options{'runSequential'} eq "T" );
my $newParameters = &Galacticus::Constraints::Parameters::Convert_Parameters_To_Galacticus($options{'parameter'});

# Expand any environment variable names in the scratch directory.
while ( $options{'scratchPath'} =~ m/\$([_A-Z]+)/ ) {
    my $environmentVariableName  = $1;
    my $environmentVariableValue = $ENV{$environmentVariableName};
    $options{'scratchPath'} =~ s/\$$environmentVariableName/$environmentVariableValue/g;
}

# Ensure scratch and work directories exist.
system("mkdir -p ".$options{'workPath'}."/mcmc")
    unless ( -e $options{'workPath'}."/mcmc" );
system("mkdir -p ".$options{'scratchPath'})
    unless ( -e $options{'scratchPath'} );

# Report.
if ( $options{'report'} eq "T" ) {
    print "Report from galacticusLikelihood.pl:\n";
    print "  MPI rank is : ".$options{'mpiRank'}."\n";
    print "  Output in   : ".$options{'scratchPath'}."/newParameters_".$options{'mpiRank'}.".xml\n";
    print "  Parameters  : \n";
    print Dumper($newParameters);
}

# Scale the CPU limit by the number of threads.
$options{'cpuLimit'} = $options{'cpuLimit'}*$options{'threads'}
   if ( exists($options{'cpuLimit'}) );

# Extract project directory.
my $projectDirectory = $options{'workPath'};
my $galacticusFile   = "galacticusLikelihood_".$options{'mpiRank'}.".hdf5";

# Bad log likelihood (highly improbable) which we will return in failure conditions.
my $badLogLikelihood         = -1.0e30;
my $badLogLikelihoodVariance =  0.0;

# Initialize a list of temporary files to remove after we're finished.
my @temporaryFiles;

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($options{'compilation'},$options{'baseParameters'},$options{'adjustMasses'},$options{'adjustOutputs'});
my @constraints = @{$constraintsRef};

# Remove any old semaphore file.
my $semaphoreName = "galacticus";
$semaphoreName = $parameters->{'treeEvolveThreadLockName'}->{'value'}
    if ( exists($parameters->{'treeEvolveThreadLockName'}) );
unlink("/dev/shm/sem.".$semaphoreName)
    if ( -e "/dev/shm/sem.".$semaphoreName );

# Set an output file name.
$parameters->{'galacticusOutputFileName'}->{'value'} = $options{'scratchPath'}."/".$galacticusFile;
push(@temporaryFiles,$parameters->{'galacticusOutputFileName'}->{'value'});

# Set state file names.
my $stateFileRoot;
if ( $options{'saveState'} eq "T" ) {
    $stateFileRoot = $options{'scratchPath'}."/".$galacticusFile;
    $stateFileRoot =~ s/\.hdf5//;
    $parameters->{'stateFileRoot'}->{'value'} = $stateFileRoot;
    push(@temporaryFiles,$stateFileRoot."*");
}

# Set a random number seed.
$parameters->{'randomSeed'}->{'value'} = int(rand(10000))+1
    if ( $options{'randomize'} eq "T" );

# Build core dump controls.
my $coredump     = $options{'coreDump'} eq "T" ? "YES"                    : "NO";
my $coreDumpSize = $options{'coreDump'} eq "T" ? $options{'coreDumpSize'} : "0" ;

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
		/$options{'temperature'}
	    ),
	    ($options{'adjustMasses'} eq "T" ? $options{'treesPerDecadeMinimum'} : 0)
	);
    $temperatureEffective = $parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}/$treesPerDecadeEffective;
    $parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $treesPerDecadeEffective;
}

# If fixed sets of trees are to be used, create them as necessary, and store to a file.
if ( $options{'useFixedTrees'} eq "T" ) {
    # Record the required set of output redshifts.
    my $outputRedshifts = $parameters->{'outputRedshifts' }              ->{'value'};
    # Record and remove any analyses.
    my $savedAnalyses   = $parameters->{'mergerTreeOutput'}->{'analyses'}->{'value'};
    delete($parameters->{'mergerTreeOutput'}->{'analyses'});
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
    if ( $options{'fixedTreesInScratch'} eq "T" ) {
	$fixedTreeDirectory = $options{'scratchPath'}."/";
    } else {
	$fixedTreeDirectory = $options{'workPath'   }."/";
    }
    my $fixedTreeFile      = $fixedTreeDirectory .       "fixedTrees".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".hdf5";
    my $buildFixedTreeFile = $options{'workPath'}."/trees/fixedTrees".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".hdf5";
    system("mkdir -p ".$options{'workPath'}."/trees");
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
		    $parameters->{'outputRedshifts'             } ->{'value'} = "10000.0";
		    $parameters->{'mergerTreePruneBaryons'      } ->{'value'} = "false";
		    my $newOperator;
		    $newOperator                                  ->{'value'} = "export";
		    $newOperator->{'outputFileName'              }->{'value'} = $buildFixedTreeFile;
		    $newOperator->{'mergerTreeExportOutputFormat'}->{'value'} = "galacticus";
		    if ( exists($parameters->{'mergerTreeOperatorMethod'}) ) {
			my @operators;
			if ( $parameters->{'mergerTreeOperatorMethod'}->{'value'} eq "sequence" ) {
			    @operators = &List::ExtraUtils::as_array($parameters->{'mergerTreeOperatorMethod'}->{'mergerTreeOperatorMethod'});
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
		    open(pHndl,">".$options{'workPath'}."/trees/treeBuildParameters".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".xml");
		    print pHndl $treeXML->XMLout($parameters);
		    close pHndl;		  
		    my $treeCommand;
		    $treeCommand .= "ulimit -t ".$options{'cpuLimit'}."; "
			if ( exists($options{'cpuLimit'}) );
		    $treeCommand .= "ulimit -v ".$options{'memoryLimit'}."; "
			if ( exists($options{'memoryLimit'}) );
		    $treeCommand .= "ulimit -c ".$coreDumpSize."; export GFORTRAN_ERROR_DUMPCORE=".$coredump."; ulimit -a; date; ".$options{'executable'}." ".$options{'workPath'}."/trees/treeBuildParameters".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".xml";
		    my $treeLog = $options{'workPath'}."/trees/treeBuildParameters".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".log";
		    &System::Redirect::tofile($treeCommand,$treeLog);
		    unless ( $? == 0 ) {
			system("mv ".$treeLog." ".$treeLog.".failed.".$$);
			die("galacticusLikelihood.pl: Galacticus model failed");		
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
    $parameters->{'mergerTreeConstructMethod'}              ->{'value'} = "read";
    $parameters->{'mergerTreeReadFileName'   }              ->{'value'} = $fixedTreeFile;
    $parameters->{'outputRedshifts'          }              ->{'value'} = $outputRedshifts;
    $parameters->{'mergerTreeOutput'         }->{'analyses'}->{'value'} = $savedAnalyses;
    $parameters->{'mergerTreeOperatorMethod' }                          = $savedMergerTreeOperator
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
	my $valueIndex;
	if ( $newParameterName =~ m/^(.*)\{(\d+)\}$/ ) {
	    $newParameterName = $1;
	    $valueIndex = $2;
	}
	foreach ( split(/::/,$newParameterName) ) {
	    # Check if the parameter name contains an array reference.
	    if ( $_ =~ m/^(.*)\[(\d+)\]$/ ) {
		# Parameter name contains array reference. Step through to the relevant parameter in the list. If the parameter is
		# not an array, allow this only if the array index given is zero.
		if ( reftype($parameter->{$1}) eq "ARRAY" ) {
		    $parameter->{$1}->[$2]->{'value'} = undef()
			unless ( scalar(@{$parameter->{$1}}) > $2 );
		    $parameter = $parameter->{$1}->[$2];
		} else {
		    die('galacticusLikelihood.pl: attempt to access non-existant array')
			unless ( $2 == 0 );
		    $parameter->{$1}->{'value'} = undef()
			unless ( exists($parameter->{$1}) );
		    $parameter = $parameter->{$1};
		}
	    } else {
		# Parameter does not contain an array reference - so simply step through to the named parameter.
		$parameter->{$_}->{'value'} = undef()
		    unless ( exists($parameter->{$_}) );
		$parameter = $parameter->{$_};
	    }
	}
	# Test if the parameter name contains a value index.
	if ( defined($valueIndex) ) {
	    # A value index is given - set the relevant entry.
	    my @values = split(" ",$parameter->{'value'});
	    $values[$valueIndex] = $newParameters->{$newParameterName."{".$valueIndex."}"};
	    $parameter->{'value'} = join(" ",@values);
	} else {
	    # No value index is given - simply set the value of the parameter.
	    $parameter->{'value'} = $newParameters->{$newParameterName};
	}    
    }
}

# Write the modified parameters to file.
&Galacticus::Constraints::Parameters::Output($parameters,$options{'scratchPath'}."/galacticusLikelihoodParameters".$options{'mpiRank'}.".xml");
push(@temporaryFiles,$options{'scratchPath'}."/galacticusLikelihoodParameters".$options{'mpiRank'}.".xml");

# Run the Galacticus model.
my $glcCommand;
$glcCommand .= "ulimit -t ".$options{'cpuLimit'}."; "
    if ( exists($options{'cpuLimit'}) );
$glcCommand .= "ulimit -v ".$options{'memoryLimit'}."; "
    if ( exists($options{'memoryLimit'}) );
$glcCommand .= "export OMP_NUM_THREADS=".$options{'threads'}."; ";
if ( exists($options{'environment'}) ) {
    $glcCommand .= "export ".$_."; "
	foreach ( &List::ExtraUtils::as_array($options{'environment'}) );
}
$glcCommand .= "ulimit -c ".$coreDumpSize."; export GFORTRAN_ERROR_DUMPCORE=".$coredump."; ulimit -a; date;";
# Determine if walltime should be limited using the "timelimit" tool (http://devel.ringlet.net/sysutils/timelimit/).
if (
    $options{'wallTimeLimit'} eq "T"
    &&
    &File::Which::which('timelimit')
    &&
    exists($options{'cpuLimit'})
    ) {
    my $wallTimeLimit = int(1.1*$options{'cpuLimit'});
    $glcCommand .= " timelimit -t ".$wallTimeLimit." -T 30";
}
$glcCommand .= " ".$options{'executable'}." ".$options{'scratchPath'}."/galacticusLikelihoodParameters".$options{'mpiRank'}.".xml";
my $logFile = $options{'scratchPath'}."/galacticusLikelihoodParameters".$options{'mpiRank'}.".log";
push(@temporaryFiles,$logFile);
# my $timeGalacticusStart = [gettimeofday];
&System::Redirect::tofile($glcCommand,$logFile);
unless ( $? == 0 ) {
    # Since a failure occurred, post to the semaphore to avoid blocking other jobs.
    &semaphorePost($options{'threads'},$semaphoreName);
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
	&reportFailure($options{'scratchPath'},$logFile,$stateFileRoot,$cpuTimeExceeded);
	# Display the final likelihood.
	&outputLikelihood($badLogLikelihood,$badLogLikelihoodVariance);
	print "galacticusLikelihood.pl: Galacticus model failed to complete - second attempt";
	unlink(@temporaryFiles)
	    if ( $options{'cleanUp'} eq "T" && scalar(@temporaryFiles) > 0 );
	exit;
    } else {
	# Job died for some other reason, try it again. Issue a failure.
	&reportFailure($options{'scratchPath'},$logFile,$stateFileRoot,$cpuTimeExceeded);
	# Try running the model again - in case this was a random error.
	if ( $options{'rerunOnError'} eq "T" ) {
	    print "ERROR: Galacticus model failed to complete - retrying\n";
	    &System::Redirect::tofile($glcCommand,$logFile);
	    unless ( $? == 0 ) {
		# Since a failure occurred, post to the semaphore to avoid blocking other jobs.
		&semaphorePost($options{'threads'},$semaphoreName);
		# Display the final likelihood.
		&outputLikelihood($badLogLikelihood,$badLogLikelihoodVariance);
		print "galacticusLikelihood.pl: Galacticus model failed to complete - second attempt";
		unlink(@temporaryFiles)
		    if ( $options{'cleanUp'} eq "T" && scalar(@temporaryFiles) > 0 );
		exit;
	    }
	} else {
		# Since a failure occurred, post to the semaphore to avoid blocking other jobs.
		&semaphorePost($options{'threads'},$semaphoreName);
		# Display the final likelihood.
		&outputLikelihood($badLogLikelihood,$badLogLikelihoodVariance);
		print "galacticusLikelihood.pl: Galacticus model failed to complete";
		unlink(@temporaryFiles)
		    if ( $options{'cleanUp'} eq "T" && scalar(@temporaryFiles) > 0 );
		exit;
	}
    }
}
# my $timeGalacticusElapsed = tv_interval($timeGalacticusStart,[gettimeofday]);
# print "%% timing : Galacticus : ".$timeGalacticusElapsed."\n";

# Free model lock if necessary.
$modelLock->unlock()
    if ( $options{'runSequential'} eq "T" );

# Store raw XML parameter file in the model.
&storeXML($parameters->{'galacticusOutputFileName'}->{'value'},$options{'scratchPath'}."/galacticusLikelihoodParameters".$options{'mpiRank'}.".xml");

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
    $analysisCommand   .= " ".$options{'scratchPath'}."/".$galacticusFile." --outputFile ".$options{'scratchPath'}."/likelihood".$options{'mpiRank'}.":".$i.".xml --quiet 1";
    $analysisCommand .= " --temperature ".$temperatureEffective;
    $analysisCommand .= " --modelDiscrepancies ".$projectDirectory."/modelDiscrepancy"
	if ( -e $projectDirectory."/modelDiscrepancy" );
    if ( exists($options{'massLimit'}) ) {
	foreach ( keys(%{$options{'massLimit'}}) ) {
	    $analysisCommand .= " --".$_."Minimum ".$options{'massLimit'}->{$_};
	}
    }
    unless ( $options{'store'} eq "F" ) {
	my $resultFile = $options{'scratchPath'}."/results".$options{'mpiRank'}.":".$i.".hdf5";
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
    my $descriptor           = shift(@threads);
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
	&reportFailure($options{'scratchPath'},$logFile,$stateFileRoot,0);
	# Display the final likelihood.
	&outputLikelihood($badLogLikelihood,$badLogLikelihoodVariance);
	print "galacticusLikelihood.pl: analysis code failed\n";
	unlink(@temporaryFiles)
	    if ( $options{'cleanUp'} eq "T" && scalar(@temporaryFiles) > 0 );
	my $ignoredResults = $_->{'thread'}->join()
	    foreach ( @threads );
	exit;
    }
    # Read the likelihood.
    ++$i;
    my $waitTime = 0; # Wait for the file to appear if necessary.
    while ( $waitTime < 10 && ! -e $options{'scratchPath'}."/likelihood".$options{'mpiRank'}.":".$i.".xml" ) {
	sleep(1);
	++$waitTime;
    }
    sleep(3);
    my $likelihood = $xml->XMLin($options{'scratchPath'}."/likelihood".$options{'mpiRank'}.":".$i.".xml");
    if ( $likelihood->{'logLikelihood'} =~ m/nan/i ) {
	# Issue a failure.
        print "ERROR: Likelihood is NaN\n";
	&reportFailure($options{'scratchPath'},$logFile,$stateFileRoot,0);
	# Display the final likelihood.
	&outputLikelihood($badLogLikelihood,$badLogLikelihoodVariance);
	print "galacticusLikelihood.pl: likelihood calculation failed";
	system("rm ".join(" ",@temporaryFiles))
	    if ( $options{'cleanUp'} eq "T" && scalar(@temporaryFiles) > 0 );
	my $ignoredResults = $_->{'thread'}->join()
	    foreach ( @threads );
	exit;
    }
    if ( $likelihood->{'logLikelihood'} =~ m/inf/i ) {
	# Issue a failure.
        print "ERROR: Likelihood is infinity\n";
	&reportFailure($options{'scratchPath'},$logFile,$stateFileRoot,0);
	# Display the final likelihood.
	&outputLikelihood($badLogLikelihood,$badLogLikelihoodVariance);
	print "galacticusLikelihood.pl: likelihood calculation failed";
	system("rm ".join(" ",@temporaryFiles))
	    if ( $options{'cleanUp'} eq "T" && scalar(@temporaryFiles) > 0 );
	my $ignoredResults = $_->{'thread'}->join()
	    foreach ( @threads );
	exit;
    }
    $likelihood->{'logLikelihoodVariance'} = 0.0
	if ( $likelihood->{'logLikelihoodVariance'} =~ m/nan/i || $likelihood->{'logLikelihoodVariance'} =~ m/infinity/i );
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
    unlink($options{'scratchPath'}."/likelihood".$options{'mpiRank'}.":".$i.".xml")
	if ( $options{'cleanUp'} eq "T" );
}

# Remove or store the model.
if ( $options{'store'} eq "F" ) {
    system("rm ".join(" ",@temporaryFiles))
	if ( $options{'cleanUp'} eq "T" && scalar(@temporaryFiles) > 0 );
} else {
    my $storeDirectory = $options{'workPath'}."/mcmc/store/model_".$options{'mpiRank'}."_".$options{'store'};
    system("mkdir -p ".$storeDirectory);
    foreach my $file ( @temporaryFiles ) {
	system("mv ".$file." ".$storeDirectory."/");
    }
}
# Display the final likelihood.
&outputLikelihood($logLikelihood,$logLikelihoodVariance);

# Script timing.
# my $timeElapsed = tv_interval($timeStart,[gettimeofday]);
# print "%% timing : ".$self." : ".$timeElapsed."\n";
exit;

sub outputLikelihood {
    # Output the log likelihood.
    my $logLikelihood         = shift;
    my $logLikelihoodVariance = shift;
    open(oHndl,">".$options{'likelihoodFile'});
    print oHndl $logLikelihood        ."\n";
    print oHndl $logLikelihoodVariance."\n";
    close(oHndl);
}

sub reportFailure {
    # Handle failures of model or analysis.
    my $scratchPath      = shift;
    my $logFile          = shift;
    my $stateFileRoot    = shift;
    my $cpuTimeExceeded  = shift;
    if ( exists($options{'failPath'}) ) {
	system("mkdir -p ".$options{'failPath'});
	my $failArchiveName = $options{'failPath'}."/archive.tar.bz2";
	$failArchiveName = $options{'failPath'}."/archiveXCPU.tar.bz2"
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
		open(my $eHndl,">".$options{'failPath'}."/galacticusLikelihoodError.txt");
		print $eHndl "A Galacticus model failed while seeking constraints.\n";
		print $eHndl "Failed model is in: ".$failArchiveName."\n";
		close($eHndl);
	    }
	}
	if ( ! -e $failArchiveName ) {
	    system("touch ".$failArchiveName);
	    my $tarCommand = "tar cvfj ".$failArchiveName." ".$scratchPath."/galacticusLikelihoodParameters".$options{'mpiRank'}.".xml ".$logFile." ".$scratchPath."/".$galacticusFile;
	    opendir(dHndl,".");
	    while ( my $fileName = readdir(dHndl) ) {
		$tarCommand .= " ".$fileName
		    if ( $fileName eq "core.".$$ );
	    }
	    closedir(dHndl);
	    $tarCommand .= " ".$stateFileRoot.".*state*"
		if ( $options{'saveState'} eq "T" );
	    system($tarCommand);
	}
    }
    if ( exists($options{'failPath'}) ) {
	system("mkdir -p ".$options{'failPath'});
	my $count = 0;
	my $failCountFileName = $options{'failPath'}."/failCount.txt";
	$failCountFileName    = $options{'failPath'}."/failCountXCPU.txt"
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
    if ( -e "/dev/shm/sem.".$semaphoreName ) {
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
}

sub systemThread {
    # Launch a system command.
    my $command = shift();    
    system($command);
    my $status = $?;
    return $status;
}

sub storeXML {
    my $galacticusFileName = shift();
    my $parametersFileName = shift();
    my $galacticusModel = new PDL::IO::HDF5(">".$galacticusFileName);
    my $parametersGroup = $galacticusModel->group('Parameters');
    my $parametersRaw   = read_file($parametersFileName);
    $parametersGroup->attrSet('rawXML' => $parametersRaw);    
}
