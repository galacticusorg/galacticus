# Launch models on PBS system.

package PBS;
use strict;
use warnings;
use Data::Dumper;
use Sys::CPU;
use File::Which;
require Galacticus::Launch::Hooks;
require Galacticus::Launch::PostProcess;
require Galacticus::Options;
require System::Redirect;
require List::ExtraUtils;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     pbs => {
	 validate       => \&Validate        ,
	 outputFileName => \&Output_File_Name,
	 launch         => \&Launch
     }
    );


sub Validate {
    # Validate the launch script.
    my $launchScript = shift();
    # Set defaults.
    my %defaults = 
	(
	 mpiLaunch               => "yes"            ,
	 mpiRun                  => "mpirun --bynode",
	 maxJobsInQueue          => -1               ,
	 postSubmitSleepDuration => 10               ,
	 jobWaitSleepDuration    => 60               ,
	 analyze                 => "yes"
	);
    # Attempt to detect MPI implementation.
    my $mpiIs = &mpiDetect();
    if ( $mpiIs eq "OpenMPI" ) {
	$defaults{'mpiLaunch'} = "yes"            ;
	$defaults{'mpiRun'   } = "mpirun --bynode";
    }    
    # Apply defaults.
    foreach ( keys(%defaults) ) {
	$launchScript->{'pbs'}->{$_} = $defaults{$_}
	unless ( exists($launchScript->{'pbs'}->{$_}) );
    }
    # Forbid splitting models if analysis is requested on node.
    die("PBS::Validate: can not analyze models on PBS if models are being split")
	if ( $launchScript->{'pbs'}->{'analyze'} eq "yes" && $launchScript->{'splitModels'} > 1 );
}

sub Output_File_Name {
    # Return the output file name for Galacticus models.
    my $outputFileName = shift();
    my $launchScript   = shift();
    if ( exists($launchScript->{'pbs'}->{'scratchPath'}) ) {
	$outputFileName = 
	    $launchScript->{'pbs'}->{'scratchPath'}.
	    "/model_".$launchScript->{'modelCounter'}."_".$$."/galacticus.hdf5";
    }
    return $outputFileName;
}

sub Launch {
    # Launch models on local machine.
    my @jobs         = @{shift()};
    my $launchScript =   shift() ;
    # Iterate over jobs.
    my @modelQueue;
    foreach my $job ( @jobs ) {
	# Create the PBS submission script.
	my $pbsScript = $job->{'directory'}."/pbsLaunch.sh";
	open(my $pbsFile,">".$pbsScript);
	print $pbsFile "#!/bin/bash\n";
	print $pbsFile "#PBS -N Galacticus_".$job->{'label'}."\n";
	if ( exists($launchScript->{'config'}->{'contact'}->{'email'}) ) {
	    if ( $launchScript->{'config'}->{'contact'}->{'email'} =~ m/\@/ && exists($launchScript->{'emailReport'}) && $launchScript->{'emailReport'} eq "yes" ) {
		print $pbsFile "#PBS -M ".$launchScript->{'config'}->{'contact'}->{'email'}."\n";
		print $pbsFile "#PBS -m bea\n";
	    }
	}
	print $pbsFile "#PBS -l walltime=".$launchScript->{'pbs'}->{'wallTime'}."\n"
	    if ( exists($launchScript->{'pbs'}->{'wallTime'}) );
	print $pbsFile "#PBS -l mem=".$launchScript->{'pbs'}->{'memory'}."\n"
	    if ( exists($launchScript->{'pbs'}->{'memory'}) );
	if ( exists($launchScript->{'pbs'}->{'ompThreads'}) ) {
	    print $pbsFile "#PBS -l nodes=1:ppn=".$launchScript->{'pbs'}->{'ompThreads'}."\n";
	} else {
	    print $pbsFile "#PBS -l nodes=1:ppn=1\n";
	}
	print $pbsFile "#PBS -j oe\n";
	my $pwd = "";
	unless ( $job->{'directory'} =~ m/^\// ) {
	    $pwd = `pwd`;
	    chomp($pwd);
	    $pwd .= "/";
	}
	print $pbsFile "#PBS -o ".$pwd.$job->{'directory'}."/galacticus.log\n";
	print $pbsFile "#PBS -q ".$launchScript->{'pbs'}->{'queue'}."\n"
	    if ( exists($launchScript->{'pbs'}->{'queue'}) );
	print $pbsFile "#PBS -V\n";
	print $pbsFile "cd \$PBS_O_WORKDIR\n";
	if ( exists($launchScript->{'pbs'}->{'environment'}) ) {
	    foreach my $environment ( @{$launchScript->{'pbs'}->{'environment'}} ) {
		print $pbsFile "export ".$environment."\n";
	    }
	}
	my $coreDump = "no";
	$coreDump = $launchScript->{'pbs'}->{'coreDump'}
	if ( exists($launchScript->{'pbs'}->{'coreDump'}) );
	if ( $coreDump eq "yes" ) {
	    print $pbsFile "ulimit -c unlimited\n";
	    print $pbsFile "export GFORTRAN_ERROR_DUMPCORE=YES\n";
	} else {
	    print $pbsFile "ulimit -c 0\n";
	    print $pbsFile "export GFORTRAN_ERROR_DUMPCORE=NO\n";
	}
	print $pbsFile "ulimit -t unlimited\n";
	print $pbsFile "export OMP_NUM_THREADS=".$launchScript->{'pbs'}->{'ompThreads'}."\n"
	    if ( exists($launchScript->{'pbs'}->{'ompThreads'}) );
	if ( exists($launchScript->{'pbs'}->{'scratchPath'}) ) {
	    my $scratchPath = $launchScript->{'pbs'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/";
	    print $pbsFile "mkdir -p ".$scratchPath."\n";
	}
	if ( exists($launchScript->{'pbs'}->{'preCommand'}) ) {
	    foreach my $command ( &ExtraUtils::as_array($launchScript->{'pbs'}->{'preCommand'}) ) {
		print $pbsFile $command."\n";
	    }
	}
	print $pbsFile $launchScript->{'pbs'}->{'mpiRun'}." -np 1 "
	    if ( $launchScript->{'pbs'}->{'mpiLaunch'} eq "yes" );
	print $pbsFile "./Galacticus.exe ".$job->{'directory'}."/parameters.xml\n";
	if ( exists($launchScript->{'pbs'}->{'scratchPath'}) ) {
	    print $pbsFile "mv ".$launchScript->{'pbs'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus.hdf5 ".$job->{'directory'}."/galacticus.hdf5\n";
	    if ( $launchScript->{'useStateFile'} eq "yes" ) {
		print $pbsFile "mv ".$launchScript->{'pbs'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus_".$job->{'modelCounter'}."_".$$.".state* ".$job->{'directory'}."/\n";
		print $pbsFile "mv ".$launchScript->{'pbs'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus_".$job->{'modelCounter'}."_".$$.".fgsl.state* ".$job->{'directory'}."/\n";
	    }
	}
	print $pbsFile "mv core* ".$job->{'directory'}."/\n";
	print $pbsFile $job->{'analysis'}
	    if ( defined($job->{'analysis'}) && $launchScript->{'pbs'}->{'analyze'} eq "yes" );
	close($pbsFile);
	# Enter the model in the queue to be run.
	push(
	    @modelQueue,
	    {
		script => $pbsScript,
		job    => $job
	    }
	    );	
    }
    # Submit jobs.
    print " -> waiting for PBS jobs to finish...\n"
	if ( $launchScript->{'verbosity'} > 0 );
    my %pbsJobs;
    while ( scalar(keys %pbsJobs) > 0 || scalar(@modelQueue) > 0 ) {
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
		print " -> PBS job ".$jobID." has finished; post-processing....\n";
		&PostProcess::Analyze($pbsJobs{$jobID}->{'job'},$launchScript)
		    unless ( $launchScript->{'pbs'}->{'analyze'} eq "yes" );
		&PostProcess::CleanUp($pbsJobs{$jobID}->{'job'},$launchScript);
		# Remove the job ID from the list of active PBS jobs.
		delete($pbsJobs{$jobID});
	    }
	}
	# If there are jobs remaining to be submitted, and not too many are already queued, launch one.
	if ( 
	    scalar(@modelQueue) > 0 
	    &&
	    (
	     scalar(keys %pbsJobs) < $launchScript->{'pbs'}->{'maxJobsInQueue'} 
	     ||
	     $launchScript->{'pbs'}->{'maxJobsInQueue'} < 0 
	    )
	    ) {
	    my $pbsJob             = shift(@modelQueue);
	    print " -> submitting script: ".$pbsJob->{'script'}."\n";
	    open(my $pHndl,"qsub ".$pbsJob->{'script'}."|");
	    my $jobID = "";
	    while ( my $line = <$pHndl> ) {
		if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
	    }
	    close($pHndl);
	    # Add job number to active job hash
	    unless ( $jobID eq "" ) {
		$pbsJobs{$jobID}->{'job'} = $pbsJob->{'job'};
	    }
	    sleep $launchScript->{'pbs'}->{'postSubmitSleepDuration'};
	} else {
	    sleep $launchScript->{'pbs'}->{'jobWaitSleepDuration'   };
	}
    }
}

sub mpiDetect {
    # Attempt to detect MPI implementation.
    my $haveMpiRun = which("mpirun");
    if ( defined($haveMpiRun) ) {
	my $mpiVersion = `$haveMpiRun -V 2>&1`;
	return "OpenMPI"
	    if ( $mpiVersion =~ m/Open MPI/ );
    }
}

sub SubmitJobs {
    # Submit jobs to PBS and wait for them to finish.
    my %arguments = %{shift()};
    my @pbsStack  = @_;
    my %pbsJobs;
    # Find the appropriate PBS section.
    my $pbsConfig = &Options::Config("pbs");
    # Determine maximum number allowed in queue at once.
    my $jobMaximum = 10;
    $jobMaximum = $pbsConfig->{'jobMaximum'}
       if ( exists($pbsConfig->{'jbMaximum'}) );
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
	# If fewer than maximum number of jobs are in the queue, pop one off the stack.
	if ( scalar(@pbsStack) > 0 && scalar(keys %pbsJobs) < $jobMaximum ) {
	    my $newJob = pop(@pbsStack);
	    my $batchScript;
	    if ( ref($newJob) ) {
		# Check for a dependency.
		if ( exists($newJob->{'dependsOn'}) && ! -e $newJob->{'dependsOn'} ) {
		    unshift(@pbsStack,$newJob);
		    sleep 5;
		    last;
		}
		# Create the batch script.
		my $ppn = 1;
		$ppn = $newJob->{'ppn'}
		    if ( exists($newJob->{'ppn'}) );
		open(my $scriptFile,">".$newJob->{'launchFile'});
		print $scriptFile "#!/bin/bash\n";
		print $scriptFile "#PBS -N ".$newJob->{'label'}."\n";
		print $scriptFile "#PBS -l walltime=".$newJob->{'walltime'}."\n"
		    if ( exists($newJob->{'wallTime'}) );
		print $scriptFile "#PBS -l mem=".$newJob->{'mem'}."\n"
		    if ( exists($newJob->{'mem'}) );
		print $scriptFile "#PBS -l nodes=1:ppn=".$ppn."\n";
		print $scriptFile "#PBS -j oe\n";
		print $scriptFile "#PBS -o ".$newJob->{'logFile'}."\n";
		print $scriptFile "#PBS -V\n";
		print $scriptFile "cd \$PBS_O_WORKDIR\n";
		print $scriptFile "export LD_LIBRARY_PATH=/home/abenson/Galacticus/Tools/lib:/home/abenson/Galacticus/Tools/lib64:\$LD_LIBRARY_PATH\n";
		print $scriptFile "export ".$_."\n"
		    foreach ( &ExtraUtils::as_array($pbsConfig->{'environment'}) );
		print $scriptFile "ulimit -t unlimited\n";
		print $scriptFile "ulimit -c unlimited\n";
		print $scriptFile "export OMP_NUM_THREADS=".$ppn."\n";
		close($scriptFile);
	    } else {
		# Batch script is already created - we're given its name.
		$batchScript = $newJob;
	    }
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


1;
