# Launch models on SLURM system.

package Galacticus::Launch::Slurm;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Cwd;
use File::Which;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PostProcess;
use System::Redirect;
use List::ExtraUtils;

# Insert hooks for our functions.
%Galacticus::Launch::Hooks::moduleHooks = 
    (
     %Galacticus::Launch::Hooks::moduleHooks,
     slurm => {
	       validate       => \&Validate        ,
               outputFileName => \&Output_File_Name,
	       launch         => \&Launch          ,
	       jobArrayLaunch => \&SubmitJobs
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
	$launchScript->{'slurm'}->{$_} = $defaults{$_}
	unless ( exists($launchScript->{'slurm'}->{$_}) );
    }
    # Forbid splitting models if analysis is requested on node.
    die("SLURM::Validate: can not analyze models on SLURM if models are being split")
	if ( $launchScript->{'slurm'}->{'analyze'} eq "yes" && $launchScript->{'splitModels'} > 1 );
}

sub Output_File_Name {
    # Return the output file name for Galacticus models.
    my $outputFileName = shift();
    my $launchScript   = shift();
    if ( exists($launchScript->{'slurm'}->{'scratchPath'}) ) {
	$outputFileName = 
	    $launchScript->{'slurm'}->{'scratchPath'}.
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
	# Create the SLURM submission script.
	my $slurmScript = $job->{'directory'}."/slurmLaunch.sh";
	open(my $slurmFile,">".$slurmScript);
	print $slurmFile "#!/bin/bash\n";
	print $slurmFile "#SBATCH -J Galacticus_".$job->{'label'}."\n";
	my $currentDirectory = cwd();
	print $slurmFile "#SBATCH --chdir=".$currentDirectory."\n";
	if ( exists($launchScript->{'config'}->{'contact'}->{'email'}) ) {
	    if ( $launchScript->{'config'}->{'contact'}->{'email'} =~ m/\@/ && exists($launchScript->{'emailReport'}) && $launchScript->{'emailReport'} eq "yes" ) {
		print $slurmFile "#SBATCH --mail-user=".$launchScript->{'config'}->{'contact'}->{'email'}."\n";
		print $slurmFile "#SBATCH --mail-type=end\n";
	    }
	}
	print $slurmFile "#SBATCH --time=".$launchScript->{'slurm'}->{'wallTime'}."\n"
	    if ( exists($launchScript->{'slurm'}->{'wallTime'}) );
	print $slurmFile "#SBATCH -l mem=".$launchScript->{'slurm'}->{'memory'}."\n"
	    if ( exists($launchScript->{'slurm'}->{'memory'}) );
	print $slurmFile "#SBATCH --nodes=1\n";
	if ( exists($launchScript->{'slurm'}->{'ompThreads'}) ) {
	    print $slurmFile "#SBATCH --cpus-per-task=".$launchScript->{'slurm'}->{'ompThreads'}."\n";
	} else {
	    print $slurmFile "#SBATCH --cpus-per-task=1\n";
	}
	my $pwd = "";
	unless ( $job->{'directory'} =~ m/^\// ) {
	    $pwd = `pwd`;
	    chomp($pwd);
	    $pwd .= "/";
	}
	print $slurmFile "#SBATCH -o ".$pwd.$job->{'directory'}."/galacticus.log\n";
	print $slurmFile "#SBATCH -A ".$launchScript->{'slurm'}->{'account'}."\n"
	    if ( exists($launchScript->{'slurm'}->{'account'}) );
	if ( exists($launchScript->{'slurm'}->{'environment'}) ) {
	    foreach my $environment ( @{$launchScript->{'slurm'}->{'environment'}} ) {
		print $slurmFile "export ".$environment."\n";
	    }
	}
	if ( exists($launchScript->{'slurm'}->{'module'}) ) {
	    foreach my $module ( @{$launchScript->{'slurm'}->{'module'}} ) {
		print $slurmFile "module load ".$module."\n";
	    }
	}
	my $coreDump = "no";
	$coreDump = $launchScript->{'slurm'}->{'coreDump'}
	if ( exists($launchScript->{'slurm'}->{'coreDump'}) );
	if ( $coreDump eq "yes" ) {
	    print $slurmFile "ulimit -c unlimited\n";
	    print $slurmFile "export GFORTRAN_ERROR_DUMPCORE=YES\n";
	} else {
	    print $slurmFile "ulimit -c 0\n";
	    print $slurmFile "export GFORTRAN_ERROR_DUMPCORE=NO\n";
	}
	print $slurmFile "ulimit -t unlimited\n";
	print $slurmFile "pwd\n";
	print $slurmFile "export OMP_NUM_THREADS=".$launchScript->{'slurm'}->{'ompThreads'}."\n"
	    if ( exists($launchScript->{'slurm'}->{'ompThreads'}) );
	if ( exists($launchScript->{'slurm'}->{'scratchPath'}) ) {
	    my $scratchPath = $launchScript->{'slurm'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/";
	    print $slurmFile "mkdir -p ".$scratchPath."\n";
	}
	if ( exists($launchScript->{'slurm'}->{'preCommand'}) ) {
	    foreach my $command ( &List::ExtraUtils::as_array($launchScript->{'slurm'}->{'preCommand'}) ) {
		print $slurmFile $command."\n";
	    }
	}
	print $slurmFile $launchScript->{'slurm'}->{'mpiRun'}." -np 1 "
	    if ( $launchScript->{'slurm'}->{'mpiLaunch'} eq "yes" );
	print $slurmFile $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$job->{'directory'}."/parameters.xml\n";
	if ( exists($launchScript->{'slurm'}->{'scratchPath'}) ) {
	    print $slurmFile "mv ".$launchScript->{'slurm'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus.hdf5 ".$job->{'directory'}."/galacticus.hdf5\n";
	    if ( $launchScript->{'useStateFile'} eq "yes" ) {
		print $slurmFile "mv ".$launchScript->{'slurm'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus_".$job->{'modelCounter'}."_".$$.".state* ".$job->{'directory'}."/\n";
	    }
	}
	print $slurmFile "mv core* ".$job->{'directory'}."/\n";
	print $slurmFile $job->{'analysis'}
	    if ( defined($job->{'analysis'}) && $launchScript->{'slurm'}->{'analyze'} eq "yes" );
	close($slurmFile);
	# Enter the model in the queue to be run.
	push(
	    @modelQueue,
	    {
		script => $slurmScript,
		job    => $job
	    }
	    );	
    }
    # Submit jobs.
    print " -> waiting for SLURM jobs to finish...\n"
	if ( $launchScript->{'verbosity'} > 0 );
    my %slurmJobs;
    while ( scalar(keys %slurmJobs) > 0 || scalar(@modelQueue) > 0 ) {
	# Find all SLURM jobs that are running.
	my %runningSLURMJobs;
	undef(%runningSLURMJobs);
	open(pHndl,"squeue|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^\s+(\S+)/ ) {$runningSLURMJobs{$1} = 1};
	}
	close(pHndl);
	foreach my $jobID ( keys(%slurmJobs) ) {
	    unless ( exists($runningSLURMJobs{$jobID}) ) {
		print " -> SLURM job ".$jobID." has finished; post-processing....\n";
		&Galacticus::Launch::PostProcess::Analyze($slurmJobs{$jobID}->{'job'},$launchScript)
		    unless ( $launchScript->{'slurm'}->{'analyze'} eq "yes" );
		&Galacticus::Launch::PostProcess::CleanUp($slurmJobs{$jobID}->{'job'},$launchScript);
		# Remove the job ID from the list of active SLURM jobs.
		delete($slurmJobs{$jobID});
	    }
	}
	# If there are jobs remaining to be submitted, and not too many are already queued, launch one.
	if ( 
	    scalar(@modelQueue) > 0 
	    &&
	    (
	     scalar(keys %slurmJobs) < $launchScript->{'slurm'}->{'maxJobsInQueue'} 
	     ||
	     $launchScript->{'slurm'}->{'maxJobsInQueue'} < 0 
	    )
	    ) {
	    my $slurmJob             = shift(@modelQueue);
	    print " -> submitting script: ".$slurmJob->{'script'}."\n";
	    open(my $pHndl,"sbatch ".$slurmJob->{'script'}."|");
	    my $jobID = "";
	    while ( my $line = <$pHndl> ) {
		if ( $line =~ m/^Submitted batch job (\d+)/ ) {$jobID = $1};
	    }
	    close($pHndl);
	    # Add job number to active job hash
	    unless ( $jobID eq "" ) {
		print "    -> job ID: ".$jobID."\n";
		$slurmJobs{$jobID}->{'job'} = $slurmJob->{'job'};
	    }
	    sleep $launchScript->{'slurm'}->{'postSubmitSleepDuration'};
	} else {
	    sleep $launchScript->{'slurm'}->{'jobWaitSleepDuration'   };
	}
    }
}

sub mpiDetect {
    # Attempt to detect MPI implementation.
    my $haveMpiRun = &File::Which::which("mpirun");
    if ( defined($haveMpiRun) ) {
	my $mpiVersion = `$haveMpiRun -V 2>&1`;
	return "OpenMPI"
	    if ( $mpiVersion =~ m/Open MPI/ );
    }
}

sub SubmitJobs {
    # Submit jobs to SLURM and wait for them to finish.
    my %arguments = %{shift()};
    my @jobStack  = @_;
    my %slurmJobs;
    # Find the appropriate PBS section.
    my $slurmConfig = &Galacticus::Options::Config("slurm");
    # Determine sleep times between jobs.
    my $submitSleepDuration = exists($arguments{'submitSleepDuration'}) ? $arguments{'submitSleepDuration'} :  5;
    my $waitSleepDuration   = exists($arguments{'waitSleepDuration'  }) ? $arguments{'waitSleepDuration'  } : 60;
    # Determine maximum number allowed in queue at once.
    my $jobMaximum = 10;
    $jobMaximum = $slurmConfig->{'jobMaximum'}
       if ( exists($slurmConfig->{'jobMaximum'}) );
    $jobMaximum = $arguments{'slurmJobMaximum'}
       if ( exists($arguments{'slurmJobMaximum'}) );
    my $limitQueuedOnly = exists($arguments{'limitQueuedOnly'}) ? $arguments{'limitQueuedOnly'} : 0;
    # Submit jobs and wait.
    print "Waiting for SLURM jobs to finish...\n";
    my @activeStates = ( "R", "PD" );
    while ( scalar(keys %slurmJobs) > 0 || scalar(@jobStack) > 0 ) {
	# Find all SLURM jobs that are running.
	my %runningSLURMJobs;
	undef(%runningSLURMJobs);
	my $jobID;
	my $squeueStatus = open(my $squeue,"squeue |");
	unless ( defined($squeueStatus) ) {
	    # If squeue pipe failed, sleep for a while, then start over.
	    print "Could not open pipe to 'squeue' - will sleep and then try again\n";
	    sleep $waitSleepDuration;
	    next;
	}
	my $line = <$squeue>;
	while ( my $line = <$squeue> ) {
	    my @columns = split(" ",$line);
	    $runningSLURMJobs{$columns[0]} = $columns[4]
		if ( grep {$_ eq $columns[4]} @activeStates );
	}
	close($squeue);
	foreach my $jobID ( keys(%slurmJobs) ) {
	    unless ( exists($runningSLURMJobs{$jobID}) ) {
		print "SLURM job ".$jobID." has finished.\n";
		# Call any "on completion" function.
		if ( exists($slurmJobs{$jobID}->{'onCompletion'}) ) {
		    my $exitStatus = 0;
		    if ( exists($slurmJobs{$jobID}->{'tracejob'}) && $slurmJobs{$jobID}->{'tracejob'} eq "yes" ) {
			my $sacct = &File::Which::which("sacct");
			if ( defined($sacct) ) {
			    my $traceStatus = open(my $trace,$sacct." -b -j ".$jobID." |");
			    if ( defined($traceStatus) ) {
				while ( my $line = <$trace> ) {
				    my @columns = split(" ",$line);
				    if ( $columns[0] eq $jobID ) {
					if ( $columns[2] =~ m/(\d+):\d+/ ) {
					    $exitStatus = $1;
					}
				    }
				}
				close($trace);
			    } else {
				$exitStatus = -1;
			    }
			}
		    }
		    &{$slurmJobs{$jobID}->{'onCompletion'}->{'function'}}(@{$slurmJobs{$jobID}->{'onCompletion'}->{'arguments'}},$jobID,$exitStatus,\@jobStack,$slurmJobs{$jobID});
		}
		# Remove the job ID from the list of active SLURM jobs.
		delete($slurmJobs{$jobID});		
	    }
	}
	# If fewer than maximum number of jobs are in the queue, pop one off the stack.
	my $limitingJobs = $limitQueuedOnly ? scalar(grep {($runningSLURMJobs{$_} eq "PD") && exists($slurmJobs{$_})} keys(%runningSLURMJobs)) : scalar(keys(%slurmJobs));
	if ( scalar(@jobStack) > 0 && $limitingJobs < $jobMaximum ) {
	    my $newJob = pop(@jobStack);
	    my $batchScript;
	    if ( ref($newJob) ) {
		# Check for a dependency.
		if ( exists($newJob->{'dependsOn'}) && ! -e $newJob->{'dependsOn'} ) {
		    unshift(@jobStack,$newJob);
		    sleep 5;
		    last;
		}
		# Create the batch script.
		my $resourceModel = exists($newJob->{'resourceModel'}) ? $newJob->{'resourceModel'} : "nodes";
		my $ppn = 1;
		if ( exists($newJob->{'ppn'}) ) {
		    $ppn = $newJob->{'ppn'};
		} elsif ( exists($slurmConfig->{'ppn'}) ) {
		    $ppn = $slurmConfig->{'ppn'};
		}
		my $nodes = 1;
		if ( exists($newJob->{'nodes'}) ) {
		    $nodes = $newJob->{'nodes'};
		} elsif ( exists($slurmConfig->{'nodes'}) ) {
		    $nodes = $slurmConfig->{'nodes'};
		}
		my $ompThreads = 1;
		if ( exists($newJob->{'ompThreads'}) ) {
		    $ompThreads = $newJob->{'ompThreads'};
		} elsif ( exists($slurmConfig->{'ompThreads'}) ) {
		    $ompThreads = $slurmConfig->{'ompThreads'};
		}
		open(my $scriptFile,">".$newJob->{'launchFile'});
		print $scriptFile "#!/bin/bash\n";
		print $scriptFile "#SBATCH --job-name=\"".$newJob->{'label'}."\"\n";
		if ( exists($newJob->{'walltime'}) ) {
		    print $scriptFile "#SBATCH --time=".$newJob->{'walltime'}."\n"
		} elsif ( exists($slurmConfig->{'timeMaximum'}) ) {
		    print $scriptFile "#SBATCH --time=".$slurmConfig->{'timeMaximum'}."\n"
		}
		if ( exists($newJob->{'queue'}) ) {
		    print $scriptFile "#SBATCH --partition=".$newJob->{'queue'}."\n";
		} elsif ( exists($arguments{'queue'}) ) {
		    print $scriptFile "#SBATCH --partition=".$arguments{'queue'}."\n";
		} elsif ( exists($slurmConfig->{'queue'}) ) {
		    print $scriptFile "#SBATCH --partition=".$slurmConfig->{'queue'}."\n";
		}
		if ( $resourceModel eq "nodes" ) {
		    print $scriptFile "#SBATCH --nodes=".$nodes."\n";
		    print $scriptFile "#SBATCH --ntasks-per-node=".$ppn."\n";		    
		    print $scriptFile "#SBATCH --mem=".$newJob->{'mem'}."\n"
			if ( exists($newJob->{'mem'}) );
		} else {
		    die("Galacticus::Launch::SLURM::SubmitJobs: unknown resource model");
		}
		print $scriptFile "#SBATCH --error=".$newJob->{'logFile'}."\n";
		print $scriptFile "#SBATCH --output=".$newJob->{'logFile'}."\n";
		# Find the working directory.
		print $scriptFile "if [ ! -z \${SLURM_SUBMIT_DIR+x} ]; then\n";
		print $scriptFile " cd \$SLURM_SUBMIT_DIR\n";
		print $scriptFile "fi\n";
		print $scriptFile "export "     .$_."\n"
		    foreach ( &List::ExtraUtils::as_array($slurmConfig->{'environment'}) );
		print $scriptFile "module load ".$_."\n"
		    foreach ( &List::ExtraUtils::as_array($slurmConfig->{'module'     }) );
		print $scriptFile "ulimit -t unlimited\n";
		print $scriptFile "ulimit -c unlimited\n";
		print $scriptFile "export OMP_NUM_THREADS=".$ompThreads."\n";
		print $scriptFile (exists($newJob->{'mpi'}) && $newJob->{'mpi'} eq "yes" ? "mpirun -np ".($ppn*$nodes)." " : "").$newJob->{'command'}."\n";
		print $scriptFile "exit\n";
		close($scriptFile);
	    } else {
		# Batch script is already created - we're given its name.
		$batchScript = $newJob;
		undef($newJob);
		$newJob->{'launchFile'} = $batchScript;
	    }
	    # Submit the SLURM job.
	    my $sbatchStatus = undef();
	    my $sbatch;
	    while ( ! defined($sbatchStatus) ) {
		$sbatchStatus = open($sbatch,"sbatch ".$newJob->{'launchFile'}."|");
		unless ( defined($sbatchStatus) ) {
		    print "Could not open pipe to 'sbatch' - will sleep and then try again\n";
		    sleep $waitSleepDuration;
		}
	    }
	    my $jobID = "";
	    while ( my $line = <$sbatch> ) {
	    	if ( $line =~ m/^Submitted batch job (\d+)/ ) {$jobID = $1};
	    }
	    close($sbatch);	    
	    # Add the job number to the active job hash.
	    $slurmJobs{$jobID} = $newJob
		unless ( $jobID eq "" );
	    sleep $submitSleepDuration;
	} else {
	    # Wait.
	    sleep $waitSleepDuration;
	}
    }
}

1;
