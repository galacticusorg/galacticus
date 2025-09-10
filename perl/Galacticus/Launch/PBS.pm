# Launch models on PBS system.

package Galacticus::Launch::PBS;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use File::Which;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PostProcess;
use Galacticus::Options;
use System::Redirect;
use List::ExtraUtils;

# Insert hooks for our functions.
%Galacticus::Launch::Hooks::moduleHooks = 
    (
     %Galacticus::Launch::Hooks::moduleHooks,
     pbs => {
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
	 mpiLaunch               => "yes"                                                                     ,
	 mpiRun                  => "mpirun --map-by node --mca mpi_preconnect_mpi 1 -hostfile \$PBS_NODEFILE",
	 mpiProcesses            => 1                                                                         ,
	 maxJobsInQueue          => -1                                                                        ,
	 postSubmitSleepDuration => 10                                                                        ,
	 jobWaitSleepDuration    => 60                                                                        ,
	 analyze                 => "yes"
	);
    # Attempt to detect MPI implementation.
    my $mpiIs = &mpiDetect();
    if ( $mpiIs eq "OpenMPI" ) {
	$defaults{'mpiLaunch'} = "yes"                                                                     ;
	$defaults{'mpiRun'   } = "mpirun --map-by node --mca mpi_preconnect_mpi 1 -hostfile \$PBS_NODEFILE";
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
    my %options      = %{shift()}
        if ( scalar(@_) > 0 );
    # Find localized configuration.
    my $pbsConfig = &Galacticus::Options::Config("pbs");
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
	my $nodes;
	my $ppn  ;
	if ( $launchScript->{'pbs'}->{'mpiLaunch'} eq "yes" ) {
	    if ( exists($launchScript->{'pbs'}->{'mpiNodes'}) ) {
		$nodes = $launchScript->{'pbs'}->{'mpiNodes'};
		$ppn   = $launchScript->{'pbs'}->{'mpiProcesses'}*$launchScript->{'pbs'}->{'ompThreads'}/$launchScript->{'pbs'}->{'mpiNodes'};
	    } else {
		$nodes = 1;
		$ppn   = $launchScript->{'pbs'}->{'mpiProcesses'}*$launchScript->{'pbs'}->{'ompThreads'}
	    }
	} else {
	    $nodes = 1;
	    $ppn   = $launchScript->{'pbs'}->{'ompThreads'};
	}
	print $pbsFile "#PBS -l nodes=".$nodes.":ppn=".$ppn."\n";
	print $pbsFile "#PBS -j oe\n";
	my $pwd = "";
	unless ( $job->{'directory'} =~ m/^\// ) {
	    $pwd = `pwd`;
	    chomp($pwd);
	    $pwd .= "/";
	}
	print $pbsFile "#PBS -o ".$pwd.$job->{'directory'}."/galacticus.log\n";
	if ( exists($launchScript->{'pbs'}->{'queue'}) ) {
	    print $pbsFile "#PBS -q ".$launchScript->{'pbs'}->{'queue'}."\n";
	} elsif ( exists($pbsConfig->{'queue'}) ) {
	    print $pbsFile "#PBS -q ".$pbsConfig->{'queue'}."\n";
	}
	print $pbsFile "#PBS -V\n";
	# Find the working directory - we support either PBS or SLURM environment variables here.
	print $pbsFile "if [ ! -z \${PBS_O_WORKDIR+x} ]; then\n";
	print $pbsFile " cd \"\${PBS_O_WORKDIR}\"\n";
	print $pbsFile "elif [ ! -z \${SLURM_SUBMIT_DIR+x} ]; then\n";
	print $pbsFile " cd \"\${SLURM_SUBMIT_DIR}\"\n";
	print $pbsFile "fi\n";
	print $pbsFile "export ".$_."\n"
	    foreach ( &List::ExtraUtils::as_array($pbsConfig->{'environment'}) );
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
	    foreach my $command ( &List::ExtraUtils::as_array($launchScript->{'pbs'}->{'preCommand'}) ) {
		print $pbsFile $command."\n";
	    }
	}
	print $pbsFile $launchScript->{'pbs'}->{'mpiRun'}." -np ".(exists($launchScript->{'pbs'}->{'mpiProcesses'}) ? $launchScript->{'pbs'}->{'mpiProcesses'} : "1")." "
	    if ( $launchScript->{'pbs'}->{'mpiLaunch'} eq "yes" );
	print $pbsFile $ENV{'GALACTICUS_EXEC_PATH'}."/".($launchScript->{'pbs'}->{'executable'} ? $launchScript->{'pbs'}->{'executable'} : "Galacticus.exe")." ".$job->{'directory'}."/parameters.xml\n";
	if ( exists($launchScript->{'pbs'}->{'scratchPath'}) ) {
	    print $pbsFile "mv ".$launchScript->{'pbs'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus.hdf5 ".$job->{'directory'}."/galacticus.hdf5\n";
	    if ( $launchScript->{'useStateFile'} eq "yes" ) {
		print $pbsFile "mv ".$launchScript->{'pbs'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus_".$job->{'modelCounter'}."_".$$.".state* ".$job->{'directory'}."/\n";
	    }
	}
	print $pbsFile "if compgen -G \"core.*\" > /dev/null; then\n";
	print $pbsFile " mv core* ".$job->{'directory'}."/\n";
	print $pbsFile "fi\n";
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
    # Determine the maximum number of jobs in the queue.
    my $jobMaximum = $launchScript->{'pbs'}->{'maxJobsInQueue'};
    $jobMaximum    = $options{'pbsJobMaximum'}
        if ( exists($options{'pbsJobMaximum'}) );
    # Determine sleep durations.
    my $postSubmitSleepDuration = $launchScript->{'pbs'}->{'postSubmitSleepDuration'};
    $postSubmitSleepDuration = $options{'submitSleepDuration'}
        if ( exists($options{'submitSleepDuration'}) );
    my $jobWaitSleepDuration    = $launchScript->{'pbs'}->{'jobWaitSleepDuration'   };
    $jobWaitSleepDuration    = $options{'waitSleepDuration'  }
        if ( exists($options{'waitSleepDuration'  }) );    
    # Submit jobs.
    print " -> waiting for PBS jobs to finish...\n"
	if ( $launchScript->{'verbosity'} > 0 );
    my %pbsJobs;
    while ( scalar(keys %pbsJobs) > 0 || scalar(@modelQueue) > 0 ) {
	# Find all PBS jobs that are running.
	my %runningPBSJobs;
	undef(%runningPBSJobs);
	my $jobID;
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^Job\sId:\s+(\S+)/ ) {$jobID = $1};
	    if ( $line =~ m/^\s*job_state\s*=\s*[^C]\s*$/ ) {$runningPBSJobs{$jobID} = 1};
	}
	close(pHndl);
	my $jobFinished = 0;
	foreach my $jobID ( keys(%pbsJobs) ) {
	    unless ( exists($runningPBSJobs{$jobID}) ) {
		print " -> PBS job ".$jobID." has finished; post-processing....\n";
		$jobFinished = 1;
		# Extract job timing information, and report on it.
		my $exitStatus = 0;
		my $startAt    = "unknown";
		my $stopAt     = "unknown";
		my $traceJob = &File::Which::which("tracejob");
		my $sacct    = &File::Which::which("sacct"   );
		if ( defined($traceJob) ) {
		    open(my $trace,$traceJob." -q -n 100 ".$jobID." |");
		    while ( my $line = <$trace> ) {
			if ( $line =~ m/Exit_status=(\d+)/ ) {
			    $exitStatus = $1;
			}
			if ( $line =~ m/^([0-9\/]+\s[0-9:]+)\s+S\s+Job Run/ ) {
			    $startAt = $1;
			}
			if ( $line =~ m/^([0-9\/]+\s[0-9:]+)\s+S\s+Exit_status/ ) {
			    $stopAt  = $1;
			}
		    }
		    close($trace);
		} elsif ( defined($sacct) ) {
		    open(my $trace,$sacct." -b -j ".$jobID." |");
		    while ( my $line = <$trace> ) {
			my @columns = split(" ",$line);
			if ( $columns[0] eq $jobID ) {
			    if ( $columns[2] =~ m/(\d+):\d+/ ) {
				$exitStatus = $1;
			    }
			}
		    }
		    close($trace);
		}
		print "  -> Job metadata: Exit status : Start time : End time = ".$pbsJobs{$jobID}->{'job'}->{'directory'}." : ".$exitStatus." : ".$startAt." : ".$stopAt."\n";
		# Perform analysis.
		&Galacticus::Launch::PostProcess::Analyze($pbsJobs{$jobID}->{'job'},$launchScript)
		    unless ( $launchScript->{'pbs'}->{'analyze'} eq "yes" );
		# Clean up.
		&Galacticus::Launch::PostProcess::CleanUp($pbsJobs{$jobID}->{'job'},$launchScript);
		# Remove the job ID from the list of active PBS jobs.
		delete($pbsJobs{$jobID});
	    }
	}
	# If there are jobs remaining to be submitted, and not too many are already queued, launch one.
	if ( 
	    scalar(@modelQueue) > 0 
	    &&
	    (
	     scalar(keys %pbsJobs) < $jobMaximum
	     ||
	     $jobMaximum < 0 
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
	    sleep $postSubmitSleepDuration;
	} else {
	    sleep $jobWaitSleepDuration
		unless ( $jobFinished == 1 );
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
    # Submit jobs to PBS and wait for them to finish.
    my %arguments = %{shift()};
    my @pbsStack  = @_;
    my %pbsJobs;
    # Find the appropriate PBS section.
    my $pbsConfig = &Galacticus::Options::Config("pbs");
    # Determine sleep times between jobs.
    my $submitSleepDuration = exists($arguments{'submitSleepDuration'}) ? $arguments{'submitSleepDuration'} :  5;
    my $waitSleepDuration   = exists($arguments{'waitSleepDuration'  }) ? $arguments{'waitSleepDuration'  } : 60;
    # Determine maximum number allowed in queue at once.
    my $jobMaximum = 10;
    $jobMaximum = $pbsConfig->{'jobMaximum'}
       if ( exists($pbsConfig->{'jobMaximum'}) );
    $jobMaximum = $arguments{'pbsJobMaximum'}
       if ( exists($arguments{'pbsJobMaximum'}) );
    my $limitQueuedOnly = exists($arguments{'limitQueuedOnly'}) ? $arguments{'limitQueuedOnly'} : 0;
    # Submit jobs and wait.
    print "Waiting for PBS jobs to finish...\n";
    while ( scalar(keys %pbsJobs) > 0 || scalar(@pbsStack) > 0 ) {
	# Find all PBS jobs that are running.
	my %runningPBSJobs;
	undef(%runningPBSJobs);
	my $jobID;
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^(Job\sId|Job):\s+(\S+)/ ) {$jobID = $2};
	    if ( $line =~ m/job_state\s*=\s*([RQHE])/ && $jobID ) {
		$runningPBSJobs{$jobID} = $1;
		undef($jobID);
	    }
	}
	close(pHndl);
	foreach my $jobID ( keys(%pbsJobs) ) {
	    unless ( exists($runningPBSJobs{$jobID}) ) {
		print "PBS job ".$jobID." has finished.\n";
		# Call any "on completion" function.
		if ( exists($pbsJobs{$jobID}->{'onCompletion'}) ) {
		    my $exitStatus = 0;
		    my $startAt;
		    my $stopAt;
		    if ( exists($pbsJobs{$jobID}->{'tracejob'}) && $pbsJobs{$jobID}->{'tracejob'} eq "yes" ) {
			my $traceJob = &File::Which::which("tracejob");
			my $sacct    = &File::Which::which("sacct"   );
			if ( defined($traceJob) ) {
			    open(my $trace,$traceJob." -q -n 100 ".$jobID." |");
			    while ( my $line = <$trace> ) {
				if ( $line =~ m/Exit_status=(\d+)/ ) {
				    $exitStatus = $1;
				}
				if ( $line =~ m/^([0-9\/]+\s[0-9:]+)\s+S\s+Job Run/ ) {
				    $startAt = $1;
				}
				if ( $line =~ m/^([0-9\/]+\s[0-9:]+)\s+S\s+Exit_status/ ) {
				    $stopAt  = $1;
				}
			    }
			    close($trace);
			} elsif ( defined($sacct) ) {
			    open(my $trace,$sacct." -b -j ".$jobID." |");
			    while ( my $line = <$trace> ) {
				my @columns = split(" ",$line);
				if ( $columns[0] eq $jobID ) {
				    if ( $columns[2] =~ m/(\d+):\d+/ ) {
					$exitStatus = $1;
				    }
				}
			    }
			    close($trace);
			}
		    }
		    &{$pbsJobs{$jobID}->{'onCompletion'}->{'function'}}(@{$pbsJobs{$jobID}->{'onCompletion'}->{'arguments'}},$jobID,$exitStatus,\@pbsStack,$pbsJobs{$jobID},$startAt,$stopAt);
		}
		# Remove the job ID from the list of active PBS jobs.
		delete($pbsJobs{$jobID});		
	    }
	}
	# If fewer than maximum number of jobs are in the queue, pop one off the stack.
	my $limitingJobs = $limitQueuedOnly ? scalar(grep {($runningPBSJobs{$_} eq "Q" || $runningPBSJobs{$_} eq "H") && exists($pbsJobs{$_})} keys(%runningPBSJobs)) : scalar(keys(%pbsJobs));
	if ( scalar(@pbsStack) > 0 && $limitingJobs < $jobMaximum ) {
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
		my $resourceModel = exists($newJob->{'resourceModel'}) ? $newJob->{'resourceModel'} : "nodes";
		my $nodes    = 1;
		my $ppn      = 1;
		$ppn   = $pbsConfig->{'ppn'  }
		    if ( exists($pbsConfig->{'ppn'  }) );
		$ppn   = $arguments  {'ppn'  }
		    if ( exists($arguments  {'ppn'  }) );
		$nodes = $arguments  {'nodes'}
		    if ( exists($arguments{'nodes'}) );
		$ppn   = $newJob   ->{'ppn'  }
		    if ( exists($newJob   ->{'ppn'  }) );
		$nodes = $newJob   ->{'nodes'}
		    if ( exists($newJob   ->{'nodes'}) );
		my $mpiProcs = $nodes*$ppn;
		$mpiProcs = $arguments{'mpiProcs'}
		    if ( exists($arguments{'mpiProcs'}) );
		$mpiProcs = $newJob->{'mpiProcs'}
		    if ( exists($newJob->{'mpiProcs'}) );
		open(my $scriptFile,">".$newJob->{'launchFile'});
		print $scriptFile "#!/bin/bash\n";
		print $scriptFile "#PBS -N ".$newJob->{'label'}."\n";
		print $scriptFile "#PBS -l walltime=".$newJob->{'walltime'}."\n"
		    if ( exists($newJob->{'walltime'}) );
		if ( exists($arguments{'queue'}) ) {
		    print $scriptFile "#PBS -q ".$arguments{'queue'}."\n";
		} elsif ( exists($pbsConfig->{'queue'}) ) {
		    print $scriptFile "#PBS -q ".$pbsConfig->{'queue'}."\n";
		}
		if ( $resourceModel eq "nodes" ) {
		    print $scriptFile "#PBS -l nodes=".$nodes.":ppn=".$ppn.(exists($newJob->{'model'}) && defined($newJob->{'model'}) ? ":model=".$newJob->{'model'} : "")."\n";
		    print $scriptFile "#PBS -l mem=".$newJob->{'mem'}."\n"
			if ( exists($newJob->{'mem'}) );
		} elsif ( $resourceModel eq "select" ) {
		    print $scriptFile "#PBS -l select=".$nodes.":ncpus=".$ppn.(exists($newJob->{'model'}) && defined($newJob->{'model'}) ? ":model=".$newJob->{'model'} : "")."\n";
		} else {
		    die("Galacticus::Launch::PBS::SubmitJobs: unknown resource model");
		}
		print $scriptFile "#PBS -j oe\n";
		(my $logFileName = $newJob->{'logFile'}) =~ s/:/_/g;
		print $scriptFile "#PBS -o ".$logFileName."\n";
		print $scriptFile "#PBS -V\n";
		# Find the working directory - we support either PBS or SLURM environment variables here.
		print $scriptFile "if [ ! -z \${PBS_O_WORKDIR+x} ]; then\n";
		print $scriptFile " cd \"\${PBS_O_WORKDIR}\"\n";
		print $scriptFile "elif [ ! -z \${SLURM_SUBMIT_DIR+x} ]; then\n";
		print $scriptFile " cd \"\${SLURM_SUBMIT_DIR}\"\n";
		print $scriptFile "fi\n";
		print $scriptFile "export ".$_."\n"
		    foreach ( &List::ExtraUtils::as_array($pbsConfig->{'environment'}) );
		print $scriptFile "ulimit -c unlimited\n";
		my $coreDump = "no";
		$coreDump = $arguments{'coreDump'}
		    if ( exists($arguments{'coreDump'}) );
		$coreDump = $newJob->{'coreDump'}
		    if ( exists($newJob->{'coreDump'}) );
		if ( $coreDump eq "yes" ) {
		    print $scriptFile "ulimit -c unlimited\n";
		    print $scriptFile "export GFORTRAN_ERROR_DUMPCORE=YES\n";
		} else {
		    print $scriptFile "ulimit -c 0\n";
		    print $scriptFile "export GFORTRAN_ERROR_DUMPCORE=NO\n";
		}
		print $scriptFile "ulimit -t unlimited\n";
		my $mpi = (exists($arguments{'mpi'}) && $arguments{'mpi'} eq "yes") || (exists($newJob->{'mpi'}) && $newJob->{'mpi'} eq "yes");
		print $scriptFile "export OMP_NUM_THREADS=".($mpi ? 1 : $ppn)."\n";
		print $scriptFile ($mpi ? "mpirun --map-by node --mca mpi_preconnect_mpi 1 --mca btl ^ib -np ".$mpiProcs." " : "").$newJob->{'command'}."\n";
		print $scriptFile "exit\n";
		close($scriptFile);
	    } else {
		# Batch script is already created - we're given its name.
		$batchScript = $newJob;
		undef($newJob);
		$newJob->{'launchFile'} = $batchScript;
	    }	    
	    # Submit the PBS job.
	    open(pHndl,"qsub ".$newJob->{'launchFile'}."|");
	    my $jobID = "";
	    while ( my $line = <pHndl> ) {
	    	if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
	    }
	    close(pHndl);	    
	    # Add the job number to the active job hash.
	    $pbsJobs{$jobID} = $newJob
		unless ( $jobID eq "" );
	    sleep $submitSleepDuration;
	} else {
	    # Wait.
	    sleep $waitSleepDuration;
	}
    }
}

1;
