# Launch models on SLURM system.

package Slurm;
use strict;
use warnings;
use Data::Dumper;
use Sys::CPU;
use Cwd;
require File::Which;
require Galacticus::Launch::Hooks;
require Galacticus::Launch::PostProcess;
require System::Redirect;
require List::ExtraUtils;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     slurm => {
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
	print $slurmFile "#SBATCH --workdir=".$currentDirectory."\n";
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
	    foreach my $command ( &ExtraUtils::as_array($launchScript->{'slurm'}->{'preCommand'}) ) {
		print $slurmFile $command."\n";
	    }
	}
	print $slurmFile $launchScript->{'slurm'}->{'mpiRun'}." -np 1 "
	    if ( $launchScript->{'slurm'}->{'mpiLaunch'} eq "yes" );
	print $slurmFile "./Galacticus.exe ".$job->{'directory'}."/parameters.xml\n";
	if ( exists($launchScript->{'slurm'}->{'scratchPath'}) ) {
	    print $slurmFile "mv ".$launchScript->{'slurm'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus.hdf5 ".$job->{'directory'}."/galacticus.hdf5\n";
	    if ( $launchScript->{'useStateFile'} eq "yes" ) {
		print $slurmFile "mv ".$launchScript->{'slurm'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus_".$job->{'modelCounter'}."_".$$.".state* ".$job->{'directory'}."/\n";
		print $slurmFile "mv ".$launchScript->{'slurm'}->{'scratchPath'}."/model_".$job->{'modelCounter'}."_".$$."/galacticus_".$job->{'modelCounter'}."_".$$.".fgsl.state* ".$job->{'directory'}."/\n";
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
		&PostProcess::Analyze($slurmJobs{$jobID}->{'job'},$launchScript)
		    unless ( $launchScript->{'slurm'}->{'analyze'} eq "yes" );
		&PostProcess::CleanUp($slurmJobs{$jobID}->{'job'},$launchScript);
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

1;
