# Launch models on PBS system using a single job.

package MonolithicPBS;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use Data::Dumper;
use Sys::CPU;
require Galacticus::Launch::Hooks;
require Galacticus::Launch::PostProcess;
require Galacticus::Launch::PBS;
require System::Redirect;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     monolithicPBS => {
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
	 mpiRun               => "mpirun"             ,
	 nodes                => 1                    ,
	 threadsPerNode       => Sys::CPU::cpu_count(),
	 ompThreads           => Sys::CPU::cpu_count(),
	 analyze              => "yes"                ,
	 jobWaitSleepDuration => 60
	);
    # Attempt to detect MPI implementation.
    my $mpiIs;
    if ( exists($launchScript->{'monolithicPBS'}->{'mpiImplementation'}) ) {
	$mpiIs = $launchScript->{'monolithicPBS'}->{'mpiImplementation'};
    } else {
	$mpiIs = &PBS::mpiDetect();
    }
    if ( $mpiIs eq "OpenMPI" ) {
	$defaults{'mpiRun'} = "mpirun --bynode";
    } elsif ( $mpiIs eq "SGI MPT" ) {
	$defaults{'mpiRun'} = "mpiexec omplace";
    }    
    # Apply defaults.
    foreach ( keys(%defaults) ) {
	$launchScript->{'monolithicPBS'}->{$_} = $defaults{$_}
	unless ( exists($launchScript->{'monolithicPBS'}->{$_}) );
    }
    # Verify that number of OpenMP threads is a factor of the number of threads per node.
    die("MonolithicPBS::Validate: number of OpenMP threads must be a factor of the number of threads per node")
	unless ( $launchScript->{'monolithicPBS'}->{'threadsPerNode'} % $launchScript->{'monolithicPBS'}->{'ompThreads'} == 0 );
    # Determine the number of MPI threads to launch.
    $launchScript->{'monolithicPBS'}->{'mpiThreads'} = 
	$launchScript->{'monolithicPBS'}->{'nodes'         }
       *$launchScript->{'monolithicPBS'}->{'threadsPerNode'}
       /$launchScript->{'monolithicPBS'}->{'ompThreads'    };
    # Attempt to compile a simple code which will return MPI rank.
    system("mkdir -p ".$launchScript->{'modelRootDirectory'});
    open(my $rankCode,">".$launchScript->{'modelRootDirectory'}."/mpiRank.c");
    print $rankCode <<CODE;
#include <stdio.h>
#include <mpi.h>
int main (argc, argv)
     int argc;
     char *argv[];
{
  int rank;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  printf("%d\\n",rank);
  MPI_Finalize();
  return 0;
}
CODE
    close($rankCode);
    system("mpicc ".$launchScript->{'modelRootDirectory'}."/mpiRank.c -o ".$launchScript->{'modelRootDirectory'}."/mpiRank");
    die("MonolithicPBS::Validate: failed to compile MPI rank code")
        unless ( $? == 0 );
}

sub Output_File_Name {
    # Return the output file name for Galacticus models.
    my $outputFileName = shift();
    my $launchScript   = shift();
    if ( exists($launchScript->{'monolithicPBS'}->{'scratchPath'}) ) {
	$outputFileName = 
	    $launchScript->{'monolithicPBS'}->{'scratchPath'}.
	    "/model_".$launchScript->{'modelCounter'}."_".$$."/galacticus.hdf5";
    }
    return $outputFileName;
}

sub Launch {
    # Launch models on local machine.
    my @jobs         = @{shift()};
    my $launchScript =   shift() ;
    # Create the monolithic PBS job script.
    my $singleJobScript;
    $singleJobScript .= "#!/bin/bash\n";
    $singleJobScript .=  "#PBS -N Galacticus\n";
    if ( exists($launchScript->{'config'}->{'contact'}->{'email'}) ) {
	if ( $launchScript->{'config'}->{'contact'}->{'email'} =~ m/\@/ && $launchScript->{'emailReport'} eq "yes" ) {
	    $singleJobScript .=  "#PBS -M ".$launchScript->{'config'}->{'contact'}->{'email'}."\n";
	    $singleJobScript .=  "#PBS -m bea\n";
	}
    }
    $singleJobScript .=  "#PBS -l walltime=".$launchScript->{'monolithicPBS'}->{'wallTime'}."\n"
	if ( exists($launchScript->{'monolithicPBS'}->{'wallTime'}) );
    $singleJobScript .=  "#PBS -l mem=".$launchScript->{'monolithicPBS'}->{'memory'}."\n"
	if ( exists($launchScript->{'monolithicPBS'}->{'memory'}) );
    $singleJobScript .=  "#PBS -l nodes=".$launchScript->{'monolithicPBS'}->{'nodes'}.":ppn=".$launchScript->{'monolithicPBS'}->{'threadsPerNode'}."\n";
    $singleJobScript .=  "#PBS -j oe\n";
    $singleJobScript .=  "#PBS -o ".$launchScript->{'modelRootDirectory'}."/galacticus.log\n";
    $singleJobScript .=  "#PBS -q ".$launchScript->{'monolithicPBS'}->{'queue'}."\n"
	if ( exists($launchScript->{'monolithicPBS'}->{'queue'}) );
    $singleJobScript .= "#PBS -V\n";
    $singleJobScript .= "cd \$PBS_O_WORKDIR\n";
    $singleJobScript .= "chmod u=wrx ".$launchScript->{'modelRootDirectory'}."/launchGalacticus.sh\n";
    $singleJobScript .= $launchScript->{'monolithicPBS'}->{'mpiRun'}." -np ".$launchScript->{'monolithicPBS'}->{'mpiThreads'}." ".$launchScript->{'modelRootDirectory'}."/launchGalacticus.sh\n";
    my $pbsScriptFile = $launchScript->{'modelRootDirectory'}."/pbsLaunch.sh";
    open(my $pbsHndl,">".$pbsScriptFile);
    print $pbsHndl $singleJobScript;
    close($pbsHndl);
    # Iterate over jobs.
    my $singleJobRank      = -1;
    my $singleLaunchScript = "#!/usr/bin/env sh\n";
    $singleLaunchScript   .= "MPI_RANK=`".$launchScript->{'modelRootDirectory'}."/mpiRank`\n";
    foreach my $job ( @jobs ) {
	# Create the PBS submission script.
	my $launchFileName = $job->{'directory'}."/launch.sh";
	open(my $launchFile,">".$launchFileName);
	print $launchFile "#!/bin/bash\n";
	if ( exists($launchScript->{'monolithicPBS'}->{'environment'}) ) {
	    foreach my $environment ( @{$launchScript->{'monolithicPBS'}->{'environment'}} ) {
		print $launchFile "export ".$environment."\n";
	    }
	}
	print $launchFile "export GFORTRAN_ERROR_DUMPCORE=YES\n";
	print $launchFile "ulimit -t unlimited\n";
	print $launchFile "ulimit -c unlimited\n";
	print $launchFile "export OMP_NUM_THREADS=".$launchScript->{'monolithicPBS'}->{'ompThreads'}."\n";
	if ( exists($launchScript->{'monolithicPBS'}->{'scratchPath'}) ) {
	    my $scratchPath = $launchScript->{'monolithicPBS'}->{'scratchPath'}."/model_".$launchScript->{'modelCounter'}."_".$$."/";
	    print $launchFile "mkdir -p ".$scratchPath."\n";
	}
	print $launchFile $galacticusPath."/Galacticus.exe ".$job->{'directory'}."/parameters.xml\n";
	if ( exists($launchScript->{'monolithicPBS'}->{'scratchPath'}) ) {
	    print $launchFile "mv ".$launchScript->{'monolithicPBS'}->{'scratchPath'}."galacticus.hdf5 ".$job->{'directory'}."/galacticus.hdf5\n";
	    if ( $launchScript->{'useStateFile'} eq "yes" ) {
		print $launchFile "mv ".$launchScript->{'monolithicPBS'}->{'scratchPath'}."galacticus_".$launchScript->{'modelCounter'}."_".$$.".state* ".$job->{'directory'}."/\n";
		print $launchFile "mv ".$launchScript->{'monolithicPBS'}->{'scratchPath'}."galacticus_".$launchScript->{'modelCounter'}."_".$$.".fgsl.state* ".$job->{'directory'}."/\n";
	    }
	}
	print $launchFile "mv core* ".$job->{'directory'}."/\n";
	print $launchFile $job->{'analysis'}
	    if ( $launchScript->{'monolithicPBS'}->{'analyze'} eq "yes" );
	close($launchFile);
	# Enter the model in the queue to be run.
	$singleJobRank       = ($singleJobRank+1) % $launchScript->{'monolithicPBS'}->{'mpiThreads'};
	$singleLaunchScript .= "if [ \$MPI_RANK -eq ".$singleJobRank." ]; then\n";
	$singleLaunchScript .= " source ".$launchFileName."\n";
	$singleLaunchScript .= "fi\n";
    }
    $singleLaunchScript .= "exit\n";
    open(my $launchHndl,">".$launchScript->{'modelRootDirectory'}."/launchGalacticus.sh");
    print $launchHndl $singleLaunchScript;
    close($launchHndl);
    # Submit jobs.
    print " -> submitting script: ".$pbsScriptFile."\n";
    open(my $pHndl,"qsub ".$pbsScriptFile."|");
    my $jobID;
    while ( my $line = <$pHndl> ) {
	if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
    }
    close($pHndl);
    die ("MonolithicPBS: failed to capture job ID")
	unless ( defined($jobID) );
    print " -> waiting for PBS jobs to finish...\n"
	if ( $launchScript->{'verbosity'} > 0 );
    my $jobFinished = 0;
    while ( $jobFinished == 0 ) {
	# Find all PBS jobs that are running.
	$jobFinished  = 1;
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^Job\sId:\s+(\S+)/ ) {
		my $runningJobID = $1;
		$jobFinished = 0
		    if ( $runningJobID eq $jobID );
	    };
	}
	close(pHndl);
	sleep $launchScript->{'monolithicPBS'}->{'jobWaitSleepDuration'}
	    unless ( $jobFinished == 1 );
    }
    # Post-process.
    foreach my $job ( @jobs ) {
	&PostProcess::Analyze($job,$launchScript)
	    unless ( $launchScript->{'monolithicPBS'}->{'analyze'} eq "yes" );
	&PostProcess::CleanUp($job,$launchScript);
    }
}

1;
