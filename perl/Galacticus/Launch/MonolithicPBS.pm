# Launch models on PBS system using a single job.

package Galacticus::Launch::MonolithicPBS;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use System::CPU;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PostProcess;
use Galacticus::Launch::PBS;
use System::Redirect;
use List::ExtraUtils;

# Insert hooks for our functions.
%Galacticus::Launch::Hooks::moduleHooks = 
    (
     %Galacticus::Launch::Hooks::moduleHooks,
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
	 mpiRun               => "mpirun"               ,
	 mpiOptions           => ""                     ,
	 nodes                => 1                      ,
	 threadsPerNode       => System::CPU::get_ncpu(),
	 ompThreads           => System::CPU::get_ncpu(),
	 shell                => "bash"                 ,
	 analyze              => "yes"                  ,
	 jobWaitSleepDuration => 60
	);
    # Attempt to detect MPI implementation.
    my $mpiIs;
    if ( exists($launchScript->{'monolithicPBS'}->{'mpiImplementation'}) ) {
	$mpiIs = $launchScript->{'monolithicPBS'}->{'mpiImplementation'};
    } else {
	$mpiIs = &Galacticus::Launch::PBS::mpiDetect();
    }
    if ( $mpiIs eq "OpenMPI" ) {
	    $defaults{'mpiRun'         } = "mpirun";
	    $defaults{'mpiOptions'     } = "--bynode";
	    $defaults{'resourceRequest'} = "#PBS -l nodes=%%NODES%%:ppn=%%THREADS%%";
    } elsif ( $mpiIs eq "SGI MPT" ) {
	    $defaults{'mpiRun'         } = "mpiexec";
	    $defaults{'mpiOptions'     } = "omplace";
	    $defaults{'resourceRequest'} = "#PBS -l select=%%NODES%%:ncpus=%%THREADS%%:mpiprocs=%%THREADS%%";
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
  int rank, size;
  char command[1024];
  if ( argc != 2 )
    {
      printf("usage: %s <executable>\\n",argv[0]);
      return 1;
    }
  else
    {
      MPI_Init (&argc, &argv);
      MPI_Comm_rank (MPI_COMM_WORLD, &rank);
      MPI_Comm_size (MPI_COMM_WORLD, &size);
      sprintf(command,"%s %d %d",argv[1],rank,size);
      printf("mpiRank running: %s\\n",command);
      system(command);
      MPI_Finalize();
      return 0;
    }
}
CODE
    close($rankCode);
    my $compileCommand;
    $compileCommand .= &Set_Environment($launchScript);
    $compileCommand .= "mpicc ".$launchScript->{'modelRootDirectory'}."/mpiRank.c -o ".$launchScript->{'modelRootDirectory'}."/mpiRank";
    $compileCommand .= " ".join(" ",map {"-I".$_} &List::ExtraUtils::as_array($launchScript->{'monolithicPBS'}->{'includePath'}))
        if ( exists($launchScript->{'monolithicPBS'}->{'includePath'}) );
    $compileCommand .= " ".join(" ",map {"-L".$_} &List::ExtraUtils::as_array($launchScript->{'monolithicPBS'}->{'libraryPath'}))
        if ( exists($launchScript->{'monolithicPBS'}->{'libraryPath'}) );
    system($compileCommand);
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
    $singleJobScript .= "#!/usr/bin/env ".$launchScript->{'monolithicPBS'}->{'shell'}."\n";
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
    my $resources = $launchScript->{'monolithicPBS'}->{'resourceRequest'};
    $resources =~ s/\%\%NODES\%\%/$launchScript->{'monolithicPBS'}->{'nodes'}/g;
    $resources =~ s/\%\%THREADS\%\%/$launchScript->{'monolithicPBS'}->{'threadsPerNode'}/g;
    $singleJobScript .= $resources."\n";
    $singleJobScript .=  "#PBS -j oe\n";
    $singleJobScript .=  "#PBS -o ".$launchScript->{'modelRootDirectory'}."/galacticus.log\n";
    $singleJobScript .=  "#PBS -q ".$launchScript->{'monolithicPBS'}->{'queue'}."\n"
	if ( exists($launchScript->{'monolithicPBS'}->{'queue'}) );
    $singleJobScript .= "#PBS -V\n";
    $singleJobScript .= "cd \$PBS_O_WORKDIR\n";
    $singleJobScript .= "chmod u=wrx ".$launchScript->{'modelRootDirectory'}."/launchGalacticus.sh\n";
    $singleJobScript .= join("\n",&List::ExtraUtils::as_array($launchScript->{'monolithicPBS'}->{'pbsCommand'}))."\n"
	if ( exists($launchScript->{'monolithicPBS'}->{'pbsCommand'}) );
    $singleJobScript .= "setenv MPI_DSM_DISTRIBUTE 0\n";
    $singleJobScript .= "setenv KMP_AFFINITY disabled\n";
    $singleJobScript .= "setenv OMP_NUM_THREADS ".$launchScript->{'monolithicPBS'}->{'ompThreads'}."\n";
    $singleJobScript .= $launchScript->{'monolithicPBS'}->{'mpiRun'}." -np ".$launchScript->{'monolithicPBS'}->{'mpiThreads'}." ".$launchScript->{'monolithicPBS'}->{'mpiOptions'}." ".$launchScript->{'modelRootDirectory'}."/mpiRank ".$launchScript->{'modelRootDirectory'}."/launchGalacticus.sh\n";
    my $pbsScriptFile = $launchScript->{'modelRootDirectory'}."/pbsLaunch.sh";
    open(my $pbsHndl,">".$pbsScriptFile);
    print $pbsHndl $singleJobScript;
    close($pbsHndl);
    # Iterate over jobs.
    my $singleJobRank      = -1;
    my $singleLaunchScript = "#!/usr/bin/env ".$launchScript->{'monolithicPBS'}->{'shell'}."\n";
    $singleLaunchScript   .= "MPI_RANK=\$1\n";
    foreach my $job ( @jobs ) {
	# Create the PBS submission script.
	my $launchFileName = $job->{'directory'}."/launch.sh";
	open(my $launchFile,">".$launchFileName);
	print $launchFile "#!/usr/bin/env ".$launchScript->{'monolithicPBS'}->{'shell'}."\n";
	print $launchFile &Set_Environment($launchScript);
	print $launchFile "export GFORTRAN_ERROR_DUMPCORE=YES\n";
	print $launchFile "ulimit -t unlimited\n";
	print $launchFile "ulimit -c unlimited\n";
	print $launchFile "export OMP_NUM_THREADS=".$launchScript->{'monolithicPBS'}->{'ompThreads'}."\n";
	if ( exists($launchScript->{'monolithicPBS'}->{'scratchPath'}) ) {
	    my $scratchPath = $launchScript->{'monolithicPBS'}->{'scratchPath'}."/model_".$launchScript->{'modelCounter'}."_".$$."/";
	    print $launchFile "mkdir -p ".$scratchPath."\n";
	}
	print $launchFile $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$job->{'directory'}."/parameters.xml\n";
	if ( exists($launchScript->{'monolithicPBS'}->{'scratchPath'}) ) {
	    print $launchFile "mv ".$launchScript->{'monolithicPBS'}->{'scratchPath'}."galacticus.hdf5 ".$job->{'directory'}."/galacticus.hdf5\n";
	    if ( $launchScript->{'useStateFile'} eq "yes" ) {
		print $launchFile "mv ".$launchScript->{'monolithicPBS'}->{'scratchPath'}."galacticus_".$launchScript->{'modelCounter'}."_".$$.".state* ".$job->{'directory'}."/\n";
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
	&Galacticus::Launch::PostProcess::Analyze($job,$launchScript)
	    unless ( $launchScript->{'monolithicPBS'}->{'analyze'} eq "yes" );
	&Galacticus::Launch::PostProcess::CleanUp($job,$launchScript);
    }
}

sub Set_Environment {
    my $launchScript = shift;
    my $environment = "";
    if ( exists($launchScript->{'monolithicPBS'}->{'environment'}) ) {
	if ( $launchScript->{'monolithicPBS'}->{'shell'} eq "bash" ) {
	    $environment = join("\n",map {"export ".$_} &List::ExtraUtils::as_array($launchScript->{'monolithicPBS'}->{'environment'}))."\n";
	} elsif ( $launchScript->{'monolithicPBS'}->{'shell'} eq "csh" ) {
	    $environment = join("\n",map {local $_ = $_; s/(.*)=(.*)/setenv $1 $2/; $_} &List::ExtraUtils::as_array($launchScript->{'monolithicPBS'}->{'environment'}))."\n";
	}
    }
    return $environment;
}

1;
