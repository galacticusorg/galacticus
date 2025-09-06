# Launch models on local system.

package Galacticus::Launch::Local;
use strict;
use warnings;
use Config;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use System::CPU;
use Clone qw(clone);
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PostProcess;
use System::Redirect;

our $useThreads;
BEGIN {
    if ($Config{usethreads}) {
        # We have threads
        $useThreads = 1;
	require threads;
    } else {
	$useThreads = 0;
    }
}

# Insert hooks for our functions.
%Galacticus::Launch::Hooks::moduleHooks = 
    (
     %Galacticus::Launch::Hooks::moduleHooks,
     local => {
	 validate       => \&Validate        ,
	 outputFileName => \&Output_File_Name,
	 launch         => \&Launch          ,
	 jobArrayLaunch => \&jobArrayLaunch
     }
    );

sub Validate {
    # Validate the launch script.
    my $launchScript = shift();
    # Set defaults.
    my %defaults = 
	(
	 threadCount => 1               ,
	 ompThreads  => "maximum"       ,
	 executable  => "Galacticus.exe"
	);
    foreach ( keys(%defaults) ) {
	$launchScript->{'local'}->{$_} = $defaults{$_}
	unless ( exists($launchScript->{'local'}->{$_}) );
    }
    # Determine how many threads to launch.
    $launchScript->{'local'}->{'threadCount'} = System::CPU::get_ncpu() 
	if ( $launchScript->{'local'}->{'threadCount'} eq "maximum" );
    $launchScript->{'local'}->{'ompThreads'} = System::CPU::get_ncpu() 
	if ( $launchScript->{'local'}->{'ompThreads'} eq "maximum" );
}

sub Output_File_Name {
    # Return the output file name for Galacticus models.
    my $outputFileName = shift;
    return $outputFileName;
}

sub Launch {
    # Launch models on local machine.
    my @jobs         = @{shift()};
    my $launchScript =   shift() ;
    my %options      = %{shift()}
        if ( scalar(@_) > 0 );
    if ( $useThreads ) {
	# Determine number of threads to use.
	my $threadCount = $launchScript->{'local'}->{'threadCount'};
	$threadCount = $options{'threadMaximum'}
        if ( exists($options{'threadMaximum'}) );
	# Launch model threads.
	my @threads;
	for(my $iThread=0;$iThread<$threadCount;++$iThread) {
	    print " -> launching thread ".$iThread." of ".$launchScript->{'local'}->{'threadCount'}."\n"
		if ( $launchScript->{'verbosity'} > 0 );
	    push(
		@threads,
		threads->create(\&Launch_Models, $iThread, $threadCount, \@jobs, $launchScript, \%options)
		);
	}
	# Wait for threads to finish.
	foreach my $thread ( @threads ) {
	    my @returnData = $thread->join();
	}
    } else {
	&Launch_Models(0,1,\@jobs,$launchScript,\%options);
    }
}

sub Launch_Models {
    my $iThread      =   shift() ;
    my $threadCount  =   shift() ;
    my @jobs         = @{shift()};
    my $launchScript =   shift() ;
    my %options      = %{shift()}
        if ( scalar(@_) > 0 );
    my $ompThreads = $launchScript->{'local'}->{'ompThreads'};
    $ompThreads = $options{'ompThreads'}
        if ( exists($options{'ompThreads'}) );
    my $verbosity = $launchScript->{'verbosity'};
    $verbosity = $options{'verbosity'}
    if ( exists($options{'verbosity'}) );
    my $executable = $launchScript->{'local'}->{'executable'};
    $executable = $options{'executable'}
        if ( exists($options{'executable'}) );
    for(my $i=0;$i<scalar(@jobs);++$i) {
	if ( ( $i % $threadCount ) == $iThread ) {
 	    print " -> thread ".$iThread." running job: ".$jobs[$i]->{'label'}."\n"
		if ( $verbosity > 0 );
	    &System::Redirect::tofile(
		    "ulimit -t unlimited;"        .
		    "ulimit -c unlimited;"        .
		    "export GFORTRAN_ERROR_DUMPCORE=YES;".
		    "export OMP_NUM_THREADS=".$ompThreads.";".
		    $ENV{'GALACTICUS_EXEC_PATH'}."/".$executable." ".$jobs[$i]->{'directory'}."/parameters.xml",
		    $jobs[$i]->{'directory'}."/galacticus.log"
		);
	    if ( $? == 0 ) {
		&Galacticus::Launch::PostProcess::Analyze($jobs[$i],$launchScript);
	    } else {
		&Galacticus::Launch::PostProcess::Failed ($jobs[$i],$launchScript);
	    }
	    &Galacticus::Launch::PostProcess::CleanUp    ($jobs[$i],$launchScript);
	}
    }
}

sub jobArrayLaunch {
    # Launch an array of jobs as threads on the local machine.
    my %arguments = %{shift()};
    my @jobStack  = @_;
    # Find the appropriate PBS section.
    my $localConfig = &Galacticus::Options::Config("local");
    # Determine maximum number allowed in queue at once.
    my $jobMaximum = System::CPU::get_ncpu() ;
    $jobMaximum = $localConfig->{'jobMaximum'}
       if ( exists($localConfig->{'jobMaximum'}) );
    $jobMaximum = $arguments{'threadMaximum'}
       if ( exists($arguments{'threadMaximum'}) );
    # Submit jobs and wait.
    print "Launching jobs locally...\n";
    my @jobsToRun;
    my %activeThreads;
    my $jobNumber = -1;
    while ( scalar(@jobStack) > 0 || scalar(keys(%activeThreads)) > 0 ) {
	# Record of whether any work was done.
	my $workDone = 0;
	# Launch a job if possible.
	if ( scalar(@jobStack) > 0 && scalar(keys(%activeThreads)) < $jobMaximum ) {
	    # Pop a job from the job stack.
	    my $newJob = pop(@jobStack);
	    # Check for a dependency.
	    if ( exists($newJob->{'dependsOn'}) && ! -e $newJob->{'dependsOn'} ) {
		unshift(@jobStack,$newJob);
		sleep 5;
		last;
	    }
	    # Create the script.
	    open(my $scriptFile,">".$newJob->{'launchFile'});
	    print $scriptFile "#!/bin/bash\n";
	    # Find the working directory.
	    print $scriptFile "cd ".$ENV{'GALACTICUS_EXEC_PATH'}."\n";
	    print $scriptFile "export ".$_."\n"
		foreach ( &List::ExtraUtils::as_array($localConfig->{'environment'}) );
	    print $scriptFile "ulimit -t ".$newJob->{'local'}->{'cpuTimeLimit'}."\n"
		if ( exists($newJob->{'local'}->{'cpuTimeLimit'}) );
	    print $scriptFile "ulimit -c unlimited\n";
	    my $ompThreads;
	    if ( exists($newJob->{'local'}->{'ompThreads'}) ) {
		$ompThreads = $newJob->{'local'}->{'ompThreads'};
	    } elsif ( exists($localConfig->{'ompThreads'}) ) {
		$ompThreads = $localConfig->{'ompThreads'};
	    }
	    $ompThreads = $arguments{'ompThreads'}
	        if ( exists($arguments{'ompThreads'}) );
	    print $scriptFile "export OMP_NUM_THREADS=".$ompThreads."\n"
		if ( defined($ompThreads) );
	    print $scriptFile $newJob->{'command'}."\n";
	    print $scriptFile "exit\n";
	    close($scriptFile);
	    # Launch the thread.
	    ++$jobNumber;
	    print "Launch job ".$jobNumber."\n";
	    if ( exists($arguments{'allowThreads'}) && $arguments{'allowThreads'} == 1 ) {
		&PDL::no_clone_skip_warning();
		if ( $^V lt v5.14.0 ) { 
		    $activeThreads{$jobNumber} = 
		    {
			thread => threads->create( sub {system("chmod u+x ".$newJob->{'launchFile'}."; ".$newJob->{'launchFile'}." >& ".$newJob->{'logFile'}); return $?;}),
			job    => $newJob,
			jobID  => $jobNumber
		    };
		} else {
		    $activeThreads{$jobNumber} = 
		    {
			thread => threads->create({'context' => 'scalar'}, sub {system("chmod u+x ".$newJob->{'launchFile'}."; ".$newJob->{'launchFile'}." >& ".$newJob->{'logFile'}); return $?;}),
			job    => $newJob,
			jobID  => $jobNumber
		    };
		}
	    } else {
		$activeThreads{$jobNumber} = 
		{
		    job    => $newJob,
		    jobID  => $jobNumber
		};
		push(@jobsToRun,$newJob);
		if ( scalar(@jobsToRun) == $jobMaximum || scalar(@jobStack) == 0 ) {
		    # Maximum number of jobs have been accumulated (or there are no more left on the stack). Launch them all and
		    # then clear the list.
		    system($ENV{'GALACTICUS_EXEC_PATH'}."/scripts/aux/localLaunchWrapper.pl ".join(" ",map {$_->{'launchFile'}." ".$_->{'logFile'}} @jobsToRun));
		    @jobsToRun = ();
		}
	    }
	    $workDone = 1;
	} else {
	    # Check for threads which can be rejoined.
	    foreach ( keys(%activeThreads) ) {
		if ( exists($arguments{'allowThreads'}) && $arguments{'allowThreads'} == 1 ) {
		    if ( $activeThreads{$_}->{'thread'}->is_joinable() ) {
			print "Job ".$activeThreads{$_}->{'jobID'}." has finished\n";
			my $exitStatus = $activeThreads{$_}->{'thread'}->join();		  
			&{$activeThreads{$_}->{'job'}->{'onCompletion'}->{'function'}}(@{$activeThreads{$_}->{'job'}->{'onCompletion'}->{'arguments'}},$activeThreads{$_}->{'jobID'},$exitStatus)
			    if ( exists($activeThreads{$_}->{'job'}->{'onCompletion'}) );
			delete($activeThreads{$_});
			$workDone = 1;
		    }
		} else {
		    print "Job ".$activeThreads{$_}->{'jobID'}." has finished\n";
		    my $exitStatus = 0;
		    &{$activeThreads{$_}->{'job'}->{'onCompletion'}->{'function'}}(@{$activeThreads{$_}->{'job'}->{'onCompletion'}->{'arguments'}},$activeThreads{$_}->{'jobID'},$exitStatus)
			if ( exists($activeThreads{$_}->{'job'}->{'onCompletion'}) );
			delete($activeThreads{$_});
		    $workDone = 1;
		}
	   }
	    # Pause if no work was done.
	    sleep(2)
		if ( $workDone == 0);
	}
    }
}

1;
