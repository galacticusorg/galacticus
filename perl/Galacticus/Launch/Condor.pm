# Launch models on Condor system.

package Condor;
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

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     condor => {
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
	 galacticusDirectory     => "/home/condor/Galacticus/v0.9.3"                        ,
	 universe                => "vanilla"                                               ,
	 environment             => ""                                                      ,
	 requirement             => []                                                      ,
	 transferFiles           => [ $galacticusPath."/Galacticus.exe", "parameters.xml"  ],
	 postSubmitSleepDuration =>  5                                                      ,
	 jobWaitSleepDuration    => 10
	);
    foreach ( keys(%defaults) ) {
	if ( UNIVERSAL::isa($defaults{$_},"ARRAY") ) {
	    push(
		@{$launchScript->{'condor'}->{$_}},
		@{$defaults                  {$_}}
		);
	} else {
	    $launchScript->{'condor'}->{$_} = $defaults{$_}
	        unless ( exists($launchScript->{'condor'}->{$_}) );
	}
    }
}

sub Output_File_Name {
    # Return the output file name for Galacticus models.
    return "galacticus.hdf5";
}

sub Launch {
    # Launch models on Condor cluster.
    my @jobs         = @{shift()};
    my $launchScript =   shift() ;
    # Set pwd for files to transfer.
    my $pwd = `pwd`;
    chomp($pwd);
    $_ =~ s/{PWD}/$pwd/
	foreach ( @{$launchScript->{'condor'}->{'transferFiles'}} );
    # Iterate over jobs.
    my %condorJobs;
    foreach my $job ( @jobs ) {
	# Create the script that Condor will execute.
	my $condorScript = "condor_run_".$launchScript->{'modelCounter'}."_".$$.".csh";
	open(my $oHndl,">".$condorScript);
	print $oHndl "#!/bin/csh\n";
	print $oHndl "ln -sf ".$launchScript->{'condor'}->{'galacticusDirectory'}."/aux\n";
	print $oHndl "ln -sf ".$launchScript->{'condor'}->{'galacticusDirectory'}."/data\n";
	print $oHndl "ln -sf ".$launchScript->{'condor'}->{'galacticusDirectory'}."/perl\n";
	print $oHndl "ln -sf ".$launchScript->{'condor'}->{'galacticusDirectory'}."/scripts\n";
	print $oHndl "ln -sf ".$launchScript->{'condor'}->{'galacticusDirectory'}."/work\n";
	print $oHndl "exec ./Galacticus.exe parameters.xml\n";
	print $oHndl "rm -f aux data perl scripts work\n";
	print $oHndl "exit\n";
	close($oHndl);
	# Creat a submit file.
	my $condorSubmit = "condor_submit_".$launchScript->{'modelCounter'}."_".$$.".txt";
	open($oHndl,">".$condorSubmit);
	print $oHndl "Executable              = ".$condorScript."\n";
	print $oHndl "Error                   = condor.err\n";
	print $oHndl "Output                  = condor.out\n";
	print $oHndl "Log                     = condor.log\n";
	print $oHndl "InitialDir              = ".$job->{'directory'}."\n";
	print $oHndl "Should_Transfer_Files   = YES\n";
	print $oHndl "When_To_Transfer_Output = ON_EXIT\n";
	print $oHndl "Transfer_Input_Files    = ".join(",",@{$launchScript->{'condor'}->{'transferFiles'}})."\n";
	print $oHndl "Universe                = ".$launchScript->{'condor'}->{'universe'}."\n";
	print $oHndl "Allow_Startup_Script    = True\n" 
	    if ( $launchScript->{'condor'}->{'universe'} eq "standard" );
	print $oHndl "Environment             = ".$launchScript->{'condor'}->{'environment'}."\n"
	    unless ( $launchScript->{'condor'}->{'environment'} eq "" );
	print $oHndl "Requirements            = (".join(") && (",@{$launchScript->{'condor'}->{'requirement'}}).")\n"
	    if ( scalar(@{$launchScript->{'condor'}->{'requirement'}}) > 0 );
	print $oHndl "+RequiresWholeMachine   = True\n"
	    if ( $launchScript->{'condor'}->{'wholeMachine'} eq "true" );
	print $oHndl "Queue\n";
	close($oHndl);
	# Submit the job - capture the job number.
	open(my $pHndl,"condor_submit -verbose ".$condorSubmit."|");
	my $jobID = "";
	while ( my $line = <$pHndl> ) {
	    if ( $line =~ m/^\*\*\s+Proc\s+([\d\.]+):/ ) {$jobID = $1};
	}
	close($pHndl);
	# Add job number to active job hash
	unless ( $jobID eq "" ) {
	    $condorJobs{$jobID}->{'job'} = $job;
	    @{$condorJobs{$jobID}->{'temporaryFiles'}} = [$condorScript,$condorSubmit];
	}
	sleep 10;	
    }
    # Wait for Condor models to finish.
    print " -> waiting for Condor jobs to finish...\n";
    while ( scalar(keys %condorJobs) > 0 ) {
	# Find all Condor jobs that are running.
	my %runningCondorJobs;
	undef(%runningCondorJobs);
	my $condorQuerySuccess = 0;
	open(my $pHndl,"condor_q|");
	while ( my $line = <$pHndl> ) {
	    if ( $line =~ m/^\s*(\d+\.\d+)\s/ ) {$runningCondorJobs{$1} = 1};
	    if ( $line =~ m/^\d+\s+jobs;/ ) {$condorQuerySuccess = 1};
	}
	close($pHndl);
	if ( $condorQuerySuccess == 1 ) {
	    foreach my $jobID ( keys(%condorJobs) ) {
		unless ( exists($runningCondorJobs{$jobID}) ) {
		    print " -> Condor job ".$jobID." has finished; post-processing....\n";
		    &PostProcess::Analyze($condorJobs{$jobID}->{'job'},$launchScript);
		    &PostProcess::CleanUp($condorJobs{$jobID}->{'job'},$launchScript);
		    # Remove any temporary files associated with this job.
		    unlink(@{$condorJobs{$jobID}->{'temporaryFiles'}});
		    # Remove the job ID from the list of active Condor jobs.
		    delete($condorJobs{$jobID});
		}
	    }
	    sleep $launchScript->{'condor'}->{'postSubmitSleepDuration'};
	} else {
	    sleep $launchScript->{'condor'}->{'jobWaitSleepDuration'   };
	}
    }
}

1;
