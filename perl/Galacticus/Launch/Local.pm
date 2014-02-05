# Launch models on local system.

package Local;
use strict;
use warnings;
use Data::Dumper;
use Sys::CPU;
use threads;
require Galacticus::Launch::Hooks;
require Galacticus::Launch::PostProcess;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     local => {
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
	 threadCount => 1        ,
	 ompThreads  => "maximum"
	);
    foreach ( keys(%defaults) ) {
	$launchScript->{'local'}->{$_} = $defaults{$_}
	unless ( exists($launchScript->{'local'}->{$_}) );
    }
    # Determine how many threads to launch.
    $launchScript->{'local'}->{'threadCount'} = Sys::CPU::cpu_count() 
	if ( $launchScript->{'local'}->{'threadCount'} eq "maximum" );
    $launchScript->{'local'}->{'ompThreads'} = Sys::CPU::cpu_count() 
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
    # Launch model threads.
    my @threads;
    for(my $iThread=0;$iThread<$launchScript->{'local'}->{'threadCount'};++$iThread) {
	print " -> launching thread ".$iThread." of ".$launchScript->{'local'}->{'threadCount'}."\n"
		if ( $launchScript->{'verbosity'} > 0 );
 	push(
	    @threads,
	    threads->create(\&Launch_Models, $iThread, \@jobs, $launchScript)
	    );
   }
    # Wait for threads to finish.
    foreach my $thread ( @threads ) {
 	my @returnData = $thread->join();
    }
}

sub Launch_Models {
    my $iThread      =   shift() ;
    my @jobs         = @{shift()};
    my $launchScript =   shift() ;
    for(my $i=0;$i<scalar(@jobs);++$i) {
	if ( ( $i % $launchScript->{'local'}->{'threadCount'} ) == $iThread ) {
 	    print " -> thread ".$iThread." running job: ".$jobs[$i]->{'label'}."\n"
		if ( $launchScript->{'verbosity'} > 0 );
	    system(
		    "ulimit -t unlimited;"        .
		    "ulimit -c unlimited;"        .
		    "export GFORTRAN_ERROR_DUMPCORE=YES;".
		    "export OMP_NUM_THREADS=".$launchScript->{'local'}->{'ompThreads'}.";".
		    "Galacticus.exe ".$jobs[$i]->{'directory'}."/parameters.xml &> ".
		    $jobs[$i]->{'directory'}."/galacticus.log"
		);
	    if ( $? == 0 ) {
		&PostProcess::Analyze($jobs[$i],$launchScript);
	    } else {
		&PostProcess::Failed ($jobs[$i],$launchScript);
	    }
	    &PostProcess::CleanUp    ($jobs[$i],$launchScript);
	}
    }
}

1;
