#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use threads;

# A simple wrapper script which launches multiple threads to process the input list of scripts.
# Andrew Benson (31-May-2016)

# Get the list of scripts (and log files).
my @tasks = @ARGV;

# Create an array of threads.
my @threads;

# Launch a thread for each given script.
while ( scalar(@tasks) > 0 ) {
    my $scriptFileName = shift(@tasks);
    my $logFileName    = shift(@tasks);
    my $command        = "chmod u+x ".$scriptFileName."; ".$scriptFileName." >& ".$logFileName;
    if ( $^V lt v5.14.0 ) { 
	push(@threads,threads->create( sub {system($command);}));
    } else {
	push(@threads,threads->create({'context' => 'void'}, sub {system($command);}));
    }
}
    
# Wait for threads to finish.
foreach my $thread ( @threads ) {
    $thread->join();
}

exit;
