#!/usr/bin/env perl
use strict;
use warnings;
use threads;

# A simple wrapper script which launches multiple threads to process the input list of files through Grasil in parallel.
# Andrew Benson (19-July-2011)

# Get a list of the file root names to run through Grasil.
my @grasilFilesRoots = @ARGV;

# Create an array of threads.
my @grasilThreads;

# Get the current working directory.
my $pwd = `pwd`;
chomp($pwd);

if ( scalar(@grasilFilesRoots) > 1 ) {

    # Launch a thread for each given name.
    foreach my $grasilFilesRoot ( @grasilFilesRoots ) {
	$grasilThreads[++$#grasilThreads] = threads->create({'void' => 1},sub {system("cd `dirname ".$grasilFilesRoot."`; ".$pwd."/aux/Grasil/grasil `basename ".$grasilFilesRoot."` &> `basename ".$grasilFilesRoot.".log`")});
    }
    
    # Wait for threads to finish.
    foreach my $grasilThread ( @grasilThreads ) {
	$grasilThread->join();
    }

} else {

    # Run a single instance without launching a thread.
    system("cd `dirname ".$grasilFilesRoots[0]."`; ".$pwd."/aux/Grasil/grasil `basename ".$grasilFilesRoots[0]."` &> `basename ".$grasilFilesRoots[0].".log`");

}

exit;
