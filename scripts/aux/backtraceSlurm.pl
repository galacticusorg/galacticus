#!/usr/bin/env perl
use strict;
use warnings;

# Obtain stack traces from all threads of all running Galacticus processes that are part of a given Slurm job and output to file.
# Andrew Benson (14-May-2025)

# Get job ID.
die("Usage: backtraceSlurm.pl <slurmJobID> <outputFile>")
    unless ( scalar(@ARGV) == 2 );
my $jobID      = $ARGV[0];
my $outputFile = $ARGV[1];

# Parse the node list.
my $nodeList;
open(my $job,"scontrol show job ".$jobID." |");
while ( my $line = <$job> ) {
    if ( $line =~ m/^\s*NodeList=(.+)/ ) {
	$nodeList = $1;
    }
}
close($job);
die("unable to find NodeList")
    unless ( $nodeList );

# Extract cannonical node names.
my @nodeNames;
open(my $nodeListPipe,"scontrol show hostnames ".$nodeList." |");
while ( my $line = <$nodeListPipe> ) {
    chomp($line);
    push(@nodeNames,$line);
}
close($nodeListPipe);

# Open output file.
open(my $output,">",$outputFile);

# Iterate over nodes.
foreach my $nodeName ( @nodeNames ) {
    print $output "Node: ".$nodeName."\n";
    open(my $srun,"srun --jobid=".$jobID." -w ".$nodeName." --pty bash -c 'pgrep -u ".$ENV{'USER'}." \"Galacticus.exe\" | xargs -i{} gdb --pid {} -batch -ex \"info threads\" -ex \"thread apply all where\"' 2>&1 |");
    while ( my $line = <$srun> ) {
	print $output $line;
    }    
    close($srun);
    # Reset the terminal as, for some reason, the above `srun` command seems to leave it in a weird state.
    system("reset");
    print $output "\n\n";
}
close($output);
