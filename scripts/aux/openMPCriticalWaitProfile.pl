#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Extract metadata on OpenMP critical section wait times from a Galacticus run.
# Andrew Benson (30-November-2016)

# Get the name of the file to profile.
die "Usage: openMPCriticalWaitProfile.pl <modelFileName>"
    unless ( scalar(@ARGV) == 1 );
my $modelFileName = $ARGV[0];
# Extract meta-data from file.
my $modelFile                = new PDL::IO::HDF5($modelFileName);
my $metaDataGroup            = $modelFile    ->group  ('metaData'                )       ;
my $openMPGroup              = $metaDataGroup->group  ('openMP'                  )       ;
my $criticalSectionNames     = $openMPGroup  ->dataset('criticalSectionNames'    )->get();
my $criticalSectionWaitTimes = $openMPGroup  ->dataset('criticalSectionWaitTimes')->get();
# Find total critical section wait time.
my $waitTimeTotal = $criticalSectionWaitTimes->sum();
# Report.
my $i                  = 0    ;
my $waitTimePercentage = 100.0;
my $waitTimeRanks      = $criticalSectionWaitTimes->qsorti();
print "Total wait time for critical section across all threads = ".$waitTimeTotal." s\n";
print "Significant contributing critical sections are:\n";
while ( $i < 10 || $waitTimePercentage > 10.0 ) {
    ++$i;
    $waitTimePercentage     = 100.0*$criticalSectionWaitTimes->($waitTimeRanks)->((-$i))/$waitTimeTotal;
    (my $criticalSectionName = $criticalSectionNames->(:,($waitTimeRanks->((-$i))))) =~ s/'(\S+)\s*'/$1/;
    print " -> ".sprintf("%8.3f",$waitTimePercentage)."% : ".$criticalSectionName."\n";
}
exit 0;
