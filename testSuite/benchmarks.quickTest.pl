#!/usr/bin/env perl
use strict;
use warnings;

# Time running our benchmark.
system("/usr/bin/time -o testSuite/benchmark-outputs/quickTest.log --format=\"%S %U\" Galacticus.exe testSuite/parameters/benchmark-quickTest.xml");
if ( $? == 0 ) {
    open(my $logFile,"testSuite/benchmark-outputs/quickTest.log");
    my $benchmarkTimes = <$logFile>;
    close($logFile);
    chomp($benchmarkTimes);
    my @times = split(" ",$benchmarkTimes);
    my $benchmarkTime = $times[0]+$times[1];
    print "BENCHMARK quickTest \"Quick test of Galacticus galaxy evolvution\" ".$benchmarkTime." 0.0 \"s\"\n";
} else {
    print "FAIL: failed to run quickTest benchmark\n";
}

# Run our benchmark in Callgrind and stage the results for archiving.
system("mkdir -p nightlies-archive; valgrind --tool=callgrind --callgrind-out-file=nightlies-archive/benchmarks.quickTest.callgrind Galacticus.exe testSuite/parameters/benchmark-quickTest.xml");

exit 0;
