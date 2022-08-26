#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Stats::Basic;
use JSON::PP;

# Run models to benchmark performance of a Milky Way model.
# Andrew Benson (10-August-2022)

# Make output directory.
system("mkdir -p outputs/");

# Run the benchmark model multiple times.
my $runTimes = pdl [];
for(my $i=0;$i<11;++$i) {
    # Run the model.
    system("cd ..; /usr/bin/time --format=\"\%e\" --output=testSuite/outputs/benchmark_milkyWay.log ./Galacticus.exe testSuite/parameters/benchmark_milkyWay.xml");
    unless ( $? == 0 ) {
	print "FAIL: Milky Way benchmark model failed to run\n";
	exit;
    }
    # Extract timing data.
    open(my $logFile,"outputs/benchmark_milkyWay.log");
    my $runTime = <$logFile>;
    close($logFile);
    chomp($runTime);
    $runTimes = $runTimes->append($runTime)
	unless ( $i == 0 );
}

# Find average and standard deviation of run times.
my $runTimeAverage           = $runTimes->average();
my $runTimeStandardDeviation = $runTimes->stdv   ();
print "Benchmark results: ".$runTimeAverage." ± ".$runTimeStandardDeviation." s";

# Generate JSON report.
my @output =
    (
     {
	 name  => "Milky Way model - Wall Time",
	 unit  => "seconds"                              ,
	 value => $runTimeAverage          ->sclr()      ,
	 range => $runTimeStandardDeviation->sclr()
     }
    );

my $json = JSON::PP->new()->pretty()->encode(\@output);
open(my $reportFile,">","outputs/benchmark_milkyWay.json");
print $reportFile $json;
close($reportFile);

print "SUCCESS: Milky Way benchmark model\n";

exit;
