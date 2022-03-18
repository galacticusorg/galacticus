#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use File::Find;
use Galacticus::Launch::PBS;
use Galacticus::Options;
use Data::Dumper;

# Run a set of short Galacticus models which test cases that have failed before,
# in order to catch regressions.
# Andrew Benson (01-Feb-2011)

# Indicate that this test can manage its own jobs.
# selfManage: true

# Get PBS configuration and determine number of threads.
my $pbsConfig = &Galacticus::Options::Config("pbs");
my $ppn       = exists($pbsConfig->{'ppn'}) ? $pbsConfig->{'ppn'} : 1;

# Find all regression parameter files and run them.
my $outputDirectory = "outputs/regressions";
system("mkdir -p ".$outputDirectory);
my @regressionDirs = ( "regressions" );
my @pbsJobs;
find(\&runRegressions,@regressionDirs);
my %options =
    (
     pbsJobMaximum       => 64,
     submitSleepDuration =>  1,
     waitSleepDuration   => 10
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

print "\n\n:--> Running regression test cases...\n";
&Galacticus::Launch::PBS::SubmitJobs(\%options,@pbsJobs)
    if ( scalar(@pbsJobs) > 0 );

exit;

sub runRegressions {
    # Run a regression case.
    my $fileName = $_;
    chomp($fileName);

    # Test if this is a parameter file to run.
    if ( $fileName =~ m/^(.*)\.xml$/ ) {
	my $name = $1;
	my %job =
	    (
	     launchFile   => "outputs/regressions/".$name.".pbs",
	     label        => $name,
	     logFile      => "outputs/regressions/".$name.".log",
	     command      => "cd ..; pwd; ./Galacticus.exe testSuite/".$File::Find::dir."/".$fileName,
	     ppn          => $ppn,
	     onCompletion =>
	     {
		 function  => \&testRegressions,
		 arguments => [ $name, "outputs/regressions/".$name.".log" ]
	     }
	    );
	push(@pbsJobs,\%job);
    }

    # Test if this is a script to run.
    if ( $fileName =~ m/^(.*)\.pl$/ ) {
	my $name = $1;
	my %job =
	    (
	     launchFile   => "outputs/regressions/".$name.".pbs",
	     label        => $name,
	     logFile      => "outputs/regressions/".$name.".log",
	     command      => "cd ..; pwd; testSuite/".$File::Find::dir."/".$fileName,
	     ppn          => $ppn,
	     onCompletion =>
	     {
		 function  => \&testRegressions,
		 arguments => [ $name, "outputs/regressions/".$name.".log" ]
	     }
	    );
	push(@pbsJobs,\%job);
    }
}

sub testRegressions {
    # Test a finished regressions model for success.
    my $name    = shift();
    my $logFile = shift();
    my $jobID   = shift();
    my $status  = shift();
    my $failed  = 0;
    if ( $status == 0 ) {
	# Job completed, check the output log.
	system("grep -q -i -e fatal -e aborted ".$logFile);
	$failed = 1
	    if ( $? == 0 );
	system("grep -q FAIL ".$logFile);
	$failed = 1
	    if ( $? == 0 );
    } else {
	# Job failed to complete.
	$failed = 1;
    }
    if ( $failed ) {
	# Job failed .
	print "FAILED: regression test case: ".$name."\n";
	system("cat ".$logFile);
    } else {
	print "SUCCESS: regression test case: ".$name."\n";
    }
}
