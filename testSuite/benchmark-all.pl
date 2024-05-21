#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Date::Format;
use XML::Simple;
use MIME::Lite;
use Net::SMTP::SSL;
use Data::Dumper;
use File::Slurp qw( slurp );
use File::Find;
use Term::ReadKey;
use System::Redirect;
use Galacticus::Launch::PBS;
use List::ExtraUtils;
use List::Uniq ':all';
use Scalar::Util 'reftype';

# Run a suite of benchmarks on the Galacticus code.
# Andrew Benson (21-May-2019).

# Get command line arguments.
die("benchmark-all.pl <statusPath>")
    unless ( scalar(@ARGV) == 1 );
my $statusPath = $ARGV[0];

# Read in any configuration options.
my $config = &Galacticus::Options::LoadConfig();

# Identify e-mail options for this host.
my $emailConfig;
my $smtpPassword;
if ( exists($config->{'email'}->{'host'}->{$ENV{'HOSTNAME'}}) ) {
    $emailConfig = $config->{'email'}->{'host'}->{$ENV{'HOSTNAME'}};
} elsif ( exists($config->{'email'}->{'host'}->{'default'}) ) {
    $emailConfig = $config->{'email'}->{'host'}->{'default'};
} else {
    $emailConfig->{'method'} = "sendmail";
}
if ( $emailConfig->{'method'} eq "smtp" && exists($emailConfig->{'passwordFrom'}) ) {
    # Get any password now.
    if ( $emailConfig->{'passwordFrom'} eq "input" ) {
	print "Please enter your e-mail SMTP password:\n";
	$smtpPassword = &getPassword;
    }
}

my $repoRevisionHash;
open(my $gitHndl,"git rev-parse HEAD|");
while ( my $line = <$gitHndl> ) {
    $repoRevisionHash = $line;
    chomp($repoRevisionHash);
}
close($gitHndl);

# Determine the current repo branch.
open(my $gitBranch,"git branch |");
my $repoBranchName = <$gitBranch>;
chomp($repoBranchName);
$repoBranchName =~ s/^\*\s*//;
close($gitBranch);

# Array of benchmark results.
my @benchmarks;

# Open a log file.
my $logFile = "testSuite/allBenchmarks.log";
open(lHndl,">".$logFile);
 
# Write header to log file.
print lHndl ":-> Running benchmarks:\n";
print lHndl "    -> Host:\t".$ENV{'HOSTNAME'}."\n";
print lHndl "    -> Time:\t".time2str("%a %b %e %T (%Z) %Y", time)."\n";

# Create a directory for test suite outputs.
system("rm -rf testSuite/benchmark-outputs");
system("mkdir -p testSuite/benchmark-outputs");

# Stack to be used for PBS jobs.
my @jobStack;

# Set options for PBS launch.
my %pbsOptions =
    (
     pbsJobMaximum       => 100,
     submitSleepDuration =>   1,
     waitSleepDuration   =>  10
    );

# Define a list of benchmarks to run. Each entry must give the name of the executable or script to run.
my @benchmarksToRun = (
    {
	name    => "stellarLuminosities"                ,
	compile => "benchmarks.stellar_luminosities.exe",
	run     => "benchmarks.stellar_luminosities.exe"
    },
    {
	name    => "quickTest"                        ,
	compile => "Galacticus.exe"                   ,
	ppn     => 16                                 ,
	run     => "testSuite/benchmarks.quickTest.pl"
    }
    );

# Build all benchmarks.
my @compileList = map {exists($_->{'compile'}) ? $_->{'compile'} : ()} @benchmarksToRun;
@compileList = uniq(sort(@compileList));
my $compileCommand = "rm -rf work/build; rm *.exe;make -j16 ".join(" ",@compileList)."\n";
my %benchmarkBuildJob   =
    (
     launchFile   => "testSuite/benchmark-outputs/compileBenchmarks.pbs",
     label        => "testSuite-compileBenchmarks"                      ,
     logFile      => "testSuite/benchmark-outputs/compileBenchmarks.log",
     command      =>  $compileCommand                                   ,
     ppn          => 16                                                 ,
     tracejob     => "yes"                                              ,
     onCompletion => 
     {
	 function  => \&benchmarkCompileFailure,
	 arguments => [ "testSuite/benchmark-outputs/compileBenchmarks.log", "Benchmark code compilation" ]
     }
    );
push(@jobStack,\%benchmarkBuildJob);
&Galacticus::Launch::PBS::SubmitJobs(\%pbsOptions,@jobStack);
unlink("testSuite/benchmark-outputs/compileBenchmarks.pbs");

# Launch all benchmarks.
my @launchFiles;
@jobStack = ();
foreach my $benchmark ( @benchmarksToRun ) {
    # Generate the job.
    if ( exists($benchmark->{'name'}) ) {
	my $label = $benchmark->{'name'};
	my $ppn = exists($benchmark->{'ppn'}) ? $benchmark->{'ppn'} : 1;
	my $launchFile = "testSuite/benchmark-outputs/".$label.".pbs";
	push(@launchFiles,$launchFile);
	$benchmark->{'expect'} = "success";
	my %job =
	    (
	     launchFile   => $launchFile               ,
	     label        => "testSuite-".$label       ,
	     logFile      => "testSuite/benchmark-outputs/".$label.".log",
	     ppn          => $ppn                      ,
	     tracejob     => "yes"                     ,
	     onCompletion => 
	     {
		 function  => \&benchmarkRecord,
		 arguments => [ "testSuite/benchmark-outputs/".$label.".log", "Benchmark code: ".$benchmark->{'name'}, $benchmark->{'expect'} ]
	     }
	    );
	$job{'command'} = $benchmark->{'run'};
	push(@jobStack,\%job);
    }
}
&Galacticus::Launch::PBS::SubmitJobs(\%pbsOptions,@jobStack);
unlink(@launchFiles);

# Close the log file.
close(lHndl);

# Scan the log file for FAILED.
my $lineNumber = 0;
my @failLines;
open(lHndl,$logFile);
while ( my $line = <lHndl> ) {
    ++$lineNumber;
    if ( $line =~ m/FAILED/ ) {
	push(@failLines,$lineNumber);
    }
    if ( $line =~ m/SKIPPED/ ) {
	push(@failLines,$lineNumber);
    }
}
close(lHndl);
open(lHndl,">>".$logFile);
my $emailSubject = "Galacticus benchmarks log";
my $exitStatus;
if ( scalar(@failLines) == 0 ) {
    print lHndl "\n\n:-> All benchmarks were successful.\n";
    print       "All benchmarks were successful.\n";
    $emailSubject .= " [success]";
    $exitStatus = 0;
} else {
    print lHndl "\n\n:-> Failures found. See following lines in log file:\n\t".join("\n\t",@failLines)."\n";
    print "Failure(s) found - see ".$logFile." for details.\n";
    $emailSubject .= " [FAILURE]";
    $exitStatus = 1;
}
close(lHndl);

# Update benchmarks file.
if ( -e $statusPath."/benchmarks.xml" ) {
    my $xml = new XML::Simple();
    my $status = $xml->XMLin($statusPath."/benchmarks.xml");
    foreach my $branchName ( keys(%{$status->{'galacticus'}}) ) {
	if ( reftype($status->{'galacticus'}->{$branchName}) eq "HASH" ) {
	    my $content = $status->{'galacticus'}->{$branchName};
	    delete($status->{'galacticus'}->{$branchName});
	    $status->{'galacticus'}->{$branchName}->[0] = $content;
	}
    }
    foreach my $benchmark ( @benchmarks ) {
	my $label = $benchmark->{'label'};
	delete($benchmark->{'label'});	
	my @benchmarkList;
	if ( exists($status->{'galacticus'}->{$repoBranchName}->[0]->{$label}) ) {
	    @benchmarkList = &List::ExtraUtils::as_array($status->{'galacticus'}->{$repoBranchName}->[0]->{$label});
	    delete($status->{'galacticus'}->{$repoBranchName}->[0]->{$label});
	}
	push(@benchmarkList,$benchmark);
	@{$status->{'galacticus'}->{$repoBranchName}->[0]->{$label}} = @benchmarkList;
    }
    open(my $output,">",$statusPath."/benchmarks.xml");
    print $output $xml->XMLout($status,rootName => "benchmarks");
    close($output);
}

# If we have an e-mail address to send the log to, then do so.
if ( defined($config->{'contact'}->{'email'}) ) {
    if ( $config->{'contact'}->{'email'} =~ m/\@/ ) {
	# Get e-mail configuration.
	my $sendMethod = $emailConfig->{'method'};
	# Construct the message.
	my $message  = "Galacticus benchmark log is attached.\n";
	my $msg = MIME::Lite->new(
	    From    => '',
	    To      => $config->{'contact'}->{'email'},
	    Subject => $emailSubject,
	    Type    => 'TEXT',
	    Data    => $message
	    );
	system("bzip2 -f ".$logFile);
	$msg->attach(
	    Type     => "application/x-bzip",
	    Path     => $logFile.".bz2",
	    Filename => "allBenchmarks.log.bz2"
	    );
	if ( $sendMethod eq "sendmail" ) {
	    $msg->send;
	}
	elsif ( $sendMethod eq "smtp" ) {
	    my $smtp; 
	    $smtp = Net::SMTP::SSL->new($config->{'email'}->{'host'}, Port=>465) or die "Can't connect";
	    $smtp->auth($config->{'email'}->{'user'},$smtpPassword) or die "Can't authenticate:".$smtp->message();
	    $smtp->mail( $config->{'contact'}->{'email'}) or die "Error:".$smtp->message();
	    $smtp->to( $config->{'contact'}->{'email'}) or die "Error:".$smtp->message();
	    $smtp->data() or die "Error:".$smtp->message();
	    $smtp->datasend($msg->as_string) or die "Error:".$smtp->message();
	    $smtp->dataend() or die "Error:".$smtp->message();
	    $smtp->quit() or die "Error:".$smtp->message();
	}
    }
}

exit $exitStatus;

sub getPassword {
    # Read a password from standard input while echoing asterisks to the screen.
    ReadMode('noecho');
    ReadMode('raw');
    my $password = '';
    while (1) {
	my $c;
	1 until defined($c = ReadKey(-1));
	last if $c eq "\n";
	print "*";
	$password .= $c;
    }
    ReadMode('restore');
    print "\n";
    return $password;
}

sub benchmarkRecord {
    # Callback function which checks for failure of jobs run in PBS and records benchmarking results.
    my $logFile    = shift();
    my $jobMessage = shift();
    my $expect     = shift();
    my $jobID      = shift();
    my $exitStatus = shift();
    my $result;
    # Check for failure message in log file.
    system("grep -q FAIL ".$logFile);
    if ( $? == 0 ) {
	$result = "FAILED: ".$jobMessage."\n";
    } else {
	$result = "SUCCESS: ".$jobMessage."\n";
    }
    # Report success or failure.
    print lHndl $result;
    if ( $result =~ m/^FAILED:/ ) {
	# Job failed.
	print lHndl "Job output follows:\n";
	print lHndl slurp($logFile);
    } else {
	# Scan for benchmark results.
	open(my $log,$logFile);
	while ( my $line = <$log> ) {
	    if ( $line =~ m/^BENCHMARK\s+([a-zA-Z_]+)\s+"([^"]+)"\s+([\d\.]+)\s+([\d\.]+)\s+"([^"]+)"/ ) {
		my $benchmark =
		{
		    label       => $1                                   ,
		    description => $2                                   ,
		    time        => $3                                   ,
		    uncertainty => $4                                   ,
		    units       => $5                                   ,
		    hash        => $repoRevisionHash                    ,
		    timestamp   => time2str("%a %b %e %T (%Z) %Y", time)
		};
		push(@benchmarks,$benchmark);
	    }
	}
	close($log);
    }
}

sub benchmarkCompileFailure {
    # Callback function which checks for failure of compile jobs run in PBS.
    my $logFile     = shift();
    my $jobMessage  = shift();
    my $jobID       = shift();
    my $errorStatus = shift();
    # Check for failure message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q FAIL ".$logFile);
	$errorStatus = 1
	    if ( $? == 0 );	
    }
    # Check for compiler error message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q Error: ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = "Compiler errors issued\n".$jobMessage;
	}
    }
    # Check for compiler warning message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q Warning: ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = "Compiler warnings issued\n".$jobMessage;
	}
    }
    # Check for make error message in log file.
    if ( $errorStatus == 0 ) {
	system("grep -q make: ".$logFile);
	if ( $? == 0 ) {
	    $errorStatus = 1;
	    $jobMessage = "Make errors issued\n".$jobMessage;
	}
    }
    # Report success or failure.
    if ( $errorStatus == 0 ) {
	# Job succeeded.
	print lHndl "SUCCESS: ".$jobMessage."\n";
	unlink($logFile);
    } else {
	# Job failed.
	print lHndl "FAILED: ".$jobMessage."\n";
	print lHndl "Job output follows:\n";
	print lHndl slurp($logFile);
    }
}
