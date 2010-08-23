#!/usr/bin/env perl
use lib "./perl";
use System::Redirect;

# Run a suite of tests on the Galacticus code.
# Andrew Benson (19-Aug-2010).

# Open a log file.
$logFile = "allTests.log";
open(lHndl,">".$logFile);

# Define a list of executables to run.
@executablesToRun = (
    "tests.IO.HDF5.exe" # Tests of HDF5 IO routines.
    );

# Run all executables.
foreach $executable ( @executablesToRun ) {
    print lHndl "Running test: ".$executable."\n";
    system("make ".$executable);
    if ( $? == 0 ) {
	# Run the test and copy any output to our log file.
	&SystemRedirect::tofile($executable,"allTests.tmp");
	open(iHndl,"allTests.tmp");
	while ( $line = <iHndl> ) {
	    print lHndl $line;
	}
	close(iHndl);
	unlink("allTests.tmp");
    } else {
	# Build failed, report an error in the log file.
	print lHndl "FAILED: building ".$executable." failed\n";
    }
}

# Close the log file.
close(lHndl);

# Scan the log file for FAILED.
$hasFailures = 0;
open(lHndl,$logFile);
while ( $line = <lHndl> ) {
    if ( $line =~ m/FAILED/ ) {++$hasFailures};
}
close(lHndl);
print $hasFailures." failure(s) found - see ".$logFile." for details.\n" unless ( $hasFailures == 0 );

exit;
