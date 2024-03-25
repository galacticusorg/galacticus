#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use System::Redirect;
use List::Uniq ':all';

# Run a Galacticus model to test allowed parameter functionality.
# Andrew Benson (16-June-2017)

# Run the test model.
&System::Redirect::tofile("mkdir -p outputs; cd ..; ./Galacticus.exe testSuite/parameters/test-allowed-parameters.xml","outputs/test-allowed-parameters.log");

# Parse the log file looking for disallowed parameters.
my @disallowed;
open(my $log,"outputs/test-allowed-parameters.log");
while ( my $line = <$log> ) {
    if ( $line =~ m/unrecognized parameter \[([a-zA-Z0-9\-]+) in [a-zA-Z\/]+\]/ ) {
	push(@disallowed,$1);
    }
}
close($log);

# Extract a sorted list of unique disallowed parameters.
@disallowed = uniq(sort(@disallowed));

# Check against expectations.
my $status = join(":",@disallowed) eq "scaleCutOff" ? "success" : "FAILURE";
print "  -> unrecognized parameters: ".join(", ",@disallowed)."\n"
    if ( $status eq "FAILURE" );
print "Test allowed parameters functionality: ".$status."\n";

exit 0;
