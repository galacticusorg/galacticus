#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";

# Postprocess output from the linker to remove irrelevant warnings.
# Andrew Benson (12-April-2024)

# Get the name of the reference file.
die("Usage: postprocessLinker.pl")
    unless ( scalar(@ARGV) == 0 );

my $status = 0;
while ( my $line = <STDIN> ) {
    my $dropBuffer = 0;
    # Drop common warnings emitted by the linker during static lniking that we can safely ignore.
    $dropBuffer = 1
	if (
	    $line =~ m/warning: the use of `mktemp' is dangerous, better use `mkstemp'/
	    ||
	    $line =~ m/warning: Using 'dlopen' in statically linked applications requires at runtime the shared libraries from the glibc version used for linking/
	    ||
	    $line =~ m/warning: ([^:]+): requires executable stack \(because the \.note\.GNU\-stack section is executable\)/
	);
    $line = ""
	if ( $dropBuffer );
    print $line;
    $status = 1
	if ( $line =~ m/^Error:/ );
}

exit $status;
