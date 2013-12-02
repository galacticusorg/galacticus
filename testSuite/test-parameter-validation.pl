#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;

# Validate a set of test parameter files.
# Andrew Benson (01-December-2013)

# Validate the default paramter file.
system("cd ..; scripts/aux/validateParameters.pl parameters.xml");
if ( $? == 0 ) {
    print "PASSED: validation of default parameter file\n";
} else {
    print "FAILED: validation of default parameter file\n";
}

# Find all validation parameter files and run them.
my @validationDirs = ( "parameters/validation" );
find(\&runValidations,@validationDirs);

exit;

sub runValidations {
    # Run a validation case.
    my $fileName = $_;
    chomp($fileName);
    # Test if this is a parameter fil to run.
    if ( $fileName =~ m/\-(valid|invalid)\.xml$/ ) {
	my $validity = $1;
	system("cd ../../..; scripts/aux/validateParameters.pl testSuite/parameters/validation/".$fileName);
	if (
	    ( $? == 0 && $validity eq   "valid" )
	    ||
	    ( $? != 0 && $validity eq "invalid" )
	    ) {
	    print "PASSED: validation of '".$fileName."'\n";
	} else {	
	    print "FAILED: validation of '".$fileName."'\n";
	}
    }
}
