#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;

# Test migration of parameter files.
# Andrew Benson (27-January-2023)

# Migrate the test parameter file.
system("mkdir outputs; cd ..; scripts/aux/parametersMigrate.pl testSuite/parameters/parameterMigration.xml testSuite/outputs/parameterMigrated.xml");
if ( $? == 0 ) {
    print "PASSED: migration of parameter file\n";
} else {
    print "FAILED: migration of parameter file\n";
}

# Parse the migrated parameter file.
my $xml        = new XML::Simple();
my $parameters = $xml->XMLin("outputs/parameterMigrated.xml");

# Check expected state.
print "FAILED: missing parameter 'massDestructionAbsolute'\n"
    unless ( exists($parameters->{'nodeOperator'}->{'nodeOperator'}->{'massDestructionAbsolute'}) );
print "FAILED: unremoved parameter spheroidVerySimpleTrackLuminosities''\n"
    if ( exists($parameters->{'spheroidVerySimpleTrackLuminosities'}) );

exit;
