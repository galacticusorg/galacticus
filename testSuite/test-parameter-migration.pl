#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Data::Dumper;
use List::ExtraUtils;

# Test migration of parameter files.
# Andrew Benson (27-January-2023)

# Migrate the test parameter file.
system("mkdir -p outputs; cd ..; scripts/aux/parametersMigrate.pl testSuite/parameters/parameterMigration.xml testSuite/outputs/parameterMigrated.xml --lastModifiedRevision 6eab8997cd73cb0a474228ade542d133890ad138^");
if ( $? == 0 ) {
    print "PASSED: migration of parameter file\n";
} else {
    print "FAILED: migration of parameter file\n";
}

# Parse the migrated parameter file.
my $xml        = new XML::Simple();
my $parameters = $xml->XMLin("outputs/parameterMigrated.xml");

# Check expected state.
my @nodeOperators = &List::ExtraUtils::as_array($parameters->{'nodeOperator'}->{'nodeOperator'});
print "FAILED: missing parameter 'massDestructionAbsolute'\n"
    unless ( exists($nodeOperators[0]->{'massDestructionAbsolute'}) );
print "FAILED: unremoved parameter 'spheroidVerySimpleTrackLuminosities'\n"
    if ( exists($parameters->{'spheroidVerySimpleTrackLuminosities'}) );

exit;
