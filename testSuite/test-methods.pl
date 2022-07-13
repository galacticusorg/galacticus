#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Options;

# Run a set of short Galacticus models spanning a full range of method options to ensure
# that they at least run to completion.
# Andrew Benson (04-Sep-2010)

# Read in any configuration options.
my $config = &Galacticus::Options::LoadConfig();

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     })
    if ( defined($queueManager) );

# Set default options.
my %options =
    (
     'pbsJobMaximum' => (defined($queueConfig) && exists($queueConfig->{'jobMaximum'})) ? $queueConfig->{'jobMaximum'} : 100,
    );

# Get any command line options.
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Determine dynamic data path.
my $dataDynamicPath = exists($ENV{'GALACTICUS_DYNAMIC_DATA_PATH'}) ? $ENV{'GALACTICUS_DYNAMIC_DATA_PATH'} : $ENV{'GALACTICUS_DATA_PATH'};
# Remove automatically generated files to force them to be regenerated (and, therefore, to test that the generating code works
# correctly).
#
# FSPS stellar population synthesis code and associated file.
system("rm -f  ".$dataDynamicPath."/dynamic/stellarPopulations/SSP_Spectra_Conroy-et-al_v2.5_imfSalpeter.hdf5");
system("rm -rf ".$dataDynamicPath."/dynamic/FSPS_v2.5");
# Core files (older than 7 days).
system("find ../ -name '".$_."' -ctime +7 -exec rm {} \\;")
    foreach ( "core.*", "vgcore.*" );
# Noninstantaneous recycling files (older than 14 days).
system("find ".$dataDynamicPath."/dynamic/stellarPopulations -name '".$_."' -ctime +14 -exec rm {} \\;")
    foreach ( "yield*.hdf5", "recycledFraction*.hdf5", "energyOutput*.hdf5" );
# CAMB transfer function files (older than 14 days).
system("find ".$dataDynamicPath."/dynamic/largeScaleStructure -name 'transfer_function_CAMB_*.xml' -ctime +14 -exec rm {} \\;");

# Simply run the models.
system("cd ..; scripts/aux/launch.pl testSuite/test-methods.xml ".join(" ",map {"--".$_." ".$options{$_}} keys(%options)));

# Check for failed models.
system("grep -q -i -e fatal -e aborted outputs/test-methods/galacticus_*/galacticus.log");
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i -e fatal -e aborted outputs/test-methods/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS!\n";
}

exit;
