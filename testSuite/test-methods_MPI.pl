#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Options;

# Run a set of short Galacticus models under MPI spanning a range of method options to ensure that they at least run to
# completion.
# Andrew Benson (14-Jun-2018)

# Read in any configuration options.
my $config;
if ( -e "galacticusConfig.xml" ) {
    my $xml = new XML::Simple;
    $config = $xml->XMLin("galacticusConfig.xml");
}

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Set default options.
my %options =
    (
     'pbsJobMaximum' => exists($queueConfig->{'jobMaximum'}) ? $queueConfig->{'jobMaximum'} : 100,
    );

# Get any command line options.
&Galacticus::Options::Parse_Options(\@ARGV,\%options);


# Simply run the models.
system("cd ..; scripts/aux/launch.pl testSuite/test-methods_MPI.xml ".join(" ",map {"--".$_." ".$options{$_}} keys(%options)));

# Check for failed models.
system("grep -q -i fatal outputs/test-methods_MPI/galacticus_*/galacticus.log");
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i fatal outputs/test-methods_MPI/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS!\n";
}

exit;
