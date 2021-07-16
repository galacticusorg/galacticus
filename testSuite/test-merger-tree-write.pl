#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use File::Which;
use Galacticus::Options;

# Export trees from Galacticus and check that they are written correctly.
# Andrew Benson (12-October-2012)

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

# Run models.
system("cd ..; mkdir -p testSuite/outputs/test-merger-tree-write; scripts/aux/launch.pl testSuite/parameters/test-merger-tree-write.xml ".join(" ",map {"--".$_." ".$options{$_}} keys(%options))."; scripts/aux/launch.pl testSuite/parameters/test-merger-tree-write-secondary.xml ".join(" ",map {"--".$_." ".$options{$_}} keys(%options)));

# Check for failed models.
system("grep -q -i fatal outputs/test-merger-tree-write/galacticus_*/galacticus.log");
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i fatal outputs/test-merger-tree-write/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS: model run\n";
}

# Validate the IRATE-format output.
my $validator = &File::Which::which('iratevalidate');
if ( $validator ) {
    system("iratevalidate outputs/test-merger-tree-write/exportedTreesIRATE.hdf5");
    die("FAILED: IRATE-format file ouput by Galacticus did not validate")
	unless ( $? == 0 );
} else {
    print "SKIP: iratevalidate is not installed - validation of IRATE-format file will be skipped\n";
}

exit;
