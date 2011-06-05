#!/usr/bin/env perl
use File::Find;

# Run a set of short Galacticus models which test cases that have failed before,
# in order to catch regressions.
# Andrew Benson (01-Feb-2011)

# Find all regression parameter files and run them.
$outputDirectory = "outputs/regressions";
system("mkdir -p ".$outputDirectory);
@regressionDirs = ( "regressions" );
find(\&runRegressions,@regressionDirs);

exit;

sub runRegressions {
    # Run a regression case.
    $fileName = $_;
    chomp($fileName);

    # Test if this is a script to run.
    if ( $fileName =~ m/\.xml$/ ) {
	print "\n\n:--> Running regression test case: ".$fileName."\n";
	system("cd ../..; Galacticus.exe testSuite/".$File::Find::dir."/".$fileName);       
	print "FAILED: regression test case: ".$fileName."\n" unless ( $? == 0 );
    }
}
