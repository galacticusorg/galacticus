#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test the `nodePropertyExtractorOutputSelector` class.
# Andrew Benson (21-February-2024)

# Make output directory.
system("mkdir -p outputs/");

# Run the model.
system("cd ..; ./Galacticus.exe testSuite/parameters/outputSelector.xml");
unless ( $? == 0 ) {
    print "FAIL: output selector model failed to run\n";
    exit;
}

# Check for correctly selected outputs.
my $model   = new PDL::IO::HDF5("outputs/outputSelector.hdf5");
my $outputs = $model->group('Outputs');
my $success = 1;
foreach my $outputName ( $outputs->groups() ) {
    my $output   = $outputs->group  ($outputName);
    my $nodeData = $output ->group  ('nodeData' );
    (my $expansionFactor) = $output->attrGet('outputExpansionFactor');
    my $redshift = 1.0/$expansionFactor-1.0;
    my $found = grep {$_ eq "mergerTreeIndex"} $nodeData->datasets();
    if ( abs($redshift-0.04) < 1.0e-3 || abs($redshift-2.34) < 1.0e-3 ) {
	# Output property should be present at this redshift.
	$success = 0
	    unless ( $found );
	print $outputName."\tz=".$redshift."\t".($found ? "found (expected)" : "not found (unexpected)")."\n";
    } else {
	# Output property should not be present at this redshift.
	$success = 0
	    if     ( $found );
	print $outputName."\tz=".$redshift."\t".($found ? "found (unexpected)" : "not found (expected)")."\n";
    }
}
my $status = $success ? "SUCCESS" : "FAILED";
print $status.": nodePropertyExtractorOutputSelector\n";
exit;
