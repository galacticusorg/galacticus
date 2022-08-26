#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use JSON::PP;

# Run models to validate a Milky Way model.
# Andrew Benson (10-August-2022)

# Make output directory.
system("mkdir -p outputs/");

# Run the validate model.
system("cd ..; ./Galacticus.exe testSuite/parameters/validate_milkyWay.xml");
unless ( $? == 0 ) {
    print "FAIL: Milky Way validation model failed to run\n";
    exit;
}

# Extract and validate the likelihoods.
my @output;
my $model    = new PDL::IO::HDF5("outputs/validate_milkyWay.hdf5");
my $analyses = $model->group('analyses');
foreach my $analysisName ( $analyses->groups() ) {
    my $analysis = $analyses->group($analysisName);
    (my $logLikelihood) = $analysis->attrGet('logLikelihood');
    print $analysisName."\t".$logLikelihood."\n";
    push(
	@output,
	{
	 name  => "Milky Way model - Likelihood - ".$analysisName,
	 unit  => "|logℒ|"                                      ,
	 value => abs($logLikelihood->sclr())
 	}
	);
}
my $json = JSON::PP->new()->pretty()->encode(\@output);
open(my $reportFile,">","outputs/validate_milkyWay.json");
print $reportFile $json;
close($reportFile);

print "SUCCESS: Milky Way validation model\n";

exit;
