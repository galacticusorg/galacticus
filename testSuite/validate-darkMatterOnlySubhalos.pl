#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use JSON::PP;

# Run models to validate a dark matter only subhalo evolution model.
# Andrew Benson (05-August-2022)

# Make output directory.
system("mkdir -p outputs/");

# Run the validate model.
system("cd ..; ./Galacticus.exe testSuite/parameters/validate_darkMatterOnlySubHalos.xml");
unless ( $? == 0 ) {
    print "FAIL: dark matter-only subhalos validation model failed to run\n";
    exit;
}

# Extract and validate the likelihoods.
my @output;
my $model    = new PDL::IO::HDF5("outputs/validate_darkMatterOnlySubHalos.hdf5");
my $analyses = $model->group('analyses');
foreach my $analysisName ( $analyses->groups() ) {
    my $analysis = $analyses->group($analysisName);
    (my $logLikelihood) = $analysis->attrGet('logLikelihood');
    print $analysisName."\t".$logLikelihood."\n";
    push(
	@output,
	{
	 name  => "Dark Matter Only Subhalos - Likelihood - ".$analysisName,
	 unit  => "|logℒ|"                                                ,
	 value => abs($logLikelihood->sclr())
 	}
	);
}
my $json = JSON::PP->new()->pretty()->encode(\@output);
open(my $reportFile,">","outputs/validate_darkMatterOnlySubhalos.json");
print $reportFile $json;
close($reportFile);

print "SUCCESS: dark matter-only subhalos validation model\n";

exit;
