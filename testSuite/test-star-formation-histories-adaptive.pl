#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use File::Slurp qw(slurp);
use System::Redirect;

# Check internal self-consistency of adaptive star formation histories.
# Andrew Benson (25-March-2021)

# Run the model.
system("mkdir -p outputs/test-star-formation-histories-adapative");
&System::Redirect::tofile("cd ..; ./Galacticus.exe testSuite/parameters/test-star-formation-histories-adaptive.xml","outputs/test-star-formation-histories-adapative/galacticus.log");

# Check for failed models.
system("grep -q -i -e fatal -e \"Galacticus experienced an error in the GSL library\" outputs/test-star-formation-histories-adapative/galacticus.log");
if ( $? == 0 ) {
    print "FAILED: model run:\n";
    system("cat outputs/test-star-formation-histories-adapative/galacticus.log");
} else {
    print "SUCCESS: model run\n";
}

# Read the model data and check for consistency.
my $model                  = new PDL::IO::HDF5("outputs/test-star-formation-histories-adapative/galacticus.hdf5");
my $outputs                = $model            ->group  ('Outputs'               )                            ;
my $starFormationHistories = $model            ->group  ('starFormationHistories')                            ;
my $stellarPopulation      = $model            ->group  ('Parameters'            )->group('stellarPopulation');
(my $recycledFraction)     = $stellarPopulation->attrGet('recycledFraction'      )                            ;
foreach my $outputName ( sort($outputs->groups()) ) {
    print $outputName.":\n";
    # Read stellar masses of disk and bulge, along with node and tree indices.
    my $nodeData                 = $outputs ->group  ($outputName          )->group('nodeData');
    my $nodesMassStellarDisk     = $nodeData->dataset('diskMassStellar'    )->get  (          );
    my $nodesMassStellarSpheroid = $nodeData->dataset('spheroidMassStellar')->get  (          );
    # Read corresponding star formation histories.
    my $starFormationHistoryData = $starFormationHistories->group($outputName);
    my $sfhDiskStarFormationHistory     = $starFormationHistoryData->dataset('diskStarFormationHistory'    )->get();
    my $sfhSpheroidStarFormationHistory = $starFormationHistoryData->dataset('spheroidStarFormationHistory')->get();
    # Iterate over nodes, computing the integrated star formation history.
    my $nodesSFHIntegratedDisk     = pdl zeroes(nelem($nodesMassStellarDisk    ));
    my $nodesSFHIntegratedSpheroid = pdl zeroes(nelem($nodesMassStellarSpheroid));
    for(my $i=0;$i<nelem($nodesMassStellarDisk);++$i) {
	$nodesSFHIntegratedDisk    ->(($i)) .= $sfhDiskStarFormationHistory    (($i),:,:)->sum()*(1.0-$recycledFraction);
	$nodesSFHIntegratedSpheroid->(($i)) .= $sfhSpheroidStarFormationHistory(($i),:,:)->sum()*(1.0-$recycledFraction);
    }
    my $tolerance     = pdl 1.0e-3;
    my $statusDisk     = any(abs($nodesSFHIntegratedDisk    -$nodesMassStellarDisk    ) > $tolerance*$nodesMassStellarDisk    ) ? "FAILED" : "SUCCESS";
    my $statusSpheroid = any(abs($nodesSFHIntegratedSpheroid-$nodesMassStellarSpheroid) > $tolerance*$nodesMassStellarSpheroid) ? "FAILED" : "SUCCESS";
    print " -> ".$statusDisk    .": disk stellar mass\n"    ;
    print " -> ".$statusSpheroid.": spheroid stellar mass\n";
}

exit 0;
