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
my $outputs                = $model            ->group  ('Outputs'               )                                  ;
my $starFormationHistories = $model            ->group  ('starFormationHistories')                                  ;
my $stellarPopulation      = $model            ->group  ('Parameters'            )->group('stellarPopulationMethod');
(my $recycledFraction)     = $stellarPopulation->attrGet('recycledFraction'      )                                  ;
foreach my $outputName ( sort($outputs->groups()) ) {
    print $outputName.":\n";
    # Read stellar masses of disk and bulge, along with node and tree indices.
    my $nodeData                 = $outputs ->group  ($outputName          )->group('nodeData');
    my $nodesNodeIndex           = $nodeData->dataset('nodeIndex'          )->get  (          );
    my $nodesTreeIndex           = $nodeData->dataset('mergerTreeIndex'    )->get  (          );
    my $nodesMassStellarDisk     = $nodeData->dataset('diskMassStellar'    )->get  (          );
    my $nodesMassStellarSpheroid = $nodeData->dataset('spheroidMassStellar')->get  (          );
    # Read corresponding star formation histories.
    my $starFormationHistoryData = $starFormationHistories->group($outputName);
    my $sfhDiskNodeIndex                = $starFormationHistoryData->dataset('diskNodeIndex'               )->get();
    my $sfhDiskTreeIndex                = $starFormationHistoryData->dataset('diskTreeIndex'               )->get();
    my $sfhDiskStarFormationHistory     = $starFormationHistoryData->dataset('diskStarFormationHistory'    )->get();
    my $sfhSpheroidNodeIndex            = $starFormationHistoryData->dataset('spheroidNodeIndex'           )->get();
    my $sfhSpheroidTreeIndex            = $starFormationHistoryData->dataset('spheroidTreeIndex'           )->get();
    my $sfhSpheroidStarFormationHistory = $starFormationHistoryData->dataset('spheroidStarFormationHistory')->get();
    # Iterate over nodes, computing the integrated star formation history.
    my $nodesSFHIntegratedDisk     = pdl zeroes(nelem($nodesNodeIndex));
    my $nodesSFHIntegratedSpheroid = pdl zeroes(nelem($nodesNodeIndex));
    for(my $i=0;$i<nelem($nodesNodeIndex);++$i) {
	my $disk     = which(($sfhDiskNodeIndex     == $nodesNodeIndex->(($i))) & ($sfhDiskTreeIndex     == $nodesTreeIndex->(($i))));
	my $spheroid = which(($sfhSpheroidNodeIndex == $nodesNodeIndex->(($i))) & ($sfhSpheroidTreeIndex == $nodesTreeIndex->(($i))));
	if ( nelem($disk) == 1 ) {
	    $nodesSFHIntegratedDisk    ->(($i)) .= $sfhDiskStarFormationHistory    (($disk    ->((0))),:,:)->sum()*(1.0-$recycledFraction);
	} else {
	    die("FAILED: no matching disk star formation history found");
	}
	if ( nelem($spheroid) == 1 ) {
	    $nodesSFHIntegratedSpheroid->(($i)) .= $sfhSpheroidStarFormationHistory(($spheroid->((0))),:,:)->sum()*(1.0-$recycledFraction);
	} else {
	    die("FAILED: no matching spheroid star formation history found");
	}
    }
    my $tolerance     = pdl 1.0e-3;
    my $statusDisk     = any(abs($nodesSFHIntegratedDisk    -$nodesMassStellarDisk    ) > $tolerance*$nodesMassStellarDisk    ) ? "FAILED" : "SUCCESS";
    my $statusSpheroid = any(abs($nodesSFHIntegratedSpheroid-$nodesMassStellarSpheroid) > $tolerance*$nodesMassStellarSpheroid) ? "FAILED" : "SUCCESS";
    print " -> ".$statusDisk    .": disk stellar mass\n"    ;
    print " -> ".$statusSpheroid.": spheroid stellar mass\n";
}

exit 0;
