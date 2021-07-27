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

# Check calculations of stellar mass-weighted ages.
# Andrew Benson (23-July-2021)

# These models used fixed timescales for star formation, and have no inflow or outflow, resulting in analytically solvable
# (exponential) star formation histories. The mean age can therefore be computed directly.

{
    # Simple model - a disk of gas evolves in isolation forming stars.
    # Compute expectations.
    my $timescaleDisk    = pdl  1.0   ;
    my $timeStart        = pdl  1.0   ;
    my $timeEnd          = pdl 13.8   ;
    my $massDiskStart    = pdl  1.0e11;
    my $massIntegral     = $massDiskStart*(1.0-exp(-($timeEnd-$timeStart)/$timescaleDisk));
    my $massTimeIntegral = $massDiskStart*($timeStart+$timescaleDisk-($timeEnd+$timescaleDisk)*exp(-($timeEnd-$timeStart)/$timescaleDisk));
    my $ageDiskTarget    = $timeEnd-$massTimeIntegral/$massIntegral;
    # Run the simple model.
    &System::Redirect::tofile("cd ..; ./Galacticus.exe testSuite/parameters/stellarMassWeightedAgesSimple.xml","outputs/stellarMassWeightedAgesSimple.log");
    unless ( $? == 0 ) {
	print "FAILED: simple model run:\n";
	system("cat outputs/stellarMassWeightedAgesSimple.log");
    } else {
	print "SUCCESS: simple model run\n";
    }
    # Read the model data and check for consistency.
    my $model     = new PDL::IO::HDF5("outputs/stellarMassWeightedAgesSimple.hdf5");
    my $outputs   = $model   ->group  ('Outputs'                   )       ;
    my $output    = $outputs ->group  ('Output1'                   )       ;
    my $nodeData  = $output  ->group  ('nodeData'                  )       ;
    my $indexNode = $nodeData->dataset('nodeIndex'                 )->get();
    my $ageDisk   = $nodeData->dataset('diskAgeStellarMassWeighted')->get();
    if ( abs($ageDisk->((0))-$ageDiskTarget) < 1.0e-3 ) {
	print "SUCCESS: simple model age\n";
    } else {
	print "FAILED: simple model age: ".$ageDisk->((0))." ≇ ".$ageDiskTarget."\n";
    }
}

{
    # Merging model - two disks of gas evolve in isolation up to t=6 Gyr, then merge, forming a spheroid,
    # Compute expectations.
    my $timescaleDisk      = pdl  1.00   ;
    my $timescaleSpheroid  = pdl  0.75   ;
    my $timeStartCentral   = pdl  1.00   ;
    my $timeStartSatellite = pdl  3.00   ;
    my $timeMerge          = pdl  6.00   ;
    my $timeEnd            = pdl 13.80   ;
    my $massDiskStart      = pdl  1.00e11;
    my $massIntegralCentralPreMerge       = $massDiskStart*(1.0-exp(-($timeMerge-$timeStartCentral)/$timescaleDisk));
    my $massTimeIntegralCentralPreMerge   = $massDiskStart*($timeStartCentral+$timescaleDisk-($timeMerge+$timescaleDisk)*exp(-($timeMerge-$timeStartCentral)/$timescaleDisk));
    my $massGasCentralMerge               = $massDiskStart*exp(-($timeMerge-$timeStartCentral)/$timescaleDisk);
    my $massIntegralSatellitePreMerge     = $massDiskStart*(1.0-exp(-($timeMerge-$timeStartSatellite)/$timescaleDisk));
    my $massTimeIntegralSatellitePreMerge = $massDiskStart*($timeStartSatellite+$timescaleDisk-($timeMerge+$timescaleDisk)*exp(-($timeMerge-$timeStartSatellite)/$timescaleDisk));
    my $massGasSatelliteMerge             = $massDiskStart*exp(-($timeMerge-$timeStartSatellite)/$timescaleDisk);
    my $massSpheroidStart                 = $massGasCentralMerge+$massGasSatelliteMerge;
    my $massIntegralCentralPostMerge      = $massSpheroidStart*(1.0-exp(-($timeEnd-$timeMerge)/$timescaleSpheroid));
    my $massTimeIntegralCentralPostMerge  = $massSpheroidStart*($timeMerge+$timescaleSpheroid-($timeEnd+$timescaleSpheroid)*exp(-($timeEnd-$timeMerge)/$timescaleSpheroid));
    my $massIntegral                      = $massIntegralCentralPreMerge+$massIntegralSatellitePreMerge+$massIntegralCentralPostMerge;
    my $massTimeIntegral                  = $massTimeIntegralCentralPreMerge+$massTimeIntegralSatellitePreMerge+$massTimeIntegralCentralPostMerge;
    my $ageSpheroidTarget                 = $timeEnd-$massTimeIntegral/$massIntegral;    
    # Run the merging model.
    &System::Redirect::tofile("cd ..; ./Galacticus.exe testSuite/parameters/stellarMassWeightedAgesMerging.xml","outputs/stellarMassWeightedAgesMerging.log");
    unless ( $? == 0 ) {
	print "FAILED: merging model run:\n";
	system("cat outputs/stellarMassWeightedAgesMerging.log");
    } else {
	print "SUCCESS: merging model run\n";
    }
    # Read the model data and check for consistency.
    my $model       = new PDL::IO::HDF5("outputs/stellarMassWeightedAgesMerging.hdf5");
    my $outputs     = $model   ->group  ('Outputs'                       )       ;
    my $output      = $outputs ->group  ('Output1'                       )       ;
    my $nodeData    = $output  ->group  ('nodeData'                      )       ;
    my $indexNode   = $nodeData->dataset('nodeIndex'                     )->get();
    my $ageSpheroid = $nodeData->dataset('spheroidAgeStellarMassWeighted')->get();
    if ( abs($ageSpheroid->((0))-$ageSpheroidTarget) < 1.0e-3 ) {
	print "SUCCESS: merging model age\n";
    } else {
	print "FAILED: merging model age: ".$ageSpheroid->((0))." ≇ ".$ageSpheroidTarget."\n";
    }
}

exit 0;
