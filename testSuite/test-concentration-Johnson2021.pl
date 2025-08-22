#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Stats::Basic;
use System::Redirect;
use File::Slurp qw(slurp);

# Run models that check that the results of Johnson, Benson & Grin (2021; https://ui.adsabs.harvard.edu/abs/2021ApJ...908...33J)
# applying their concentration model to merger trees can be reproduced.
# Andrew Benson (23-September-2024)

# Make output directory.
system("mkdir -p outputs/");

# Specify the tests to run. Mean and tolerance targets are taken from Table 1 of Johnson, Benson & Grin (2021).
my @tests =
    (
     {
	 # Model with no subsampling.
	 suffix           => "",
	 mean             => 1.084,
	 meanTolerance    => 0.070,
	 scatter          => 0.134,
	 scatterTolerance => 0.040
     },
     {
	 # Model with subsampling below 10⁸M☉.
	 suffix           => "Subsample1e8",
	 mean             => 1.084,
	 meanTolerance    => 0.070,
	 scatter          => 0.134,
	 scatterTolerance => 0.040
     },
     {
	 # Model with subsampling below 10⁷M☉.
	 suffix           => "Subsample1e7",
	 mean             => 1.084,
	 meanTolerance    => 0.070,
	 scatter          => 0.134,
	 scatterTolerance => 0.040
     }
    );

# Iterate over tests.
foreach my $test ( @tests ) {

    # Run the model.
    system("cd ..; ./Galacticus.exe testSuite/parameters/concentrationDistributionJohnson2021".ucfirst($test->{'suffix'}).".xml");
    unless ( $? == 0 ) {
	print "FAILED: Johnson2021 ".$test->{'suffix'}." concentration model failed to run\n";
	exit;
    }
    
    # Select halos in the mass range used for the results in Table 2 of Benson, Ludlow, & Cole (2019), and compute the mean and
    # scatter in concentration.
    my $model          = new PDL::IO::HDF5("outputs/concentrationDistributionJohnson2021".ucfirst($test->{'suffix'}).".hdf5");
    my $outputs        = $model   ->group  ('Outputs'                )       ;
    my $output         = $outputs ->group  ('Output1'                )       ;
    my $nodeData       = $output  ->group  ('nodeData'               )       ;
    my $nodeIsIsolated = $nodeData->dataset('nodeIsIsolated'         )->get();
    my $massHalo       = $nodeData->dataset('massHaloEnclosedCurrent')->get();
    my $concentration  = $nodeData->dataset('concentration'          )->get();
    # Compute N-body measurement uncertainties to both halo masses and concentrations.
    my $massParticle             = pdl 1.6e5; # COCO simulation particle mass.
    my $countParticles           = $massHalo/$massParticle;
    my $alpha                    = -0.20+1.46*log10($concentration)-0.25*log10($concentration)**2;
    my $b                        = -0.54;
    my $concentrationUncertainty = 10.0**($alpha+$b*log10($countParticles)); # Concentration uncertainty model from equation 5 of Benson, Ludlow, & Cole (2019).
    my $massUncertainty          = 0.135*(1000.0/$countParticles)**(1.0/3.0); # Mass uncertainty model from Trenti et al. (2010; http://adsabs.harvard.edu/abs/2012ApJ...747..100S).
    # Construct perturbed halo masses.
    my $massMeasured = $massHalo*exp($massUncertainty*grandom(nelem($massHalo)));
    # Select halos for inclusion.
    my $selection = which(($nodeIsIsolated == 1) & ($massMeasured->log10() >= 9.402) & ($massMeasured->log10() < 9.902));
    # Construct perturbed concentrations.
    my $concentrationsMeasured = $concentration->log10()->($selection)+$concentrationUncertainty->($selection)*grandom(nelem($selection));
    # Evaluate the mean and scatter and check these match the target values.
    my $mean    = $concentrationsMeasured->average();
    my $scatter = $concentrationsMeasured->stdv();
    print "Mean   : ".sprintf("%5.3f",$mean   )." (target = ".$test->{'mean'   }." ± ".$test->{'meanTolerance'   }.")\n";
    print "Scatter: ".sprintf("%5.3f",$scatter)." (target = ".$test->{'scatter'}." ± ".$test->{'scatterTolerance'}.")\n";
    my $status  =
	(
	 abs($mean   -$test->{'mean'   }) < $test->{'meanTolerance'   }
	 &&
	 abs($scatter-$test->{'scatter'}) < $test->{'scatterTolerance'}
	)
	? "SUCCESS"
	: "FAILED";
    print $status.": Johnson2021 ".$test->{'suffix'}." concentration\n";
}

exit;
