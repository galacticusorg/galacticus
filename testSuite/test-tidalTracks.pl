#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test that subhalo tidal track evolution by validating against the fitting function of Errani & Navarro (2021;
# https://ui.adsabs.harvard.edu/abs/2021MNRAS.505...18E).
# Andrew Benson (17-December-2021)

# Make output directory.
system("mkdir -p outputs/");

# Run the tidal tracks model.
system("cd ..; ./Galacticus.exe testSuite/parameters/tidalTracks.xml");
unless ( $? == 0 ) {
    print "FAIL: tidal track model failed to run\n";
    exit;
}

my $velocityMaximumInitial                                                ;
my $radiusMaximumInitial                                                  ;
my $offsetMaximum          = 0.0                                         ;
my $model                  = new PDL::IO::HDF5("outputs/tidalTracks.hdf5");
my $outputs                = $model->group('Outputs')                     ;
for(my $i=1;$i<602;++$i) {
    my $output = $outputs->group('Output'  .$i);
    my $nodes  = $output ->group('nodeData'   );
    my $properties;
    $properties->{$_} = $nodes->dataset($_)->get()
	foreach ( 'darkMatterProfileDMORadiusVelocityMaximum', 'darkMatterProfileDMOVelocityMaximum' );
    # Validate number of nodes.
    unless ( nelem($properties->{'satelliteBoundMass'}) == 1 ) {
	print "FAILED: expected only 1 node in the output\n";
	exit;
    }
    # Capture initial values.
    $velocityMaximumInitial = $properties->{'darkMatterProfileDMOVelocityMaximum'      }->((0))
	unless ( defined($velocityMaximumInitial) );
    $radiusMaximumInitial   = $properties->{'darkMatterProfileDMORadiusVelocityMaximum'}->((0))
	unless ( defined($radiusMaximumInitial  ) );
    # Evaluate fit (equation 5 of Errani & Navarro, 20201; https://ui.adsabs.harvard.edu/abs/2021MNRAS.505...18E).
    my $alpha                       = pdl 0.40;
    my $beta                        = pdl 0.65;
    my $radiusMaximumScaleFree      = $properties->{'darkMatterProfileDMORadiusVelocityMaximum'}->((0))/$radiusMaximumInitial  ;
    my $velocityMaximumScaleFree    = $properties->{'darkMatterProfileDMOVelocityMaximum'      }->((0))/$velocityMaximumInitial;
    my $velocityMaximumScaleFreeFit = 2.0**$alpha*$radiusMaximumScaleFree**$beta/(1.0+$radiusMaximumScaleFree**2)**$alpha;
    my $offset                      = abs(log10($velocityMaximumScaleFree/$velocityMaximumScaleFreeFit));
    $offsetMaximum                  = $offset
	if ( $offset > $offsetMaximum );
}

my $status = $offsetMaximum < 0.01 ? "SUCCESS" : "FAILED";
print $status.": subhalo tidal tracks\n";

exit;
