#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;

# Run models that test that subhalo tidal track evolution by validating against the fitting function of Errani & Navarro (2021;
# https://ui.adsabs.harvard.edu/abs/2021MNRAS.505...18E).
# Andrew Benson (17-December-2021)

# Make output directory.
system("mkdir -p outputs/");

# Iterate over models to run.
my @testCases =
    (
     {
      	 label     => "nonMonotonic",
      	 gamma     => 1.00          ,
      	 fitMetric => 0.03
     },
     {
	 label     => "monotonic",
	 gamma     => 1.0000     ,
	 fitMetric => 0.0062
     },
     {
      	 label     => "monotonic",
     	 gamma     => 0.5000     ,
	 fitMetric => 0.0056
     },
     {
      	 label     => "monotonic",
      	 gamma     => 0.000      ,
	 fitMetric => 0.028
     }
    );
foreach my $testCase ( @testCases ) {

    # Run the tidal tracks model.
    system("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/tidalTracks".ucfirst($testCase->{'label'})."_gamma".sprintf("%3.1f",$testCase->{'gamma'}).".xml");
    unless ( $? == 0 ) {
	print "FAIL: tidal track model '".$testCase->{'label'}." gamma=".sprintf("%3.1f",$testCase->{'gamma'})."' failed to run\n";
	exit;
   }

    my $velocityMaximumInitial                                                                                                                                ;
    my $radiusMaximumInitial                                                                                                                                  ;
    my $offsetMaximum          = 0.0                                                                                                                          ;
    my $model                  = new PDL::IO::HDF5("outputs/tidalTracks".ucfirst($testCase->{'label'})."_gamma".sprintf("%3.1f",$testCase->{'gamma'}).".hdf5");
    my $outputs                = $model->group('Outputs')                                                                                                     ;
    my $track;
    $track->{'model'}->{'t'} = pdl [];
    $track->{'model'}->{'r'} = pdl [];
    $track->{'model'}->{'v'} = pdl [];
    foreach my $outputName ( $outputs->groups() ) {
	my $output   = $outputs->group  ($outputName );
	my $nodeData = $output ->group  ('nodeData'  );
	my @propertyNames = ( 'time', 'nodeIsIsolated', 'darkMatterProfileDMOVelocityMaximum', 'darkMatterProfileDMORadiusVelocityMaximum' );
	my $snapshot;
	$snapshot->{$_} = $nodeData->dataset($_)->get()
	    foreach ( @propertyNames );
	my $selection = which($snapshot->{'nodeIsIsolated'} == 0);
	$track->{'model'}->{'t'} = $track->{'model'}->{'t'}->append($snapshot->{'time'                                     }->($selection));
	$track->{'model'}->{'r'} = $track->{'model'}->{'r'}->append($snapshot->{'darkMatterProfileDMORadiusVelocityMaximum'}->($selection));
	$track->{'model'}->{'v'} = $track->{'model'}->{'v'}->append($snapshot->{'darkMatterProfileDMOVelocityMaximum'      }->($selection));
    }
    my $order = $track->{'model'}->{'t'}->qsorti();
    my $r0    = $track->{'model'}->{'r'}->($order)->((0))->copy();
    my $v0    = $track->{'model'}->{'v'}->($order)->((0))->copy();
    $track->{'model'}->{'r'} /= $r0;
    $track->{'model'}->{'v'} /= $v0;
    # Construct N-body tidal track.
    if ( $testCase->{'gamma'} == 1.0 ) {
	# Tidal track from Errani & Navarro (2021).
	my $la    = pdl sequence(1000)/999.0*6.0-6.0;
	my $a     = 10.0**$la;
	my $alpha = pdl 0.40;
	my $beta  = pdl 0.65;
	$track->{'nBody'}->{'r'} = $a;
	$track->{'nBody'}->{'v'} = 2.0**$alpha*$a**$beta/(1.0+$a**2)**$alpha;
    } else {
	# Tidal tracks from Penarrubia et al. (2010).
	my $mur ;
	my $etar;
	my $muv ;
	my $etav;
	if ( $testCase->{'gamma'} == 0.0 ) {
	    $mur   = pdl -1.30;
	    $etar  = pdl +0.05;
	    $muv   = pdl +0.40;
	    $etav  = pdl +0.37;
	} elsif ( $testCase->{'gamma'} == 0.5 ) {
	    $mur   = pdl -0.40;
	    $etar  = pdl +0.27;
	    $muv   = pdl +0.40;
	    $etav  = pdl +0.35;
	} elsif ( $testCase->{'gamma'} == 1.5 ) {
	    $mur   = pdl +0.00;
	    $etar  = pdl +0.48;
	    $muv   = pdl +0.40;
	    $etav  = pdl +0.24;
	} else {
	    die("unknown gamma");
	}
	my $la    = pdl sequence(1000)/999.0*6.0-6.0;
	my $a     = 10.0**$la;
	my $r     = 2.0**$mur*$a**$etar/(1.0+$a)**$mur;
	my $v     = 2.0**$muv*$a**$etav/(1.0+$a)**$muv;
	my $order = $r->qsorti();
	$track->{'nBody'}->{'r'} = $r->($order);
	$track->{'nBody'}->{'v'} = $v->($order);
    }
    my $selection = which(($track->{'model'}->{'v'} > 0.1) & ($track->{'model'}->{'r'} > 0.1));
    my $fitMetric = 0.0;
    for(my $i=0;$i<nelem($selection);++$i) {
	my $distanceSquared =
	    +($track->{'nBody'}->{'r'}->log10()-$track->{'model'}->{'r'}->($selection)->(($i))->log10())**2
	    +($track->{'nBody'}->{'v'}->log10()-$track->{'model'}->{'v'}->($selection)->(($i))->log10())**2;
	$fitMetric += $distanceSquared->minimum();
    }
 
    my $status = $fitMetric < $testCase->{'fitMetric'} ? "SUCCESS" : "FAILED";
    print $status.": subhalo tidal tracks '".$testCase->{'label'}." gamma=".sprintf("%3.1f",$testCase->{'gamma'})."'\n";

}

exit;
