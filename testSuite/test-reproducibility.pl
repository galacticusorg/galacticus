#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test reproducibility.
# Andrew Benson (27-June-2015)

# Build Galacticus.
system("cd ..; make Galacticus.exe");

# Make output directory.
system("mkdir outputs/reproducibility");

# Define tests.
my @tests =
    (
     {
	 name           => "cooling"                                         ,
	 parameters     => "testSuite/parameters/reproducibility/cooling.xml",
	 outputFileName => "testSuite/outputs/reproducibility/cooling.hdf5"  ,
	 assertions     =>
	     [
	      {
		  name              => "hot halo mass"            ,
		  output            => 1                          ,
		  property          => "hotHaloMass"              ,
		  values            => pdl ( 7.97045921944598e10 ),
		  toleranceRelative => 1.0e-12
	      }
	     ]
     },
     {
	 # Closed box chemical evolution model - values below are computed from the analytic solution for this model.
	 name           => "closedBox"                                         ,
	 parameters     => "testSuite/parameters/reproducibility/closedBox.xml",
	 outputFileName => "testSuite/outputs/reproducibility/closedBox.hdf5"  ,
	 assertions     =>
	     [
	      {
		  name              => "gas mass"                   ,
		  output            => 1                            ,
		  property          => "diskMassGas"                ,
		  values            => pdl ( 9.0717953e9  )         ,
		  toleranceRelative => 1.0e-2
	      },
	      {
		  name              => "stellar mass"               ,
		  output            => 1                            ,
		  property          => "diskMassStellar"            ,
		  values            => pdl ( 9.0928205e10  )        ,
		  toleranceRelative => 1.0e-2
	      },
	      {
		  name              => "gas metals"                 ,
		  output            => 1                            ,
		  property          => "diskAbundancesGasMetals"    ,
		  values            => pdl ( 9.0717953e8 )          ,
		  toleranceRelative => 1.0e-2
	      },
	      {
		  name              => "stellar metals"             ,
		  output            => 1                            ,
		  property          => "diskAbundancesStellarMetals",
		  values            => pdl ( 2.8814957E9 )          ,
		  toleranceRelative => 1.0e-2
	      }
	     ]
     }
    );

# Run tests.
foreach my $test ( @tests ) {
    system("cd ..; Galacticus.exe ".$test->{'parameters'});
    if ( $? == 0 ) {
	my $model   = new PDL::IO::HDF5("../".$test->{'outputFileName'});
	my $outputs = $model->group('Outputs');
	foreach my $assertion ( @{$test->{'assertions'}} ) {
	    my $output       = $outputs ->group  ('Output'  .$assertion->{'output'  });
	    my $nodeData     = $output  ->group  ('nodeData'                         );
	    my $property     = $nodeData->dataset(           $assertion->{'property'});
	    my $values       = $property->get    (                                   );
	    my $difference   = abs($values-$assertion->{'values'});
	    my $allowedError = $assertion->{'toleranceRelative'}*abs($values);
	    if ( all($difference < $allowedError) ) {
		print "SUCCESS: assertion '".$assertion->{'name'}."' of reproducibility test '".$test->{'name'}."' passed\n";
	    } else {
		print    "FAIL: assertion '".$assertion->{'name'}."' of reproducibility test '".$test->{'name'}."' failed\n";
	    }
	}
    } else {
	print "FAIL: reproducibility test '".$test->{'name'}."' model failed to run\n";
    }
}

exit;
