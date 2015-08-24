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
		for(my $i=0;$i<nelem($difference);++$i) {
		    print "\t".$values->(($i))."\t".$assertion->{'values'}->(($i))."\t".$difference->(($i))."\t".$allowedError->(($i))."\n";
		}
	    }
	}
    } else {
	print "FAIL: reproducibility test '".$test->{'name'}."' model failed to run\n";
    }
}

exit;
