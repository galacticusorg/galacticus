#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;

# Run models that test reproducibility.
# Andrew Benson (27-June-2015)

# Make output directory.
system("mkdir -p outputs/reproducibility");

# Define constants.
my $gravitationalConstant = pdl 4.3011827419096073e-9;

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
	 name           => "adiabaticContraction"                                         ,
	 parameters     => "testSuite/parameters/reproducibility/adiabaticContraction.xml",
	 outputFileName => "testSuite/outputs/reproducibility/adiabaticContraction.hdf5"  ,
	 assertions     =>
	     [
	      {
		  name              => "spheroid radius"          ,
		  output            => 1                          ,
		  property          => "spheroidRadius"           ,
		  values            => pdl ( 0.00360653416232037 ),
		  toleranceRelative => 1.0e-12
	      },
	      {
		  name              => "spheroid angular momentum"                                                                ,
		  output            => 1                                                                                          ,
		  expression        => "(%[spheroidRadius]*%[spheroidVelocity]*%[spheroidMassStellar])/%[spheroidAngularMomentum]",
		  values            => pdl ( 0.5 )                                                                                ,
		  toleranceRelative => 4.0e-6
	      },
	      {
		  name              => "rotation curve"                                                        ,
		  output            => 1                                                                       ,
		  expression        => "%[rotationCurve:spheroidRadius:all:all:loaded:1.0]/%[spheroidVelocity]",
		  values            => pdl ( 1.0 )                                                             ,
		  toleranceRelative => 3.0e-6
	      },
	      {
		  # The adiabatic contraction calculation solves the following system:
		  #  Mt(ri)*fi*ri - Mt(rf)*ff*rf - Mb*rf = 0
		  # where:
		  #  * Mt(r) is the total mass of the node (assuming baryons trace dark matter) within radius r;
		  #  * ri is the initial radius in the dark matter halo;
		  #  * ri is the final radius in the dark matter halo;
		  #  * fi is the initial mass fraction distributed as the dark matter;
		  #  * ff is the final mass fraction distributed as the dark matter;
		  #  * Mb is the mass of baryons within rf.
		  #		  
		  # The second two terms are equal to (Vf rf)^2/G where Vf is the final circular velocity at radius rf.  In an
		  # isothermal halo Mt(r) = Mv*r/rv where Mv is the virial mass of the halo, and rv is the virial radius. Using
		  # G*Mv/rv=Vv^2 with Vv being the virial velocity of the halo the system is reduced to:		  
		  #  fi*(Vv*ri)^2 - (Vf*rf)^2 = 0,
		  # or,
		  #  sqrt(fi)*Vv*ri/Vf/rf = 1.		  
		  # Finally, ri = rv Mt(rf)/Mv, and Mt(rf)*ff = VDMf^2*rf/G where VDMf^2 is the dark matter contrbution to the
		  # final rotation curve.  In this calculation, ff=(0.3-0.05)/0.3=0.833333, and fi=ff+1.0e10/1.0e12=0.8433333, 
		  name              => "initial specific angular momentum"              ,
		  output            => 1                                                ,
		  expression        =>
		      "+sqrt(0.84333333)"                                               .
		      "*%[nodeVirialVelocity]"                                          .
		      "/%[spheroidVelocity]"                                            .
		      "*("                                                              .
		      "  +%[nodeVirialRadius]"                                          .
		      "  *("                                                            .
		      "    +%[rotationCurve:spheroidRadius:darkHalo:dark:loaded:1.0]**2".
		      "    *%[spheroidRadius]"                                          .
		      "    /$gravitationalConstant"                                     .
		      "    /0.83333333"                                                 .
		      "   )"                                                            .
		      "  /%[basicMass]"                                                 .
		      " )".
		      "/%[spheroidRadius]"                                              ,
		  values            => pdl ( 1.0 )                                      ,
		  toleranceRelative => 3.0e-3
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
	    my $values;
	    if ( exists($assertion->{'property'}) ) {
		my $property = $nodeData->dataset($assertion->{'property'});
		$values      = $property->get    (                        );
	    } elsif ( exists($assertion->{'expression'}) ) {
		my %propertyNames;
		while ( my @matches = $assertion->{'expression'} =~ m/%\[([^%]+)\]/ ) {
		    my $propertyName = $matches[0];
		    $assertion->{'expression'} =~ s/%\[([^%]+)\]/\$properties{'$1'}/;
		    ++$propertyNames{$propertyName};
		}
		my %properties;
		foreach ( keys(%propertyNames) ) {
		    $properties{$_} = $nodeData->dataset($_)->get();
		}
		$values = eval($assertion->{'expression'});
	    } else {
		die('assertion gives no property or expression to test')
	    }
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
