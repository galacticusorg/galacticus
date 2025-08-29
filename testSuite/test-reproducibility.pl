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
                  values            => pdl ( 7.99799857410098e10 ),
     		  toleranceRelative => 4.0e-6
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
     		  values            => pdl ( 2.8814957e9 )          ,
     		  toleranceRelative => 1.0e-2
     	      }
     	     ]
     },
     {
     	 # Leaky box chemical evolution model - values below are computed from the analytic solution for this model.
     	 name           => "leakyBox"                                         ,
     	 parameters     => "testSuite/parameters/reproducibility/leakyBox.xml",
     	 outputFileName => "testSuite/outputs/reproducibility/leakyBox.hdf5"  ,
     	 assertions     =>
     	     [
     	      {
     		  name              => "gas mass"                   ,
     		  output            => 1                            ,
     		  property          => "diskMassGas"                ,
     		  values            => pdl ( 4.0762204e9  )         ,
     		  toleranceRelative => 1.1e-2
     	      },
     	      {
     		  name              => "stellar mass"               ,
     		  output            => 1                            ,
     		  property          => "diskMassStellar"            ,
     		  values            => pdl ( 3.5971417e10  )        ,
     		  toleranceRelative => 1.0e-2
     	      },
     	      {
     		  name              => "gas metals"                 ,
     		  output            => 1                            ,
     		  property          => "diskAbundancesGasMetals"    ,
     		  values            => pdl ( 2.03811e8   )          ,
     		  toleranceRelative => 1.0e-2
     	      },
     	      {
     		  name              => "stellar metals"             ,
     		  output            => 1                            ,
     		  property          => "diskAbundancesStellarMetals",
     		  values            => pdl ( 4.85624e8   )          ,
     		  toleranceRelative => 1.0e-2
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
		  values            => pdl ( 0.00360702914995244 ),
		  toleranceRelative => 1.4e-5
	      },
	      {
		  name              => "spheroid angular momentum"                                                                ,
		  output            => 1                                                                                          ,
		  expression        => "(%[spheroidRadius]*%[spheroidVelocity]*%[spheroidMassStellar])/%[spheroidAngularMomentum]",
		  values            => pdl ( 0.5 )                                                                                ,
		  toleranceRelative => 2.0e-4
	      },
	      {
		  name              => "rotation curve"                         ,
		  output            => 1                                        ,
		  expression        => "%[rotationCurve{0}]/%[spheroidVelocity]",
		  values            => pdl ( 1.0 )                              ,
		  toleranceRelative => 2.0e-4
	      },
	      {
		  # The adiabatic contraction calculation solves the following system:
		  #  Mₜ(r₀) f₀ r₀ - Mₜ(r₁) f₁ r₁ - Mᵦ r₁ = 0
		  # where:
		  #  * Mₜ(r) is the total mass of the node (assuming baryons trace dark matter) within radius r;
		  #  * r₀ is the initial radius in the dark matter halo;
		  #  * r₁ is the final radius in the dark matter halo;
		  #  * f₀ is the initial mass fraction distributed as the dark matter;
		  #  * f₁ is the final mass fraction distributed as the dark matter;
		  #  * Mᵦ is the mass of baryons within r₁.
		  #		  
		  # The second two terms are equal to (V₁ r₁)²/G where V₁ is the final circular velocity at radius r₁.  In an
		  # isothermal halo Mₜ(r) = Mᵥ r/rᵥ where Mᵥ is the virial mass of the halo, and rᵥ is the virial radius. Using
		  # G Mᵥ/rᵥ=Vᵥ² with Vᵥ being the virial velocity of the halo the system is reduced to:		  
		  #  f₀ (Vᵥ r₀)² - (V₁ r₁)² = 0,
		  # or,
		  #  √f₀ Vᵥ r₀/V₁/r₁ = 1.		  
		  # Finally, r₀ = rᵥ Mₜ(r₁)/Mᵥ, and Mₜ(r₁) f₁ = Vᵪ₁² r₁/G where Vᵪ₁² is the dark matter contrbution to the
		  # final rotation curve.  In this calculation, f₁=(0.3-0.05)/0.3=0.833333, and f₀=f₁+10¹⁰/10¹²=0.8433333, 
		  name              => "initial specific angular momentum",
		  output            => 1                                  ,
		  expression        =>
		      "+sqrt(0.84333333)"                                 .
		      "*%[darkMatterOnlyVelocityVirial]"                  .
		      "/%[spheroidVelocity]"                              .
		      "*("                                                .
		      "  +%[darkMatterOnlyRadiusVirial]"                  .
		      "  *("                                              .
		      "    +%[rotationCurve{1}]**2"                       .
		      "    *%[spheroidRadius]"                            .
		      "    /$gravitationalConstant"                       .
		      "    /0.83333333"                                   .
		      "   )"                                              .
		      "  /%[basicMass]"                                   .
		      " )".
		      "/%[spheroidRadius]"                                ,
		  values            => pdl ( 1.0 )                        ,
		  toleranceRelative => 3.0e-3
	      }
	     ]
     }
    );
    
# Run tests.
foreach my $test ( @tests ) {
    print "Running test \"".$test->{'name'}."\": Galacticus.exe ".$test->{'parameters'}."\n";
    system("cd ..; ./Galacticus.exe ".$test->{'parameters'});
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
		while ( my @matches = $assertion->{'expression'} =~ m/%\[([^%]+?)(\{(\d+)\})??\]/ ) {
		    my $propertyName  = $matches[0];
		    my $propertyIndex = $matches[2];
		    my $replacement = "\$properties{'$1'}";
		    $replacement .= "->slice('(".$propertyIndex.")','')"
			if ( defined($propertyIndex) );
		    $assertion->{'expression'} =~ s/%\[[^%]+\]/$replacement/;
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
