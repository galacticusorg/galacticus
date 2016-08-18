#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/../perl';
use PDL;
use PDL::NiceSlice;
use Galacticus::HDF5;

# Run a standard hot halo model to test mass conservation.
# Andrew Benson (23-May-2016)

# Run the model.
system("cd ..; scripts/aux/launch.pl testSuite/parameters/test-mass-conservation-standard.xml");

# Extract masses and check they add up to the expected value.
my $galacticus;
$galacticus->{'file' } = "outputs/test-mass-conservation-standard/galacticus_0:1/galacticus.hdf5";
$galacticus->{'store'} = 0;
$galacticus->{'tree' } = "all";
&Galacticus::HDF5::Get_Parameters($galacticus    );
&Galacticus::HDF5::Count_Trees   ($galacticus    );
&Galacticus::HDF5::Select_Output ($galacticus,0.0);
&Galacticus::HDF5::Get_Dataset($galacticus,['mergerTreeWeight','diskMassStellar','diskMassGas','spheroidMassStellar','spheroidMassGas','hotHaloMass','hotHaloOutflowedMass','nodeIsIsolated','basicMass','hotHaloUnaccretedMass','mergerTreeIndex']);
my $properties = $galacticus->{'dataSets'  };
my $parameters = $galacticus->{'parameters'};
# Find centrals.
my $centrals   = which($properties->{'nodeIsIsolated'} == 1);
my $satellites = which($properties->{'nodeIsIsolated'} == 0);
# Find masses in satellites
my $massSatellites = pdl zeroes(nelem($centrals));
for(my $i=0;$i<nelem($centrals);++$i) {
    my $inTree = which($properties->{'mergerTreeIndex'}->($satellites) == $properties->{'mergerTreeIndex'}->($centrals)->(($i)));
    $massSatellites->(($i)) .=
	    (
	     +$properties->{'diskMassStellar'      }->($satellites)->($inTree)->sum()
	     +$properties->{'diskMassGas'          }->($satellites)->($inTree)->sum()	
	     +$properties->{'spheroidMassStellar'  }->($satellites)->($inTree)->sum()
	     +$properties->{'spheroidMassGas'      }->($satellites)->($inTree)->sum()	
	     +$properties->{'hotHaloMass'          }->($satellites)->($inTree)->sum()
	     +$properties->{'hotHaloOutflowedMass' }->($satellites)->($inTree)->sum()
	     +$properties->{'hotHaloUnaccretedMass'}->($satellites)->($inTree)->sum()
	    )
	    /(
		+$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
		/$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
	    )
	    /$properties->{'basicMass'}->($centrals)->(($i));
}
$properties->{'diskMassStellar'      } /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'diskMassGas'          } /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'spheroidMassStellar'  } /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'spheroidMassGas'      } /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'hotHaloMass'          } /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'hotHaloOutflowedMass' } /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'hotHaloUnaccretedMass'} /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
my $massTotal =
    +$properties    ->{'diskMassStellar'      }->($centrals)
    +$properties    ->{'diskMassGas'          }->($centrals)
    +$properties    ->{'spheroidMassStellar'  }->($centrals)
    +$properties    ->{'spheroidMassGas'      }->($centrals)
    +$properties    ->{'hotHaloMass'          }->($centrals)
    +$properties    ->{'hotHaloOutflowedMass' }->($centrals)
    +$properties    ->{'hotHaloUnaccretedMass'}->($centrals)
    +$massSatellites;
# Check that all masses are unity.
if ( any(abs($massTotal-1.0) > 1.0e-4) ) {
    print "FAILED: mass conservation failure -> ".($massTotal-1.0)."\n";
} else {
    print "SUCCESS: mass conservation\n";
}
exit;
