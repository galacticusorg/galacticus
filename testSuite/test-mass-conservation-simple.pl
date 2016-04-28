#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "../";
 $ENV{"GALACTICUS_ROOT_V094"} = getcwd()."/../";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
require Galacticus::HDF5;

# Run a simple model to test mass conservation.
# Andrew Benson (28-April-2016)

# Run the model.
system("cd ..; scripts/aux/launch.pl testSuite/test-mass-conservation.xml");

# Extract masses and check they add up to the expected value.
my $galacticus;
$galacticus->{'file' } = "outputs/test-mass-conservation/galacticus_0:1/galacticus.hdf5";
$galacticus->{'store'} = 0;
$galacticus->{'tree' } = "all";
&HDF5::Get_Parameters($galacticus    );
&HDF5::Count_Trees   ($galacticus    );
&HDF5::Select_Output ($galacticus,0.0);
&HDF5::Get_Dataset($galacticus,['mergerTreeWeight','diskMassStellar','diskMassGas','hotHaloMass','hotHaloOutflowedMass','nodeIsIsolated','basicMassBertschinger','hotHaloUnaccretedMass','mergerTreeIndex']);
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
	     +$properties->{'diskMassStellar'}->($satellites)->($inTree)->sum()
	     +$properties->{'diskMassGas'}->($satellites)->($inTree)->sum()	
	     +$properties->{'hotHaloMass'}->($satellites)->($inTree)->sum()
	     +$properties->{'hotHaloOutflowedMass'}->($satellites)->($inTree)->sum()
	     +$properties->{'hotHaloUnaccretedMass'}->($satellites)->($inTree)->sum()
	    )
	    /(
		+$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
		/$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
	    )
	    /$properties->{'basicMassBertschinger'}->($centrals)->(($i));
}
$properties->{'diskMassStellar'} /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'};
$properties->{'diskMassGas'} /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'  };
$properties->{'hotHaloMass'} /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'  };
$properties->{'hotHaloOutflowedMass'} /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'  };
$properties->{'hotHaloUnaccretedMass'} /= 
    (
     +$parameters->{'cosmologyParametersMethod'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParametersMethod'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'  };
my $massTotal =
    +$properties    ->{'diskMassStellar'      }->($centrals)
    +$properties    ->{'diskMassGas'          }->($centrals)
    +$properties    ->{'hotHaloMass'          }->($centrals)
    +$properties    ->{'hotHaloOutflowedMass' }->($centrals)
    +$properties    ->{'hotHaloUnaccretedMass'}->($centrals)
    +$massSatellites;
# Check that all masses are unity.
if ( any(abs($massTotal-1.0) > 1.0e-6) ) {
    print "FAILED: mass conservation failure -> ".($massTotal-1.0)."\n";
} else {
    print "SUCCESS: mass conservation\n";
}
exit;
