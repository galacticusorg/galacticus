#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use Galacticus::HDF5;
use Data::Dumper;

# Run a standard hot halo model to test mass conservation.
# Andrew Benson (23-May-2016)

# Run the model.
system("cd ..; scripts/aux/launch.pl testSuite/parameters/test-mass-conservation-standard.xml");

# Check for failed models.
system('grep -q -i -e fatal -e "Galacticus experienced" outputs/test-mass-conservation-standard/galacticus_*/galacticus.log');
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i -e fatal -e "Galacticus experienced" outputs/test-mass-conservation-standard/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS: model run\n";
}

# Extract masses and check they add up to the expected value.
my $galacticus;
$galacticus->{'file' } = "outputs/test-mass-conservation-standard/galacticus_0:1/galacticus.hdf5";
$galacticus->{'store'} = 0;
$galacticus->{'tree' } = "all";
&Galacticus::HDF5::Get_Parameters($galacticus    );
&Galacticus::HDF5::Count_Trees   ($galacticus    );
&Galacticus::HDF5::Select_Output ($galacticus,0.0);
&Galacticus::HDF5::Get_Dataset($galacticus,['mergerTreeWeight','blackHoleMass','diskMassStellar','diskMassGas','spheroidMassStellar','spheroidMassGas','hotHaloMass','hotHaloOutflowedMass','nodeIsIsolated','basicMass','hotHaloUnaccretedMass','mergerTreeIndex']);
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
	     +$properties->{'blackHoleMass'        }->($satellites)->($inTree)->sum()
	    )
	    /(
		+$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
		/$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
	    )
	    /$properties->{'basicMass'}->($centrals)->(($i));
}
$properties->{'diskMassStellar'      } /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'diskMassGas'          } /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'spheroidMassStellar'  } /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'spheroidMassGas'      } /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'hotHaloMass'          } /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'hotHaloOutflowedMass' } /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'hotHaloUnaccretedMass'} /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMass'};
$properties->{'blackHoleMass'        } /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
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
    +$properties    ->{'blackHoleMass'        }->($centrals)
    +$massSatellites;
# Check that all masses are unity.
if ( any(abs($massTotal-1.0) > 1.0e-3) ) {
    print "FAILED: mass conservation [standard] failure -> ".($massTotal-1.0)."\n";
} else {
    print "SUCCESS: mass conservation\n";
}
exit;
