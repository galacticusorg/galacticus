#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use Galacticus::HDF5;
use Galacticus::Options;

# Run a simple model to test mass conservation.
# Andrew Benson (28-April-2016)

# Read in any configuration options.
my $config;
if ( -e "galacticusConfig.xml" ) {
    my $xml = new XML::Simple;
    $config = $xml->XMLin("galacticusConfig.xml");
}

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Set default options.
my %options =
    (
     'pbsJobMaximum' => exists($queueConfig->{'jobMaximum'}) ? $queueConfig->{'jobMaximum'} : 100,
    );

# Get any command line options.
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Run the model.
system("cd ..; scripts/aux/launch.pl testSuite/parameters/test-mass-conservation-simple.xml ".join(" ",map {"--".$_." ".$options{$_}} keys(%options)));

# Check for failed models.
system("grep -q -i fatal outputs/test-mass-conservation-simple/galacticus_*/galacticus.log");
if ( $? == 0 ) {
    # Failures were found. Output their reports.
    my @failures = split(" ",`grep -l -i fatal outputs/test-mass-conservation-simple/galacticus_*/galacticus.log`);
    foreach my $failure ( @failures ) {
	print "FAILED: log from ".$failure.":\n";
	system("cat ".$failure);
    }
} else {
    print "SUCCESS: model run\n";
}

# Extract masses and check they add up to the expected value.
my $galacticus;
$galacticus->{'file' } = "outputs/test-mass-conservation-simple/galacticus_0:1/galacticus.hdf5";
$galacticus->{'store'} = 0;
$galacticus->{'tree' } = "all";
&Galacticus::HDF5::Get_Parameters($galacticus    );
&Galacticus::HDF5::Count_Trees   ($galacticus    );
&Galacticus::HDF5::Select_Output ($galacticus,0.0);
&Galacticus::HDF5::Get_Dataset($galacticus,['mergerTreeWeight','diskMassStellar','diskMassGas','hotHaloMass','hotHaloOutflowedMass','nodeIsIsolated','basicMassBertschinger','hotHaloUnaccretedMass','mergerTreeIndex']);
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
	     +$properties->{'hotHaloMass'          }->($satellites)->($inTree)->sum()
	     +$properties->{'hotHaloOutflowedMass' }->($satellites)->($inTree)->sum()
	     +$properties->{'hotHaloUnaccretedMass'}->($satellites)->($inTree)->sum()
	    )
	    /(
		+$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
		/$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
	    )
	    /$properties->{'basicMassBertschinger'}->($centrals)->(($i));
}
$properties->{'diskMassStellar'} /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'};
$properties->{'diskMassGas'} /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'  };
$properties->{'hotHaloMass'} /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'  };
$properties->{'hotHaloOutflowedMass'} /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
    )
    *$properties->{'basicMassBertschinger'  };
$properties->{'hotHaloUnaccretedMass'} /= 
    (
     +$parameters->{'cosmologyParameters'}->{'OmegaBaryon'}->{'value'}
     /$parameters->{'cosmologyParameters'}->{'OmegaMatter'}->{'value'}
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
    print "FAILED: mass conservation [simple] failure -> ".($massTotal-1.0)."\n";
} else {
    print "SUCCESS: mass conservation\n";
}
exit;
