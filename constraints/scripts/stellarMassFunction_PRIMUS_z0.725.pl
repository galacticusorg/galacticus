#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use XML::Simple;
require Galacticus::Options;
require Galacticus::Constraints::MassFunctions;

# Compute likelihood (and make a plot) for a Galacticus model given the PRIMUS z=0.725 stellar mass function data from
# Moustakas et al. (2013; http://adsabs.harvard.edu/abs/2013ApJ...767...50M),

# Data structure to hold the specification for our mass function.
my $massFunctionConfig;

# Get name of input and output files.
die("stellarMassFunction_PRIMUS_z0.725.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
$massFunctionConfig->{'self'          } = $0;
$massFunctionConfig->{'galacticusFile'} = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments;
&Options::Parse_Options(\@ARGV,\%arguments);

# Specify the properties of this mass function.
my $entry                                    = 5;
$massFunctionConfig->{'redshift'           } = pdl 0.725;
$massFunctionConfig->{'analysisLabel'      } = "primusStellarMassFunctionZ0.725";
$massFunctionConfig->{'discrepancyFileName'} = "discrepancy".ucfirst($massFunctionConfig->{'analysisLabel'}).".hdf5";
$massFunctionConfig->{'massType'           } = "stellarMass";
$massFunctionConfig->{'massErrorRandomDex' } = 0.07;
$massFunctionConfig->{'xRange'             } = "5.0e8:2.0e12";
$massFunctionConfig->{'yRange'             } = "5.0e-8:1.0e-1";
$massFunctionConfig->{'xLabel'             } = "\$M_\\star\$ [\$M_\\odot\$]";
$massFunctionConfig->{'yLabel'             } = "\${\\rm d}n/{\\rm d}\\log M_\\star\$ [Mpc\$^{-3}\$]";
$massFunctionConfig->{'title'              } = "Stellar mass function at \$z\\approx 0.725\$";

# Read the observed data.
my $xml          = new XML::Simple;
my $observations = $xml->XMLin("data/observations/massFunctionsStellar/Stellar_Mass_Function_PRIMUS_2013.xml");
my $columns      = ${$observations->{'stellarMassFunction'}}[$entry]->{'columns'};
$massFunctionConfig->{'x'                         } = pdl @{$columns->{'stellarMass' }->{'datum'}};
$massFunctionConfig->{'y'                         } = pdl @{$columns->{'massFunction'}->{'datum'}};
$massFunctionConfig->{'yUpperError'               } = pdl @{$columns->{'upperError'  }->{'datum'}};
$massFunctionConfig->{'yLowerError'               } = pdl @{$columns->{'lowerError'  }->{'datum'}};
$massFunctionConfig->{'yIsPer'                    } = "log10";
$massFunctionConfig->{'xScaling'                  } = $columns->{'stellarMass' }->{'scaling'};
$massFunctionConfig->{'yScaling'                  } = $columns->{'massFunction'}->{'scaling'};
$massFunctionConfig->{'observationLabel'          } = $observations->{'label'};
$massFunctionConfig->{'hubbleConstantObserved'    } = $columns->{'stellarMass' }->{'hubble'        };
$massFunctionConfig->{'massHubbleExponent'        } = $columns->{'stellarMass' }->{'hubbleExponent'};
$massFunctionConfig->{'massFunctionHubbleExponent'} = $columns->{'massFunction'}->{'hubbleExponent'};

# Construct the mass function.
&MassFunctions::Construct(\%arguments,$massFunctionConfig);

exit;
