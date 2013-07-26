#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
require Galacticus::Options;
require Galacticus::Constraints::MassFunctions;

# Compute likelihood (and make a plot) for a Galacticus model given the HI mass function data from Martin et al. (2010;
# http://adsabs.harvard.edu/abs/2010ApJ...723.1359M),

# Data structure to hold the specification for our mass function.
my $massFunctionConfig;

# Get name of input and output files.
die("hiMassFunction_ALFALFA_z0.0.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
$massFunctionConfig->{'self'          } = $0;
$massFunctionConfig->{'galacticusFile'} = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments;
&Options::Parse_Options(\@ARGV,\%arguments);

# Specify the properties of this mass function.
my $entry                                    = 0;
$massFunctionConfig->{'redshift'           } = pdl 0.000;
$massFunctionConfig->{'analysisLabel'      } = "alfalfaHiMassFunctionZ0.00";
$massFunctionConfig->{'discrepancyFileName'} = "discrepancy".ucfirst($massFunctionConfig->{'analysisLabel'}).".hdf5";
$massFunctionConfig->{'massType'           } = "hiGasMass";
$massFunctionConfig->{'massErrorRandomDex' } = 0.08;
$massFunctionConfig->{'xRange'             } = "1.0e6:3.0e11";
$massFunctionConfig->{'yRange'             } = "1.0e-6:1.0e0";
$massFunctionConfig->{'xLabel'             } = "\$M_{\\rm HI}\$ [\$M_\\odot\$]";
$massFunctionConfig->{'yLabel'             } = "\${\\rm d}n/{\\rm d}\\log M_{\\rm HI}\$ [Mpc\$^{-3}\$]";
$massFunctionConfig->{'title'              } = "HI mass function at \$z\\approx 0.00\$";

# Read the observed data.
my $xml                          = new XML::Simple;
my $observations                 = $xml->XMLin("data/observations/massFunctionsHI/HI_Mass_Function_ALFALFA_2010.xml");
my $columns                      = $observations->{'massFunctions'}->{'columns'};
$massFunctionConfig->{'x'                         } = pdl @{$observations->{'massFunctions'}->{'columns'}->{'mass'        }->{'data'}};
$massFunctionConfig->{'y'                         } = pdl @{$observations->{'massFunctions'}->{'columns'}->{'massFunction'}->{'data'}};
$massFunctionConfig->{'yUpperError'               } = pdl @{$observations->{'massFunctions'}->{'columns'}->{'upperError'  }->{'data'}};
$massFunctionConfig->{'yLowerError'               } = pdl @{$observations->{'massFunctions'}->{'columns'}->{'lowerError'  }->{'data'}};
$massFunctionConfig->{'yIsPer'                    } = "log10";
$massFunctionConfig->{'xScaling'                  } = $columns->{'mass'        }->{'scaling'};
$massFunctionConfig->{'yScaling'                  } = $columns->{'massFunction'}->{'scaling'};
$massFunctionConfig->{'observationLabel'          } = $observations->{'massFunctions'}->{'label'};
$massFunctionConfig->{'hubbleConstantObserved'    } = $observations->{'massFunctions'}->{'columns'}->{'mass'        }->{'hubble'        };
$massFunctionConfig->{'massHubbleExponent'        } = $observations->{'massFunctions'}->{'columns'}->{'mass'        }->{'hubbleExponent'};
$massFunctionConfig->{'massFunctionHubbleExponent'} = $observations->{'massFunctions'}->{'columns'}->{'massFunction'}->{'hubbleExponent'};

# Construct the mass function.
&MassFunctions::Construct(\%arguments,$massFunctionConfig);

exit;
