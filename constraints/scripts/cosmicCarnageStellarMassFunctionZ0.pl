#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::IO::HDF5;
require Galacticus::Options;
require Galacticus::Constraints::MassFunctions;

# Compute likelihood (and make a plot) for a Galacticus model given the stellar mass function data for z=0 used in the Cosmic
# CARNage workshop.

# Data structure to hold the specification for our mass function.
my $massFunctionConfig;

# Get name of input and output files.
die("cosmicCarnageStellarMassFunctionZ0.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
$massFunctionConfig->{'self'          } = $0;
$massFunctionConfig->{'galacticusFile'} = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Specify the properties of this mass function.
my $entry                                             = 0;
$massFunctionConfig->{'redshift'                    } = pdl 0.0;
$massFunctionConfig->{'analysisLabel'               } = "carnageStellarMassFunctionZ0";
$massFunctionConfig->{'massType'                    } = "stellarMass";
$massFunctionConfig->{'massErrorRandomDex'          } = 0.08;
$massFunctionConfig->{'xRange'                      } = "1.0e7:1.0e12";
$massFunctionConfig->{'yRange'                      } = "3.0e-6:3.0e-1";
$massFunctionConfig->{'xLabel'                      } = "\$M_\\star\$ [\$M_\\odot\$]";
$massFunctionConfig->{'yLabel'                      } = "\${\\rm d}n/{\\rm d}\\log M_\\star\$ [Mpc\$^{-3}\$]";
$massFunctionConfig->{'title'                       } = "Stellar mass function at \$z=0\$";

# Read the observed data.
my $observations                                      = new PDL::IO::HDF5($galacticusPath."data/observations/cosmicCarnageWorkshop/massFunctionStellarZ0.hdf5");
$massFunctionConfig->{'x'                           }  = $observations->dataset('mass'                )->get    (                  );
$massFunctionConfig->{'y'                           }  = $observations->dataset('massFunctionObserved')->get    (                  );
$massFunctionConfig->{'yIsPer'                      }  = "ln";
$massFunctionConfig->{'xScaling'                    }  = "linear";
$massFunctionConfig->{'yScaling'                    }  = "linear";
$massFunctionConfig->{'covariance'                  }  = $observations->dataset('covariance'          )->get    (                  );
$massFunctionConfig->{'cosmologyScalingMass'        } = "none";
$massFunctionConfig->{'cosmologyScalingMassFunction'} = "none";
$massFunctionConfig->{'hubbleConstantObserved'      } = 67.770000;
$massFunctionConfig->{'omegaMatterObserved'         } =  0.307115;
$massFunctionConfig->{'omegaDarkEnergyObserved'     } =  0.692885;
$massFunctionConfig->{'observationLabel'            } = "Constraint";

# Limit comparison mass range.
$massFunctionConfig->{'constraintMassMinimum'       } = 2.0e8;

# Construct the mass function.
&MassFunctions::Construct(\%arguments,$massFunctionConfig);

exit;
