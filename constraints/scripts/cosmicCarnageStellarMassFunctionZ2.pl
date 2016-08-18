#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Galacticus::Path;
use PDL;
use PDL::IO::HDF5;
use Galacticus::Options;
use Galacticus::Constraints::MassFunctions;

# Compute likelihood (and make a plot) for a Galacticus model given the stellar mass function data for z=2 used in the Cosmic
# CARNage workshop.

# Data structure to hold the specification for our mass function.
my $massFunctionConfig;

# Get name of input and output files.
die("cosmicCarnageStellarMassFunctionZ2.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
$massFunctionConfig->{'self'          } = $0;
$massFunctionConfig->{'galacticusFile'} = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Specify the properties of this mass function.
my $entry                                             = 0;
$massFunctionConfig->{'redshift'                    } = pdl 2.0;
$massFunctionConfig->{'analysisLabel'               } = "carnageStellarMassFunctionZ2";
$massFunctionConfig->{'massType'                    } = "stellarMass";
$massFunctionConfig->{'massErrorRandomDex'          } = 0.24;
$massFunctionConfig->{'xRange'                      } = "3.0e8:1.0e12";
$massFunctionConfig->{'yRange'                      } = "1.0e-6:3.0e-1";
$massFunctionConfig->{'xLabel'                      } = "\$M_\\star\$ [\$M_\\odot\$]";
$massFunctionConfig->{'yLabel'                      } = "\${\\rm d}n/{\\rm d}\\log M_\\star\$ [Mpc\$^{-3}\$]";
$massFunctionConfig->{'title'                       } = "Stellar mass function at \$z=2\$";

# Read the observed data.
my $observations                                      = new PDL::IO::HDF5(&galacticusPath()."data/observations/cosmicCarnageWorkshop/massFunctionStellarZ2.hdf5");
$massFunctionConfig->{'x'                           } = $observations->dataset('mass'                )->get    (                  );
$massFunctionConfig->{'y'                           } = $observations->dataset('massFunctionObserved')->get    (                  );
$massFunctionConfig->{'yIsPer'                      } = "ln";
$massFunctionConfig->{'xScaling'                    } = "linear";
$massFunctionConfig->{'yScaling'                    } = "linear";
$massFunctionConfig->{'covariance'                  } = $observations->dataset('covariance'          )->get    (                  );
$massFunctionConfig->{'cosmologyScalingMass'        } = "none";
$massFunctionConfig->{'cosmologyScalingMassFunction'} = "none";
$massFunctionConfig->{'hubbleConstantObserved'      } = 67.770000;
$massFunctionConfig->{'omegaMatterObserved'         } =  0.307115;
$massFunctionConfig->{'omegaDarkEnergyObserved'     } =  0.692885;
$massFunctionConfig->{'observationLabel'            } = "Constraint";

# Construct the mass function.
&Galacticus::Constraints::MassFunctions::Construct(\%arguments,$massFunctionConfig);

exit;
