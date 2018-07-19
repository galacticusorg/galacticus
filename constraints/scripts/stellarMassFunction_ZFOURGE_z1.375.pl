#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use Galacticus::Options;
use Galacticus::Constraints::MassFunctions;

# Compute likelihood (and make a plot) for a Galacticus model given the stellar mass function data from Tmoczak (2014;
# http://adsabs.harvard.edu/abs/2014ApJ...783...85T).

# Data structure to hold the specification for our mass function.
my $massFunctionConfig;

# Get name of input and output files.
die("stellarMassFunction_ZFOURGE_z1.375.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
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
my $entry                                    = 0;
$massFunctionConfig->{'redshift'           } = pdl 1.375;
$massFunctionConfig->{'analysisLabel'      } = "zfourgeStellarMassFunctionZ1.375";
$massFunctionConfig->{'discrepancyFileName'} = "discrepancy".ucfirst($massFunctionConfig->{'analysisLabel'}).".hdf5";
$massFunctionConfig->{'massType'           } = "massStellar";
$massFunctionConfig->{'massErrorRandomDex' } = 0.2;
$massFunctionConfig->{'xRange'             } = "3.0e8:1.0e12";
$massFunctionConfig->{'yRange'             } = "3.0e-7:1.0e-1";
$massFunctionConfig->{'xLabel'             } = "\$M_\\star\$ [\$M_\\odot\$]";
$massFunctionConfig->{'yLabel'             } = "\$\\mathrm{d}n/\\mathrm{d}\\log M_\\star\$ [Mpc\$^{-3}\$]";
$massFunctionConfig->{'title'              } = "Stellar mass function at \$z\\approx 1.375\$";
$massFunctionConfig->{'systematicParameter'} = "zfourgeStellarMassFunctionZ1.375MassSystematic";
$massFunctionConfig->{'systematicOrder'    } = 1;
$massFunctionConfig->{'systematicZeroPoint'} = 11.3;

# Read the observed data.
my $observations                                      = new PDL::IO::HDF5(&galacticusPath()."data/observations/massFunctionsStellar/Stellar_Mass_Function_ZFOURGE_2014_z1.25_1.50.hdf5");
$massFunctionConfig ->{'x'                           }  = $observations->dataset('mass'                )->get    (                  );
$massFunctionConfig ->{'y'                           }  = $observations->dataset('massFunctionObserved')->get    (                  );
$massFunctionConfig ->{'yIsPer'                      }  = "ln";
$massFunctionConfig ->{'xScaling'                    }  = "linear";
$massFunctionConfig ->{'yScaling'                    }  = "linear";
$massFunctionConfig ->{'covariance'                  }  = $observations->dataset('covariance'          )->get    (                  );
$massFunctionConfig ->{'errorModel'                  }  = "logNormal"                                                                ;
($massFunctionConfig->{'observationLabel'            }) = $observations                                 ->attrGet('label'           );
($massFunctionConfig->{'hubbleConstantObserved'      }) = $observations->group  ('Parameters'          )->attrGet('H_0'             );
($massFunctionConfig->{'omegaMatterObserved'         }) = $observations->group  ('Parameters'          )->attrGet('Omega_Matter'    );
($massFunctionConfig->{'omegaDarkEnergyObserved'     }) = $observations->group  ('Parameters'          )->attrGet('Omega_DE'        );
($massFunctionConfig->{'cosmologyScalingMass'        }) = $observations->dataset('mass'                )->attrGet('cosmologyScaling');
($massFunctionConfig->{'cosmologyScalingMassFunction'}) = $observations->dataset('massFunction'        )->attrGet('cosmologyScaling');

# Construct the mass function.
&Galacticus::Constraints::MassFunctions::Construct(\%arguments,$massFunctionConfig);

exit;
