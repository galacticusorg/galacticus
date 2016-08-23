#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Galacticus::Path;
use PDL;
use PDL::IO::HDF5;
use Galacticus::Options;
use Galacticus::Constraints::FractionFunctions;

# Compute likelihood (and make a plot) for a Galacticus model given the passive fraction
# function data for z=2 used in the Cosmic CARNage workshop.

# Data structure to hold the specification for our mass function.
my $fractionFunctionConfig;

# Get name of input and output files.
die("cosmicCarnagePassiveFractionZ2.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
$fractionFunctionConfig->{'self'          } = $0;
$fractionFunctionConfig->{'galacticusFile'} = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Specify the properties of this fraction function.
my $entry                                                     = 0;
$fractionFunctionConfig->{'redshift'                        } = pdl 2.0;
$fractionFunctionConfig->{'analysisLabel'                   } = "carnagePassiveFractionFunctionZ2";
$fractionFunctionConfig->{'xRange'                          } = "1.0e9:1.0e12";
$fractionFunctionConfig->{'yRange'                          } = "1.0e-2:1.5";
$fractionFunctionConfig->{'xLabel'                          } = "\$M_\\star\$ [\$M_\\odot\$]";
$fractionFunctionConfig->{'yLabel'                          } = "\$f_{\\rm passive}\$ []";
$fractionFunctionConfig->{'title'                           } = "Passive fraction function at \$z=2\$";

# Read the observed data.
my $observations = new PDL::IO::HDF5(&galacticusPath()."data/observations/cosmicCarnageWorkshop/passiveFractionZ2.hdf5");
$fractionFunctionConfig->{'x'                               }  = $observations->dataset('mass'            )->get    (                  );
$fractionFunctionConfig->{'y'                               }  = $observations->dataset('fractionObserved')->get    (                  );
$fractionFunctionConfig->{'xScaling'                        }  = "linear";
$fractionFunctionConfig->{'yScaling'                        }  = "linear";
$fractionFunctionConfig->{'covariance'                      }  = $observations->dataset('covariance'      )->get    (                  );
$fractionFunctionConfig->{'cosmologyScalingMass'            } = "none";
$fractionFunctionConfig->{'cosmologyScalingFractionFunction'} = "none";
$fractionFunctionConfig->{'hubbleConstantObserved'          } = 67.770000;
$fractionFunctionConfig->{'omegaMatterObserved'             } =  0.307115;
$fractionFunctionConfig->{'omegaDarkEnergyObserved'         } =  0.692885;
$fractionFunctionConfig->{'observationLabel'                } = "Constraint";

# Construct the mass function.
&Galacticus::Constraints::FractionFunctions::Construct(\%arguments,$fractionFunctionConfig);

exit;
