#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Galacticus::Path;
use PDL;
use PDL::IO::HDF5;
use Galacticus::Options;
use Galacticus::Constraints::LuminosityFunctions;

# Compute likelihood (and make a plot) for a Galacticus model given the z=1.47 Halpha luminosity function from Sobral et al. (2013;
# http://adsabs.harvard.edu/abs/2013MNRAS.428.1128S).

# Data structure to hold the specification for our mass function.
my $luminosityFunctionConfig;

# Get name of input and output files.
die("luminosityFunction_HiZELS_Halpha_z1.47.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
$luminosityFunctionConfig->{'self'          } = $0;
$luminosityFunctionConfig->{'galacticusFile'} = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Specify the properties of this mass function.
$luminosityFunctionConfig->{'magnitudes'         } = "no";
$luminosityFunctionConfig->{'redshift'           } = pdl 1.47;
$luminosityFunctionConfig->{'analysisLabel'      } = "luminosityFunctionHalphaSobral2013HiZELSZ3";
$luminosityFunctionConfig->{'discrepancyFileName'} = "discrepancy".ucfirst($luminosityFunctionConfig->{'analysisLabel'}).".hdf5";
$luminosityFunctionConfig->{'xRange'             } = "2.0e40:5.0e43";
$luminosityFunctionConfig->{'yRange'             } = "2.0e-6:1.0e-1";
$luminosityFunctionConfig->{'xLabel'             } = "\$L_{\\mathrm{H}\\alpha}\$ [ergs s\$^{-1}\$]";
$luminosityFunctionConfig->{'yLabel'             } = "\$\\mathrm{d}n/\\mathrm{d}\\log L_{\\mathrm{H}\\alpha}\$ [Mpc\$^{-3}\$]";
$luminosityFunctionConfig->{'title'              } = "HiZELS H\$\\alpha\$ luminosity function at \$z \\approx 1.47\$";

# Read the observed data.
my $observations                                  = new PDL::IO::HDF5(&galacticusPath()."data/observations/luminosityFunctions/hAlphaLuminosityFunctionSobral2013HiZELSZ1.47.hdf5");
$luminosityFunctionConfig ->{'x'               }  = $observations->dataset('luminosity'             )->get();
$luminosityFunctionConfig ->{'y'               }  = $observations->dataset('luminosityFunction'     )->get();
$luminosityFunctionConfig ->{'error'           }  = $observations->dataset('luminosityFunctionError')->get();
($luminosityFunctionConfig->{'observationLabel'}) = "Sobral et al. (2013)";

# Construct the mass function.
&Galacticus::Constraints::LuminosityFunctions::Construct(\%arguments,$luminosityFunctionConfig);

exit;
