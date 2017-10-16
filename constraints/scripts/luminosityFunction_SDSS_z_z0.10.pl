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

# Compute likelihood (and make a plot) for a Galacticus model given the z-band luminosity function from Montero-Dorta & Prada
# (2009; http://adsabs.harvard.edu/abs/2009MNRAS.399.1106M).

# Data structure to hold the specification for our mass function.
my $luminosityFunctionConfig;

# Get name of input and output files.
die("luminosityFunction_SDSS_z_z0.10.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
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
$luminosityFunctionConfig->{'magnitudes'         } = "yes";
$luminosityFunctionConfig->{'redshift'           } = pdl 0.01;
$luminosityFunctionConfig->{'analysisLabel'      } = "luminosityFunctionMonteroDorta2009SDSSz";
$luminosityFunctionConfig->{'discrepancyFileName'} = "discrepancy".ucfirst($luminosityFunctionConfig->{'analysisLabel'}).".hdf5";
$luminosityFunctionConfig->{'xRange'             } = "-24.3:-16.9";
$luminosityFunctionConfig->{'yRange'             } = "1.0e-8:1.0e-1";
$luminosityFunctionConfig->{'xLabel'             } = "\$M_\\mathrm{0.1z}\$";
$luminosityFunctionConfig->{'yLabel'             } = "\$\\mathrm{d}n/\\mathrm{d}M_\\mathrm{0.1z}\$ [Mpc\$^{-3}\$]";
$luminosityFunctionConfig->{'title'              } = "SDSS z-band luminosity function at \$z = 0.1\$";

# Read the observed data.
my $observations                                  = new PDL::IO::HDF5(&galacticusPath()."data/observations/luminosityFunctions/zLuminosityFunctionMonteroDorta2009SDSS.hdf5");
$luminosityFunctionConfig ->{'x'               }  = $observations->dataset('magnitudeAbsolute'      )->get();
$luminosityFunctionConfig ->{'y'               }  = $observations->dataset('luminosityFunction'     )->get();
$luminosityFunctionConfig ->{'error'           }  = $observations->dataset('luminosityFunctionError')->get();
($luminosityFunctionConfig->{'observationLabel'}) = "Montero-Dorta \& Prada (2009)";


# Construct the mass function.
&Galacticus::Constraints::LuminosityFunctions::Construct(\%arguments,$luminosityFunctionConfig);

exit;
