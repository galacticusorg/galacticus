#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use Galacticus::Options;
use Galacticus::Constraints::LuminosityFunctions;

# Compute likelihood (and make a plot) for a Galacticus model given the i-band luminosity function from Montero-Dorta & Prada
# (2009; http://adsabs.harvard.edu/abs/2009MNRAS.399.1106M).

# Data structure to hold the specification for our mass function.
my $luminosityFunctionConfig;

# Get name of input and output files.
die("luminosityFunction_SDSS_i_z0.10.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
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
$luminosityFunctionConfig->{'analysisLabel'      } = "luminosityFunctionMonteroDorta2009SDSSi";
$luminosityFunctionConfig->{'discrepancyFileName'} = "discrepancy".ucfirst($luminosityFunctionConfig->{'analysisLabel'}).".hdf5";
$luminosityFunctionConfig->{'xRange'             } = "-24.2:-16.5";
$luminosityFunctionConfig->{'yRange'             } = "1.0e-8:1.0e-1";
$luminosityFunctionConfig->{'xLabel'             } = "\$M_\\mathrm{0.1i}\$";
$luminosityFunctionConfig->{'yLabel'             } = "\$\\mathrm{d}n/\\mathrm{d}M_\\mathrm{0.1i}\$ [Mpc\$^{-3}\$]";
$luminosityFunctionConfig->{'title'              } = "SDSS i-band luminosity function at \$z = 0.1\$";

# Read the observed data.
my $observations                                  = new PDL::IO::HDF5(&galacticusPath()."data/observations/luminosityFunctions/iLuminosityFunctionMonteroDorta2009SDSS.hdf5");
$luminosityFunctionConfig ->{'x'               }  = $observations->dataset('magnitudeAbsolute'      )->get();
$luminosityFunctionConfig ->{'y'               }  = $observations->dataset('luminosityFunction'     )->get();
$luminosityFunctionConfig ->{'error'           }  = $observations->dataset('luminosityFunctionError')->get();
($luminosityFunctionConfig->{'observationLabel'}) = "Montero-Dorta \& Prada (2009)";


# Construct the mass function.
&Galacticus::Constraints::LuminosityFunctions::Construct(\%arguments,$luminosityFunctionConfig);

exit;
