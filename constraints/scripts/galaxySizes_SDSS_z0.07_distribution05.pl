#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use Astro::Cosmology;
use XML::Simple;
use Data::Dumper;
use Galacticus::Options;
use Galacticus::HDF5;
use Galacticus::Constraints::Covariances;
use Galacticus::Constraints::DiscrepancyModels;
use Galacticus::Constraints::GalaxySizes;
use Stats::Histograms;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Compute likelihood (and make a plot) for a Galacticus model given the disk galaxy size data from Shen et al. (2003;
# http://adsabs.harvard.edu/abs/2003MNRAS.343..978S).
# Andrew Benson (18-July-2013)

# Get name of input and output files.
die("galaxySizes_SDSS_z0.07_distribution05.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named options.
my $iArg = -1;
my %options =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options            (\@ARGV                  ,\%options);
&Galacticus::Constraints::GalaxySizes::SDSS2003($galacticusFileName,"05",\%options);
exit 0;
