#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
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
use Galacticus::Constraints::Concentrations;
use Stats::Histograms;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Compute likelihood (and make a plot) for a Galacticus model given the concentratio distribution in a mass bin from the COCO CDM
# simulation (Wojciech et al.; 2016; MNRAS; 457; 3492; http://adsabs.harvard.edu/abs/2016MNRAS.457.3492H).
# Andrew Benson (20-April-2018)

# Get name of input and output files.
die("concentrations_CDM_COCO_distribution5.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named options.
my $iArg = -1;
my %options =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options              (\@ARGV                  ,\%options);
&Galacticus::Constraints::Concentrations::COCOCDM($galacticusFileName,"05",\%options);
exit 0;
