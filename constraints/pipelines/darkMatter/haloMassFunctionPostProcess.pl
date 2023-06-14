#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use PDL::Math;
use XML::Simple;
use Data::Dumper;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;

# Generate a halo mass function using the optimal parameters.
# Andrew Benson (22-September-2020)

# Get arguments.
die('Usage: haloMassFunctionPostProcess.pl <outputDirectory> [options...]')
    unless ( scalar(@ARGV) >= 1 );
my $outputDirectory = $ARGV[0];
my %options;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Define simulations to process.
my @simulations =
(
 {
     label               => "MilkyWay",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo023", "Halo088", "Halo119", "Halo188", "Halo247", "Halo268", "Halo270", "Halo288", "Halo327", "Halo349", "Halo364", "Halo374", "Halo414", "Halo415", "Halo416", "Halo440", "Halo460", "Halo469", "Halo490", "Halo530", "Halo558", "Halo567", "Halo570", "Halo606", "Halo628", "Halo641", "Halo675", "Halo718", "Halo738", "Halo749", "Halo797", "Halo800", "Halo825", "Halo829", "Halo852", "Halo878", "Halo881", "Halo925", "Halo926", "Halo937", "Halo939", "Halo967", "Halo9749", "Halo9829", "Halo990", "Halo004", "Halo113" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0bf7f4",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0bf7f4",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM3",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 3keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM6.5",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 6.5keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM10",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 10keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM3_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 3keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0bf7f4",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM4_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 4keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0bf7f4",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM5_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 5keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0bf7f4",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM6_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 6keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0bf7f4",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM6.5_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 6.5keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0bf7f4",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_WDM10_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way WDM 10keV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0bf7f4",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_25.9e-22eV",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 25.9e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_69.4e-22eV",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 69.4e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_113e-22eV",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 113e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_151e-22eV",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 151e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_185e-22eV",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 185e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_490e-22eV",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 185e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_25.9e-22eV_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 25.9e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_69.4e-22eV_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 69.4e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_113e-22eV_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 113e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_151e-22eV_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 151e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_185e-22eV_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 185e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_fdm_490e-22eV_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way FDM 185e-22 eV zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#0ba0f7",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1GeV_envelope",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1GeV (envelope) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1GeV_halfmode",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1GeV (halfmode) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1e-2GeV_envelope",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1e-2GeV (envelope) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1e-2GeV_halfmode",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1e-2GeV (halfmode) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1e-4GeV_envelope",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1e-4GeV (envelope) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1e-4GeV_halfmode",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1e-4GeV (halfmode) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004", "Halo113", "Halo023" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 2.81981e5,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1GeV_envelope_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1GeV (envelope) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1GeV_halfmode_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1GeV (halfmode) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1e-2GeV_envelope_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1e-2GeV (envelope) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1e-2GeV_halfmode_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1e-2GeV (halfmode) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1e-4GeV_envelope_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1e-4GeV (envelope) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "MilkyWay_IDM1e-4GeV_halfmode_hires",
     description         => "Halo mass function for non-backsplash z=0 halos from Milky Way, 1e-4GeV (halfmode) IDM zoom-in simulations.",
     subpath             => "ZoomIns",
     realizations        => [ "Halo004" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000 ],
     color               => "#f70bdf",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "LMC",
     description         => "Halo mass function for non-backsplash z=0 halos from LMC zoom-in simulations.",
     subpath             => "ZoomIns",
     # Halo479 and Halo888 have been excluded due to contamination by low-res particles 
     realizations        => [ "Halo032", "Halo059", "Halo0662", "Halo083", "Halo088", "Halo097", "Halo104", "Halo110", "Halo202", "Halo208", "Halo218", "Halo296", "Halo301", "Halo303", "Halo340", "Halo374", "Halo380", "Halo391", "Halo405", "Halo440", "Halo463", "Halo4662", "Halo511", "Halo524", "Halo539", "Halo567", "Halo575", "Halo602", "Halo697", "Halo711", "Halo721", "Halo767", "Halo802", "Halo824", "Halo850", "Halo853", "Halo914", "Halo932", "Halo933" ],
     simulationReference => "Nadler et al.",
     simulationURL       => "https://www",
     hubbleConstant      => 0.7,
     massParticle        => 3.52476e4,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     color               => "#e6e6fa",
     plotModify          => \&zoomInsPlotModify,
     validate            => \&zoomInsValidate
 },
 {
     label               => "VSMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the VSMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/vsmdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 6.2e6,
     boxLength           => 160.0,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000, 0.67110, 0.50250, 0.33000, 0.24750 ],
     color               => "#46019b"
 },
 {
     label               => "SMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the SMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/smdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 9.63e7,
     boxLength           => 400.0,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000, 0.66430, 0.50000, 0.33100, 0.24800 ],
     color               => "#007efe"
 },
 {
     label               => "MDPL2",
     description         => "Halo mass function for non-backsplash z=0 halos from the MDPL2 simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/mdpl2/",
     hubbleConstant      => 0.6777,
     massParticle        => 1.51e9,
     boxLength           => 1000.0,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000, 0.67120, 0.50320, 0.33030, 0.24230 ],
     color               => "#00bb00"
 },
 {
     label               => "BigMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the BigMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/bigmdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 2.359e10,
     boxLength           => 2500.0,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000, 0.67040, 0.50000, 0.31800 ],
     color               => "#fef601"
 },
 {
     label               => "HugeMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the HugeMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/hugemdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 7.9e10,
     boxLength           => 4000.0,
     countHaloMinimum    => 0,
     expansionFactors    => [ 1.00000, 0.67120, 0.50320, 0.33030 ],
     color               => "#dd0000"
 }
);

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Get an XML object.
my $xml = new XML::Simple();

# Specify the minimum number of particles used in fitting halo mass functions.
my $countParticlesMinimum = 100;

# Iterate over simulations to get the model mass function.
my @jobs;
foreach my $simulation ( @simulations ) {

    # Compute darker color.
    $simulation->{'colorDark'} =
	sprintf(
	    "#%02lx%02lx%02lx",
	    hex(substr($simulation->{'color'},1,2))/2,
	    hex(substr($simulation->{'color'},3,2))/2,
	    hex(substr($simulation->{'color'},5,2))/2,
	);

    # Add a single, null realization if the simulation has none.
    $simulation->{'realizations'} = [ "" ]
	unless ( exists($simulation->{'realizations'}) );

    # Convert particle mass to Solar masses.
    $simulation->{'massParticle'} /= $simulation->{'hubbleConstant'};

    # Convert box length to Mpc.
    $simulation->{'boxLength'} /= $simulation->{'hubbleConstant'}
        if ( exists($simulation->{'boxLength'}) );

    # Iterate over realizations.
    foreach my $realization ( @{$simulation->{'realizations'}} ) {
	my $realizationLabel = $realization eq "" ? "" : "_".$realization;
	
	# Iterate over expansion factors.
	foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	    my $redshift      = 1.0/$expansionFactor-1.0;
	    my $redshiftLabel = sprintf("z%5.3f",$redshift);

	    # Validate parameters.
	    my $parameters = $xml->XMLin($outputDirectory."/haloMassFunctionBase_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".xml");
	    &{$simulation->{'validate'}}($simulation,$parameters,$realization,$expansionFactor)
		if ( exists($simulation->{'validate'}) );

	    # Generate a job to create our mass function fit.
	    {
		my $job;
		$job->{'command'   } =
		    "Galacticus.exe ".$outputDirectory."/haloMassFunctionBase_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".xml";
		$job->{'launchFile'} = $outputDirectory."/haloMassFunction_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".sh" ;
		$job->{'logFile'   } = $outputDirectory."/haloMassFunction_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".log";
		$job->{'label'     } =                   "haloMassFunction_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel       ;
		$job->{'ppn'       } = 1;
		$job->{'nodes'     } = 1;
		$job->{'mpi'       } = "no";
		push(@jobs,$job)
		    unless ( -e $outputDirectory."/haloMassFunction".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.":MPI0000.hdf5" );
	    }

	    # Modify parameters to use the Despali et al. (2015) mass function.
	    $parameters->{'outputFileName'}->{'value'} = $outputDirectory."/haloMassFunction_Despali2015".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".hdf5";
	    $parameters->{'haloMassFunction'}->{'value'} = "despali2015";
	    delete($parameters->{'haloMassFunction'}->{$_})
	    	foreach ( "errorFractionalMaximum", "toleranceRelative", "haloMassFunction" );
	    open(my $parameterFile,">",$outputDirectory."/haloMassFunctionBase_Despali2015_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".xml");
	    print $parameterFile $xml->XMLout($parameters,RootName => "parameters");
	    close($parameterFile);
	    # Generate a job to create Despali et al. (2015) mass functions.
	    {
	    	my $job;
	    	$job->{'command'   } =
	    	    "Galacticus.exe ".$outputDirectory."/haloMassFunctionBase_Despali2015_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".xml";
	    	$job->{'launchFile'} = $outputDirectory."/haloMassFunction_Despali2015_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".sh" ;
	    	$job->{'logFile'   } = $outputDirectory."/haloMassFunction_Despali2015_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".log";
	    	$job->{'label'     } =                   "haloMassFunction_Despali2015_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel       ;
	    	$job->{'ppn'       } = 1;
	    	$job->{'nodes'     } = 1;
	    	$job->{'mpi'       } = "no";
	    	push(@jobs,$job)
	    	    unless ( -e $outputDirectory."/haloMassFunction_Despali2015".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.":MPI0000.hdf5" );
	    }
	}
    }
}
&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs)
    if ( scalar(@jobs) > 0 );

# Read all data.
my $iMax = -1;
my @expansionFactors;
foreach my $simulation ( @simulations ) {
    # Iterate over realizations.
    foreach my $realization ( @{$simulation->{'realizations'}} ) {
	my $realizationLabel = $realization eq "" ? "" : "_".$realization;
	# Iterate over expansion factors.
	my $i = -1;
	foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	    ++$i;
	    $expansionFactors[$i] =     $expansionFactor    ;
	    my $redshift          = 1.0/$expansionFactor-1.0;
	    my $redshiftLabel     = sprintf("z%5.3f",$redshift);

	    # Read the resulting halo mass function.
	    my $model                                                                                = new PDL::IO::HDF5($outputDirectory."/haloMassFunction".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.":MPI0000.hdf5");
	    my $outputs                                                                              = $model           ->group  ('Outputs'                       )       ;
	    my $output                                                                               = $outputs         ->group  ('Output1'                       )       ;
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'mass'              } = $output          ->dataset('haloMass'                      )->get();
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunction'      } = $output          ->dataset('haloMassFunctionLnM'           )->get();
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionBinned'} = $output          ->dataset('haloMassFunctionLnMBinAveraged')->get();
	    
	    # Read the target dataset.
	    my $target                                                                               = new PDL::IO::HDF5($ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".hdf5");
	    my $targetSimulation                                                                     = $target          ->group  ('simulation0001'                )       ;
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'             } = $targetSimulation->dataset('mass'                          )->get();
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'     } = $targetSimulation->dataset('massFunction'                  )->get();
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'count'            } = $targetSimulation->dataset('count'                         )->get();	    

	    # Construct the target dataset errors.
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}                    =       $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}                            ->copy()     ;
	    my $nonZeroTarget                                                                                           = which($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'count'       }                                     > 0);
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}->($nonZeroTarget) /=       $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'count'       }->double()->($nonZeroTarget)->sqrt()     ;

	    # Interpolate mass function.
	    ($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'})  
		= interpolate(
		$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'              }->log(),
		$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'mass'              }->log(),
		$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionBinned'}->log()
		);
	    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'} .= exp($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'});
	}
	$iMax = $i
	    if ( $i > $iMax );
    }
}

# Begin creating the plots.
for(my $i=0;$i<=$iMax;++$i) {
    my $expansionFactor = $expansionFactors[$i];
    my $redshift      = 1.0/$expansionFactor-1.0;
    my $redshiftLabel = sprintf("z%5.3f",$redshift);
    unless ( -e $outputDirectory."/haloMassFunction_".$redshiftLabel.".pdf" ) {
	## Halo mass function.
	my $plot;
	my $gnuPlot;
	my $plotFileTeX = $outputDirectory."/haloMassFunction_".$redshiftLabel.".tex";
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
	print $gnuPlot "set output '".$plotFileTeX."'\n";
	print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
	print $gnuPlot "set ylabel '\$ \\mathrm{d} n / \\mathrm{d} \\log M \$ [Mpc\$^{-3}\$]\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.475,0.45\n";
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [8.0e6:1.0e16]\n";
	print $gnuPlot "set yrange [1.0e-9:1.0e2]\n";
	print $gnuPlot "set pointsize 1.0\n";
	foreach my $simulation ( @simulations ) {
	    next
		unless ( scalar(@{$simulation->{'expansionFactors'}}) > $i );
	    my $massHaloMinimum = $simulation->{'massParticle'}*$countParticlesMinimum;
	    my $xLimit          = pdl [ $massHaloMinimum, $massHaloMinimum ];
	    my $yLimit          = pdl [ 1.0e-9          , 1.0e2            ];
	    (my $title          = $simulation->{'label'}) =~ s/_/ /g;
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                             ,
		$xLimit                                                            ,
		$yLimit                                                            ,
		style        => "line"                                             ,
		weight       => [1,1]                                              ,
		linePattern  => 3                                                  ,
		color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		);
 	    # Iterate over realizations.
	    my $massTarget;
	    my $massFunctionTarget;
	    my $massFunctionTargetError;
	    my $massModel;
	    my $massFunctionModel;
	    my $massFunctionModelBinned;
	    foreach my $realization ( @{$simulation->{'realizations'}} ) {
		my $realizationLabel = $realization eq "" ? "" : "_".$realization;
		my $countHalos =
		    (
		     +$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'     }
		     /$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}
		    )**2;
		my $selectConversion = which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'        } > $massHaloMinimum)
			&
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'} > 0.0             )
		    );
		my $conversion = average($countHalos->($selectConversion)/$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}->($selectConversion));
		my $countHalosModel = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'}*$conversion;
		(my $selected, my $unselected) =
		    which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'            } >  $massHaloMinimum)
			&
			($simulation                                                      ->{'countHaloMinimum'} <= $countHalosModel)
		    );
		my $selectedModel =
		    which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'mass'            } >  $massHaloMinimum)
		    );
		unless ( defined($massTarget) ) {
		    $massTarget               = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'              }->($selected     )->copy()   ;
		    $massFunctionTarget       = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'      }->($selected     )->copy()   ;
		    $massFunctionTargetError  = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError' }->($selected     )->copy()**2;
		    $massModel                = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'mass'              }->($selectedModel)->copy()   ;
		    $massFunctionModel        = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunction'      }->($selectedModel)->copy()   ;
		    $massFunctionModelBinned  = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionBinned'}->($selectedModel)->copy()   ;
		} else {
		    $massFunctionTarget      += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'      }->($selected     )           ;
		    $massFunctionTargetError += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError' }->($selected     )        **2;
		    $massFunctionModel       += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunction'      }->($selectedModel)           ;
		    $massFunctionModelBinned += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionBinned'}->($selectedModel)           ;
		}
	    }
	    $massFunctionTargetError .= $massFunctionTargetError->sqrt();
	    $massFunctionTarget      /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionTargetError /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionModel       /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionModelBinned /= scalar(@{$simulation->{'realizations'}});
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                             ,
		$massTarget                                                        ,
		$massFunctionTarget                                                ,
		errorUp      => $massFunctionTargetError                           ,
		errorDown    => $massFunctionTargetError                           ,
		style        => "point"                                            ,
		weight       => [2,1]                                              ,
		symbol       => [6,6]                                              ,
		pointSize    => 0.5                                                ,
		color        => [$simulation->{'colorDark'},$simulation->{'color'}],
		%options
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                                  ,
		$massModel,
		$massFunctionModel,
		style        => "line"                                                                                  ,
		weight       => [3,1]                                                                                   ,
		color        => [$simulation->{'colorDark'},$simulation->{'color'}],
		title        => $title
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                                    ,
		$massModel,
		$massFunctionModelBinned,
		style        => "point"                                                                                   ,
		weight       => [2,1]                                                                                     ,
		symbol       => [2,2]                                                                                     ,
		pointSize    => 0.3                                                                                       ,
		color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		);
	}  
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
    }

    # Individual plots.
    foreach my $simulation ( @simulations ) {
	next
	    unless ( scalar(@{$simulation->{'expansionFactors'}}) > $i );
	# Iterate over realizations.
	foreach my $realization ( @{$simulation->{'realizations'}} ) {
	    my $realizationLabel = $realization eq "" ? "" : "_".$realization;
	    unless ( -e $outputDirectory."/haloMassFunction_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".pdf" ) {
		# If a custom plot modifier function is defined, call it.
		my $plotOptions;
		($plotOptions->{'title'   } = $simulation->{'label'}) =~ s/_/ /g;
		$plotOptions ->{'xMinimum'} = 8.0e+06;
		$plotOptions ->{'xMaximum'} = 1.0e+16;
		$plotOptions ->{'yMinimum'} = 1.0e-09;
		$plotOptions ->{'yMaximum'} = 1.0e+02;
		$plotOptions ->{'xKey'    } = 0.475;
		$plotOptions ->{'yKey'    } = 0.450;
		&{$simulation->{'plotModify'}}($simulation,$plotOptions,$realization,$simulation->{'expansionFactors'}->[$i])
		    if ( exists($simulation->{'plotModify'}) );
		## Halo mass function.
		my $plot;
		my $gnuPlot;
		my $plotFileTeX = $outputDirectory."/haloMassFunction_".$simulation->{'label'}.$realizationLabel."_".$redshiftLabel.".tex";
		open($gnuPlot,"|gnuplot");
		print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
		print $gnuPlot "set output '".$plotFileTeX."'\n";
		print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
		print $gnuPlot "set ylabel '\$ \\mathrm{d} n / \\mathrm{d} \\log M \$ [Mpc\$^{-3}\$]\n";
		print $gnuPlot "set lmargin screen 0.15\n";
		print $gnuPlot "set rmargin screen 0.95\n";
		print $gnuPlot "set bmargin screen 0.15\n";
		print $gnuPlot "set tmargin screen 0.95\n";
		print $gnuPlot "set key spacing 1.2\n";
		print $gnuPlot "set key at screen ".sprintf("%5.3f",$plotOptions->{'xKey'}).",".sprintf("%5.3f",$plotOptions->{'yKey'})."\n";
		print $gnuPlot "set logscale xy\n";
		print $gnuPlot "set mxtics 10\n";
		print $gnuPlot "set mytics 10\n";
		print $gnuPlot "set format x '\$10^{\%L}\$'\n";
		print $gnuPlot "set format y '\$10^{\%L}\$'\n";
		print $gnuPlot "set xrange [".sprintf("%12.6e",$plotOptions->{'xMinimum'}).":".sprintf("%12.6e",$plotOptions->{'xMaximum'})."]\n";
		print $gnuPlot "set yrange [".sprintf("%12.6e",$plotOptions->{'yMinimum'}).":".sprintf("%12.6e",$plotOptions->{'yMaximum'})."]\n";
		print $gnuPlot "set pointsize 1.0\n";
		my $massHaloMinimum = $simulation->{'massParticle'}*$countParticlesMinimum;
		my $xLimit          = pdl [ $massHaloMinimum, $massHaloMinimum ];
		my $yLimit          = pdl [ 1.0e-9          , 1.0e2            ];
		&GnuPlot::PrettyPlots::Prepare_Dataset(
		    \$plot                                                                                                  ,
		    $xLimit                                                                                                 ,
		    $yLimit                                                                                                 ,
		    style        => "line"                                                                                  ,
		    weight       => [1,1]                                                                                   ,
		    linePattern  => 3                                                                                       ,
		    color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		    );
		my $countHalos =
		    (
		     +$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'     }
		     /$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}
		    )**2;
		my $selectConversion = which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'        } > $massHaloMinimum)
			&
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'} > 0.0             )
		    );
		my $conversion = average($countHalos->($selectConversion)/$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}->($selectConversion));
		my $countHalosModel = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'}*$conversion;
		(my $selected, my $unselected) =
		    which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'            } >  $massHaloMinimum)
			&
			($simulation                                                      ->{'countHaloMinimum'} <= $countHalosModel)
		    );
		my $selectedModel =
		    which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'mass'            } >  $massHaloMinimum)
		    );
		&GnuPlot::PrettyPlots::Prepare_Dataset(
		    \$plot                                                                                                               ,
		    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'        }->($selected)                     ,
		    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}->($selected)                     ,
		    errorUp      => $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}->($selected),
		    errorDown    => $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}->($selected),
		    style        => "point"                                                                                              ,
		    weight       => [2,1]                                                                                                ,
		    symbol       => [6,6]                                                                                                ,
		    pointSize    => 0.5                                                                                                  ,
		    color        => [$simulation->{'colorDark'},$simulation->{'color'}]                                                  ,
		    title        => $plotOptions->{'title'}
		    );
		&GnuPlot::PrettyPlots::Prepare_Dataset(
		    \$plot                                                                                                  ,
		    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'mass'        }->($selectedModel)    ,
		    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunction'}->($selectedModel)    ,
		    style        => "line"                                                                                  ,
		    weight       => [3,1]                                                                                   ,
		    color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		    );
		&GnuPlot::PrettyPlots::Prepare_Dataset(
		    \$plot                                                                                                    ,
		    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'mass'              }->($selectedModel),
		    $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionBinned'}->($selectedModel),
		    style        => "point"                                                                                   ,
		    weight       => [2,1]                                                                                     ,
		    symbol       => [2,2]                                                                                     ,
		    pointSize    => 0.3                                                                                       ,
		    color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		    );		
		&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
		close($gnuPlot);
		&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
	    }
	}
    }

    # Realization plots.
    foreach my $simulation ( @simulations ) {
	next
	    unless ( scalar(@{$simulation->{'expansionFactors'}}) > $i );
	next
	    unless ( scalar(@{$simulation->{'realizations'    }}) > 1  );
	# Skip if plot exists.
	next
	    if ( -e $outputDirectory."/haloMassFunctionRealizations_".$simulation->{'label'}."_".$redshiftLabel.".pdf" );
	# If a custom plot modifier function is defined, call it.
	my $plotOptions;
	($plotOptions->{'title'   } = $simulation->{'label'}) =~ s/_/ /g;
	$plotOptions ->{'xMinimum'} = 8.0e+06;
	$plotOptions ->{'xMaximum'} = 1.0e+16;
	$plotOptions ->{'yMinimum'} = 1.0e-09;
	$plotOptions ->{'yMaximum'} = 1.0e+02;
	$plotOptions ->{'xKey'    } = 0.475;
	$plotOptions ->{'yKey'    } = 0.450;
	&{$simulation->{'plotModify'}}($simulation,$plotOptions,undef(),$simulation->{'expansionFactors'}->[$i])
	    if ( exists($simulation->{'plotModify'}) );
	## Halo mass function.
	my $plot;
	my $gnuPlot;
	my $plotFileTeX = $outputDirectory."/haloMassFunctionRealizations_".$simulation->{'label'}."_".$redshiftLabel.".tex";
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
	print $gnuPlot "set output '".$plotFileTeX."'\n";
	print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
	print $gnuPlot "set ylabel '\$ \\mathrm{d} n / \\mathrm{d} \\log M \$ [Mpc\$^{-3}\$]\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen ".sprintf("%5.3f",$plotOptions->{'xKey'}).",".sprintf("%5.3f",$plotOptions->{'yKey'})."\n";
	print $gnuPlot "set logscale x\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [".sprintf("%12.6e",$plotOptions->{'xMinimum'}).":".sprintf("%12.6e",$plotOptions->{'xMaximum'})."]\n";
	print $gnuPlot "set yrange [0.0:2.0]\n";
	print $gnuPlot "set pointsize 1.0\n";
	my $massHaloMinimum = $simulation->{'massParticle'}*$countParticlesMinimum;
	my $xLimit          = pdl [ $massHaloMinimum, $massHaloMinimum ];
	my $yLimit          = pdl [ 0.0             , 2.0              ];
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	    \$plot                                                                                                  ,
	    $xLimit                                                                                                 ,
	    $yLimit                                                                                                 ,
	    style        => "line"                                                                                  ,
	    weight       => [1,1]                                                                                   ,
	    linePattern  => 3                                                                                       ,
	    color        => [$simulation->{'colorDark'},$simulation->{'color'}]
	    );
	# Iterate over realizations.
	my $massFunctionTarget;
	my $massFunctionTargetError;
	my $massFunctionModel;
	foreach my $realization ( @{$simulation->{'realizations'}} ) {
	    my $realizationLabel = $realization eq "" ? "" : "_".$realization;
	    my $countHalos =
		(
		 +$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'     }
		 /$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}
		)**2;
	    my $selectConversion = which(
		($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'        } > $massHaloMinimum)
		&
		($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'} > 0.0             )
		);
	    my $conversion = average($countHalos->($selectConversion)/$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}->($selectConversion));
	    my $countHalosModel = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'}*$conversion;
	    (my $selected, my $unselected) =
		which(
		    ($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'            } >  $massHaloMinimum)
		    &
		    ($simulation                                                      ->{'countHaloMinimum'} <= $countHalosModel)
		);
	    my $selectedModel =
		which(
		    ($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'mass'            } >  $massHaloMinimum)
		);
	    unless ( defined($massFunctionModel) ) {
		$massFunctionTarget       = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'      }->($selected     )->copy()   ;
		$massFunctionTargetError  = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError' }->($selected     )->copy()**2;
		$massFunctionModel        = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunction'      }->($selectedModel)->copy()   ;
	    } else {
		$massFunctionTarget      += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'      }->($selected     )           ;
		$massFunctionTargetError += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError' }->($selected     )        **2;
		$massFunctionModel       += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunction'      }->($selectedModel)           ;
	    }
	}
	$massFunctionTargetError .= $massFunctionTargetError->sqrt();
	$massFunctionTarget      /= scalar(@{$simulation->{'realizations'}});
	$massFunctionTargetError /= scalar(@{$simulation->{'realizations'}});
	$massFunctionModel       /= scalar(@{$simulation->{'realizations'}});
	# Iterate over realizations.
	foreach my $realization ( @{$simulation->{'realizations'}} ) {
	    my $redshiftLabelLocal   = sprintf("%5.3f",1.0/$simulation->{'expansionFactors'}->[$i]-1.0);
	    my $target                   = new PDL::IO::HDF5($ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$simulation->{'label'}."_".$realization."_z".$redshiftLabelLocal.".hdf5");
	    my $targetSimulation         = $target          ->group  ('simulation0001'        );
	    (my $massEnvironment       ) = $targetSimulation->attrGet('massEnvironment'       );
	    (my $overdensityEnvironment) = $targetSimulation->attrGet('overdensityEnvironment');
	    my $realizationLabel = $realization eq "" ? "" : "_".$realization;
	    my $countHalos =
		(
		 +$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'     }
		 /$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}
		)**2;
	    my $selectConversion = which(
		($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'        } > $massHaloMinimum)
		&
		($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'} > 0.0             )
		);
	    my $conversion = average($countHalos->($selectConversion)/$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}->($selectConversion));
	    my $countHalosModel = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'}*$conversion;
	    (my $selected, my $unselected) =
		which(
		    ($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'            } >  $massHaloMinimum)
		    &
		    ($simulation                                                      ->{'countHaloMinimum'} <= $countHalosModel)
		);
	    my $selectedModel =
		which(
		    ($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'mass'            } >  $massHaloMinimum)
		);
	    my $overdensityLow  = -1.50;
	    my $overdensityHigh = +1.10;
	    my $fraction;
	    if ( $overdensityEnvironment      < $overdensityLow  ) {
		$fraction = 0.0;
	    } elsif ( $overdensityEnvironment > $overdensityHigh ) {
		$fraction = 1.0;
	    } else {
		$fraction = sclr(($overdensityEnvironment-$overdensityLow)/($overdensityHigh-$overdensityLow));
	    }
	    my @colorLowDensity  = ( 237.0, 0.83, 0.94 );
	    my @colorHighDensity = (  13.0, 0.83, 0.94 ); 
	    my $color = &GnuPlot::PrettyPlots::Color_Gradient($fraction,\@colorLowDensity,\@colorHighDensity);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                                               ,
		$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'        }->($selected)                     ,
		$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}->($selected)                     
		/$massFunctionTarget,
		style        => "line"                                                                                  ,
		style        => "point"                                                                                              ,
		weight       => [2,1]                                                                                                ,
		symbol       => [6,6]                                                                                                ,
		pointSize    => 0.5                                                                                                  ,
		color        => [$color,$color]
	    	);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                                  ,
		$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'mass'        }->($selectedModel)    ,
	 	$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunction'}->($selectedModel)
		/$massFunctionModel,
		style        => "line"                                                                                  ,
		weight       => [3,1]                                                                                   ,
		color        => [$color,$color]
		);
	}
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
    }
    
    unless ( -e $outputDirectory."/haloMassFunctionResidualsFractional_".$redshiftLabel.".pdf" ) {
	## Halo mass function residuals, fractional.
	my $plot;
	my $gnuPlot;
	my $plotFileTeX = $outputDirectory."/haloMassFunctionResidualsFractional_".$redshiftLabel.".tex";
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
	print $gnuPlot "set output '".$plotFileTeX."'\n";
	print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
	print $gnuPlot "set ylabel '\$ [(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{model}-(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{target}]/(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{model} \$\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.65,0.49\n";
	print $gnuPlot "set logscale x\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [8.0e6:1.0e16]\n";
	print $gnuPlot "set yrange [-1.0:+1.0]\n";
	print $gnuPlot "set pointsize 1.0\n";
	foreach my $simulation ( @simulations ) {
	    next
		unless ( scalar(@{$simulation->{'expansionFactors'}}) > $i );
	    my $massHaloMinimum = $simulation->{'massParticle'}*$countParticlesMinimum;
	    my $xLimit          = pdl [ $massHaloMinimum, $massHaloMinimum ];
	    my $yLimit          = pdl [ -1.0            , 1.0              ];
	    # Iterate over realizations.
	    my $massTarget;
	    my $massFunctionTarget;
	    my $massFunctionTargetError;
	    my $massModel;
	    my $massFunctionModel;
	    my $massFunctionModelBinned;
	    foreach my $realization ( @{$simulation->{'realizations'}} ) {
		my $realizationLabel = $realization eq "" ? "" : "_".$realization;
		my $countHalos =
		    (
		     +$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'     }
		     /$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}
		    )**2;
		my $selectConversion = which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'        } > $massHaloMinimum)
		    );
		my $conversion = average($countHalos->($selectConversion)/$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}->($selectConversion));
		my $countHalosModel = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'}*$conversion;
		(my $selected, my $unselected) =
		    which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'            } >  $massHaloMinimum)
		    );
		unless ( defined($massTarget) ) {
		    $massTarget               = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'                    }->($selected)   ;
		    $massFunctionTarget       = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'            }->($selected)   ;
		    $massFunctionTargetError  = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'       }->($selected)**2;
		    $massModel                = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'mass'                    }->($selected)   ;
		    $massFunctionModel        = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunction'            }->($selected)   ;
		    $massFunctionModelBinned  = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionInterpolated'}->($selected)   ;
		} else {
		    $massFunctionTarget      += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'            }->($selected)   ;
		    $massFunctionTargetError += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'       }->($selected)**2;
		    $massFunctionModel       += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunction'            }->($selected)   ;
		    $massFunctionModelBinned += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionInterpolated'}->($selected)   ;
		}
	    }
	    $massFunctionTargetError .= $massFunctionTargetError->sqrt();
	    $massFunctionTarget      /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionTargetError /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionModel       /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionModelBinned /= scalar(@{$simulation->{'realizations'}});
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                             ,
		$xLimit                                                            ,
		$yLimit                                                            ,
		style        => "line"                                             ,
		weight       => [1,1]                                              ,
		linePattern  => 3                                                  ,
		color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                                        ,
		$massTarget,
		(
		 +$massFunctionModelBinned
		 -$massFunctionTarget
		)
		/$massFunctionModelBinned,
		style        => "point"                                            ,
		weight       => [2,1]                                              ,
		symbol       => [1,1]                                              ,
		pointSize    => 0.2                                                ,
		color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		);
	}
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
    }

    unless ( -e $outputDirectory."/haloMassFunctionResidualsNormalized_".$redshiftLabel.".pdf" ) {
	## Halo mass function residuals, normalized.
	my $plot;
	my $gnuPlot;
	my $plotFileTeX = $outputDirectory."/haloMassFunctionResidualsNormalized_".$redshiftLabel.".tex";
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
	print $gnuPlot "set output '".$plotFileTeX."'\n";
	print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
	print $gnuPlot "set ylabel '\$ |(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{model}-(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{target}|/\\sigma_\\mathrm{model} \$\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.65,0.49\n";
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [8.0e6:1.0e16]\n";
	print $gnuPlot "set yrange [0.1:1.0e3]\n";
	print $gnuPlot "set pointsize 1.0\n";
	foreach my $simulation ( @simulations ) {
	    next
		unless ( scalar(@{$simulation->{'expansionFactors'}}) > $i );
	    my $massHaloMinimum = $simulation->{'massParticle'}*$countParticlesMinimum;
	    my $xLimit          = pdl [ $massHaloMinimum, $massHaloMinimum ];
	    my $yLimit          = pdl [ 0.1             , 1.0e3            ];
	    # Iterate over realizations.
	    my $massTarget;
	    my $massFunctionTarget;
	    my $massFunctionModelError;
	    my $massModel;
	    my $massFunctionModel;
	    my $massFunctionModelBinned;
	    foreach my $realization ( @{$simulation->{'realizations'}} ) {
		my $realizationLabel = $realization eq "" ? "" : "_".$realization;
		my $countHalos =
		    (
		     +$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'     }
		     /$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunctionError'}
		    )**2;
		my $selectConversion = which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'        } > $massHaloMinimum)
			&
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'} > 0.0             )
		    );
		my $conversion = average($countHalos->($selectConversion)/$simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'}->($selectConversion));
		my $countHalosModel = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model'}->{'massFunctionInterpolated'}*$conversion;
		(my $selected, my $unselected) =
		    which(
			($simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'            } >  $massHaloMinimum)
		    );
		unless ( defined($massTarget) ) {
		    $massTarget               = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'mass'                    }->($selected)   ;
		    $massFunctionTarget       = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'            }->($selected)   ;
		    $massModel                = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'mass'                    }->($selected)   ;
		    $massFunctionModel        = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunction'            }->($selected)   ;
		    $massFunctionModelBinned  = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionInterpolated'}->($selected)   ;
		    $massFunctionModelError   = $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionInterpolated'}->($selected)**2/$countHalosModel->($selected);
		} else {
		    $massFunctionTarget      += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'target'}->{'massFunction'            }->($selected)   ;
		    $massFunctionModel       += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunction'            }->($selected)   ;
		    $massFunctionModelBinned += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionInterpolated'}->($selected)   ;
		    $massFunctionModelError  += $simulation->{$realizationLabel}->{'redshifts'}->[$i]->{'model' }->{'massFunctionInterpolated'}->($selected)**2/$countHalosModel->($selected);
		}
	    }
	    $massFunctionModelError  .= $massFunctionModelError->sqrt();
	    $massFunctionTarget      /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionModelError  /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionModel       /= scalar(@{$simulation->{'realizations'}});
	    $massFunctionModelBinned /= scalar(@{$simulation->{'realizations'}});
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                             ,
		$xLimit                                                            ,
		$yLimit                                                            ,
		style        => "line"                                             ,
		weight       => [1,1]                                              ,
		linePattern  => 3                                                  ,
		color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                                           ,
		$massTarget,
		abs(
		    +$massFunctionTarget
		    -$massFunctionModelBinned
		)
		/$massFunctionModelError,
		style        => "point"                                            ,
		weight       => [2,1]                                              ,
		symbol       => [1,1]                                              ,
		pointSize    => 0.2                                                ,
		color        => [$simulation->{'colorDark'},$simulation->{'color'}]
		);
	}
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
    }

}

exit 0;

sub zoomInsPlotModify {
    # Modify plotting arguments for the halo mass function.
    my $simulation      = shift();
    my $plotOptions     = shift();
    my $realization     = shift();
    my $expansionFactor = shift();
    my $redshiftLabel   = sprintf("%5.3f",1.0/$expansionFactor-1.0);
    # Appending environment to title.
    if ( defined($realization) ) {
	my $target                   = new PDL::IO::HDF5($ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$simulation->{'label'}."_".$realization."_z".$redshiftLabel.".hdf5");                                    
	my $targetSimulation         = $target          ->group  ('simulation0001'        );                                                                                                                                  
	(my $massEnvironment       ) = $targetSimulation->attrGet('massEnvironment'       );                                                                                                                                  
	(my $overdensityEnvironment) = $targetSimulation->attrGet('overdensityEnvironment');
	$plotOptions->{'title'   } = "\\\\tiny ".$plotOptions->{'title'}.sprintf("; ".$realization."; \$z = ".$redshiftLabel."; \\\\log_{10}(M_\\\\mathrm{env}/\\\\mathrm{M}_\\\\odot)=%5.2f; \\\\delta_\\\\mathrm{env}=%+6.3f",log10($massEnvironment),$overdensityEnvironment);
    }
    # Set plotting limits.
    $plotOptions->{'xMinimum'} = 8.0e+06;
    $plotOptions->{'xMaximum'} = 3.0e+12;
    $plotOptions->{'yMinimum'} = 3.0e-04;
    $plotOptions->{'yMaximum'} = 1.0e+02;
    # Set key location.
    $plotOptions->{'xKey'    } = 0.975;
    $plotOptions->{'yKey'    } = 0.900;
}

sub zoomInsValidate {
    # Validate parameters for the halo mass function.
    my $simulation      = shift();
    my $parameters      = shift();
    my $realization     = shift();
    my $expansionFactor = shift();
    my $redshiftLabel   = sprintf("%5.3f",1.0/$expansionFactor-1.0);
    # Check environment properties match those in the halo mass function file
    if ( defined($realization) ) {
	my $target                   = new PDL::IO::HDF5($ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$simulation->{'label'}."_".$realization."_z".$redshiftLabel.".hdf5");
	my $targetSimulation         = $target          ->group  ('simulation0001'        );                                                                                                                                  
	(my $massEnvironment       ) = $targetSimulation->attrGet('massEnvironment'       );                                                                                                                                  
	(my $overdensityEnvironment) = $targetSimulation->attrGet('overdensityEnvironment');
	if ( abs($massEnvironment       -$parameters->{'haloEnvironment'}->{'massEnvironment'}->{'value'}) > 1.0e-3*$massEnvironment ) {
	    print "environment mass in parameter file does not match that in halo mass function file for ".$simulation->{'label'}." ".$realization." z=".$redshiftLabel."\n";
	    print "           parameter file: ".$parameters->{'haloEnvironment'}->{'massEnvironment'}->{'value'}."\n";
	    print "  halo mass function file: ".                                   $massEnvironment             ."\n";
	    exit 1;
	}
	if ( abs($overdensityEnvironment-$parameters->{'haloEnvironment'}->{'overdensity'    }->{'value'}) > 1.0e-3                  ) {
	    print "environment overdensity in parameter file does not match that in halo mass function file for ".$simulation->{'label'}." ".$realization." z=".$redshiftLabel."\n";
	    print "           parameter file: ".$parameters->{'haloEnvironment'}->{'overdensity'}->{'value'}."\n";
	    print "  halo mass function file: ".                                   $overdensityEnvironment  ."\n";
	    exit 1;
	}
    }
}
