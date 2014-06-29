#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::Fit::Polynomial;
use PDL::IO::HDF5;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Find relation between mass and maximum distance for the sample of Bernardi et al. (2013;
# http://adsabs.harvard.edu/abs/2013MNRAS.436..697B).  Andrew Benson (24-April-2014)

# Tabulated data provided by Mariangela Bernardi.
my $log10MassStellar = pdl # Stellar masses are in units of M_Solar.
    (
     8.655397415,
     8.758553505,
     8.856801987,
     8.954696655,
     9.052925110,
     9.151350975,
     9.255033493,
     9.351667404,
     9.454163551,
     9.553510666,
     9.651016235,
     9.750605583,
     9.852743149,
     9.954152107,
     10.053371429,
     10.153795242,
     10.253108025,
     10.352433205,
     10.452750206,
     10.551837921,
     10.651668549,
     10.751879692,
     10.849768639,
     10.950572014,
     11.049318314,
     11.147963524,
     11.247805595,
     11.346734047,
     11.446274757,
     11.545268059,
     11.644490242,
     11.742430687,
     11.840696335,
     11.941007614,
     12.043179512,
     12.135819435
    );
my $volumeMaximum = pdl # Comoving volumes are in units of 10^9 Mpc^3.
    (
     0.002602330,
     0.002865336,
     0.003689888,
     0.005141599,
     0.006540476,
     0.008913630,
     0.011576199,
     0.015452851,
     0.020556541,
     0.027066843,
     0.035638396,
     0.047748027,
     0.063546862,
     0.082955151,
     0.106538044,
     0.134253938,
     0.161168912,
     0.196770098,
     0.232025141,
     0.281275746,
     0.349784664,
     0.444777516,
     0.562816249,
     0.717034707,
     0.912019113,
     1.143364514,
     1.425395784,
     1.767075533,
     2.164316387,
     2.631458218,
     3.202964305,
     3.810579128,
     4.508823099,
     5.098926730,
     5.946447509,
     6.690833525
    );

# Solid angle of the sample.
system($galacticusPath."scripts/aux/mangleRansack.pl ".$galacticusPath."constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/sdss_dr72safe0_res6d.pol ".$galacticusPath."constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/solidAngle.hdf5 0");
my $solidAngleFile = new PDL::IO::HDF5($galacticusPath."constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/solidAngle.hdf5");
my $solidAngles    = $solidAngleFile->dataset('solidAngle')->get();
my $solidAngle     = $solidAngles->sum();

# Convert volumes to maximum distance.
my $log10DistanceMaximum = log10((3.0*$volumeMaximum*1.0e9/$solidAngle)**(1.0/3.0));

# Fit a polynomial.
(my $fit, my $coeffs) = fitpoly1d($log10MassStellar,$log10DistanceMaximum,6);
my $fitMass     = sequence(1000)*3.6/999.0+8.6;
my $fitDistance = pdl zeroes(nelem($fitMass));
my $closing     = "";
for(my $i=0;$i<nelem($coeffs);++$i) {
    $fitDistance += $coeffs->(($i))*($fitMass**$i);
    print $coeffs->(($i))."d0";
    if ( $i == nelem($coeffs)-1 ) {
	print $closing."\n";
    } else {
	print "+logarithmicMass*(";
	$closing .= ")";
    }
}

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = $galacticusPath."constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/massDistanceRelation.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Maximum Distance vs. Limiting Stellar Mass for Bernardi et al. (2013) Sample'\n";
print $gnuPlot "set xlabel 'Limiting stellar mass; \$M_\\star [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Maximum distance; \$D_{\\rm max}\$ [Mpc]'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.16,0.79\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [3.0e8:2.0e12]\n";
print $gnuPlot "set yrange [100.0:2500.0]\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$log10MassStellar,
    10.0**$log10DistanceMaximum,
    style      => "point",
    weight     => [5,3],
    symbol     => [6,7],
    color      => $PrettyPlots::colorPairs{'mediumSeaGreen'},
    title      => 'Bernardi et al. (2013)'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$fitMass    ,
    10.0**$fitDistance,
    style      => "line",
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{'redYellow'},
    title      => 'fit'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
