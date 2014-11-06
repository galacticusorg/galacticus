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
use PDL::IO::Misc;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Fit stellar mass completeness limits from the ULTRAVISTA survey (Muzzin et al. 2013). Fits the tabulated results in the data file
# downloaded from ULTRAVISTA web site.
# Andrew Benson (13-August-2014)

# Specify work directory.
my $workDirectory = $galacticusPath."constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/";

# Get the data file if we do not already have it.
system("wget http://www.strw.leidenuniv.nl/galaxyevolution/ULTRAVISTA/Mstar_redshift_completeness_emp_uvista_v4.1_95.dat -O ".$workDirectory."Mstar_redshift_completeness_emp_uvista_v4.1_95.dat")
    unless ( -e $workDirectory."Mstar_redshift_completeness_emp_uvista_v4.1_95.dat" );

# Read the data file.
(my $redshift, my $mass) = rcols($workDirectory."Mstar_redshift_completeness_emp_uvista_v4.1_95.dat",0,1);

# Find unique masses.
my $unique = which($redshift < 3.0); #$mass->uniqind();

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = $workDirectory."massRedshiftRelation.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Redshift vs. Limiting Stellar Mass for Muzzin et al. (2013) Sample'\n";
print $gnuPlot "set xlabel 'Limiting stellar mass; \$M_\\star [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Redshift; \$z\$'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.74,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale x\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [3.0e8:3.0e11]\n";
print $gnuPlot "set yrange [0.2:4.0]\n";
print $gnuPlot "set pointsize 2.0\n";
# Fit limiting redshifts.
(my $fit, my $coeffs) = fitpoly1d($mass->($unique),$redshift->($unique),6);
my $fitMass     = sequence(1000)*5.0/999.0+7.0;
my $fitRedshift = pdl zeroes(nelem($fitMass));
my $closing     = "";
for(my $i=0;$i<nelem($coeffs);++$i) {
    $fitRedshift += $coeffs->(($i))*($fitMass**$i);
    print $coeffs->(($i))."d0";
    if ( $i == nelem($coeffs)-1 ) {
	print $closing."\n";
    } else {
	print "+logarithmicMass*(";
	$closing .= ")";
    }
}
$fitRedshift .= $fitRedshift/(1.0-exp(($fitMass-11.24)/0.02));
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$fitMass,
    $fitRedshift,
    style      => "line",
    weight     => [3,1],
    color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[0]},
    title      => 'fit'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$mass,
    $redshift,
    style      => "line",
    weight     => [5,3],
    linePattern => 3,
    color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[0]},
    title      => 'observed'
    );    
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
