#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::Fit::Polynomial;
use Astro::Cosmology;
require Stats::Percentiles;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Determine the relation between stellar mass and Spitzer IRAC 4.5um apparent magnitude using semi-analytic models from the
# Millennium Database. Specifically, the Henriques2012a models (http://adsabs.harvard.edu/abs/2012MNRAS.421.2904H) are used -
# these correspond to Guo et al. (2011; http://adsabs.harvard.edu/abs/2011MNRAS.413..101G) which use a Chabrier IMF. We therefore
# include a correction to the Salpeter IMF assumed by Caputi et al. (2011).
# Andrew Benson (21-August-2012)

# Define a working directory.
my $workDirectory = $galacticusPath."constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/massLuminosityWork/";
system("mkdir -p ".$workDirectory);

# Define little Hubble parameter for the Millennium Simulation.
my $hubble = pdl 0.7;
# Define Millennium Simulation cosmology.
my $cosmology = Astro::Cosmology->new(
    omega_matter => 0.3,
    omega_lambda => 0.7,
    H0           => 70.0
    );
 
# Define the database URL.
my $databaseURL = "http://gavo.mpa-garching.mpg.de/MyMillennium?action=doQuery&SQL=";

# Define magnitude range for sample.
my $chabrierToSalpeter = pdl 1.8;
my $magnitudeLimit     = pdl 24.0;
my $magnitudeWidth     = pdl  0.2;
my $magnitudeMinimum   = $magnitudeLimit-0.5*$magnitudeWidth-2.5*log10($chabrierToSalpeter);
my $magnitudeMaximum   = $magnitudeLimit+0.5*$magnitudeWidth-2.5*log10($chabrierToSalpeter);

# Define the SQL query.
my $sqlQuery = "select light.z_app, mass.stellarmass, light.i2 from Guo2010a..MR as mass, Henriques2012a.wmap1.BC03_001 as light where mass.galaxyid = light.galid and light.i2 < ".$magnitudeMaximum." and light.i2 > ".$magnitudeMinimum." and light.z_app > 2.5 and light.z_app < 5.5";

# Download the data.
system("wget --http-user=abenson --http-passwd=n70KZIVc \"".$databaseURL.$sqlQuery."\" -O ".$workDirectory."millenniumDB.csv")
    unless ( -e $workDirectory."millenniumDB.csv" );

# Read the data.
(my $redshift, my $mass, my $magnitude) = rcols($workDirectory."millenniumDB.csv",{EXCLUDE => '/(^#|stellarMass)/', COLSEP => ","});

# Convert mass to logarithmic and scale to Solar units.
my $logMass = log10(1.0e10*$mass/$hubble);

# Define weights.
my $weight = pdl ones(nelem($logMass));

# Define redshift bins.
my $redshiftMinimum = pdl 2.5;
my $redshiftMaximum = pdl 5.5;
my $redshiftCount   = 10;
my $redshiftBins    = sequence($redshiftCount)*($redshiftMaximum-$redshiftMinimum)/($redshiftCount-1)+$redshiftMinimum;

# Find the median stellar mass as a function of redshift.
my $percentiles  = pdl ( 50.0 );
my $medianMassV  = &Percentiles::BinnedPercentiles($redshiftBins,$redshift,$logMass,$weight,$percentiles);
my $medianMass = $medianMassV->(:,(0));

# Generate a fit to the data.
my $nonZero = which($medianMass > 0.0);
(my $fit, my $coeffs) = fitpoly1d($medianMass->index($nonZero),$redshiftBins->index($nonZero),2);
my $fitMass     = sequence(1000)*5.0/999.0+8.0;
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

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = $workDirectory."massLuminosityRelation.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Redshift vs. Limiting Stellar Mass for Caputi et al. (2011) Sample'\n";
print $gnuPlot "set xlabel 'Limiting stellar mass; \$M_\\star [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Redshift; \$z\$'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.19,0.8\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale x\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [6.0e9:6.0e10]\n";
print $gnuPlot "set yrange [2.5:5.5]\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$fitMass,
    $fitRedshift,
    style      => "line",
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{'mediumSeaGreen'},
    title      => 'Fit'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$medianMass,
    $redshiftBins,
    style      => "point",
    weight     => [5,3],
    symbol     => [6,7],
    color      => $PrettyPlots::colorPairs{'redYellow'},
    title      => 'Guo et al. (2011) SAM'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
