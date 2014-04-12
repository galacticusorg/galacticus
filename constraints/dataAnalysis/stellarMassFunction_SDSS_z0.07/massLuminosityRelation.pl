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

# Determine the relation between stellar mass and SDSS r-band absolute magnitude using semi-analytic models from the Millennium
# Database. Specifically, the DeLucia2006a models are used - these correspond to De Lucia & Blaizot (2007;
# http://adsabs.harvard.edu/abs/2007MNRAS.375....2D) which use a Chabrier IMF, consistent with the IMF used by Li & White (2009).
# Andrew Benson (10-July-2012)

# Define a working directory.
my $workDirectory = $galacticusPath."constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/massLuminosityWork/";
system("mkdir -p ".$workDirectory);

# Define little Hubble parameter for the Millennium Simulation.
my $hubble = pdl 0.7;

# Define the limiting apparent magnitude of the Li & White (2009) sample.
my $apparentMagnitudeLimit = pdl 17.6;

# Define Millennium Simulation cosmology.
my $cosmology = Astro::Cosmology->new(
    omega_matter => 0.3,
    omega_lambda => 0.7,
    H0           => 70.0
    );
 
# Define the database URL.
my $databaseURL = "http://gavo.mpa-garching.mpg.de/Millennium/?action=doQuery&SQL=";

# Get list of snapshot numbers and corresponding redshifts.
my $sqlQuery = "select snapnum, z from millimil..Snapshots";
system("wget \"".$databaseURL.$sqlQuery."\" -O ".$workDirectory."millenniumDB_snapshots.csv")
    unless ( -e $workDirectory."millenniumDB_snapshots.csv" );
(my $snapshotNumbers, my $redshifts) = rcols($workDirectory."millenniumDB_snapshots.csv",{EXCLUDE => '/(^#|snapnum)/', COLSEP => ","});

# Initialize vectors to store results.
my $redshiftTable     = pdl ( 0.0 );
my $limitingMassTable = pdl ( 6.0 );

# Loop over redshifts.
for(my $iSnapshot=nelem($snapshotNumbers)-1;$iSnapshot>=0;--$iSnapshot) {
    my $redshift = $redshifts->(($iSnapshot));

    # Only consider reasonably low redshifts.
    if ( $redshift <= 0.5 && $redshift > 0.0 ) {

	# Compute the limiting absolute magnitude at this redshift.
	my $limitingMagnitude = $cosmology->absolute_magnitude($apparentMagnitudeLimit,$redshift);

	# Define the SQL query.
	my $sqlQuery = "select mass.stellarMass, light.r_sdss from millimil..DeLucia2006a as mass, millimil..DeLucia2006a_sdss2mass as light where mass.snapnum = ".$iSnapshot." and mass.galaxyID = light.galaxyID";

	# Download the data.
	system("wget \"".$databaseURL.$sqlQuery."\" -O ".$workDirectory."millenniumDB_".$iSnapshot.".csv")
	    unless ( -e $workDirectory."millenniumDB_".$iSnapshot.".csv" );

	# Read the data.
	(my $mass, my $magnitude) = rcols($workDirectory."millenniumDB_".$iSnapshot.".csv",{EXCLUDE => '/(^#|stellarMass)/', COLSEP => ","});

	# Convert mass to logarithmic scale and Solar units.
	my $logMass = log10(1.0e10*$mass/$hubble);

	# Define weights.
	my $weight = pdl ones(nelem($logMass));

	# Define mass bins.
	my $logMassMinimum = pdl  8.0;
	my $logMassMaximum = pdl 13.0;
	my $logMassCount   = 10;
	my $logMassBins    = sequence($logMassCount)*($logMassMaximum-$logMassMinimum)/($logMassCount-1)+$logMassMinimum;

	# Find the median magnitude as a function of mass.
	my $percentiles = pdl ( 50.0 );
	my $median      = &Percentiles::BinnedPercentiles($logMassBins,$logMass,$magnitude,$weight,$percentiles);

	# Fit a polynomial to the results.
	my $nonZeroBins = which($median->(:,(0)) < 0.0);
	(my $fit, my $coeffs) = fitpoly1d($median->(:,(0))->index($nonZeroBins),$logMassBins->index($nonZeroBins),3);

	# Compute the limiting mass at this redshift.
	my $limitingMass = 
	    +$coeffs->((0))
	    +$coeffs->((1))*$limitingMagnitude
	    +$coeffs->((2))*$limitingMagnitude**2;

	# Store the results.
	$redshiftTable     = $redshiftTable    ->append($redshift    );
	$limitingMassTable = $limitingMassTable->append($limitingMass);

    }
}

# Generate a fit to the data.
(my $fit, my $coeffs) = fitpoly1d($limitingMassTable,$redshiftTable,5);
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
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Redshift vs. Limiting Stellar Mass for Li \\& White (2009) Sample'\n";
print $gnuPlot "set xlabel 'Limiting stellar mass; \$M_\\star [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Redshift; \$z\$'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.76\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale x\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e8:1.0e13]\n";
print $gnuPlot "set yrange [0.0:0.5]\n";
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
    10.0**$limitingMassTable,
    $redshiftTable,
    style      => "point",
    weight     => [5,3],
    symbol     => [6,7],
    color      => $PrettyPlots::colorPairs{'redYellow'},
    title      => 'De Lucia (2006) SAM'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
