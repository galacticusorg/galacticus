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
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::Fit::Polynomial;
use XML::Simple;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Construct a mass error model for the ALFALFA survey.
# Andrew Benson (10-July-2013)

# Define working directory.
my $workDirectory = "constraints/dataAnalysis/hiMassFunction_ALFALFA_z0.00/";

# Load the data.
my $xml                  = new XML::Simple;
my $data                 = $xml->XMLin($workDirectory."alfalfaHIMassErrorModel.xml");
my $logarithmicMass      = pdl @{$data->{'mass'     }->{'datum'}};
my $logarithmicMassError = pdl @{$data->{'massError'}->{'datum'}};

# Generate a fit to the data. Functional form chosen by simple eyeball examination of the data.
my $aBest          = -1.0;
my $bBest = 5.8;
my $cBest = 0.5;
my $fitMeasureBest = 1.0e30;
for(my $a=0.05;$a<=0.15;$a+=0.005) {
    for(my $b=5.00;$b<=7.00;$b+=0.005) {
	for(my $c=0.25;$c<=0.75;$c+=0.005) {
	    my $fit        = $a+exp(-($logarithmicMass-$b)/$c);
	    my $fitMeasure = sum(($fit-$logarithmicMassError)**2);
	    if ( $fitMeasure < $fitMeasureBest ) {
		$aBest          = $a;
		$bBest          = $b;
		$cBest          = $c;
		$fitMeasureBest = $fitMeasure;
	    }
	}
    }
}
print "Parameters of best fitting model:\n";
print "  a = ".$aBest."\n";
print "  b = ".$bBest."\n";
print "  c = ".$cBest."\n";
my $massFine = pdl sequence(1000)*7.0/1000.0+5.0;
my $fitFine  = $aBest+exp(-($massFine-$bBest)/$cBest);

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = $workDirectory."alfalfaHIMassErrorModel.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'HI mass error model for ALFALFA 40\\\% (\$\\alpha.40\$) sample'\n";
print $gnuPlot "set xlabel 'HI gas mass; \$M_{\\rm HI} [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Error in \$\\log_{10}\$ of HI gas mass; \$\\sigma\$ []'\n";
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
print $gnuPlot "set xrange [3.1e5:3.1e11]\n";
print $gnuPlot "set yrange [0.0:0.8]\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$logarithmicMass,
    $logarithmicMassError,
    style      => "point",
    weight     => [5,3],
    symbol     => [6,7],
    pointSize  => 2.0,
    color      => $PrettyPlots::colorPairs{'peachPuff'},
    title      => 'Haynes et al. (2011)'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$massFine,
    $fitFine,
    style      => "line",
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{'mediumSeaGreen'},
    title      => 'Fit'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
