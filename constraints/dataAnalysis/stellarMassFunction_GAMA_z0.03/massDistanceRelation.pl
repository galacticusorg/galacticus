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
use PDL::Fit::Polynomial;
use XML::Simple;
use Data::Dumper;
require Stats::Percentiles;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Determine the relation between stellar mass and limiting distance for the GAMA stellar mass function.
# Andrew Benson (04-June-2014)

# Define a working directory.
my $workDirectory = $galacticusPath."constraints/dataAnalysis/stellarMassFunction_GAMA_z0.03/";

# Define survey solid angle.
my $squareDegreesPerSteradian = pdl 3282.80635;
my $solidAngle                = pdl 143.0/$squareDegreesPerSteradian;

# Read the stellar mass function.
my $xml = new XML::Simple();
my $stellarMassFunction = $xml->XMLin($galacticusPath."data/observations/massFunctionsStellar/Stellar_Mass_Function_GAMA_2012.xml");
my $logarithmicMass = pdl $stellarMassFunction->{'massFunction'}->{'columns'}->{'mass'        }->{'datum'};
my $massFunction    = pdl $stellarMassFunction->{'massFunction'}->{'columns'}->{'massFunction'}->{'datum'};
my $number          = pdl $stellarMassFunction->{'massFunction'}->{'columns'}->{'number'      }->{'datum'};
my $binWidth        = pdl 0.2*ones(nelem($massFunction));
$binWidth->(0:1)   .= 0.5;

# Compute the effective volume.
my $volume = $number/$massFunction/$binWidth;

# Define the fields.
my @fields =
    (
     {
	 label      => "G09/G15",
	 depth      => 19.4,
	 fieldCount => 2
     },
     {
	 label      => "G12",
	 depth      => 19.8,
	 fieldCount => 1
     }
    );

# Find relative volumes of the fields.
my $volumeTotal = 0.0;
foreach ( @fields ) {
    $_->{'distance'} = 10.0**(0.4*($_->{'depth'}-19.0));
    $_->{'volume'  } = $_->{'distance'}**3;
    $volumeTotal += $_->{'fieldCount'}*$_->{'volume'};
}
foreach ( @fields ) {
    $_->{'volume'  } /= $volumeTotal;
}

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = $workDirectory."massDistanceRelation.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Limiting Distance vs. Stellar Mass for GAMA Survey'\n";
print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Limiting distance; \$D_{\\rm max}\$ [Mpc]'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.575,0.26\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e6:1.0e12]\n";
print $gnuPlot "set yrange [20.0:500]\n";
print $gnuPlot "set pointsize 2.0\n";

# Iterate over fields.
my $color = -1;
foreach my $field ( @fields ) {
    ++$color;

    # Find the maximum distance.
    my $distance            = (3.0*$volume*$field->{'volume'}/($solidAngle/3.0))**(1.0/3.0);
    my $logarithmicDistance = log10($distance);

    # Fit a polynomial to the results.
    my $fitTo             = which($logarithmicMass < 9.0);
    (my $fit, my $coeffs) = fitpoly1d($logarithmicMass->($fitTo),$logarithmicDistance->($fitTo),2);

    # Generate a fit to the data.
    my $fitMass     = sequence(1000)*6.0/999.0+6.0;
    my $fitDistance = pdl zeroes(nelem($fitMass));
    my $closing     = "";
    print $field->{'label'}.":\t";
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

    &PrettyPlots::Prepare_Dataset(
	\$plot,
	10.0**$fitMass,
	10.0**$fitDistance,
	style      => "line",
	weight     => [5,3],
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence2'}}[$color]},
	title      => 'Fit ['.$field->{'label'}.']'
	);
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	10.0**$logarithmicMass,
	10.0**$logarithmicDistance,
	style      => "point",
	weight     => [5,3],
	symbol     => [6,7],
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[$color]},
	title      => 'GAMA ['.$field->{'label'}.']'
	);

}

&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
