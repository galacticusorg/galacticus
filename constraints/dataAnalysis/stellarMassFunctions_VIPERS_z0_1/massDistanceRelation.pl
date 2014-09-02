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
use Astro::Cosmology;
require Stats::Percentiles;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Determine the relation between stellar mass and limiting distance for the VIPERS stellar mass functions.
# Andrew Benson (04-June-2014)

# Define a working directory.
my $workDirectory = $galacticusPath."constraints/dataAnalysis/stellarMassFunctions_VIPERS_z0_1/";

# Define survey solid angle (computed from mangle polygons).
my $solidAngle = pdl 0.003137; #0.004143678167878;

# Read the stellar mass function.
my $xml = new XML::Simple();
my $stellarMassFunction = $xml->XMLin($galacticusPath."data/observations/massFunctionsStellar/Stellar_Mass_Functions_VIPERS_2013.xml");

# Construct cosmology.
my $cosmologyObserved 
    = Astro::Cosmology->new(
    omega_matter => $stellarMassFunction->{'cosmology'}->{'omegaMatter'    }, 
    omega_lambda => $stellarMassFunction->{'cosmology'}->{'omegaDarkEnergy'}, 
    h0           => $stellarMassFunction->{'cosmology'}->{'hubble'         }
    );

# Estimated sampling rates for each redshift.
my $samplingRate = pdl [ 0.44, 0.44, 0.4 ];

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = $workDirectory."massDistanceRelation.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Limiting Distance vs. Stellar Mass for VIPERS Survey'\n";
print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Limiting distance; \$D_{\\rm max}\$ [Mpc]'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.35,0.7\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [3.0e9:3.0e12]\n";
print $gnuPlot "set yrange [2000:4000]\n";
print $gnuPlot "set pointsize 2.0\n";

# Iterate over redshifts.
my $redshiftBin = -1;
for my $columns ( @{$stellarMassFunction->{'massFunction'}->{'columns'}} ) {
    ++$redshiftBin;

    # Extract mass function.
    my $logarithmicMass = pdl $columns->{'mass'        }->{'datum'};
    my $massFunction    = pdl $columns->{'massFunction'}->{'datum'};
    my $number          = pdl $columns->{'number'      }->{'datum'};
    my $binWidth        = pdl 0.2*ones(nelem($massFunction));
    my $redshiftMinimum = pdl $columns->{'redshiftLow' };
    my $redshiftMaximum = pdl $columns->{'redshiftHigh'};
    my $redshiftLabel   = sprintf("\$%3.1f<z<%3.1f\$",$redshiftMinimum->sclr(),$redshiftMaximum->sclr());

    # Find distance to minimum and maximum redshifts.
    my $distanceMinimum = $cosmologyObserved->comoving_distance($redshiftMinimum);
    my $distanceMaximum = $cosmologyObserved->comoving_distance($redshiftMaximum);

    # Compute the effective volume.
    my $volume = $number/$samplingRate->(($redshiftBin))/$massFunction/$binWidth;
 
    # Find the maximum distance.
    my $distance            = (3.0*$volume/$solidAngle+$distanceMinimum**3)**(1.0/3.0);
    my $logarithmicDistance = log10($distance);

    # Fit a polynomial to the results.
    my $fitTo             = which(($massFunction > 0.0) & ($distance < 0.98*$distanceMaximum) & ($logarithmicMass < 12));
    (my $fit, my $coeffs) = fitpoly1d($logarithmicMass->($fitTo),$logarithmicDistance->($fitTo),2);

    # Generate a fit to the data.
    my $fitMass     = sequence(1000)*9.0/999.0+4.0;
    my $fitDistance = pdl zeroes(nelem($fitMass));
    my $closing     = "";
    print $redshiftMinimum." < z < ".$redshiftMaximum.": ";
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
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence2'}}[$redshiftBin]},
	title      => 'Fit ['.$redshiftLabel.']'
	);
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	10.0**$logarithmicMass,
	10.0**$logarithmicDistance,
	style      => "point",
	weight     => [5,3],
	symbol     => [6,7],
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[$redshiftBin]},
	title      => 'VIPERS ['.$redshiftLabel.']'
	);
}

&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
