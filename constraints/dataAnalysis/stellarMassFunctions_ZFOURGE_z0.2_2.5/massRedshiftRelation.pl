#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::Fit::Polynomial;
use PDL::IO::Misc;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Fit stellar mass completeness limits from the ZFOURGE survey (Tomczak et al. 2014). Fits the tabulated results in the data file
# given by R. Quadri.
# Andrew Benson (11-August-2014)

# Read data.
my @fields = 
    (
     {
	 name => "NMBS"
     },
     {
	 name => "ZFOURGE"
     }
    );
(my $redshift, my $massNMBS, my $massZFOURGE) = rcols("zfourge-SMF-supplemental/masslimits.dat",0,1,2);
$fields[0]->{'mass'    } = $massNMBS;
$fields[1]->{'mass'    } = $massZFOURGE;
$fields[0]->{'redshift'} = $redshift;
$fields[1]->{'redshift'} = $redshift;

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = "constraints/dataAnalysis/stellarMassFunctions_ZFOURGE_z0.2_2.5/massRedshiftRelation.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Redshift vs. Limiting Stellar Mass for Tomczak et al. (2014) Sample'\n";
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
print $gnuPlot "set xrange [3.0e7:3.0e11]\n";
print $gnuPlot "set yrange [0.2:3.1]\n";
print $gnuPlot "set pointsize 2.0\n";

# Fit limiting redshifts.
my $iField = -1;
foreach my $field ( @fields ) {
    ++$iField;
    (my $fit, my $coeffs) = fitpoly1d($field->{'mass'},$field->{'redshift'},4);
    my $fitMass     = sequence(1000)*5.0/999.0+7.0;
    my $fitRedshift = pdl zeroes(nelem($fitMass));
    my $closing     = "";
    print $field->{'name'}." : ";
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
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	10.0**$fitMass,
	$fitRedshift,
	style      => "line",
	weight     => [3,1],
	color      => $GnuPlot::PrettyPlots::colorPairs{${$GnuPlot::PrettyPlots::colorPairSequences{'sequence1'}}[$iField]},
	title      => $field->{'name'}.' [fit]'
    );
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	10.0**$field->{'mass'},
	$field->{'redshift'},
	style      => "line",
	weight     => [5,3],
	linePattern => 3,
	color      => $GnuPlot::PrettyPlots::colorPairs{${$GnuPlot::PrettyPlots::colorPairSequences{'sequence1'}}[$iField]},
	title      => $field->{'name'}.' [observed]'
	);    
}
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
