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
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Fit stellar mass completeness limits from the PRIMUS survey (Moustakas et al. 2013).
# Andrew Benson (24-April-2014)

# Completeness limit data for "All" galaxies taken from Table 2 of Moustakas et al. (2013).
my @fields =
    (
     {
	 name     => "COSMOS"                                          ,
	 redshift => pdl (0.250, 0.350,  0.450,  0.575,  0.725,  0.900),
	 mass     => pdl (8.730, 9.140,  9.510,  9.920, 10.330, 10.710)
     },
     {
	 name     => "XMM-SXDS"                                        ,
	 redshift => pdl (0.250, 0.350,  0.450,  0.575,  0.725,  0.900),
	 mass     => pdl (8.860, 9.230,  9.580,  9.970, 10.380, 10.780)
     },
     {
	 name     => "XMM-CFHTLS"                                      ,
	 redshift => pdl (0.250, 0.350,  0.450,  0.575,  0.725,  0.900),
	 mass     => pdl (8.950, 9.230,  9.510,  9.870, 10.310, 10.830)
     },
     {
	 name     => "CDFS"                                            ,
	 redshift => pdl (0.250, 0.350,  0.450,  0.575,  0.725,  0.900),
	 mass     => pdl (9.620, 9.870, 10.100, 10.370, 10.650, 10.940 )
     },
     {
	 name     => "ELAIS-S1"                                        ,
	 redshift => pdl (0.250, 0.350,  0.450,  0.575,  0.725,  0.900),
	 mass     => pdl (9.700, 9.990, 10.260, 10.560, 10.870, 11.170)
     }
    );

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = "redshifts.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Redshift vs. Limiting Stellar Mass for Moustakas et al. (2013) Sample'\n";
print $gnuPlot "set xlabel 'Limiting stellar mass; \$M_\\star [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Redshift; \$z\$'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.34,0.56\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale x\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [3.0e8:1.0e11]\n";
print $gnuPlot "set yrange [0.0:1.0]\n";
print $gnuPlot "set pointsize 2.0\n";

# Fit limiting redshifts.
my $iField = -1;
foreach my $field ( @fields ) {
    ++$iField;
    (my $fit, my $coeffs) = fitpoly1d($field->{'mass'},$field->{'redshift'},3);
    my $fitMass     = sequence(1000)*5.0/999.0+8.0;
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
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	10.0**$fitMass,
	$fitRedshift,
	style      => "line",
	weight     => [5,3],
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[$iField]},
	title      => $field->{'name'}.' [fit]'
    );
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	10.0**$field->{'mass'},
	$field->{'redshift'},
	style      => "point",
	weight     => [5,3],
	symbol     => [6,7],
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[$iField]},
	title      => $field->{'name'}.' [PRIMUS]'
	);
    
}
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
