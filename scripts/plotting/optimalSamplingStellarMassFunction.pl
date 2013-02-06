#!/usr/bin/env perl
use strict;
use warnings;
use lib "./perl";
use PDL;
use XML::Simple;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use Data::Dumper;

# Make a plot of the tree sampling function.
# Andrew Benson (20-November-2010)

# Compile the sampling code.
system("make Optimal_Sampling.Stellar_Mass_Function.exe");

# Create the tree timing file.
my $treeTiming;
$treeTiming->{'fit'}->{'coefficient'} = [-19.64, 2.57, -0.066];
my $xmlOut = new XML::Simple (NoAttr=>1, RootName=>"timing");
open(oHndl,">treeTimingFit.xml");
print oHndl $xmlOut->XMLout($treeTiming);
close(oHndl);

# Run the sampling code.
system("Optimal_Sampling.Stellar_Mass_Function.exe");

# Read data from the tree weighting file.
my $haloMass = pdl [];
my $gamma    = pdl [];
my $haloMF   = pdl [];
my $timeTree = pdl [];
open(iHndl,"treeSampling.data");
while ( my $line = <iHndl> ) {
    $line =~ s/^\s*//;
    $line =~ s/\s*$//;
    my @columns = split(/\s+/,$line);
    $haloMass = $haloMass->append($columns[1]);
    $gamma    = $gamma   ->append($columns[2]);
    $haloMF   = $haloMF  ->append($columns[3]);
    $timeTree = $timeTree->append($columns[4]);
}
close(iHndl);

# Normalize to unit time.
my $deltaLogM  = log($haloMass->index(1)/$haloMass->index(0));
my $gammaTime  = sum($gamma *$timeTree*$deltaLogM);
my $haloMFTime = sum($haloMF*$timeTree*$deltaLogM);
$gamma        /= $gammaTime;
$haloMF       /= $haloMFTime;

# Declare variables for GnuPlot;
my ($gnuPlot, $plot);

# Open a pipe to GnuPlot.
my $plotFileEPS = "plots/optimalSamplingStellarMassFunction.eps";
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.4,0.2\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e10:1.0e15]\n";
print $gnuPlot "set yrange [1.0e-7:3.0]\n";
print $gnuPlot "set title 'Optimal sampling density'\n";
print $gnuPlot "set xlabel '\$M_{\\rm halo} [M_\\odot]\$'\n";
print $gnuPlot "set ylabel '\$\\gamma(M_{\\rm halo})\$'\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(\$plot,
			      $haloMass,$haloMF,
			      style     => "line",
			      weight    => [5,3],
			      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
			      title     => 'Halo mass function weight'
    );
&PrettyPlots::Prepare_Dataset(\$plot,
			      $haloMass,$gamma,
			      style     => "line",
			      weight    => [5,3],
			      color     => $PrettyPlots::colorPairs{'redYellow'},
			      title     => 'Optimal weight'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

# Clean up.
unlink("treeTimingFit.xml","Optimal_Sampling.Stellar_Mass_Function.exe","stellarMassFunction.data","treeSampling.data");

exit;
