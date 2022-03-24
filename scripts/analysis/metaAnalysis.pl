#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use GnuPlot::LaTeX;
use GnuPlot::PrettyPlots;
use Data::Dumper;

# Create plots of collected metadata on merger tree ODE evolver.
# Andrew Benson (05-June-2011)

# Get arguments.
die ("Usage: metaAnalysis.pl <modelFile> <outputFolder>")
    unless ( scalar(@ARGV) == 2 );
my $modelFile    = $ARGV[0];
my $outputFolder = $ARGV[1];

# Open the file containing the meta-data.
my $HDFfile = new PDL::IO::HDF5($modelFile);

# Read the timestep histogram.
my $foundGroup      = 0;
my @groupsAvailable = $HDFfile->groups();
foreach my $groupAvailable ( @groupsAvailable ) {
    $foundGroup = 1 if ( $groupAvailable eq "metaData" );
}
die ("metaAnalysis.pl: metaData group does not exist") unless ( $foundGroup == 1 );
$foundGroup      = 0;
@groupsAvailable = $HDFfile->group("metaData")->groups();
foreach my $groupAvailable ( @groupsAvailable ) {
    $foundGroup = 1 if ( $groupAvailable eq "evolverProfiler" );
}
die ("metaAnalysis.pl: metaData/evolverProfile group does not exist") unless ( $foundGroup == 1 );
my $timeSteps        = $HDFfile->group("metaData/evolverProfiler")->dataset("timeStep"        )->get();
my $timeStepCount    = $HDFfile->group("metaData/evolverProfiler")->dataset("timeStepCount"   )->get();
my $propertyNames    = $HDFfile->group("metaData/evolverProfiler")->dataset("propertyNames"   )->get();
my $propertyHitCount = $HDFfile->group("metaData/evolverProfiler")->dataset("propertyHitCount")->get();
my $propertyHitRate  = 100.0*$propertyHitCount/$propertyHitCount->sum();
my $propertyIndex    = pdl 1..nelem($propertyHitRate);

# Convert count into frequency.
$timeStepCount = float($timeStepCount)/$timeStepCount->sum();

# Declare variables for GnuPlot;
my ($gnuPlot, $outputFile, $outputFileEPS, $plot);

# Open a pipe to GnuPlot and make a plot of halo mass to stellar mass ratio vs. stellar mass.
$outputFile = $outputFolder."/timesteps.pdf";
($outputFileEPS = $outputFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$outputFileEPS."'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.2,0.2\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e-9:10.0]\n";
print $gnuPlot "set yrange [1.0e-6:1.0]\n";
print $gnuPlot "set title 'Timestep distribution'\n";
print $gnuPlot "set xlabel 'Timestep [Gyr]'\n";
print $gnuPlot "set ylabel 'Frequency'\n";
&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
			      $timeSteps,$timeStepCount,
			      style => "point", symbol => [6,7], weight => [5,3],
			      color => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'});
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($outputFileEPS);

# Open a pipe to GnuPlot and make a plot of halo mass to stellar mass ratio vs. stellar mass.
$outputFile = $outputFolder."/hitRates.pdf";
($outputFileEPS = $outputFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$outputFileEPS."'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.2,0.2\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
my $xExtent = nelem($propertyHitCount)+1;
print $gnuPlot "set xrange [0:".$xExtent."]\n";
my $yExtent = maximum($propertyHitRate)*1.05;
print $gnuPlot "set yrange [0:".$yExtent."]\n";
print $gnuPlot "set title 'Hit Frequency'\n";
print $gnuPlot "set xlabel 'Property'\n";
print $gnuPlot "set ylabel 'Hit frequency [\\%]'\n";
&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
			      $propertyIndex,$propertyHitRate,
			      style => "boxes", symbol => [6,7], weight => [5,3],
			      color => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'});
for (my $iProperty=0;$iProperty<nelem($propertyHitCount);++$iProperty) {
    print $gnuPlot "set label '".$propertyNames->atstr($iProperty)."' at ".$propertyIndex->index($iProperty).",0 front rotate by 90\n";
}
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($outputFileEPS);

exit;
