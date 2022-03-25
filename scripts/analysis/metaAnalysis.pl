#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
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
my $timeSteps                  = $HDFfile->group("metaData/evolverProfiler")->dataset("timeStep"                  )->get();
my $timeStepCount              = $HDFfile->group("metaData/evolverProfiler")->dataset("timeStepCount"             )->get();
my $evaluationCount            = $HDFfile->group("metaData/evolverProfiler")->dataset("evaluationCount"           )->get();
my $timeCPU                    = $HDFfile->group("metaData/evolverProfiler")->dataset("timeCPU"                   )->get();
my $timeStepCountInterrupted   = $HDFfile->group("metaData/evolverProfiler")->dataset("timeStepCountInterrupted"  )->get();
my $evaluationCountInterrupted = $HDFfile->group("metaData/evolverProfiler")->dataset("evaluationCountInterrupted")->get();
my $timeCPUInterrupted         = $HDFfile->group("metaData/evolverProfiler")->dataset("timeCPUInterrupted"        )->get();
my $propertyNames              = $HDFfile->group("metaData/evolverProfiler")->dataset("propertyNames"             )->get();
my $propertyHitCount           = $HDFfile->group("metaData/evolverProfiler")->dataset("propertyHitCount"          )->get();
my $propertyHitRate            = 100.0*$propertyHitCount/$propertyHitCount->sum();
my $propertyIndex              = pdl 1..nelem($propertyHitRate);

# Convert counts into frequencies.
my $timeStepFrequency                    = float(  $timeStepCount           )/sum(  $timeStepCount           );
my $timeStepWeightedFrequency            = float($evaluationCount           )/sum($evaluationCount           );
my $timeCPUFrequency                     = float(  $timeCPU                 )/sum(  $timeCPU                 );
my $timeStepInterruptedFrequency         = float(  $timeStepCountInterrupted)/sum(  $timeStepCountInterrupted);
my $timeStepInterruptedWeightedFrequency = float($evaluationCountInterrupted)/sum($evaluationCountInterrupted);
my $timeCPUInterruptedFrequency          = float(  $timeCPUInterrupted      )/sum(  $timeCPUInterrupted      );

# Find cumulative frequencies.
my $timeStepFrequencyCumulative                    = $timeStepFrequency                   ->cumusumover()/$timeStepFrequency                   ->sum();
my $timeStepWeightedFrequencyCumulative            = $timeStepWeightedFrequency           ->cumusumover()/$timeStepWeightedFrequency           ->sum();
my $timeCPUFrequencyCumulative                     = $timeCPUFrequency                    ->cumusumover()/$timeCPUFrequency                    ->sum();
my $timeStepInterruptedFrequencyCumulative         = $timeStepInterruptedFrequency        ->cumusumover()/$timeStepInterruptedFrequency        ->sum();
my $timeStepInterruptedWeightedFrequencyCumulative = $timeStepInterruptedWeightedFrequency->cumusumover()/$timeStepInterruptedWeightedFrequency->sum();
my $timeCPUInterruptedFrequencyCumulative          = $timeCPUInterruptedFrequency         ->cumusumover()/$timeCPUInterruptedFrequency         ->sum();

# Find ranges for timestep plot.
my $frequencies        = $timeStepFrequency->append($timeStepWeightedFrequency)->append($timeStepInterruptedFrequencyCumulative)->append($timeStepInterruptedWeightedFrequencyCumulative);
my $frequenciesNonZero = which($frequencies > 0.0);
my $timestepMinimum    = 10.0**(floor(log10(        $timeSteps  ->(( 0)               ) )));
my $timestepMaximum    = 10.0**(ceil (log10(        $timeSteps  ->((-1)               ) )));
my $frequencyMinimum   = 10.0**(floor(log10(minimum($frequencies->($frequenciesNonZero)))));
my $frequencyMaximum   = 10.0**(ceil (log10(maximum($frequencies->($frequenciesNonZero)))));

# Declare variables for GnuPlot;
my ($gnuPlot, $outputFile, $outputFileTeX, $plot);

# Open a pipe to GnuPlot and make a plot of timestep frequency.
$outputFileTeX = $outputFolder."/timesteps.tex";
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 6in,4in\n";
print $gnuPlot "set output '".$outputFileTeX."'\n";
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
print $gnuPlot "set xrange [". $timestepMinimum.":". $timestepMaximum."]\n";
print $gnuPlot "set yrange [".$frequencyMinimum.":".$frequencyMaximum."]\n";
print $gnuPlot "set title offset 0,-0.8 'Timestep distribution'\n";
print $gnuPlot "set xlabel 'Timestep [Gyr]'\n";
print $gnuPlot "set ylabel 'Frequency'\n";
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepFrequency,
     style  => "point",
     symbol => [6,7],
     weight => [3,1],
     title  => "Steps",
     color  => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepWeightedFrequency,
     style  => "point",
     symbol => [6,7],
     weight => [3,1],
     title  => "Evaluations",
     color  => $GnuPlot::PrettyPlots::colorPairs{'mediumSeaGreen'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeCPUFrequency,
     style  => "point",
     symbol => [6,7],
     weight => [3,1],
     title  => "CPU",
     color  => $GnuPlot::PrettyPlots::colorPairs{'indianRed'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepInterruptedFrequency,
     style  => "point",
     symbol => [6,6],
     weight => [3,1],
     title  => "Steps (interrupted)",
     color  => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepInterruptedWeightedFrequency,
     style  => "point",
     symbol => [6,6],
     weight => [3,1],
     title  => "Evaluations (interrupted)",
     color  => $GnuPlot::PrettyPlots::colorPairs{'mediumSeaGreen'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepInterruptedWeightedFrequency,
     style  => "point",
     symbol => [6,6],
     weight => [3,1],
     title  => "Evaluations (interrupted)",
     color  => $GnuPlot::PrettyPlots::colorPairs{'mediumSeaGreen'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeCPUInterruptedFrequency,
     style  => "point",
     symbol => [6,6],
     weight => [3,1],
     title  => "CPU (interrupted)",
     color  => $GnuPlot::PrettyPlots::colorPairs{'indianRed'}
    );
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($outputFileTeX);

# Open a pipe to GnuPlot and make a plot of cumulative timestep frequency.
$outputFileTeX = $outputFolder."/timestepsCumulative.tex";
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 6in,4in\n";
print $gnuPlot "set output '".$outputFileTeX."'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.2,0.6\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale x\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [".$timestepMinimum.":".$timestepMaximum."]\n";
print $gnuPlot "set yrange [-0.02:1.02]\n";
print $gnuPlot "set title offset 0,-0.8 'Cumulative timestep distribution'\n";
print $gnuPlot "set xlabel 'Timestep [Gyr]'\n";
print $gnuPlot "set ylabel 'Cumulative fraction'\n";
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepFrequencyCumulative,
     style  => "point",
     symbol => [6,7],
     weight => [3,1],
     title  => "Steps",
     color  => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepWeightedFrequencyCumulative,
     style  => "point",
     symbol => [6,7],
     weight => [3,1],
     title  => "Evaluations",
     color  => $GnuPlot::PrettyPlots::colorPairs{'mediumSeaGreen'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeCPUFrequencyCumulative,
     style  => "point",
     symbol => [6,7],
     weight => [3,1],
     title  => "CPU",
     color  => $GnuPlot::PrettyPlots::colorPairs{'indianRed'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepInterruptedFrequencyCumulative,
     style  => "point",
     symbol => [6,6],
     weight => [3,1],
     title  => "Steps (interrupted)",
     color  => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeStepInterruptedWeightedFrequencyCumulative,
     style  => "point",
     symbol => [6,6],
     weight => [3,1],
     title  => "Evaluations (interrupted)",
     color  => $GnuPlot::PrettyPlots::colorPairs{'mediumSeaGreen'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset
    (
     \$plot,
     $timeSteps,
     $timeCPUInterruptedFrequencyCumulative,
     style  => "point",
     symbol => [6,6],
     weight => [3,1],
     title  => "CPU (interrupted)",
     color  => $GnuPlot::PrettyPlots::colorPairs{'indianRed'}
    );
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($outputFileTeX);

# Open a pipe to GnuPlot and make a plot of property hit rates.
$outputFileTeX = $outputFolder."/hitRates.tex";
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 6in,4in\n";
print $gnuPlot "set output '".$outputFileTeX."'\n";
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
print $gnuPlot "set title offset 0,-0.8 'Hit Frequency'\n";
print $gnuPlot "set xlabel 'Property'\n";
print $gnuPlot "set ylabel 'Hit frequency [\\%]'\n";
&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
			      $propertyIndex,$propertyHitRate,
			      style => "boxes", symbol => [6,7], weight => [3,1],
			      color => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'});
for (my $iProperty=0;$iProperty<nelem($propertyHitCount);++$iProperty) {
    print $gnuPlot "set label '".$propertyNames->atstr($iProperty)."' at ".$propertyIndex->index($iProperty).",0 front rotate by 90\n";
}
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($outputFileTeX);

exit;
