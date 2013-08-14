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
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
require Stats::Histograms;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;

# Get name of input and output files.
die("Plot_bJ_Luminosity_Function.pl <galacticusFile> <outputDir/File> [<showFit>]")
    unless ( scalar(@ARGV) == 2 || scalar(@ARGV) == 3 );
my $self           = $0;
my $galacticusFile = $ARGV[0];
my $outputTo       = $ARGV[1];
my $showFit;
if ( scalar(@ARGV) == 3 ) {
    $showFit    = $ARGV[2];
    if ( lc($showFit) eq "showfit"   ) {$showFit = 1};
    if ( lc($showFit) eq "noshowfit" ) {$showFit = 0};
} else {
    $showFit = 0;
}

# Check if output location is file or directory.
my $outputFile;
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/bJ_Luminosity_Function.pdf";
}
(my $fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/$1/;

# Create data structure to read the results.
my $dataSet;
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);

# Read the XML data file.
my $xml     = new XML::Simple;
my $data    = $xml->XMLin($galacticusPath."data/observations/luminosityFunctions/bJ_Luminosity_Function_2dFGRS_Norberg_2002.xml");
my $columns = $data->{'luminosityFunction'}->{'columns'};
my $xBins   = pdl @{$columns->{'magnitude'}->{'data'}};
my $x       = pdl @{$columns->{'magnitude'}->{'data'}};
my $y       = pdl @{$columns->{'luminosityFunction'}->{'data'}};
my $error   = pdl @{$columns->{'error'}->{'data'}};
$xBins   = $xBins-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'});
$x       = $x-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'});
$y       = $y*($dataSet->{'parameters'}->{'H_0'}/$columns->{'luminosityFunction'}->{'hubble'})**3;
$error   = $error*($dataSet->{'parameters'}->{'H_0'}/$columns->{'luminosityFunction'}->{'hubble'})**3;

# Reverse the order of the vectors.
$xBins = $xBins(-1:0);
$x     = $x    (-1:0);
$y     = $y    (-1:0);
$error = $error(-1:0);

# Read galaxy data and construct luminosity function.
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['mergerTreeWeight','magnitudeTotal:bJ:rest:z0.0000:vega','magnitudeTotal:bJ:rest:z0.0000:dustAtlas:vega']);
my $dataSets      = $dataSet->{'dataSets'};
my $magnitude     = $dataSets->{'magnitudeTotal:bJ:rest:z0.0000:dustAtlas:vega'};
my $magnitudeFree = $dataSets->{'magnitudeTotal:bJ:rest:z0.0000:vega'};
my $weight        = $dataSets->{'mergerTreeWeight'};
delete($dataSet->{'dataSets'});
(my $yGalacticus    , my $errorGalacticus    ) = &Histograms::Histogram($xBins,$magnitude    ,$weight,differential => 1);
(my $yGalacticusFree, my $errorGalacticusFree) = &Histograms::Histogram($xBins,$magnitudeFree,$weight,differential => 1);

# Compute chi^2.
my $chiSquared = sum(($yGalacticus-$y)**2/($errorGalacticus**2+$error**2));
my $degreesOfFreedom = nelem($y);
if ( $showFit == 1 ) {
    my %fitData;
    $fitData{'name'} = "Norberg et al. (2002) bJ-band luminosity function";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make plot of stellar mass function.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot"); # 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'b\$_{\\rm J}\$-band Luminosity Function at \$z=0\$'\n";
print $gnuPlot "set xlabel 'b\$_{\\rm J}\$-band absolute magnitude; \$M_{\\rm b_J}\$'\n";
print $gnuPlot "set ylabel 'Comoving number density; \${\\rm d}n/{\\rm d}M_{\\rm b_J}(M_{\\rm b_J}) [\\hbox{Mpc}^{-3}]\$'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale y\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [-24.5:-14.5]\n";
print $gnuPlot "set yrange [1.0e-8:3.0e-2]\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$y,
    errorUp    => $error,
    errorDown  => $error,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[0]},
    title      => $data->{'luminosityFunction'}->{'label'}.' [observed]'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$yGalacticusFree,
    errorUp    => $errorGalacticusFree,
    errorDown  => $errorGalacticusFree,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{'peachPuff'},
    title      => 'Galacticus (no dust extinction)'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$yGalacticus,
    errorUp    => $errorGalacticus,
    errorDown  => $errorGalacticus,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{'redYellow'},
    title      => 'Galacticus'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);
&MetaData::Write($plotFile,$galacticusFile,$self);

exit;
