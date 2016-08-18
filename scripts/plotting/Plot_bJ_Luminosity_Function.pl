#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Galacticus::Path;
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
use Stats::Histograms;
use Galacticus::HDF5;
use Galacticus::Magnitudes;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use XMP::MetaData;

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
&Galacticus::HDF5::Get_Parameters($dataSet);
&Galacticus::HDF5::Count_Trees($dataSet);
&Galacticus::HDF5::Select_Output($dataSet,0.0);

# Read the XML data file.
my $xml     = new XML::Simple;
my $data    = $xml->XMLin(&galacticusPath()."data/observations/luminosityFunctions/bJ_Luminosity_Function_2dFGRS_Norberg_2002.xml");
my $columns = $data->{'luminosityFunction'}->{'columns'};
my $xBins   = pdl @{$columns->{'magnitude'}->{'data'}};
my $x       = pdl @{$columns->{'magnitude'}->{'data'}};
my $y       = pdl @{$columns->{'luminosityFunction'}->{'data'}};
my $error   = pdl @{$columns->{'error'}->{'data'}};
$xBins   = $xBins-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant'}->{'value'});
$x       = $x    -5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant'}->{'value'});
$y       = $y    *($dataSet->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant'}->{'value'}/$columns->{'luminosityFunction'}->{'hubble'})**3;
$error   = $error*($dataSet->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant'}->{'value'}/$columns->{'luminosityFunction'}->{'hubble'})**3;

# Reverse the order of the vectors.
$xBins = $xBins(-1:0);
$x     = $x    (-1:0);
$y     = $y    (-1:0);
$error = $error(-1:0);

# Read galaxy data and construct luminosity function.
$dataSet->{'tree'} = "all";
&Galacticus::HDF5::Get_Dataset($dataSet,['mergerTreeWeight','magnitudeTotal:bJ:rest:z0.0000:vega','magnitudeTotal:bJ:rest:z0.0000:dustAtlas:vega']);
my $dataSets      = $dataSet->{'dataSets'};
my $magnitude     = $dataSets->{'magnitudeTotal:bJ:rest:z0.0000:dustAtlas:vega'};
my $magnitudeFree = $dataSets->{'magnitudeTotal:bJ:rest:z0.0000:vega'};
my $weight        = $dataSets->{'mergerTreeWeight'};
delete($dataSet->{'dataSets'});
(my $yGalacticus    , my $errorGalacticus    ) = &Stats::Histograms::Histogram($xBins,$magnitude    ,$weight,differential => 1);
(my $yGalacticusFree, my $errorGalacticusFree) = &Stats::Histograms::Histogram($xBins,$magnitudeFree,$weight,differential => 1);

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
(my $plotFileTeX = $plotFile) =~ s/\.pdf$/.tex/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
print $gnuPlot "set output '".$plotFileTeX."'\n";
print $gnuPlot "set title offset 0,-1 'b\$_\\mathrm{J}\$-band Luminosity Function at \$z\\approx0\$'\n";
print $gnuPlot "set xlabel '\$M_\\mathrm{b_J}\$'\n";
print $gnuPlot "set ylabel '\$\\mathrm{d}n/\\mathrm{d}M_\\mathrm{b_J}(M_\\mathrm{b_J})\\,\\,[\\hbox{Mpc}^{-3}]\$'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.5,0.2\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale y\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [-24.5:-14.5]\n";
print $gnuPlot "set yrange [1.0e-8:3.0e-2]\n";
print $gnuPlot "set pointsize 2.0\n";
&GnuPlot::PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$y,
    errorUp    => $error,
    errorDown  => $error,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    pointSize  => 0.5,
    color      => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
    title      => $data->{'luminosityFunction'}->{'label'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$yGalacticus,
    errorUp    => $errorGalacticus,
    errorDown  => $errorGalacticus,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    pointSize  => 0.5,
    color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
    title      => 'Galacticus'
    );
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
&XMP::MetaData::Write($plotFile,$galacticusFile,$self);

exit;
