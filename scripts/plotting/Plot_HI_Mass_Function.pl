#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
require Galacticus::HDF5;
require Stats::Histograms;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_HI_Gas_Mass_Function.pl <galacticusFile> <outputDir/File> [<showFit>]")};
$self           = $0;
$galacticusFile = $ARGV[0];
$outputTo       = $ARGV[1];
if ( $#ARGV == 2 ) {
    $showFit    = $ARGV[2];
    if ( lc($showFit) eq "showfit"   ) {$showFit = 1};
    if ( lc($showFit) eq "noshowfit" ) {$showFit = 0};
} else {
    $showFit = 0;
}

# Check if output location is file or directory.
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/HI_Gas_Mass_Function.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);

# Read the XML data file.
$xml = new XML::Simple;
$data = $xml->XMLin($galacticusPath."data/observations/massFunctionsHI/HI_Mass_Function_Zwaan_2005.xml");
$columns = $data->{'massFunction'}->{'columns'};
$xBins = pdl @{$columns->{'mass'}->{'data'}};
$x = pdl @{$columns->{'mass'}->{'data'}};
$y = pdl @{$columns->{'massFunction'}->{'data'}};
$errorUp = pdl @{$columns->{'upperError'}->{'data'}};
$errorDown = pdl @{$columns->{'lowerError'}->{'data'}};
$errorUp   = (+10.0**($y+$errorUp)  -10.0**$y)*($dataSet->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**3;
$errorDown = (-10.0**($y-$errorDown)+10.0**$y)*($dataSet->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**3;
$x         = (10.0**$x             )*($dataSet->{'parameters'}->{'H_0'}/$columns->{'mass'}->{'hubble'})**2;
$y         = (10.0**$y             )*($dataSet->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**3;
$xBins     = $xBins+log10(           ($dataSet->{'parameters'}->{'H_0'}/$columns->{'mass'}->{'hubble'})**2);

# Read galaxy data and construct mass function.
$yGalacticus = zeroes nelem($xBins);
$errorGalacticus = zeroes nelem($xBins);
$binStep = $xBins->index(1)-$xBins->index(0);
$binMin = $xBins->index(0)-0.5*$binStep;
$binMax = $xBins->index(nelem($xBins)-1)+0.5*$binStep;
# Factor to convert cold gas mass to HI mass from Power, Baugh & Lacey (2009; http://adsabs.harvard.edu/abs/2009arXiv0908.1396P).
$gasMassToHIMassFactor = pdl 0.54;
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['mergerTreeWeight','diskMassGas','spheroidMassGas']);
$dataSets           = $dataSet->{'dataSets'};
$logarithmicMassGas = log10(($dataSets->{'diskMassGas'}+$dataSets->{'spheroidMassGas'})*$gasMassToHIMassFactor);
$weight             = $dataSets->{'mergerTreeWeight'};
delete($dataSet->{'dataSets'});
($yGalacticus,$errorGalacticus) = &Histograms::Histogram($xBins,$logarithmicMassGas,$weight,differential => 1);

# Compute chi^2.
$chiSquared = sum(($yGalacticus-$y)**2/($errorGalacticus**2+(0.5*($errorUp-$errorDown))**2));
$degreesOfFreedom = nelem($y);
if ( $showFit == 1 ) {
    $fitData{'name'} = "Zwaan et al. (2005) HI gas mass function";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make plot of stellar mass function.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'HI Gas Mass Function at \$z=0\$'\n";
print $gnuPlot "set xlabel 'Galaxy HI mass; \$M_{\\rm HI} [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Comoving number density; \${\\rm d}n/{\\rm d}\\ln M_{\\rm HI}(M_{\\rm HI}) [\\hbox{Mpc}^{-3}]\$'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e7:1.0e11]\n";
print $gnuPlot "set yrange [1.0e-6:1.0e0]\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$y,
    errorUp    => $errorUp,
    errorDown  => $errorDown,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[0]},
    title      => $data->{'massFunction'}->{'label'}.' [observed]'
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
