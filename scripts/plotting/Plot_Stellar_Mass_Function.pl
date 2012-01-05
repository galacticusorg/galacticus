#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
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

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Stellar_Mass_Function.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
    $outputFile = $outputTo."/Stellar_Mass_Function.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.1);

# Read the XML data file.
$xml       = new XML::Simple;
$data      = $xml->XMLin($galacticusPath."data/Stellar_Mass_Function_Li_White_2009.xml");
$columns   = $data->{'stellarMassFunction'}->{'columns'};
$xBins     = pdl @{$columns->{'stellarMass' }->{'datum'}};
$x         = pdl @{$columns->{'stellarMass' }->{'datum'}};
$y         = pdl @{$columns->{'massFunction'}->{'datum'}};
$errorUp   = pdl @{$columns->{'upperError'  }->{'datum'}};
$errorDown = pdl @{$columns->{'lowerError'  }->{'datum'}};
$errorUp   = (+10.0**($y+$errorUp  )-10.0**$y)*($dataSet->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;
$errorDown = (-10.0**($y-$errorDown)+10.0**$y)*($dataSet->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;
$xBins     = $xBins+log10(           ($dataSet->{'parameters'}->{'H_0'}/$columns->{'stellarMass' }->{'hubble'})**$columns->{'stellarMass' }->{'hubbleExponent'});
$x         = (10.0**$x             )*($dataSet->{'parameters'}->{'H_0'}/$columns->{'stellarMass' }->{'hubble'})**$columns->{'stellarMass' }->{'hubbleExponent'} ;
$y         = (10.0**$y             )*($dataSet->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;

# Read galaxy data and construct mass function.
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['volumeWeight','diskStellarMass','spheroidStellarMass']);
$dataSets               = $dataSet->{'dataSets'};
$logarithmicStellarMass = log10(($dataSets->{'diskStellarMass'}+$dataSets->{'spheroidStellarMass'}));
$weight                 = $dataSets->{'volumeWeight'};
delete($dataSet->{'dataSets'});
($yGalacticus,$errorGalacticus) = &Histograms::Histogram($xBins,$logarithmicStellarMass,$weight,differential => 1);

# Compute chi^2.
$chiSquared = sum(($yGalacticus-$y)**2/($errorGalacticus**2+(0.5*($errorUp-$errorDown))**2));
$degreesOfFreedom = nelem($y);
if ( $showFit == 1 ) {
    $fitData{'name'} = "Li and White (2009) stellar mass function";
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
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Stellar Mass Function at \$z\\approx0.1\$'\n";
print $gnuPlot "set xlabel 'Galaxy stellar mass; \$M_\\star [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Comoving number density; \${\\rm d}n/{\\rm d}\\ln M_\\star(M_\\star) [\\hbox{Mpc}^{-3}]\$'\n";
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
print $gnuPlot "set xrange [1.0e8:3.0e12]\n";
print $gnuPlot "set yrange [1.0e-8:1.0e-1]\n";
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
	 title      => $data->{'stellarMassFunction'}->{'label'}.' [observed]'
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

exit;
