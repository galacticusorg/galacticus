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
use XML::Simple;
use Astro::Cosmology;
use Math::SigFigs;
use Data::Dumper;
require Stats::Means;
require Galacticus::HDF5;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;

# Get name of input and output files.
die("Plot_Star_Formation_History.pl <galacticusFile> <outputDir/File> [<showFit>]")
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

# Get parameters for the Galacticus model.
my $dataBlock;
$dataBlock->{'file'} = $galacticusFile;
&HDF5::Get_Parameters($dataBlock);

# Check if output location is file or directory.
my $outputFile;
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/Star_Formation_History.pdf";
}
(my $fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/$1/;

# Extract global data
&HDF5::Get_History($dataBlock,['historyExpansion','historyStarFormationRate']);
my $history        = $dataBlock->{'history'};
my $time           = $history->{'historyTime'};
my $redshift       = 1.0/$history->{'historyExpansion'}-1.0;
my $SFR            = $history->{'historyStarFormationRate'}/1.0e9;
my $stellarDensity = $history->{'historyStellarDensity'};

# Determine IMF correction factor.
my $imfCorrection = 1.0;
if ( $dataBlock->{'parameters'}->{'imfSelectionFixed'} eq "Chabrier" ) {
    $imfCorrection = 0.6;
}

# Define redshift bins.
my $redshiftPoints = pdl 16;
my $redshiftMin    = pdl 0.0;
my $redshiftMax    = pdl 8.0;
my $redshiftBin    = pdl ($redshiftMax-$redshiftMin)/$redshiftPoints;
my $redshiftBins   = pdl (0..$redshiftPoints-1)*$redshiftBin+$redshiftMin+0.5*$redshiftBin;

# Read the XML data file.
my $xml = new XML::Simple;
my $data = $xml->XMLin($galacticusPath."data/observations/starFormationRate/Star_Formation_Rate_Data.xml");
my $iDataset = -1;
my $chiSquared = 0.0;
my $degreesOfFreedom = 0;
my @dataSets;
foreach my $dataSet ( @{$data->{'starFormationRate'}} ) {
    my $columns = $dataSet->{'columns'};
    my $x = pdl @{$columns->{'redshift'}->{'data'}};
    my $xLowerError = pdl @{$columns->{'redshiftErrorDown'}->{'data'}};
    my $xUpperError = pdl @{$columns->{'redshiftErrorUp'}->{'data'}};
    $xLowerError = $x-$xLowerError;
    $xUpperError = $x+$xUpperError;
    my $y = pdl @{$columns->{'sfr'}->{'data'}};
    my $yLowerError = pdl @{$columns->{'sfrErrorDown'}->{'data'}};
    my $yUpperError = pdl @{$columns->{'sfrErrorUp'}->{'data'}};
    $yUpperError = (10.0**($y+$yUpperError));
    $yLowerError = (10.0**($y-$yLowerError));
    $y = (10.0**$y);

    # Compute cosmology corrections.
    my $cosmologyData = Astro::Cosmology->new(
	omega_matter => $columns->{'sfr'}->{'omega'},
	omega_lambda => $columns->{'sfr'}->{'lambda'},
	H0           => $columns->{'sfr'}->{'hubble'}
	);
    my $cosmologyGalacticus = Astro::Cosmology->new(
	omega_matter => $dataBlock->{'parameters'}->{'Omega_Matter'},
	omega_lambda => $dataBlock->{'parameters'}->{'Omega_DE'},
	H0           => $dataBlock->{'parameters'}->{'H_0'}
	);

    my $volumeElementData            = $cosmologyData      ->differential_comoving_volume($x);
    my $volumeElementGalacticus      = $cosmologyGalacticus->differential_comoving_volume($x);
    my $luminosityDistanceData       = $cosmologyData      ->luminosity_distance         ($x);
    my $luminosityDistanceGalacticus = $cosmologyGalacticus->luminosity_distance         ($x);
    my $cosmologyCorrection = ($volumeElementData/$volumeElementGalacticus)*($luminosityDistanceGalacticus/$luminosityDistanceData)**2;
    
    # Apply cosmology corrections.
    $y           = $y          *$cosmologyCorrection;
    $yUpperError = $yUpperError*$cosmologyCorrection;
    $yLowerError = $yLowerError*$cosmologyCorrection;

    # Apply IMF correction.
    $y           = $y          *$imfCorrection;
    $yUpperError = $yUpperError*$imfCorrection;
    $yLowerError = $yLowerError*$imfCorrection;

    # Store the dataset.
    ++$iDataset;
    $dataSets[$iDataset]->{'x'}           =               $x;
    $dataSets[$iDataset]->{'xLowerError'} = -$xLowerError+$x;
    $dataSets[$iDataset]->{'xUpperError'} = +$xUpperError-$x;
    $dataSets[$iDataset]->{'y'}           =               $y;
    $dataSets[$iDataset]->{'yLowerError'} = -$yLowerError+$y;
    $dataSets[$iDataset]->{'yUpperError'} = +$yUpperError-$y;
    $dataSets[$iDataset]->{'label'}       = $dataSet->{'label'};

    # Compute a binned mean star formation rate.
    my $e      = sqrt($yUpperError**2+$yLowerError**2);
    my $weight = 1.0/($yUpperError**2+$yLowerError**2);
    (my $yBinned, my $yBinnedError, my $ySigma, my $ySigmaError)
	= &Means::BinnedMean($redshiftBins,$x,$y,$weight);
    (my $eBinned, my $eBinnedError, my $eSigma, my $eSigmaError)
	= &Means::BinnedMean($redshiftBins,$x,$e,$weight);
    my $sigmaMax = which($ySigma > $eBinned);
    $eBinned->index($sigmaMax) .= $ySigma->index($sigmaMax);
    my $empty = which($yBinnedError == 0.0);
    $eBinned->index($empty) .= 1.0e30;

    # Interpolate model to data points and compute chi^2.
    (my $sfrInterpolated, my $error) = interpolate($redshiftBins,$redshift,$SFR);
    $chiSquared += sum((($yBinned-$sfrInterpolated)/$eBinned)**2);
    $degreesOfFreedom += nelem($yBinned)-nelem($empty);
}

# Display chi^2 information if requested.
if ( $showFit == 1 ) {
    my %fitData;
    $fitData{'name'} = "Volume averaged star formation rate history.";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make plot of redshift evolution.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Star Formation Rate History'\n";
print $gnuPlot "set xlabel 'Redshift; \$z\$'\n";
print $gnuPlot "set ylabel 'Comoving star formation rate density; \$\\dot{\\rho}(z) [M_\\odot/\\hbox{yr}/\\hbox{Mpc}^{-3}]\$'\n";
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
print $gnuPlot "set xrange [-0.25:9.0]\n";
print $gnuPlot "set yrange [0.005:1.0]\n";
print $gnuPlot "set pointsize 2.0\n";
for($iDataset=0;$iDataset<scalar(@dataSets);++$iDataset) {
    &PrettyPlots::Prepare_Dataset(
	 \$plot,
	 $dataSets[$iDataset]->{'x'},$dataSets[$iDataset]->{'y'},
	 errorRight => $dataSets[$iDataset]->{'xUpperError'},
	 errorLeft  => $dataSets[$iDataset]->{'xLowerError'},
	 errorUp    => $dataSets[$iDataset]->{'yUpperError'},
	 errorDown  => $dataSets[$iDataset]->{'yLowerError'},
	 style      => "point",
	 symbol     => [6,7],
	 weight     => [5,3],
	 color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[$iDataset]},
	 title      => $dataSets[$iDataset]->{'label'}.' [observed]'
	);
}
my $nonZeroPoints = which($SFR > 0.0);
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $redshift->index($nonZeroPoints),$SFR->index($nonZeroPoints),
    style      => "line",
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{'redYellow'},
    title      => 'Galacticus'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);
&MetaData::Write($plotFile,$galacticusFile,$self);

exit;
