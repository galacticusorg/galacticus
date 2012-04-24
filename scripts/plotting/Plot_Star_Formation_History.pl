#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
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
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Star_Formation_History.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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

# Get parameters for the Galacticus model.
$dataBlock->{'file'} = $galacticusFile;
&HDF5::Get_Parameters($dataBlock);

# Check if output location is file or directory.
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/Star_Formation_History.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Extract global data
&HDF5::Get_History($dataBlock,['historyExpansion','historyStarFormationRate']);
$history        = $dataBlock->{'history'};
$time           = $history->{'historyTime'};
$redshift       = 1.0/$history->{'historyExpansion'}-1.0;
$SFR            = $history->{'historyStarFormationRate'}/1.0e9;
$stellarDensity = $history->{'historyStellarDensity'};

# Determine IMF correction factor.
$imfCorrection = 1.0;
if ( $dataBlock->{'parameters'}->{'imfSelectionFixed'} eq "Chabrier" ) {
    $imfCorrection = 0.6;
}

# Define redshift bins.
$redshiftPoints = pdl 16;
$redshiftMin    = pdl 0.0;
$redshiftMax    = pdl 8.0;
$redshiftBin    = pdl ($redshiftMax-$redshiftMin)/$redshiftPoints;
$redshiftBins   = pdl (0..$redshiftPoints-1)*$redshiftBin+$redshiftMin+0.5*$redshiftBin;

# Read the XML data file.
$xml = new XML::Simple;
$data = $xml->XMLin($galacticusPath."data/Star_Formation_Rate_Data.xml");
$iDataset = -1;
$chiSquared = 0.0;
$degreesOfFreedom = 0;
foreach $dataSet ( @{$data->{'starFormationRate'}} ) {
    $columns = $dataSet->{'columns'};
    $x = pdl @{$columns->{'redshift'}->{'data'}};
    $xLowerError = pdl @{$columns->{'redshiftErrorDown'}->{'data'}};
    $xUpperError = pdl @{$columns->{'redshiftErrorUp'}->{'data'}};
    $xLowerError = $x-$xLowerError;
    $xUpperError = $x+$xUpperError;
    $y = pdl @{$columns->{'sfr'}->{'data'}};
    $yLowerError = pdl @{$columns->{'sfrErrorDown'}->{'data'}};
    $yUpperError = pdl @{$columns->{'sfrErrorUp'}->{'data'}};
    $yUpperError = (10.0**($y+$yUpperError));
    $yLowerError = (10.0**($y-$yLowerError));
    $y = (10.0**$y);

    # Compute cosmology corrections.
    $cosmologyData = Astro::Cosmology->new(
	omega_matter => $columns->{'sfr'}->{'omega'},
	omega_lambda => $columns->{'sfr'}->{'lambda'},
	H0           => $columns->{'sfr'}->{'hubble'}
	);
    $cosmologyGalacticus = Astro::Cosmology->new(
	omega_matter => $dataBlock->{'parameters'}->{'Omega_Matter'},
	omega_lambda => $dataBlock->{'parameters'}->{'Omega_DE'},
	H0           => $dataBlock->{'parameters'}->{'H_0'}
	);

    $volumeElementData            = $cosmologyData      ->differential_comoving_volume($x);
    $volumeElementGalacticus      = $cosmologyGalacticus->differential_comoving_volume($x);
    $luminosityDistanceData       = $cosmologyData      ->luminosity_distance         ($x);
    $luminosityDistanceGalacticus = $cosmologyGalacticus->luminosity_distance         ($x);
    $cosmologyCorrection = ($volumeElementData/$volumeElementGalacticus)*($luminosityDistanceGalacticus/$luminosityDistanceData)**2;
    
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
    $e      = sqrt($yUpperError**2+$yLowerError**2);
    $weight = 1.0/($yUpperError**2+$yLowerError**2);
    ($yBinned,$yBinnedError,$ySigma,$ySigmaError)
	= &Means::BinnedMean($redshiftBins,$x,$y,$weight);
    ($eBinned,$eBinnedError,$eSigma,$eSigmaError)
	= &Means::BinnedMean($redshiftBins,$x,$e,$weight);
    $sigmaMax = which($ySigma > $eBinned);
    $eBinned->index($sigmaMax) .= $ySigma->index($sigmaMax);
    $empty = which($yBinnedError == 0.0);
    $eBinned->index($empty) .= 1.0e30;

    # Interpolate model to data points and compute chi^2.
    ($sfrInterpolated,$error) = interpolate($redshiftBins,$redshift,$SFR);
    $chiSquared += sum((($yBinned-$sfrInterpolated)/$eBinned)**2);
    $degreesOfFreedom += nelem($yBinned)-nelem($empty);
}

# Display chi^2 information if requested.
if ( $showFit == 1 ) {
    $fitData{'name'} = "Volume averaged star formation rate history.";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make plot of redshift evolution.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot > /dev/null 2&>1");
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
for($iDataset=0;$iDataset<=$#dataSets;++$iDataset) {
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
