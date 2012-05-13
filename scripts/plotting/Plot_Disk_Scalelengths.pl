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
use LaTeX::Encode;
require Stats::Histograms;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require System::Redirect;
require XMP::MetaData;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Disk_Scalelengths.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
my $outputDir;
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
    $outputDir = ".";
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/Disk_Scalelengths.pdf";
    $outputDir = $outputTo;
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataBlock->{'file'} = $galacticusFile;
$dataBlock->{'store'} = 0;
&HDF5::Get_Parameters($dataBlock);
&HDF5::Count_Trees($dataBlock);
&HDF5::Select_Output($dataBlock,0.0);
$dataBlock->{'tree'} = "all";
&HDF5::Get_Dataset($dataBlock,['volumeWeight','diskScaleLength','magnitudeTotal:RGO_I:rest:z0.0000:dustAtlas[faceOn]:vega','bulgeToTotalLuminosity:RGO_I:rest:z0.0000:dustAtlas']);
$dataSets = $dataBlock->{'dataSets'};
$scaleLength = $dataSets->{'diskScaleLength'};
$magnitude = $dataSets->{'magnitudeTotal:RGO_I:rest:z0.0000:dustAtlas[faceOn]:vega'};
$morphology = $dataSets->{'bulgeToTotalLuminosity:RGO_I:rest:z0.0000:dustAtlas'};
$weight = $dataSets->{'volumeWeight'};
delete($dataBlock->{'dataSets'});
  
# Initialize chi^2 accumulator.
$chiSquared = 0.0;
$degreesOfFreedom = 0;

# Read the XML data file.
undef(@tmpFiles);
$xml = new XML::Simple;
$data = $xml->XMLin($galacticusPath."data/observations/galaxySizes/Disk_Sizes_Dejong_2000.xml");
my $i = -1;
my @leafFiles;
my @plotFiles;
foreach $dataSet ( @{$data->{'sizeDistribution'}} ) {
    $columns = $dataSet->{'columns'};
    $x = pdl @{$columns->{'scaleLength'}->{'data'}};
    $y = pdl @{$columns->{'distribution'}->{'data'}};
    $yUpperLimit = pdl @{$columns->{'distributionUpperLimit'}->{'data'}};
    $yUpperError = pdl @{$columns->{'distributionErrorUp'}->{'data'}};
    $yLowerError = pdl @{$columns->{'distributionErrorDown'}->{'data'}};
    $x = $x*($columns->{'scaleLength'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'});
    $yUpperError = +$y*(10.0**$yUpperError)-$y;
    $yLowerError = -$y/(10.0**$yLowerError)+$y;

    my $zeroPoints = which($y <= 0.0);
    my $yP           = $y          ->copy();
    my $yUpperErrorP = $yUpperError->copy();
    my $yLowerErrorP = $yLowerError->copy();
    $yP          ->index($zeroPoints) .=      $yUpperLimit->index($zeroPoints);
    $yUpperErrorP->index($zeroPoints) .=  0.0                                 ;
    $yLowerErrorP->index($zeroPoints) .= -0.7*$yUpperLimit->index($zeroPoints);
    $magnitudeMinimum = $dataSet->{'magnitudeRange'}->{'minimum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'});;
    $magnitudeMaximum = $dataSet->{'magnitudeRange'}->{'maximum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'});;
    $morphologyMinimum = 0.05; # Morpology constraints are approximate, based on De Jong & Lacey's T-type selection of 3 < T < 8.
    $morphologyMaximum = 0.30;

    # Select Galacticus galaxies.
    $logScaleLengthSelected = where(3.0+log10($scaleLength),$magnitude > $magnitudeMinimum & $magnitude <= $magnitudeMaximum & $morphology >= $morphologyMinimum & $morphology <= $morphologyMaximum);
    $weightSelected         = where($weight                ,$magnitude > $magnitudeMinimum & $magnitude <= $magnitudeMaximum & $morphology >= $morphologyMinimum & $morphology <= $morphologyMaximum);
    $xBins = log10($x);
    ($yGalacticus,$errorGalacticus) = &Histograms::Histogram($xBins,$logScaleLengthSelected,$weightSelected
							     ,differential => 1);
    $yGalacticus         /= ($magnitudeMaximum-$magnitudeMinimum);
    $errorGalacticus     /= ($magnitudeMaximum-$magnitudeMinimum);
    $errorUpGalacticus    = $yGalacticus+$errorGalacticus;
    $errorDownGalacticus  = $yGalacticus-$errorGalacticus;

    # Compute chi^2.
    $chiSquaredList = where(($yGalacticus-$y)**2/($errorGalacticus**2+(0.5*($yUpperError-$yLowerError))**2),$y > 0.0 & $errorGalacticus > 0.0);
    $chiSquaredRange = sum($chiSquaredList);
    $degreesOfFreedomRange = nelem($chiSquaredList);
    $chiSquared += $chiSquaredRange;
    $degreesOfFreedom += $degreesOfFreedomRange;

    # Make plot of disk size distribution.
    ++$i;
    my $plot;
    my $gnuPlot;
    (my $plotFile = $outputFile) =~ s/\.pdf/_$i.pdf/;
    (my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
    open($gnuPlot,"|gnuplot > /dev/null 2>&1");
    print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
    print $gnuPlot "set output '".$plotFileEPS."'\n";
    my $title = $dataSet->{'description'};
    $title =~ s/([\-\d]+\s*<\s*)M(\s*<\s*[\-\d]+)/\$$1M_{\\rm I,0}$2\$/;
    print $gnuPlot "set title '".$title."'\n";
    print $gnuPlot "set xlabel 'Disk scale length; \$r_{\\rm disk}\$ [kpc]'\n";
    print $gnuPlot "set ylabel 'Comoving number density; \${\\rm d}^2n/{\\rm d}\\log_{10}r_{\\rm disk}/{\\rm d}M_{\\rm I,0} [\\hbox{Mpc}^{-3} \\hbox{mag}^{-1}]\$'\n";
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
    print $gnuPlot "set xrange [0.2:30]\n";
    print $gnuPlot "set yrange [1.0e-5:1.0e0]\n";
    print $gnuPlot "set pointsize 2.0\n";
    my $label = latex_encode($data->{'label'});
    $label =~ s/\\/\\\\/g;
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x,$yP,
	errorDown  => $yLowerErrorP,
	errorUp    => $yUpperErrorP,
	style      => "point",
	symbol     => [6,7],
	weight     => [5,3],
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[0]},
	title      => $label.' [observed]'
	);
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x,$yGalacticus,
	errorDown  => $errorDownGalacticus,
	errorUp    => $errorUpGalacticus,
	style      => "point",
	symbol     => [6,7],
	weight     => [5,3],
	color      => $PrettyPlots::colorPairs{'redYellow'},
	title      => 'Galacticus'
	);
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileEPS);
    (my $leafName = $plotFile) =~ s/^.*\/([^\/]+)$/$1/;
    push(@leafFiles,$leafName);
    push(@plotFiles,$plotFile);
}
&SystemRedirect::tofile("rm -f ".$outputFile."; cd ".$outputDir."; pdfmerge ".join(" ",@leafFiles)." tmp.pdf; cd -; mv ".$outputDir."/tmp.pdf ".$outputFile,"/dev/null");
unlink(@plotFiles);
&MetaData::Write($outputFile,$galacticusFile,$self);

# Display chi^2 information
if ( $showFit == 1 ) {
    $fitData{'name'} = "Dejong & Lacey (2000) disk scale length distributions";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

exit;
