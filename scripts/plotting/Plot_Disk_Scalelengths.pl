#!/usr/bin/env perl
use lib "./perl";
use PDL;
use XML::Simple;
use Galacticus::HDF5;
use Galacticus::Magnitudes;
use Math::SigFigs;
use Stats::Histograms;
use Data::Dumper;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Disk_Scalelengths.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
    $outputFile = $outputTo."/Disk_Scalelengths.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['volumeWeight','diskScaleLength','magnitudeTotal:RGO_I:rest:z0.0000:dustAtlas[faceOn]:vega']);
$dataSets = $dataSet->{'dataSets'};
$scaleLength = $dataSets->{'diskScaleLength'};
$magnitude = $dataSets->{'magnitudeTotal:RGO_I:rest:z0.0000:dustAtlas[faceOn]:vega'};
$weight = $dataSets->{'volumeWeight'};
delete($dataSet->{'dataSets'});

# Open a pipe to GnuPlot.
open(gnuPlot,"|gnuplot");
print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
print gnuPlot "set output \"tmp.ps\"\n";

# Initialize chi^2 accumulator.
$chiSquared = 0.0;
$degreesOfFreedom = 0;

# Read the XML data file.
undef(@tmpFiles);
$xml = new XML::Simple;
$data = $xml->XMLin("data/Disk_Sizes_Dejong_2000.xml");
foreach $dataSet ( @{$data->{'sizeDistribution'}} ) {
    $columns = $dataSet->{'columns'};
    $x = pdl @{$columns->{'scaleLength'}->{'data'}};
    $y = pdl @{$columns->{'distribution'}->{'data'}};
    $yUpperLimit = pdl @{$columns->{'distributionUpperLimit'}->{'data'}};
    $yUpperError = pdl @{$columns->{'distributionErrorUp'}->{'data'}};
    $yLowerError = pdl @{$columns->{'distributionErrorDown'}->{'data'}};
    $x = $x*($columns->{'scaleLength'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
    $yUpperError = $y*(10.0**$yUpperError);
    $yLowerError = $y/(10.0**$yLowerError);
    $yUpperLimitArrow = pdl -0.3*$yUpperLimit;
    $magnitudeMinimum = $dataSet->{'magnitudeRange'}->{'minimum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});;
    $magnitudeMaximum = $dataSet->{'magnitudeRange'}->{'maximum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});;

    # Select Galacticus galaxies.
    $logScaleLengthSelected = where(3.0+log10($scaleLength),$magnitude > $magnitudeMinimum & $magnitude <= $magnitudeMaximum);
    $weightSelected         = where($weight                ,$magnitude > $magnitudeMinimum & $magnitude <= $magnitudeMaximum);
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

    # Make the plot.
    print gnuPlot "set xlabel \"r_{disk} [kpc]\"\n";
    print gnuPlot "set ylabel \"d^2n/dlog_{10}r_{disk}/dM_{I,0} [Mpc^{-3} mag^{-1}]\"\n";
    print gnuPlot "set title \"".$dataSet->{'description'}."\"\n";
    print gnuPlot "set logscale xy\n";
    print gnuPlot "set xrange [0.1:20]\n";
    print gnuPlot "set mxtics 2\n";
    print gnuPlot "set mytics 2\n";
    print gnuPlot "set pointsize 1.0\n";
    print gnuPlot "unset label\n";
    print gnuPlot "set label \"{/Symbol c}^2=".FormatSigFigs($chiSquaredRange,4)." [".$degreesOfFreedomRange."]\" at screen 0.6, screen 0.2\n";
    $plotCommand = "plot ";
    $join = "";
    if ( any($y > 0.0) ) {$plotCommand .= "'-' with errorbars lt 1 pt 6 title\"".$data->{'label'}."\""; $join = ", "};
    if ( any($y <= 0.0) ) {$plotCommand .= $join."'-' with vectors lt 1 title \"".$data->{'label'}." upper limit\""; $join = ", "};
    if ( any($yGalacticus > 0.0) ) {$plotCommand .= $join."'-' with errorbars lt 2 pt 6 title 'Galacticus'"};
    print gnuPlot $plotCommand."\n";
    if ( any($y > 0.0) ) {
	for($i=0;$i<nelem($x);++$i) {
	    if ( $y->index($i) > 0.0 ) {print gnuPlot $x->index($i)." ".$y->index($i)." ".$yLowerError->index($i)." ".$yUpperError->index($i)."\n"};
	}
	print gnuPlot "e\n";
    }
    if ( any($y <= 0.0) ) {
	for($i=0;$i<nelem($x);++$i) {
	    if ( $y->index($i) <= 0.0 ) {print gnuPlot $x->index($i)." ".$yUpperLimit->index($i)." 0.0 ".$yUpperLimitArrow->index($i)."\n"};
	}
	print gnuPlot "e\n";
    }
    if ( any($yGalacticus > 0.0) ) {
	for($i=0;$i<nelem($x);++$i) {
	    print gnuPlot $x->index($i)." ".$yGalacticus->index($i)." ".$errorDownGalacticus->index($i)." ".$errorUpGalacticus->index($i)."\n";
	}
	print gnuPlot "e\n";
    }
}

# Close the pipe to GnuPlot.
close(gnuPlot);

# Convert to PDF.
system("ps2pdf tmp.ps ".$outputFile);

# Clean up files.
unlink("tmp.ps",@tmpFiles);

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
