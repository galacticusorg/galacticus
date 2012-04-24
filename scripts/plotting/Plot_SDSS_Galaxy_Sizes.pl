#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
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
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require Stats::Histograms;
require XMP::MetaData;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_SDSS_Galaxy_Sizes.pl <galacticusFile> <outputDir/File> [<showFit>]")};
my $self           = $0;
my $galacticusFile = $ARGV[0];
my $outputTo       = $ARGV[1];
my $showFit;
if ( $#ARGV == 2 ) {
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
    $outputFile = $outputTo."/SDSS_Galaxy_Sizes.pdf";
}
(my $fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/$1/;

# Read tabulation of half-radii vs. disk/spheroid properties and extract to PDLs.
my $xml                 = new XML::Simple;
my $halfRadiiData       = $xml->XMLin($galacticusPath."data/Half_Radii_Exponential_Hernquist.xml");
my $spheroidRadii       = pdl @{$halfRadiiData->{'spheroidRadius'}->{'data'}};
my $spheroidMasses      = pdl @{$halfRadiiData->{'spheroidMass'}->{'data'}};
my $spheroidRadiiIndex  = pdl 0..nelem($spheroidRadii )-1;
my $spheroidMassesIndex = pdl 0..nelem($spheroidMasses)-1;
my $halfRadiiTable      = zeroes(nelem($spheroidRadii),nelem($spheroidMasses));
for(my $iRadius=0;$iRadius<nelem($spheroidRadii);++$iRadius) {
    $halfRadiiTable(($iRadius),:) .= pdl @{${$halfRadiiData->{'halfRadius'}}[$iRadius]->{'data'}};
}

# Create data structure to read the results.
my $dataBlock;
$dataBlock->{'file'} = $galacticusFile;
$dataBlock->{'store'} = 0;
&HDF5::Get_Parameters($dataBlock);
&HDF5::Count_Trees($dataBlock);
&HDF5::Select_Output($dataBlock,0.1);
$dataBlock->{'tree'} = "all";
&HDF5::Get_Dataset($dataBlock,['volumeWeight'
			      ,'diskStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'
			      ,'spheroidStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'
			      ,'magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB'
			      ,'diskScaleLength'
			      ,'spheroidScaleLength'
		   ]);
my $dataSets = $dataBlock->{'dataSets'};
my $spheroidScaleLength = $dataSets->{'spheroidScaleLength'}/$dataSets->{'diskScaleLength'};
my $spheroidLuminosity  = $dataSets->{'spheroidStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'}/$dataSets->{'diskStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'};
my $indexScaleLength = interpol($spheroidScaleLength,$spheroidRadii,$spheroidRadiiIndex);
my $indexLuminosity  = interpol($spheroidLuminosity ,$spheroidMasses,$spheroidMassesIndex);
my $radius           = $halfRadiiTable->interpND(transpose(cat($indexScaleLength,$indexLuminosity)));
$radius          *= 1000.0*$dataSets->{'diskScaleLength'};
my $morphology       = $dataSets->{'spheroidStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'}/($dataSets->{'diskStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'}+$dataSets->{'spheroidStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'});
my $magnitude        = $dataSets->{'magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB'};
my $weight           = $dataSets->{'volumeWeight'};

# Initialize chi^2 accumulator.
my $chiSquared = 0.0;
my $degreesOfFreedom = 0;

# Open a pipe to GnuPlot.
open(gnuPlot,"|gnuplot > /dev/null 2&>1");
print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
print gnuPlot "set output \"tmp.ps\"\n";

# Read the XML data file.
my @tmpFiles;
undef(@tmpFiles);
my $data = $xml->XMLin($galacticusPath."data/SDSS_Galaxy_Sizes_Shen_2003.xml");
foreach my $dataSet ( @{$data->{'sizeDistribution'}} ) {
    my $columns = $dataSet->{'columns'};
    my $x = pdl @{$columns->{'radius'}->{'data'}};
    my $y = pdl @{$columns->{'distribution'}->{'data'}};
    my $yError = pdl @{$columns->{'distributionError'}->{'data'}};
    $x = (10.0**$x)*($columns->{'radius'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'});
    my $yUpperError = $y+$yError;
    my $yLowerError = $y-$yError;

    # Select Galacticus galaxies and compute distribution.
    my $logRadiusSelected = where(log10($radius),
			    ($magnitude >= $dataSet->{'magnitudeRange'}->{'minimum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'}))
			    & ($magnitude < $dataSet->{'magnitudeRange'}->{'maximum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'}))
			    & ($morphology >= $dataSet->{'morphologyRange'}->{'minimum'})
			    & ($morphology < $dataSet->{'morphologyRange'}->{'maximum'})
	);
    my $weightSelected = where($weight,
			    ($magnitude >= $dataSet->{'magnitudeRange'}->{'minimum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'}))
			    & ($magnitude < $dataSet->{'magnitudeRange'}->{'maximum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'}))
			    & ($morphology >= $dataSet->{'morphologyRange'}->{'minimum'})
			    & ($morphology < $dataSet->{'morphologyRange'}->{'maximum'})
	);
    
    my $xGalacticus;
    my $yGalacticus;
    my $errorGalacticus;
    my $chiSquaredRange;
    my $degreesOfFreedomRange;
    if ( nelem($logRadiusSelected) > 0 ) {
	# Compute histogram first for chi^2 purposes.
	my $xBins = log10($x);
	($yGalacticus,$errorGalacticus) = &Histograms::Histogram($xBins,$logRadiusSelected,$weightSelected
								 ,normalized => 1, differential => 1);
	# Compute chi^2.
	my $chiSquaredList = where(($yGalacticus-$y)**2/($errorGalacticus**2+$yError**2),($y > 0.0) & ($errorGalacticus > 0.0));
	$chiSquaredRange = sum($chiSquaredList);
	$degreesOfFreedomRange = nelem($chiSquaredList);
	$chiSquared += $chiSquaredRange;
	$degreesOfFreedom += $degreesOfFreedomRange;
	# Now recompute for plotting.
	my $delta = $xBins->index(1)-$xBins->index(0);
	my $xMin = $xBins->index(0)-(int(($xBins->index(0)+0.5)/$delta)+1)*$delta;
	my $xMax = $xBins->index(0)+(int((1.5-$xBins->index(0))/$delta)+1)*$delta;
	my $nBin = int(($xMax-$xMin)/$delta)+1;
	my $xBinsGalacticus = $xMin+$delta*sequence($nBin);
	$xGalacticus = 10.0**$xBinsGalacticus;
	($yGalacticus,$errorGalacticus) = &Histograms::Histogram($xBinsGalacticus,$logRadiusSelected,$weightSelected
								 ,normalized => 1, differential => 1);
    } else {
	$xGalacticus     = pdl zeroes(nelem($x));
	$yGalacticus     = pdl zeroes(nelem($x));
	$errorGalacticus = pdl zeroes(nelem($x));
    }
    
    # Make the plot.
    print gnuPlot "set xlabel \"r_{disk} [kpc]\"\n";
    print gnuPlot "set ylabel \"dn/dlog_{10}r_{disk}/n\"\n";
    print gnuPlot "set title \"".$dataSet->{'description'}."\"\n";
    print gnuPlot "set logscale xy\n";
    print gnuPlot "set xrange [0.3:30]\n";
    print gnuPlot "set mxtics 2\n";
    print gnuPlot "set mytics 2\n";
    print gnuPlot "set pointsize 1.0\n";
    print gnuPlot "unset label\n";
    if ( nelem($logRadiusSelected) > 0 ) {print gnuPlot "set label \"{/Symbol c}^2=".FormatSigFigs($chiSquaredRange,4)." [".$degreesOfFreedomRange."]\" at screen 0.6, screen 0.2\n"};
    print gnuPlot "plot '-' with errorbars lt 1 pt 6 title \"".$data->{'label'}."\"";
    if ( any($yGalacticus > 0.0) ) {
	print gnuPlot ", '-' with errorbars lt 2 pt 6 title \"Galacticus\"";
    }
    print gnuPlot "\n";
    for(my $i=0;$i<nelem($x);++$i) {
	if ( $y->index($i) > 0.0 ) {print gnuPlot $x->index($i)." ".$y->index($i)." ".$yLowerError->index($i)." ".$yUpperError->index($i)."\n"};
    }
    print gnuPlot "e\n";
    if ( any($yGalacticus > 0.0)  ) {
	for(my $i=0;$i<nelem($x);++$i) {
	    if ( $yGalacticus->index($i) > 0.0 ) {print gnuPlot $xGalacticus->index($i)." ".$yGalacticus->index($i)." ".$errorGalacticus->index($i)."\n"};
	}
	print gnuPlot "e\n";
    }
}

# Close the pipe to GnuPlot.
close(gnuPlot);

# Convert to PDF.
system("ps2pdf tmp.ps ".$outputFile);
&MetaData::Write($outputFile,$galacticusFile,$self);

# Clean up files.
unlink("tmp.ps",@tmpFiles);

# Display chi^2 information
if ( $showFit == 1 ) {
    my %fitData;
    $fitData{'name'} = "Shen et al. (2003) galaxy half-light radius distributions";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

exit;
