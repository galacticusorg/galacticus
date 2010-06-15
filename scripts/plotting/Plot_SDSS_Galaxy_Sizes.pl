#!/usr/bin/env perl
use lib "./perl";
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Galacticus::HDF5;
use Galacticus::Magnitudes;
use Math::SigFigs;
use Stats::Histograms;
use Data::Dumper;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_SDSS_Galaxy_Sizes.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
    $outputFile = $outputTo."/SDSS_Galaxy_Sizes.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Read tabulation of half-radii vs. disk/spheroid properties and extract to PDLs.
$xml                 = new XML::Simple;
$halfRadiiData       = $xml->XMLin("data/Half_Radii_Exponential_Hernquist.xml");
$spheroidRadii       = pdl @{$halfRadiiData->{'spheroidRadius'}->{'data'}};
$spheroidMasses      = pdl @{$halfRadiiData->{'spheroidMass'}->{'data'}};
$spheroidRadiiIndex  = pdl 0..nelem($spheroidRadii )-1;
$spheroidMassesIndex = pdl 0..nelem($spheroidMasses)-1;
$halfRadiiTable      = zeroes(nelem($spheroidRadii),nelem($spheroidMasses));
for($iRadius=0;$iRadius<nelem($spheroidRadii);++$iRadius) {
    $halfRadiiTable(($iRadius),:) .= pdl @{${$halfRadiiData->{'halfRadius'}}[$iRadius]->{'data'}};
}

# Create data structure to read the results.
$dataSet{'file'} = $galacticusFile;
$dataSet{'store'} = 0;
&HDF5::Get_Parameters(\%dataSet);
&HDF5::Count_Trees(\%dataSet);
&HDF5::Select_Output(\%dataSet,0.1);
$dataSet{'tree'} = "all";
&HDF5::Get_Dataset(\%dataSet,['volumeWeight'
			      ,'diskStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'
			      ,'spheroidStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'
			      ,'magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB'
			      ,'diskScaleLength'
			      ,'spheroidScaleLength'
		   ]);
$dataSets = \%{$dataSet{'dataSets'}};
$spheroidScaleLength = ${$dataSets->{'spheroidScaleLength'}}/${$dataSets->{'diskScaleLength'}};
$spheroidLuminosity  = ${$dataSets->{'spheroidStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'}}/${$dataSets->{'diskStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'}};
$indexScaleLength = interpol($spheroidScaleLength,$spheroidRadii,$spheroidRadiiIndex);
$indexLuminosity  = interpol($spheroidLuminosity ,$spheroidMasses,$spheroidMassesIndex);
$radius           = $halfRadiiTable->interpND(transpose(cat($indexScaleLength,$indexLuminosity)));
$radius          *= 1000.0*${$dataSets->{'diskScaleLength'}};
$morphology       = ${$dataSets->{'spheroidStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'}}/(${$dataSets->{'diskStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'}}+${$dataSets->{'spheroidStellarLuminosity:SDSS_r:observed:z0.1000:dustAtlas'}});
$magnitude        = ${$dataSets->{'magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB'}};
$weight           = ${$dataSets->{'volumeWeight'}};

# Initialize chi^2 accumulator.
$chiSquared = 0.0;
$degreesOfFreedom = 0;

# Open a pipe to GnuPlot.
open(gnuPlot,"|gnuplot");
print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
print gnuPlot "set output \"tmp.ps\"\n";

# Read the XML data file.
undef(@tmpFiles);
$xml = new XML::Simple;
$data = $xml->XMLin("data/SDSS_Galaxy_Sizes_Shen_2003.xml");
foreach $dataSet ( @{$data->{'sizeDistribution'}} ) {
    $columns = $dataSet->{'columns'};
    $x = pdl @{$columns->{'radius'}->{'data'}};
    $y = pdl @{$columns->{'distribution'}->{'data'}};
    $yError = pdl @{$columns->{'distributionError'}->{'data'}};
    $x = (10.0**$x)*($columns->{'radius'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
    $yUpperError = $y+$yError;
    $yLowerError = $y-$yError;

    # Select Galacticus galaxies and compute distribution.
    $logRadiusSelected = where(log10($radius),
			    $magnitude >= $dataSet->{'magnitudeRange'}->{'minimum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'})
			    & $magnitude < $dataSet->{'magnitudeRange'}->{'maximum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'})
			    & $morphology >= $dataSet->{'morphologyRange'}->{'minimum'}
			    & $morphology < $dataSet->{'morphologyRange'}->{'maximum'}
	);
    $weightSelected = where($weight,
			    $magnitude >= $dataSet->{'magnitudeRange'}->{'minimum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'})
			    & $magnitude < $dataSet->{'magnitudeRange'}->{'maximum'}-5.0*log10($dataSet->{'magnitudeRange'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'})
			    & $morphology >= $dataSet->{'morphologyRange'}->{'minimum'}
			    & $morphology < $dataSet->{'morphologyRange'}->{'maximum'}
	);
    
    if ( nelem($logRadiusSelected) > 0 ) {
	$xBins = log10($x);
	($yGalacticus,$errorGalacticus) = &Histograms::Histogram($xBins,$logRadiusSelected,$weightSelected
								 ,normalized => 1, differential => 1);
	# Compute chi^2.
	$chiSquaredList = where(($yGalacticus-$y)**2/($errorGalacticus**2+$yError**2),$y > 0.0 & $errorGalacticus > 0.0);
	$chiSquaredRange = sum($chiSquaredList);
	$degreesOfFreedomRange = nelem($chiSquaredList);
	$chiSquared += $chiSquaredRange;
	$degreesOfFreedom += $degreesOfFreedomRange;
    } else {
	$yGalacticus = zeroes(nelem($xBins));
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
    for($i=0;$i<nelem($x);++$i) {
	if ( $y->index($i) > 0.0 ) {print gnuPlot $x->index($i)." ".$y->index($i)." ".$yLowerError->index($i)." ".$yUpperError->index($i)."\n"};
    }
    print gnuPlot "e\n";
    if ( any($yGalacticus > 0.0)  ) {
	for($i=0;$i<nelem($x);++$i) {
	    if ( $yGalacticus->index($i) > 0.0 ) {print gnuPlot $x->index($i)." ".$yGalacticus->index($i)." ".$errorGalacticus->index($i)."\n"};
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
    $fitData{'name'} = "Shen et al. (2003) galaxy half-light radius distributions";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

exit;
