#!/usr/bin/env perl
use lib "./perl";
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Galacticus::HDF5;
use Galacticus::Magnitudes;
use Math::SigFigs;
use Data::Dumper;
use Stats::Histograms;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Morphological_Luminosity_Function.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
    $outputFile = $outputTo."/Morphological_Luminosity_Function.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataSet{'file'} = $galacticusFile;
$dataSet{'store'} = 0;
&HDF5::Get_Parameters(\%dataSet);
&HDF5::Count_Trees(\%dataSet);
&HDF5::Select_Output(\%dataSet,0.0);

# Read galaxy data.
$dataSet{'tree'} = "all";
&HDF5::Get_Dataset(\%dataSet,[
		       'volumeWeight',
		       'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega',
		       'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'
		   ]);
$dataSets  = \%{$dataSet{'dataSets'}};

# Read the XML data file.
$xml     = new XML::Simple;
$data    = $xml->XMLin("data/Morphological_Luminosity_Functions_2MASS_Devereux_2009.xml");

# Estimate bulge-to-total ratio ranges for each morphological class.
foreach $morphology ( @{$data->{'morphology'}} ) {
	# Get the luminosity function.
	$x  = pdl @{$morphology->{'magnitude'         }->{'datum'}};
	$y  = pdl @{$morphology->{'luminosityFunction'}->{'datum'}};
	$xMin    = where($x-5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'})-0.25                                                         ,$x == -23.25);
	$xMax    = where($x-5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'})+0.25                                                         ,$x == -23.25);
	$ySelect = where($y*($dataSet{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'},$x == -23.25);
	$bulgeToTotal->{$morphology->{"class"}}->{"abundance"}              = $ySelect->index(0);
	$bulgeToTotal->{$morphology->{"class"}}->{"magnitude"}->{"minimum"} = $xMin   ->index(0);
	$bulgeToTotal->{$morphology->{"class"}}->{"magnitude"}->{"maximum"} = $xMax   ->index(0);
}

# Morphological classes (ordered).
@classes = ("Sc-Scd","Sb-Sbc","Sa-Sab","Lenticular","Elliptical");

# Structure giving bulge-to-total ranges for each morphological class.
$bulgeToTotal->{"Total"}->{"minimum"} = 0.00;
$bulgeToTotal->{"Total"}->{"maximum"} = 1.00;

# Compute cumulative morphological fractions.
$previousFraction = 0.0;
foreach $class ( @classes ) {
    unless ( $class eq "Total" ) {
	$bulgeToTotal->{$class}->{"abundance"} /= $bulgeToTotal->{"Total"}->{"abundance"};
	$bulgeToTotal->{$class}->{"abundance"} += $previousFraction;
	$previousFraction = $bulgeToTotal->{$class}->{"abundance"};
    }
}

# Compute cumulative fraction of model galaxies by bulge-to-total ratio.
$weight = where(${$dataSets->{'volumeWeight'}},
		${$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'}} >= $bulgeToTotal->{"Total"}->{"magnitude"}->{"minimum"} &
		${$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'}} < $bulgeToTotal->{"Total"}->{"magnitude"}->{"maximum"}
		);
$ratio  = where(${$dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'}},
		${$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'}} >= $bulgeToTotal->{"Total"}->{"magnitude"}->{"minimum"} &
		${$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'}} < $bulgeToTotal->{"Total"}->{"magnitude"}->{"maximum"}
		);

$indices          = $ratio->qsorti;
$totalWeight      = $weight->sum;
$orderedRatios    = $ratio->index($indices);
$cumulativeWeight = ($weight->index($indices)->cumusumover)/$totalWeight;

# Interpolate to get corresponding model bulge-to-total ratios.
$previousFraction = 0.0;
foreach $class ( @classes ) {
    ($ratioBoundary,$error) = interpolate($bulgeToTotal->{$class}->{"abundance"},$cumulativeWeight,$orderedRatios);
    $bulgeToTotal->{$class}->{"minimum"} = $previousFraction;
    $bulgeToTotal->{$class}->{"maximum"} = $ratioBoundary;
    $previousFraction = $ratioBoundary;
}

# Open a pipe to GnuPlot.
open(gnuPlot,"|gnuplot");
print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
print gnuPlot "set output \"tmp.ps\"\n";

# Loop over all morphologies and make plots.
foreach $morphology ( @{$data->{'morphology'}} ) {
    if ( exists($bulgeToTotal->{$morphology->{'class'}}) ) {
	# Get the luminosity function.
	$x     = pdl @{$morphology->{'magnitude'              }->{'datum'}};
	$y     = pdl @{$morphology->{'luminosityFunction'     }->{'datum'}};
	$error = pdl @{$morphology->{'luminosityFunctionError'}->{'datum'}};
	
	# Scale for Hubble parameter.
	$x     += -5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
	$y     *= ($dataSet{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'};
	$error *= ($dataSet{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'};

	# Construct Galacticus luminosity function.
	$magnitude = where(${$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'}},
			     ${$dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'}} >= $bulgeToTotal->{$morphology->{"class"}}->{"minimum"}
			   & ${$dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'}} <= $bulgeToTotal->{$morphology->{"class"}}->{"maximum"}
	    );
	$weight    = where(${$dataSets->{'volumeWeight'}},
			     ${$dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'}} >= $bulgeToTotal->{$morphology->{"class"}}->{"minimum"}
			   & ${$dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'}} <= $bulgeToTotal->{$morphology->{"class"}}->{"maximum"}
	    );

	($yGalacticus,$errorGalacticus) = &Histograms::Histogram($x,$magnitude,$weight,differential => 1);
	
        # Compute chi^2.
	$thisChiSquared        = sum(($yGalacticus-$y)**2/($errorGalacticus**2+$error**2));
	$thisDegreesOfFreedom  = nelem($y);
	$chiSquared           += $thisChiSquared;
	$degreesOfFreedom     += $thisDegreesOfFreedom;

	# Make the plot.
	print gnuPlot "set xlabel \"M_{K,vega}\"\n";
	print gnuPlot "set ylabel \"dn/dlogM_{K,vega} [Mpc^{-3}]\"\n";
	print gnuPlot "set title '".$morphology->{'class'}." K Luminosity Function'\n";
	print gnuPlot "unset label\n";
	print gnuPlot "set label \"{/Symbol c}^2=".FormatSigFigs($thisChiSquared,4)." [".$thisDegreesOfFreedom."]\" at screen 0.6, screen 0.2\n";
	print gnuPlot "set key left\n";
	print gnuPlot "set logscale y\n";
	print gnuPlot "set mxtics 2\n";
	print gnuPlot "set mytics 10\n";
	print gnuPlot "set format y '10^{\%L}'\n";
	print gnuPlot "set pointsize 1.0\n";
	print gnuPlot "plot '-' with errorbars lt 1 pt 6 title '".$data->{'label'}."', '-' with errorbars lt 2 pt 4 title \"Galacticus\"\n";
	for ($i=0;$i<nelem($x);++$i) {
	    print gnuPlot $x->index($i)." ".$y->index($i)." ".$error->index($i)."\n";
	}
	print gnuPlot "e\n";
	for ($i=0;$i<nelem($x);++$i) {
	    print gnuPlot $x->index($i)." ".$yGalacticus->index($i)." ".$errorGalacticus->index($i)."\n";
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

# Output fit data.
if ( $showFit == 1 ) {
    $fitData{'name'} = "Cole et al. (2001) K-band luminosity function";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

exit;
