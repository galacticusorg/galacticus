#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
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
use File::Which;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require Stats::Histograms;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;

# Option to control whether the morphological classification should be coarse-grained into just three bins.
my $coarseGrain = 0;

# Get name of input and output files.
die("Plot_Morphological_Luminosity_Function.pl <galacticusFile> <outputDir/File> [<showFit>]")
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
my $outputDir;
my $outputFile;
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
    $outputDir = ".";
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/Morphological_Luminosity_Function.pdf";
    $outputDir = $outputTo;
}
(my $fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/$1/;

# Create data structure to read the results.
my $dataSet;
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);

# Read galaxy data.
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,[
		       'mergerTreeWeight',
		       'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega',
		       'bulgeToTotalLuminosities:2MASS_Ks:observed:z0.0000:dustAtlas'
		   ]);
my $dataSets  = $dataSet->{'dataSets'};

# Read the XML data file.
my $xml     = new XML::Simple;
my $data    = $xml->XMLin($galacticusPath."data/observations/luminosityFunctions/Morphological_Luminosity_Functions_2MASS_Devereux_2009.xml");

# Estimate bulge-to-total ratio ranges for each morphological class.
my $bulgeToTotal;
foreach my $morphology ( @{$data->{'morphology'}} ) {
	# Get the luminosity function.
	my $x  = pdl @{$morphology->{'magnitude'         }->{'datum'}};
	my $y  = pdl @{$morphology->{'luminosityFunction'}->{'datum'}};
	my $xMin    = where($x-5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'})-0.25                                                         ,$x == -23.25);
	my $xMax    = where($x-5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'})+0.25                                                         ,$x == -23.25);
	my $ySelect = where($y*($dataSet->{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'},$x == -23.25);
	$bulgeToTotal->{$morphology->{"class"}}->{"abundance"}              = $ySelect->index(0);
	$bulgeToTotal->{$morphology->{"class"}}->{"magnitude"}->{"minimum"} = $xMin   ->index(0);
	$bulgeToTotal->{$morphology->{"class"}}->{"magnitude"}->{"maximum"} = $xMax   ->index(0);
}

# Morphological classes (ordered).
my @classes = ("Sc-Scd","Sb-Sbc","Sa-Sab","Lenticular","Elliptical");

# Structure giving bulge-to-total ranges for each morphological class.
$bulgeToTotal->{"Total"}->{"minimum"} = 0.00;
$bulgeToTotal->{"Total"}->{"maximum"} = 1.00;

# Compute cumulative morphological fractions.
my $previousFraction = 0.0;
foreach my $class ( @classes ) {
    unless ( $class eq "Total" ) {
	$bulgeToTotal->{$class}->{"abundance"} /= $bulgeToTotal->{"Total"}->{"abundance"};
	$bulgeToTotal->{$class}->{"abundance"} += $previousFraction;
	$previousFraction = $bulgeToTotal->{$class}->{"abundance"};
    }
}

# Compute cumulative fraction of model galaxies by bulge-to-total ratio.
my $weight = where($dataSets->{'mergerTreeWeight'},
		($dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'} >= -23.50) &
		($dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'} <  -23.00)
		);
my $ratio  = where($dataSets->{'bulgeToTotalLuminosities:2MASS_Ks:observed:z0.0000:dustAtlas'},
		($dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'} >= -23.50) &
		($dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'} <  -23.00)
		);
die("Plot_Morphological_Luminosity_Function.pl: no galaxies found in normalization magnitude range")
    if ( nelem($ratio) == 0 );
my $indices          = $ratio->qsorti;
my $totalWeight      = $weight->sum;
my $orderedRatios    = $ratio->index($indices);
my $cumulativeWeight = ($weight->index($indices)->cumusumover)/$totalWeight;

# Interpolate to get corresponding model bulge-to-total ratios.
$previousFraction = 0.0;
foreach my $class ( @classes ) {
    (my $ratioBoundary, my $error) = interpolate($bulgeToTotal->{$class}->{"abundance"},$cumulativeWeight,$orderedRatios);
    $bulgeToTotal->{$class}->{"minimum"} = $previousFraction;
    $bulgeToTotal->{$class}->{"maximum"} = $ratioBoundary;
    $previousFraction = $ratioBoundary;
}

# Open a pipe to GnuPlot.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;

# Loop over all morphologies and make plots.
my $iPlot = 0;
my @leafFiles;
my @plotFiles;
my $chiSquared = 0;
my $degreesOfFreedom = 0;
foreach my $morphology ( @{$data->{'morphology'}} ) {

    next
	if ( $coarseGrain == 1 && ( $morphology->{'class'} eq "Lenticular" || $morphology->{'class'} eq "Sc-Scd" ) );

    if ( exists($bulgeToTotal->{$morphology->{'class'}}) ) {
	# Get the luminosity function.
	my $x     = pdl @{$morphology->{'magnitude'              }->{'datum'}};
	my $y     = pdl @{$morphology->{'luminosityFunction'     }->{'datum'}};
	my $error = pdl @{$morphology->{'luminosityFunctionError'}->{'datum'}};
	
	if ( $coarseGrain == 1 && $morphology->{'class'} eq "Elliptical" ) {
	    my $m1     = ${$data->{'morphology'}}[$iPlot+1];
	    my $y1     = pdl @{$m1->{'luminosityFunction'     }->{'datum'}};
	    my $error1 = pdl @{$m1->{'luminosityFunctionError'}->{'datum'}};
	    $y     .= $y+$y1;
	    $error .= sqrt($error**2+$error1**2);
        }
	if ( $coarseGrain == 1 && $morphology->{'class'} eq "Sb-Sbc" ) {
	    my $m1     = ${$data->{'morphology'}}[$iPlot+1];
	    my $y1     = pdl @{$m1->{'luminosityFunction'     }->{'datum'}};
	    my $error1 = pdl @{$m1->{'luminosityFunctionError'}->{'datum'}};
	    $y     .= $y+$y1;
	    $error .= sqrt($error**2+$error1**2);
        }
	
	# Scale for Hubble parameter.
	$x     += -5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'});
	$y     *= ($dataSet->{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'};
	$error *= ($dataSet->{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'};

	# Construct Galacticus luminosity function.
	my $bulgeToTotalMinimum = $bulgeToTotal->{$morphology->{"class"}}->{"minimum"};
	my $bulgeToTotalMaximum = $bulgeToTotal->{$morphology->{"class"}}->{"maximum"};
	if ( $coarseGrain == 1 && $morphology->{'class'} eq "Elliptical" ) {
	    $bulgeToTotalMinimum = $bulgeToTotal->{"Lenticular"}->{"minimum"};
	}
	if ( $coarseGrain == 1 && $morphology->{'class'} eq "Sb-Sbc" ) {
	    $bulgeToTotalMinimum = $bulgeToTotal->{"Sc-Scd"}->{"minimum"};
	}
	my $magnitude = where($dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'},
			     ($dataSets->{'bulgeToTotalLuminosities:2MASS_Ks:observed:z0.0000:dustAtlas'} >= $bulgeToTotalMinimum)
			   & ($dataSets->{'bulgeToTotalLuminosities:2MASS_Ks:observed:z0.0000:dustAtlas'} <= $bulgeToTotalMaximum)
	    );
	my $weight    = where($dataSets->{'mergerTreeWeight'},
			     ($dataSets->{'bulgeToTotalLuminosities:2MASS_Ks:observed:z0.0000:dustAtlas'} >= $bulgeToTotalMinimum)
			   & ($dataSets->{'bulgeToTotalLuminosities:2MASS_Ks:observed:z0.0000:dustAtlas'} <= $bulgeToTotalMaximum)
	    );

	(my $yGalacticus, my $errorGalacticus) = &Histograms::Histogram($x,$magnitude,$weight,differential => 1);
	
        # Compute chi^2.
	my $thisChiSquared        = sum(($yGalacticus-$y)**2/($errorGalacticus**2+$error**2));
	my $thisDegreesOfFreedom  = nelem($y);
	$chiSquared           += $thisChiSquared;
	$degreesOfFreedom     += $thisDegreesOfFreedom;

	# Make the plot.
	++$iPlot;
	(my $thisPlot = $plotFile) =~ s/\.pdf/$iPlot.pdf/;
	(my $thisPlotEPS = $thisPlot) =~ s/\.pdf$/.eps/;
	open(my $gnuPlot,"|gnuplot");# 1>/dev/null 2>&1");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$thisPlotEPS."'\n";
	print $gnuPlot "set xlabel '\$M_{\\rm K,vega}\$'\n";
	print $gnuPlot "set ylabel '\${\\rm d}n/{\\rm d}\log M_{\\rm K,vega} [\\hbox{Mpc}^{-3}]\$'\n";
	my $label = $morphology->{'class'};
	$label .= " + Lenticular" if ( $coarseGrain == 1 && $morphology->{'class'} eq "Elliptical" );
	$label = "Sb-Scd" if ( $coarseGrain == 1 && $morphology->{'class'} eq "Sb-Sbc" );
	print $gnuPlot "set title '".$label." K Luminosity Function'\n";
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
	print $gnuPlot "set xrange [-26.0:-19.0]\n";
	print $gnuPlot "set yrange [3.0e-6:3.0e-2]\n";
	print $gnuPlot "set pointsize 2.0\n";
	&PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $x,$y,
	    errorUp    => $error,
	    errorDown  => $error,
	    style      => "point",
	    symbol     => [6,7],
	    weight     => [5,3],
	    color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[0]},
	    title      => $data->{'label'}.' [observed]'
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

        # Close the pipe to GnuPlot.
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($thisPlotEPS);
	(my $leafName = $thisPlot) =~ s/^.*\/([^\/]+)$/$1/;
	push(@leafFiles,$leafName);
	push(@plotFiles,$thisPlot);
    }
}
die("Plot_Morphological_Luminosity_Function.pl: 'pdfmerge' tool is required")
    unless ( which("pdfmerge") );
&SystemRedirect::tofile("rm -f ".$outputFile."; cd ".$outputDir."; pdfmerge ".join(" ",@leafFiles)." tmp.pdf; cd -; mv ".$outputDir."/tmp.pdf ".$outputFile,"/dev/null");
&MetaData::Write($outputFile,$galacticusFile,$self);
unlink(@plotFiles);

# Output fit data.
if ( $showFit == 1 ) {
    my %fitData;
    $fitData{'name'} = "Devereuc et al. (2009) morphologically segregated K-band luminosity functions";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

exit;
