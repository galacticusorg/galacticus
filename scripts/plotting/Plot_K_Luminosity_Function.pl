#!/usr/bin/env perl
use lib "./perl";
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Graphics::GnuplotIF;
use Galacticus::HDF5;
use Galacticus::Magnitudes;
use Math::SigFigs;
use Data::Dumper;
use Stats::Histograms;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_K_Luminosity_Function.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
    $outputFile = $outputTo."/K_Luminosity_Function.pdf";
}

# Create data structure to read the results.
$dataSet{'file'} = $galacticusFile;
$dataSet{'store'} = 0;
&HDF5::Get_Parameters(\%dataSet);
&HDF5::Count_Trees(\%dataSet);
&HDF5::Select_Output(\%dataSet,0.0);

# Read the XML data file.
$xml     = new XML::Simple;
$data    = $xml->XMLin("data/K_Luminosity_Function_Cole_2001.xml");
$columns = $data->{'luminosityFunction'}->{'columns'};
$xBins   = pdl @{$columns->{'magnitude'}->{'data'}};
$x       = pdl @{$columns->{'magnitude'}->{'data'}};
$y       = pdl @{$columns->{'luminosityFunction'}->{'data'}};
$error   = pdl @{$columns->{'error'}->{'data'}};
$xBins   = $xBins-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
$x       = $x-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
$y       = $y*($dataSet{'parameters'}->{'H_0'}/$columns->{'luminosityFunction'}->{'hubble'})**3;
$error   = $error*($dataSet{'parameters'}->{'H_0'}/$columns->{'luminosityFunction'}->{'hubble'})**3;

# Reverse the order of the vectors.
$xBins = $xBins(-1:0);
$x = $x(-1:0);
$y = $y(-1:0);
$error = $error(-1:0);

# Read galaxy data and construct mass function.
$xGalacticus = $xBins;
$dataSet{'tree'} = "all";
&HDF5::Get_Dataset(\%dataSet,['volumeWeight','magnitudeTotal:UKIRT_K:rest:z0.0000:dustAtlas:vega']);
$dataSets  = \%{$dataSet{'dataSets'}};
$magnitude = ${$dataSets->{'magnitudeTotal:UKIRT_K:rest:z0.0000:dustAtlas:vega'}};
$weight    = ${$dataSets->{'volumeWeight'}};
($yGalacticus,$errorGalacticus) = &Histograms::Histogram($xGalacticus,$magnitude,$weight,differential => 1);

# Compute chi^2.
$chiSquared = sum(($yGalacticus-$y)**2/($errorGalacticus**2+$error**2));
$degreesOfFreedom = nelem($y);
if ( $showFit == 1 ) {
    $fitData{'name'} = "Cole et al. (2001) K-band luminosity function";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make the plot.
$plot1  = Graphics::GnuplotIF->new();
$plot1->gnuplot_hardcopy( '| ps2pdf - '.$outputFile, 
			  'postscript enhanced', 
			  'color lw 3 solid' );
$plot1->gnuplot_set_xlabel("M_{K,vega}");
$plot1->gnuplot_set_ylabel("dn/dlogM_{K,vega} [Mpc^{-3}]");
$plot1->gnuplot_set_title("K Luminosity Function");
$plot1->gnuplot_cmd("set label \"{/Symbol c}^2=".FormatSigFigs($chiSquared,4)." [".$degreesOfFreedom."]\" at screen 0.6, screen 0.2");
$plot1->gnuplot_cmd("set key left");
$plot1->gnuplot_cmd("set logscale y");
$plot1->gnuplot_cmd("set mxtics 2");
$plot1->gnuplot_cmd("set mytics 10");
$plot1->gnuplot_cmd("set format y \"10^{\%L}\"");
$plot1->gnuplot_cmd("set pointsize 1.0");
$plot1->gnuplot_cmd("plot '-' with errorbars lt 1 pt 6 title \"".$data->{'luminosityFunction'}->{'label'}."\", '-' with errorbars lt 2 pt 4 title \"Galacticus\"");
for ($i=0;$i<nelem($x);++$i) {
   $plot1->gnuplot_cmd($x->index($i)." ".$y->index($i)." ".$error->index($i));
}
$plot1->gnuplot_cmd("e");
for ($i=0;$i<nelem($xGalacticus);++$i) {
   $plot1->gnuplot_cmd($xGalacticus->index($i)." ".$yGalacticus->index($i)." ".$errorGalacticus->index($i));
}
$plot1->gnuplot_cmd("e");

exit;
