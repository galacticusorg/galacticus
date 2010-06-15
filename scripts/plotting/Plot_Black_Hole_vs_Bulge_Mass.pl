#!/usr/bin/env perl
use lib "./perl";
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Graphics::GnuplotIF;
use Galacticus::HDF5;
use Galacticus::Magnitudes;
use Math::SigFigs;
use Stats::Means;
use Data::Dumper;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Black_Hole_vs_Bulge_Mass.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
    $outputFile = $outputTo."/Black_Hole_vs_Bulge_Mass.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Define mass bins.
$logSpheroidMassPoints = pdl 5;
$logSpheroidMassMin    = pdl 9.8;
$logSpheroidMassMax    = pdl 11.5;
$logSpheroidMassBin    = pdl ($logSpheroidMassMax-$logSpheroidMassMin)/$logSpheroidMassPoints;
$logSpheroidMassBins   = pdl (0..$logSpheroidMassPoints-1)*$logSpheroidMassBin+$logSpheroidMassMin+0.5*$logSpheroidMassBin;

# Create data structure to read the results.
$dataSet{'file'} = $galacticusFile;
$dataSet{'store'} = 0;
&HDF5::Get_Parameters(\%dataSet);
&HDF5::Count_Trees(\%dataSet);
&HDF5::Select_Output(\%dataSet,0.0);
$dataSet{'tree'} = "all";
&HDF5::Get_Dataset(\%dataSet,['volumeWeight','spheroidStellarMass','blackHoleMass']);
$dataSets         = \%{$dataSet{'dataSets'}};
$spheroidMass     = ${$dataSets->{'spheroidStellarMass'}};
$blackHoleMass    = ${$dataSets->{'blackHoleMass'}};
unless (exists($dataSets->{'blackHoleMass'})) {
    if ( $showFit == 1 ) {
	$fitData{'name'} = "Haering & Rix (2003) black hole vs. bulge mass relation";
	$fitData{'chiSquared'} = "undefined";
	$fitData{'degreesOfFreedom'} = "undefined";
	$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
	print $xmlOutput->XMLout(\%fitData);
    }
} else {
    $weight           = ${$dataSets->{'volumeWeight'}};
    $logSpheroidMass  = where(log10($spheroidMass ),$spheroidMass > 0.0 & $blackHoleMass > 0.0);
    $logBlackHoleMass = where(log10($blackHoleMass),$spheroidMass > 0.0 & $blackHoleMass > 0.0);
    if ( nelem($logSpheroidMass) > 0 ) {
	($logBlackHoleMassMeanGalacticus,$logBlackHoleMassMeanErrorGalacticus,$logBlackHoleMassSigmaGalacticus,$logBlackHoleMassSigmaErrorGalacticus)
	    = &Means::BinnedMean($logSpheroidMassBins,$logSpheroidMass,$logBlackHoleMass,$weight);
    }

# Read the XML data file.
    $xml = new XML::Simple;
    $data = $xml->XMLin("data/Black_Hole_Mass_vs_Galaxy_Properties.xml");
    $columns = $data->{'blackHoleData'}->{'columns'};
    $x = pdl @{$columns->{'bulgeMass'}->{'data'}};
    $y = pdl @{$columns->{'blackHoleMass'}->{'data'}};
    $yErrorUp = pdl @{$columns->{'blackHoleMassErrorUp'}->{'data'}};
    $yErrorDown = pdl @{$columns->{'blackHoleMassErrorDown'}->{'data'}};
    $yErrorUp = $y+$yErrorUp;
    $yErrorDown = $y-$yErrorDown;
    $x = $x*($columns->{'bulgeMass'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'})**2;
    $y = $y*($columns->{'blackHoleMass'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
    $yErrorUp = $yErrorUp*($columns->{'blackHoleMass'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
    $yErrorDown = $yErrorDown*($columns->{'blackHoleMass'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
    $logSpheroidMass       = log10($x);
    $logBlackHoleMass      = log10($y);
    $logError = 0.5*($yErrorUp-$yErrorDown)/$y/log(10.0);
    $weights = 1.0/$logError**2;
($logBlackHoleMassMean,$logBlackHoleMassMeanError,$logBlackHoleMassSigma,$logBlackHoleMassSigmaError)
    = &Means::BinnedMean($logSpheroidMassBins,$logSpheroidMass,$logBlackHoleMass,$weights);

# Compute chi^2.
if ( nelem($logSpheroidMass) > 0 ) {
    $degreesOfFreedom = 2*nelem($logBlackHoleMassMean);
    $chiSquared = sum((($logBlackHoleMassMean-$logBlackHoleMassMeanGalacticus)**2)/($logBlackHoleMassMeanError**2+$logBlackHoleMassMeanErrorGalacticus**2))
	+sum((($logBlackHoleMassSigma-$logBlackHoleMassSigmaGalacticus)**2)/($logBlackHoleMassSigmaError**2+$logBlackHoleMassSigmaErrorGalacticus**2));
} else {
    $chiSquared = 0.0;
    $degreesOfFreedom = 0;
}

if ( $showFit == 1 ) {
    $fitData{'name'} = "Haering & Rix (2003) black hole vs. bulge mass relation";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make the plot.
$plot1  = Graphics::GnuplotIF->new();
$plot1->gnuplot_hardcopy( '| ps2pdf - '.$outputFile, 
			  'postscript enhanced', 
			  'color lw 3 solid' );
$plot1->gnuplot_set_xlabel("M_{bulge} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]");
$plot1->gnuplot_set_ylabel("M_{black hole} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]");
$plot1->gnuplot_set_title("Black hole mass vs. bulge mass");
$plot1->gnuplot_cmd("set label \"{/Symbol c}^2=".FormatSigFigs($chiSquared,4)." [".$degreesOfFreedom."]\" at screen 0.6, screen 0.2");
$plot1->gnuplot_cmd("set xrange [1.0e8:1.0e13]");
$plot1->gnuplot_cmd("set yrange [1.0e6:1.0e10]");
$plot1->gnuplot_cmd("set key left");
$plot1->gnuplot_cmd("set logscale xy");
$plot1->gnuplot_cmd("set mxtics 10");
$plot1->gnuplot_cmd("set mytics 10");
$plot1->gnuplot_cmd("set format x \"10^{\%L}\"");
$plot1->gnuplot_cmd("set format y \"10^{\%L}\"");
$plot1->gnuplot_cmd("set pointsize 1.0");
$plotCommand = "plot '-' with errorbars pt 6 title \"".$data->{'blackHoleData'}->{'label'}."\"";
if ( nelem($logSpheroidMass) > 0 ) {$plotCommand .= ", '-' pt 4 title \"Galacticus\""};
$plot1->gnuplot_cmd($plotCommand);
for ($i=0;$i<nelem($x);++$i) {
    $plot1->gnuplot_cmd($x->index($i)." ".$y->index($i)." ".$yErrorDown->index($i)." ".$yErrorUp->index($i));
}
$plot1->gnuplot_cmd("e");
if ( nelem($logSpheroidMass) > 0 ) {
    for ($i=0;$i<nelem($spheroidMass);++$i) {
	$plot1->gnuplot_cmd($spheroidMass->index($i)." ".$blackHoleMass->index($i));
    }
    $plot1->gnuplot_cmd("e");
}
}

exit;
