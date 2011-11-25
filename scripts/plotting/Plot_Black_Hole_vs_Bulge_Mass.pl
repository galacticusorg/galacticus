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
use PDL::NiceSlice;
use XML::Simple;
use Graphics::GnuplotIF;
require GnuPlot::LaTeX;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
use Math::SigFigs;
require Stats::Means;
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
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['volumeWeight','spheroidStellarMass','blackHoleMass']);
$dataSets         = $dataSet->{'dataSets'};
$volumeWeight     = where($dataSets->{'volumeWeight'}       ,$dataSets->{'spheroidStellarMass'} > 3.0e8);
$spheroidMass     = where($dataSets->{'spheroidStellarMass'},$dataSets->{'spheroidStellarMass'} > 3.0e8);
$blackHoleMass    = where($dataSets->{'blackHoleMass'}      ,$dataSets->{'spheroidStellarMass'} > 3.0e8);
unless (exists($dataSets->{'blackHoleMass'})) {
    if ( $showFit == 1 ) {
	$fitData{'name'} = "Haering & Rix (2003) black hole vs. bulge mass relation";
	$fitData{'chiSquared'} = "undefined";
	$fitData{'degreesOfFreedom'} = "undefined";
	$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
	print $xmlOutput->XMLout(\%fitData);
    }
} else {
    $weight           = where($volumeWeight        ,$spheroidMass > 0.0 & $blackHoleMass > 0.0);
    $logSpheroidMass  = where(log10($spheroidMass ),$spheroidMass > 0.0 & $blackHoleMass > 0.0);
    $logBlackHoleMass = where(log10($blackHoleMass),$spheroidMass > 0.0 & $blackHoleMass > 0.0);
    if ( nelem($logSpheroidMass) > 0 ) {
	($logBlackHoleMassMeanGalacticus,$logBlackHoleMassMeanErrorGalacticus,$logBlackHoleMassSigmaGalacticus,$logBlackHoleMassSigmaErrorGalacticus)
	    = &Means::BinnedMean($logSpheroidMassBins,$logSpheroidMass,$logBlackHoleMass,$weight);
    }

    # Define constants.
    $solarMass = pdl 1.98892e30;
    $kilo      = pdl 1.0e3;

    # Read the XML data file.
    $x      = pdl [];
    $y      = pdl [];
    $yError = pdl [];
    $xml = new XML::Simple;
    $data = $xml->XMLin($galacticusPath."data/Black_Hole_Mass_vs_Galaxy_Properties_Feoli_Mancini_2009.xml", KeyAttr => "");
    foreach $parameter ( @{$data->{'cosmology'}->{'parameter'}} ) {
	$cosmology{$parameter->{'name'}} = $parameter->{'value'};
    }
    foreach $system ( @{$data->{'galaxies'}->{'system'}} ) {
	$x      = $x     ->append($system->{'spheroidMass'      });
	$y      = $y     ->append($system->{'blackHoleMass'     });
	$yError = $yError->append($system->{'blackHoleMassError'});
    }
    $x               .= $x     *$data->{'units'}->{'mass'    }->{'unitsInSI'}/$solarMass;
    $y               .= $y     *$data->{'units'}->{'velocity'}->{'unitsInSI'}/$kilo;
    $yError          .= $yError*$data->{'units'}->{'velocity'}->{'unitsInSI'}/$kilo;
    if ( exists($cosmology{'H_0'}) ) {
	$x      .= $x     *($dataSet->{'parameters'}->{'H_0'}/$cosmology{'H_0'})**$data->{'units'}->{'mass'    }->{'hubbleExponent'};
	$y      .= $y     *($dataSet->{'parameters'}->{'H_0'}/$cosmology{'H_0'})**$data->{'units'}->{'velocity'}->{'hubbleExponent'};
	$yError .= $yError*($dataSet->{'parameters'}->{'H_0'}/$cosmology{'H_0'})**$data->{'units'}->{'velocity'}->{'hubbleExponent'};
    }
    $xError           = $x*(10.0**0.18-1.0);
    $logSpheroidMass  = log10($x);
    $logBlackHoleMass = log10($y);
    $logError         = $yError/$y/log(10.0);
    $weights          = 1.0/$logError**2;
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
	$fitData{'name'} = "Feoli & Mancini (2009) black hole vs. bulge mass relation";
	$fitData{'chiSquared'} = $chiSquared;
	$fitData{'degreesOfFreedom'} = $degreesOfFreedom;
	$fitData{'fileName'} = $fileName;
	$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
	print $xmlOutput->XMLout(\%fitData);
    }
    
    # Make the plot.
    ($outputFileEPS = $outputFile) =~ s/\.pdf$/.eps/;
    open(gnuPlot,"|gnuplot");
    print gnuPlot "set terminal epslatex color lw 3 solid 7\n";
    print gnuPlot "set output '".$outputFileEPS."'\n";
    print gnuPlot "set xlabel '\$M_{\\star,\\rm bulge} [M_\\odot]\$'\n";
    print gnuPlot "set ylabel '\$M_\\bullet [M_\\odot]\$'\n";
    print gnuPlot "set title 'Black hole mass vs. bulge mass'\n";
    print gnuPlot "set label '\$\\chi^2=".FormatSigFigs($chiSquared,4)."\$ [".$degreesOfFreedom."]' at screen 0.73, screen 0.2\n";
    print gnuPlot "set key left\n";
    print gnuPlot "set logscale xy\n";
    print gnuPlot "set mxtics 10\n";
    print gnuPlot "set mytics 10\n";
    print gnuPlot "set format x '\$10^{\%L}\$'\n";
    print gnuPlot "set format y '\$10^{\%L}\$'\n";
    print gnuPlot "set pointsize 1.0\n";
    $plotCommand = "plot '-' with xyerrorbars pt 6 title \"".$data->{'label'}."\"";
    if ( nelem($logSpheroidMass) > 0 ) {$plotCommand .= ", '-' pt 4 title \"Galacticus\""};
    print gnuPlot $plotCommand."\n";
    for ($i=0;$i<nelem($x);++$i) {
	print gnuPlot $x->index($i)." ".$y->index($i)." ".$xError->index($i)." ".$yError->index($i)."\n";
    }
    print gnuPlot "e\n";
    if ( nelem($logSpheroidMass) > 0 ) {
	for ($i=0;$i<nelem($spheroidMass);++$i) {
	    print gnuPlot $spheroidMass->index($i)." ".$blackHoleMass->index($i)."\n";
	}
	print gnuPlot "e\n";
    }
}
close(gnuPlot);
&LaTeX::GnuPlot2PDF($outputFileEPS);

exit;
