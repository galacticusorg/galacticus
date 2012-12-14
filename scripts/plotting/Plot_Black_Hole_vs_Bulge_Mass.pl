#!/usr/bin/env perl
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
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
require Stats::Means;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require XMP::MetaData;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Black_Hole_vs_Bulge_Mass.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/Black_Hole_vs_Bulge_Mass.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Define mass bins.
$logSpheroidMassPoints = pdl 20;
$logSpheroidMassMin    = pdl 9.0;
$logSpheroidMassMax    = pdl 12.0;
$logSpheroidMassBin    = pdl ($logSpheroidMassMax-$logSpheroidMassMin)/$logSpheroidMassPoints;
$logSpheroidMassBins   = pdl (0..$logSpheroidMassPoints-1)*$logSpheroidMassBin+$logSpheroidMassMin+0.5*$logSpheroidMassBin;

# Create data structure to read the results.
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['volumeWeight','spheroidMassStellar','blackHoleMass']);
$dataSets         = $dataSet->{'dataSets'};
$volumeWeight     = where($dataSets->{'volumeWeight'}       ,$dataSets->{'spheroidMassStellar'} > 3.0e8);
$spheroidMass     = where($dataSets->{'spheroidMassStellar'},$dataSets->{'spheroidMassStellar'} > 3.0e8);
$blackHoleMass    = where($dataSets->{'blackHoleMass'}      ,$dataSets->{'spheroidMassStellar'} > 3.0e8);
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
    $data = $xml->XMLin($galacticusPath."data/observations/blackHoles/Black_Hole_Mass_vs_Galaxy_Properties_Feoli_Mancini_2009.xml", KeyAttr => "");
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
	$nonzero = which($logBlackHoleMassMeanGalacticus > 0.0);
	$degreesOfFreedom = nelem($nonzero);
	$chiSquared = sum((($logBlackHoleMassMean->index($nonzero)-$logBlackHoleMassMeanGalacticus->index($nonzero))**2)/($logBlackHoleMassMeanError->index($nonzero)**2+$logBlackHoleMassMeanErrorGalacticus->index($nonzero)**2));
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
    
    # Make plot of stellar mass function.
    my $plot;
    my $gnuPlot;
    my $plotFile = $outputFile;
    (my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
    open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
    print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
    print $gnuPlot "set output '".$plotFileEPS."'\n";
    print $gnuPlot "set title 'Black Hole Mass vs. Bulge Stellar Mass \$z=0\$'\n";
    print $gnuPlot "set xlabel 'Bulge stellar mass; \$M_{\\rm bulge} [M_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Black hole mass; \$M_\\bullet [M_\\odot]\$'\n";
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
    print $gnuPlot "set xrange [3.0e8:3.0e12]\n";
    print $gnuPlot "set yrange [1.0e5:1.0e10]\n";
    print $gnuPlot "set pointsize 2.0\n";
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x,$y,
	errorLeft  => $xError,
	errorRight => $xError,
	errorUp    => $yError,
	errorDown  => $yError,
	style      => "point",
	symbol     => [6,7],
	weight     => [5,3],
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[0]},
	title      => $data->{'label'}.' [observed]'
	);
    $spheroidMassBins   =  10.0** $logSpheroidMassBins;
    $blackHoleMassMean  =  10.0** $logBlackHoleMassMeanGalacticus;
    $blackHoleMassUpper = +10.0**($logBlackHoleMassMeanGalacticus+$logBlackHoleMassSigmaGalacticus)-$blackHoleMassMean;
    $blackHoleMassLower = -10.0**($logBlackHoleMassMeanGalacticus-$logBlackHoleMassSigmaGalacticus)+$blackHoleMassMean;
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$spheroidMassBins,$blackHoleMassMean,
	errorDown  => $blackHoleMassLower,
	errorUp    => $blackHoleMassUpper,
	style      => "point",
	symbol     => [6,7],
	weight     => [5,3],
	color      => $PrettyPlots::colorPairs{'redYellow'},
	title      => 'Galacticus'
	);
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileEPS);
    &MetaData::Write($plotFile,$galacticusFile,$self);
}

exit;
