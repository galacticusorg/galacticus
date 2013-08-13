#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
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
die("Plot_Black_Hole_vs_Bulge_Mass.pl <galacticusFile> <outputDir/File> [<showFit>]")
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
my $outputFile;
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/Black_Hole_vs_Bulge_Mass.pdf";
}
(my $fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/$1/;

# Define mass bins.
my $logSpheroidMassPoints = pdl 20;
my $logSpheroidMassMin    = pdl 9.0;
my $logSpheroidMassMax    = pdl 12.0;
my $logSpheroidMassBin    = pdl ($logSpheroidMassMax-$logSpheroidMassMin)/$logSpheroidMassPoints;
my $logSpheroidMassBins   = pdl (0..$logSpheroidMassPoints-1)*$logSpheroidMassBin+$logSpheroidMassMin+0.5*$logSpheroidMassBin;

# Create data structure to read the results.
my $dataSet;
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['mergerTreeWeight','spheroidMassStellar','blackHoleMass']);
my $dataSets         = $dataSet->{'dataSets'};
my $mergerTreeWeight     = where($dataSets->{'mergerTreeWeight'}       ,$dataSets->{'spheroidMassStellar'} > 3.0e8);
my $spheroidMass     = where($dataSets->{'spheroidMassStellar'},$dataSets->{'spheroidMassStellar'} > 3.0e8);
my $blackHoleMass    = where($dataSets->{'blackHoleMass'}      ,$dataSets->{'spheroidMassStellar'} > 3.0e8);
unless (exists($dataSets->{'blackHoleMass'})) {
    if ( $showFit == 1 ) {
	my %fitData;
	$fitData{'name'} = "Haering & Rix (2003) black hole vs. bulge mass relation";
	$fitData{'chiSquared'} = "undefined";
	$fitData{'degreesOfFreedom'} = "undefined";
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
	print $xmlOutput->XMLout(\%fitData);
    }
} else {
    my $weight                     = where($mergerTreeWeight    ,($spheroidMass > 0.0) & ($blackHoleMass > 0.0));
    my $logSpheroidMassGalacticus  = where(log10($spheroidMass ),($spheroidMass > 0.0) & ($blackHoleMass > 0.0));
    my $logBlackHoleMassGalacticus = where(log10($blackHoleMass),($spheroidMass > 0.0) & ($blackHoleMass > 0.0));
    my $logBlackHoleMassMeanGalacticus;
    my $logBlackHoleMassMeanErrorGalacticus;
    my $logBlackHoleMassSigmaGalacticus;
    if ( nelem($logSpheroidMassGalacticus) > 0 ) {
	($logBlackHoleMassMeanGalacticus, $logBlackHoleMassMeanErrorGalacticus, $logBlackHoleMassSigmaGalacticus, my $logBlackHoleMassSigmaErrorGalacticus)
	    = &Means::BinnedMean($logSpheroidMassBins,$logSpheroidMassGalacticus,$logBlackHoleMassGalacticus,$weight);
    } else {
	$logBlackHoleMassMeanGalacticus      = pdl zeroes(nelem($logSpheroidMassBins));
	$logBlackHoleMassMeanErrorGalacticus = pdl zeroes(nelem($logSpheroidMassBins));
    }

    # Define constants.
    my $solarMass = pdl 1.98892e30;
    my $kilo      = pdl 1.0e3;

    # Read the XML data file.
    my $x      = pdl [];
    my $y      = pdl [];
    my $yError = pdl [];
    my $xml = new XML::Simple;
    my $data = $xml->XMLin($galacticusPath."data/observations/blackHoles/Black_Hole_Mass_vs_Galaxy_Properties_Feoli_Mancini_2009.xml", KeyAttr => "");
    my %cosmology;
    foreach my $parameter ( @{$data->{'cosmology'}->{'parameter'}} ) {
	$cosmology{$parameter->{'name'}} = $parameter->{'value'};
    }
    foreach my $system ( @{$data->{'galaxies'}->{'system'}} ) {
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
    my $xError           = $x*(10.0**0.18-1.0);
    my $logSpheroidMass  = log10($x);
    my $logBlackHoleMass = log10($y);
    my $logError         = $yError/$y/log(10.0);
    my $weights          = 1.0/$logError**2;
    (my $logBlackHoleMassMean, my $logBlackHoleMassMeanError, my $logBlackHoleMassSigma, my $logBlackHoleMassSigmaError)
	= &Means::BinnedMean($logSpheroidMassBins,$logSpheroidMass,$logBlackHoleMass,$weights);

    # Compute chi^2.
    my $degreesOfFreedom;
    my $chiSquared;
    if ( nelem($logSpheroidMass) > 0 ) {
	my $nonzero = which($logBlackHoleMassMeanGalacticus > 0.0);
	$degreesOfFreedom = nelem($nonzero);
	$chiSquared = sum((($logBlackHoleMassMean->index($nonzero)-$logBlackHoleMassMeanGalacticus->index($nonzero))**2)/($logBlackHoleMassMeanError->index($nonzero)**2+$logBlackHoleMassMeanErrorGalacticus->index($nonzero)**2));
    } else {
	$chiSquared = 0.0;
	$degreesOfFreedom = 0;
    }
    
    if ( $showFit == 1 ) {
	my %fitData;
	$fitData{'name'} = "Feoli & Mancini (2009) black hole vs. bulge mass relation";
	$fitData{'chiSquared'} = $chiSquared;
	$fitData{'degreesOfFreedom'} = $degreesOfFreedom;
	$fitData{'fileName'} = $fileName;
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
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
    my $spheroidMassBins   =  10.0** $logSpheroidMassBins;
    my $blackHoleMassMean  =  10.0** $logBlackHoleMassMeanGalacticus;
    my $blackHoleMassUpper = +10.0**($logBlackHoleMassMeanGalacticus+$logBlackHoleMassSigmaGalacticus)-$blackHoleMassMean;
    my $blackHoleMassLower = -10.0**($logBlackHoleMassMeanGalacticus-$logBlackHoleMassSigmaGalacticus)+$blackHoleMassMean;
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
