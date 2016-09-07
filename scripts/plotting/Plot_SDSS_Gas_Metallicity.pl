#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Math::SigFigs;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require Stats::Percentiles;
require XMP::MetaData;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Get name of input and output files.
die("Plot_SDSS_Gas_Metallicity.pl <galacticusFile> <outputDir/File> [<showFit>]")
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
    $outputFile = $outputTo."/SDSS_Gas_Metallicity.pdf";
}
(my $fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/$1/;

# Specify Solar metallicity and oxygen abundance.
my $solarMetallicity     = pdl 0.0189;
my $solarOxygenAbundance = pdl 4.8977e-4;

# Minimum gas fraction for galaxies to be considered in this plot.
my $gasFractionMinimum   = pdl 0.1;

# Initialize chi^2 accumulator.
my $chiSquared      = 0.0;
my $degreesOfFreedom = 0;

# Create data structure to read the results.
my $dataBlock;
$dataBlock->{'file'}  = $galacticusFile;
$dataBlock->{'store'} = 0;
&HDF5::Get_Parameters($dataBlock    );
&HDF5::Count_Trees   ($dataBlock    );
&HDF5::Select_Output ($dataBlock,0.1);
$dataBlock->{'tree'} = "all";
&HDF5::Get_Dataset($dataBlock,['mergerTreeWeight'
			      ,'magnitudeTotal:SDSS_g:observed:z0.1000:dustAtlas[faceOn]:AB'
			      ,'magnitudeTotal:SDSS_z:observed:z0.1000:AB'
			      ,'diskMassStellar'
			      ,'spheroidMassStellar'
			      ,'diskMassGas'
			      ,'spheroidMassGas'
			      ,'diskAbundancesGasMetals'
			      ,'spheroidAbundancesGasMetals'
		   ]);
my $dataSets = $dataBlock->{'dataSets'};
my $gasFraction    = ($dataSets->{'diskMassGas'}+$dataSets->{'spheroidMassGas'})/($dataSets->{'diskMassGas'}+$dataSets->{'spheroidMassGas'}+$dataSets->{'diskStellarMass'}+$dataSets->{'spheroidMassStellar'});
my $gasMetallicity = where(12.0+log10(($dataSets->{'diskAbundancesGasMetals'}+$dataSets->{'spheroidAbundancesGasMetals'})/($dataSets->{'diskMassGas'}+$dataSets->{'spheroidMassGas'}))-log10($solarMetallicity)+log10($solarOxygenAbundance),$gasFraction > $gasFractionMinimum);

# Read the XML data file.
my @tmpFiles;
my $xml = new XML::Simple;
my $data = $xml->XMLin($galacticusPath."data/observations/abundances/Gas_Phase_Metallicities_SDSS_Tremonti_2004.xml");
my $iDataset = 0;
foreach my $dataSet ( @{$data->{'gasMetallicity'}} ) {
    ++$iDataset;
    my $columns = $dataSet->{'columns'};
    my $x = pdl @{$columns->{'magnitude'}->{'data'}};
    $x = $x-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataBlock->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant'}->{'value'});

    # Compute the distribution of Galacticus galaxies.
    my $filter = $columns->{'magnitude'}->{'filter'};
    my $dust   = $columns->{'magnitude'}->{'dust'  };
    my $dustLabel;
    if ( $dust eq "corrected" ) {$dustLabel = ""};
    if ( $dust eq "face-on" )   {$dustLabel = ":dustAtlas[faceOn]"};
 
    my $property    = "magnitudeTotal:".$filter.":observed:z0.1000".$dustLabel.":AB";
    my $magnitude   = where($dataSets->{$property}     ,$gasFraction > $gasFractionMinimum);
    my $weight      = where($dataSets->{'mergerTreeWeight'},$gasFraction > $gasFractionMinimum);
    my $percentiles = pdl [2.5,16.0,50.0,84.0,97.5];
    my $results     = &Percentiles::BinnedPercentiles(
	$x,
	$magnitude,
	$gasMetallicity,
	$weight,
	$percentiles,
	);
    # Make the plot.
    my $plot;
    (my $plotFileTeX = $outputFile) =~ s/\.pdf$/_$filter.tex/;
    open(my $gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set xlabel \"".$columns->{'magnitude'}->{'label'}."\"\n";
    print $gnuPlot "set ylabel \"12+[O/H]\"\n";
    print $gnuPlot "set title offset 0,-1 \"".$dataSet->{'title'}."\"\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.5,0.2\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set mxtics 2\n";
    print $gnuPlot "set mytics 2\n";
    print $gnuPlot "set xrange [".($iDataset == 1 ? "-23:-16" : "-24:-17")."]\n";
    print $gnuPlot "set yrange [7:10]\n";
    print $gnuPlot "set pointsize 1.0\n";
    my $iPercentile           = 0;
    my $chiSquaredRange       = 0.0;
    my $degreesOfFreedomRange = 0;
    my $observationalError = pdl $dataSet->{'distributionError'};
    foreach my $percentile ( @{$columns->{'distributionPercentile'}} ) {
	++$iPercentile;
	# Compute chi^2.
	my $yData               = pdl @{$percentile->{'data'}};
	my $yGalacticus         = $results(:,($iPercentile-1));
	$chiSquaredRange       += sum((($yData-$yGalacticus)/$observationalError)**2);
	$degreesOfFreedomRange += nelem($yData);
    }
    # Accumulate chi^2.
    $chiSquared       += $chiSquaredRange;
    $degreesOfFreedom += $degreesOfFreedomRange;
    # Plot datasets.
    my $yLow1  = pdl @{$columns->{'distributionPercentile'}->[0]->{'data'}};
    my $yHigh1 = pdl @{$columns->{'distributionPercentile'}->[4]->{'data'}};
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x,$yLow1,y2 => $yHigh1,
	style      => "filledCurve",
	symbol     => [6,7],
	weight     => [5,3],
	pointSize  => 0.5,
	color      => $PrettyPlots::colorPairs{'lightSkyBlue'}
	);
    my $yLow2  = pdl @{$columns->{'distributionPercentile'}->[1]->{'data'}};
    my $yHigh2 = pdl @{$columns->{'distributionPercentile'}->[3]->{'data'}};
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x,$yLow2,y2 => $yHigh2,
	style      => "filledCurve",
	symbol     => [6,7],
	weight     => [5,3],
	pointSize  => 0.5,
	color      => $PrettyPlots::colorPairs{'cornflowerBlue'}
	);
    my $yMedian = pdl @{$columns->{'distributionPercentile'}->[2]->{'data'}};
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x,$yMedian,
	style      => "line",
	symbol     => [6,7],
	weight     => [5,3],
	pointSize  => 0.5,
	color      => $PrettyPlots::colorPairs{'darkSlateBlue'}
	);
    my $nonZero1         = which(($results(:,(0)) > 0.0) & ($results(:,(4)) > 0.0));
    my $yGalacticusLow1  = $results($nonZero1,(0));
    my $yGalacticusHigh1 = $results($nonZero1,(4));
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x->($nonZero1),$yGalacticusLow1,y2 => $yGalacticusHigh1,
	style        => "filledCurve",
	symbol       => [6,7],
	weight       => [5,3],
	pointSize    => 0.5,
	color        => $PrettyPlots::colorPairs{'salmon'},
	transparency => 0.5
	);
    my $nonZero2         = which(($results(:,(1)) > 0.0) & ($results(:,(3)) > 0.0));
    my $yGalacticusLow2  = $results($nonZero2,(1));
    my $yGalacticusHigh2 = $results($nonZero2,(3));
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x->($nonZero2),$yGalacticusLow2,y2 => $yGalacticusHigh2,
	style        => "filledCurve",
	symbol       => [6,7],
	weight       => [5,3],
	pointSize    => 0.5,
	color        => $PrettyPlots::colorPairs{'orange'},
	transparency => 0.5
	);
    my $nonZero3         = which($results(:,(2)) > 0.0);
    my $yGalacticusMedian = $results($nonZero3,(2));
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x->($nonZero3),$yGalacticusMedian,
	style        => "line",
	symbol       => [6,7],
	weight       => [5,3],
	pointSize    => 0.5,
	color        => $PrettyPlots::colorPairs{'redYellow'},
	transparency => 0.5
	);
    # Finalize plotting.
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    # Close the pipe to GnuPlot.
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileTeX);
    &MetaData::Write($outputFile,$galacticusFile,$self);    
}

# Display chi^2 information
if ( $showFit == 1 ) {
    my %fitData;
    $fitData{'name'} = "Tremonti et al. (2004) SDSS gas-phase metallicity distributions";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

exit;
