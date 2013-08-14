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
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use Data::Dumper;
require Galacticus::HDF5;
require Stats::Histograms;

# Plot various statistics of galaxies computed with N-body and Monte Carlo merger trees for comparison.
# Andrew Benson (4-October-2010).

# Specify the models to plot, including labels and model number.
my @models = (
    { label => "N-body; subhalo merger times + targets", number => 1 },
    { label => "N-body; subhalo merger times"          , number => 2 },
    { label => "N-body; analytic merger times"         , number => 3 },
    { label => "Monte Carlo"                           , number => 4 }
    );

# Specify the list of properties to read from each model.
my @properties = ( "mergerTreeWeight", "diskStellarMass", "spheroidStellarMass", "diskScaleLength", "spheroidScaleLength", "diskGasMass", "spheroidGasMass", "blackHoleMass", "diskStellarMetals", "spheroidStellarMetals", "diskGasMetals", "spheroidGasMetals", "nodeMass", "nodeIsIsolated", "timeToMerge" );

# Specify a list of plots to make.
my @plots = (
	  { fileName => "tests/nBodyMergerTrees/plots/diskStellarMass.pdf"    , xMin => 3.0, xMax => 12.0, xBin => 0.5, xType => "log", x => "diskStellarMass"    , xLabel => "M_{stars, disk} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]"    , yLabel => "df/d ln M_{stars,disk}"    , title => "Disk stellar mass"     },
	  { fileName => "tests/nBodyMergerTrees/plots/spheroidStellarMass.pdf", xMin => 3.0, xMax => 12.0, xBin => 0.5, xType => "log", x => "spheroidStellarMass", xLabel => "M_{stars, spheroid} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]", yLabel => "df/d ln M_{stars,spheroid}", title => "Spheroid stellar mass" },
	  { fileName => "tests/nBodyMergerTrees/plots/diskGasMass.pdf"        , xMin => 3.0, xMax => 12.0, xBin => 0.5, xType => "log", x => "diskGasMass"        , xLabel => "M_{stars, disk} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]"    , yLabel => "df/d ln M_{stars,disk}"    , title => "Disk gas mass"         },
	  { fileName => "tests/nBodyMergerTrees/plots/spheroidGasMass.pdf"    , xMin => 3.0, xMax => 12.0, xBin => 0.5, xType => "log", x => "spheroidGasMass"    , xLabel => "M_{stars, spheroid} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]", yLabel => "df/d ln M_{stars,spheroid}", title => "Spheroid gas mass"     },
	  { fileName => "tests/nBodyMergerTrees/plots/blackHoleMass.pdf"      , xMin => 1.0, xMax => 10.0, xBin => 0.5, xType => "log", x => "blackHoleMass"      , xLabel => "M_{BH} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]"             , yLabel => "df/d ln M_{BH}"            , title => "Black hole mass"       },
	  { fileName => "tests/nBodyMergerTrees/plots/diskScaleLength.pdf"    , xMin => -4.0, xMax => 0.0, xBin => 0.25, xType => "log", x => "diskScaleLength"    , xLabel => "R_{disk} [Mpc]"    , yLabel => "df/d ln R_{disk}"    , title => "Disk scale length"    },
	  { fileName => "tests/nBodyMergerTrees/plots/spheroidScaleLength.pdf", xMin => -4.0, xMax => 0.0, xBin => 0.25, xType => "log", x => "spheroidScaleLength", xLabel => "r_{spheroid} [Mpc]", yLabel => "df/d ln r_{spheroid}", title => "Spheroid scale length"},
	  { fileName => "tests/nBodyMergerTrees/plots/diskStellarMetals.pdf"    , xMin => 2.0, xMax => 11.0, xBin => 0.5, xType => "log", x => "diskStellarMetals"    , xLabel => "M_{stars, disk} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]"    , yLabel => "df/d ln M_{stars,disk}"    , title => "Disk stellar metals"     },
	  { fileName => "tests/nBodyMergerTrees/plots/spheroidStellarMetals.pdf", xMin => 2.0, xMax => 11.0, xBin => 0.5, xType => "log", x => "spheroidStellarMetals", xLabel => "M_{stars, spheroid} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]", yLabel => "df/d ln M_{stars,spheroid}", title => "Spheroid stellar metals" },
	  { fileName => "tests/nBodyMergerTrees/plots/diskGasMetals.pdf"        , xMin => 2.0, xMax => 11.0, xBin => 0.5, xType => "log", x => "diskGasMetals"        , xLabel => "M_{stars, disk} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]"    , yLabel => "df/d ln M_{stars,disk}"    , title => "Disk gas metals"         },
	  { fileName => "tests/nBodyMergerTrees/plots/spheroidGasMetals.pdf"    , xMin => 2.0, xMax => 11.0, xBin => 0.5, xType => "log", x => "spheroidGasMetals"    , xLabel => "M_{stars, spheroid} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]", yLabel => "df/d ln M_{stars,spheroid}", title => "Spheroid gas metals"     },
	  { fileName => "tests/nBodyMergerTrees/plots/nodeMass.pdf"             , xMin => 9.0, xMax => 14.0, xBin => 0.5, xType => "log", x => "nodeMass"             , xLabel => "M_{node} [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]", yLabel => "df/d ln M_{node}", title => "Node mass"     },
	  { fileName => "tests/nBodyMergerTrees/plots/timeToMerge.pdf"          , xMin => -3.0, xMax => 2.0, xBin => 0.5, xType => "log", x => "timeToMerge"          , xLabel => "t_{merge} [Gyr]", yLabel => "df/d ln t_{merge}", title => "Time to merge"     }
	  );

# Read properties from each model.
my $iModel = -1;
my @modelValues;
foreach my $model ( @models ) {
    ++$iModel;
    my $dataSet;
    $dataSet->{'file' } = "tests/nBodyMergerTrees/models/galacticus_0:".$model->{'number'}."/galacticus.hdf5";
    $dataSet->{'store'} = 0;
    $dataSet->{'tree' } = "all";
    &HDF5::Select_Output($dataSet,0.0);
    &HDF5::Get_Dataset($dataSet,\@properties);
    my $dataSets         = $dataSet->{'dataSets'};
    foreach my $property ( @properties ) {
	${$modelValues[$iModel]}{$property} = $dataSets->{$property};
    }
}

# Make plots.
foreach my $plot ( @plots ) {

    # Open a pipe to GnuPlot.
    open(gnuPlot,"|gnuplot");
    print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
    print gnuPlot "set output \"tmp.ps\"\n";
    
    # Make the plot.
    print gnuPlot "set xlabel \"".$plot->{'xLabel'}."\"\n";
    print gnuPlot "set ylabel \"".$plot->{'yLabel'}."\"\n";
    print gnuPlot "set title '".$plot->{'title'}."'\n";
    print gnuPlot "set key right\n";
    print gnuPlot "set logscale y\n";
    print gnuPlot "set mytics 10\n";
    print gnuPlot "set format y '10^{\%L}'\n";
    if ( $plot->{'xType'} eq "log" ) {
	print gnuPlot "set logscale x\n";
	print gnuPlot "set mxtics 10\n";
	print gnuPlot "set format x '10^{\%L}'\n";
    } else {
	print gnuPlot "set mxtics 2\n";
    }
    print gnuPlot "set pointsize 1.0\n";
    print gnuPlot "plot";
    my $joiner = " ";
    foreach my $model ( @models ) {
	print gnuPlot $joiner."'-' with yerrorbars title '".$model->{'label'}."'";
	$joiner = ", ";
    }
    print gnuPlot "\n";
    foreach my $modelValue ( @modelValues ) {
	# Construct the bins to use.
	my $nBins = int(($plot->{'xMax'}-$plot->{'xMin'})/$plot->{'xBin'});
	my $bins  = pdl (0..$nBins-1)*$plot->{'xBin'}+$plot->{'xMin'}+0.5*$plot->{'xBin'};
	# Construct the x values to use.
	my $x;
	my $xBins;
	if ( $plot->{'xType'} eq "log" ) {
	    $x     = log10($modelValue->{$plot->{'x'}});
	    $xBins = 10.0**$bins;
	} else {
	    $x     = $modelValue->{$plot->{'x'}};
	    $xBins = $bins;
	}
	my $w = $modelValue->{'mergerTreeWeight'};
	(my $y, my $error) = &Histograms::Histogram($bins,$x,$w,differential => 1,normalized => 1);
	for (my $i=0;$i<nelem($xBins);++$i) {
	    print gnuPlot $xBins->index($i)." ".$y->index($i)." ".$error->index($i)."\n";
	}
	print gnuPlot "e\n";
    }
    
    # Close the pipe to GnuPlot.
    close(gnuPlot);
    
    # Convert to PDF.
    system("mkdir -p `dirname ".$plot->{'fileName'}."`;ps2pdf tmp.ps ".$plot->{'fileName'});
     
    # Clean up files.
    unlink("tmp.ps");

}

exit;
