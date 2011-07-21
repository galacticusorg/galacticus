#!/usr/bin/env perl
use lib "./perl";
use PDL;
use PDL::NiceSlice;
use Galacticus::HDF5;
use Stats::Percentiles;
use Data::Dumper;

# Plot the distribution of spheroid sizes in elliptical galaxies from any number of specified models.
# Andrew Benson (11-November-2010)

# Get parameters.
if ( $#ARGV < 2 ) {die("Plot_Spheroid_Sizes.pl <galacticusFile> <label> [<galacticusFile2> <label2> ....] <outputFile>")};

# Get the name of the output file.
$outputFile = $ARGV[$#ARGV];

# Specify files to read.
for($i=0;$i<$#ARGV;$i+=2) {
    $models{$ARGV[$i+1]} = $ARGV[$i];
}

# Define constants.
$kilo                     = pdl 1000.0;
$effectiveRadiusHernquist = pdl 1.8153;

# Define morphological selection.
$ellipticalCut = pdl 0.9;

# Specify bins in which to compute medians.
$logarithmicMassMinimum  = pdl  9.0;
$logarithmicMassMaximum  = pdl 11.5;
$logarithmicMassBinWidth = pdl  0.2;
$binCount                = pdl int(($logarithmicMassMaximum-$logarithmicMassMinimum)/$logarithmicMassBinWidth+1.0);
$logarithmicMassBins     = sequence($binCount->list())*$logarithmicMassBinWidth+$logarithmicMassMinimum;
$massBins                = 10.0**$logarithmicMassBins;

# Compute the radius fitting function.
$alpha     = pdl 0.56;
$beta      = pdl 3.47e-6;
$radiusFit = $beta*($massBins**$alpha);

# Specify percentiles to compute;
$percentiles = pdl [ 10.0, 50.0, 90.0 ];

# Open a pipe to GnuPlot.
open(gnuPlot,"|gnuplot");
print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
print gnuPlot "set output \"tmp.ps\"\n";

# Make the plot.
print gnuPlot "set xlabel \"M_* [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]\"\n";
print gnuPlot "set ylabel \"r_{1/2} [kpc]\"\n";
print gnuPlot "set title 'Elliptical size vs. mass relation at z=0'\n";
print gnuPlot "set key right\n";
print gnuPlot "set logscale xy\n";
print gnuPlot "set mxtics 10\n";
print gnuPlot "set format x '10^{\%L}'\n";
print gnuPlot "set mytics 10\n";
print gnuPlot "set format y '10^{\%L}'\n";
print gnuPlot "set xrange [3.0e8:1.0e12]\n";
print gnuPlot "set yrange [0.3:30.0]\n";
print gnuPlot "set pointsize 1.0\n";

# Initialize plotting command.
$plotCommand = "plot '-' title 'Observed' with lines";
$joiner      = ", ";
for ($i=0;$i<nelem($massBins);++$i) {
    $plotData .= $massBins->index($i)." ".$radiusFit->index($i)."\n";
}
$plotData .= "e\n";

# Loop over input files.
foreach $model ( keys(%models) ) {
    # Get the input file name;
    $inputFile = $models{$model};
    # Create data structure to read the results.
    undef($dataSet);
    undef($dataSets);
    $dataSet->{'file'}  = $inputFile;
    $dataSet->{'store'} = 0;
    &HDF5::Get_Parameters($dataSet    );
    &HDF5::Count_Trees   ($dataSet    );
    &HDF5::Select_Output ($dataSet,0.0);
    $dataSet->{'tree'}  = "all";
    &HDF5::Get_Dataset($dataSet,['volumeWeight','diskStellarMass','spheroidStellarMass','spheroidScaleLength']);
    $dataSets         = $dataSet->{'dataSets'};

    # Compute the morphology.
    $morphology       = $dataSets->{'spheroidStellarMass'}/($dataSets->{'diskStellarMass'}+$dataSets->{'spheroidStellarMass'});

    # Select galaxies that are ellipticals.
    $volumeWeight     = where($dataSets->{'volumeWeight'       },$morphology > $ellipticalCut);
    $diskMass         = where($dataSets->{'diskStellarMass'    },$morphology > $ellipticalCut);
    $spheroidMass     = where($dataSets->{'spheroidStellarMass'},$morphology > $ellipticalCut);
    $spheroidRadius   = where($dataSets->{'spheroidScaleLength'},$morphology > $ellipticalCut);

    # Compute the total stellar mass and its logarithm.
    $mass             = $diskMass+$spheroidMass;
    $logarithmicMass  = log10($mass);

    # Convert spheroid radius to projected effective radius in kpc.
    $spheroidRadius  *= $kilo*$effectiveRadiusHernquist;

    # Compute the median range.
    $results = &Percentiles::BinnedPercentiles(
        $logarithmicMassBins,
        $logarithmicMass,
        $spheroidRadius,
        $volumeWeight,
        $percentiles
        );

    # Add this plot to the plot command.
    $plotCommand .= $joiner."'-' title '".$model."' with errorbars";
    $joiner = ", ";

    # Add data to the plot data.
    for ($i=0;$i<nelem($massBins);++$i) {
	$plotData .= $massBins->index($i)." ".$results(($i),(1))." ".$results(($i),(0))." ".$results(($i),(2))."\n";
    }
    $plotData .= "e\n";
    
}

# Send the plot commands to GnuPlot.
print gnuPlot $plotCommand."\n";
print gnuPlot $plotData;

# Close the pipe to GnuPlot.
close(gnuPlot);

# Convert to PDF.
system("ps2pdf tmp.ps ".$outputFile);

# Clean up files.
unlink("tmp.ps");

exit;
