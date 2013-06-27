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
use PDL::Basic;
use Graphics::GnuplotIF;
use XML::Simple;
require Galacticus::HDF5;
require Stats::Histograms;

# Run a test of the progenitor mass function construction algorithms.
# Andrew Benson (12-Feb-2010)

# Write startup message.
print "Progenitor_Mass_Function.pl is running: will create a plot of progenitor mass functions from Galacticus...\n";

# Specify name of Galacticus output file to be used.
my $galacticusOutput = $galacticusPath."tests/progenitorMassFunction/progenitorMassFunctionTest.hdf5";

# Hubble constant.
my $h0 = 0.73;

# Specify range and bin sizes in logarithmic mass ratio for conditional mass functions.
my $lgMmin  = -6.0;
my $lgMmax  =  1.0;
my $lgMstep =  0.2;

# Specify number of outputs (including z=0).
my $outputCount = 5;

# Specify range of M2 (final mass) bins.
my @rootBins = ( 12.0, 13.5, 15.0 );
my $rootBinStep = 0.301029996;
my $rootBinCount = $#rootBins+1;

# Create array of masses for conditional mass function.
my $dummy = pdl[1..10];
(my $lgMval, my $hist) = hist($dummy,$lgMmin,$lgMmax,$lgMstep);
my $Mval = 10.0**$lgMval;

# Initialize PDLs to store progenitor mass function and accumulated weights.
my $progenitorMF  = zeroes($rootBinCount,$outputCount-1,nelem($hist));
my $summedWeights = zeroes($rootBinCount,$outputCount-1);

# Read data from the Millennium Simulation.
my $xml = new XML::Simple;
my $data = $xml->XMLin($galacticusPath."data/darkMatter/Progenitor_Mass_Function_Millennium_Simulation.xml");
my %millenniumData;
foreach my $massFunction ( @{$data->{'massFunction'}} ) {
    my $label = $massFunction->{'rootMass'}.":".$massFunction->{'redshift'};
    @{${$millenniumData{$label}}{'Mass'}} = @{$massFunction->{'massRatio'}};
    @{${$millenniumData{$label}}{'Fcmf'}} = @{$massFunction->{'Fcmf'}};
}

# Run Galacticus to generate the data.
print "  -> Running Galacticus to generate merger trees...\n";
system($galacticusPath."Galacticus.exe ".$galacticusPath."tests/progenitorMassFunction/Progenitor_Mass_Function_Parameters.xml");

# Create data structure to read the results.
my $dataSet;
$dataSet->{'file'} = $galacticusOutput;

# Get a count of the number of trees present.
&HDF5::Count_Trees($dataSet);
my $treesCount = scalar(@{$dataSet->{'mergerTreesAvailable'}});
print "  -> Found ".$treesCount." trees: processing.......\n";

# Loop through trees.
for (my $iTree=1;$iTree<=$treesCount;$iTree+=1) {
    $dataSet->{'tree'} = $iTree;
    # Loop over outputs.
    for (my $iOutput=$outputCount;$iOutput>0;--$iOutput) {
	$dataSet->{'output'} = $iOutput;
	# Read the node masses and which nodes are isolated.
	&HDF5::Get_Dataset($dataSet,['basicMass','nodeIsIsolated','mergerTreeWeight']);
	my $dataSets = $dataSet->{'dataSets'};
	# Get a list of isolated node masses.
	my $isolatedNodeMass = where($dataSets->{'basicMass'},$dataSets->{'nodeIsIsolated'} == 1);
	$isolatedNodeMass = $isolatedNodeMass*$h0; # Put into "h" units as were used by Cole et al.
	my $weight = $isolatedNodeMass;
	$isolatedNodeMass = log10($isolatedNodeMass);

	# Clean up.
	delete($dataSets->{'nodeIsIsolated'});
	delete($dataSets->{'basicMass'});
	# Compute fractional node masses and take the log.
	my $logRootNodeMass;
	my $rootNodeMass;
	if ( $iOutput == $outputCount ) {
	    $logRootNodeMass = $isolatedNodeMass->index(0);
	    $rootNodeMass = $weight->index(0);
	}
	$isolatedNodeMass = $isolatedNodeMass-$logRootNodeMass;
	$weight = $weight/$rootNodeMass;

	if ( $iOutput < $outputCount ) {
	    # Determine in which tree bin this should lie.
	    my $rootBin = -1;
	    for(my $iRootBin=0;$iRootBin<$rootBinCount;++$iRootBin) {
		if ( $logRootNodeMass-$rootBins[$iRootBin] > -0.5*$rootBinStep && $logRootNodeMass-$rootBins[$iRootBin] <= 0.5*$rootBinStep ) {$rootBin=$iRootBin};
	    }
	    if ( $rootBin >= 0 && $rootBin < $rootBinCount) {
		# Build a histogram.
		(my $hist, my $histErrors) = &Histograms::Histogram($lgMval,$isolatedNodeMass,$weight);
		# Accumulate.
		my $vWeight = $dataSets->{'mergerTreeWeight'}->index(0);
		$progenitorMF->(($rootBin),($iOutput-1),:) += $hist*$vWeight;
		$summedWeights->(($rootBin),($iOutput-1)) += $vWeight;
	    }
	}
	delete($dataSets->{'mergerTreeWeight'});
    }
}

# Divide through by weights and bin width.
for (my $rootBin=0;$rootBin<$rootBinCount;++$rootBin) {
    for (my $iOutput=$outputCount-1;$iOutput>0;--$iOutput) {
	if ( $summedWeights->(($rootBin),($iOutput-1)) > 0.0 ) {$progenitorMF->(($rootBin),($iOutput-1),:) /= $lgMstep*$summedWeights->(($rootBin),($iOutput-1))};
    }
}

# Create the plot.
print "  -> Creating the plot...\n";
my $plotName = $galacticusPath."tests/progenitorMassFunction/progenitorMassFunction.pdf";
my $plot1  = Graphics::GnuplotIF->new();
$plot1->gnuplot_hardcopy( '| ps2pdf - '.$plotName, 
			  'postscript enhanced', 
			  'color lw 3' );
$plot1->gnuplot_cmd("set multiplot");
$plot1->gnuplot_cmd("set size 0.33,0.25");
$plot1->gnuplot_cmd("set label \"M_1/M_2\" at screen 0.5, screen -0.08 center");
$plot1->gnuplot_cmd("set label \"df_{CMF}/dlog_{10}M_1\" at screen -0.06, screen 0.5 center rotate");
for (my $rootBin=0;$rootBin<$rootBinCount;++$rootBin) {
    for (my $iOutput=$outputCount-1;$iOutput>0;--$iOutput) {
	my $ox = 0.33*$rootBin;
	my $oy = 1.0-0.25*$iOutput;
	$plot1->gnuplot_cmd("set origin ".$ox.",".$oy);
	$plot1->gnuplot_cmd("set lmargin 0");
	$plot1->gnuplot_cmd("set rmargin 0");
	$plot1->gnuplot_cmd("set tmargin 0");
	$plot1->gnuplot_cmd("set bmargin 0");
	$plot1->gnuplot_cmd("set logscale xy");
	$plot1->gnuplot_cmd("unset key");
	if ( $rootBin == 0 ) {
	    $plot1->gnuplot_cmd("set format y \"10^{\%L}\"");
	} else {
	    $plot1->gnuplot_cmd("set format y \"\"");
	}
	if ( $iOutput == $outputCount-1 ) {
	    $plot1->gnuplot_cmd("set format x \"10^{\%L}\"");
	} else {
	    $plot1->gnuplot_cmd("set format x \"\"");
	}
	$plot1->gnuplot_cmd("set yrange [0.01:3.0]");
	my @x = list($Mval);
	my @y = list($progenitorMF->(($rootBin),($iOutput-1),:));
	my $M2;
	if ( $rootBin == 0 ) {
	    $M2 = "1.000e+12";
	    $plot1->gnuplot_cmd("set xrange [0.003:2.0]");
	}
	if ( $rootBin == 1 ) {
	    $M2 = "3.160e+13";
	    $plot1->gnuplot_cmd("set xrange [0.0001:2.0]");
	}
	if ( $rootBin == 2 ) {
	    $M2 = "1.000e+15";
	    $plot1->gnuplot_cmd("set xrange [0.000003:2.0]");
	}
	my $z1;
	if ( $iOutput == 4 ) {$z1 = "0.5"};
	if ( $iOutput == 3 ) {$z1 = "1.0"};
	if ( $iOutput == 2 ) {$z1 = "2.0"};
	if ( $iOutput == 1 ) {$z1 = "4.0"};
	my $label = $M2.":".$z1;
	my @xM = @{${$millenniumData{$label}}{'Mass'}};
	my @yM = @{${$millenniumData{$label}}{'Fcmf'}};
	$plot1->gnuplot_plot_many( \@x, \@y, \@xM, \@yM );
    }
}
$plot1->gnuplot_cmd("unset multiplot");
print "  -> Plot is available in: ".$plotName."\n";

# Remove the Galacticus output.
unlink($galacticusOutput);

exit;
