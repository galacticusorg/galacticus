#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
use XML::Simple;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Visualize the posterior distribution from a BIE statelog file.
# Andrew Benson (10-June-2012)

# Get file name to process.
die("Usage: mcmcVisualize.pl fileRoot configFileName [options]")
    unless ( scalar(@ARGV) > 1 );
my $fileRoot       = $ARGV[0];
my $configFileName = $ARGV[1];

# Create a hash of named arguments.
my %arguments;
my $iArg = 0;
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	my $argument = $1;
	my $value = $ARGV[$iArg+1];
	if ( $ARGV[$iArg] =~ m/^\"/ ) {
	    until ( $value =~ m/\"$/ ) {
		++$iArg;
		$value .= $ARGV[$iArg];
	    }
	}
	++$iArg;
	if ( $argument eq "range" ) {
	    my @rangeSpecifiers = split(/:/,$value);
	    my %rangeData = 
		(
		 name  => $rangeSpecifiers[0],
		 lower => $rangeSpecifiers[1],
		 upper => $rangeSpecifiers[2],
		);		
	    push(@{$arguments{$argument}},\%rangeData);
	} else {
	    $arguments{$argument} = $value;
	    $arguments{$argument} =~ s/^\"(.*)\"$/$1/;
	}
    }
}

# Extract options.
my $workDirectory = ".";
$workDirectory = $arguments{'workDirectory'}
    if ( exists($arguments{'workDirectory'}) );
my $nGrid = 100;
$nGrid = $arguments{'ngrid'}
    if ( exists($arguments{'ngrid'}) );
my $xScale = "linear";
$xScale = $arguments{'xScale'}
    if ( exists($arguments{'xScale'}) );
my $xLabel = "x";
$xLabel = $arguments{'xLabel'}
    if ( exists($arguments{'xLabel'}) );
my $yScale = "linear";
$yScale = $arguments{'yScale'}
    if ( exists($arguments{'yScale'}) );
my $yLabel = "y";
$yLabel = $arguments{'yLabel'}
    if ( exists($arguments{'yLabel'}) );
my $zScale = "linear";
$zScale = $arguments{'zScale'}
    if ( exists($arguments{'zScale'}) );
my $zLabel = "\${\\rm d}p/{\\rm d}x\$";
$zLabel = "\${\\rm d}p/{\\rm d}\\log x\$"
    if ( $xScale eq "log" );
$zLabel = $arguments{'zLabel'}
    if ( exists($arguments{'zLabel'}) );
my $plotFile = "mcmcVisualize.pdf";
$plotFile = $arguments{'output'}
    if ( exists($arguments{'output'}) );
my $title = "Marginalized posterior";
$title = $arguments{'title'}
    if ( exists($arguments{'title'}) );
my $colorbox = 1;
$colorbox = $arguments{'colorbox'}
    if ( exists($arguments{'colorbox'}) );
my $textSize = 7;
$textSize = $arguments{'textSize'}
    if ( exists($arguments{'textSize'}) );
my $labelStyle = "normalsize";
$labelStyle = $arguments{'labelStyle'}
    if ( exists($arguments{'labelStyle'}) );

# Validate.
die("mcmcVisualize.pl: if outliers are specified chainCount must also be specified")
    if ( exists($arguments{'outliers'}) && ! exists($arguments{'chainCount'}) );

# Parse label options.
my $labelX = "x";
my $labelY = "y";
if ( exists($arguments{'labels'}) ) {
    $labelX = "";
    $labelY = "";
    my @labelChoices = split(/,/,$arguments{'labels'});
    foreach my $labelChoice ( @labelChoices ) {
	$labelX = "x"
	    if ( $labelChoice eq "x" );
	$labelY = "y"
	    if ( $labelChoice eq "y" );
	$labelY = "y2"
	    if ( $labelChoice eq "y2" );
    }
}

# Open the config file and parse the available parameter names.
my $xml    = new XML::Simple;
my $config = $xml->XMLin($configFileName, KeyAttr => []);
my @properties;
foreach my $parameter ( @{$config->{'parameters'}->{'parameter'}} ) {
    push(@properties,$parameter->{'name'})
	 if ( exists($parameter->{'prior'}) );
}

# Find which columns to plot.
my $xColumn = 5;
my $yColumn = 6;
my $dimensions = 1;
if ( exists($arguments{'xProperty'}) ) {
    for(my $i=0;$i<scalar(@properties);++$i) {
	$xColumn = $i+5
	    if ( $properties[$i] eq $arguments{'xProperty'} );
    } 
 }
if ( exists($arguments{'yProperty'}) ) {
    $dimensions = 2;
    for(my $i=0;$i<scalar(@properties);++$i) {
	$yColumn = $i+5
	    if ( $properties[$i] eq $arguments{'yProperty'} );
    } 
}
if ( exists($arguments{'range'}) ) {
    foreach my $range ( @{$arguments{'range'}} ) {
	for(my $i=0;$i<scalar(@properties);++$i) {
	    $range->{'column'} = $i
		if ( $properties[$i] eq $range->{'name'} );
	}
    } 
}

# Extract only entries within a given range if necessary.
my $kdeFileName = $workDirectory."/tmp.kde";
my @outlierChains;
@outlierChains = split(/,/,$arguments{'outliers'})
    if ( exists($arguments{'outliers'}) );
open(oHndl,">".$kdeFileName);
for(my $iChain=0; -e $fileRoot."_".sprintf("%4.4d",$iChain).".log";++$iChain) {
    my $skip = 0;
    if ( scalar(@outlierChains) > 0 ) {
	foreach ( @outlierChains ) {
	    $skip = 1
		if ( $iChain == $_ );
	}
    }
    next
	if ( $skip == 1 );
    open(iHndl,$fileRoot."_".sprintf("%4.4d",$iChain).".log");
    my $i = -1;
    while ( my $line = <iHndl> ) {
	my $lineCopy = $line;
	$lineCopy =~ s/^\s*//;
	$lineCopy =~ s/\s*$//;
	my @columns = split(/\s+/,$lineCopy);
	my $include = 1;
	if ( exists($arguments{'range'}) ) {
	    foreach my $range ( @{$arguments{'range'}} ) {
		$include = 0
		    if (
			$columns[$range->{'column'}] < $range->{'lower'}
			||
			$columns[$range->{'column'}] > $range->{'upper'}
		    );
	    }
	}
	print oHndl $line
	    if ( $include == 1 );
    }
    close(iHndl);
}
close(oHndl);

# Run script to do the kernel density estimation.
my $command = "python constraints/visualization/kernelDensityEstimation.py ".$kdeFileName." ".$xColumn;
$command .= " ".$yColumn
    if ( $dimensions == 2 );
$command .= " --ngood=".$arguments{'ngood'}
    if ( exists($arguments{'ngood'}) );
$command .= " --ngrid=".$nGrid." --output=".$workDirectory."/kde.txt";
system($command);
unlink($kdeFileName)
    if ( exists($arguments{'range'}) );

# Read the file.
my $x = pdl [];
my $y = pdl [];
my $p = pdl [];
open(iHndl,$workDirectory."/kde.txt");
while ( my $line = <iHndl> ) {
    chomp($line);
    my @columns = split(/\s+/,$line);
    $x = $x->append($columns[0]);
    if ( $dimensions == 1 ) {
	$p = $p->append($columns[1]);
    } else {
	$y = $y->append($columns[1]);
	$p = $p->append($columns[2]);
    }
}
close(iHndl);

# Determine contour levels for confidence regions.
my $pSorted          = $p->qsort();
my $total            = $p->sum ();
my $pCumulant        = $pSorted->cumusumover()/$total;
my $confidence       = pdl ( 0.682689492137,  0.954499736104, 0.997300203937 );
my $excludedFraction = 1.0-$confidence;
(my $levels, my $error) = interpolate($excludedFraction,$pCumulant,$pSorted);

# Convert axes with logarithmic scaling.
$x = exp($x)
    if ( $xScale eq "log" );
$y = exp($y)
    if ( $yScale eq "log" );

# Find ranges.
my $xMin = $x->min();
my $xMax = $x->max();
my $yMin = $y->min();
my $yMax = $y->max();
my $pMin = 0.0;
$pMin = 1.0e-3
    if ( $zScale eq "log" );
my $pMax = 1.02*$p->max();
my $pMaxShort = $pMax*0.9;

# Remove temporary file.
unlink($workDirectory."/kde.txt");

# Plot the posterior.
# Make plot of redshift evolution.
my $plot;
my $gnuPlot;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid ".$textSize."\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title '".$title."'\n"
    unless ( $title eq "none" );
print $gnuPlot "set xlabel '{\\".$labelStyle." ".$xLabel."}'\n"
    unless ( $labelX eq "" );
unless ( $labelY eq "" ) {
    if ( $dimensions == 1 ) {
	print $gnuPlot "set ylabel '{\\".$labelStyle." ".$zLabel."}'\n";
    } else {
	print $gnuPlot "set ylabel '{\\".$labelStyle." ".$yLabel."}'";
	print $gnuPlot " offset graph 1.2,0"
	    if ( $labelY eq "y2" );
	print $gnuPlot "\n";
    }
}
print $gnuPlot "set ytics offset graph 1.1,0\n"
    if ( $labelY eq "y2" );
print $gnuPlot "unset colorbox\n"
    if ( $colorbox == 0 );
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
if ( $xScale eq "log" ) {
    print $gnuPlot "set logscale x\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n"
	unless ( $labelX eq "" );
}
if ( ( $dimensions == 1 && $zScale eq "log" ) || ( $dimensions == 2 && $yScale eq "log" ) ) {
    print $gnuPlot "set logscale y\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n"
	unless ( $labelY eq "" );
}
print $gnuPlot "set xrange [".$xMin.":".$xMax."]\n";
if ( $dimensions == 1 ) {
    print $gnuPlot "set yrange [".$pMin.":".$pMax."]\n";
} else {
    print $gnuPlot "set yrange [".$yMin.":".$yMax."]\n";
}
print $gnuPlot "set format x ''\n"
    if ( $labelX eq "" );
print $gnuPlot "set format y ''\n"
    if ( $labelY eq "" );
print $gnuPlot "set pointsize 2.0\n";
if ( $dimensions == 1 ) {
    for(my $i=2;$i>=0;--$i) {
	my $selection = which($p <= $levels->(($i)));
	my $pLevel = $p->copy();
	$pLevel->index($selection) .= -1.0;
	&PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $x,$pLevel,
	    style       => "filledCurve",
	    weight      => [5,3],
	    color       => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[$i]},
	    filledCurve => "x1"
	    );
    }
    &PrettyPlots::Prepare_Dataset(
	 \$plot,
	 $x,$p,
	 style       => "line",
	 weight      => [5,3],
	 color       => $PrettyPlots::colorPairs{'redYellow'},
    );
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    # Output data to file if required.
    if ( exists($arguments{'data'}) ) {
	my $xml = new XML::Simple(RootName => "data", NoAttr => 1);
	my $data;
	@{$data->{'x'         }} = $x         ->list();
	@{$data->{'p'         }} = $p         ->list();
	@{$data->{'confidence'}} = $levels    ->list();
	@{$data->{'levels'    }} = $confidence->list();
	$data  ->{'xProperty' }  = $arguments{'xProperty'};
	$data  ->{'xLabel'    }  = $xLabel;
	$data  ->{'pLabel'    }  = $zLabel;
	$data  ->{'xScale'    }  = $xScale;
	open(my $oHndl,">".$arguments{'data'});
	print $oHndl $xml->XMLout($data);
	close($oHndl);
    }
} else {
    # Generate confidence contours.
    open(ppHndl,"|gnuplot");
    print ppHndl "set table 'contour.dat'\n";
    print ppHndl "unset surface\n";
    print ppHndl "set contour base; set cntrparam levels discrete ".join(",",$levels->list())."\n";
    print ppHndl "splot '-'\n";
    my $k = -1;
    for(my $i=0;$i<$nGrid;++$i) {
	for(my $j=0;$j<$nGrid;++$j) {
	    ++$k;
	    print ppHndl $x->(($k))." ".$y->(($k))." ".$p->(($k))."\n";
	}
	print ppHndl "\n"
	    unless ( $i == $nGrid-1 );
    }
    print ppHndl "e\n";
    print ppHndl "unset table\n";
    close(ppHndl);
    system("awk \"NF<2{printf\\\"\\n\\\"}{print}\" <contour.dat >contour1.dat");
    print $gnuPlot "set pm3d map\n";
    print $gnuPlot "set pm3d explicit\n";
    print $gnuPlot "set pm3d corners2color c1\n";
    print $gnuPlot "set palette rgbformulae 21,22,23 negative\n";
    print $gnuPlot "set log cb\n"
	if ( $zScale eq "log" );
    print $gnuPlot "splot '-' with pm3d notitle, 'contour1.dat' with line lt -1 lc rgbcolor \"#FF44FF\" notitle\n";
    $k = -1;
    for(my $i=0;$i<$nGrid;++$i) {
	for(my $j=0;$j<$nGrid;++$j) {
	    ++$k;
	    print $gnuPlot $x->(($k))." ".$y->(($k))." ".$p->(($k))."\n";
	}
	print $gnuPlot "\n"
	    unless ( $i == $nGrid-1 );
    }
    print $gnuPlot "e\n";
    # Output data to file if required.
    if ( exists($arguments{'data'}) ) {
	my $xml = new XML::Simple(RootName => "data", NoAttr => 1);
	my $data;
	@{$data->{'x'         }} = $x         ->list();
	@{$data->{'y'         }} = $y         ->list();
	@{$data->{'p'         }} = $p         ->list();
	@{$data->{'confidence'}} = $levels    ->list();
	@{$data->{'levels'    }} = $confidence->list();
	$data  ->{'xProperty' }  = $arguments{'xProperty'};
	$data  ->{'yProperty' }  = $arguments{'yProperty'};
	$data  ->{'xLabel'    }  = $xLabel;
	$data  ->{'yLabel'    }  = $yLabel;
	$data  ->{'xScale'    }  = $xScale;
	$data  ->{'yScale'    }  = $yScale;
	open(my $oHndl,">".$arguments{'data'});
	print $oHndl $xml->XMLout($data);
	close($oHndl);
    }
}
close($gnuPlot);
unlink("contour.dat","contour1.dat");
&LaTeX::GnuPlot2PDF($plotFileEPS, margin => -1);

exit;
