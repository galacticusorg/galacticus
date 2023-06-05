#!/usr/bin/env perl
use Cwd;
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}.'/perl';
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
use XML::Simple;
use List::Uniq ':all';
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Visualize the posterior distribution from an MCMC simulation log file.
# Andrew Benson (10-June-2012)

# Get file name to process.
die("Usage: mcmcVisualize.pl parameterFileName fileRoot1 fileRoot2.... [options]")
    unless ( scalar(@ARGV) > 1 );
my $parameterFileName = $ARGV[0];
my @fileRoots;
for(my $i=1;$i<scalar(@ARGV);++$i) {
    last
	if ( $ARGV[$i] =~ m/^\-\-/ );
    push(@fileRoots,$ARGV[$i]);
}

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
my $xScale = "identity";
$xScale = $arguments{'xScale'}
if ( exists($arguments{'xScale'}) );
die("xScale must be 'identity' or 'logarithm'")
    unless ( $xScale eq "identity" || $xScale eq "logarithm" );
my $xLabel = "x";
$xLabel = $arguments{'xLabel'}
if ( exists($arguments{'xLabel'}) );
my $yScale = "identity";
$yScale = $arguments{'yScale'}
    if ( exists($arguments{'yScale'}) );
die("yScale must be 'identity' or 'logarithm'")
    unless ( $yScale eq "identity" || $yScale eq "logarithm" );
my $yLabel = "y";
$yLabel = $arguments{'yLabel'}
    if ( exists($arguments{'yLabel'}) );
my $zScale = "identity";
$zScale = $arguments{'zScale'}
if ( exists($arguments{'zScale'}) );
die("zScale must be 'identity' or 'logarithm'")
    unless ( $zScale eq "identity" || $zScale eq "logarithm" );
my $zLabel = "\${\\rm d}p/{\\rm d}x\$";
$zLabel = "\${\\rm d}p/{\\rm d}\\log x\$"
    if ( $xScale eq "logarithm" );
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
my $plotSize = "5cm,5cm";
$plotSize = $arguments{'plotSize'}
    if ( exists($arguments{'plotSize'}) );
my $textSize = 7;
$textSize = $arguments{'textSize'}
    if ( exists($arguments{'textSize'}) );
my $labelStyle = "normalsize";
$labelStyle = $arguments{'labelStyle'}
    if ( exists($arguments{'labelStyle'}) );
my $showLabels = "yes";
$showLabels = $arguments{'showLabels'}
    if ( exists($arguments{'showLabels'}) );
my $lineWeight = "5,3";
$lineWeight = $arguments{'lineWeight'}
    if ( exists($arguments{'lineWeight'}) );
my $firstParameterColumn = exists($arguments{'oldChainFormat'}) && $arguments{'oldChainFormat'} eq "yes" ? 5 : 6;

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
my $xml        = new XML::Simple;
my $parameters = $xml->XMLin($parameterFileName);
my @properties;
foreach my $parameter ( @{$parameters->{'posteriorSampleSimulation'}->{'modelParameter'}} ) {
    push(@properties,$parameter->{'name'}->{'value'})
	 if ( $parameter->{'value'} eq "active" );
}

# Find which columns to plot.
my $xColumn;
my $yColumn;
my $dimensions = 1;
if ( exists($arguments{'xProperty'}) ) {
    for(my $i=0;$i<scalar(@properties);++$i) {
	print $i."\t". $properties[$i]."\n";
	$xColumn = $i+$firstParameterColumn
	    if ( $properties[$i] eq $arguments{'xProperty'} );
    } 
    die("x-parameter '".$arguments{'xProperty'}."' not found")
	unless ( defined($xColumn) )
} else {
    $xColumn = $firstParameterColumn;
}
if ( exists($arguments{'yProperty'}) ) {
    $dimensions = 2;
    for(my $i=0;$i<scalar(@properties);++$i) {
	$yColumn = $i+$firstParameterColumn
	    if ( $properties[$i] eq $arguments{'yProperty'} );
    } 
    die("y-parameter '".$arguments{'yProperty'}."' not found")
	unless ( defined($yColumn) )
} else {
    $yColumn = $firstParameterColumn+1;
}

if ( exists($arguments{'range'}) ) {
    foreach my $range ( @{$arguments{'range'}} ) {
	for(my $i=0;$i<scalar(@properties);++$i) {
	    $range->{'column'} = $i+$firstParameterColumn
		if ( $properties[$i] eq $range->{'name'} );
	}
    } 0
}

# Colors for contours.
my @contourColors = ( "FF44FF", "AAAAAA" );

# Plot the posterior.
my $plot;
my $gnuPlot;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid size ".$plotSize." ".$textSize."\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title '".$title."'\n"
    unless ( $title eq "none" );
print $gnuPlot "set xlabel '{\\".$labelStyle." ".$xLabel."}'\n"
    unless ( $labelX eq "" || $showLabels eq "no" );
unless ( $labelY eq "" || $showLabels eq "no" ) {
    if ( $dimensions == 1 ) {
	print $gnuPlot "set ylabel '{\\".$labelStyle." ".$zLabel."}'\n";
    } else {
	print $gnuPlot "set ylabel '{\\".$labelStyle." ".$yLabel."}'";
	print $gnuPlot " offset graph 1.4,0"
	    if ( $labelY eq "y2" );
	print $gnuPlot "\n";
    }
}
print $gnuPlot "set ytics offset graph 1.2,0\n"
    if ( $labelY eq "y2" && $showLabels eq "yes");
print $gnuPlot "unset colorbox\n"
    if ( $colorbox == 0 );
if ( $showLabels eq "yes" ) {
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
} else {
    print $gnuPlot "set lmargin screen 0\n";
    print $gnuPlot "set rmargin screen 0.99\n";
    print $gnuPlot "set bmargin screen 0\n";
    print $gnuPlot "set tmargin screen 0.99\n";
}
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
if ( $xScale eq "logarithm" ) {
    print $gnuPlot "set logscale x\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n"
	unless ( $labelX eq "" );
}
if ( ( $dimensions == 1 && $zScale eq "logarithm" ) || ( $dimensions == 2 && $yScale eq "logarithm" ) ) {
    print $gnuPlot "set logscale y\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n"
	unless ( $labelY eq "" );
}
print $gnuPlot "set format x ''\n"
    if ( $labelX eq "" );
print $gnuPlot "set format y ''\n"
    if ( $labelY eq "" );
if ( $showLabels eq "no" ) {
    print $gnuPlot "set xtics format ''\n";
    print $gnuPlot "set ytics format ''\n";
}
print $gnuPlot "set pointsize 2.0\n";

# Temporary files.
my @tempFiles;

print $gnuPlot "set multiplot\n"
    if ( $dimensions == 2 );

# Iterate over chain sets.
my $iFileRoot = -1;
foreach my $fileRoot ( @fileRoots ) {
    ++$iFileRoot;
    
    # Extract only converged entries and entries within a given range if necessary.
    my $kdeFileName = $plotFile.".tmp.kde";
    push(@tempFiles,$kdeFileName);
    my @outlierChains;
    @outlierChains = split(/,/,$arguments{'outliers'})
	if ( exists($arguments{'outliers'}) );
    # If sampling from end of chains, determine chain lengths.
    my $iChainBegin = 0;
    if ( exists($arguments{'ngood'}) ) {
	my $chainLength = 0;
	open(iHndl,$fileRoot."_0000.log");
	while ( my $line = <iHndl> ) {
	    ++$chainLength;
	}
	close(iHndl);
	$iChainBegin = $chainLength+1-$arguments{'ngood'};
    }
    # Initialize maximum likelihood state.
    my $likelihoodMaximum           = pdl -1.0e60;
    my @likelihoodMaximumParameters;
    # Read chains.
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
	    ++$i;
	    my $lineCopy = $line;
	    $lineCopy =~ s/^\s*//;
	    $lineCopy =~ s/\s*$//;
	    my @columns = split(/\s+/,$lineCopy);
	    my $include = 1;
	    $include = 0
		if ( $columns[3] eq "F" && ( ! exists($arguments{'useUnconverged'}) || $arguments{'useUnconverged'} eq "no" ) );
	    $include = 0
		if ( exists($arguments{'minimumLogL'}) && $columns[4] < $arguments{'minimumLogL'} );
	    $include = 0
		if ( $i < $iChainBegin );
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
	    # Proceed only if this state is to be included.
	    if ( $include == 1 ) {
		# Check for maximum likelihood parameters.
		if ( $columns[4] > $likelihoodMaximum ) {
		    $likelihoodMaximum = $columns[4];
		    @likelihoodMaximumParameters = @columns;
		}
		# Take the log if necessary.
		$columns[$xColumn] = log($columns[$xColumn])
		    if ( $xScale eq "logarithm" );
		$columns[$yColumn] = log($columns[$yColumn])
		    if ( $yScale eq "logarithm" );
		# Output the state.
		print oHndl join("\t",@columns)."\n";
	    }
	}
	close(iHndl);
    }
    close(oHndl);
    
    # Run script to do the kernel density estimation.
    my $command = "python constraints/visualization/kernelDensityEstimation.py ".$kdeFileName." ".$xColumn;
    $command .= " ".$yColumn
	if ( $dimensions == 2 );
    $command .= " --ngrid=".$nGrid." --output=".$plotFile.".kde";
    push(@tempFiles,$plotFile.".kde");
    system($command);

    # Read the file.
    my $x = pdl [];
    my $y = pdl [];
    my $p = pdl [];
    open(iHndl,$plotFile.".kde");
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

    # Find mode of the distribution.
    my $modeX = $x->(($p->maximum_ind()))->sclr();
    
    # Determine contour levels for confidence regions.
    my $pSorted          = $p->qsort();
    my $total            = $p->sum  ();
    my $pCumulant        = $pSorted->cumusumover()/$total;
    my $confidence       = pdl ( 0.682689492137,  0.954499736104, 0.997300203937 );
    my $excludedFraction = 1.0-$confidence;
    (my $levels, my $error) = interpolate($excludedFraction,$pCumulant,$pSorted);

    # Find span of confidence regions.
    my @spanLowerX;
    my @spanLowerY;
    my @spanUpperX;
    my @spanUpperY;
    for(my $i=0;$i<nelem($confidence);++$i) {
	my $selected = which($p > $levels->(($i)));
	if ( nelem($selected) > 0 ) {
	    my $indexX   = $x->($selected)->qsorti();
	    my $indexY   = $y->($selected)->qsorti()
		if ( $dimensions == 2 );
	    my $lowerX = $x->($selected)->($indexX)->(( 0))->copy();
	    my $upperX = $x->($selected)->($indexX)->((-1))->copy();
	    if ( $xScale eq "logarithm" ) {
		$lowerX .= exp($lowerX);
		$upperX .= exp($upperX);
	    }
	    push(@spanLowerX,$lowerX->sclr());
	    push(@spanUpperX,$upperX->sclr());
	    if ( $dimensions == 2 ) {
		my $lowerY = $y->($selected)->($indexY)->(( 0))->copy();
		my $upperY = $y->($selected)->($indexY)->((-1))->copy();
		if ( $yScale eq "logarithm" ) {
		    $lowerY .= exp($lowerY);
		    $upperY .= exp($upperY);
		}
		push(@spanLowerY,$lowerY->sclr());
		push(@spanUpperY,$upperY->sclr());
	    }
	} else {
	    push(@spanLowerX,"undef");
	    push(@spanLowerY,"undef");
	    push(@spanUpperX,"undef");
	    push(@spanUpperY,"undef");
	}
    }

    # Convert axes with logarithmic scaling.
    $x = exp($x)
	if ( $xScale eq "logarithm" );
    $y = exp($y)
	if ( $yScale eq "logarithm" );

    # Find ranges.
    my $pMin;
    my $pMax;
    if ( $iFileRoot == 0 ) {
	my $xMin = $x->min();
	my $xMax = $x->max();
	if ( exists($arguments{'xRange'}) && $arguments{'xRange'} =~ m/(.*):(.*)/ ) {
	    $xMin = $1;
	    $xMax = $2;
	}
	my $yMin = $y->min();
	my $yMax = $y->max();
	if ( exists($arguments{'yRange'}) && $arguments{'yRange'} =~ m/(.*):(.*)/ ) {
	    $yMin = $1;
	    $yMax = $2;
	}
	$pMin = 0.0;
	$pMin = 1.0e-3
	    if ( $zScale eq "logarithm" );
	$pMax = 1.02*$p->max();
	if ( exists($arguments{'pRange'}) && $arguments{'pRange'} =~ m/(.*):(.*)/ ) {
	    $pMin = $1;
	    $pMax = $2;
	}
	print $gnuPlot "set xrange [".$xMin.":".$xMax."]\n";
	if ( $dimensions == 1 ) {
	    print $gnuPlot "set yrange [".$pMin.":".$pMax."]\n";
	} else {
	    print $gnuPlot "set yrange [".$yMin.":".$yMax."]\n";
	}
    }

    # Remove temporary file.
    unlink($workDirectory."/kde.txt");

    if ( $dimensions == 1 ) {
	my @lineWeight = split(",",$lineWeight);
	for(my $i=2;$i>=0;--$i) {
	    my $selection = which($p <= $levels->(($i)));
	    my $pLevel = $p->copy();
	    $pLevel->index($selection) .= -1.0;
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot,
		$x,$pLevel,
		style       => "filledCurve",
		weight      => \@lineWeight,
		color       => $GnuPlot::PrettyPlots::colorPairs{${$GnuPlot::PrettyPlots::colorPairSequences{'sequence1'}}[$i]},
		filledCurve => "x1"
		);
	}
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $x,$p,
	    style       => "line",
	    weight      => \@lineWeight,
	    color       => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	    );
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	# Output data to file if required.
	if ( exists($arguments{'data'}) ) {
	    my $xml = new XML::Simple(RootName => "data", NoAttr => 1);
	    my $data;
	    @{$data->{'x'                      }} = $x         ->list();
	    @{$data->{'p'                      }} = $p         ->list();
	    @{$data->{'confidence'             }} = $levels    ->list();
	    @{$data->{'levels'                 }} = $confidence->list();
	    @{$data->{'spanLowerX'             }} = @spanLowerX;
	    @{$data->{'spanUpperX'             }} = @spanUpperX;
	    $data  ->{'modeValueX'             }  = $modeX;
	    $data  ->{'maximumLikelihoodValueX'}  = $likelihoodMaximumParameters[$xColumn];
	    $data  ->{'xProperty'              }  = $arguments{'xProperty'};
	    $data  ->{'xLabel'                 }  = $xLabel;
	    $data  ->{'pLabel'                 }  = $zLabel;
	    $data  ->{'xScale'                 }  = $xScale;
	    open(my $oHndl,">".$arguments{'data'});
	    print $oHndl $xml->XMLout($data);
	    close($oHndl);
	}
    } else {
	# Generate confidence contours.
	open(ppHndl,"|gnuplot");
	print ppHndl "set table '".$plotFile.".contour.dat'\n";
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
	push(@tempFiles,$plotFile.".contour.dat",$plotFile.".contour1.dat");
	system("awk \"NF<2{printf\\\"\\n\\\"}{print}\" <".$plotFile.".contour.dat >".$plotFile.".contour1.dat");
	print $gnuPlot "set pm3d map\n";
	print $gnuPlot "set pm3d explicit\n";
	print $gnuPlot "set pm3d corners2color c1\n";
	print $gnuPlot "set palette rgbformulae 21,22,23 negative\n";
	print $gnuPlot "set log cb\n"
	    if ( $zScale eq "logarithm" );
	if ( $iFileRoot == 0 ) {
	    print $gnuPlot "set cbrange [".$pMin.":".$pMax."]\n";
	} elsif ( $iFileRoot == 1 ) {
	    print $gnuPlot "unset colorbox\n"
		unless ( $colorbox == 0 );
	}
	if ( $iFileRoot == 0 ) {
	    print $gnuPlot "splot '-' with pm3d notitle, '".$plotFile.".contour1.dat' with line lt -1 lc rgbcolor \"#FF44FF\" notitle\n";
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
 	} else {
	    print $gnuPlot "splot 'contour".$iFileRoot.".dat' with line lt -1 lc rgbcolor \"#".$contourColors[$iFileRoot]."\" notitle\n";
	}
	# Output data to file if required.
	if ( exists($arguments{'data'}) ) {
	    my $xml = new XML::Simple(RootName => "data", NoAttr => 1);
	    my $data;
	    @{$data->{'x'                      }} = $x         ->list();
	    @{$data->{'y'                      }} = $y         ->list();
	    @{$data->{'p'                      }} = $p         ->list();
	    @{$data->{'confidence'             }} = $levels    ->list();
	    @{$data->{'levels'                 }} = $confidence->list();	    
	    @{$data->{'spanLowerX'             }} = @spanLowerX;
	    @{$data->{'spanUpperX'             }} = @spanUpperX;
	    @{$data->{'spanLowerY'             }} = @spanLowerY;
	    @{$data->{'spanUpperY'             }} = @spanUpperY;
	    $data  ->{'maximumLikelihoodValueX'}  = $likelihoodMaximumParameters[$xColumn];
	    $data  ->{'maximumLikelihoodValueY'}  = $likelihoodMaximumParameters[$yColumn];
	    $data  ->{'xProperty'              }  = $arguments{'xProperty'};
	    $data  ->{'yProperty'              }  = $arguments{'yProperty'};
	    $data  ->{'xLabel'                 }  = $xLabel;
	    $data  ->{'yLabel'                 }  = $yLabel;
	    $data  ->{'xScale'                 }  = $xScale;
	    $data  ->{'yScale'                 }  = $yScale;
	    open(my $oHndl,">".$arguments{'data'});
	    print $oHndl $xml->XMLout($data);
	    close($oHndl);
	}
    }
}

print $gnuPlot "unset multiplot\n"
    if ( $dimensions == 2 );

close($gnuPlot);
unlink(uniq(sort(@tempFiles)));
&GnuPlot::LaTeX::GnuPlot2PDF($plotFileEPS, margin => 1);

exit;
