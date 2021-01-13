#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use Cwd;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use Data::Dumper;
use Galacticus::Options;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Make plots of all analyses that were output to a Galacticus model file.
# Andrew Benson (13-February-2019)

# Get arguments and options.
die("Usage: analysesPlot.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName  = $ARGV[0];
my %options;
$options{'outputDirectory'} = $galacticusFileName =~ m/^(.*)\/[^\/]+$/ ? $1 : getcwd();
&Galacticus::Options::Parse_Options(\@ARGV,\%options);
# Open the Galacticus file.
my $galacticusFile = new PDL::IO::HDF5($galacticusFileName);
# Iterate over all analyses.
my $analysesGroup = $galacticusFile->group('analyses');
my @analyses      =  $analysesGroup->groups();
foreach my $analysisName ( @analyses ) {
    print " -> ".$analysisName."\n";    # Extract all available attributes from this analysis.
    my $analysisGroup  = $analysesGroup->group($analysisName);
    my @attributeNames = $analysisGroup->attrs();
    my $attributes =
    {
	 xAxisLabel  => "x",
	 yAxisLabel  => "y",
	 xAxisIsLog  =>  0 ,
	 yAxisIsLog  =>  1 ,
	 targetLabel => "data"
    };
    ($attributes->{$_}) = $analysisGroup->attrGet($_)
	foreach ( @attributeNames );
    # Skip cases for which we have no "type" specified.
    unless ( exists($attributes->{'type'}) ) {
	print "Warning: analysis '".$analysisName."' has no 'type' attribute, so it can not be processed.\n";
    } elsif ( $attributes->{'type'} eq "function1D" ) {
    	# Simple 1D function - will be shown as an x-y scatter plot.
	# Validate attributes.
	my @datasetNames =
	    (
	     {name => 'xDataset'         , required => 1},
	     {name => 'yDataset'         , required => 1},
	     {name => 'yDatasetTarget'   , required => 0},
	     {name => 'yCovariance'      , required => 0},
	     {name => 'yCovarianceTarget', required => 0},
	     {name => 'yErrorLower'      , required => 0},
	     {name => 'yErrorUpper'      , required => 0},
	     {name => 'yErrorLowerTarget', required => 0},
	     {name => 'yErrorUpperTarget', required => 0},
	    );
	foreach my $attribute ( @datasetNames ) {
	    die("Error: attribute '".$attribute->{'name'}."' is missing from analysis '".$analysisName."' but is required.")
		unless ( ! $attribute->{'required'} || grep {$_ eq $attribute->{'name'}} keys(%{$attributes}) );
	}
	# Read the datasets.
	my $data;
	foreach my $dataset ( @datasetNames ) {
	    if ( grep {$_ eq $dataset->{'name'}} keys(%{$attributes}) ) {
		(my $analysisDatasetName) = $analysisGroup->attrGet($dataset->{'name'});
		$data->{$dataset->{'name'}} = $analysisGroup->dataset($analysisDatasetName)->get();
	    }
	}
	&function1DPlot($data,$attributes,$analysisName);
    } elsif ( $attributes->{'type'} eq "function1DSequence" ) {
    	# Sequence of simple 1D functions - will be shown as multiple x-y scatter plot.
	# Validate attributes.
	my @datasetNames =
	    (
	     {name => 'xDataset'         , required => 1},
	     {name => 'yDataset'         , required => 1},
	     {name => 'yDatasetTarget'   , required => 1},
	     {name => 'yCovariance'      , required => 0},
	     {name => 'yCovarianceTarget', required => 0},
	     {name => 'yErrorLower'      , required => 0},
	     {name => 'yErrorUpper'      , required => 0},
	     {name => 'yErrorLowerTarget', required => 0},
	     {name => 'yErrorUpperTarget', required => 0},
	    );
	foreach my $attribute ( @datasetNames ) {
	    die("Error: attribute '".$attribute->{'name'}."' is missing from analysis '".$analysisName."' but is required.")
		unless ( ! $attribute->{'required'} || grep {$_ eq $attribute->{'name'}} keys(%{$attributes}) );
	}
	# Read the datasets.
	my $data;
	foreach my $dataset ( @datasetNames ) {
	    if ( grep {$_ eq $dataset->{'name'}} keys(%{$attributes}) ) {
		(my $analysisDatasetName) = $analysisGroup->attrGet($dataset->{'name'});
		$data->{$dataset->{'name'}} = $analysisGroup->dataset($analysisDatasetName)->get();
	    }
	}
	# Iterate over each individual plot in the sequence.
	for(my $i=0;$i<$data->{'yDataset'}->dim(1);++$i) {
	    my $dataSequence                              ;
	    my $size1D       = $data->{'xDataset'}->dim(0);
	    $dataSequence->{'xDataset'} = $data->{'xDataset'};
	    foreach my $dataset ( 'yDataset', 'yDatasetTarget', 'yErrorLower', 'yErrorUpper', 'yErrorLowerTarget', 'yErrorUpperTarget' ) {
		$dataSequence->{$dataset} = $data->{$dataset}->(:,($i))
		    if ( exists($data->{$dataset}) );
	    }
	    foreach my $dataset ( 'yCovariance', 'yCovarianceTarget' ) {
		$dataSequence->{$dataset} = $data->{$dataset}->($size1D*$i:$size1D*($i+1)-1,$size1D*$i:$size1D*($i+1)-1)
		    if ( exists($data->{$dataset}) );
	    }
	    &function1DPlot($dataSequence,$attributes,$analysisName.":".$i);
	}
    }
}

exit 0;

sub function1DPlot {
    # Plot a 1-D function.
    my $data         = shift();
    my $attributes   = shift();
    my $analysisName = shift();
    # Fill in missing target dataset.
    my $haveTarget   = exists($data->{'yDatasetTarget'});
    $data->{'yDatasetTarget'} = $data->{'yDataset'}->copy()
	unless ( $haveTarget );
    # Determine non-zero entries.
    my $nonZero       = which( $data->{'yDataset'      } > 0.0)                                     ;
    my $nonZeroTarget = which(                                     $data->{'yDatasetTarget'} > 0.0) ;
    my $nonZeroBoth   = which(($data->{'yDataset'      } > 0.0) | ($data->{'yDatasetTarget'} > 0.0));
    # Determine plot ranges.
    my $yErrorLower      ;
    my $yErrorUpper      ;
    my $yErrorLowerTarget;
    my $yErrorUpperTarget;
    if ( exists($data->{'yErrorLower'      }) ) {
	$yErrorLower       = $data->{'yErrorLower'      }                       ;
    }
    if ( exists($data->{'yErrorUpper'      }) ) {
	$yErrorUpper       = $data->{'yErrorUpper'      }                       ;
    }
    if ( exists($data->{'yCovariance'      }) ) {
	$yErrorUpper       = $data->{'yCovariance'      }->diagonal(0,1)->sqrt();
	$yErrorLower       = $data->{'yCovariance'      }->diagonal(0,1)->sqrt();
    }
    if ( exists($data->{'yErrorLowerTarget'}) ) {
	$yErrorLowerTarget = $data->{'yErrorLowerTarget'}                       ;
    }
    if ( exists($data->{'yErrorUpperTarget'}) ) {
	$yErrorUpperTarget = $data->{'yErrorUpperTarget'}                       ;
    }
    if ( exists($data->{'yCovarianceTarget'}) ) {
	$yErrorUpperTarget = $data->{'yCovarianceTarget'}->diagonal(0,1)->sqrt();
	$yErrorLowerTarget = $data->{'yCovarianceTarget'}->diagonal(0,1)->sqrt();
    }
    my $yLower       = $data->{'yDataset'      }->copy();
    my $yUpper       = $data->{'yDataset'      }->copy();
    my $yLowerTarget = $data->{'yDatasetTarget'}->copy();
    my $yUpperTarget = $data->{'yDatasetTarget'}->copy();
    $yLower       -= $yErrorLower       if ( defined($yErrorLower      ) );
    $yUpper       += $yErrorUpper       if ( defined($yErrorUpper      ) );
    $yLowerTarget -= $yErrorLowerTarget if ( defined($yErrorLowerTarget) );
    $yUpperTarget += $yErrorUpperTarget if ( defined($yErrorUpperTarget) );
    if ( $attributes->{'yAxisIsLog'} ) {
	my $negativeY                      = which($yLower       <= 1.0e-10*$data->{'yDataset'      });
	my $negativeYTarget                = which($yLowerTarget <= 1.0e-10*$data->{'yDatasetTarget'});
	$yLower      ->($negativeY      ) .= $data->{'yDataset'      }->($negativeY      );
	$yLowerTarget->($negativeYTarget) .= $data->{'yDatasetTarget'}->($negativeYTarget);
    }
    my $yUpperBoth   = $yUpper->($nonZero)->append($yUpperTarget->($nonZeroTarget));
    my $yLowerBoth   = $yLower->($nonZero)->append($yLowerTarget->($nonZeroTarget));
    my $xMinimum = $data->{'xDataset'}->($nonZeroBoth)->minimum();
    my $xMaximum = $data->{'xDataset'}->($nonZeroBoth)->maximum();
    my $yMinimum = $yLowerBoth                        ->minimum();
    my $yMaximum = $yUpperBoth                        ->maximum();
    if ( $attributes->{'xAxisIsLog'} ) {
	$xMinimum /= 1.05;
	$xMaximum *= 1.05;
	# Ensure that we span at least one integer.
	my $xMinimumLog = log10($xMinimum);
	my $xMaximumLog = log10($xMaximum);
	if ( int($xMaximumLog) == int($xMinimumLog) ) {
	    if ( $xMinimumLog-int($xMinimumLog) < int($xMaximumLog)+1-$xMaximumLog ) {
		$xMinimum = 10.0** int($xMinimumLog)   ;
	    } else {
		$xMaximum = 10.0**(int($xMaximumLog)+1);
	    }
	}		
    } else {
	my $xRange  = $xMaximum-$xMinimum;
	$xRange    .= $xMinimum
	    if ( $xRange == 0.0 );
	$xMinimum  -= 0.05*$xRange;
	$xMaximum  += 0.05*$xRange;
    }
    if ( $attributes->{'yAxisIsLog'} ) {
	$yMinimum /= 1.05;
	$yMaximum *= 1.05;
	# Ensure that we span at least one integer.
	my $yMinimumLog = log10($yMinimum);
	my $yMaximumLog = log10($yMaximum);
	if ( int($yMaximumLog) == int($yMinimumLog) ) {
	    if ( $yMinimumLog-int($yMinimumLog) < int($yMaximumLog)+1-$yMaximumLog ) {
		$yMinimum = 10.0** int($yMinimumLog)   ;
	    } else {
		$yMaximum = 10.0**(int($yMaximumLog)+1);
	    }
	}		
    } else {
	my $yRange  = $yMaximum-$yMinimum;
	$yRange    .= $yMinimum
	    if ( $yRange == 0.0 );
	$yMinimum  -= 0.05*$yRange;
	$yMaximum  += 0.05*$yRange;
    }
    # Override axes.
    $xMinimum .= $options{'xMinimum'}
        if ( exists($options{'xMinimum'}) );
    $xMaximum .= $options{'xMaximum'}
        if ( exists($options{'xMaximum'}) );
    $yMinimum .= $options{'yMinimum'}
        if ( exists($options{'yMinimum'}) );
    $yMaximum .= $options{'yMaximum'}
        if ( exists($options{'yMaximum'}) );
    # Determine a suitable pointsize.
    my $pointSize = "1.0";
    $pointSize = 0.25
	if ( nelem($data->{'xDataset'}) > 100 );
    # Construct a plot.
    my $plot;
    my $gnuPlot;
    my $plotFileTeX = $options{'outputDirectory'}."/".$analysisName.".tex";
    open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set title offset 0,-0.8 '\\tiny ".$attributes->{'description'}."'\n"
	if ( exists($attributes->{'description'}) );
    print $gnuPlot "set xlabel '".$attributes->{'xAxisLabel'}."'\n";
    print $gnuPlot "set ylabel '".$attributes->{'yAxisLabel'}."'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.20\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.1\n";
    print $gnuPlot "set key bmargin box\n";
    foreach my $axis ( 'x', 'y' ) {
	if ( $attributes->{$axis.'AxisIsLog'} ) {
	    print $gnuPlot "set logscale ".$axis."\n";
	    print $gnuPlot "set m".$axis."tics 10\n";
	    print $gnuPlot "set format ".$axis." '\$10^{\%L}\$'\n";
	}
    }
    print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
    print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
    print $gnuPlot "set pointsize ".$pointSize."\n";
    (my $targetLabel = $attributes->{'targetLabel'}) =~ s/&/\\&/g;
    my %plotOptionsTarget =
	(
	 style      => "point",
	 symbol     => [6,7],
	 weight     => [2,1],
	 color      => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	 title      => "\\\\footnotesize ".$targetLabel
	);
    $plotOptionsTarget{'errorDown'} = $yErrorLowerTarget->($nonZeroTarget)
	if ( defined($yErrorLowerTarget) );
    $plotOptionsTarget{'errorUp'  } = $yErrorUpperTarget->($nonZeroTarget)
	if ( defined($yErrorUpperTarget) );
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $data->{'xDataset'      }->($nonZeroTarget),
	    $data->{'yDatasetTarget'}->($nonZeroTarget),
	    %plotOptionsTarget
	)
	if ( $haveTarget );
    my %plotOptions =
	(
	 style      => "point",
	 symbol     => [6,7],
	 weight     => [2,1],
	 color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	 title      => "\\\\footnotesize Galacticus"
	);
    $plotOptions{'errorDown'} = $yErrorLower->($nonZero)
	if ( defined($yErrorLower) );
    $plotOptions{'errorUp'  } = $yErrorUpper->($nonZero)
	if ( defined($yErrorUpper) );
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $data->{'xDataset'}->($nonZero),
	    $data->{'yDataset'}->($nonZero),
	    %plotOptions
	);
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
}
