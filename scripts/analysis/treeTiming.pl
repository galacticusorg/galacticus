#!/usr/bin/env perl
use strict;
use warnings;
use lib "./perl";
use PDL;
use PDL::IO::HDF5;
use GnuPlot::LaTeX;
use GnuPlot::PrettyPlots;
use Data::Dumper;
use List::Util qw(first);
use PDL::Fit::Polynomial;
use PDL::NiceSlice;
use XML::Simple;
use Stats::Percentiles;

# Analyze tree timing data and construct a suitable fitting function for tree processing times.
# Andrew Benson (02-September-2011)

# Get arguments.
die ("Usage: treeTiming.pl <modelFile> [options.....]") unless ( $#ARGV >= 0 );
my $modelFile = $ARGV[0];

# Create a hash of named arguments.
my $iArg = -1;
my %arguments;
while ( $iArg < scalar(@ARGV)-1 ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	my $argumentName = $1;
	if ( $iArg == scalar(@ARGV)-1 ) {
	    $arguments{$argumentName} = 0;
	} else {
	    if ( $ARGV[$iArg+1] =~ m/^\-\-(.*)/ ) {
		$arguments{$argumentName} = 0;
	    } else {
		$arguments{$argumentName} = $ARGV[$iArg+1];
		++$iArg;
	    }
	}
    }
}

# Initalize datasets.
my $treeMasses = pdl [];
my $treeTimes  = pdl [];

# Read existing data if requested.
if ( exists($arguments{'accumulate'}) ) {
    die("treeTiming.xml: an output file must be given when accumulating") unless ( exists($arguments{'outputFile'}) );
    if ( -e $arguments{'outputFile'} ) {
	my $xml = new XML::Simple;
	my $currentData = $xml->XMLin($arguments{'outputFile'});
	foreach my $tree ( @{$currentData->{'tree'}} ) {
	    $treeMasses = $treeMasses->append($tree->{'mass'});
	    $treeTimes  = $treeTimes ->append($tree->{'time'});
	}
    }
}

# Open the file containing the meta-data.
my $HDFfile = new PDL::IO::HDF5($modelFile);

# Read the timing data.
my @groupsAvailable;
@groupsAvailable = $HDFfile->groups();
die ("treeTiming.pl: metaData group does not exist") unless ( defined(first { $_ eq "metaData" } @groupsAvailable) );
@groupsAvailable = $HDFfile->group("metaData")->groups();
die ("treeTiming.pl: treeTiming group does not exist") unless ( defined(first { $_ eq "treeTiming" } @groupsAvailable) );
$treeMasses            = $treeMasses->append($HDFfile->dataset("metaData/treeTiming/treeMasses"        )->get());
my $treeConstructTimes =                     $HDFfile->dataset("metaData/treeTiming/treeConstructTimes")->get() ;
my $treeEvolveTimes    =                     $HDFfile->dataset("metaData/treeTiming/treeEvolveTimes"   )->get() ;
$treeTimes             = $treeTimes ->append($treeConstructTimes+$treeEvolveTimes                              );

# Get logarithmic quantities.
my $logTreeMasses = log10($treeMasses);
my $logTreeTimes  = log10($treeTimes );

# Get the median tree times.
my $percentiles    = pdl [ 50.0 ];
my $logTreeMassBin = sequence(40)*0.5;
my $weight         = pdl ones(nelem($logTreeMasses));
my $medians        = &Percentiles::BinnedPercentiles($logTreeMassBin,$logTreeMasses,$logTreeTimes,$weight,$percentiles);
my $mediansFlat    = $medians(:,(0));
my $nonzero        = which($mediansFlat != 0.0);

# Find a polynomial fit to the timing data.
(my $yfit, my $coeffs) = fitpoly1d $logTreeMassBin->index($nonzero), $mediansFlat->index($nonzero), 3;

# Create the fit.
my $logTreeMass = sequence(700)/100+9;
my $treeMass    = 10.0**$logTreeMass;
my $treeTime    = 10.0**(
    $coeffs->index(0)+
    $coeffs->index(1)*$logTreeMass+
    $coeffs->index(2)*$logTreeMass**2
    );

# Write the data to file - keeping up to maxPoints datapoints drawn from the most recently added.
if ( exists($arguments{'outputFile'}) ) {
    my $startIndex;
    if ( exists($arguments{'maxPoints'}) ) {
	$startIndex = nelem($treeMasses)-$arguments{'maxPoints'};
	$startIndex = 0 if ( $startIndex < 0 );
    } else {
	$startIndex = 0;
    }
    my $outputData;
    for(my $i=$startIndex;$i<nelem($treeMasses);++$i) {
	push(
	    @{$outputData->{'tree'}},
	    {
		mass => $treeMasses->index($i)->list(),
		time => $treeTimes ->index($i)->list()
	    }
	    );
    }
    $outputData->{'fit'} = {
	coefficient => [$coeffs->list()]
    };
    system("mkdir -p `dirname ".$arguments{'outputFile'}."`");
    my $xmlOut = new XML::Simple(RootName=>"timing",NoAttr => 1);
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOut->XMLout($outputData);
    close(oHndl);
}

# Make a plot of tree timing data.
if ( exists($arguments{'plotFile'}) ) {
    my ($gnuPlot, $outputFile, $outputFileEPS, $plot);
    $outputFile = $arguments{'plotFile'};
    ($outputFileEPS = $outputFile) =~ s/\.pdf$/.eps/;
    open($gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
    print $gnuPlot "set output '".$outputFileEPS."'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.2,0.2\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [1.0e9:1.0e16]\n";
    print $gnuPlot "set yrange [1.0e-1:1.0e5]\n";
    print $gnuPlot "set title 'Tree processing times'\n";
    print $gnuPlot "set xlabel 'Tree mass \$[M_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Processing time [s]'\n";
    &PrettyPlots::Prepare_Dataset(\$plot,
				  $treeMasses,$treeTimes,
				  style => "point", symbol => [6,7], weight => [5,3],
				  color => $PrettyPlots::colorPairs{'cornflowerBlue'});
    &PrettyPlots::Prepare_Dataset(\$plot,
				  $treeMass,$treeTime,
				  style => "line", weight => [5,3],
				  color => $PrettyPlots::colorPairs{'redYellow'});
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($outputFileEPS);
}

exit;
