#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
 $ENV{"GALACTICUS_ROOT_V093"} = "/";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
require Stats::Means;
require Stats::Histograms;

# Test Perl modules.
# Andrew Benson (7-July-2013)

# Check individual compilation.
system("find . -type f -name \"*.pl\" | xargs -n 1 perl -c 1>perl.tmp 2>&1; sed -r s/failed/FAILED/g perl.tmp; rm perl.tmp");
system("find perl -type f -name \"*.pm\" | xargs -n 1 perl -c 1>perl.tmp 2>&1; sed -r s/failed/FAILED/g perl.tmp; rm perl.tmp");

# Statistics
#  Create bins.
my $bins = pdl sequence(3);

## Stats::Means
my $xm = pdl ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2 );
my $ym = pdl ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 );
my $wm = pdl ( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1,  1,  1,  1,   1,  1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 );
my $eps      = pdl 1.0e-3;
my $mExpect  = pdl ( 4.5       , 14.5       , 24.5        );
my $meExpect = pdl ( 0.90829511,  0.90829511,  0.90829511 );
my $dExpect  = pdl ( 3.0276503 ,  3.0276503 ,  3.0276503  );
my $deExpect = pdl ( 1.6371658 ,  1.6371658 ,  1.6371658  );
(my $m, my $me, my $d, my $de) = &Means::BinnedMean($bins,$xm,$ym,$wm);
print "FAILED: Stats::Means::BinnedMean fails to compute correct means\n"
    unless ( all(abs($m - $mExpect) < $eps*$mExpect ) );
print "FAILED: Stats::Means::BinnedMean fails to compute correct errors on means\n"
    unless ( all(abs($me-$meExpect) < $eps*$meExpect) );
print "FAILED: Stats::Means::BinnedMean fails to compute correct standard deviation\n"
    unless ( all(abs($d - $dExpect) < $eps*$dExpect ) );
print "FAILED: Stats::Means::BinnedMean fails to compute correct errors on standard deviation\n"
    unless ( all(abs($de-$deExpect) < $eps*$deExpect) );

## Stats::Histograms
### Test binning.
my $xBins = pdl sequence(3);
my $yBins = pdl sequence(3);
my $x     = pdl 4.0*random(100);
my $y     = pdl 3.0*random(100);
my $w     = pdl random(100);
(my $h, my $e, my $c) = &Histograms::Histogram2D($xBins,$yBins,$x,$y,$w);
foreach my $i ( $xBins->list() ) {
    foreach my $j ( $yBins->list() ) {
	print "FAILED: Stats::Histograms::Histogram2D fails to compute correct weight in bin #1\n"
	    unless (
		&agree(
		     $h->(($i),($j)),
		     sum(
			 where(
			     $w,
			     ($x >= $i-0.5) & ($x < $i+0.5) & ($y >= $j-0.5) & ($y < $j+0.5)
			 )
		     ),
		     1.0e-3,
		     absolute => 1
		)
	    );
    }
}
### Test differential.
my $x2Bins = pdl 2*sequence(3);
my $y2Bins = pdl 2*sequence(3);
my $x2     = pdl 8.0*random(100);
my $y2     = pdl 6.0*random(100);
my $w2     = pdl random(100);
(my $h2, my $e2, my $c2) = &Histograms::Histogram2D($x2Bins,$y2Bins,$x2,$y2,$w2,differential => "xy");
for(my $i=0;$i<3;++$i) {
    for(my $j=0;$j<3;++$j) {
	print "FAILED: Stats::Histograms::Histogram2D fails to compute correct weight in bin #2\n"
	    unless (
		&agree(
		     $h2->(($i),($j)),
		     sum(
			 where(
			     $w2,
			     ($x2 >= 2.0*$i-1.0) & ($x2 < 2.0*$i+1.0) & ($y2 >= 2.0*$j-1.0) & ($y2 < 2.0*$j+1.0)
			 )
		     )
		     /4.0,
		     1.0e-3,
		     absolute => 1
		)
	    );
    }
}
### Test normalization.
#### xy
(my $hxyh, my $exyh, my $cxyh) = &Histograms::Histogram2D($xBins,$yBins,$x,$y,$w,normalized => "xy",normalizeBy => "histogram");
print "FAILED: Stats::Histograms::Histogram2D fails to normalize [xy:histogram]\n"
    unless ( &agree(sum($hxyh),1.0,1.0e-6) );
(my $hxyw, my $exyw, my $cxyw) = &Histograms::Histogram2D($xBins,$yBins,$x,$y,$w,normalized => "xy",normalizeBy => "weights"  );
print "FAILED: Stats::Histograms::Histogram2D fails to normalize [xy:weights]\n"
    unless (
	&agree(
	     sum($hxyw) 
	     , 
	     sum(
		 where(
		     $w,
		     ($x < 2.5)
		     &
		     ($y < 2.5)
		 )
	     )
	     /
	     sum($w)
	     ,
	     1.0e-6
	)
    );
#### x
(my $hxh, my $exh, my $cxh) = &Histograms::Histogram2D($xBins,$yBins,$x,$y,$w,normalized => "x", normalizeBy => "histogram");
(my $hxw, my $exw, my $cxw) = &Histograms::Histogram2D($xBins,$yBins,$x,$y,$w,normalized => "x", normalizeBy => "weights"  );
for(my $j=0;$j<3;++$j) {
    print "FAILED: Stats::Histograms::Histogram2D fails to normalize [x:histogram]\n"
	unless ( &agree(sum($hxh->(:,($j))),1.0,1.0e-6) );
    print "FAILED: Stats::Histograms::Histogram2D fails to normalize [x:weights]\n"
	unless ( &agree(sum($hxw->(:,($j))),sum(where($w,($x > -0.5) & ($x < 2.5) & ($y > $j-0.5) & ($y < $j+0.5)))/sum(where($w,($y > $j-0.5) & ($y < $j+0.5))),1.0e-6) );
}
#### y
(my $hyh, my $eyh, my $cyh) = &Histograms::Histogram2D($xBins,$yBins,$x,$y,$w,normalized => "y", normalizeBy => "histogram");
(my $hyw, my $eyw, my $cyw) = &Histograms::Histogram2D($xBins,$yBins,$x,$y,$w,normalized => "y", normalizeBy => "weights"  );
for(my $i=0;$i<3;++$i) {
    print "FAILED: Stats::Histograms::Histogram2D fails to normalize [y:histogram]\n"
	unless ( &agree(sum($hyh->(($i),:)),1.0,1.0e-6) );
    print "FAILED: Stats::Histograms::Histogram2D fails to normalize [y:weights]\n"
	unless ( &agree(sum($hyw->(($i),:)),sum(where($w,($x > $i-0.5) & ($x < $i+0.5) & ($y > -0.5) & ($y < 2.5)))/sum(where($w,($x > $i-0.5) & ($x < $i+0.5))),1.0e-6) );
}
### Test Gaussian smoothing
my $xs = pdl ( 1.0 );
my $ys = pdl ( 1.0 );
my $sx = pdl ( 1.0 );
my $sy = pdl ( 1.0 );
my $ws = pdl ( 1.0 );
(my $hs, my $es, my $cs) = &Histograms::Histogram2D($xBins,$yBins,$xs,$ys,$ws,gaussianSmoothX => $sx,gaussianSmoothY => $sy);
my $hst = pdl
    [
     [0.058433556, 0.092564571, 0.058433556],
     [0.092564571,   0.1466315, 0.092564571],
     [0.058433556, 0.092564571, 0.058433556]
    ];
print "FAILED: Stats::Histograms::Histogram2D fails to smooth\n"
    unless ( &agree($hs,$hst,1.0e-6) );
print "FAILED: Stats::Histograms::Histogram2D covariance fails for gaussian smoothing\n"
    unless ( &agree(sqrt($cs->diagonal(0,1)),$es->flat(),1.0e-6) );
#### Normalized in x
(my $hsxh, my $esxh, my $csxh) = &Histograms::Histogram2D($xBins,$yBins,$xs,$ys,$ws,gaussianSmoothX => $sx,gaussianSmoothY => $sy, normalized => "x", normalizeBy => "histogram");
(my $hsxw, my $esxw, my $csxw) = &Histograms::Histogram2D($xBins,$yBins,$xs,$ys,$ws,gaussianSmoothX => $sx,gaussianSmoothY => $sy, normalized => "x", normalizeBy => "weights");
for(my $j=0;$j<3;++$j) {
    print "FAILED: Stats::Histograms::Histogram2D fails to normalize [x:histogram]\n"
	unless ( &agree(sum($hsxh->(:,($j))),1.0000000000,1.0e-6) );
    print "FAILED: Stats::Histograms::Histogram2D fails to normalize [x:weights]\n"
	unless ( &agree(sum($hsxw->(:,($j))),0.8663855075,1.0e-6) );
}
####  Note that covariance should be zero as we have a single point and we've normalized it.
my $csxht = pdl zeroes(9);
print "FAILED: Stats::Histograms::Histogram2D covariance fails for gaussian smoothing\n"
    unless ( &agree(sqrt($csxh->diagonal(0,1)),$csxht,1.0e-6,absolute => 1) );
#### Normalized in y
(my $hsyh, my $esyh, my $csyh) = &Histograms::Histogram2D($xBins,$yBins,$xs,$ys,$ws,gaussianSmoothX => $sx,gaussianSmoothY => $sy, normalized => "y", normalizeBy => "histogram");
(my $hsyw, my $esyw, my $csyw) = &Histograms::Histogram2D($xBins,$yBins,$xs,$ys,$ws,gaussianSmoothX => $sx,gaussianSmoothY => $sy, normalized => "y", normalizeBy => "weights");
for(my $i=0;$i<3;++$i) {
    print "FAILED: Stats::Histograms::Histogram2D fails to normalize [y:histogram]\n"
	unless ( &agree(sum($hsyh->(($i),:)),1.0000000000,1.0e-6) );
    print "FAILED: Stats::Histograms::Histogram2D fails to normalize [y:weights]\n"
	unless ( &agree(sum($hsyw->(($i),:)),0.8663855075,1.0e-6) );
}
####  Note that covariance should be zero as we have a single point and we've normalized it.
my $csyht = pdl zeroes(9);
print "FAILED: Stats::Histograms::Histogram2D covariance fails for gaussian smoothing\n"
    unless ( &agree(sqrt($csyh->diagonal(0,1)),$csyht,1.0e-6,absolute => 1) );

exit;

sub agree {
    # Return true if two values agree to within the given fractional tolerance.
    my $a         = shift;
    my $b         = shift;
    my $tolerance = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};
    $options{'absolute'} = 0
	unless ( exists($options{'absolute'}) );
    if ( UNIVERSAL::isa($a, 'PDL') ) {
	my $scale = 0.5*(abs($a)+abs($b));
	$scale .= 1.0
	    if ( $options{'absolute'} == 1 );
	if ( all(abs($a-$b) < $tolerance*$scale) ) {
	    return 1;
	} else {
	    return 0;
	}
    } else {
	my $scale = 0.5*(abs($a)+abs($b));
	$scale = 1.0
	    if ( $options{'absolute'} == 1 );
	if ( abs($a-$b) < $tolerance*$scale ) {
	    return 1;
	} else {
	    return 0;
	}
    }
}
