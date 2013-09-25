#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V091"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use DateTime;
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::IO::HDF5;
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Construct a set of random points that lie within the angular mask of the UKIDSS UDS sample used by Caputi et al. (2011;
# http://adsabs.harvard.edu/abs/2011MNRAS.413..162C). The mask is defined by a set of boundaries, plus a set of rectangular
# cut-outs. The boundaries are hard-coded into this script. The rectangles required are found by randomly placing rectangles of an
# initial minimum size into the field and growing them to be as large as possible while not containing any galaxies. Once enough
# suqares have been placed, random points are generated and rejected if they lie outside of the bounaries or inside any cut-out
# rectangle.
# Andrew Benson (23-August-2012)

# Define a working directory.
my $workDirectory = $galacticusPath."constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/";

# Read galaxy positions.
(my $RA, my $dec) = rcols($workDirectory."data/surveyGeometry.txt",{EXCLUDE => '/#/'});

# Specify minimum rectangle edge size (in degrees).
my $rectangleSizeMinimum = pdl 0.01;

# Initialize an empty list of rectangles.
my @rectangles;

# Define boundaries of the survey region. These are specifed by two points which define a line and a direction indicating whether
# the boundary excludes the region above, below, left, or right of the line.
my @boundaries;
push(
    @boundaries,
    {
	x1 => 34.460,
	y1 => -4.648,
	x2 => 34.750,
	y2 => -4.648,
	type => "above"
    },
    {
	x1 => 34.000,
	y1 => -4.908,
	x2 => 34.460,
	y2 => -4.648,
	type => "above"
    },
    {
    	x1 => 34.900,
    	y1 => -4.920,
    	x2 => 34.750,
    	y2 => -4.648,
    	type => "above"
    },
    {
    	x1 => 34.440,
    	y1 => -5.518,
    	x2 => 34.905,
    	y2 => -5.255,
    	type => "below"
    },
    {
    	x1 => 34.210,
    	y1 => -5.518,
    	x2 => 34.440,
    	y2 => -5.518,
    	type => "below"
    },
    {
    	x1 => 34.210,
    	y1 => -5.518,
    	x2 => 34.000,
    	y2 => -5.150,
    	type => "below"
    },
    {
    	x1 => 34.000,
    	y1 => -4.908,
    	x2 => 34.000,
    	y2 => -5.150,
    	type => "left"
    },
    {
    	x1 => 34.900,
    	y1 => -4.920,
    	x2 => 34.905,
    	y2 => -5.255,
    	type => "right"
    }
    );

# Initialize count of failed rectangle placement attempts and specify the maximum number of failed attempts allowed before we
# decide that we've filled all rectangles possible.
my $failCount   =   0;
my $failMaximum = 100;

# Do placement of rectangles.
while ( $failCount < $failMaximum ) {

    # Choose initial rectangle size.
    my $sizeX = pdl $rectangleSizeMinimum;
    my $sizeY = pdl $rectangleSizeMinimum;

    # Choose initial location.
    my $i                 =    0; # Number of iterations of rectangle placement completed.
    my $iterationsMaximum = 1000; # Maximum number of iterations to try before failing.
    my $good              =    0; # Is this rectangle "good" (i.e. within boundaries and not overlapping any galaxies)?
    my $x;
    my $y;
    while ( $i < $iterationsMaximum && $good == 0 ) {
	++$i;                        # Increment iteration count.
	$x = pdl random(1)*0.9+34.0; # Select x and y positions of the rectangle center.
	$y = pdl random(1)*1.0- 5.6;
	# Test if rectangle is "good".
	$good = &rectangleIsAcceptable($x,$y,$sizeX,$sizeY);
    }

    # Is the rectangle good...?
    if ( $good == 1 ) {
	# It is, attempt to grow it to fill the available space.
	my $maxSizeReached = 0;
	my $growFactor     = 0.01;
	while ( $maxSizeReached < 4 ) { # Only stop when all four sides of the rectangle can grow no further.
	    $maxSizeReached = 0;
	    # Grow left-edge.
	    $x          -= 0.5*$growFactor*$sizeX;
	    $sizeX      *= 1.0+$growFactor;
	    $good        = &rectangleIsAcceptable($x,$y,$sizeX,$sizeY);
	    if ( $good == 0 ) {
		$maxSizeReached += 1;
		$sizeX         /= 1.0+$growFactor;
		$x              += 0.5*$growFactor*$sizeX;
	    }
	    # Grow right-edge.
	    $x          += 0.5*$growFactor*$sizeX;
	    $sizeX      *= 1.0+$growFactor;
	    $good        = &rectangleIsAcceptable($x,$y,$sizeX,$sizeY);
	    if ( $good == 0 ) {
		$maxSizeReached += 1;
		$sizeX          /= 1.0+$growFactor;
		$x              -= 0.5*$growFactor*$sizeX;
	    }
	    # Grow bottom-edge.
	    $y          -= 0.5*$growFactor*$sizeY;
	    $sizeY      *= 1.0+$growFactor;
	    $good        = &rectangleIsAcceptable($x,$y,$sizeX,$sizeY);
	    if ( $good == 0 ) {
		$maxSizeReached += 1;
		$sizeY          /= 1.0+$growFactor;
		$y              += 0.5*$growFactor*$sizeY;
	    }
	    # Grow top-edge.
	    $y          += 0.5*$growFactor*$sizeY;
	    $sizeY      *= 1.0+$growFactor;
	    $good        = &rectangleIsAcceptable($x,$y,$sizeX,$sizeY);
	    if ( $good == 0 ) {
		$maxSizeReached += 1;
		$sizeY          /= 1.0+$growFactor;
		$y              -= 0.5*$growFactor*$sizeY;
	    }
	}
	
	# Add the rectangle to the rectangle list.
	push(
	    @rectangles,
	    {
		x     => $x,
		y     => $y,
		sizeX => $sizeX,
		sizeY => $sizeY,
	    }
	    );
	my $rectangleCount = scalar(@rectangles);
	print "Placed rectangle ".$rectangleCount." at (".$x.",".$y.") with size ".$sizeX."x".$sizeY."\n";
    } else {
	++$failCount;
    }   
}

# Generate random points.
# Specify range of RA and dec that encompasses the entire survey.
my $cosThetaMin = pdl cos(1.6517967);
my $cosThetaMax = pdl cos(1.6675065);
my $phiMin      = pdl 0.5931317;
my $phiMax      = pdl 0.60971425;
# Generate points at random within this region.
my $randomCount = 1000000;
my $phi         = pdl random($randomCount)*($phiMax     -$phiMin     )+$phiMin;
my $cosTheta    = pdl random($randomCount)*($cosThetaMax-$cosThetaMin)+$cosThetaMin;
my $theta       = acos($cosTheta);
my $Pi          = pdl 3.1415927;
my $randomRA    = $phi*180.0/$Pi;
my $randomDec   = 90.0-$theta*180.0/$Pi;
# Find the subset of points that are acceptable (i.e. within bounds and not inside any rectangle).
my $randomAccept     = &pointsAreAcceptable($randomRA,$randomDec);
# Compute the solid angle of the survey using the random points to do a Monte-Carlo integration.
my $solidAngleTotal  = ($phiMax-$phiMin)*($cosThetaMin-$cosThetaMax);
my $solidAngleSurvey = $solidAngleTotal*nelem($randomAccept)/$randomCount;
print "Survey solid angle is: ".$solidAngleSurvey."\n";

# Rotate random points such that pole of coordinate system lies approximately in center of survey field.
# Approximate center of the field
my $theta0     = pdl 95.05*$Pi/180.0;
my $phi0       = pdl 34.45*$Pi/180.0;
# Define a coordinate system with c1 vector pointing to the field center.
my $c1         = pdl ( sin($theta0)*cos($phi0), sin($theta0)*sin($phi0), cos($theta0) );
my $c2         = pdl ( sin($phi0), -cos($phi0), 0.0 );
my $c3         = crossp($c1,$c2);
# Construct Cartesian coordinate representation of random points.
my $x          = pdl [ sin($theta)*cos($phi), sin($theta)*sin($phi), cos($theta) ];
# Project onto new coordinate vectors.
my $x1         = transpose($x) x transpose($c1);
my $x2         = transpose($x) x transpose($c2);
my $x3         = transpose($x) x transpose($c3);
# Convert to theta and phi in this new coordinate system.
my $thetaPrime = transpose(acos ($x1    ))->(:,(0));
my $phiPrime   = transpose(atan2($x3,$x2))->(:,(0));

# Write random points to file.
my $hdfFile = new PDL::IO::HDF5(">".$workDirectory."data/surveyGeometryRandoms.hdf5");
$hdfFile->dataset('theta')->set($thetaPrime->index($randomAccept));
$hdfFile->dataset('phi'  )->set($phiPrime  ->index($randomAccept));
my $dt = DateTime->now->set_time_zone('local');
(my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
$hdfFile->attrSet(
    createdBy   => "Galacticus; constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/surveyGeometryRandoms.pl",
    description => "Random points in survey geometry of UKIDSS UDS for Caputi et al. (2011; http://adsabs.harvard.edu/abs/2011MNRAS.413..162C)",
    timestamp   => $now
    );

# Plot the results.
my $plot;
my $gnuPlot;
my $plotFile = $workDirectory."surveyMask.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Survey Mask for Caputi et al. (2011) Sample'\n";
print $gnuPlot "set xlabel 'right ascension [\$^\\circ\$]'\n";
print $gnuPlot "set ylabel 'declination [\${^\\circ}\$]'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.19,0.8\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set xrange [33.98:34.93]\n";
print $gnuPlot "set yrange [-5.53:-4.63]\n";
print $gnuPlot "set pointsize 2.0\n";
# Add survey boundaries.
foreach my $boundary ( @boundaries ) {
    my $sx = pdl ( $boundary->{'x1'}, $boundary->{'x2'} );
    my $sy = pdl ( $boundary->{'y1'}, $boundary->{'y2'} );
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$sx,
	$sy,
	style      => "line",
	weight     => [1,1],
	color      => $PrettyPlots::colorPairs{'indianRed'},
	);
}
# Add rectangles.
foreach my $rectangle ( @rectangles ) {
    my $sx = pdl (
	$rectangle->{'x'}-0.5*$rectangle->{'sizeX'}, 
	$rectangle->{'x'}+0.5*$rectangle->{'sizeX'}
	);
    my $syl = pdl (
	$rectangle->{'y'}-0.5*$rectangle->{'sizeY'}, 
	$rectangle->{'y'}-0.5*$rectangle->{'sizeY'} 
	);
    my $syu = pdl (
	$rectangle->{'y'}+0.5*$rectangle->{'sizeY'}, 
	$rectangle->{'y'}+0.5*$rectangle->{'sizeY'} 
	);
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$sx->((0),:),
	$syl->((0),:),,
	y2         => $syu->((0),:),,
	style      => "filledCurve",
	weight     => [1,1],
	color      => $PrettyPlots::colorPairs{'cornflowerBlue'},
	);
}
# Add galaxies.
my $uniformRandom = random(nelem($RA));
my $selection     = which($uniformRandom <= 1.0);
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $RA ->index($selection),
    $dec->index($selection),
    style      => "point",
    weight     => [1,1],
    symbol     => [0,0],
    color      => $PrettyPlots::colorPairs{'mediumSeaGreen'},
    );
# Add random points.
my $uniformRandomRandom = random(nelem($randomAccept));
my $selectionRandom     = which($uniformRandomRandom <= 0.1);
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $randomRA ->index($randomAccept->index($selectionRandom)),
    $randomDec->index($randomAccept->index($selectionRandom)),
    style      => "point",
    weight     => [1,1],
    symbol     => [0,0],
    color      => $PrettyPlots::colorPairs{'redYellow'},
    );
# Finalize the plot.
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;

sub rectangleIsAcceptable {
    my $x     = shift;
    my $y     = shift;
    my $sizeX = shift;
    my $sizeY = shift;

    return 0
	if ( any( (abs($RA-$x) < 0.5*$sizeX) & (abs($dec-$y) < 0.5*$sizeY) ) );

    my $acceptable = 1;
    foreach my $boundary ( @boundaries ) {
	if ( $boundary->{'type'} eq "above" || $boundary->{'type'} eq "below" ) {
	    my $yLimit;
	    $yLimit = $y+0.5*$sizeY
		if ( $boundary->{'type'} eq "above" );
	    $yLimit = $y-0.5*$sizeY
		if ( $boundary->{'type'} eq "below" );
	    foreach my $xLimit ( $x ) { #$x-0.5*$sizeX, $x+0.5*$sizeX ) {
		my $boundaryY = ($boundary->{'y2'}-$boundary->{'y1'})*($xLimit-$boundary->{'x1'})/($boundary->{'x2'}-$boundary->{'x1'})+$boundary->{'y1'};
		$acceptable = 0
		    if ( $boundaryY < $yLimit && $boundary->{'type'} eq "above" );
		$acceptable = 0
		    if ( $boundaryY > $yLimit && $boundary->{'type'} eq "below" );
	    }
	}
	if ( $boundary->{'type'} eq "left" || $boundary->{'type'} eq "right" ) {
	    my $xLimit;
	    $xLimit = $x+0.5*$sizeX
		if ( $boundary->{'type'} eq "right" );
	    $xLimit = $x-0.5*$sizeX
		if ( $boundary->{'type'} eq "left" );
	    foreach my $yLimit ( $y ) { #( $y-0.5*$sizeY, $y+0.5*$sizeY ) {
		my $boundaryX = ($boundary->{'x2'}-$boundary->{'x1'})*($yLimit-$boundary->{'y1'})/($boundary->{'y2'}-$boundary->{'y1'})+$boundary->{'x1'};
		$acceptable = 0
		    if ( $boundaryX < $xLimit && $boundary->{'type'} eq "right" );
		$acceptable = 0
		    if ( $boundaryX > $xLimit && $boundary->{'type'} eq "left" );
	    }
	}
    }	
    return 0
	if ( $acceptable == 0 );

    foreach my $rectangle ( @rectangles ) {
	return 0
	    if (
		$rectangle->{'x'}+0.5*$rectangle->{'sizeX'} > $x &&
		$rectangle->{'x'}-0.5*$rectangle->{'sizeX'} < $x &&
		$rectangle->{'y'}+0.5*$rectangle->{'sizeY'} > $y &&
		$rectangle->{'y'}-0.5*$rectangle->{'sizeY'} < $y
	    );
    }

    return 1;
}

sub pointsAreAcceptable {
    my $pRA  = shift;
    my $pDec = shift;

    my $acceptable = pdl sequence(nelem($pRA));

    foreach my $boundary ( @boundaries ) {
	my $inBoundary;
	$inBoundary = which($pDec < ($boundary->{'y2'}-$boundary->{'y1'})*($pRA-$boundary->{'x1'})/($boundary->{'x2'}-$boundary->{'x1'})+$boundary->{'y1'})
	    if ( $boundary->{'type'} eq "above" );
	$inBoundary = which($pDec > ($boundary->{'y2'}-$boundary->{'y1'})*($pRA-$boundary->{'x1'})/($boundary->{'x2'}-$boundary->{'x1'})+$boundary->{'y1'})
	    if ( $boundary->{'type'} eq "below" );
	$inBoundary = which($pRA < ($boundary->{'x2'}-$boundary->{'x1'})*($pDec-$boundary->{'y1'})/($boundary->{'y2'}-$boundary->{'y1'})+$boundary->{'x1'})
	    if ( $boundary->{'type'} eq "right" );
	$inBoundary = which($pRA > ($boundary->{'x2'}-$boundary->{'x1'})*($pDec-$boundary->{'y1'})/($boundary->{'y2'}-$boundary->{'y1'})+$boundary->{'x1'})
	    if ( $boundary->{'type'} eq "left" );
	if ( defined($inBoundary) ) {
	    my $stillAcceptable = setops($acceptable,'AND',$inBoundary);
	    $acceptable = pdl $stillAcceptable;
	}
    }	

    foreach my $rectangle ( @rectangles ) {
	my $outsideRectangle = 
	    which(
		($pRA  < $rectangle->{'x'}-0.5*$rectangle->{'sizeX'}) |
		($pRA  > $rectangle->{'x'}+0.5*$rectangle->{'sizeX'}) |
		($pDec < $rectangle->{'y'}-0.5*$rectangle->{'sizeY'}) |
		($pDec > $rectangle->{'y'}+0.5*$rectangle->{'sizeY'})
	    );
	my $stillAcceptable = setops($acceptable,'AND',$outsideRectangle);
	$acceptable = pdl $stillAcceptable;
    }
    
    return $acceptable;
}
