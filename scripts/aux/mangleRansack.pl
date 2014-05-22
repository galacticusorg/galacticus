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
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
require Mangle;

# Generate a set of random points in a survey window using mangle (http://space.mit.edu/~molly/mangle/).
# Andrew Benson (18-April-2014)

# Get arguments.
die("Usage: mangleRansack.pl <file1> [<file2>]....... <outputFile> <nPoints>")
    unless ( scalar(@ARGV) > 1 );
my @files          = @ARGV;
my $pointCount     = pop(@files);
my $outputFileName = pop(@files);

# Build mangle.
&Mangle::Build;

# Determine the solid angle of each file.
my $solidAngles = pdl [];
foreach my $file ( @files ) {
    my $solidAngle;
    open(my $pipe,"aux/mangle/bin/harmonize ".$file." /dev/null |");
    while ( my $line = <$pipe> ) {
	if ( $line =~ m/^area of \(weighted\) region is ([0-9\.]+) str/ ) {
	    $solidAngle = $1;
	}
    }
    close($pipe);
    die("mangleRansack.pl: unable to determine solid angle for file ".$file)
	unless ( defined($solidAngle) );
    $solidAngles = $solidAngles->append($solidAngle);
}

# Determine fractions of total solid angle in each window.
my $fraction = $solidAngles/sum($solidAngles);

# Determine number of points in each window.
my $pointCounts = long($pointCount*$fraction);

# Generate random points for each window.
my $theta;
my $phi;
if ( $pointCount > 0 ) {
    my $azimuth   = pdl [];
    my $elevation = pdl [];
    for(my $i=0;$i<scalar(@files);++$i) {
	system("aux/mangle2.2/bin/ransack -r".$pointCounts->(($i))->sclr()." ".$files[$i]." ".$files[$i].".ran");
	die("mangleRansack.pl: unable to generate random points for ".$files[$i])
	    unless ( $? == 0 );
	# Read points.
	(my $thisAzimuth, my $thisElevation) = rcols($files[$i].".ran",0,1,{LINES => '1:-1'});
	$azimuth   = $azimuth  ->append($thisAzimuth  );
	$elevation = $elevation->append($thisElevation);
    }
    # Convert to spherical coordinates.
    $theta = PI*(0.5-$elevation/180.0);
    $phi   = PI*     $azimuth  /180.0 ;
}

# Output.
my $outputFile = new PDL::IO::HDF5(">".$outputFileName);
$outputFile->attrSet(files => join(" : ",@files));
$outputFile->dataset('solidAngle')->set($solidAngles);
if ( $pointCount > 0 ) {
    $outputFile->dataset('theta')->set($theta);
    $outputFile->dataset('phi'  )->set($phi  );
}

exit;
