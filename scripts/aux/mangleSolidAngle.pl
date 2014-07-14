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
    my $solidAngle = pdl 0.0;
    my @subFiles;
    if ( $file =~ m/^[\+\-]/ ) {
	@subFiles = split(/\:/,$file);
	print scalar(@subFiles)."\n";
    } else {
	push(@subFiles,"+".$file);
    }
    foreach my $subFile ( @subFiles ) {
	print $subFile."\n";
	if ( $subFile =~ m/^([\+\-])(.*)/ ) {
	    my $sign     = $1;
	    my $fileName = $2;
	    my $multiplier;
	    if ( $sign eq "+" ) {
		$multiplier = pdl +1.0;
	    } else {
		$multiplier = pdl -1.0;
	    }
	    my $subSolidAngle;
	    open(my $pipe,"aux/mangle/bin/harmonize ".$fileName." /dev/null |");
	    while ( my $line = <$pipe> ) {
		if ( $line =~ m/^area of \(weighted\) region is ([0-9\.]+) str/ ) {
		    $subSolidAngle = $1;
		}
	    }
	    close($pipe);
	    die("mangleRansack.pl: unable to determine solid angle for file ".$file)
		unless ( defined($solidAngle) );
	    $solidAngle += $multiplier*$subSolidAngle;
	}
    }
    $solidAngles = $solidAngles->append($solidAngle);
}

# Output.
my $outputFile = new PDL::IO::HDF5(">".$outputFileName);
$outputFile->attrSet(files => join(" : ",@files));
$outputFile->dataset('solidAngle')->set($solidAngles);

exit;
