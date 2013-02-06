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
use PDL::IO::HDF5;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require Galacticus::HDF5;
require XMP::MetaData;

# Plot the spectrum of a galaxy computed by Galacticus+Grasil.
# Andrew Benson (18-January-2013)

# Get ourself.
my $self             = $0;

# Read arguments.
die("Usage: plotGrasilSpectrum.pl <galacticusFile> <output> <mergerTreeIndex> <nodeIndex> <inclination> <outputFile>")
    unless ( scalar(@ARGV) == 6 );
my $galacticusFile   = $ARGV[0];
my $outputNumber     = $ARGV[1];
my $mergerTreeIndex  = $ARGV[2];
my $nodeIndex        = $ARGV[3];
my $inclination      = $ARGV[4];
my $outputFile       = $ARGV[5];

# Access the Galacticus file.
my $galacticus       = new PDL::IO::HDF5($galacticusFile);
   
# Access the SED group.
my $groupName        = "grasilSEDs/Output".$outputNumber."/mergerTree".$mergerTreeIndex."/node".$nodeIndex;
my $sedGroup         = $galacticus->group($groupName);

# Read wavelengths.
my $wavelengths      = $sedGroup->dataset('wavelength')->get();
$wavelengths        /= 1.0e4; # Convert to microns.

# Read inclinations.
my $inclinations     = $sedGroup->dataset('inclination')->get();

# Read SED.
my $seds             = $sedGroup->dataset('SED')->get();

# Interpolate the SED to the requested inclination.
my $inclinationIndex = interpol($inclination,$inclinations,sequence(nelem($inclinations)));
my $index            = pdl zeroes(2,nelem($wavelengths));
$index->((0),:)     .= $inclinationIndex;
$index->((1),:)     .= sequence(nelem($wavelengths));
my $sed              = $seds->interpND($index);

# Make plot of stellar mass function.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Galacticus+Grasil spectrum for output ".$outputNumber.", tree ".$mergerTreeIndex.", node ".$nodeIndex."' offset 0,-0.25\n";
print $gnuPlot "set xlabel 'Wavelength; [\$\\mu\$m]'\n";
print $gnuPlot "set ylabel 'Luminosity; [\$10^{30}\$ erg s\$^{-1}\$ \\AA\$^{-1}\$]'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e-2:1.0e5]\n";
print $gnuPlot "set yrange [1.0e1:1.0e11]\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $wavelengths,
    $sed,
    style      => "line",
    weight     => [3,1],
    color      => $PrettyPlots::colorPairs{'redYellow'},
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS, margin => 2);
&MetaData::Write($plotFile,$galacticusFile,$self);

exit;
