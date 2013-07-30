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
use XML::Simple;
use Data::Dumper;
use PDL;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Make a plot of the specified transfer function file.
# Andrew Benson (27-Jan-2009)

die "Usage: cooling_function.pl <coolingFunctionFile> [<pdfFile>]"
    unless ( scalar(@ARGV) == 1 || scalar(@ARGV) == 2 );
my $coolingFunctionFile = $ARGV[0];
my $pdfFile;
if ( scalar(@ARGV) == 2 ) {
    $pdfFile         = $ARGV[1];
} else {
    ($pdfFile = $coolingFunctionFile) =~ s/^data\/(.+)\.xml$/plots\/$1\.pdf/;
}

# Read the XML data file.
my $xml = new XML::Simple;
my $data = $xml->XMLin($coolingFunctionFile);
my @coolingFunctionsArray = @{$data -> {'coolingFunction'}};

# Create plot titles.
my $lZMin = +100.0;
my $lZMax = -100.0;
my @plotTitles;
for (my $iCoolingFunction=0;$iCoolingFunction<scalar(@coolingFunctionsArray);++$iCoolingFunction) {
    my $metallicity = ${$coolingFunctionsArray[$iCoolingFunction]}{'metallicity'};
    if ( $metallicity <= -999.0 ) {
	push(@plotTitles,"Primordial");
    } else {
	push(@plotTitles,"\$\\log_{10}(Z/Z_\\odot) = ".$metallicity."\$");
	if ( $metallicity > $lZMax ) {$lZMax=$metallicity};
	if ( $metallicity < $lZMin ) {$lZMin=$metallicity};
    }
}

# Make a plot of the cooling function.
my $plot;
my $gnuPlot;
my $plotFile = $pdfFile;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Atomic CIE Cooling Function (Cloudy 13.02), colored by \$\\log_{10}(Z/Z_\\odot)\$'\n";
print $gnuPlot "set xlabel 'Temperature [K]'\n";
print $gnuPlot "set ylabel '\$\\Lambda(T)\$ [erg cm\$^3\$ s\$^{-1}\$]'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key off\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set palette rgbformulae 33,13,10\n";
print $gnuPlot "set cbrange [".$lZMin.":".$lZMax."]\n";
print $gnuPlot "set xrange [316.0:1.0e9]\n";
print $gnuPlot "set yrange [3.0e-30:1.0e-19]\n";
print $gnuPlot "set pointsize 2.0\n";
# Loop over cooling functions.
for (my $iCoolingFunction=0;$iCoolingFunction<scalar(@coolingFunctionsArray);++$iCoolingFunction) {
    # Get the data for this cooling function.
    my $x = pdl @{$coolingFunctionsArray[$iCoolingFunction]->{'temperature'}->{'datum'}};
    my $y = pdl @{$coolingFunctionsArray[$iCoolingFunction]->{'coolingRate'}->{'datum'}};
    my $cFrac = $iCoolingFunction/(scalar(@coolingFunctionsArray)-1);
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$x,$y,
	style  => "line",
	weight => [3,1],
	color  => ["palette frac ".$cFrac,"palette frac ".$cFrac]
	);
}
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
