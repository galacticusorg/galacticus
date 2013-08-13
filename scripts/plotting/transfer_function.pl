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
use XML::Simple;
use Data::Dumper;
use PDL;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Make a plot of the specified transfer function file.
# Andrew Benson (15-Dec-2009)

die "Usage: transfer_function.pl <transferFunctionFile> [<pdfFile>]"
    unless ( scalar(@ARGV) == 1 || scalar(@ARGV) == 2 );
my $transferFunctionFile = $ARGV[0];
my $pdfFile;
if ( scalar(@ARGV) == 2 ) {
    $pdfFile          = $ARGV[1];
} else {
    ($pdfFile = $transferFunctionFile) =~ s/^data\/(.+)\.xml$/plots\/$1\.pdf/;
}

# Read the XML data file.
my $xml = new XML::Simple;
my $data = $xml->XMLin($transferFunctionFile);
my @dataArray = @{$data -> {'datum'}};

# Extract the data.
my $x = pdl [];
my $y = pdl [];
foreach my $datum ( @dataArray ) {
    my @columns = split(/\s+/,$datum);
    $x = $x->append($columns[0]);
    $y = $y->append($columns[1]);
}

# Make a plot of the transfer function.
my $plot;
my $gnuPlot;
my $plotFile = $pdfFile;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Transfer Function' offset screen 0,-0.02\n";
print $gnuPlot "set xlabel 'Wavenumber [Mpc\$^{-1}\$]'\n";
print $gnuPlot "set ylabel 'Transfer function'\n";
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
my $xMinimum = minimum($x)*0.9;
my $xMaximum = maximum($x)/0.9;
my $yMinimum = minimum($y)*0.9;
my $yMaximum = maximum($y)/0.9;
print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$y,
    style  => "line",
    weight => [3,1],
    color  => $PrettyPlots::colorPairs{'redYellow'}
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS, margin => 2);

exit;
