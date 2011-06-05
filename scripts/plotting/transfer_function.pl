#!/usr/bin/env perl

# Make a plot of the specified transfer function file.
# Andrew Benson (15-Dec-2009)

use XML::Simple;
use Data::Dumper;
use Graphics::GnuplotIF;
use POSIX;

unless ( $#ARGV == 0 || $#ARGV == 1 ) {die "Usage: transfer_function.pl <transferFunctionFile> [<pdfFile>]"};
$transferFunctionFile = $ARGV[0];
if ( $#ARGV == 1 ) {
    $pdfFile          = $ARGV[1];
} else {
    ($pdfFile = $transferFunctionFile) =~ s/^data\/(.+)\.xml$/plots\/\1\.pdf/;
}

# Read the XML data file.
$xml = new XML::Simple;
$data = $xml->XMLin($transferFunctionFile);
@dataArray = @{$data -> {'datum'}};

# Extract the data.
foreach $datum ( @dataArray ) {
    @columns = split(/\s+/,$datum);
    $x[++$#x] = $columns[0];
    $y[++$#y] = $columns[1];
}

# Make the plot.
$plot1  = Graphics::GnuplotIF->new();
$plot1->gnuplot_hardcopy( '| ps2pdf - '.$pdfFile, 
			  'postscript enhanced', 
			  'color lw 3' );
$plot1->gnuplot_set_xlabel("Wavenumber [Mpc^{-1}]");
$plot1->gnuplot_set_ylabel("Transfer function");
$plot1->gnuplot_set_plot_titles( "T(k)" );
$plot1->gnuplot_set_title("Transfer Function");
$plot1->gnuplot_cmd("set logscale xy");
$plot1->gnuplot_cmd("set format y \"10^{\%L}\"");
$plot1->gnuplot_cmd("set format x \"10^{\%L}\"");
$plot1->gnuplot_plot_xy( \@x, \@y );

exit;
