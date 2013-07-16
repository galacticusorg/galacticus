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
use Graphics::GnuplotIF;
use POSIX;

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
my @x;
my @y;
foreach my $datum ( @dataArray ) {
    my @columns = split(/\s+/,$datum);
    push(@x,$columns[0]);
    push(@y,$columns[1]);
}

# Make the plot.
my $plot1  = Graphics::GnuplotIF->new();
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
