#!/usr/bin/env perl

# Make a plot of the specified transfer function file.
# Andrew Benson (27-Jan-2009)

use XML::Simple;
use Data::Dumper;
use Graphics::GnuplotIF;
use POSIX;

unless ( $#ARGV == 0 || $#ARGV == 1 ) {die "Usage: cooling_function.pl <coolingFunctionFile> [<pdfFile>]"};
$coolingFunctionFile = $ARGV[0];
if ( $#ARGV == 1 ) {
    $pdfFile         = $ARGV[1];
} else {
    ($pdfFile = $coolingFunctionFile) =~ s/^data\/(.+)\.xml$/plots\/\1\.pdf/;
}

# Initialize the plot.
$plot1  = Graphics::GnuplotIF->new();
$plot1->gnuplot_hardcopy( '| ps2pdf - '.$pdfFile, 
			  'postscript enhanced', 
			  'color lw 3' );
$plot1->gnuplot_set_xlabel("Temperature [K]");
$plot1->gnuplot_set_ylabel("{/Symbol L}(T) [erg cm^3 s^{-1}]");
$plot1->gnuplot_set_title("Atomic CIE Cooling Function (Cloudy 08.00), colored by log_{10}(Z/Z_{{/=12 O}&{/*-.66 O}{/=12 \267}})");
$plot1->gnuplot_cmd("set logscale xy");
$plot1->gnuplot_cmd("set format y '10^{\%L}'");
$plot1->gnuplot_cmd("set format x '10^{\%L}'");
$plot1->gnuplot_cmd("set palette rgbformulae 33,13,10");
$plot1->gnuplot_cmd("set key off");

# Read the XML data file.
$xml = new XML::Simple;
$data = $xml->XMLin($coolingFunctionFile);
@coolingFunctionsArray = @{$data -> {'coolingFunction'}};

# Create plot titles.
$lZMin = 100.0;
$lZMax = -100.0;
for ($iCoolingFunction=0;$iCoolingFunction<=$#coolingFunctionsArray;++$iCoolingFunction) {
    $metallicity = ${$coolingFunctionsArray[$iCoolingFunction]}{'metallicity'};
    if ( $metallicity <= -999.0 ) {
	$plotTitles[++$#plotTitles] = "Primordial";
    } else {
	$plotTitles[++$#plotTitles] = "log_{10}(Z/Z_{{/=12 O}&{/*-.66 O}{/=12 \267}}) = ".$metallicity;
	if ( $metallicity > $lZMax ) {$lZMax=$metallicity};
	if ( $metallicity < $lZMin ) {$lZMin=$metallicity};
    }
}

$plot1->gnuplot_set_plot_titles( @plotTitles );
$plot1->gnuplot_cmd("set cbrange [".$lZMin.":".$lZMax."]");

# Loop over cooling functions.
for ($iCoolingFunction=0;$iCoolingFunction<=$#coolingFunctionsArray;++$iCoolingFunction) {
    
    # Get the data for this cooling function.
    @{$x[$iCoolingFunction]} = @{$coolingFunctionsArray[$iCoolingFunction]->{'temperature'}->{'datum'}};
    @{$y[$iCoolingFunction]} = @{$coolingFunctionsArray[$iCoolingFunction]->{'coolingRate'}->{'datum'}};
    ${$hashArray[$iCoolingFunction]}{'x_values'} = \@{$x[$iCoolingFunction]};
    ${$hashArray[$iCoolingFunction]}{'y_values'} = \@{$y[$iCoolingFunction]};
    $cFrac = $iCoolingFunction/$#coolingFunctionsArray;
    ${$hashArray[$iCoolingFunction]}{'style_spec'} = "lines linetype 1 palette frac ".$cFrac;

    $refArray[$iCoolingFunction] = \%{$hashArray[$iCoolingFunction]};
}

# Make the plot.
$plot1->gnuplot_plot_many_style( @refArray );
    
exit;
