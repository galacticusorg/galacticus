#!/usr/bin/env perl
use lib "./perl";
use PDL;
use XML::Simple;
use Graphics::GnuplotIF;
use Galacticus::HDF5;
use Math::SigFigs;
use Data::Dumper;
use Stats::Histograms;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Stellar_Mass_Function.pl <galacticusFile> <outputDir/File> [<showFit>]")};
$galacticusFile = $ARGV[0];
$outputTo       = $ARGV[1];
if ( $#ARGV == 2 ) {
    $showFit    = $ARGV[2];
    if ( lc($showFit) eq "showfit"   ) {$showFit = 1};
    if ( lc($showFit) eq "noshowfit" ) {$showFit = 0};
} else {
    $showFit = 0;
}

# Check if output location is file or directory.
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/Stellar_Mass_Function.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataSet{'file'} = $galacticusFile;
$dataSet{'store'} = 0;
&HDF5::Get_Parameters(\%dataSet);
&HDF5::Count_Trees(\%dataSet);
&HDF5::Select_Output(\%dataSet,0.0);

# Read the XML data file.
$xml       = new XML::Simple;
$data      = $xml->XMLin("data/Stellar_Mass_Function_Li_White_2009.xml");
$columns   = $data->{'stellarMassFunction'}->{'columns'};
$xBins     = pdl @{$columns->{'stellarMass' }->{'datum'}};
$x         = pdl @{$columns->{'stellarMass' }->{'datum'}};
$y         = pdl @{$columns->{'massFunction'}->{'datum'}};
$errorUp   = pdl @{$columns->{'upperError'  }->{'datum'}};
$errorDown = pdl @{$columns->{'lowerError'  }->{'datum'}};
$errorUp   = (10.0**($y+$errorUp  ))*($dataSet{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;
$errorDown = (10.0**($y-$errorDown))*($dataSet{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;
$xBins     = $xBins+log10(           ($dataSet{'parameters'}->{'H_0'}/$columns->{'stellarMass' }->{'hubble'})**$columns->{'stellarMass' }->{'hubbleExponent'});
$x         = (10.0**$x             )*($dataSet{'parameters'}->{'H_0'}/$columns->{'stellarMass' }->{'hubble'})**$columns->{'stellarMass' }->{'hubbleExponent'} ;
$y         = (10.0**$y             )*($dataSet{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;

# Read galaxy data and construct mass function.
$dataSet{'tree'} = "all";
&HDF5::Get_Dataset(\%dataSet,['volumeWeight','diskStellarMass','spheroidStellarMass']);
$dataSets               = \%{$dataSet{'dataSets'}};
$logarithmicStellarMass = log10((${$dataSets->{'diskStellarMass'}}+${$dataSets->{'spheroidStellarMass'}}));
$weight                 = ${$dataSets->{'volumeWeight'}};
delete($dataSet{'dataSets'});
($yGalacticus,$errorGalacticus) = &Histograms::Histogram($xBins,$logarithmicStellarMass,$weight,differential => 1);

# Compute chi^2.
$chiSquared = sum(($yGalacticus-$y)**2/($errorGalacticus**2+(0.5*($errorUp-$errorDown))**2));
$degreesOfFreedom = nelem($y);
if ( $showFit == 1 ) {
    $fitData{'name'} = "Li and White (2009) stellar mass function";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make the plot.
$plot1  = Graphics::GnuplotIF->new();
$plot1->gnuplot_hardcopy( '| ps2pdf - '.$outputFile, 
			  'postscript enhanced', 
			  'color lw 3 solid' );
$plot1->gnuplot_set_xlabel("M_* [M_{{/=12 O}&{/*-.66 O}{/=12 \267}}]");
$plot1->gnuplot_set_ylabel("dn/dlog_{10}M_* [Mpc^{-3}]");
$plot1->gnuplot_set_title("Stellar Mass Function");
$plot1->gnuplot_cmd("set label \"{/Symbol c}^2=".FormatSigFigs($chiSquared,4)." [".$degreesOfFreedom."]\" at screen 0.2, screen 0.2");
$plot1->gnuplot_cmd("set logscale xy");
$plot1->gnuplot_cmd("set mxtics 10");
$plot1->gnuplot_cmd("set mytics 10");
$plot1->gnuplot_cmd("set format y \"10^{\%L}\"");
$plot1->gnuplot_cmd("set format x \"10^{\%L}\"");
$plot1->gnuplot_cmd("set pointsize 1.0");
$plot1->gnuplot_cmd("plot '-' with errorbars lt 1 pt 6 title \"".$data->{'stellarMassFunction'}->{'label'}."\", '-' with errorbars lt 2 pt 4 title \"Galacticus\"");
for ($i=0;$i<nelem($x);++$i) {
   $plot1->gnuplot_cmd($x->index($i)." ".$y->index($i)." ".$errorDown->index($i)." ".$errorUp->index($i));
}
$plot1->gnuplot_cmd("e");
for ($i=0;$i<nelem($x);++$i) {
   $plot1->gnuplot_cmd($x->index($i)." ".$yGalacticus->index($i)." ".$errorGalacticus->index($i));
}
$plot1->gnuplot_cmd("e");

exit;
