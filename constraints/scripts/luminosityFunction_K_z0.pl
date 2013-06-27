#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use XML::Simple;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require Stats::Histograms;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;

# Get name of input and output files.
die("luminosityFunction_K_z0.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $self           = $0;
my $galacticusFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments;
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
my $outputFile;
}

# Create data structure to read the results.
my $model;
$model->{'file' } = $galacticusFile;
$model->{'store'} = 0;
&HDF5::Get_Parameters($model    );
&HDF5::Count_Trees   ($model    );
&HDF5::Select_Output ($model,0.0);

# Read the XML data file.
my $xml     = new XML::Simple;
my $data    = $xml->XMLin($galacticusPath."data/observations/luminosityFunctions/K_Luminosity_Function_Cole_2001.xml");
my $columns = $data->{'luminosityFunction'}->{'columns'};
my $xBins   = pdl @{$columns->{'magnitude'}->{'data'}};
my $x       = pdl @{$columns->{'magnitude'}->{'data'}};
my $y       = pdl @{$columns->{'luminosityFunction'}->{'data'}};
my $error   = pdl @{$columns->{'error'}->{'data'}};
$xBins      = $xBins-5.0*log10($columns->{'magnitude'}->{'hubble'}/$model->{'parameters'}->{'H_0'});
$x          = $x-5.0*log10($columns->{'magnitude'}->{'hubble'}/$model->{'parameters'}->{'H_0'});
$y          = $y*($model->{'parameters'}->{'H_0'}/$columns->{'luminosityFunction'}->{'hubble'})**3;
$error      = $error*($model->{'parameters'}->{'H_0'}/$columns->{'luminosityFunction'}->{'hubble'})**3;

# Reverse the order of the vectors.
$xBins = $xBins(-1:0);
$x     = $x    (-1:0);
$y     = $y    (-1:0);
$error = $error(-1:0);

# Read galaxy data and construct mass function.
my $xGalacticus = $xBins;
$model->{'tree'} = "all";
&HDF5::Get_Dataset($model,['mergerTreeWeight','magnitudeTotal:UKIRT_K:rest:z0.0000:dustAtlas:vega']);
my $dataSets  = $model->{'dataSets'};
my $magnitude = $dataSets->{'magnitudeTotal:UKIRT_K:rest:z0.0000:dustAtlas:vega'};
my $weight    = $dataSets->{'mergerTreeWeight'};
(my $yGalacticus, my $errorGalacticus) = &Histograms::Histogram($xGalacticus,$magnitude,$weight,differential => 1);

# Output the results to file if requested.
if ( exists($arguments{'resultFile'}) ) {
    my $results;
    @{$results->{'x'    }} = $xGalacticus    ->list();
    @{$results->{'y'    }} = $yGalacticus    ->list();
    @{$results->{'error'}} = $errorGalacticus->list();
    my $xmlOut = new XML::Simple (RootName=>"results", NoAttr => 1);;
    # Output the parameters to file.
    open(pHndl,">".$arguments{'resultFile'});
    print pHndl $xmlOut->XMLout($results);
    close pHndl;
}

# Output accuracy to file if requested.
if ( exists($arguments{'accuracyFile'}) ) {
    my $results;
    @{$results->{'x'         }} = $xGalacticus    ->list();
    @{$results->{'yModel'    }} = $yGalacticus    ->list();
    @{$results->{'yData'     }} = $y              ->list();
    @{$results->{'errorModel'}} = $errorGalacticus->list();
    @{$results->{'errorData' }} = $error          ->list();
    my $xmlOut = new XML::Simple (RootName=>"accuracy", NoAttr => 1);;
    # Output the parameters to file.
    open(pHndl,">".$arguments{'accuracyFile'});
    print pHndl $xmlOut->XMLout($results);
    close pHndl;
}

# Compute chi^2.
my $chiSquared       = sum(($yGalacticus-$y)**2/($errorGalacticus**2+$error**2));
my $degreesOfFreedom = nelem($y);
my $logLikelihood    = -0.5*$chiSquared;
my $constraint;
$constraint->{'logLikelihood'} = $logLikelihood;

# Output the constraint.
if ( exists($arguments{'outputFile'}) ) {
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

# Create a plot of the luminosity function.
if ( exists($arguments{'plotFile'}) ) {
    
    # Declare variables for GnuPlot;
    my ($gnuPlot, $plotFileEPS, $plot);
    
    # Open a pipe to GnuPlot.
    ($plotFileEPS =  $arguments{'plotFile'}) =~ s/\.pdf$/.eps/;
    open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
    print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
    print $gnuPlot "set output '".$plotFileEPS."'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.4,0.2\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale y\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [-28:-18]\n";
    print $gnuPlot "set yrange [1.0e-8:1.0e-2]\n";
    print $gnuPlot "set title 'K-band luminosity function at \$z=0\$'\n";
    print $gnuPlot "set xlabel '\$M_{\\rm K,Vega}\$ [dust-extinguished]'\n";
    print $gnuPlot "set ylabel '\${\\rm d}n/{\\rm d}M_{\\rm K,Vega}\$'\n";
    &PrettyPlots::Prepare_Dataset(\$plot,
				  $x,$y,
				  errorUp   => $error,
				  errorDown => $error,
				  style     => "point",
				  symbol    => [6,7], 
				  weight    => [5,3],
				  color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
				  title     => $data->{'luminosityFunction'}->{'label'}
	);
    &PrettyPlots::Prepare_Dataset(\$plot,
				  $xGalacticus,$yGalacticus,
				  errorUp   => $errorGalacticus,
				  errorDown => $errorGalacticus,
				  style     => "point",
				  symbol    => [6,7], 
				  weight    => [5,3],
				  color     => $PrettyPlots::colorPairs{'redYellow'},
				  title     => "Galacticus"
	);
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileEPS);
    &MetaData::Write($arguments{'plotFile'},$galacticusFile,$self);
}

exit;
