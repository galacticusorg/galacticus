#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V091"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
require Galacticus::HDF5;
require Galacticus::StellarMass;
require Stats::Histograms;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;

# Compute likelihood (and make a plot) for a Galacticus model given the stellar mass function data from Li & White (2009;
# http://adsabs.harvard.edu/abs/2011ApJ...728...46W).

# Get name of input and output files.
die("stellarMassFunction_z0.07.pl <galacticusFile> <outputFile> <plotFile>") unless ( $#ARGV == 1 || $#ARGV == 2 );
my $self           = $0;
my $galacticusFile = $ARGV[0];
my $outputFile     = $ARGV[1];
my $plotFile;
$plotFile          = $ARGV[2] if ( $#ARGV == 2 );

# Read the XML data file.
my $xml       = new XML::Simple;
my $data      = $xml->XMLin($galacticusPath."data/observations/massFunctionsStellar/Stellar_Mass_Function_Li_White_2009.xml");
my $columns   = $data->{'stellarMassFunction'}->{'columns'};
my $xBins     = pdl @{$columns->{'stellarMass' }->{'datum'}};
my $x         = pdl @{$columns->{'stellarMass' }->{'datum'}};
my $y         = pdl @{$columns->{'massFunction'}->{'datum'}};
my $errorUp   = pdl @{$columns->{'upperError'  }->{'datum'}};
my $errorDown = pdl @{$columns->{'lowerError'  }->{'datum'}};

# Create data structure to read the results.
my $galacticus;
$galacticus->{'file' } = $galacticusFile;
$galacticus->{'store'} = 0;
&HDF5::Get_Parameters($galacticus     );
&HDF5::Count_Trees   ($galacticus     );
&HDF5::Select_Output ($galacticus,0.07);

# Convert to Hubble constant of Galacticus model.
$errorUp   = (10.0**($y+$errorUp  ))*($galacticus->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;
$errorDown = (10.0**($y-$errorDown))*($galacticus->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;
$xBins     = $xBins+log10(           ($galacticus->{'parameters'}->{'H_0'}/$columns->{'stellarMass' }->{'hubble'})**$columns->{'stellarMass' }->{'hubbleExponent'});
$x         = (10.0**$x             )*($galacticus->{'parameters'}->{'H_0'}/$columns->{'stellarMass' }->{'hubble'})**$columns->{'stellarMass' }->{'hubbleExponent'} ;
$y         = (10.0**$y             )*($galacticus->{'parameters'}->{'H_0'}/$columns->{'massFunction'}->{'hubble'})**$columns->{'massFunction'}->{'hubbleExponent'} ;
$errorUp   .= +$errorUp  -$y;
$errorDown .= -$errorDown+$y;

# Read galaxy data and construct mass function.
$galacticus->{'tree'} = "all";
&HDF5::Get_Dataset($galacticus,['volumeWeight','stellarMass']);
my $dataSets               = $galacticus->{'dataSets'};
my $logarithmicStellarMass = log10($dataSets->{'stellarMass' });
my $weight                 =       $dataSets->{'volumeWeight'};
delete($galacticus->{'dataSets'});

# Construct the mass function. Assume 0.2 dex errors on stellar masses. This is approximate, but motivated by the discussion of
# Bell et al. (2003; ApJS; 149; 289-312).
my $sigma = pdl ones(nelem($logarithmicStellarMass))*0.2;
(my $yGalacticus,my $errorGalacticus) = &Histograms::Histogram($xBins,$logarithmicStellarMass,$weight,differential => 1,gaussianSmooth => $sigma);

# Estimate reasonable errors for the stellar mass function. The errors quoted by Li & White were derived from mock catalogs drawn
# from the Millennium Simulation. Unfortunately, they didn't quote the covariance matrix, nor do they account for possible
# systematics due to stellar population synthesis code or IMF uncertainties. Therefore we simply prevent the errors from dropping
# below 20% (as suggested by Yang et al. (2012; arXiv:1110.1420).
my $error              = 0.5*($errorUp+$errorDown);
my $errorMinimum       = 0.2*$y;
my $limit              = which($error < $errorMinimum);
$error->index($limit) .= $errorMinimum->index($limit);

# Compute chi^2.
my $chiSquared       = sum(($yGalacticus-$y)**2/($errorGalacticus**2+$error**2));
my $degreesOfFreedom = nelem($y);
my $logLikelihood    = -0.5*$chiSquared;
my $constraint;
$constraint->{'logLikelihood'} = $logLikelihood;

# Output the constraint.
unless ( $outputFile eq "none" ) {
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$outputFile);
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

# Create a plot of the luminosity function.
if ( defined($plotFile) ) {
    
    # Declare variables for GnuPlot;
    my ($gnuPlot, $plotFileEPS, $plot);
    
    # Open a pipe to GnuPlot.
    ($plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
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
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [1.0e8:1.0e13]\n";
    print $gnuPlot "set yrange [1.0e-9:1.0e-1]\n";
    print $gnuPlot "set title 'Stellar mass function at \$z=0.07\$'\n";
    print $gnuPlot "set xlabel '\$M_\\star\$ [\$M_\\odot\$]'\n";
    print $gnuPlot "set ylabel '\${\\rm d}n/{\\rm d}\\log_{10}M_\\star\$'\n";
    &PrettyPlots::Prepare_Dataset(\$plot,
				  $x,$y,
				  errorUp   => $errorUp,
				  errorDown => $errorDown,
				  style     => "point",
				  symbol    => [6,7], 
				  weight    => [5,3],
				  color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
				  title     => $data->{'stellarMassFunction'}->{'label'}
	);
    &PrettyPlots::Prepare_Dataset(\$plot,
				  $x,$yGalacticus,
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
    &MetaData::Write($plotFile,$galacticusFile,$self);
}

exit;
