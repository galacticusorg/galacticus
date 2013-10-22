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
use Astro::Cosmology;
use Math::SigFigs;
use Data::Dumper;
require Galacticus::HDF5;
require Galacticus::StellarMass;
require Stats::Histograms;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;

# Compute likelihood (and make plots) for a Galacticus model given the stellar mass function data from Caputi et al. (2011;
# http://adsabs.harvard.edu/abs/2011MNRAS.413..162C).

# Get name of input and output files.
die("stellarMassFunction_Caputi_2011.pl <galacticusFile> <outputFile> <plotFile>") unless ( $#ARGV == 1 || $#ARGV == 2 );
my $self           = $0;
my $galacticusFile = $ARGV[0];
my $outputFile     = $ARGV[1];
my $plotFile;
$plotFile          = $ARGV[2] if ( $#ARGV == 2 );

# Read the XML data file.
my $xml       = new XML::Simple;
my $data      = $xml->XMLin($galacticusPath."data/observations/massFunctionsStellar/Stellar_Mass_Functions_Caputi_2011.xml");

# Initialize the constraint data.
my $constraint;

# Loop over redshifts.
foreach my $dataset ( @{$data->{'dataset'}} ) {

    # Determine the redshift.
    my $redshift = ($dataset->{'redshiftLow'}+$dataset->{'redshiftHigh'})/2.0;

    # Extract data arrays.
    my $xBins     = pdl @{$dataset->{'stellarMass' }->{'datum'}};
    my $x         = pdl @{$dataset->{'stellarMass' }->{'datum'}};
    my $y         = pdl @{$dataset->{'massFunction'}->{'datum'}};
    my $errorUp   = pdl @{$dataset->{'upperError'  }->{'datum'}};
    my $errorDown = pdl @{$dataset->{'lowerError'  }->{'datum'}};

    # Create data structure to read the results.
    my $galacticus;
    $galacticus->{'file' } = $galacticusFile;
    $galacticus->{'store'} = 0;
    &HDF5::Get_Parameters($galacticus          );
    &HDF5::Count_Trees   ($galacticus          );
    &HDF5::Select_Output ($galacticus,$redshift);

    # Compute cosmology corrections.
    my $cosmologyData = Astro::Cosmology->new(
	omega_matter => $data->{'cosmology'}->{'omega'},
	omega_lambda => $data->{'cosmology'}->{'lambda'},
	H0           => $data->{'cosmology'}->{'hubble'}
	);
    my $cosmologyGalacticus = Astro::Cosmology->new(
	omega_matter => $galacticus->{'parameters'}->{'Omega_Matter'},
	omega_lambda => $galacticus->{'parameters'}->{'Omega_DE'},
	H0           => $galacticus->{'parameters'}->{'H_0'}
	);

    my $volumeElementData            = $cosmologyData      ->differential_comoving_volume($redshift);
    my $volumeElementGalacticus      = $cosmologyGalacticus->differential_comoving_volume($redshift);
    my $luminosityDistanceData       = $cosmologyData      ->luminosity_distance         ($redshift);
    my $luminosityDistanceGalacticus = $cosmologyGalacticus->luminosity_distance         ($redshift);
    my $massCorrection               = ($luminosityDistanceGalacticus/$luminosityDistanceData)**2;
    my $volumeCorrection             = $volumeElementData/$volumeElementGalacticus;
    
    # Convert to Hubble constant of Galacticus model.
    $errorUp   = (10.0**($y+$errorUp  ))*      $volumeCorrection ;
    $errorDown = (10.0**($y-$errorDown))*      $volumeCorrection ;
    $xBins     =        $xBins          +log10($massCorrection  );
    $x         = (10.0**$x             )*      $massCorrection   ;
    $y         = (10.0**$y             )*      $volumeCorrection ;
    $errorUp   .= +$errorUp  -$y;
    $errorDown .= -$errorDown+$y;

    # Read galaxy data and construct mass function.
    $galacticus->{'tree'} = "all";
    &HDF5::Get_Dataset($galacticus,['volumeWeight','stellarMass']);
    my $dataSets               = $galacticus->{'dataSets'};
    my $logarithmicStellarMass = log10($dataSets->{'stellarMass' });
    my $weight                 =       $dataSets->{'volumeWeight'};
    delete($galacticus->{'dataSets'});

    # Construct the mass function. Assume 0.3 dex errors on stellar masses which is the dispersion in Delta M/M* reported by
    # Caputi et al. (2011).
    my $sigma = pdl ones(nelem($logarithmicStellarMass))*0.3;
    (my $yGalacticus,my $errorGalacticus) = &Histograms::Histogram($xBins,$logarithmicStellarMass,$weight,differential => 1,gaussianSmooth => $sigma);

    # Estimate reasonable errors for the stellar mass function. No covariance matrix is available.
    my $error              = 0.5*($errorUp+$errorDown);

    # Compute chi^2.
    my $chiSquared       = sum(($yGalacticus-$y)**2/($errorGalacticus**2+$error**2));
    my $degreesOfFreedom = nelem($y);
    my $logLikelihood    = -0.5*$chiSquared;
    $constraint->{'logLikelihood'} += $logLikelihood;
    
    # Create a plot of the luminosity function.
    if ( defined($plotFile) ) {
	
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	
	# Open a pipe to GnuPlot.
	(my $pRedshift = $redshift) =~ s/\./p/;
	(my $thisFile = $plotFile) =~ s/\.pdf$/_z$pRedshift.pdf/;
	($plotFileEPS = $thisFile) =~ s/\.pdf$/.eps/;
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
	print $gnuPlot "set xrange [1.0e10:1.0e12]\n";
	print $gnuPlot "set yrange [1.0e-7:1.0e-3]\n";
	print $gnuPlot "set title 'Stellar mass function at \$z=".$redshift."\$'\n";
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
				      title     => $data->{'label'}
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
	&MetaData::Write($thisFile,$galacticusFile,$self);
    }

}
    
# Output the constraint.
unless ( $outputFile eq "none" ) {
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$outputFile);
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

exit;
