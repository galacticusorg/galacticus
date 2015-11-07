#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::LinearAlgebra;
require Galacticus::Options;
require Galacticus::HDF5;
require Galacticus::StellarMass;
require Galacticus::Constraints::Covariances;

# Compute likelihood (and make a plot) for a Galacticus model given the stellar mass-gas-phase metallicity relation data for z=0 from
# Andrews & Martini (2013).

# Get name of input and output files.
die("stellarMassMetallicityRelation_Andrews_2013.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Define Galacticus unit system.
my $massSolar        = pdl 1.989e30;
my $metallicitySolar = pdl 10.0**8.86;

# Read observational data.
my  $observations                   = new PDL::IO::HDF5($galacticusPath."data/observations/abundances/gasPhaseMetallicityAndrews2013.hdf5");
my  $massStellarObserved            = $observations->dataset('mass'                 )->get    (           );
my  $metallicityMeanObserved        = $observations->dataset('metallicity'          )->get    (           );
my  $metallicityCovarianceObserved  = $observations->dataset('metallicityCovariance')->get    (           );
(my $massUnits                    ) = $observations->dataset('mass'                 )->attrGet("unitsInSI");
(my $metallicityUnits             ) = $observations->dataset('metallicity'          )->attrGet("unitsInSI");
(my $massScaling                  ) = $observations->dataset('mass'                 )->attrGet("scaling"  );
(my $metallicityScaling           ) = $observations->dataset('metallicity'          )->attrGet("scaling"  );

# Extract errors.
my $massCount                = nelem($massStellarObserved);
my $metallicityErrorObserved = sqrt($metallicityCovarianceObserved->diagonal(0,1));

# Convert mass to linear scalings.
if      ( $massScaling eq "linear" ) {
    # No conversion needed.
} elsif ( $massScaling eq "log10"  ) {
    # Convert from log10.
    $massStellarObserved .= 10.0**$massStellarObserved;
} else {
    die('stellarMassMetallicityRelation_Andrews_2013.pl: unknown scaling');
}
# Convert metallicity to linear scaling.
if      ( $metallicityScaling eq "linear" ) {
    # No conversion needed.
} elsif ( $metallicityScaling eq "log10"  ) {
    # Convert from log10.
    $metallicityMeanObserved       .= 10.0**$metallicityMeanObserved;
    $metallicityErrorObserved      .= $metallicityErrorObserved     *log(10.0)   *      $metallicityMeanObserved                          ;
    $metallicityCovarianceObserved .= $metallicityCovarianceObserved*log(10.0)**2*outer($metallicityMeanObserved,$metallicityMeanObserved);
} else {
    die('stellarMassMetallicityRelation_Andrews_2013.pl: unknown scaling');
}
# Convert to preferred unit system.
$massStellarObserved           *=  $massUnits       /$massSolar           ;
$metallicityMeanObserved       *=  $metallicityUnits/$metallicitySolar    ;
$metallicityErrorObserved      *=  $metallicityUnits/$metallicitySolar    ;
$metallicityCovarianceObserved *= ($metallicityUnits/$metallicitySolar)**2;

# Read model metallicity relation.
my $model                        = new PDL::IO::HDF5($galacticusFileName);
my $analysisGroup                = $model           ->group  ('analysis'                         )       ;
my $metallicityGroup             = $analysisGroup   ->group  ('sdssGasMetallicityZ0.07'          )       ;
my $massStellarModel             = $metallicityGroup->dataset('mass'                             )->get();
my $metallicityDistributionModel = $metallicityGroup->dataset('metallicityDistribution'          )->get();
my $metallicityCovarianceModel   = $metallicityGroup->dataset('metallicityDistributionCovariance')->get();
my $metallicityMeanModel         = $metallicityDistributionModel->((0),:);
my $metallicityMeanErrorModel    = sqrt($metallicityCovarianceModel->diagonal(0,1));

# Compute likelihood if required.
if ( exists($arguments{'outputFile'}) ) {
    # Construct the full covariance matrix, which is the covariance matrix of the observations
    # plus that of the model.
    my $fullCovariance = $metallicityCovarianceObserved+$metallicityCovarianceModel;
    # Compute the likelihood.
    my $constraint;
    my $logDeterminant;
    my $offsets;
    my $inverseCovariance;    
    my $logLikelihood = &Covariances::ComputeLikelihood($metallicityMeanObserved,$metallicityMeanModel,$fullCovariance, determinant => \$logDeterminant, inverseCovariance => \$inverseCovariance, offsets => \$offsets, quiet => $arguments{'quiet'}, inversionMethod => "eigendecomposition");
    $constraint->{'logLikelihood'} = $logLikelihood;
    # Find the Jacobian of the log-likelihood with respect to the model mass function.
    my $jacobian = pdl zeroes(1,nelem($metallicityMeanModel));
    for(my $i=0;$i<nelem($metallicityMeanModel);++$i) {
	$jacobian->((0),($i)) .= sum($inverseCovariance->(($i),:)*$offsets);
    }
    # Compute the variance in the log-likelihood due to errors in the model.
    my $logLikelihoodVariance = transpose($jacobian) x $metallicityCovarianceModel x $jacobian;
    $constraint->{'logLikelihoodVariance'} = $logLikelihoodVariance->sclr();
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

# Make a plot if requested.
if ( exists($arguments{'plotFile'}) ) {
    require GnuPlot::PrettyPlots;
    require GnuPlot::LaTeX;
    # Declare variables for GnuPlot;
    my ($gnuPlot, $plotFileEPS, $plot);
    # Open a pipe to GnuPlot.
    ($plotFileEPS = $arguments{'plotFile'}) =~ s/\.pdf$/.eps/;
    open($gnuPlot,"|gnuplot ");
    print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
    print $gnuPlot "set output '".$plotFileEPS."'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.45,0.2\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set ytics ('7.6' 3.98107e+07, '7.8' 6.30957e+07, '8.0' 100000000, '8.2' 1.58489e+08, '8.4' 2.51189e+08, '8.6' 3.98107e+08, '8.8' 6.30957e+08, '9.0' 1e+09, '9.2' 1.58489e+09, '9.4' 2.51189e+09)\n";
    print $gnuPlot "set xrange [1.0e7:1.0e11]\n";
    print $gnuPlot "set yrange [4.0e7:1.0e9]\n";
    print $gnuPlot "set title offset 0,-0.5 'Gas-phase mass-metallicity relation'\n";
    print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star\\ [{\\rm M}_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Metallicity; \$12 + \\log_{10}(\\hbox{O/H})\$'\n";
    print $gnuPlot "set pointsize 1.0\n";
    &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massStellarObserved,
     				  $metallicityMeanObserved,
				  errorUp   => $metallicityErrorObserved,
				  errorDown => $metallicityErrorObserved,
				  style     => "point",
				  symbol    => [6,7],
     				  weight    => [3,1],
     				  color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
     				  title     => "Andrews \\\\& Martini (2013)"
     	);
    &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massStellarModel,
     				  $metallicityMeanModel,
				  errorUp   => $metallicityMeanErrorModel,
				  errorDown => $metallicityMeanErrorModel,
				  style     => "point",
				  symbol    => [6,7],
     				  weight    => [3,1],
     				  color     => $PrettyPlots::colorPairs{'mediumSeaGreen'},
     				  title     => "Galacticus"
     	);
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
}

exit;
