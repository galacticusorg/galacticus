#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
use Galacticus::HDF5;
use Galacticus::Options;
use Galacticus::Constraints::Covariances;
use Stats::Histograms;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Evaluate likelihood of a model compared to the spin distribution measured in the Millennium Simulation by Bett et al. (2007).
# Andrew Benson (01-August-2017)

# Get command line arguments.
die("Usage: spinDistributionFunctionBett2007.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName  =$ARGV[0];
# Create a hash of named arguments.
my %options;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);
# Read spin distribution function from model.
my $galacticusFile                     = new PDL::IO::HDF5($galacticusFileName);         
my $analyses                           = $galacticusFile->group  ('analyses'                          )       ;
my $spinAnalysis                       = $analyses      ->group  ('spinDistribution'                  )       ;
my $spinAnalysisSpin                   = $spinAnalysis  ->dataset('spin'                              )->get();
my $spinAnalysisDistribution           = $spinAnalysis  ->dataset('spinDistributionFunction'          )->get();
my $spinAnalysisDistributionCovariance = $spinAnalysis  ->dataset('spinDistributionFunctionCovariance')->get();
# Read the N-body spin distribution.
my $nbodyFile                  = new PDL::IO::HDF5(&galacticusPath().'data/darkMatter/bett2007HaloSpinDistribution.hdf5');
my $spinNBody                  = $nbodyFile->dataset('spinParameter'    )->get();
my $spinDistributionNBody      = $nbodyFile->dataset('distribution'     )->get();
my $spinDistributionErrorNBody = $nbodyFile->dataset('distributionError')->get();
# Convert N-body distribution to per natural log.
$spinDistributionNBody      /= log(10.0);
$spinDistributionErrorNBody /= log(10.0);
# Compute the likelihood.
if ( exists($options{'outputFile'}) ) {
    # For bins in the N-body distribution with zero halos, set an appropriate error.
    my $nbodyHaloCount = pdl 1503922.0; # Total number of halos in the TREEclean sample of Bett et al. (2007).
    my $nbodyBinWidth  = log($spinNBody->((1))/$spinNBody->((0)));
    my $nbodyZeroBins  = which($spinDistributionNBody <= 0.0);
    $spinDistributionErrorNBody->($nbodyZeroBins) .= 1.0/$nbodyHaloCount/$nbodyBinWidth;
    # Compute combined covariance.
    my $fullCovariance = $spinAnalysisDistributionCovariance+stretcher($spinDistributionErrorNBody**2);
    # Compute the likelihood.
    my $constraint;
    my $logDeterminant;
    my $offsets;
    my $inverseCovariance;
    my $jacobian;
    my $logLikelihood =
	&Galacticus::Constraints::Covariances::ComputeLikelihood
	(
	 $spinDistributionNBody                               ,
	 $spinAnalysisDistribution                            ,
	 $fullCovariance                                      ,
	 jacobian                     => \$jacobian           ,
	 offsets                      => \$offsets            ,
	 inversionMethod              => "eigendecomposition" ,
	 productMethod                => "linearSolver"
	);
    $constraint->{'label'        } = "spinDistributionFunctionBett2007";
    $constraint->{'logLikelihood'} = $logLikelihood;
    # Compute the variance in the log-likelihood due to errors in the model.
    my $logLikelihoodVariance = transpose($jacobian) x $spinAnalysisDistributionCovariance x $jacobian;
    $constraint->{'logLikelihoodVariance'} = $logLikelihoodVariance->sclr();
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$options{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}
# Plot the distribution.
if ( exists($options{'plotFile'}) ) {
    my $plot;
    my $gnuPlot;
    (my $plotFileTeX = $options{'plotFile'}) =~ s/\.pdf$/.tex/;
    open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set title offset 0,-1 'Spin distribution at \$z=0\$'\n";
    print $gnuPlot "set xlabel 'Spin parameter; \$\\lambda\\, []\$'\n";
    print $gnuPlot "set ylabel 'Distribution; \$\\mathrm{d}f/\\mathrm{d}\\ln \\lambda\\, []\$'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.175,0.80\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [1.0e-3:0.3]\n";
    print $gnuPlot "set yrange [1.0e-3:1.0e0]\n";
    print $gnuPlot "set pointsize 1.0\n";
     &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	$spinNBody,
	$spinDistributionNBody,
	errorDown  => $spinDistributionErrorNBody,
	errorUp    => $spinDistributionErrorNBody,
	style      => "point",
	symbol     => [6,7],
	weight     => [2,1],
	pointSize  => 0.2,
	color      => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	title      => "Bett et al. (2007)"
	);
   &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	$spinAnalysisSpin,
	$spinAnalysisDistribution+$spinAnalysisDistributionCovariance->diagonal(0,1)->sqrt(),
	y2         => $spinAnalysisDistribution-$spinAnalysisDistributionCovariance->diagonal(0,1)->sqrt(),
	style      => "filledCurve",
	weight     => [2,1],
	color      => $GnuPlot::PrettyPlots::colorPairs{'mediumSeaGreen'}
	);
   &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	$spinAnalysisSpin,
	$spinAnalysisDistribution,
	style      => "line",
	weight     => [2,1],
	color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	title      => "Galacticus"
	);
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
}

exit 0;
