#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Fit::Polynomial;
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
use Galacticus::HDF5;
use Galacticus::Options;
use Galacticus::Path;
use Galacticus::Constraints::Covariances;
use Stats::Histograms;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Evaluate likelihood of a model compared to the concentration-mass relation measured in the Copernicus Complexio CDM simulation by Ludlow et al. (2016).
# Andrew Benson (11-August-2017)

# Get command line arguments.
die("Usage: concentrationMassRelationCDMLudlow2016.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName  =$ARGV[0];
# Create a hash of named arguments.
my %options;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);
# Read concentration-mass relation from model.
my $galacticusFile                                          = new PDL::IO::HDF5($galacticusFileName);         
my $analyses                                                = $galacticusFile       ->group  ('analyses'                                  )       ;
my $concentrationAnalysis                                   = $analyses             ->group  ('concentrationHaloMassRelationCDMLudlow2016')       ;
my $concentrationAnalysisMassHalo                           = $concentrationAnalysis->dataset('massHalo'                                  )->get();
my $concentrationAnalysisConcentrationLogarithmic           = $concentrationAnalysis->dataset('log10Concentration'                        )->get();
my $concentrationAnalysisConcentrationLogarithmicCovariance = $concentrationAnalysis->dataset('log10ConcentrationCovariance'              )->get();
# Read the N-body concentration-mass relation.
my $nbodyFile          = new PDL::IO::HDF5(&galacticusPath().'data/darkMatter/concentrationMassRelationCDMLudlow2016.hdf5');
my $massHaloNBody      = $nbodyFile->dataset('massHalo'         )->get();
my $concentrationNBody = $nbodyFile->dataset('concentrationMean')->get();
# Convert N-body concentration to logarithmic form and estimate errors on them. We do this by fitting a 3rd order polynomial to
# the N-body points and considering the residuals to be error estimates. The squares of the residuals are then fit with a 1st
# order polynomial to predict the variance.
my  $concentrationLogarithmicNBody          = log10($concentrationNBody);
(my $concentrationLogarithmicNBodyFit)      = fitpoly1d(log10($massHaloNBody),$concentrationLogarithmicNBody,3);
my $concentrationLogarithmicNBodyResidual   = $concentrationLogarithmicNBody-$concentrationLogarithmicNBodyFit;
(my $concentrationLogarithmicNBodyVariance) = fitpoly1d(log10($massHaloNBody),log10($concentrationLogarithmicNBodyResidual**2),2);
# Compute the likelihood.
if ( exists($options{'outputFile'}) ) {
    # Compute combined covariance.
    my $fullCovariance = $concentrationAnalysisConcentrationLogarithmicCovariance+stretcher(10.0**$concentrationLogarithmicNBodyVariance);
    # Compute the likelihood.
    my $constraint;
    my $logDeterminant;
    my $offsets;
    my $inverseCovariance;
    my $jacobian;
    my $logLikelihood =
	&Galacticus::Constraints::Covariances::ComputeLikelihood
	(
	 $concentrationLogarithmicNBody                       ,
	 $concentrationAnalysisConcentrationLogarithmic       ,
	 $fullCovariance                                      ,
	 jacobian                     => \$jacobian           ,
	 offsets                      => \$offsets            ,
	 inversionMethod              => "eigendecomposition" ,
	 productMethod                => "linearSolver"
	);
    $constraint->{'label'        } = "concentrationMassRelationCDMLudlow2016";
    $constraint->{'logLikelihood'} = $logLikelihood;
    # Compute the variance in the log-likelihood due to errors in the model.
    my $logLikelihoodVariance = transpose($jacobian) x $concentrationAnalysisConcentrationLogarithmicCovariance x $jacobian;
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
    print $gnuPlot "set title offset 0,-1 'Concentration-mass relation at \$z=0\$'\n";
    print $gnuPlot "set xlabel 'Halo mass; \$M_\\mathrm{200c}\\, [\\mathrm{M}_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Mean log concentration; \$\\langle \\log_{10} c_\\mathrm{200c}\\rangle\\, []\$'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.5,0.80\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale x\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [1.0e8:1.0e13]\n";
    print $gnuPlot "set yrange [0.8:1.2]\n";
    print $gnuPlot "set pointsize 1.0\n";
     &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	$massHaloNBody,
	$concentrationLogarithmicNBody,
	errorDown  => 10.0**(0.5*$concentrationLogarithmicNBodyVariance),
	errorUp    => 10.0**(0.5*$concentrationLogarithmicNBodyVariance),
	style      => "point",
	symbol     => [6,7],
	weight     => [2,1],
	color      => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	title      => "Ludlow et al. (2016)"
	);
   &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	$concentrationAnalysisMassHalo,
	$concentrationAnalysisConcentrationLogarithmic,
	errorUp   => $concentrationAnalysisConcentrationLogarithmicCovariance->diagonal(0,1)->sqrt(),
	errorDown => $concentrationAnalysisConcentrationLogarithmicCovariance->diagonal(0,1)->sqrt(),
	style      => "point",
	weight     => [3,1],
	symbol     => [6,7],
	color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	title      => "Galacticus"
	);
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
}

exit 0;
