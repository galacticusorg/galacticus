#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Galacticus::Path;
use PDL;
use PDL::NiceSlice;
use PDL::Math;
use PDL::IO::HDF5;
use PDL::LinearAlgebra;
use Galacticus::Options;
use Galacticus::Constraints::Covariances;
use Galacticus::HDF5;

# Compute likelihood (and make a plot) for a Galacticus model given the mass-metallicity relation data for z=0 from Blanc et
# al. (2017).

# Get name of input and output files.
die("massMetallicityRelation_Blanc2017.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Read observational data.
my  $observations                   = new PDL::IO::HDF5(&galacticusPath()."data/observations/abundances/massMetallicityRelationBlanc2017.hdf5");
my  $massStellarObserved            = $observations->dataset('massStellar'               )->get();
my  $metallicityObserved            = $observations->dataset('abundanceOxygenMean'       )->get();
my  $metallicity16PercentCIObserved = $observations->dataset('abundanceOxygen16PercentCI')->get();
my  $metallicity84PercentCIObserved = $observations->dataset('abundanceOxygen84PercentCI')->get();

# Estimate error on observed mean by taken 16-84% range as measure of the dispersion in the measurements, and assuming the same
# number of galaxies per bin.
my $galaxyCountTotal               = pdl 43306.0;                                   # Blanc et al. (2017)
my $galaxyCountBin                 = $galaxyCountTotal/nelem($massStellarObserved);
my $errorObserved                  = 0.5*($metallicity84PercentCIObserved+$metallicity16PercentCIObserved)/sqrt($galaxyCountBin);

# Read model data.
my $model;
$model->{'file'}    = $galacticusFileName;
&Galacticus::HDF5::Open_File     ($model);
&Galacticus::HDF5::Get_Parameters($model);
my $analysis        = $model   ->       {'hdf5File'                 }->group ('analyses'                )
                                                                     ->group ('massMetallicityBlanc2017');
my $massStellarModel = $analysis->dataset('massStellar'              )->get  (                          );
my $metallicityModel = $analysis->dataset('metallicityMean'          )->get  (                          );
my $covarianceModel  = $analysis->dataset('metallicityMeanCovariance')->get  (                          );
my $errorModel       = sqrt($covarianceModel->diagonal(0,1));

# Compute confidence intervals on data.

# Output the results to file if requested.
if ( exists($arguments{'resultFile'}) ) {
    my $resultsFile = new PDL::IO::HDF5(">".$arguments{'resultFile'});
    $resultsFile->dataset('x'             )->set($massStellarObserved);
    $resultsFile->dataset('y'             )->set($metallicityModel   );
    $resultsFile->dataset('covariance'    )->set($covarianceModel    );
    $resultsFile->dataset('yData'         )->set($metallicityObserved);
    $resultsFile->dataset('covarianceData')->set($errorObserved      );
}

# Compute the likelihood:
if ( exists($arguments{'outputFile'}) ) {
    # Build the full covariance matrix.
    my $covarianceObserved = stretcher($errorObserved);
    my $fullCovariance     = $covarianceModel+$covarianceObserved;
    # Compute the likelihood.
    my $offsets;
    my $jacobian;
    my $constraint;
    my $logLikelihood =
	&Galacticus::Constraints::Covariances::ComputeLikelihood
	(
	 $metallicityModel                           ,
	 $metallicityObserved                        ,
	 $fullCovariance                             ,
	 jacobian              => \$jacobian         ,
	 offsets               => \$offsets          ,
	 quiet                 => $arguments{'quiet'},
	 productMethod         => "linearSolver" 
	);
    $constraint->{'label'} = "massMetallicityRelationZ0.0Blanc2017";
    $constraint->{'logLikelihood'} = $logLikelihood;
    # Compute the variance in the log-likelihood due to errors in the model.
    my $logLikelihoodVariance = transpose($jacobian) x $covarianceModel x $jacobian;
    $constraint->{'logLikelihoodVariance'} = $logLikelihoodVariance->sclr();
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}
 
# Make a plot if requested.
if ( exists($arguments{'plotFile'}) ) {
    use GnuPlot::PrettyPlots;
    use GnuPlot::LaTeX;
    # Declare variables for GnuPlot;
    my ($gnuPlot, $plotFileTeX, $plot);
    # Open a pipe to GnuPlot.
    ($plotFileTeX = $arguments{'plotFile'}) =~ s/\.pdf$/.tex/;
    open($gnuPlot,"|gnuplot ");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.45,0.2\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale x\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [3.0e8:2.0e11]\n";
    print $gnuPlot "set yrange [8.1:9.0]\n";
    print $gnuPlot "set title offset 0,-0.9 'Mass-metallicity relation'\n";
    print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star\\ [\\mathrm{M}_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Metallicity; \$12+\\log_{10}(\\mathrm{O}/\\mathrm{H})\$\n";
    print $gnuPlot "set pointsize 0.5\n";
    &GnuPlot::PrettyPlots::Prepare_Dataset
	(
	 \$plot,
	 $massStellarObserved,
	 $metallicityObserved,
	 errorDown => $errorObserved,
	 errorUp   => $errorObserved,
	 style     => "point",
	 symbol    => [6,7], 
	 weight    => [1,1],
	 color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	 title     => "Blanc et al. (2017)"
    	);
    # Plot model.
    &GnuPlot::PrettyPlots::Prepare_Dataset
	(
	 \$plot,
	 $massStellarModel,
	 $metallicityModel,
	 errorDown => $errorModel,
	 errorUp   => $errorModel,
	 style     => "point",
	 symbol    => [6,7], 
	 weight    => [1,1],
	 color     => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	 title     => "Galacticus"
    	);
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX,margin => 1);
}

exit 0;
