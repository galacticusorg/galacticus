#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Galacticus::Path;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::LinearAlgebra;
use Galacticus::Options;
use Galacticus::HDF5;
use Galacticus::StellarMass;
use Galacticus::Constraints::Covariances;
use Stats::Means;

# Compute likelihood (and make a plot) for a Galacticus model given the black hole-bulge mass relation at z~0 from Kormendy & Ho
# (2013).

# Get name of input and output files.
die("blackHoleBulgeMassRelation_KormendyHo2013.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Define Galacticus unit system.
my $massSolar = pdl 1.989e30;

# Read observational data.
my  $observations                = new PDL::IO::HDF5(&galacticusPath()."data/observations/blackHoles/blackHoleMassVsBulgeMass_KormendyHo2013.hdf5");
my  $massBulgeObserved           = $observations->dataset('massBulge'         )->get    (           );
my  $massBulgeErrorObserved      = $observations->dataset('massBulgeError'    )->get    (           );
my  $massBlackHoleObserved       = $observations->dataset('massBlackHole'     )->get    (           );
my  $massBlackHoleErrorObserved  = $observations->dataset('massBlackHoleError')->get    (           );
(my $massBulgeUnits            ) = $observations->dataset('massBulge'         )->attrGet("unitsInSI");
(my $massBlackHoleUnits        ) = $observations->dataset('massBlackHole'     )->attrGet("unitsInSI");
(my $massBulgeScaling          ) = $observations->dataset('massBulge'         )->attrGet("scaling"  );
(my $massBlackHoleScaling      ) = $observations->dataset('massBlackHole'     )->attrGet("scaling"  );

# Convert mass to linear scalings.
if      ( $massBulgeScaling eq "linear" ) {
    # No conversion needed.
} elsif ( $massBulgeScaling eq "log10"  ) {
    # Convert from log10.
    $massBulgeObserved      .= 10.0**$massBulgeObserved;
    $massBulgeErrorObserved .= $massBulgeErrorObserved*log(10.0)*$massBulgeObserved;
} else {
    die('blackHoleMassDistribution_KormendyHo2013.pl: unknown scaling');
}
# Convert massBlackHole to linear scaling.
if      ( $massBlackHoleScaling eq "linear" ) {
    # Convert to log10.
    $massBlackHoleErrorObserved .= $massBlackHoleErrorObserved/$massBlackHoleObserved/log(10.0);
    $massBlackHoleObserved      .= log10($massBlackHoleObserved);
} elsif ( $massBlackHoleScaling eq "log10"  ) {
    # No conversion needed.
} else {
    die('blackHoleMassDistribution_KormendyHo2013.pl: unknown scaling');
}
# Convert to preferred unit system.
$massBulgeObserved          *=        $massBulgeUnits    /$massSolar ;
$massBulgeErrorObserved     *=        $massBulgeUnits    /$massSolar ;
$massBlackHoleObserved      +=  log10($massBlackHoleUnits/$massSolar);
$massBlackHoleErrorObserved +=  log10($massBlackHoleUnits/$massSolar);

# Read model black hole mass distribution.
my $model                        = new PDL::IO::HDF5($galacticusFileName);
my $analysisGroup                = $model         ->group  ('analyses'               )       ;
my $blackHoleGroup               = $analysisGroup ->group  ('blackHoleBulgeRelation' )       ;
my $massBulgeModel               = $blackHoleGroup->dataset('massStellarSpheroid'    )->get();
my $massBlackHoleModel           = $blackHoleGroup->dataset('massBlackHole'          )->get();
my $massBlackHoleCovarianceModel = $blackHoleGroup->dataset('massBlackHoleCovariance')->get();

# Compute binned mean of observed black hole masses.
my $weights = 1.0/$massBlackHoleErrorObserved**2;
(my $massBlackHoleMeanObserved, my $massBlackHoleMeanErrorObserved) = &Stats::Means::BinnedMean(log10($massBulgeModel),log10($massBulgeObserved),$massBlackHoleObserved,$weights);

# Compute likelihood if required.
if ( exists($arguments{'outputFile'}) ) {
     # Build the full covariance matrix.
    my $covarianceObserved = stretcher($massBlackHoleMeanErrorObserved**2);
    my $fullCovariance     = $massBlackHoleCovarianceModel+$covarianceObserved;
    # Compute the likelihood.
    my $offsets;
    my $jacobian;
    my $constraint;
    my $logLikelihood =
	&Galacticus::Constraints::Covariances::ComputeLikelihood
	(
	 $massBlackHoleModel                         ,
	 $massBlackHoleMeanObserved                  ,
	 $fullCovariance                             ,
	 jacobian              => \$jacobian         ,
	 offsets               => \$offsets          ,
	 quiet                 => $arguments{'quiet'},
	 productMethod         => "linearSolver" 
	);
    $constraint->{'logLikelihood'} = $logLikelihood;
    $constraint->{'label'        } = "blackHoleBulgeMassRelation_KormendyHo2013";
    # Compute the variance in the log-likelihood due to errors in the model.
    my $logLikelihoodVariance = transpose($jacobian) x $massBlackHoleCovarianceModel x $jacobian;
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
    print $gnuPlot "set key at screen 0.25,0.75\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [3.0e8:3.0e12]\n";
    print $gnuPlot "set yrange [3.0e5:1.0e11]\n";
    print $gnuPlot "set title offset 0,-0.9 'Black hole mass distribution'\n";
    print $gnuPlot "set xlabel 'Bulge stellar mass; \$M_{\\star, \\mathrm{bulge}}\\ [{\\rm M}_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Black hole mass; \$M_\\bullet\\ [{\\rm M}_\\odot]\$'\n";
    print $gnuPlot "set pointsize 0.5\n";
    # Plot points.
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	$massBulgeObserved,
	10.0**$massBlackHoleObserved,
	errorUp    => +10.0**(+$massBlackHoleErrorObserved+$massBlackHoleObserved)-10.0**$massBlackHoleObserved,
	errorDown  => -10.0**(-$massBlackHoleErrorObserved+$massBlackHoleObserved)+10.0**$massBlackHoleObserved,
	errorLeft  => $massBulgeErrorObserved,
	errorRight => $massBulgeErrorObserved,
	style      => "point",
	symbol     => [6,7],
	weight     => [1,1],
	pointSize  => 0.1,
	color      => $GnuPlot::PrettyPlots::colorPairs{'peachPuff'},
	title      => "Kormendy \\\\& Ho (2013); individual"
     	);
     &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	$massBulgeModel,
	10.0**$massBlackHoleMeanObserved,
	errorUp   => +10.0**(+$massBlackHoleMeanErrorObserved+$massBlackHoleMeanObserved)-10.0**$massBlackHoleMeanObserved,
	errorDown => -10.0**(-$massBlackHoleMeanErrorObserved+$massBlackHoleMeanObserved)+10.0**$massBlackHoleMeanObserved,
	style     => "point",
	symbol    => [6,7],
	weight    => [3,1],
	color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	title      => "Kormendy \\\\& Ho (2013); mean"
     	);
      &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot,
	$massBulgeModel,
	10.0**$massBlackHoleModel,
	errorUp   => $massBlackHoleCovarianceModel->diagonal(0,1)->sqrt(),
	errorDown => $massBlackHoleCovarianceModel->diagonal(0,1)->sqrt(),
	style     => "point",
	symbol    => [6,7],
	weight    => [3,1],
	color     => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	title      => "Galacticus"
     	);
   &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);    
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX,margin => 1);
}

exit;
