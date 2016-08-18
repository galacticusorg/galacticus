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
use Galacticus::HDF5;

# Compute likelihood (and make a plot) for a Galacticus model given the early-type relation data
# for z=0 from Kelvin et al. (2014).

# Get name of input and output files.
die("earlyTypeFraction_Kelvin2014.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Specify constants.
my $massSolar = pdl 1.9891e30;

# Read observational data.
my  $observations              = new PDL::IO::HDF5(&galacticusPath()."data/observations/morphology/earlyTypeFractionGAMA.hdf5");
my  $massStellarObserved       = $observations->dataset('mass'             )->get();
my  $countAllObserved          = $observations->dataset('countAll'         )->get();
my  $countEarlyObserved        = $observations->dataset('countEarly'       )->get();
my  $countEarlyLBSObserved     = $observations->dataset('countEarlyPlusLBS')->get();
(my $labelObserved           ) = $observations->attrGet('label'            )       ;
my  $fractionEarlyObserved     = float($countEarlyObserved   )/float($countAllObserved);
my  $fractionEarlyLBSObserved  = float($countEarlyLBSObserved)/float($countAllObserved);

# Perform scaling and unit conversion.
(my $massUnitsInSI, my $massScaling) = $observations->dataset('mass')->attrGet('unitsInSI','scaling');
if ( $massScaling eq "linear" ) {
    # Nothing to do.
} elsif ( $massScaling eq "log10" ) {
    $massStellarObserved .= 10.0**$massStellarObserved;
} else {
    die("earlyTypeFraction_Kelvin2014.pl: unknown mass scaling");
}

# Read model data.
my $model;
$model->{'file'}    = $galacticusFileName;
&Galacticus::HDF5::Open_File     ($model);
&Galacticus::HDF5::Get_Parameters($model);
my $analysis        = $model   ->       {'hdf5File'                 }->group ('analysis'                  )
                                                                     ->group ('gamaEarlyTypeFractionZ0.03');
my $massModel       = $analysis->dataset('mass'                      )->get  (                            );
my $fractionModel   = $analysis->dataset('fractionFunction'          )->get  (                            );
my $covarianceModel = $analysis->dataset('fractionFunctionCovariance')->get  (                            );
my $errorModel      = sqrt($covarianceModel->diagonal(0,1));

# Extract probability of LBS class being included in early-type class.
my $lbsInclusionProbability = pdl 0.1;
$lbsInclusionProbability .= $model->{'parameter'}->{'gamaEarlyTypeLBSProbability'}->{'value'}
    if ( exists($model->{'parameter'}->{'gamaEarlyTypeLBSProbability'}) );

# Compute confidence intervals on data. (Follows the approach described in Cameron 2011; PASA; 28; 128;
# http://adsabs.harvard.edu/abs/2011PASA...28..128C)
my $confidenceLevel       = 0.683; # i.e. 1 sigma
my $alpha                 = 1.0-$confidenceLevel;
my $cRange                = pdl [ 0.5*$alpha, 1.0-0.5*$alpha ];
my $a                     =                   +$countEarlyObserved   +1;
my $aLBS                  =                   +$countEarlyLBSObserved+1;
my $b                     = +$countAllObserved-$countEarlyObserved   +1;
my $bLBS                  = +$countAllObserved-$countEarlyLBSObserved+1;
my $errorUpperObserved    = pdl zeroes(nelem($massStellarObserved));
my $errorUpperLBSObserved = pdl zeroes(nelem($massStellarObserved));
my $errorLowerObserved    = pdl zeroes(nelem($massStellarObserved));
my $errorLowerLBSObserved = pdl zeroes(nelem($massStellarObserved));
my $intervalCount      = 1000;
for(my $i=0;$i<nelem($massStellarObserved);++$i) {
    my $p    = pdl sequence($intervalCount)/($intervalCount-1.0);
    my $f    = $p**($a   ->(($i))-1.0)*(1.0-$p)**($b   ->(($i))-1.0);
    my $fLBS = $p**($aLBS->(($i))-1.0)*(1.0-$p)**($bLBS->(($i))-1.0);
    my $c    = $f   ->cumusumover()/$f   ->sum();
    my $cLBS = $fLBS->cumusumover()/$fLBS->sum();
    (my $pRange   , my $pRangeError   ) = interpolate($cRange,$c   ,$p);
    (my $pRangeLBS, my $pRangeErrorLBS) = interpolate($cRange,$cLBS,$p);
    $errorLowerObserved   ->(($i)) .= -$pRange   ->((0))+$fractionEarlyObserved   ->(($i));
    $errorLowerLBSObserved->(($i)) .= -$pRangeLBS->((0))+$fractionEarlyLBSObserved->(($i));
    $errorUpperObserved   ->(($i)) .= +$pRange   ->((1))-$fractionEarlyObserved   ->(($i));
    $errorUpperLBSObserved->(($i)) .= +$pRangeLBS->((1))-$fractionEarlyLBSObserved->(($i));
}
$errorUpperObserved   ->($errorUpperObserved   <0.0;?) .= 0.0;
$errorUpperLBSObserved->($errorUpperLBSObserved<0.0;?) .= 0.0;
$errorLowerObserved   ->($errorLowerObserved   <0.0;?) .= 0.0;
$errorLowerLBSObserved->($errorLowerLBSObserved<0.0;?) .= 0.0;


# Output the results to file if requested.
if ( exists($arguments{'resultFile'}) ) {
    my $resultsFile = new PDL::IO::HDF5(">".$arguments{'resultFile'});
    $resultsFile->dataset('x'             )->set($massStellarObserved  );
    $resultsFile->dataset('y'             )->set($fractionModel        );
    $resultsFile->dataset('covariance'    )->set($covarianceModel      );
    $resultsFile->dataset('yData'         )->set($fractionEarlyObserved);
    $resultsFile->dataset('covarianceData')->set(stretcher(0.5*($errorUpperObserved+$errorLowerObserved)));
}

# Compute the likelihood:
if ( exists($arguments{'outputFile'}) ) {
    # Find a Cholesky decomposition of the model covariance matrix for later construction of
    # model realizations.
    my $modelCholesky = mchol($covarianceModel);
    # Compute factors needed in likelihood function.
    (my $facA    )    =lgamma(float($a         ));
    (my $facALBS )    =lgamma(float($aLBS      ));
    (my $facB    )    =lgamma(float(      $b   ));
    (my $facBLBS )    =lgamma(float(      $bLBS));
    (my $facAB   )    =lgamma(float($a   +$b   ));
    (my $facABLBS)    =lgamma(float($aLBS+$bLBS));
    # Compute a base model likelihood.
    my $logLikelihoodsBase = &logLikelihood($a,$b,$facA,$facB,$facAB,$fractionModel);
    # Compute model realizations.
    my $realizationCount   =     100000;
    my $probability        = pdl      0.0;
    my $probabilitySquared = pdl      0.0;
    srand(10);
    for(my $i=0;$i<$realizationCount;++$i) {
	# Construct model realization.
	my $normalDeviates           = grandom(nelem($fractionModel));
	my $offsets                  = $normalDeviates x $modelCholesky;
	my $fractionModelRealization = $fractionModel+$offsets->clump(2);
	$fractionModelRealization->($fractionModelRealization < 0.001;?) .= 0.001;
	$fractionModelRealization->($fractionModelRealization > 0.999;?) .= 0.999;
	# Evaluate likelihood with and without LBS class included.
	my $logLikelihoods     = &logLikelihood($a   ,$b   ,$facA   ,$facB   ,$facAB   ,$fractionModelRealization)-$logLikelihoodsBase;
	my $logLikelihoodsLBS  = &logLikelihood($aLBS,$bLBS,$facALBS,$facBLBS,$facABLBS,$fractionModelRealization)-$logLikelihoodsBase;	
	# Compute the net probability.
	my $probabilityRealization = 
	    +     $lbsInclusionProbability *exp(sum($logLikelihoodsLBS))
	    +(1.0-$lbsInclusionProbability)*exp(sum($logLikelihoods   ));
	$probability        += $probabilityRealization   ;
	$probabilitySquared += $probabilityRealization**2;
    }
    # Normalize the probability, convert to log and add the base log-likelihood back in.
    $probability            /= $realizationCount;
    my $probabilityVariance  = ($probabilitySquared-$realizationCount*$probability**2)/($realizationCount-1);
    my $logLikelihood         = log($probability)+sum($logLikelihoodsBase);
    my $logLikelihoodVariance = $probabilityVariance/$probability**2;
    my $constraint;
    $constraint->{'logLikelihood'        } = $logLikelihood        ->sclr();
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
    print $gnuPlot "set logscale x\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [5.0e8:6.0e11]\n";
    print $gnuPlot "set yrange [-0.05:1.05]\n";
    print $gnuPlot "set title offset 0,-0.5 'Early-type fraction'\n";
    print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star\\ [{\\rm M}_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Early-type fraction; \$f\\ []\$'\n";
    print $gnuPlot "set pointsize 1.0\n";
    # Plot observations. Areas of points for observations with and without LBS class included
    # are scaled by the probability for LBS inclusion. Error bars showing the 68% confidence
    # interval are displayed.
    &GnuPlot::PrettyPlots::Prepare_Dataset
	(
	 \$plot,
	 $massStellarObserved,
	 $fractionEarlyLBSObserved,
	 errorDown => $errorLowerLBSObserved,
	 errorUp   => $errorUpperLBSObserved,
	 style     => "point",
	 symbol    => [6,7], 
	 weight    => [1,1],
	 pointSize => sqrt($lbsInclusionProbability),
	 color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'}
    	);
    &GnuPlot::PrettyPlots::Prepare_Dataset
	(
	 \$plot,
	 $massStellarObserved,
	 $fractionEarlyObserved,
	 errorDown => $errorLowerObserved,
	 errorUp   => $errorUpperObserved,
	 style     => "point",
	 symbol    => [6,7], 
	 weight    => [1,1],
	 pointSize => sqrt(1.0-$lbsInclusionProbability),
	 color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	 title     => $labelObserved
    	);
    # Plot model.
    &GnuPlot::PrettyPlots::Prepare_Dataset
	(
	 \$plot,
	 $massModel,
	 $fractionModel,
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
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
}

exit;

sub logLikelihood {
    # Evaluate the log likelihood using a binomial distribution.
    my $a        = shift();
    my $b        = shift();
    my $facA     = shift();
    my $facB     = shift();
    my $facAB    = shift();
    my $fraction = shift();
    return ($a-1.0)*log($fraction)+($b-1.0)*log(1.0-$fraction)+$facAB-$facA-$facB;
}
