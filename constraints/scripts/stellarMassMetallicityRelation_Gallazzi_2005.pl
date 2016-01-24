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

# Compute likelihood (and make a plot) for a Galacticus model given the stellar mass-metallicity relation data for z=0 from
# Gallazzi et al. (2005).

# Get name of input and output files.
die("stellarMassMetallicityRelation_Gallazzi_2005.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
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
my $metallicitySolar = pdl 0.0189;

# Read observational data.
my  $observations                   = new PDL::IO::HDF5($galacticusPath."data/observations/abundances/stellarPhaseMetallicityGallazzi2005.hdf5");
my  $massStellarObserved            = $observations->dataset('mass'                   )->get    (           );
my  $metallicityObserved16          = $observations->dataset('metallicityPercentile16')->get    (           );
my  $metallicityObserved50          = $observations->dataset('metallicityPercentile50')->get    (           );
my  $metallicityObserved84          = $observations->dataset('metallicityPercentile84')->get    (           );
my  $metallicityCovarianceObserved  = $observations->dataset('metallicityCovariance'  )->get    (           );
(my $massUnits                    ) = $observations->dataset('mass'                   )->attrGet("unitsInSI");
(my $metallicityUnits             ) = $observations->dataset('metallicityPercentile50')->attrGet("unitsInSI");
(my $massScaling                  ) = $observations->dataset('mass'                   )->attrGet("scaling"  );
(my $metallicityScaling           ) = $observations->dataset('metallicityPercentile16')->attrGet("scaling"  );

# Extract errors.
my $massCount                  = nelem($massStellarObserved);
my $metallicityErrorObserved16 = sqrt($metallicityCovarianceObserved->diagonal(0,1)->(0*$massCount:1*$massCount-1));
my $metallicityErrorObserved50 = sqrt($metallicityCovarianceObserved->diagonal(0,1)->(1*$massCount:2*$massCount-1));
my $metallicityErrorObserved84 = sqrt($metallicityCovarianceObserved->diagonal(0,1)->(2*$massCount:3*$massCount-1));

# Create single metallicity vector.
my $percentilesObserved = $metallicityObserved16->append($metallicityErrorObserved50)->append($metallicityErrorObserved84);

# Convert mass to linear scalings.
if      ( $massScaling eq "linear" ) {
    # No conversion needed.
} elsif ( $massScaling eq "log10"  ) {
    # Convert from log10.
    $massStellarObserved .= 10.0**$massStellarObserved;
} else {
    die('stellarMassMetallicityRelation_Gallazzi_2005.pl: unknown scaling');
}
# Convert metallicity to linear scaling.
if      ( $metallicityScaling eq "linear" ) {
    # No conversion needed.
} elsif ( $metallicityScaling eq "log10"  ) {
    # Convert from log10.
    $metallicityObserved16         .= 10.0**$metallicityObserved16;
    $metallicityObserved50         .= 10.0**$metallicityObserved50;
    $metallicityObserved84         .= 10.0**$metallicityObserved84;
    $percentilesObserved           .= 10.0**$percentilesObserved;
    $metallicityErrorObserved16    .= $metallicityErrorObserved16*log(10.0)*$metallicityObserved16;
    $metallicityErrorObserved50    .= $metallicityErrorObserved50*log(10.0)*$metallicityObserved50;
    $metallicityErrorObserved84    .= $metallicityErrorObserved84*log(10.0)*$metallicityObserved84;
    $metallicityCovarianceObserved .= $metallicityCovarianceObserved*log(10.0)**2*outer($percentilesObserved,$percentilesObserved);
} else {
    die('stellarMassMetallicityRelation_Gallazzi_2005.pl: unknown scaling');
}
# Convert to preferred unit system.
$massStellarObserved           *=  $massUnits/$massSolar    ;
$metallicityCovarianceObserved *= ($massUnits/$massSolar)**2;
foreach ( $percentilesObserved, $metallicityObserved16, $metallicityObserved50, $metallicityObserved84, $metallicityErrorObserved16, $metallicityErrorObserved50, $metallicityErrorObserved84 ) {
    $_ *= $metallicityUnits/$metallicitySolar;
}

# Read model metallicity distribution.
my $model                        = new PDL::IO::HDF5($galacticusFileName);
my $analysisGroup                = $model           ->group  ('analysis'                         )       ;
my $metallicityGroup             = $analysisGroup   ->group  ('sdssStellarMetallicityZ0.07'      )       ;
my $massStellarModel             = $metallicityGroup->dataset('mass'                             )->get();
my $metallicityModel             = $metallicityGroup->dataset('metallicity'                      )->get();
my $metallicityDistributionModel = $metallicityGroup->dataset('metallicityDistribution'          )->get();
my $metallicityCovarianceModel   = $metallicityGroup->dataset('metallicityDistributionCovariance')->get();

# Compute the model metallicity distribution percentiles.
my $metallicityModel16 = pdl zeroes($massCount);
my $metallicityModel50 = pdl zeroes($massCount);
my $metallicityModel84 = pdl zeroes($massCount);
my $quantiles = pdl [ 0.16, 0.50, 0.84 ];
for(my $i=0;$i<$massCount;++$i) {
    my $distributionCumulative = 
	+$metallicityDistributionModel->(:,($i))->cumusumover()
	/$metallicityDistributionModel->(:,($i))->sum        ();
    (my $metallicities) = interpolate($quantiles,$distributionCumulative,$metallicityModel);
    $metallicityModel16->(($i)) .= $metallicities->((0));
    $metallicityModel50->(($i)) .= $metallicities->((1));
    $metallicityModel84->(($i)) .= $metallicities->((2));
}
my $percentilesModel = $metallicityModel16->append($metallicityModel50)->append($metallicityModel84);

# Compute the model metallicity percentiles covariance.
my $realizationCount = 1000;
my $realizations     = pdl zeroes($realizationCount,3*$massCount);
my $covariances      = pdl zeroes(3*$massCount     ,3*$massCount);
# Normalize the model metallicity distribution.
my $metallicityCount = nelem($metallicityModel);
for(my $i=0;$i<$massCount;++$i) {
    my $normalization = $metallicityDistributionModel->(:,($i))->sum();
    $metallicityDistributionModel->(:,$i                                               ) /= $normalization;
    $metallicityCovarianceModel  ->(:,$i*$metallicityCount:($i+1)*$metallicityCount-1  ) /= $normalization;
    $metallicityCovarianceModel  ->(  $i*$metallicityCount:($i+1)*$metallicityCount-1,:) /= $normalization;
}
# Get Cholesky decomposition. We do this only for the non-zero entries in the metallicity
# distribution (since the covrariance is zero for any zero entries). We reconstruct the
# covariance matrix from its eigenvectors to ensure that it is semi-positive definite.
my $nonEmpty                        = which($metallicityDistributionModel > 0.0);
my $nonEmptyCovariance              = $metallicityCovarianceModel->($nonEmpty,$nonEmpty);
(my $eigenVectors, my $eigenValues) = eigens_sym($nonEmptyCovariance);
my $nonEmptyCovarianceReconstructed = $eigenVectors x stretcher($eigenValues) x minv($eigenVectors);
my $cholesky                        = mchol($nonEmptyCovarianceReconstructed);
# Make realizations.
for(my $j=0;$j<$realizationCount;++$j) {
    # Generate random deviates.
    my $deviates = pdl grandom(nelem($nonEmpty));
    # Generate realization of metallicity distribution.
    my $realization = pdl zeroes($metallicityDistributionModel);
    $realization->flat()->($nonEmpty) .= $metallicityDistributionModel->flat()->($nonEmpty)+flat($deviates x $cholesky);
    $realization->($realization < 0.0;?) .= 0.0;
    # Compute percentiles.
    for(my $i=0;$i<$massCount;++$i) {
	my $distributionCumulative = 
	    +$realization->(:,($i))->cumusumover()
	    /$realization->(:,($i))->sum        ();
	(my $metallicities) = interpolate($quantiles,$distributionCumulative,$metallicityModel);
	$realizations->(($j),$i:$i+2*$massCount:$massCount) .= $metallicities;
    }
}
# Compute covariances.
my $means           = $realizations->average();
my $covarianceModel = pdl zeroes(3*$massCount,3*$massCount);
for(my $m=0;$m<3*$massCount;++$m) {
    for(my $n=0;$n<3*$massCount;++$n) {
	$covarianceModel->(($m),($n)) .=
	    sum(
		($realizations->(:,($m))-$means->(($m)))
		*
		($realizations->(:,($n))-$means->(($n)))
	    )
	    /($realizationCount-1);
    }
}
# Extract errors.
my $metallicityErrorModel16 = $covarianceModel->diagonal(0,1)->(0*$massCount:1*$massCount-1)->sqrt();
my $metallicityErrorModel50 = $covarianceModel->diagonal(0,1)->(1*$massCount:2*$massCount-1)->sqrt();
my $metallicityErrorModel84 = $covarianceModel->diagonal(0,1)->(2*$massCount:3*$massCount-1)->sqrt();

# Compute likelihood if required.
if ( exists($arguments{'outputFile'}) ) {
    # Construct the full covariance matrix, which is the covariance matrix of the observations
    # plus that of the model.
    my $fullCovariance                   = $metallicityCovarianceObserved+$covarianceModel;
    # Compute the likelihood.
    my $constraint;
    my $logDeterminant;
    my $offsets;
    my $inverseCovariance;    
    my $logLikelihood = &Covariances::ComputeLikelihood($percentilesObserved,$percentilesModel,$fullCovariance, determinant => \$logDeterminant, inverseCovariance => \$inverseCovariance, offsets => \$offsets, quiet => $arguments{'quiet'}, inversionMethod => "eigendecomposition");
    $constraint->{'logLikelihood'} = $logLikelihood;
    # Find the Jacobian of the log-likelihood with respect to the model mass function.
    my $jacobian = pdl zeroes(1,nelem($percentilesModel));
    for(my $i=0;$i<nelem($percentilesModel);++$i) {
	$jacobian->((0),($i)) .= sum($inverseCovariance->(($i),:)*$offsets);
    }
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
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [5.0e8:1.0e12]\n";
    print $gnuPlot "set yrange [0.05:3.0]\n";
    print $gnuPlot "set title offset 0,-0.5 'Stellar-phase mass-metallicity relation'\n";
    print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star\\ [{\\rm M}_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Metallicity; \$Z\\ [{\\rm Z}_\\odot]\$'\n";
    print $gnuPlot "set pointsize 1.0\n";
    &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massStellarObserved,
     				  $metallicityObserved50,
				  errorUp   => $metallicityErrorObserved50,
				  errorDown => $metallicityErrorObserved50,
				  style     => "point",
				  symbol    => [6,7],
     				  weight    => [3,1],
     				  color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
     				  title     => "Gallazzi et al. (2005; 50\\\\%)"
     	);
     &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massStellarObserved,
     				  $metallicityObserved16,
 				  errorUp   => $metallicityErrorObserved16,
				  errorDown => $metallicityErrorObserved16,
    				  style     => "point",
 				  symbol    => [6,7],
     				  weight    => [3,1],
    				  color     => $PrettyPlots::colorPairs{'lightSkyBlue'},
     				  title     => "Gallazzi et al. (2005; 16/84\\\\%)"
     	);
     &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massStellarObserved,
     				  $metallicityObserved84,
				  errorUp   => $metallicityErrorObserved84,
				  errorDown => $metallicityErrorObserved84,
     				  style     => "point",
				  symbol    => [6,7],
     				  weight    => [3,1],
     				  color     => $PrettyPlots::colorPairs{'lightSkyBlue'}
	 );

    &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massStellarModel,
     				  $metallicityModel50,
 				  errorUp   => $metallicityErrorModel50,
				  errorDown => $metallicityErrorModel50,
  				  style     => "point",
				  symbol    => [6,7],
     				  weight    => [3,1],
     				  color     => $PrettyPlots::colorPairs{'mediumSeaGreen'},
     				  title     => "Galacticus (50\\\\%)"
     	);
     &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massStellarModel,
     				  $metallicityModel16,
  				  errorUp   => $metallicityErrorModel16,
				  errorDown => $metallicityErrorModel16,
   				  style     => "point",
 				  symbol    => [6,7],
     				  weight    => [3,1],
    				  color     => $PrettyPlots::colorPairs{'lightSeaGreen'},
     				  title     => "Galacticus (16/84\\\\%)"
     	);
     &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massStellarModel,
     				  $metallicityModel84,
  				  errorUp   => $metallicityErrorModel84,
				  errorDown => $metallicityErrorModel84,
      				  style     => "point",
				  symbol    => [6,7],
     				  weight    => [3,1],
     				  color     => $PrettyPlots::colorPairs{'lightSeaGreen'}
	 );
    
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
}

exit;
