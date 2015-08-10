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

# Compute likelihood (and make a plot) for a Galacticus model given the black hole mass distribution compilation at z~0 from
# Kormendy & Ho (2013).

# Get name of input and output files.
die("blackHoleMassDistribution_KormendyHo2013.pl <galacticusFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Define Galacticus unit system.
my $massSolar = pdl 1.989e30;

# Read observational data.
my  $observations                = new PDL::IO::HDF5($galacticusPath."data/observations/blackHoles/blackHoleMassVsBulgeMass_KormendyHo2013.hdf5");
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
    # No conversion needed.
} elsif ( $massBlackHoleScaling eq "log10"  ) {
    # Convert from log10.
    $massBlackHoleObserved      .= 10.0**$massBlackHoleObserved;
    $massBlackHoleErrorObserved .= $massBlackHoleErrorObserved*log(10.0)*$massBlackHoleObserved;
} else {
    die('blackHoleMassDistribution_KormendyHo2013.pl: unknown scaling');
}
# Convert to preferred unit system.
$massBulgeObserved          *=  $massBulgeUnits    /$massSolar;
$massBulgeErrorObserved     *=  $massBulgeUnits    /$massSolar;
$massBlackHoleObserved      *=  $massBlackHoleUnits/$massSolar;
$massBlackHoleErrorObserved *=  $massBlackHoleUnits/$massSolar;

# Read model black hole mass distribution.
my $model                           = new PDL::IO::HDF5($galacticusFileName);
my $analysisGroup                   = $model           ->group  ('analysis'                           )       ;
my $metallicityGroup                = $analysisGroup   ->group  ('blackHoleDistributionKormendyHo'    )       ;
my $massBulgeModel                  = $metallicityGroup->dataset('mass'                               )->get();
my $massBlackHoleModel              = $metallicityGroup->dataset('blackHoleMass'                      )->get();
my $massDistributionModel           = $metallicityGroup->dataset('blackHoleMassDistribution'          )->get();
my $massDistributionCovarianceModel = $metallicityGroup->dataset('blackHoleMassDistributionCovariance')->get();

# Normalize the model distribution.
my $normalization = $massDistributionModel->sumover();
for(my $i=0;$i<nelem($massBulgeModel);++$i) {
    $massDistributionModel->(:,($i)) /= $normalization->(($i))
	if ( $normalization->(($i)) > 0.0 );
}
my $distributionMaximum = max($massDistributionModel);
my $distributionMinimum = 10.0**floor(log10(1.0e-3*$distributionMaximum));

# Compute likelihood if required.
if ( exists($arguments{'outputFile'}) ) {
    my $realizationCount = 1000;
    # Get Cholesky decomposition. We do this only for the non-zero entries in the metallicity
    # distribution (since the covrariance is zero for any zero entries). We reconstruct the
    # covariance matrix from its eigenvectors to ensure that it is semi-positive definite.
    my $nonEmpty                        = which($massDistributionModel > 1.0e-3);
    my $nonEmptyCovariance              = $massDistributionCovarianceModel->($nonEmpty,$nonEmpty);
    (my $eigenVectors, my $eigenValues) = eigens_sym($nonEmptyCovariance);    
    my $nonEmptyCovarianceReconstructed = $eigenVectors x stretcher($eigenValues) x minv($eigenVectors);
    my $cholesky                        = mchol($nonEmptyCovarianceReconstructed);
    # Make realizations.
    my $likelihood        = pdl 0.0;
    my $likelihoodSquared = pdl 0.0;
    my $logLikelihoodBase;
    for(my $j=0;$j<$realizationCount;++$j) {
	# Generate random deviates.
	my $deviates                       = pdl grandom(nelem($nonEmpty));
	# Generate realization of metallicity distribution.
	my $realization                    = pdl zeroes($massDistributionModel);
	$realization->flat()->($nonEmpty) .= $massDistributionModel->flat()->($nonEmpty)+flat($deviates x $cholesky);
	$realization->($realization < 0.0;?) .= 0.0;
	# Iterate over observed galaxies.
	my $logLikelihoodRealization = pdl 0.0;
	for(my $i=0;$i<nelem($massBulgeObserved);++$i) {
	    my $massBulge           = $massBulgeObserved         ->(($i));
	    my $massBlackHole       = $massBlackHoleObserved     ->(($i));
	    my $massBulgeError      = $massBulgeErrorObserved    ->(($i));
	    my $massBlackHoleError  = $massBlackHoleErrorObserved->(($i));
	    my $offsetMassBulge     = ($massBulgeModel    -$massBulge    )/$massBulgeError    ;
	    my $offsetMassBlackHole = ($massBlackHoleModel-$massBlackHole)/$massBlackHoleError;
	    my $probabilityObserved =
		outer
		(
		 exp(-0.5*$offsetMassBlackHole**2),
		 exp(-0.5*$offsetMassBulge    **2)
		);
	    $probabilityObserved      /= $probabilityObserved->sum();
	    my $probabilityModel       = sum($probabilityObserved*$realization);
	    if ( $probabilityModel > 0.0 ) {
		$logLikelihoodRealization += log($probabilityModel);
	    } else {
		$logLikelihoodRealization -= 1.0e30;
	    }
	}
	if ( $j == 0 ) {
	    $logLikelihoodBase  =          $logLikelihoodRealization                     ;
	} else {
	    $likelihood        += exp(     $logLikelihoodRealization-$logLikelihoodBase );
	    $likelihoodSquared += exp(2.0*($logLikelihoodRealization-$logLikelihoodBase));
	}
    }
    $likelihood            /= $realizationCount;
    my $likelihoodVariance  = ($likelihoodSquared-$realizationCount*$likelihood**2)/($realizationCount-1);
    my $logLikelihood       = log($likelihood)+$logLikelihoodBase;
    my $logLikelihoodVariance = $likelihoodVariance/$likelihood**2;
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
    print $gnuPlot "set key at screen 0.20,0.9\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [3.0e8:3.0e12]\n";
    print $gnuPlot "set yrange [3.0e5:1.0e11]\n";
    print $gnuPlot "set title offset 0,-0.5 'Black hole mass distribution'\n";
    print $gnuPlot "set xlabel 'Bulge stellar mass; \$M_{\\star, \\rm bulge}\\ [{\\rm M}_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'Black hole mass; \$M_\\bullet\\ [{\\rm M}_\\odot]\$'\n";
    print $gnuPlot "set pointsize 1.0\n";
    print $gnuPlot "set pm3d map\n";
    print $gnuPlot "set pm3d explicit\n";
    print $gnuPlot "set logscale zcb\n";
    print $gnuPlot "set cbrange [".$distributionMinimum.":".$distributionMaximum."]\n";
    print $gnuPlot "set palette rgbformulae 34,35,36\n";
    print $gnuPlot "set palette negative\n";
    print $gnuPlot "set multiplot\n";    
    # Plot color-shading.
    print $gnuPlot "splot '-' with pm3d notitle\n";
    for(my $i=0;$i<nelem($massBulgeModel);++$i) {
	for(my $j=0;$j<nelem($massBlackHoleModel);++$j) {
	    print $gnuPlot $massBulgeModel->index($i)." ".$massBlackHoleModel->index($j)." ".$massDistributionModel->(($j),($i))."\n";
	}
	print $gnuPlot "\n" unless ( $i == nelem($massBulgeModel)-1 );
    }
    print $gnuPlot "e\n";
    print $gnuPlot "unset pm3d\n";
    print $gnuPlot "unset label; unset border; unset xtics; unset ytics; unset x2tics; unset y2tics; set xlabel ''; set ylabel ''\n";
    # Plot points.
    &PrettyPlots::Prepare_Dataset(\$plot,
     				  $massBulgeObserved,
     				  $massBlackHoleObserved,
   				  errorUp  => $massBlackHoleErrorObserved,
   				  errorDown => $massBlackHoleErrorObserved,
   				  errorLeft   => $massBulgeErrorObserved,
   				  errorRight => $massBulgeErrorObserved,
   				  style     => "point",
   				  symbol    => [6,7],
     				  weight    => [3,1],
     				  color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
     				  title     => "Kormendy \\\\& Ho (2013)"
     	);
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot, multiPlot => 1);    
    print $gnuPlot "unset multiplot\n";
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
}

exit;
