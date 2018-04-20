# Contains a Perl module which implements construction of a variety of mass functions from
# Galacticus when fitting to constraints.

package Galacticus::Constraints::LuminosityFunctions;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw(PI);
use PDL::MatrixOps;
use PDL::IO::HDF5;
use Math::SigFigs;
use Data::Dumper;
use LaTeX::Encode;
use XML::Simple;
use Scalar::Util 'reftype';
use Galacticus::HDF5;
use Galacticus::Constraints::Covariances;
use Stats::Histograms;
use Stats::Percentiles;
use List::ExtraUtils;

sub Construct {
    # Construct a luminosity function from Galacticus for constraint purposes.
    my %arguments = %{$_[0]};
    my $config    =   $_[1];
    $arguments{'quiet'} = 0
	unless ( exists($arguments{'quiet'}) );
    
    # Make the label LaTeX compliant.
    $config->{'observationLabel'} = latex_encode($config->{'observationLabel'});
    $config->{'observationLabel'} =~ s/\\/\\\\/g;

    # Check array sizes match.
    foreach ( 'y', 'error' ) {
	die('LuminosityFunctions::Construct(): size of '.$_.' array differs from size of x array')
	    if ( exists($config->{$_}) && nelem($config->{$_}) != nelem($config->{'x'}) );
    }
    
    # Construct a convariance matrix if one is needed.
    unless ( exists($config->{'covariance'}) ) {
	$config->{'covariance'}                 = pdl zeroes(nelem($config->{'error'}),nelem($config->{'error'}));
	$config->{'covariance'}->diagonal(0,1) .= $config->{'error'}**2;
    }

    # Read the model data.
    my $galacticus;
    $galacticus->{'file' } = $config->{'galacticusFile'};
    $galacticus->{'store'} = 0;
    &Galacticus::HDF5::Open_File($galacticus);
    my $yGalacticus          = $galacticus->{'hdf5File'}->group('analyses')->group($config->{'analysisLabel'})->dataset('luminosityFunction'          )->get();
    my $covarianceGalacticus = $galacticus->{'hdf5File'}->group('analyses')->group($config->{'analysisLabel'})->dataset('luminosityFunctionCovariance')->get();
    my $errorGalacticus      = sqrt($covarianceGalacticus->diagonal(0,1));

    # Apply any shifts due to model discrepancy.
    if ( exists($arguments{'modelDiscrepancies'}) ) {
	# Locate the path which contains discrepancies.
	my $discrepancyPath = $arguments{'modelDiscrepancies'};
	# Scan the path for discrepancy files.
	opendir(discrepDir,$discrepancyPath);
	while ( my $discrepancy = readdir(discrepDir) ) {
	    my $discrepancyFileName = $discrepancyPath."/".$discrepancy."/".$config->{'discrepancyFileName'};
	    if ( -e $discrepancyFileName ) {
		my $discrepancyFile = new PDL::IO::HDF5($discrepancyFileName);
		my @datasets = $discrepancyFile->datasets();
		foreach my $dataset ( @datasets ) {
		    if ( $dataset eq "multiplicative" ) {
			# Read the multiplicative discrepancy
			my $multiplier         = $discrepancyFile->dataset('multiplicative')->get();
			$yGalacticus          *= $multiplier;
			$errorGalacticus      *= $multiplier;
			$covarianceGalacticus .= $covarianceGalacticus*outer($multiplier,$multiplier);
		    }		    
		    if ( $dataset eq "multiplicativeCovariance" ) {
			# Adjust the model accordingly.
			my $covarianceMultiplier  = $discrepancyFile->dataset('multiplicativeCovariance')->get();
			$covarianceGalacticus    += $covarianceMultiplier*outer($yGalacticus,$yGalacticus);
			$errorGalacticus         += sqrt($covarianceMultiplier)*$yGalacticus;
		    }		    
		    if ( $dataset eq "additive" ) {
			# Read the additive discrepancy
			my $addition  = $discrepancyFile->dataset('additive')->get();
			# Adjust the model accordingly.
			$yGalacticus += $addition;
		    }
		    if ( $dataset eq "additiveCovariance" ) {
			# Read the covariance of the discrepancy.
			my $covariance         = $discrepancyFile->dataset('additiveCovariance')->get();
			# Adjust the model discrepancy covariance accordingly.
			$covarianceGalacticus += $covariance;
		    }
		}
	    }
	}
    }
    
    # Output the results to file if requested.
    if ( exists($arguments{'resultFile'}) ) {
	my $resultsFile = new PDL::IO::HDF5(">".$arguments{'resultFile'});
	$resultsFile->dataset('x'             )->set($config->{'x'         });
	$resultsFile->dataset('yData'         )->set($config->{'y'         });
	$resultsFile->dataset('covarianceData')->set($config->{'covariance'});
	$resultsFile->dataset('y'             )->set($yGalacticus           );
	$resultsFile->dataset('covariance'    )->set($covarianceGalacticus  );
    }

    # Compute the likelihood:
    if ( exists($arguments{'outputFile'}) ) {
	my $constraint;
	$constraint->{'label'} = $config->{'analysisLabel'};
	# Find zero elements in the observed luminosity function.
	my $nonZeroObserved = which($config->{'y'} > 1.0e-30);
	# Scale the errors on the model luminosity functions to what they would be if we had a perfect match to the data. This prevents
	# high likelihoods for cases where the model vastly exceeds the data (in which case its covariance will also be far too
	# large).
	my $zeroModel                              = which($yGalacticus <= 1.0e-30);
	my $modelCovarianceScale                   = pdl ones($yGalacticus);
	$modelCovarianceScale->($nonZeroObserved) .= sqrt($yGalacticus->($nonZeroObserved)/$config->{'y'}->($nonZeroObserved));
	$modelCovarianceScale->($zeroModel      ) .= 1.0;
	# Construct the full covariance matrix, which is the covariance matrix of the observations
	# plus that of the model (scaled).
	my $fullCovariance        = $config->{'covariance'}+$covarianceGalacticus/outer($modelCovarianceScale,$modelCovarianceScale);
	# If the range of luminosityes over which to constrain is limited, set the model luminosity function outside of that range equal to
	# the observed luminosity function.
	my $yGalacticusLimited = $yGalacticus->copy();
	if ( exists($config->{'constraintLuminosityMinimum'}) ) {
	    my $noConstraint = which($config->{'x'} < $config->{'constraintLuminosityMinimum'});
	    $yGalacticusLimited->($noConstraint) .= $config->{'y'}->($noConstraint)
		if ( nelem($noConstraint) > 0 );
	}
	if ( exists($config->{'constraintLuminosityMaximum'}) ) {
	    my $noConstraint = which($config->{'x'} > $config->{'constraintLuminosityMaximum'});
	    $yGalacticusLimited->($noConstraint) .= $config->{'y'}->($noConstraint)
		if ( nelem($noConstraint) > 0 );
	}
	# If log-normal errors are requested, map y-values and covariances to logarithmic values.
	my $yGalacticusMapped             ;
	my $yDataMapped                   ;
	my $covarianceFullMapped          ;
	my $covarianceGalacticusMapped    ;
	my $isBad                      = 0;
	if ( exists($config->{'errorModel'}) && $config->{'errorModel'} eq "logNormal" ) {
	    if ( any($yGalacticusLimited->($nonZeroObserved) <= 0.0) ) { 
		$isBad = 1;
	    } else {
		# Currently we only map to log-normal errors for model points that are below the corresponding data point. This
		# penalizes models which lie far below the data in such a way that the log-likelihood continues to get worse as
		# the model falls to lower and lower values, while not allowing the model to deviate far above the data with only
		# modest decrease in log-likelihood.
		my $belowData = 
		    which
		    (
		     $yGalacticusLimited       ->($nonZeroObserved) 
		     <
		     $config            ->{'y'}->($nonZeroObserved)
		    );
		$yGalacticusMapped                         =           $yGalacticusLimited       ->($nonZeroObserved                 )->copy() ;
		$yDataMapped                               =           $config            ->{'y'}->($nonZeroObserved                 )->copy() ;
		$covarianceFullMapped                      =           $fullCovariance           ->($nonZeroObserved,$nonZeroObserved)->copy() ;		
		$covarianceGalacticusMapped                =           $covarianceGalacticus     ->($nonZeroObserved,$nonZeroObserved)->copy() ;		
		my $covarianceTransform                    = pdl ones($yGalacticusMapped);
		$covarianceTransform       ->($belowData) .= 1.0/$yDataMapped                    ->($belowData                       )         ;
		$yGalacticusMapped         ->($belowData) .= log($yGalacticusMapped              ->($belowData                       )        );
		$yDataMapped               ->($belowData) .= log($yDataMapped                    ->($belowData                       )        );
		$covarianceFullMapped                     *= outer($covarianceTransform,$covarianceTransform);
		$covarianceGalacticusMapped               *= outer($covarianceTransform,$covarianceTransform);
	    }
	} else {
	    $yGalacticusMapped          =     $yGalacticusLimited         ->($nonZeroObserved                 );
	    $yDataMapped                =     $config              ->{'y'}->($nonZeroObserved                 );
	    $covarianceFullMapped       =     $fullCovariance             ->($nonZeroObserved,$nonZeroObserved);
	    $covarianceGalacticusMapped =     $covarianceGalacticus       ->($nonZeroObserved,$nonZeroObserved);
	}
	# Check for a bad model.
	if ( $isBad ) {
	    # Model was bad, report a very low likelihood.
	    $constraint->{'logLikelihood'        } = -1.0e30;
	    $constraint->{'logLikelihoodVariance'} = +1.0e30;
	} else {
	    # Compute the likelihood.
	    my $offsets;
	    my $jacobian;
	    my $logLikelihood =
		&Galacticus::Constraints::Covariances::ComputeLikelihood
		(
		 $yGalacticusMapped                          ,
		 $yDataMapped                                ,
		 $covarianceFullMapped                       ,
		 jacobian              => \$jacobian         ,
		 offsets               => \$offsets          ,
		 quiet                 => $arguments{'quiet'},
		 productMethod         => "linearSolver" 
		);
	    $constraint->{'logLikelihood'} = $logLikelihood;
	    # Compute the variance in the log-likelihood due to errors in the model.
	    my $logLikelihoodVariance = transpose($jacobian) x $covarianceGalacticusMapped x $jacobian;
	    $constraint->{'logLikelihoodVariance'} = $logLikelihoodVariance->sclr();
	}
	# Output the constraint.
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
	open(oHndl,">".$arguments{'outputFile'});
	print oHndl $xmlOutput->XMLout($constraint);
	close(oHndl);
    }

    # Create a plot of the luminosity function.
    if ( exists($arguments{'plotFile'}) ) {
	require GnuPlot::PrettyPlots;
	require GnuPlot::LaTeX;
	require XMP::MetaData;
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
	my $isMagnitudes = exists($config->{'magnitudes'}) ? $config->{'magnitudes'} : "no";
	if ( $isMagnitudes eq "yes" ) {
	    print $gnuPlot "set logscale y\n";
	    print $gnuPlot "set key at screen 0.45,0.2\n";
	} else {
	    print $gnuPlot "set logscale xy\n";
	    print $gnuPlot "set mxtics 10\n";
	    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	    print $gnuPlot "set key at screen 0.25,0.2\n";
	}
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key left\n";
	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [".$config->{'xRange'}."]\n";
	print $gnuPlot "set yrange [".$config->{'yRange'}."]\n";
	print $gnuPlot "set title offset 0,-1 '".$config->{'title'}."'\n";
	print $gnuPlot "set xlabel '".$config->{'xLabel'}."'\n";
	print $gnuPlot "set ylabel '".$config->{'yLabel'}."'\n";
	# Plot a single model.
	&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
					       $config->{'x'},$config->{'y'},
					       errorUp   => $config->{'error'},
					       errorDown => $config->{'error'},
					       style     => "point",
					       symbol    => [6,7], 
					       weight    => [5,3],
					       pointSize => 0.5,
					       color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
					       title     => $config->{'observationLabel'}
	    );
	&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
					       $config->{'x'},
					       $yGalacticus,
					       errorUp      => $errorGalacticus,
					       errorDown    => $errorGalacticus,
					       style        => "point",
					       symbol       => [6,7], 
					       weight       => [5,3],
					       pointSize    => 0.5,
					       transparency => 0.5,
					       color        => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
					       title        => "Galacticus"
	    );
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX,margin => 1);

    }

}

1;
