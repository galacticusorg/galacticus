# Contains a Perl module which implements various useful functionality for constraints based on the distribution of galaxy sizes.

package Galacticus::Constraints::GalaxySizes;
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
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

sub SDSS2003 {
    # Perform likelihood analysis for the Shen et al. (2003) SDSS galaxy size distribution.
    my $galacticusFileName =   shift() ;
    my $distributionNumber =   shift() ;
    my %options            = %{shift()};
    # Distribution number label.
    my $distributionLabel = sprintf("%2.2i",$distributionNumber);
    # Validate distribution number.
    die("Galacticus::Constraints::GalaxySizes::SDSS2003(): distributionNumber ∈ [1..34] is required")
	if ( $distributionNumber < 1 || $distributionNumber > 34 );
    # Read data.
    my $dataCompilation = new PDL::IO::HDF5("data/observations/galaxySizes/Galaxy_Sizes_By_Mass_SDSS_Shen_2003.hdf5");
    my $dataSet         = $dataCompilation->group('distribution'.$distributionLabel);
    my $data;
    ($data->{$_}) = $dataSet->attrGet($_)
	foreach ( 'sersicIndexMaximum', 'massMinimum', 'massMaximum' );
    $data->{$_} = $dataSet->dataset($_)->get()
	foreach ( 'radius', 'radiusFunction', 'radiusFunctionError' );
    $data->{'covariance'} = stretcher($data->{'radiusFunctionError'}**2);
    # Read model.
    my $model        = new PDL::IO::HDF5($galacticusFileName);
    my $analyses     = $model   ->group('analyses'                          );
    my $distribution = $analyses->group('galaxySizesSDSS'.$distributionLabel);
    $model->{'radiusFunction'     } = $distribution->dataset('galaxySizesSDSSFunction'          )->get();
    $model->{'covariance'         } = $distribution->dataset('galaxySizesSDSSFunctionCovariance')->get();
    $model->{'radiusFunctionError'} = $model->{'covariance'}->diagonal(0,1);
    # Evaluate the model likelihood.
    if ( exists($options{'outputFile'}) ) {
	# Construct the full covariance matrix, which is the covariance matrix of the observations
	# plus that of the model.
	my $fullCovariance                   = $data->{'covariance'}+$model->{'covariance'};
	# Identify upper limits.
	my $upperLimits                      = which($data->{'radiusFunction'} < 0.0);
	my $dataRadiusFunction               =       $data->{'radiusFunction'}->copy();
	$dataRadiusFunction->($upperLimits) .= -$dataRadiusFunction->($upperLimits);
	# Where model points are zero, set them equal to a small fraction of the data. Where the data is an upper limit, this will not
	# affect the answers. Where data has a value this will ensure a low likelihood.
	my $modelRadiusFunction              = $model->{'radiusFunction'}->copy();
	my $modelZero                        = which($modelRadiusFunction <= 0.0);
	$modelRadiusFunction->($modelZero)  .= 1.0e-3*$dataRadiusFunction->($modelZero);
	# Find the peak value, and exclude from the likelihood calculation - this accounts for the fact that these distributions
	# are normalized to unity which introduces correlation between points (i.e. given N-1 points, remaining point is
	# determined by unitarity).
	my $include = pdl ones($data->{'radiusFunction'});
	# Get the index of the peak value in this size function.
	my $maxIndex = maximum_ind($data->{'radiusFunction'});
	# Mark this element as excluded.
	$include->(($maxIndex)) .= 0;
	# Create excluded copies.
	my $includeIndices              = which($include == 1.0);
	my $dataRadiusFunctionExcluded  = $dataRadiusFunction                 ->           ($includeIndices                )  ;
	my $modelRadiusFunctionExcluded = $modelRadiusFunction                ->           ($includeIndices                )  ;
	my $modelCovarianceExcluded     = $model              ->{'covariance'}->           ($includeIndices,$includeIndices)  ;
	my $fullCovarianceExcluded      = $fullCovariance                     ->           ($includeIndices,$includeIndices)  ;
	my $shiftedIndices              = $include                            ->cumusumover(                               )-1;
	my $upperLimitsExcluded         = $shiftedIndices                     ->           ($upperLimits                   )  ;
	# Apply discrepancies.
	if ( exists($options{'modelDiscrepancies'}) ) {
	    my $modelCovarianceExcludedOriginal = $modelCovarianceExcluded               ->copy();
	    my $modelErrorExcluded              = $modelCovarianceExcluded->diagonal(0,1)->copy();
	    # Limit multiplicative covariances which can be too large due to noise.
	    &Galacticus::Constraints::DiscrepancyModels::Apply_Discrepancies(
		"discrepancyGalaxySizeZ0.07_".$distributionLabel.".hdf5",
		$options                    {'modelDiscrepancies'}      ,
		$modelRadiusFunctionExcluded                            ,
		$modelErrorExcluded                                     ,
		$modelCovarianceExcluded                                ,
		limitMultiplicativeCovariance => 1.0
		);
	    # Find the change in the model covariance matrix and add to the full covariance matrix.
	    my $modelCovarianceExcludedChange  = $modelCovarianceExcluded      -$modelCovarianceExcludedOriginal;
	    $fullCovarianceExcluded           += $modelCovarianceExcludedChange                                 ;
	}
	# Compute the likelihood.
	my $constraint;
	my $logDeterminant;
	my $offsets;
	my $inverseCovariance;
	my $jacobian;
	my $logLikelihood =
	    &Galacticus::Constraints::Covariances::ComputeLikelihood
	    (
	     $dataRadiusFunctionExcluded                          ,
	     $modelRadiusFunctionExcluded                         ,
	     $fullCovarianceExcluded                              ,
	     upperLimits                  =>  $upperLimitsExcluded,
	     jacobian                     => \$jacobian           ,
	     offsets                      => \$offsets            ,
	     quiet                        =>  $options{'quiet'} ,
	     inversionMethod              => "eigendecomposition" ,
	     productMethod                => "linearSolver"
	    );
	$constraint->{'label'        } = "galaxySizeZ0.07";
	$constraint->{'logLikelihood'} = $logLikelihood;
	# Compute the variance in the log-likelihood due to errors in the model.
	my $logLikelihoodVariance = transpose($jacobian) x $modelCovarianceExcluded x $jacobian;
	$constraint->{'logLikelihoodVariance'} = $logLikelihoodVariance->sclr();
	# Output the constraint.
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
	open(oHndl,">".$options{'outputFile'});
	print oHndl $xmlOutput->XMLout($constraint);
	close(oHndl);
    }

    # Output the results to file if requested.
    if ( exists($options{'resultFile'}) ) {
	# Identify upper limits.
	my $upperLimits                      = which($data->{'radiusFunction'} < 0.0);
	my $dataRadiusFunction               =       $data->{'radiusFunction'}->flat()->copy();
	$dataRadiusFunction->($upperLimits) .= -$dataRadiusFunction->($upperLimits);
	# For each function, find the peak value, and exclude from the points.
	my $exclude = pdl ones($data->{'radiusFunction'});
	for(my $i=0;$i<$data->{'radiusFunction'}->dim(1);++$i) {
	    # Get the index of the peak value in this size function.
	    my $maxIndex = maximum_ind($data->{'radiusFunction'}->(:,($i)));
	    # Mark this element as excluded.
	    $exclude->(($maxIndex),($i)) .= 0;
	}
	# Create excluded copies.
	my $includeIndices              = which($exclude == 1.0);
	my $radiusExcluded              = $data              ->{'radius'        }->flat()->($includeIndices                );
	my $dataRadiusFunctionExcluded  = $dataRadiusFunction                            ->($includeIndices                );
	my $modelRadiusFunctionExcluded = $model             ->{'radiusFunction'}->flat()->($includeIndices                );
	my $dataCovarianceExcluded      = $data              ->{'covariance'    }        ->($includeIndices,$includeIndices);
	my $modelCovarianceExcluded     = $model             ->{'covariance'    }        ->($includeIndices,$includeIndices);
	my $errorExcluded               = sqrt($dataCovarianceExcluded->diagonal(0,1)+$modelCovarianceExcluded->diagonal(0,1));
	# Write out the results.
	my $resultsFile = new PDL::IO::HDF5(">".$options{'resultFile'});
	$resultsFile->dataset('x'             )->set($radiusExcluded             );
	$resultsFile->dataset('y'             )->set($modelRadiusFunctionExcluded);
	$resultsFile->dataset('error'         )->set($errorExcluded              );
	$resultsFile->dataset('covariance'    )->set($modelCovarianceExcluded    );
	$resultsFile->dataset('yData'         )->set($dataRadiusFunctionExcluded );
	$resultsFile->dataset('covarianceData')->set($dataCovarianceExcluded     );
    }

    # Create a plot of the radius function.
    if ( exists($options{'plotFile'}) ) {
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileTeX, $plot);
	# Open a pipe to GnuPlot.
	($plotFileTeX = $options{'plotFile'}) =~ s/\.pdf$/.tex/;
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
	print $gnuPlot "set output '".$plotFileTeX."'\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.2,0.8\n";
	print $gnuPlot "set key left\n";
	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set logscale x\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	my $xMinimum = 0.8*minimum($data->{'radius'});
	my $xMaximum = 1.2*maximum($data->{'radius'});
	my $yMinimum = 0.0;
	my $yMaximum = 1.2*maximum($data->{'radiusFunction'}->append($model->{'radiusFunction'}));
	print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
	print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
	my $plogMassMinimum = sprintf("%5.2f",log10($data->{'massMinimum'}));
	my $plogMassMaximum = sprintf("%5.2f",log10($data->{'massMaximum'}));
	print $gnuPlot "set title offset -3,-0.9 'Half-light radius distribution; ".($data->{'sersicIndexMaximum'} <= 2.5 ? "late" : "early")."-type; \$".$plogMassMinimum." < \\log_{10}(M_\\star/M_\\odot) < ".$plogMassMaximum." \$'\n";
	print $gnuPlot "set xlabel 'Petrosian half-light radius; \$r_{50}\$ [kpc]'\n";
	print $gnuPlot "set ylabel 'Fraction; \${\\rm d}F/{\\rm d}\\log_{10}r_{50}\$ [dex\$^{-1}\$]'\n";
	# Identify upper limits and construct appropriate inputs for plotting functions.
	my $upperLimits                             = which($data->{'radiusFunction'} < 0.0);
	my $dataRadiusFunction                      =       $data->{'radiusFunction'}->copy();
	$dataRadiusFunction->($upperLimits)        .= -$dataRadiusFunction->($upperLimits);
	my $errorUp                                 = $data->{'radiusFunctionError'}->copy();
	my $errorDown                               = $data->{'radiusFunctionError'}->copy();
	$errorUp           ->($upperLimits)        .= +0.0;
	$errorDown         ->($upperLimits)        .= -0.5;
	if ( nelem($upperLimits) > 0 ) {
	    my $shortArrows                             = which($dataRadiusFunction->($upperLimits) < 0.6);
	    $errorDown->($upperLimits)->($shortArrows) .= -$dataRadiusFunction->($upperLimits)->($shortArrows)+0.1;
	}
	&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
					       $data->{'radius'},$dataRadiusFunction,
					       errorUp   => $errorUp,
					       errorDown => $errorDown,
					       style     => "point",
					       symbol    => [6,7], 
					       weight    => [5,3],
					       color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
					       title     => "Shen et al. (2003)"
	    );   
	&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
					       $data->{'radius'},$model->{'radiusFunction'},
					       errorUp   => $model->{'radiusFunctionError'},
					       errorDown => $model->{'radiusFunctionError'},
					       style     => "point",
					       symbol    => [6,7], 
					       weight    => [5,3],
					       color     => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
					       title     => "Galacticus"
	    );
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX,margin => 1);
    }
}

1;
