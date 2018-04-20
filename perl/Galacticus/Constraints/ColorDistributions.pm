# Contains a Perl module which implements various useful functionality for constraints based on the distribution of galaxy colors.

package Galacticus::Constraints::ColorDistributions;
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

sub SDSS2004 {
    # Perform likelihood analysis for the Baldry et al. (2004) SDSS galaxy color distribution.
    my $galacticusFileName =   shift() ;
    my $distributionNumber =   shift() ;
    my %options            = %{shift()};
    # Distribution number label.
    my $distributionLabel = sprintf("%2.2i",$distributionNumber);
    # Validate distribution number.
    die("Galacticus::Constraints::ColorDistributions::SDSS2004(): distributionNumber ∈ [1..16] is required")
	if ( $distributionNumber < 1 || $distributionNumber > 16 );
    # Read data.
    my $dataCompilation = new PDL::IO::HDF5("data/observations/galaxyColors/colorDistributionsBaldry2004.hdf5");
    my $dataSet         = $dataCompilation->group('distribution'.$distributionLabel);
    my $data;
    ($data->{$_}) = $dataSet->attrGet($_)
	foreach ( 'magnitudeMinimum', 'magnitudeMaximum' );
    $data->{$_} = $dataSet->dataset($_)->get()
	foreach ( 'color', 'distribution', 'distributionError' );
    $data->{'covariance'} = stretcher($data->{'distributionError'}**2);
    # Read model.
    my $model        = new PDL::IO::HDF5($galacticusFileName);
    my $analyses     = $model   ->group('analyses'                                );
    my $distribution = $analyses->group('colorDistributionSDSS'.$distributionLabel);
    $model->{'distribution'     } = $distribution->dataset('colorDistributionSDSSFunction'          )->get();
    $model->{'covariance'       } = $distribution->dataset('colorDistributionSDSSFunctionCovariance')->get();
    $model->{'distributionError'} = $model->{'covariance'}->diagonal(0,1);
    # Evaluate the model likelihood.
    if ( exists($options{'outputFile'}) ) {
	# Construct the full covariance matrix, which is the covariance matrix of the observations
	# plus that of the model.
	my $fullCovariance                   = $data->{'covariance'}+$model->{'covariance'};
	# Where model points are zero, set them equal to a small fraction of the data. This will ensure a low likelihood.
	my $modelDistribution              = $model->{'distribution'}->copy();
	my $modelZero                      = which($modelDistribution <= 0.0);
	$modelDistribution->($modelZero)  .= 1.0e-3*$data->{'distribution'}->($modelZero);
	# Find the peak value, and exclude from the likelihood calculation - this accounts for the fact that these distributions
	# are normalized to unity which introduces correlation between points (i.e. given N-1 points, remaining point is
	# determined by unitarity).
	my $include  = pdl ones($data->{'distribution'});
	# Get the index of the peak value in this size function.
	my $maxIndex = maximum_ind($data->{'distribution'});
	# Mark this element as excluded.
	$include->(($maxIndex)) .= 0;
	# Create excluded copies.
	my $includeIndices            = which($include == 1.0);
	my $dataDistributionExcluded  = $data             ->{'distribution'}->($includeIndices                );
	my $modelDistributionExcluded = $modelDistribution                  ->($includeIndices                );
	my $modelCovarianceExcluded   = $model            ->{'covariance'  }->($includeIndices,$includeIndices);
	my $fullCovarianceExcluded    = $fullCovariance                     ->($includeIndices,$includeIndices);
	# Apply discrepancies.
	if ( exists($options{'modelDiscrepancies'}) ) {
	    my $modelCovarianceExcludedOriginal = $modelCovarianceExcluded               ->copy();
	    my $modelErrorExcluded              = $modelCovarianceExcluded->diagonal(0,1)->copy();
	    # Limit multiplicative covariances which can be too large due to noise.
	    &Galacticus::Constraints::DiscrepancyModels::Apply_Discrepancies(
		"discrepancyGalaxyColorZ0.07_".$distributionLabel.".hdf5",
		$options                    {'modelDiscrepancies'}       ,
		$modelDistributionExcluded                               ,
		$modelErrorExcluded                                      ,
		$modelCovarianceExcluded                                 ,
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
	     $dataDistributionExcluded                          ,
	     $modelDistributionExcluded                         ,
	     $fullCovarianceExcluded                            ,
	     jacobian                   => \$jacobian           ,
	     offsets                    => \$offsets            ,
	     quiet                      =>  $options{'quiet'}   ,
	     inversionMethod            => "eigendecomposition" ,
	     productMethod              => "linearSolver"
	    );
	$constraint->{'label'        } = "galaxyColorSDSSZ0.10";
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
	my $dataDistribution = $data->{'distribution'}->flat()->copy();
	# For each function, find the peak value, and exclude from the points.
	my $exclude = pdl ones($data->{'distribution'});
	for(my $i=0;$i<$data->{'distribution'}->dim(1);++$i) {
	    # Get the index of the peak value in this color distribution.
	    my $maxIndex = maximum_ind($data->{'distribution'}->(:,($i)));
	    # Mark this element as excluded.
	    $exclude->(($maxIndex),($i)) .= 0;
	}
	# Create excluded copies.
	my $includeIndices            = which($exclude == 1.0);
	my $colorExcluded             = $data              ->{'color'       }->flat()->($includeIndices                );
	my $dataDistributionExcluded  = $dataDistribution                            ->($includeIndices                );
	my $modelDistributionExcluded = $model             ->{'distribution'}->flat()->($includeIndices                );
	my $dataCovarianceExcluded    = $data              ->{'covariance'  }        ->($includeIndices,$includeIndices);
	my $modelCovarianceExcluded   = $model             ->{'covariance'  }        ->($includeIndices,$includeIndices);
	my $errorExcluded             = sqrt($dataCovarianceExcluded->diagonal(0,1)+$modelCovarianceExcluded->diagonal(0,1));
	# Write out the results.
	my $resultsFile = new PDL::IO::HDF5(">".$options{'resultFile'});
	$resultsFile->dataset('x'             )->set($colorExcluded            );
	$resultsFile->dataset('y'             )->set($modelDistributionExcluded);
	$resultsFile->dataset('error'         )->set($errorExcluded            );
	$resultsFile->dataset('covariance'    )->set($modelCovarianceExcluded  );
	$resultsFile->dataset('yData'         )->set($dataDistributionExcluded );
	$resultsFile->dataset('covarianceData')->set($dataCovarianceExcluded   );
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
	my $xMinimum = minimum($data->{'color'})-0.1;
	my $xMaximum = maximum($data->{'color'})+0.1;
	my $yMinimum = 0.0;
	my $yMaximum = 1.1*maximum($data->{'distribution'}->append($model->{'distribution'}));
	print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
	print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
	my $pMagnitudeMinimum = sprintf("%5.2f",$data->{'magnitudeMinimum'});
	my $pMagnitudeMaximum = sprintf("%5.2f",$data->{'magnitudeMaximum'});
	print $gnuPlot "set title offset -3,-0.9 'SDSS u\$-\$r color distribution; \$".$pMagnitudeMinimum." < M_\\mathrm{r^{0.1}} < ".$pMagnitudeMaximum." \$'\n";
	print $gnuPlot "set xlabel 'Color; \$\\mathrm{u}^{0.1}-\\mathrm{r}^{0.1}\$'\n";
	print $gnuPlot "set ylabel 'Fraction; \$\\mathrm{d}F/\\mathrm{d}(\\mathrm{u}^{0.1}-\\mathrm{r}^{0.1})\$'\n";
	&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
					       $data->{'color'},$data->{'distribution'},
					       errorUp   => $data->{'distributionError'},
					       errorDown => $data->{'distributionError'},
					       style     => "point",
					       symbol    => [6,7], 
					       weight    => [5,3],
					       color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
					       title     => "Baldry et al. (2004)"
	    );   
	&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
					       $data->{'color'},$model->{'distribution'},
					       errorUp   => $model->{'distributionError'},
					       errorDown => $model->{'distributionError'},
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
