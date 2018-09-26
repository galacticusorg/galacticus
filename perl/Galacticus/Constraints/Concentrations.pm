# Contains a Perl module which implements various useful functionality for constraints based on the distribution of halo
# concentrations.

package Galacticus::Constraints::Concentrations;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::Math;
use PDL::IO::HDF5;
use PDL::LinearAlgebra;
use Galacticus::Options;
use Galacticus::HDF5;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use LaTeX::Format;

sub COCOCDM {
    # Perform likelihood analysis for the distribution of concentrations from the COCO CDM simulations.
    my $galacticusFileName =   shift() ;
    my $distributionNumber =   shift() ;
    my %options            = %{shift()};
    # Distribution number label.
    my $distributionLabel = sprintf("%2.2i",$distributionNumber);
    # Validate distribution number.
    die("Galacticus::Constraints::Concentrations::COCOCDM(): distributionNumber ∈ [1..7] is required")
	if ( $distributionNumber < 1 || $distributionNumber > 7 );
    # Read data.
    my $dataCompilation = new PDL::IO::HDF5($ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/concentrationDistributionCocoCDM.hdf5");
    my $mass                   = $dataCompilation->dataset('mass'                  )->get();
    my $concentration          = $dataCompilation->dataset('concentration'         )->get();
    my $distribution           = $dataCompilation->dataset('distribution'          )->get();
    my $distributionCovariance = $dataCompilation->dataset('distributionCovariance')->get();
    my $data;
    $data->{'concentration'    } = $concentration                                                                                         ;
    $data->{'distribution'     } = $distribution                           ->(  :,($distributionNumber-1))                                ;
    $data->{'covariance'       } = $distributionCovariance                 ->(:,:,($distributionNumber-1))                                ;
    $data->{'massMinimum'      } = $mass                                   ->(    ($distributionNumber-1))/sqrt($mass->((1))/$mass->((0)));
    $data->{'massMaximum'      } = $mass                                   ->(    ($distributionNumber-1))*sqrt($mass->((1))/$mass->((0)));
    $data->{'distributionError'} = $data                  ->{'covariance'}->diagonal(0,1)->sqrt();
    # Read model.
    my $model        = new PDL::IO::HDF5($galacticusFileName);
    my $analyses     = $model   ->group('analyses'                                           );
    my $modelData    = $analyses->group('concentrationDistributionCDMCOCO'.$distributionLabel);
    $model->{'distribution'     } = $modelData->dataset('concentrationFunction'          )->get();
    $model->{'covariance'       } = $modelData->dataset('concentrationFunctionCovariance')->get();
    $model->{'distributionError'} = $model->{'covariance'}->diagonal(0,1);
    # Evaluate the model likelihood.
    if ( exists($options{'outputFile'}) ) {
    	# Construct the full covariance matrix, which is the covariance matrix of the observations
    	# plus that of the model.
    	my $fullCovariance                 = $data->{'covariance'}+$model->{'covariance'};
    	# Where model points are zero, set them equal to a small fraction of the data. This will ensure a low likelihood.
    	my $modelDistribution              = $model->{'distribution'}->copy();
    	my $modelZero                      = which($modelDistribution <= 0.0);
    	$modelDistribution->($modelZero)  .= 1.0e-3*$data->{'distribution'}->($modelZero);
    	# Find all non-zero entries in the data.
    	my $include                        = which($data->{'distribution'} > 0.0);
    	# Create excluded copies.
    	my $dataDistributionExcluded  = $data             ->{'distribution'}->($include         );
    	my $modelDistributionExcluded = $modelDistribution                  ->($include         );
    	my $modelCovarianceExcluded   = $model            ->{'covariance'  }->($include,$include);
    	my $fullCovarianceExcluded    = $fullCovariance                     ->($include,$include);
    	# Apply discrepancies.
    	if ( exists($options{'modelDiscrepancies'}) ) {
    	    my $modelCovarianceExcludedOriginal = $modelCovarianceExcluded               ->copy();
    	    my $modelErrorExcluded              = $modelCovarianceExcluded->diagonal(0,1)->copy();
    	    # Limit multiplicative covariances which can be too large due to noise.
    	    &Galacticus::Constraints::DiscrepancyModels::Apply_Discrepancies(
    		"discrepancyConcentrationCOCOCDM_".$distributionLabel.".hdf5",
    		$options                    {'modelDiscrepancies'}           ,
    		$modelDistributionExcluded                                   ,
    		$modelErrorExcluded                                          ,
    		$modelCovarianceExcluded                                     ,
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
    	     productMethod              => "linearSolver"       ,
	     errorTolerant              => 1
    	    );
    	$constraint->{'label'        } = "concentrationCDMCOCO";
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

    # Create a plot of the concentration function.
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
    	print $gnuPlot "set key at screen 0.25,0.8\n";
    	print $gnuPlot "set key left\n";
    	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set logscale x\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	my $nonZero  = which($data->{'distribution'} > 0.0);
    	my $xMinimum = minimum($data->{'concentration'}->($nonZero))/1.1;
    	my $xMaximum = maximum($data->{'concentration'}->($nonZero))*1.1;
	$xMaximum   .= 30.0
	    if ( $xMaximum > 30.0 );
    	my $yMinimum = 0.0;
    	my $yMaximum = 1.1*maximum($data->{'distribution'}->append($model->{'distribution'}));
    	print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
    	print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
    	print $gnuPlot "set title offset -3,-1 'COCO CDM; ".&LaTeX::Format::Number($data->{'massMinimum'},2)." \$< M_\\mathrm{200c}/\\mathrm{M}_\\odot <\$ ".&LaTeX::Format::Number($data->{'massMaximum'},2)."'\n";
    	print $gnuPlot "set xlabel 'Concentration; \$c_\\mathrm{200c}\$'\n";
    	print $gnuPlot "set ylabel 'Distribution; \$\\mathrm{d}p/\\mathrm{d}\\log c_\\mathrm{200c}\$'\n";
    	&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
    					       $data->{'concentration'},$data->{'distribution'},
    					       errorUp   => $data->{'distributionError'},
    					       errorDown => $data->{'distributionError'},
    					       style     => "point",
    					       symbol    => [6,7], 
    					       weight    => [5,3],
    					       color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
    					       title     => "COCO CDM"
    	    );   
    	&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
    					       $data->{'concentration'},$model->{'distribution'},
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
