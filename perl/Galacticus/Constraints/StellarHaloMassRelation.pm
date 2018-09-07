# Contains a Perl module which implements various useful functionality for constraints based on the stellar mass vs. halo mass
# relation.

package Galacticus::Constraints::StellarHaloMassRelation;
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
use PDL::GSL::INTERP;
use Data::Dumper;

sub COSMOS2012 {
    # Perform likelihood analysis for the Leauthaud et al. (2012) COSMOS stellar mass-halo mass relation.
    my $galacticusFileName = shift();
    my $redshiftRange      = shift();
    my %options;
    (%options) = @_
	if ( scalar(@_) > 0 );

    # Validate redshift range.
    die("Galacticus::Constraints::StellarHaloMassRelation::COSMOS2012(): redshiftRange must be 1, 2, or 3")
	if ( $redshiftRange < 1 || $redshiftRange > 3 );
    
    # Read observational data.
    my $data             = new PDL::IO::HDF5($ENV{'GALACTICUS_DATA_PATH'}."/static/observations/stellarHaloMassRelation/stellarHaloMassRelation_COSMOS_Leauthaud2012.hdf5");
    my $redshiftGroup    = $data         ->group  ('redshiftInterval'.$redshiftRange)       ;
    my $massStellarData  = $redshiftGroup->dataset('massStellar'                    )->get();
    my $massHaloMeanData = $redshiftGroup->dataset('massHaloMean'                   )->get();
    my $massHaloLowData  = $redshiftGroup->dataset('massHaloLow'                    )->get();
    my $massHaloHighData = $redshiftGroup->dataset('massHaloHigh'                   )->get();
    (my $redshiftMinimum, my $redshiftMaximum) = $redshiftGroup->attrGet('redshiftMinimum','redshiftMaximum');
    
    # Find a spline fit to the observed data, and compute the uncertainty in logarithm of halo mass.
    my $massStellarDataLogarithmic   = log($massStellarData );
    my $massHaloMeanDataLogarithmic  = log($massHaloMeanData);
    my $massHaloLowDataLogarithmic   = log($massHaloLowData );
    my $massHaloHighDataLogarithmic  = log($massHaloHighData);
    my $massHaloErrorDataLogarithmic = 0.5*($massHaloHighDataLogarithmic-$massHaloLowDataLogarithmic);
    my $spline                       = PDL::GSL::INTERP->init('cspline',$massHaloMeanDataLogarithmic,$massStellarDataLogarithmic  );
    my $splineError                  = PDL::GSL::INTERP->init('cspline',$massHaloMeanDataLogarithmic,$massHaloErrorDataLogarithmic);

    # Read model data.
    my $model                            = new PDL::IO::HDF5($galacticusFileName);
    my $analysis                         = $model   ->group ('analyses'                                            )
    	                                            ->group ('stellarHaloMassRelationLeauthaud2012z'.$redshiftRange)       ;
    my $massHaloModel                    = $analysis->dataset('massHalo'                                           )->get();
    my $massStellarModel                 = $analysis->dataset('massStellarLog10'                                   )->get();
    my $massStellarCovarianceModel       = $analysis->dataset('massStellarLog10Covariance'                         )->get();
    my $modelParameters                  = $model   ->group('Parameters');
    my $massResolution = pdl 0.0;
    if ( grep {$_ eq "mergerTreeMassResolutionMethod"} $modelParameters->groups() ) {
	my $resolutionGroup = $modelParameters->group('mergerTreeMassResolutionMethod');
	if ( grep {$_ eq "massResolution"} $resolutionGroup->attrs() ) {
	    ($massResolution) = $resolutionGroup->attrGet('massResolution');
	}
    }
    my $modelEntries                     = 
	which(
	    ($massStellarModel                          >  0.0                      ) 
	    &
	    ($massStellarCovarianceModel->diagonal(0,1) >  0.0                      )
	    &
	    ($massHaloModel                             >= $massHaloMeanData->(( 0)))
	    &
	    ($massHaloModel                             <= $massHaloMeanData->((-1)))
	    &
	    ($massHaloModel                             >  20.0*$massResolution     )
	);
    my $massHaloModelLogarithmic         = log(      $massHaloModel                            ->($modelEntries))        ;
    my $massStellarModelLogarithmic      = log(10.0)*$massStellarModel                         ->($modelEntries)         ;
    my $massStellarErrorModelLogarithmic = log(10.0)*$massStellarCovarianceModel->diagonal(0,1)->($modelEntries) ->sqrt();

    # Interpolate observational data to model points.
    my $massStellarDataLogarithmicInterpolated      =  $spline     ->eval ($massHaloModelLogarithmic);
    my $massStellarErrorDataLogarithmicInterpolated =  $spline     ->deriv($massHaloModelLogarithmic)
	                                              *$splineError->eval ($massHaloModelLogarithmic);
    # Output the results to file if requested.
    if ( exists($options{'resultFile'}) ) {
	my $resultsFile = new PDL::IO::HDF5(">".$options{'resultFile'});
	$resultsFile->dataset('x'             )->set($massHaloModelLogarithmic                   );
	$resultsFile->dataset('y'             )->set($massStellarModelLogarithmic                );
	$resultsFile->dataset('yData'         )->set($massStellarErrorDataLogarithmicInterpolated);
    }

    # Compute the likelihood:
    if ( exists($options{'outputFile'}) ) {
	my $constraint;
	if ( nelem($massStellarModelLogarithmic) > 0 ) {
	    my $logLikelihood = -0.5*sum(($massStellarModelLogarithmic-$massStellarDataLogarithmicInterpolated)**2/($massStellarErrorModelLogarithmic**2+$massStellarErrorDataLogarithmicInterpolated**2));
	    $constraint->{'logLikelihood'} = $logLikelihood;
	} else {
	    $constraint->{'logLikelihood'} = -1.0e30;
	}
	$constraint->{'logLikelihoodVariance'} = 0.0;
	$constraint->{'label'                } = "cosmosStellarHaloMassRelationZ".$redshiftRange;
	# Output the constraint.
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
	open(oHndl,">".$options{'outputFile'});
	print oHndl $xmlOutput->XMLout($constraint);
	close(oHndl);
    }
    
    # Make a plot if requested.
    if ( exists($options{'plotFile'}) ) {
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileTeX, $plot);
	# Open a pipe to GnuPlot.
	($plotFileTeX = $options{'plotFile'}) =~ s/\.pdf$/.tex/;
	open($gnuPlot,"|gnuplot ");
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
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [1.0e08:1.0e12]\n";
	print $gnuPlot "set yrange [3.0e10:3.0e14]\n";
	print $gnuPlot "set title offset 0,-1 'Stellar mass-halo mass relation at \$z=".$redshiftMinimum."\$--\$".$redshiftMaximum."\$'\n";
	print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star\\, [\\mathrm{M}_\\odot]\$'\n";
	print $gnuPlot "set ylabel 'Halo mass; \$M_\\mathrm{200b}\\, [\\mathrm{M}_\\odot]\$'\n";
	print $gnuPlot "set pointsize 1.0\n";
	# Plot data points.
	&GnuPlot::PrettyPlots::Prepare_Dataset
	    (
	     \$plot,
	     $massStellarData,
	     $massHaloLowData,
	     y2           => $massHaloHighData,
	     style        => "filledCurve",
	     weight       => [3,1],
	     color        => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	     transparency => 0.75
	    );
	&GnuPlot::PrettyPlots::Prepare_Dataset
	    (
	     \$plot,
	     $massStellarData,
	     $massHaloMeanData,
	     style     => "line",
	     weight    => [3,1],
	     color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	     title     => "Leauthaud et al. (2012)"
	    );
	# Plot model points.
	my $massStellarModelLinear    = +10.0** $massStellarModel                                                                            ;
	my $massStellarModelErrorUp   = +10.0**($massStellarModel+$massStellarCovarianceModel->diagonal(0,1)->sqrt())-$massStellarModelLinear;
	my $massStellarModelErrorDown = -10.0**($massStellarModel-$massStellarCovarianceModel->diagonal(0,1)->sqrt())+$massStellarModelLinear;
	&GnuPlot::PrettyPlots::Prepare_Dataset
	    (
	     \$plot,
	     $massStellarModelLinear,
	     $massHaloModel,
	     errorLeft  => $massStellarModelErrorDown,
	     errorRight => $massStellarModelErrorUp  ,
	     style      => "point",
	     symbol     => [6,7], 
	     weight     => [3,1],
	     pointSize  => 0.5,
	     color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	     title      => "Galacticus"
	    );
	# Finalize plot.
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX,margin => 1);
    }
}

sub COSMOS2012_FastReject {
    # Perform fast-rejection likelihood analysis for the Leauthaud et al. (2012) COSMOS stellar mass-halo mass relation.
    my $galacticusFileName = shift();
    my $redshiftRange      = shift();
    my %options;
    (%options) = @_
	if ( scalar(@_) > 0 );

    # Validate redshift range.
    die("Galacticus::Constraints::StellarHaloMassRelation::COSMOS2012_FastReject(): redshiftRange must be 1, 2, or 3")
	if ( $redshiftRange < 1 || $redshiftRange > 3 );
    
    # Read observational data.
    my $data             = new PDL::IO::HDF5($ENV{'GALACTICUS_DATA_PATH'}."/static/observations/stellarHaloMassRelation/stellarHaloMassRelation_COSMOS_Leauthaud2012.hdf5");
    my $redshiftGroup    = $data         ->group  ('redshiftInterval'.$redshiftRange)       ;
    my $massStellarData  = $redshiftGroup->dataset('massStellar'                    )->get();
    my $massHaloMeanData = $redshiftGroup->dataset('massHaloMean'                   )->get();
    my $massHaloLowData  = $redshiftGroup->dataset('massHaloLow'                    )->get();
    my $massHaloHighData = $redshiftGroup->dataset('massHaloHigh'                   )->get();
    (my $redshiftMinimum, my $redshiftMaximum) = $redshiftGroup->attrGet('redshiftMinimum','redshiftMaximum');
    
    # Find a spline fit to the observed data, and compute the uncertainty in logarithm of halo mass.
    my $massStellarDataLogarithmic   = log($massStellarData );
    my $massHaloMeanDataLogarithmic  = log($massHaloMeanData);
    my $massHaloLowDataLogarithmic   = log($massHaloLowData );
    my $massHaloHighDataLogarithmic  = log($massHaloHighData);
    my $massHaloErrorDataLogarithmic = 0.5*($massHaloHighDataLogarithmic-$massHaloLowDataLogarithmic);
    my $spline                       = PDL::GSL::INTERP->init('cspline',$massHaloMeanDataLogarithmic,$massStellarDataLogarithmic  );
    my $splineError                  = PDL::GSL::INTERP->init('cspline',$massHaloMeanDataLogarithmic,$massHaloErrorDataLogarithmic);

    # Read model data.
    my $model                            = new PDL::IO::HDF5($galacticusFileName);
    my $analysis                         = $model   ->group  ('analyses'                                            )
    	                                            ->group  ('stellarHaloMassRelationLeauthaud2012z'.$redshiftRange)       ;
    my $massHaloModel                    = $analysis->dataset('massHalo'                                            )->get();
    my $massStellarModel                 = $analysis->dataset('massStellarLog10'                                    )->get();
    my $massStellarCovarianceModel       = $analysis->dataset('massStellarLog10Covariance'                          )->get();
    my $modelParameters                  = $model   ->group  ('Parameters'                                          )       ;
    my $rawParameters;
    if ( grep {$_ eq "rawXML"} $modelParameters->attrs() ) {
	(my $rawXML) = $modelParameters->attrGet('rawXML');
	my $xmlParameters = new XML::Simple();
	$rawParameters    = $xmlParameters->XMLin($rawXML);
    } else {
	die("Galacticus::Constraints::StellarHaloMassRelation::COSMOS2012_FastReject(): raw parameters XML not found");
    }
    my $binSelect                        = $rawParameters->{'stellarHaloMassRelationCosmos2012'}->{'redshiftInterval'.$redshiftRange}->{'binSelect'}->{'value'};
    my $massHaloModelLogarithmic         = log(      $massHaloModel                            ->($binSelect))        ;
    my $massStellarModelLogarithmic      = log(10.0)*$massStellarModel                         ->($binSelect)         ;
    my $massStellarErrorModelLogarithmic = log(10.0)*$massStellarCovarianceModel->diagonal(0,1)->($binSelect) ->sqrt();
    # Interpolate observational data to model points.
    my $massStellarDataLogarithmicInterpolated      =  $spline     ->eval ($massHaloModelLogarithmic);
    my $massStellarErrorDataLogarithmicInterpolated =  $spline     ->deriv($massHaloModelLogarithmic)
	*$splineError->eval ($massHaloModelLogarithmic);
    # If available, compute the scatter in the stellar mass.
    my $sigmaLog10MassStellar;
    if ( grep {$_ eq 'stellarHaloMassSquaredRelationLeauthaud2012z'.$redshiftRange} $model->group('analyses')->groups() ) {
	my $analysisSquared              = $model          ->group  ('analyses'                                                   )
	                                                   ->group  ('stellarHaloMassSquaredRelationLeauthaud2012z'.$redshiftRange)       ;
	my $massStellarSquaredModel      = $analysisSquared->dataset('massStellarLog10Squared'                                    )->get();
	my $sigmaLog10MassStellarSquared = +$massStellarSquaredModel->(($binSelect))
	                                   -$massStellarModel       ->(($binSelect))**2;
	$sigmaLog10MassStellar           = $sigmaLog10MassStellarSquared > 0.0 ? sqrt($sigmaLog10MassStellarSquared) : 0.0;
    } else {
	# If no mean square is available, set scatter to zero (which is always considered to be viable).
	$sigmaLog10MassStellar      = pdl 0.0;
    }
    # Compute the likelihood:
    if ( exists($options{'outputFile'}) ) {
	my $constraint;
	if ( isfinite($massStellarModelLogarithmic->((0))) ) {
	    my $toleranceSigma       = $rawParameters->{'stellarHaloMassRelationCosmos2012'}->{'redshiftInterval'.$redshiftRange}->{'toleranceSigma'  }->{'value'};
	    my $factorSigma          = $rawParameters->{'stellarHaloMassRelationCosmos2012'}->{'redshiftInterval'.$redshiftRange}->{'factorSigma'     }->{'value'};
	    my $likelihoodViable     = $rawParameters->{'stellarHaloMassRelationCosmos2012'}->{'redshiftInterval'.$redshiftRange}->{'likelihoodViable'}->{'value'};
	    my $scatterMaximum       = $rawParameters->{'stellarHaloMassRelationCosmos2012'}->{'redshiftInterval'.$redshiftRange}->{'scatterMaximum'  }->{'value'};
	    my $scatterSigma         = $rawParameters->{'stellarHaloMassRelationCosmos2012'}->{'redshiftInterval'.$redshiftRange}->{'scatterSigma'    }->{'value'};
	    my $offset               = ($massStellarModelLogarithmic-$massStellarDataLogarithmicInterpolated)**2/$massStellarErrorDataLogarithmicInterpolated**2;
	    my $logLikelihoodMass    = $offset                < $toleranceSigma**2 ? 0.0 : sclr(-0.5*( $offset               -$toleranceSigma**2)*$factorSigma  **2);
	    my $logLikelihoodScatter = $sigmaLog10MassStellar < $scatterMaximum    ? 0.0 : sclr(-0.5*(($sigmaLog10MassStellar-$scatterMaximum   )/$scatterSigma)**2);
	    my $logLikelihood        = 
		(
		 $offset                < $toleranceSigma**2 
		 &&
		 $sigmaLog10MassStellar < $scatterMaximum
		)
		?
		$likelihoodViable 
		:
		+$logLikelihoodMass
		+$logLikelihoodScatter;
	    $constraint->{'logLikelihood'} = $logLikelihood;
	} else {
	    $constraint->{'logLikelihood'} =  -1.0e30;
	}
	$constraint->{'logLikelihoodVariance'} = 0.0;
	$constraint->{'label'                } = "cosmosStellarHaloMassRelationZ".$redshiftRange."FastReject";
	# Output the constraint.
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
	open(oHndl,">".$options{'outputFile'});
	print oHndl $xmlOutput->XMLout($constraint);
	close(oHndl);
    }

    # Make a plot if requested.
    if ( exists($options{'plotFile'}) ) {
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileTeX, $plot);
	# Open a pipe to GnuPlot.
	($plotFileTeX = $options{'plotFile'}) =~ s/\.pdf$/_$binSelect.tex/;
	open($gnuPlot,"|gnuplot ");
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
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [1.0e08:1.0e12]\n";
	print $gnuPlot "set yrange [3.0e10:3.0e14]\n";
	print $gnuPlot "set title offset 0,-1 'Stellar mass-halo mass relation at \$z=".$redshiftMinimum."\$--\$".$redshiftMaximum."\$'\n";
	print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star\\, [\\mathrm{M}_\\odot]\$'\n";
	print $gnuPlot "set ylabel 'Halo mass; \$M_\\mathrm{200b}\\, [\\mathrm{M}_\\odot]\$'\n";
	print $gnuPlot "set pointsize 1.0\n";
	# Plot data points.
	&GnuPlot::PrettyPlots::Prepare_Dataset
	    (
	     \$plot,
	     $massStellarData,
	     $massHaloLowData,
	     y2           => $massHaloHighData,
	     style        => "filledCurve",
	     weight       => [3,1],
	     color        => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	     transparency => 0.75
	    );
	&GnuPlot::PrettyPlots::Prepare_Dataset
	    (
	     \$plot,
	     $massStellarData,
	     $massHaloMeanData,
	     style     => "line",
	     weight    => [3,1],
	     color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	     title     => "Leauthaud et al. (2012)"
	    );
	# Plot model points.
	my $massStellarModelLinear    = +10.0** $massStellarModel                                                                            ;
	my $massStellarModelErrorUp   = +10.0**($massStellarModel+$massStellarCovarianceModel->diagonal(0,1)->sqrt())-$massStellarModelLinear;
	my $massStellarModelErrorDown = -10.0**($massStellarModel-$massStellarCovarianceModel->diagonal(0,1)->sqrt())+$massStellarModelLinear;
	&GnuPlot::PrettyPlots::Prepare_Dataset
	    (
	     \$plot,
	     $massStellarModelLinear                 ->($binSelect),
	     $massHaloModel                          ->($binSelect),
	     errorLeft  => $massStellarModelErrorDown->($binSelect),
	     errorRight => $massStellarModelErrorUp  ->($binSelect),
	     style      => "point",
	     symbol     => [6,7], 
	     weight     => [3,1],
	     pointSize  => 0.5,
	     color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	     title      => "Galacticus"
	    );
	# Finalize plot.
	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX,margin => 1);
    }
}

1;
