# Contains a Perl module which implements various useful functionality for constraints based on the stellar mass vs. halo mass
# relation.

package Galacticus::Constraints::StellarHaloMassRelation;
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
use PDL::GSL::INTERP;

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
    my $data             = new PDL::IO::HDF5("data/observations/stellarHaloMassRelation/stellarHaloMassRelation_COSMOS_Leauthaud2012.hdf5");
    my $redshiftGroup    = $data         ->group  ('redshiftInterval'.$redshiftRange)       ;
    my $massStellarData  = $redshiftGroup->dataset('massStellar'                    )->get();
    my $massHaloMeanData = $redshiftGroup->dataset('massHaloMean'                   )->get();
    my $massHaloLowData  = $redshiftGroup->dataset('massHaloLow'                    )->get();
    my $massHaloHighData = $redshiftGroup->dataset('massHaloHigh'                   )->get();
    (my $redshiftMinimum, my $redshiftMaximum) = $redshiftGroup->attrGet('redshiftMinimum','redshiftMaximum');
    
    # Find a spline fit to the observed data, and compute the uncertainty in logairhtm of halo mass.
    my $massStellarDataLogarithmic   = log($massStellarData );
    my $massHaloMeanDataLogarithmic  = log($massHaloMeanData);
    my $massHaloLowDataLogarithmic   = log($massHaloLowData );
    my $massHaloHighDataLogarithmic  = log($massHaloHighData);
    my $massHaloErrorDataLogarithmic = 0.5*($massHaloHighDataLogarithmic-$massHaloLowDataLogarithmic);
    my $spline                       = PDL::GSL::INTERP->init('cspline',$massHaloMeanDataLogarithmic,$massStellarDataLogarithmic  );
    my $splineError                  = PDL::GSL::INTERP->init('cspline',$massHaloMeanDataLogarithmic,$massHaloErrorDataLogarithmic);

    # Read model data.
    my $model = new PDL::IO::HDF5($galacticusFileName);
    my $analysis                         = $model   ->group ('analyses'                                            )
    	                                            ->group ('stellarHaloMassRelationLeauthaud2012z'.$redshiftRange)       ;
    my $massHaloModel                    = $analysis->dataset('massHalo'                                           )->get();
    my $massStellarModel                 = $analysis->dataset('massStellar'                                        )->get();
    my $massStellarCovarianceModel       = $analysis->dataset('massStellarCovariance'                              )->get();
    my $modelEntries                     = 
	which(
	    ($massStellarModel                          >  0.0                      ) 
	    &
	    ($massStellarCovarianceModel->diagonal(0,1) >  0.0                      )
	    &
	    ($massHaloModel                             >= $massHaloMeanData->(( 0)))
	    &
	    ($massHaloModel                             <= $massHaloMeanData->((-1)))
	);
    my $massHaloModelLogarithmic         = log($massHaloModel   ->($modelEntries));
    my $massStellarModelLogarithmic      = log($massStellarModel->($modelEntries));
    my $massStellarErrorModelLogarithmic = $massStellarCovarianceModel->diagonal(0,1)->($modelEntries)->sqrt()/$massStellarModel->($modelEntries);

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
	    my $logLikelihood = -0.5*($massStellarModelLogarithmic-$massStellarDataLogarithmicInterpolated)**2/($massStellarErrorModelLogarithmic**2+$massStellarErrorDataLogarithmicInterpolated**2);
	    $constraint->{'logLikelihood'} = $logLikelihood->sclr();
	} else {
	    $constraint->{'logLikelihood'} = -1.0e30;
	}
	$constraint->{'logLikelihoodVariance'} = 0.0;
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
	# Plot model points.
	my $massHaloRegion    = $massHaloLowData->append($massHaloHighData->(-1:0));
	my $massStellarRegion = $massStellarData->append($massStellarData ->(-1:0));
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
	&GnuPlot::PrettyPlots::Prepare_Dataset
	    (
	     \$plot,
	     $massStellarModel,
	     $massHaloModel,
	     errorLeft  => $massStellarCovarianceModel->diagonal(0,1)->sqrt(),
	     errorRight => $massStellarCovarianceModel->diagonal(0,1)->sqrt(),
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
