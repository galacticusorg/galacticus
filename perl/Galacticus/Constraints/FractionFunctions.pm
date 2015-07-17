# Contains a Perl module which implements construction of a variety of fraction functions from
# Galacticus when fitting to constraints.

package FractionFunctions;
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
use PDL::Constants qw(PI);
use PDL::MatrixOps;
use PDL::IO::HDF5;
use Math::SigFigs;
use Data::Dumper;
use LaTeX::Encode;
use XML::Simple;
use Scalar::Util 'reftype';
use Astro::Cosmology;
require Galacticus::HDF5;
require Galacticus::Constraints::Covariances;

sub Construct {
    # Construct a fraction function from Galacticus for constraint purposes.
    my %arguments = %{$_[0]};
    my $config    =   $_[1];
    
    # Create data structure to read the results.
    my $galacticus;
    $galacticus->{'file' } = $config->{'galacticusFile'};
    $galacticus->{'store'} = 0;

    # Make the label LaTeX compliant.
    $config->{'observationLabel'} = latex_encode($config->{'observationLabel'});
    $config->{'observationLabel'} =~ s/\\/\\\\/g;

    # Construct upper and lower error bars if not specified.
    if ( exists($config->{'covariance'}) ) {
	my $diagonalError = $config->{'covariance'}->diagonal(0,1)->copy();
	foreach ( 'yLowerError', 'yUpperError' ) {
	    $config->{$_} = sqrt($diagonalError)
		unless ( exists($config->{$_}) );
	}
    }

    # Check array sizes match.
    foreach ( 'y', 'yLowerError', 'yUpperError' ) {
	die('MassFunctions::Construct(): size of '.$_.' array differs from size of x array')
	    if ( exists($config->{$_}) && nelem($config->{$_}) != nelem($config->{'x'}) );
    }
    
    # Map from any logarithmic scaling to linear scaling.
    $config->{'x'          } = +10.0** $config->{'x'          }
        if ( $config->{'xScaling'} eq "log10" );
    $config->{'yUpperError'} = +10.0**($config->{'y'          }+$config->{'yUpperError'})-10.0**$config->{'y'          }
        if ( $config->{'yScaling'} eq "log10" );
    $config->{'yLowerError'} = -10.0**($config->{'y'          }-$config->{'LowerError' })+10.0**$config->{'y'          }
        if ( $config->{'yScaling'} eq "log10" );
    $config->{'y'          } = +10.0** $config->{'y'          }
        if ( $config->{'yScaling'} eq "log10" );

    # Construct a convariance matrix if one is needed.
    my $error = 0.5*($config->{'yUpperError'}+$config->{'yLowerError'});
    unless ( exists($config->{'covariance'}) ) {
	$config->{'covariance'}                 = pdl zeroes(nelem($error),nelem($error));
	$config->{'covariance'}->diagonal(0,1) .= $error**2;
    }

    # Get logarithmic bins in mass.
    my $xBins = log10($config->{'x'});
    my $xBinWidths;
    $xBinWidths = log10($config->{'xWidth'})
	if ( exists($config->{'xWidth'}) );

    # Model results.
    my $yGalacticus;
    my $errorGalacticus;
    my $covarianceGalacticus;

    # Determine if the model file contains a pre-computed mass function.
    &HDF5::Open_File($galacticus);
    my @analysisGroups = $galacticus->{'hdf5File'}->group('analysis')->groups();
    if ( grep {$_ eq $config->{'analysisLabel'}} @analysisGroups ) {
	$yGalacticus          = $galacticus->{'hdf5File'}->group('analysis')->group($config->{'analysisLabel'})->dataset('fraction'          )->get();
	$covarianceGalacticus = $galacticus->{'hdf5File'}->group('analysis')->group($config->{'analysisLabel'})->dataset('fractionCovariance')->get();
	$errorGalacticus      = sqrt($covarianceGalacticus->diagonal(0,1));
    }

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
		    	my $multiplier = $discrepancyFile->dataset('multiplicative')->get();
		    	$yGalacticus          *= $multiplier;
		    	$errorGalacticus      *= $multiplier;
		    	$covarianceGalacticus .= $covarianceGalacticus*outer($multiplier,$multiplier);
		    }		    
		    if ( $dataset eq "multiplicativeCovariance" ) {
		    	# Adjust the model accordingly.
		    	my $covarianceMultiplier = $discrepancyFile->dataset('multiplicativeCovariance')->get();
		    	$covarianceGalacticus   += $covarianceMultiplier*outer($yGalacticus,$yGalacticus);
		    }		    
		    if ( $dataset eq "additive" ) {
		    	# Read the additive discrepancy
		    	my $addition = $discrepancyFile->dataset('additive')->get();
		    	# Adjust the model accordingly.
		    	$yGalacticus     += $addition;
		    }
		    if ( $dataset eq "additiveCovariance" ) {
		    	# Read the covariance of the discrepancy.
		    	my $covariance = $discrepancyFile->dataset('additiveCovariance')->get();
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
	$resultsFile->dataset('x'             )->set($xBins                 );
	$resultsFile->dataset('y'             )->set($yGalacticus           );
	$resultsFile->dataset('covariance'    )->set($covarianceGalacticus  );
	$resultsFile->dataset('yData'         )->set($config->{'y'         });
	$resultsFile->dataset('covarianceData')->set($config->{'covariance'});
    }

    # Compute the likelihood:
    if ( exists($arguments{'outputFile'}) ) {
	# Construct the full covariance matrix, which is the covariance matrix of the observations
	# plus that of the model.
	my $fullCovariance = $config->{'covariance'}+$covarianceGalacticus;
	# If the range of masses over which to constrain is limited, set the model fraction function outside of that range equal to
	# the observed fraction function.
	my $yGalacticusLimited = $yGalacticus->copy();
	if ( exists($config->{'constraintMassMinimum'}) ) {
	    my $noConstraint = which($config->{'x'} < $config->{'constraintMassMinimum'});
	    $yGalacticusLimited->($noConstraint) .= $config->{'y'}->($noConstraint)
		if ( nelem($noConstraint) > 0 );
	}
	if ( exists($config->{'constraintMassMaximum'}) ) {
	    my $noConstraint = which($config->{'x'} > $config->{'constraintMassMaximum'});
	    $yGalacticusLimited->($noConstraint) .= $config->{'y'}->($noConstraint)
		if ( nelem($noConstraint) > 0 );
	}
	# Compute the likelihood.
	my $constraint;
	my $logDeterminant;
	my $offsets;
	my $inverseCovariance;
	my $logLikelihood = &Covariances::ComputeLikelihood($yGalacticusLimited,$config->{'y'},$fullCovariance, determinant => \$logDeterminant, inverseCovariance => \$inverseCovariance, offsets => \$offsets, quiet => $arguments{'quiet'});
	$constraint->{'logLikelihood'} = $logLikelihood;
	# Find the Jacobian of the log-likelihood with respect to the model mass function.
	my $jacobian = pdl zeroes(1,nelem($yGalacticus));
	for(my $i=0;$i<nelem($yGalacticus);++$i) {
	    $jacobian->((0),($i)) .= sum($inverseCovariance->(($i),:)*$offsets);
	}
	# Compute the variance in the log-likelihood due to errors in the model.
	my $logLikelihoodVariance = transpose($jacobian) x $covarianceGalacticus x $jacobian;
	$constraint->{'logLikelihoodVariance'} = $logLikelihoodVariance->sclr();
	# Output the constraint.
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
	open(oHndl,">".$arguments{'outputFile'});
	print oHndl $xmlOutput->XMLout($constraint);
	close(oHndl);
    }

    # Create a plot of the fraction function.
    if ( exists($arguments{'plotFile'}) ) {
	require GnuPlot::PrettyPlots;
	require GnuPlot::LaTeX;
	require XMP::MetaData;
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
	print $gnuPlot "set key at screen 0.4,0.2\n";
	print $gnuPlot "set key left\n";
	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [".$config->{'xRange'}."]\n";
	print $gnuPlot "set yrange [".$config->{'yRange'}."]\n";
	print $gnuPlot "set title '".$config->{'title'}."'\n";
	print $gnuPlot "set xlabel '".$config->{'xLabel'}."'\n";
	print $gnuPlot "set ylabel '".$config->{'yLabel'}."'\n";
	&PrettyPlots::Prepare_Dataset(\$plot,
				      $config->{'x'},$config->{'y'},
				      errorUp   => $error,
				      errorDown => $error,
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
				      title     => $config->{'observationLabel'}
	    );
	&PrettyPlots::Prepare_Dataset(\$plot,
				      $config->{'x'},$yGalacticus,
				      errorUp   => $errorGalacticus,
				      errorDown => $errorGalacticus,
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'redYellow'},
				      title     => "Galacticus"
	    );
	&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);

    }

}

1;
