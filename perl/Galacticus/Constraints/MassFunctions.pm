# Contains a Perl module which implements construction of a variety of mass functions from
# Galacticus when fitting to constraints.

package MassFunctions;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath  = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw(PI);
use PDL::MatrixOps;
use Math::SigFigs;
use Data::Dumper;
use LaTeX::Encode;
use XML::Simple;
use Scalar::Util 'reftype';
use Switch;
use Astro::Cosmology;
require Galacticus::HDF5;
require Galacticus::StellarMass;
require Galacticus::HIGasMass;
require Galacticus::GasMass;
require Galacticus::Constraints::Covariances;
require Stats::Histograms;

sub Construct {
    # Construct a mass function from Galacticus for constraint purposes.
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
    
    # Convert from per log10(M) to per log(M) if necessary.
    if ( $config->{'yIsPer'} eq "log10" ) {
	$config->{'y'          } /= log(10.0);
	$config->{'yUpperError'} /= log(10.0);
	$config->{'yLowerError'} /= log(10.0);
    }

    # Construct a convariance matrix if one is needed.
    my $error = 0.5*($config->{'yUpperError'}+$config->{'yLowerError'});
    unless ( exists($config->{'covariance'}) ) {
	$config->{'covariance'}                 = pdl zeroes(nelem($error),nelem($error));
	$config->{'covariance'}->diagonal(0,1) .= $error**2;
    }

    # Get logarithmic bins in mass.
    my $xBins = log10($config->{'x'});

    # Model results.
    my $yGalacticus;
    my $errorGalacticus;
    my $covarianceGalacticus;

    # Determine if the model file contains a pre-computed mass function.
    &HDF5::Open_File($galacticus);
    my $gotModelMassFunction = 0;
    my @rootGroups = $galacticus->{'hdf5File'}->groups();
    if ( grep {$_ eq "analysis"} @rootGroups ) {
    	my @analysisGroups = $galacticus->{'hdf5File'}->group('analysis')->groups();
    	if ( grep {$_ eq $config->{'analysisLabel'}} @analysisGroups ) {
    	    $gotModelMassFunction  = 1;
    	    $yGalacticus          = $galacticus->{'hdf5File'}->group('analysis')->group($config->{'analysisLabel'})->dataset('massFunction'          )->get();
    	    $covarianceGalacticus = $galacticus->{'hdf5File'}->group('analysis')->group($config->{'analysisLabel'})->dataset('massFunctionCovariance')->get();
    	    $errorGalacticus      = sqrt($covarianceGalacticus->diagonal(0,1));
    	}
    }
    # Read galaxy data and construct mass function if necessary.
    if ( $gotModelMassFunction == 0 ) {
	$galacticus->{'tree'} = "all";
	&HDF5::Get_Parameters($galacticus);
	&HDF5::Count_Trees  ($galacticus                      );
	&HDF5::Select_Output($galacticus,$config->{'redshift'});
	&HDF5::Get_Dataset  ($galacticus,['mergerTreeWeight',$config->{'massType'}]);
	my $dataSets = $galacticus->{'dataSets'};
	my $weight   = $dataSets->{'mergerTreeWeight'};
	# Find cosmological conversion factors.
	my $cosmologyObserved = Astro::Cosmology->new(omega_matter => $config->{'omegaMatterObserved'}, omega_lambda =>  $config->{'omegaDarkEnergyObserved'}, h0 =>  $config->{'hubbleConstantObserved'});
	my $cosmologyModel    = Astro::Cosmology->new(omega_matter => $galacticus->{'parameters'}->{'Omega_Matter'}, omega_lambda => $galacticus->{'parameters'}->{'Omega_DE'}, h0 => $galacticus->{'parameters'}->{'H_0'});
	my $cosmologyScalingMass        ;
	my $cosmologyScalingMassFunction;
	if ( $config->{'cosmologyScalingMass'}->atstr(0) eq 'none' ) {
	    # Nothing to do.
	} elsif ( $config->{'cosmologyScalingMass'}->atstr(0) eq 'luminosity' ) {
	    $cosmologyScalingMass = ($cosmologyObserved->luminosity_distance($config->{'redshift'})/$cosmologyModel->luminosity_distance($config->{'redshift'}))**2;
	} else {
	    die('MassFunctions::Construct: unrecognized cosmology scaling');
	}
	if ( $config->{'cosmologyScalingMassFunction'}->atstr(0) eq 'none' ) {
	    # Nothing to do.
	} elsif ( $config->{'cosmologyScalingMassFunction'}->atstr(0) eq 'inverseComovingVolume' ) {
	    $cosmologyScalingMassFunction =
		($cosmologyModel->comoving_distance($config->{'redshift'})/$cosmologyObserved->comoving_distance($config->{'redshift'}))**2
		/($cosmologyModel->h0($config->{'redshift'})/$cosmologyObserved->h0($config->{'redshift'}))**2;
	} else {
	    die('MassFunctions::Construct: unrecognized cosmology scaling');
	}
	$dataSets->{$config->{'massType'}} *= $cosmologyScalingMass        ;
	$weight                            *= $cosmologyScalingMassFunction;
	# Map masses.
	my $logarithmicMass;
	if ( exists($config->{'massMap'}) ) {
	    $logarithmicMass = &{$config->{'massMap'}}($config,$galacticus);
	} else {
	    $logarithmicMass = log10($dataSets->{$config->{'massType'}});
	}
	# Add random Gaussian errors to the masses.
	my $sigma = pdl ones(nelem($logarithmicMass));
	if ( reftype($config->{'massErrorRandomDex'}) eq "CODE" ) {
	    $sigma .= &{$config->{'massErrorRandomDex'}}($logarithmicMass,$galacticus);
	} else {
	    $sigma *= $config->{'massErrorRandomDex'};
	}
	# Add systematic shift to masses (i.e. systematic errors nuisance parameter).
	if ( exists($config->{'systematicParameter'}) ) {
	    my $systematicOffset = pdl zeroes(nelem($logarithmicMass));
	    for(my $i=0;$i<=$config->{'systematicOrder'};++$i) {
		my $parameterName = $config->{'systematicParameter'}.$i;
		$systematicOffset += 
		    $galacticus->{'parameters'}->{$parameterName}
		    *($logarithmicMass-$config->{'systematicZeroPoint'})**$i;
	    }
	    $logarithmicMass += $systematicOffset;
	}
	# Construct the mass function. 
	my %options = 
	    (
	     differential   => 1,
	     gaussianSmooth => $sigma
	    );
	($yGalacticus,$errorGalacticus,$covarianceGalacticus) = &Histograms::Histogram($xBins,$logarithmicMass,$weight,%options);
	# Convert model mass function from per log10(M) to per log(M).
	$yGalacticus          /= log(10.0);
	$errorGalacticus      /= log(10.0);
	$covarianceGalacticus /= log(10.0)**2;
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
	my $results;
	@{$results->{'x'             }} = $xBins                 ->list();
	@{$results->{'y'             }} = $yGalacticus           ->list();
	@{$results->{'error'         }} = $errorGalacticus       ->list();
	@{$results->{'covariance'    }} = $covarianceGalacticus  ->list();
	@{$results->{'yData'         }} = $config->{'y'}         ->list();
	@{$results->{'covarianceData'}} = $config->{'covariance'}->list();
	my $xmlOut = new XML::Simple (RootName=>"results", NoAttr => 1);;
	# Output the parameters to file.
	open(pHndl,">".$arguments{'resultFile'});
	print pHndl $xmlOut->XMLout($results);
	close pHndl;
    }

    # Output accuracy to file if requested.
    if ( exists($arguments{'accuracyFile'}) ) {
	my $results;
	@{$results->{'x'         }} = $xBins          ->list();
	@{$results->{'yModel'    }} = $yGalacticus    ->list();
	@{$results->{'yData'     }} = $config->{'y'}  ->list();
	@{$results->{'errorModel'}} = $errorGalacticus->list();
	@{$results->{'errorData' }} = $error          ->list();
	my $xmlOut = new XML::Simple (RootName=>"accuracy", NoAttr => 1);;
	# Output the parameters to file.
	open(pHndl,">".$arguments{'accuracyFile'});
	print pHndl $xmlOut->XMLout($results);
	close pHndl;
    }

    # Compute the likelihood:
    if ( exists($arguments{'outputFile'}) ) {
	# For any bins in which the model mass function is zero, we extrapolate from adjacent bins. This is necessary in order to
	# permit a finite likelihood to be computed. Such models will always have a very low likelihood, so the details of the
	# extrapolation should not matter.
	my $yDefined                 = pdl zeroes(nelem($yGalacticus));
	my $nonZero                  = which($yGalacticus > 0.0);
	$yDefined->($nonZero)       .= 1.0;
	my $logYGalacticus           = $yGalacticus->copy();
	$logYGalacticus->($nonZero) .= log($yGalacticus->index($nonZero));
	for(my $i=2;$i<nelem($yGalacticus);++$i) {
	    if ( $yDefined->(($i)) == 0.0 && $yDefined->(($i-1)) > 0.0 && $yDefined->(($i-2)) > 0.0 ) {
		$logYGalacticus->(($i)) .= $logYGalacticus->(($i-1))*2.0-$logYGalacticus->(($i-2));
		$yDefined      ->(($i)) .= 1.0;
	    }
	}
	for(my $i=nelem($yGalacticus)-3;$i>=0;--$i) {
	    if ( $yDefined->(($i)) == 0.0 && $yDefined->(($i+1)) > 0.0 && $yDefined->(($i+2)) > 0.0 ) {
		$logYGalacticus->(($i)) .= $logYGalacticus->(($i+1))*2.0-$logYGalacticus->(($i+2));
		$yDefined      ->(($i)) .= 1.0;
	    }
	}
	$yGalacticus .= exp($logYGalacticus);

	# Construct the full covariance matrix, which is the covariance matrix of the observations
	# plus that of the model.
	my $fullCovariance        = $config->{'covariance'}+$covarianceGalacticus;
	# Compute the likelihood.
	my $constraint;
	my $logDeterminant;
	my $logLikelihood = &Covariances::ComputeLikelihood($yGalacticus,$config->{'y'},$fullCovariance, determinant => \$logDeterminant, quiet => $arguments{'quiet'});
	$constraint->{'logLikelihood'} = $logLikelihood;
	# Output the constraint.
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
	open(oHndl,">".$arguments{'outputFile'});
	print oHndl $xmlOutput->XMLout($constraint);
	close(oHndl);
    }

    # Create a plot of the mass function.
    if ( exists($arguments{'plotFile'}) ) {
	require GnuPlot::PrettyPlots;
	require GnuPlot::LaTeX;
	require XMP::MetaData;
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	($plotFileEPS = $arguments{'plotFile'}) =~ s/\.pdf$/.eps/;
	open($gnuPlot,"|gnuplot ");#1>/dev/null 2>&1");
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
