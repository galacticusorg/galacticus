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
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::GSL::RNG;
use Data::Dumper;
require Galacticus::HDF5;
require Galacticus::Options;
require Galacticus::Constraints::Covariances;
require Galacticus::Constraints::LocalGroupDatabase;

# Compute likelihood (and make a plot) for a Galacticus model given the stellar mass function of Milky Way satellite galaxies.
# Andrew Benson (13-November-2014)

# Get name of input and output files.
die("milkyWaySatellitesStellarMassFunction.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $self               = $0;
my $galacticusFileName = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet   => 0         ,
     central => "MilkyWay"
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Define constants.
my $Pi = pdl 3.1415927;

# Define distance dataset.
my $distanceDataset = "distance".$arguments{'central'};

# Extract data from the Local Group database.
(my $names, my $masses, my $distanceModuli, my $vBandApparentMagnitudes, my $vBandAbsoluteMagnitudes, my $detectionEfficiency, my $publication, my $distances) 
    = &LocalGroup::Select(["name","massStellar","distanceModulus","magnitudeApparentV","magnitudeAbsoluteV","detectionEfficiencyHalfLight","publication",$distanceDataset], excludeCentral => 1, );
# Find the mass conversion.
my $massConversion     = $masses   ->{'meta'}->{'unitsInSI'}/1.9891e30;
# Find the distance conversion.
my $distanceConversion = $distances->{'meta'}->{'unitsInSI'}/3.0860e22;
# Estimate mass error from photometric uncertainty.
$masses->{'error'} = 0.4*log(10.0)*$masses->{'value'}*sqrt($vBandApparentMagnitudes->{'error'}**2+$distanceModuli->{'error'}**2);
my $haveAbsoluteMagnitude =which($vBandAbsoluteMagnitudes->{'error'} > 0.0);
$masses->{'error'}->($haveAbsoluteMagnitude) .= 0.4*log(10.0)*$masses->{'value'}->($haveAbsoluteMagnitude)*$vBandAbsoluteMagnitudes->{'error'}->($haveAbsoluteMagnitude);

# Convert to Solar masses.
$masses   ->{'value'    } *= $massConversion;
$masses   ->{'error'    } *= $massConversion;
$masses   ->{'errorLow' }  = $masses->{'error'};
$masses   ->{'errorHigh'}  = $masses->{'error'};
# Convert to Mpc.
$distances->{'value'    } *= $distanceConversion;
# Read model mass function.
my $galacticus;
$galacticus->{'file' } = $galacticusFileName;
&HDF5::Get_Parameters($galacticus);
my $massFunctionGroup = $galacticus->{'hdf5File'}->group('analysis')->group(lcfirst($arguments{'central'}).'MassFunction');
my $model;
($model->{'haloRadius'}, $model->{'logLikelihood'}, $model->{'haloCount'}) = $massFunctionGroup->attrGet('haloRadius','logLikelihood' ,'haloCount');
foreach ( "massStellar", "massFunctionCumulative", "massFunctionCumulativeVariance", "massFunctionCumulativePoissonVariance" ) {
    $model->{$_} = $massFunctionGroup->dataset($_)->get();
}

# Extract V-band mass-to-light ratio to assume for Local Group satellites.
my $massToLightRatio = exists($galacticus->{'parameters'}->{'localGroupSatellitesMassToLightBandV'}) ? exists($galacticus->{'parameters'}->{'localGroupSatellitesMassToLightBandV'}) : pdl 1.0;
# Scale observed masses for mass-to-light ratio.
$masses->{'value'    } *= $massToLightRatio;
$masses->{'errorLow' } *= $massToLightRatio;
$masses->{'errorHigh'} *= $massToLightRatio;
# Find members.
my $members = which(
    ($distances->{'value'} <= $model->{'haloRadius'})
    &
    ($distances->{'value'} >  0.0                   )
    &
    ($masses   ->{'value'} >  0.0                   )
    );
# Reduce to members.
$masses                 ->{'value'    } = $masses                 ->{'value'    }->($members);
$masses                 ->{'error'    } = $masses                 ->{'error'    }->($members);
$masses                 ->{'errorLow' } = $masses                 ->{'errorLow' }->($members);
$masses                 ->{'errorHigh'} = $masses                 ->{'errorHigh'}->($members);
$distanceModuli         ->{'value'    } = $distanceModuli         ->{'value'    }->($members);
$vBandApparentMagnitudes->{'value'    } = $vBandApparentMagnitudes->{'value'    }->($members);
# Determine cumulative order.
my $number              = pdl sequence(nelem($masses->{'value'}))+1;
$masses->{'order'    }  = $masses->{'value'}->qsorti();
# Find completeness, dependent on system.
if ( $arguments{'central'} eq "MilkyWay" ) {
    # Determine completeness correction for model using the model of Tollerud et al. (2008;
    # http://adsabs.harvard.edu/abs/2008ApJ...688..277T).
    my $tollerudA              = pdl 0.6;     # Fitting coefficient from Tollerud et al.
    my $tollerudB              = pdl 5.23;    # Fitting coefficient from Tollerud et al.
    my $vBandAbsoluteSolar     = pdl 4.80;    # Vega system.
    my $modernSurveySolidAngle = pdl 3.0*$Pi; # PanSTARRS survey.
    my $milkyWaySystemRadius   = pdl 0.417;   # Milky Way system radius in Mpc.
    my $vBandAbsoluteMagnitude = -2.5*log10($model->{'massStellar'}/$massToLightRatio)+$vBandAbsoluteSolar;
    my $milkyWaySystemVolume   = 4.0*$Pi*$milkyWaySystemRadius**3/3.0;
    my $volumeLimit            = 10.0**(-$tollerudA*$vBandAbsoluteMagnitude-$tollerudB);
    my $fullVolume             = which($volumeLimit > $milkyWaySystemVolume);
    my $modernSurveysOnly      = which($vBandAbsoluteMagnitude > -10.0);
    # Compute completeness.
    $model->{'completeness'}   = $volumeLimit/$milkyWaySystemVolume;
    # Limit completeness to maximum volume of the Milky Way system.
    $model->{'completeness'}->($fullVolume       ) .= 1.0;
    # Limit completeness to the survey region of the sky.
    $model->{'completeness'}->($modernSurveysOnly) *= $modernSurveySolidAngle/4.0/$Pi;
    # Handle new DES dwarfs (10-March-2015).
    my $desSurveySolidAngle = pdl 0.038785095; # sterradians, corresponding to 1600 deg^2 of survey not overlapping SDSS.
    for(my $i=0;$i<scalar(@{$publication->{'value'}});++$i) {
	$model->{'completeness'}->(($i)) .= $detectionEfficiency->{'value'}->(($i))*$desSurveySolidAngle
 	    if (
		defined($publication->{'value'}->[$i]->{'value'})
		&&
		$publication        ->{'value'}->[$i]->{'value'} eq "Koposov et al. (2015; ApJ)"
		&&
		$detectionEfficiency->{'value'}->(($i))           > 0.0 
	    );
    }   
    # Compute cumulative completeness.
    my $massFunctionDifferential = $model->{'massFunctionCumulative'}->copy();
    $massFunctionDifferential->(0:-2)        -= $model->{'massFunctionCumulative'}->(1:-1);
    $massFunctionDifferential                *= $model->{'completeness'};
    my $massFunctionIncomplete                = $massFunctionDifferential->(-1:0)->cumusumover()->(-1:0);
    my $zeroEntries                           = which($model->{'massFunctionCumulative'} <= 0.0);
    $model->{'completeness'}                 .= $massFunctionIncomplete/$model->{'massFunctionCumulative'};
    $model->{'completeness'}->($zeroEntries) .= 1.0;
} else {
    $model->{'completeness'}                  = pdl ones(nelem($model->{'massStellar'}));
}

# Create a plot of the mass function function.
if ( exists($arguments{'plotFile'}) ) {
    require GnuPlot::PrettyPlots;
    require GnuPlot::LaTeX;
    require XMP::MetaData;
    # Find label for central galaxy.
    my $centralLabel;
    $centralLabel = "Milky Way"
	if ( $arguments{'central'} eq "MilkyWay" );
    $centralLabel = "M31"
	if ( $arguments{'central'} eq "M31"      );
    # Declare variables for GnuPlot;
    my ($gnuPlot, $plotFileEPS, $plot);
    # Open a pipe to GnuPlot.
    ($plotFileEPS = $arguments{'plotFile'}) =~ s/\.pdf$/.eps/;
    open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
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
    print $gnuPlot "set xrange [100.0:3.0e10]\n";
    print $gnuPlot "set yrange [0.3:100.0]\n";
    print $gnuPlot "set title 'Cumulative stellar mass function of ".$centralLabel." satellite galaxies'\n";
    print $gnuPlot "set xlabel 'Stellar mass; \$M_\\star\$ [\$M_\\odot\$]'\n";
    print $gnuPlot "set ylabel 'Cumulative number; \$ N ( > M_\\star) \$ []'\n";
    &PrettyPlots::Prepare_Dataset(\$plot,
				  $masses              ->{'value'    }->($masses->{'order'}),
				  $number                             ->(-1:0              ),
				  errorLeft  => $masses->{'errorLow' }->($masses->{'order'}),
				  errorRight => $masses->{'errorHigh'}->($masses->{'order'}),
				  style      => "point",
				  symbol     => [6,7], 
				  weight     => [3,1],
				  color      => $PrettyPlots::colorPairs{'cornflowerBlue'},
				  title      => "Observed"
	);
    &PrettyPlots::Prepare_Dataset(\$plot,
				  +$model->{'massStellar'           },
				  +$model->{'massFunctionCumulative'}
				  *$model->{'completeness'          },
				  style      => "line",
				  weight     => [3,1],
				  color      => $PrettyPlots::colorPairs{'redYellow'},
				  title      => "Galacticus"
       );
    &PrettyPlots::Prepare_Dataset(\$plot,
				  +$model->{'massStellar'           },
				  +(
				      +     $model->{'massFunctionCumulative'        }
				      +sqrt($model->{'massFunctionCumulativeVariance'})
				  )
				  *$model->{'completeness'          },
				  style      => "line",
				  weight     => [3,1],
				  color      => $PrettyPlots::colorPairs{'redYellowFaint'}
       );
    &PrettyPlots::Prepare_Dataset(\$plot,
				  +$model->{'massStellar'           },
				  +(
				      +     $model->{'massFunctionCumulative'        }
				      -sqrt($model->{'massFunctionCumulativeVariance'})
				  )
				  *$model->{'completeness'          },
				  style      => "line",
				  weight     => [3,1],
				  color      => $PrettyPlots::colorPairs{'redYellowFaint'}
       );
    &PrettyPlots::Prepare_Dataset(\$plot,
				  +$model->{'massStellar'           },
				  +(
				      +     $model->{'massFunctionCumulative'}
				      +sqrt($model->{'massFunctionCumulative'})
				  )
				  *$model->{'completeness'          },
				  style      => "line",
				  weight     => [3,1],
				  color      => $PrettyPlots::colorPairs{'lightGoldenrod'}
       );
    &PrettyPlots::Prepare_Dataset(\$plot,
				  +$model->{'massStellar'           },
				  +(
				      +     $model->{'massFunctionCumulative'}
				      -sqrt($model->{'massFunctionCumulative'})
				  )
				  *$model->{'completeness'          },
				  style      => "line",
				  weight     => [3,1],
				  color      => $PrettyPlots::colorPairs{'lightGoldenrod'}
       );
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
}

# Output the results to file if requested.
if ( exists($arguments{'resultFile'}) ) {
    # Interpolate model mass function to observed masses.
    (my $modelMassFunction        , my $modelMassFunctionError        ) = interpolate($masses->{'value'}->($masses->{'order'}),$model->{'massStellar'},$model->{'massFunctionCumulative'               });
    (my $modelMassFunctionVariance, my $modelMassFunctionVarianceError) = interpolate($masses->{'value'}->($masses->{'order'}),$model->{'massStellar'},$model->{'massFunctionCumulativePoissonVariance'});
    my $resultsFile = new PDL::IO::HDF5(">".$arguments{'resultFile'});
    $resultsFile->dataset('x'             )->set($masses->{'value'}->($masses->{'order'}));
    $resultsFile->dataset('y'             )->set($modelMassFunction);
    $resultsFile->dataset('error'         )->set(pdl zeroes(nelem($number)));
    $resultsFile->dataset('covariance'    )->set(pdl stretcher($modelMassFunctionVariance));
    $resultsFile->dataset('yData'         )->set($number->(-1:0));
    $resultsFile->dataset('covarianceData')->set(pdl zeroes(nelem($number),nelem($number)));
}

# Compute the likelihood:
if ( exists($arguments{'outputFile'}) ) {
    my $constraint;
    my $logLikelihood;
    # Determine errors on masses in log-space.
    my $errorsLogarithmic = $masses->{'error'}/$masses->{'value'}/log(10.0);
    my $massesLogarithmic = log10($masses->{'value'});
    # Interpolate model mass function to observed masses.
    (my $modelMassFunction, my $modelMassFunctionError) = interpolate($masses->{'value'}->($masses->{'order'}),$model->{'massStellar'},$model->{'massFunctionCumulative'}*$model->{'completeness'});
    # Determine change in model mass function between each bin.
    $modelMassFunction->(0:-2) -= $modelMassFunction->(1:-1);   
    # Find zeros of the mass function.
    my $zeroMassFunction = which($modelMassFunction <= 0.0);
    # Set zero mass function to small but non-zero value.
    $modelMassFunction->($zeroMassFunction) .= 1.0e-3;
    # Determine log-likelihood of raw data.
    my $logLikelihoodRaw = sum(log($modelMassFunction)-$modelMassFunction);
    # Monte Carlo sample from the observational mass error distribution.
    my $realizationCount = 10000;
    my $likelihood       = pdl 0.0; 
    for(my $i=0;$i<$realizationCount;++$i) {
	my $massesPerturbed = 10.0**($massesLogarithmic+grandom(nelem($massesLogarithmic))*$errorsLogarithmic);
	my $orderPerturbed  = $massesPerturbed->qsorti();
	# Interpolate model mass function to observed masses.
	(my $modelMassFunction        , my $modelMassFunctionError        ) = interpolate($massesPerturbed->($orderPerturbed),$model->{'massStellar'},$model->{'massFunctionCumulative'               }*$model->{'haloCount'}   );
	(my $modelMassFunctionVariance, my $modelMassFunctionVarianceError) = interpolate($massesPerturbed->($orderPerturbed),$model->{'massStellar'},$model->{'massFunctionCumulativePoissonVariance'}*$model->{'haloCount'}**2);
	(my $modelCompleteness        , my $modelCompletenessError        ) = interpolate($massesPerturbed->($orderPerturbed),$model->{'massStellar'},$model->{'completeness'                         }                         );
	# Determine change in model mass function between each bin.
	$modelMassFunction        ->(0:-2) -= $modelMassFunction        ->(1:-1);
	$modelMassFunctionVariance->(0:-2) -= $modelMassFunctionVariance->(1:-1);
	# Find zeros of the mass function.
	my $zeroMassFunction = which(($modelMassFunction <= 0.0) | ($modelMassFunctionVariance <= 0.0));
	# Set zero mass function to small but non-zero value.
	$modelMassFunction        ->($zeroMassFunction) .= 1.0e-3;
	$modelMassFunctionVariance->($zeroMassFunction) .= 1.0e-3;
	# Sample from the distribution of model mass functions, assuming a Poisson
	# process. Since our mean mass function is a weighted sum of Poisson distributed
	# variables the exact form of the distribution is not simple. We use the approximation
	# proposed by Fay & Feuer (1997; Stat Med; 16; 791;
	# http://www.ncbi.nlm.nih.gov/pubmed/9131766).
	my $k = $modelMassFunction/$modelMassFunctionVariance;
	my $modelMassFunctionUnscaled = $modelMassFunction*$k;
	my $rng = PDL::GSL::RNG->new('taus');
	$rng->set_seed(time());
	$rng->ran_feed_poisson($modelMassFunctionUnscaled);
	my $modelMassFunctionPerturbed = $modelMassFunctionUnscaled*$modelCompleteness/$k/$model->{'haloCount'};
	# Evaluate likelihood, catching impossible cases.
	$likelihood += exp(sum(log($modelMassFunction)-$modelMassFunction)-$logLikelihoodRaw)
	    if (all($modelMassFunction > 0.0));
    }
    $logLikelihood = log($likelihood/$realizationCount)+$logLikelihoodRaw;
    # Include any log-likelihood arising from the requirement to have Magellanic-like halos in the system.
    $constraint->{'logLikelihood'        } = $logLikelihood->sclr()+$model->{'logLikelihood'}->sclr();
    $constraint->{'logLikelihoodVariance'} = 0.0;
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

exit;
