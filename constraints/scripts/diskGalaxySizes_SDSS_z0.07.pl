#!/usr/bin/env perl
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
use PDL::IO::HDF5;
use Astro::Cosmology;
use XML::Simple;
use Data::Dumper;
require Galacticus::Options;
require Galacticus::HDF5;
require Galacticus::Constraints::Covariances;
require Stats::Histograms;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Compute likelihood (and make a plot) for a Galacticus model given the disk galaxy size data from Shen et al. (2003;
# http://adsabs.harvard.edu/abs/2003MNRAS.343..978S).
# Andrew Benson (18-July-2013)

# Get name of input and output files.
die("diskGalaxySizes_SDSS_z0.07.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments;
&Options::Parse_Options(\@ARGV,\%arguments);

# Define constants.
my $megaParsec                    = pdl 3.08567758e+22;
my $sampleRedshift                =     0.07;
# Coefficient in the relation between disk scale length and halo spin times halo virial radius in the Mo, Mao & White (1998) model
# of disk formation.
my $diskRadiusCoefficient         = pdl 1.0/sqrt(2.0);
#  Conversion from exponential disk scale length to Petrosian r50. Computed by Andrew Benson (see "Math/Petrosian Radii.mw").
my $diskScaleLengthToPetrosianR50 = pdl 1.667632104; 
#  Random error (in dex) on galaxy stellar masses in the Shen et al. (2003) sample. Shen et al. (2003) quote 95% confidence
#  interval on masses of +/-40%, which corresponds to a standard deviation of 0.0806 dex assuming a normal distribution in
#  log10(stellar mass).
my $diskMassRandomError           = pdl 0.0806;
#  Random error (in dex) on galaxy Petrosian half-mass radii. This is estimated, very approximately, from the distribution of
#  fractional errors measured by Louis Abramson from the SDSS database.
my $diskRadiusRandomError         = pdl 0.0128;

# Read data.
my $xml                   = new XML::Simple;
my $dataCompilation       = $xml->XMLin("data/observations/galaxySizes/Galaxy_Sizes_By_Mass_SDSS_Shen_2003.xml");

# Create data structure to read the model results.
my $model;
$model->{'file' } = $galacticusFile;
$model->{'store'} = 0;
$model->{'tree' } = "all";
&HDF5::Get_Parameters($model                );

# Construct cosmologies for data and model.
my $cosmologyData = Astro::Cosmology->new(
    omega_matter => $dataCompilation->{'cosmology'}->{'omegaMatter'},
    omega_lambda => $dataCompilation->{'cosmology'}->{'omegaLambda'},
    H0           => $dataCompilation->{'cosmology'}->{'hubble'     }
    );
my $cosmologyGalacticus = Astro::Cosmology->new(
    omega_matter => $model->{'parameters'}->{'Omega_Matter'},
    omega_lambda => $model->{'parameters'}->{'Omega_DE'    },
    H0           => $model->{'parameters'}->{'H_0'         }
    );

# Compute cosmological correction factor for sizes.
my $angularDistanceData             = $cosmologyData      ->angular_diameter_distance($sampleRedshift);
my $angularDistanceGalacticus       = $cosmologyGalacticus->angular_diameter_distance($sampleRedshift);
my $radiusCosmologyCorrectionFactor = $angularDistanceData/$angularDistanceGalacticus;

# Compute cosmological correction factor for masses.
my $luminosityDistanceData          = $cosmologyData      ->luminosity_distance($sampleRedshift);
my $luminosityDistanceGalacticus    = $cosmologyGalacticus->luminosity_distance($sampleRedshift);
my $massCosmologyCorrectionFactor   = ($luminosityDistanceData/$luminosityDistanceGalacticus)**2;

# Count the number of late-type distributions.
my $distributionLateTypeCount = 0;
foreach my $dataSet ( @{$dataCompilation->{'distribution'}} ) {
    ++$distributionLateTypeCount
	if ( $dataSet->{'sersicIndexMaximum'} <= 2.5 );
}
# Iterate over distributions.
$dataCompilation->{'massMinimum'} = pdl zeroes($distributionLateTypeCount);
$dataCompilation->{'massMaximum'} = pdl zeroes($distributionLateTypeCount);
my $j = -1;
foreach my $dataSet ( @{$dataCompilation->{'distribution'}} ) {
    # Check Sersic index constraints - consider only disk-dominated galaxies.
    if ( $dataSet->{'sersicIndexMaximum'} <= 2.5 ) {
	++$j;
	# Extract data radii and distribution plus errors.
	my @properties = ( 'radius', 'radiusFunction', 'radiusFunctionError' );
	my $data;
	$data->{$_} = pdl @{$dataSet->{$_}->{'datum'}}
	   foreach ( @properties );
	# Construct data compilation arrays if necessary.
	unless ( exists($dataCompilation->{'radius'}) ) {
	    $dataCompilation->{'radiusCount'        } = nelem($data->{'radius'});
	    $dataCompilation->{'radius'             } = pdl zeroes(nelem($data->{'radius'}),$distributionLateTypeCount);
	    $dataCompilation->{'radiusFunction'     } = pdl zeroes(nelem($data->{'radius'}),$distributionLateTypeCount);
	    $dataCompilation->{'radiusFunctionError'} = pdl zeroes(nelem($data->{'radius'}),$distributionLateTypeCount);
	    $dataCompilation->{'covariance'         } = pdl zeroes(
		                                                   $distributionLateTypeCount*nelem($data->{'radius'}),
								   $distributionLateTypeCount*nelem($data->{'radius'})
		                                                  );
	}
	# Convert for units.
	$data->{'radius'} *= $dataSet->{'radius'}->{'unitsInSI'}/$megaParsec
	    if ( $dataSet->{'radius'}->{'units'} ne "Mpc" );
	# Convert for scalings.
	foreach (@properties ) {
	    if ( $dataSet->{$_}->{'scaling'} eq "linear" ) {
		# Nothing to do in this case.
	    } else {
		die("diskGalaxySizes_SDSS_z0.07.pl: unknown scaling for ".$_);
	    }
	}
	# Construct data covariance matrix.
	#  Construct temporary matrices used in finding covariance matrix.
	my $radiusCount                                         = nelem($data->{'radius'});
	my $radiusBinWidth                                      = log10($data->{'radius'}->((1))/$data->{'radius'}->((0)));
	my $covarianceUncorrelated                              = pdl zeroes($radiusCount,$radiusCount);
	my $jacobian                                            = pdl zeroes($radiusCount,$radiusCount);
	my $upperLimits                                         = which($data->{'radiusFunctionError'} < 0.0);
	my $dataRadiusFunction                                  =       $data->{'radiusFunction'}->copy();
	$dataRadiusFunction->($upperLimits)                    .= -$dataRadiusFunction->($upperLimits);
	$covarianceUncorrelated->diagonal(0,1)                 .= $data->{'radiusFunctionError'}**2;
	$covarianceUncorrelated->diagonal(0,1)->($upperLimits) .= 0.0;
	for(my $i=0;$i<$radiusCount;++$i) {
	    $jacobian->(($i),  : ) .= -$dataRadiusFunction*$radiusBinWidth;
	    $jacobian->(($i),($i)) += 1.0;
	}
	$data->{'covariance'}                                 = $jacobian x $covarianceUncorrelated x transpose($jacobian);
	$data->{'covariance'}->diagonal(0,1)->($upperLimits) .= $data->{'radiusFunctionError'}->($upperLimits)**2;

	# Store results in the compilation.
	$dataCompilation->{'radius'             }->(:,($j)) .= $data->{'radius'             };
	$dataCompilation->{'radiusFunction'     }->(:,($j)) .= $data->{'radiusFunction'     };
	$dataCompilation->{'radiusFunctionError'}->(:,($j)) .= $data->{'radiusFunctionError'};
	$dataCompilation->{'covariance'}->
	    (
	     $j*$dataCompilation->{'radiusCount'}:($j+1)*$dataCompilation->{'radiusCount'}-1,
	     $j*$dataCompilation->{'radiusCount'}:($j+1)*$dataCompilation->{'radiusCount'}-1
	    )
	    .= $data->{'covariance'};
	# Find minimum and maximum masses for this bin.
	$dataCompilation->{'massMinimum'}->(($j)) .= $dataSet->{'massMinimum'};
	$dataCompilation->{'massMaximum'}->(($j)) .= $dataSet->{'massMaximum'};
    }
}
$dataCompilation->{'massLogarithmic'} = log10(sqrt($dataCompilation->{'massMinimum'}*$dataCompilation->{'massMaximum'}));
my $radiusLogarithmicBins = log10($dataCompilation->{'radius'});

# Check for pre-computed size function.
my $gotModelSizeFunction = 0;
my @rootGroups = $model->{'hdf5File'}->groups();
if ( grep {$_ eq "analysis"} @rootGroups ) {
    my @analysisGroups = $model->{'hdf5File'}->group('analysis')->groups();
    if ( grep {$_ eq "sdssGalaxySizesZ0.07"} @analysisGroups ) {
	$gotModelSizeFunction  = 1;
	$model->{'radiusFunction'} = $model->{'hdf5File'}->group('analysis')->group('sdssGalaxySizesZ0.07')->dataset('sizeFunction'          )->get();
	$model->{'covariance'    } = $model->{'hdf5File'}->group('analysis')->group('sdssGalaxySizesZ0.07')->dataset('sizeFunctionCovariance')->get();
    }
}

# If no pre-computed size function was found, compute one now.
unless ( $gotModelSizeFunction == 1 ) {
    # Read required datasets from model.
    &HDF5::Count_Trees   ($model                );
    &HDF5::Select_Output ($model,$sampleRedshift);
    &HDF5::Get_Dataset
	(
	 $model,
	 [
	  'mergerTreeWeight'  ,
	  'diskMassStellar'   ,
	  'nodeVirialRadius'  ,
	  'spinSpin'          ,
	  'nodeIsOnMainBranch',
	  'basicMass'
	 ]
	);
    # Construct the model radii and masses.
    $model->{'dataSets'}->{'diskRadius'} = 
	$diskScaleLengthToPetrosianR50            *
	$diskRadiusCoefficient                    *
	$model->{'dataSets'}->{'nodeVirialRadius'}*
	$model->{'dataSets'}->{'spinSpin'        };
    $model->{'dataSets'}->{'diskRadiusLogarithmic'} = log10($model->{'dataSets'}->{'diskRadius'     }*$radiusCosmologyCorrectionFactor);
    $model->{'dataSets'}->{'massLogarithmic'      } = log10($model->{'dataSets'}->{'diskMassStellar'}*  $massCosmologyCorrectionFactor);
    # Apply mass systematics.
    my $massSystematicLogM0         = pdl 11.0;
    my $massSystematicParameterRoot = "diskGalaxySizesSDSSZ0.07MassSystematic";
    my $massSystematicOffset = pdl zeroes(nelem($model->{'dataSets'}->{'massLogarithmic'}));
    my $iMass = 0;
    while ( exists($model->{'parameters'}->{$massSystematicParameterRoot.$iMass}) ) {
	$massSystematicOffset += 
	    $model->{'parameters'}->{$massSystematicParameterRoot.$iMass}
	*($model->{'dataSets'}->{'massLogarithmic'}-$massSystematicLogM0)**$iMass;
    }
    $model->{'dataSets'}->{'massLogarithmic'} += $massSystematicOffset;
    # Apply radius systematics.
    my $radiusSystematicLogR0         = pdl 1.0;
    my $radiusSystematicParameterRoot = "diskGalaxySizesSDSSZ0.07RadiusSystematic";
    my $radiusSystematicOffset = pdl zeroes(nelem($model->{'dataSets'}->{'diskRadiusLogarithmic'}));
    my $iRadius = 0;
    while ( exists($model->{'parameters'}->{$radiusSystematicParameterRoot.$iRadius}) ) {
	$radiusSystematicOffset += 
	    $model->{'parameters'}->{$radiusSystematicParameterRoot.$iRadius}
	*($model->{'dataSets'}->{'diskRadiusLogarithmic'}-$radiusSystematicLogR0)**$iRadius;
    }
    $model->{'dataSets'}->{'diskRadiusLogarithmic'} += $radiusSystematicOffset;
    #  Build model histogram.
    my $weight  = $model->{'dataSets'}->{'mergerTreeWeight'};
    my %options = 
	(
	 differential          => "x"                                         ,
	 normalized            => "x"                                         ,
	 normalizeBy           => "weights"                                   ,
	 gaussianSmoothX       => $diskRadiusRandomError*ones(nelem($weight)) ,
	 gaussianSmoothY       => $diskMassRandomError  *ones(nelem($weight)) ,
	 covarianceModel       => "binomial"                                  ,
	 mainBranchStatus      => $model->{'dataSets'}->{'nodeIsOnMainBranch'},
	 nodeMass              => $model->{'dataSets'}->{'basicMass'         },
	 haloMassBinsPerDecade => 10                                          ,
	 haloMassBinsMinimum   => 1.0e8                                       ,
	 haloMassBinsMaximum   => 1.0e16
	);
    ($model->{'radiusFunction'},$model->{'radiusFunctionError'},$model->{'covariance'}) 
	= &Histograms::Histogram2D($radiusLogarithmicBins,$dataCompilation->{'massLogarithmic'},$model->{'dataSets'}->{'diskRadiusLogarithmic'},$model->{'dataSets'}->{'massLogarithmic'},$weight,%options);
}

# Apply any shifts due to model discrepancy.
if ( exists($arguments{'modelDiscrepancies'}) ) {
    # Locate the path which contains discrepancies.
    my $discrepancyPath = $arguments{'modelDiscrepancies'};
    # Scan the path for discrepancy files.
    opendir(discrepDir,$discrepancyPath);
    while ( my $discrepancy = readdir(discrepDir) ) {
	my $discrepancyFileName = $discrepancyPath."/".$discrepancy."/discrepancyDiskGalaxySizes_SDSS_z0.07.hdf5";
	if ( -e $discrepancyFileName ) {
	    my $discrepancyFile = new PDL::IO::HDF5($discrepancyFileName);
	    my @datasets = $discrepancyFile->datasets();
	    foreach my $dataset ( @datasets ) {
		if ( $dataset eq "multiplicative" ) {
		    # Read the multiplicative discrepancy
		    my $multiplier = $discrepancyFile->dataset('multiplicative')->get();
		    # Adjust the model accordingly.

		    my $covarianceMultiplier = pdl zeroes(nelem($multiplier),nelem($multiplier));
		    $covarianceMultiplier .= $discrepancyFile->dataset('multiplicativeCovariance')->get()
			if ( grep {$_ eq "multiplicativeCovariance"} @datasets );
		    $model->{'covariance'} .=
		 	 $model->{'covariance'}                      *outer($multiplier  ,$multiplier  )
			+                       $covarianceMultiplier*outer($model->{'y'},$model->{'y'})
			+$model->{'covariance'}*$covarianceMultiplier;
		    $model->{'radiusFunction'     } *= $multiplier;
		}
		if ( $dataset eq "additive" ) {
		    # Read the additive discrepancy
		    my $addition = $discrepancyFile->dataset('additive')->get();
		    # Adjust the model accordingly.
		    $model->{'radiusFunction'     } += $addition;
		}
		if ( $dataset eq "additiveCovariance" ) {
		    # Read the covariance of the discrepancy.
		    my $covariance = $discrepancyFile->dataset('additiveCovariance')->get();
		    # Adjust the model discrepancy covariance accordingly.
		    $model->{'covariance'} += $covariance;
		}
	    }
	}
    }
}

# Evaluate the model likelihood.
if ( exists($arguments{'outputFile'}) ) {
    # Construct the full covariance matrix, which is the covariance matrix of the observations
    # plus that of the model.
    my $fullCovariance                   = $dataCompilation->{'covariance'}+$model->{'covariance'};
    # Identify upper limits.
    my $upperLimits                      = which($dataCompilation->{'radiusFunction'} < 0.0);
    my $dataRadiusFunction               =       $dataCompilation->{'radiusFunction'}->flat()->copy();
    $dataRadiusFunction->($upperLimits) .= -$dataRadiusFunction->($upperLimits);
    # Where model points are zero, set them equal to a small fraction of the data. Where the data is an upper limit, this will not
    # affect the answers. Where data has a value this will ensure a low likelihood.
    my $modelRadiusFunction              = $model->{'radiusFunction'}->flat()->copy();
    my $modelZero                        = which($modelRadiusFunction <= 0.0);
    $modelRadiusFunction->($modelZero)  .= 1.0e-3*$dataRadiusFunction->($modelZero);
    # Compute the likelihood.
    my $constraint;
    $constraint->{'logLikelihood'} = &Covariances::ComputeLikelihood($dataRadiusFunction,$modelRadiusFunction,$fullCovariance, upperLimits => $upperLimits);
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

# Output the results to file if requested.
if ( exists($arguments{'resultFile'}) ) {
    my $results;
    for(my $j=0;$j<nelem($dataCompilation->{'massLogarithmic'});++$j) {
	push(@{$results->{'x'}},$radiusLogarithmicBins->(:,($j))->list());
	for(my $i=0;$i<$radiusLogarithmicBins->dim(0);++$i) {
	    push(@{$results->{'y'}},$dataCompilation->{'massLogarithmic'}->(($j))->sclr());
	}
    }
    @{$results->{'z'             }} = $model          ->{'radiusFunction'     }->flat()->list();
    @{$results->{'error'         }} = $model          ->{'radiusFunctionError'}->flat()->list();
    @{$results->{'covariance'    }} = $model          ->{'covariance'         }->flat()->list();
    @{$results->{'zData'         }} = $dataCompilation->{'radiusFunction'     }->flat()->list();
    @{$results->{'covarianceData'}} = $dataCompilation->{'covariance'         }->flat()->list();
    my $xmlOut = new XML::Simple (RootName=>"results", NoAttr => 1);;
    # Output the parameters to file.
    open(pHndl,">".$arguments{'resultFile'});
    print pHndl $xmlOut->XMLout($results);
    close pHndl;
}

# Output accuracy to file if requested.
if ( exists($arguments{'accuracyFile'}) ) {
    my $results;
    for(my $j=0;$j<nelem($dataCompilation->{'massLogarithmic'});++$j) {
	push(@{$results->{'x'}},$radiusLogarithmicBins->(:,($j))->list());
	for(my $i=0;$i<$radiusLogarithmicBins->dim(0);++$i) {
	    push(@{$results->{'y'}},$dataCompilation->{'massLogarithmic'}->(($j))->sclr());
	}
    }
    @{$results->{'zModel'    }} = $model          ->{'radiusFunction'     }->flat()->list();
    @{$results->{'zData'     }} = $dataCompilation->{'radiusFunction'     }->flat()->list();
    @{$results->{'errorModel'}} = $model          ->{'radiusFunctionError'}->flat()->list();
    @{$results->{'errorData' }} = $dataCompilation->{'radiusFunctionError'}->flat()->list();
    my $xmlOut = new XML::Simple (RootName=>"accuracy", NoAttr => 1);;
    # Output the parameters to file.
    open(pHndl,">".$arguments{'accuracyFile'});
    print pHndl $xmlOut->XMLout($results);
    close pHndl;
}

# Create a plot of the radius function.
if ( exists($arguments{'plotFile'}) ) {
    my @plotFiles;
    for(my $i=0;$i<$dataCompilation->{'radius'}->dim(1);++$i) {
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	(my $plotFile = $arguments{'plotFile'}) =~ s/\.pdf$/$i.pdf/;
	($plotFileEPS = $plotFile             ) =~ s/\.pdf$/.eps/;
	push(@plotFiles,$plotFile);
	open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$plotFileEPS."'\n";
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
	my $xMinimum = 0.8e3*minimum($dataCompilation->{'radius'}->(:,($i)));
	my $xMaximum = 1.2e3*maximum($dataCompilation->{'radius'}->(:,($i)));
	my $yMinimum = 0.0;
	my $yMaximum = 1.2*maximum($dataCompilation->{'radiusFunction'}->(:,($i))->append($model->{'radiusFunction'}->(:,($i))));
	print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
	print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
	my $plogMassMinimum = sprintf("%5.2f",log10($dataCompilation->{'massMinimum'}->(($i))));
	my $plogMassMaximum = sprintf("%5.2f",log10($dataCompilation->{'massMaximum'}->(($i))));
	print $gnuPlot "set title 'Petrosian half-light radius distribution; late-type; \$".$plogMassMinimum." < \\log_{10}(M_\\star/M_\\odot) < ".$plogMassMaximum." \$'\n";
	print $gnuPlot "set xlabel 'Petrosian half-light radius; \$r_{50}\$ [kpc]'\n";
	print $gnuPlot "set ylabel 'Differential fraction; \${\\rm d}F/{\\rm d}\\log_{10}r_{50}\$ [dex\$^{-1}\$]'\n";
	# Extract errors for the observed data.
	my $errorData                               = 
	    sqrt(
		$dataCompilation->{'covariance'}->(
		    $i*$dataCompilation->{'radiusCount'}:($i+1)*$dataCompilation->{'radiusCount'}-1,
		    $i*$dataCompilation->{'radiusCount'}:($i+1)*$dataCompilation->{'radiusCount'}-1
		)->diagonal(0,1)
	    );
	# Identify upper limits and construct appropriate inputs for plotting functions.
	my $upperLimits                             = which($dataCompilation->{'radiusFunction'}->(:,($i)) < 0.0);
	my $dataRadiusFunction                      =       $dataCompilation->{'radiusFunction'}->(:,($i))->copy();
	$dataRadiusFunction->($upperLimits)        .= -$dataRadiusFunction->($upperLimits);
	my $errorUp                                 = $errorData->copy();
	my $errorDown                               = $errorData->copy();
	$errorUp           ->($upperLimits)        .= 0.0;
	$errorDown         ->($upperLimits)        .= -0.5;
	my $shortArrows                             = which($dataRadiusFunction->($upperLimits) < 0.6);
	$errorDown->($upperLimits)->($shortArrows) .= -$dataRadiusFunction->($upperLimits)->($shortArrows)+0.1;
	&PrettyPlots::Prepare_Dataset(\$plot,
				      1.0e3*$dataCompilation->{'radius'}->(:,($i)),$dataRadiusFunction,
				      errorUp   => $errorUp,
				      errorDown => $errorDown,
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
				      title     => $dataCompilation->{'label'}
	    );
	my $errorModel = 
	    sqrt(
		$model->{'covariance'}->(
		    $i*$dataCompilation->{'radiusCount'}:($i+1)*$dataCompilation->{'radiusCount'}-1,
		    $i*$dataCompilation->{'radiusCount'}:($i+1)*$dataCompilation->{'radiusCount'}-1
		)
		->diagonal(0,1)
	    );
	&PrettyPlots::Prepare_Dataset(\$plot,
				      1.0e3*$dataCompilation->{'radius'}->(:,($i)),$model->{'radiusFunction'}->(:,($i)),
				      errorUp   => $errorModel,
				      errorDown => $errorModel,
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


exit;
