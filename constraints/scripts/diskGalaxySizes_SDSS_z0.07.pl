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
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
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
my %arguments =
    (
     quiet => 0
    );
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
my $dataCompilation = new PDL::IO::HDF5("data/observations/galaxySizes/Galaxy_Sizes_By_Mass_SDSS_Shen_2003.hdf5");

# Create data structure to read the model results.
my $model;
$model->{'file' } = $galacticusFile;
$model->{'store'} = 0;
$model->{'tree' } = "all";
&HDF5::Get_Parameters($model                );
# Construct cosmologies for data and model.
my $cosmologyData = Astro::Cosmology->new(
    omega_matter => &assclr($dataCompilation->group('cosmology')->attrGet('Omega_Matter')),
    omega_lambda => &assclr($dataCompilation->group('cosmology')->attrGet('Omega_DE'    )),
    H0           => &assclr($dataCompilation->group('cosmology')->attrGet('H_0'         ))
    );
my $cosmologyGalacticus = Astro::Cosmology->new(
    omega_matter => $model->{'parameters'}->{'cosmologyParametersMethod'}->{'OmegaMatter'    }->{'value'},
    omega_lambda => $model->{'parameters'}->{'cosmologyParametersMethod'}->{'OmegaDarkEnergy'}->{'value'},
    H0           => $model->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant' }->{'value'}
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
foreach my $dataSet ( $dataCompilation->groups() ) {
    if ( $dataSet =~ m/^distribution\d+$/ ) {
	++$distributionLateTypeCount
	    if ( &assclr($dataCompilation->group($dataSet)->attrGet('sersicIndexMaximum')) <= 2.5 );
    }
}
# Iterate over distributions.
my $sizeData;
$sizeData->{'label'      } = &assclr($dataCompilation->group("provenance")->attrGet('label'));
$sizeData->{'massMinimum'} = pdl zeroes($distributionLateTypeCount);
$sizeData->{'massMaximum'} = pdl zeroes($distributionLateTypeCount);
my $j = -1;
foreach my $dataSetName ( $dataCompilation->groups() ) {
    if ( $dataSetName =~ m/^distribution\d+$/ ) {
	my $dataSet = $dataCompilation->group($dataSetName);
	# Check Sersic index constraints - consider only disk-dominated galaxies.
	if ( &assclr($dataSet->attrGet('sersicIndexMaximum')) <= 2.5 ) {
	    ++$j;
	    # Extract data radii and distribution plus errors.
	    my @properties = ( 'radius', 'radiusFunction', 'radiusFunctionError' );
	    my $data;
	    $data->{$_} = $dataSet->dataset($_)->get()
		foreach ( @properties );
	    # Construct data compilation arrays if necessary.
	    unless ( exists($sizeData->{'radius'}) ) {
		$sizeData->{'radiusCount'        } = nelem($data->{'radius'});
		$sizeData->{'radius'             } = pdl zeroes(nelem($data->{'radius'}),$distributionLateTypeCount);
		$sizeData->{'radiusFunction'     } = pdl zeroes(nelem($data->{'radius'}),$distributionLateTypeCount);
		$sizeData->{'radiusFunctionError'} = pdl zeroes(nelem($data->{'radius'}),$distributionLateTypeCount);
		$sizeData->{'covariance'         } = pdl zeroes(
		    $distributionLateTypeCount*nelem($data->{'radius'}),
		    $distributionLateTypeCount*nelem($data->{'radius'})
		    );
	    }
	    # Convert for units.
	    $data->{'radius'} *= &assclr($dataSet->dataset('radius')->attrGet('unitsInSI'))/$megaParsec
		if ( &assclr($dataSet->dataset('radius')->attrGet('units')) ne "Mpc" );
	    # Convert for scalings.
	    foreach ( @properties ) {
		if ( &assclr($dataSet->dataset($_)->attrGet('scaling')) eq "linear" ) {
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
	    $sizeData->{'radius'             }->(:,($j)) .= $data->{'radius'             };
	    $sizeData->{'radiusFunction'     }->(:,($j)) .= $data->{'radiusFunction'     };
	    $sizeData->{'radiusFunctionError'}->(:,($j)) .= $data->{'radiusFunctionError'};
	    $sizeData->{'covariance'}->
		(
		 $j*$sizeData->{'radiusCount'}:($j+1)*$sizeData->{'radiusCount'}-1,
		 $j*$sizeData->{'radiusCount'}:($j+1)*$sizeData->{'radiusCount'}-1
		)
		.= $data->{'covariance'};
	    # Find minimum and maximum masses for this bin.
	    $sizeData->{'massMinimum'}->(($j)) .= &assclr($dataSet->attrGet('massMinimum'));
	    $sizeData->{'massMaximum'}->(($j)) .= &assclr($dataSet->attrGet('massMaximum'));
	}
    }
}
$sizeData->{'massLogarithmic'} = log10(sqrt($sizeData->{'massMinimum'}*$sizeData->{'massMaximum'}));
my $radiusLogarithmicBins = log10($sizeData->{'radius'});

# Check for pre-computed size function.
my $gotModelSizeFunction = 0;
my @rootGroups = $model->{'hdf5File'}->groups();
if ( grep {$_ eq "analysis"} @rootGroups ) {
    my @analysisGroups = $model->{'hdf5File'}->group('analysis')->groups();
    if ( grep {$_ eq "sdssGalaxySizesZ0.07"} @analysisGroups ) {
	$gotModelSizeFunction  = 1;
	$model->{'radiusFunction'     } = $model->{'hdf5File'}->group('analysis')->group('sdssGalaxySizesZ0.07')->dataset('sizeFunction'          )->get();
	$model->{'covariance'         } = $model->{'hdf5File'}->group('analysis')->group('sdssGalaxySizesZ0.07')->dataset('sizeFunctionCovariance')->get();
	$model->{'radiusFunctionError'} = $model->{'covariance'}->diagonal(0,1);
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
    while ( exists($model->{'parameters'}->{$massSystematicParameterRoot.$iMass}->{'value'}) ) {
	$massSystematicOffset += 
	    $model->{'parameters'}->{$massSystematicParameterRoot.$iMass}->{'value'}
	*($model->{'dataSets'}->{'massLogarithmic'}->{'value'}-$massSystematicLogM0)**$iMass;
    }
    $model->{'dataSets'}->{'massLogarithmic'} += $massSystematicOffset;
    # Apply radius systematics.
    my $radiusSystematicLogR0         = pdl 1.0;
    my $radiusSystematicParameterRoot = "diskGalaxySizesSDSSZ0.07RadiusSystematic";
    my $radiusSystematicOffset = pdl zeroes(nelem($model->{'dataSets'}->{'diskRadiusLogarithmic'}));
    my $iRadius = 0;
    while ( exists($model->{'parameters'}->{$radiusSystematicParameterRoot.$iRadius}->{'value'}) ) {
	$radiusSystematicOffset += 
	    $model->{'parameters'}->{$radiusSystematicParameterRoot.$iRadius}->{'value'}
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
	= &Histograms::Histogram2D($radiusLogarithmicBins,$sizeData->{'massLogarithmic'},$model->{'dataSets'}->{'diskRadiusLogarithmic'},$model->{'dataSets'}->{'massLogarithmic'},$weight,%options);
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
		    $model->{'covariance'    } .= $model->{'covariance'}*outer($multiplier,$multiplier);
		    $model->{'radiusFunction'} *= $multiplier;
		}
		if ( $dataset eq "multiplicativeCovariance" ) {
		    # Adjust the model accordingly.
		    my $covarianceMultiplier = $discrepancyFile->dataset('multiplicativeCovariance')->get();
		    $model->{'covariance'} += $covarianceMultiplier*outer($model->{'radiusFunction'},$model->{'radiusFunction'});
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
    my $fullCovariance                   = $sizeData->{'covariance'}+$model->{'covariance'};
    # Identify upper limits.
    my $upperLimits                      = which($sizeData->{'radiusFunction'} < 0.0);
    my $dataRadiusFunction               =       $sizeData->{'radiusFunction'}->flat()->copy();
    $dataRadiusFunction->($upperLimits) .= -$dataRadiusFunction->($upperLimits);
    # Where model points are zero, set them equal to a small fraction of the data. Where the data is an upper limit, this will not
    # affect the answers. Where data has a value this will ensure a low likelihood.
    my $modelRadiusFunction              = $model->{'radiusFunction'}->flat()->copy();
    my $modelZero                        = which($modelRadiusFunction <= 0.0);
    $modelRadiusFunction->($modelZero)  .= 1.0e-3*$dataRadiusFunction->($modelZero);
    # For each function, find the peak value, and exclude from the likelihood calculation.
    my $exclude = pdl ones($sizeData->{'radiusFunction'});
    for(my $i=0;$i<$sizeData->{'radiusFunction'}->dim(1);++$i) {
	# Get the index of the peak value in this size function.
	my $maxIndex = maximum_ind($sizeData->{'radiusFunction'}->(:,($i)));
	# Mark this element as excluded.
	$exclude->(($maxIndex),($i)) .= 0;
    }
    # Create excluded copies.
    my $includeIndices              = which($exclude == 1.0);
    my $dataRadiusFunctionExcluded  = $dataRadiusFunction                     ->           ($includeIndices                )  ;
    my $modelRadiusFunctionExcluded = $modelRadiusFunction                    ->           ($includeIndices                )  ;
    my $fullCovarianceExcluded      = $fullCovariance                         ->           ($includeIndices,$includeIndices)  ;
    my $modelCovarianceExcluded     = $model              ->    {'covariance'}->           ($includeIndices,$includeIndices)  ;
    my $shiftedIndices              = $exclude            ->flat(            )->cumusumover(                               )-1;
    my $upperLimitsExcluded         = $shiftedIndices                         ->           ($upperLimits                   )  ;
    # Compute the likelihood.
    my $constraint;
    my $logDeterminant;
    my $offsets;
    my $inverseCovariance;
    my $jacobian;
    my $logLikelihood =
	&Covariances::ComputeLikelihood
	(
	 $dataRadiusFunctionExcluded                          ,
	 $modelRadiusFunctionExcluded                         ,
	 $fullCovarianceExcluded                              ,
	 upperLimits                  =>  $upperLimitsExcluded,
	 jacobian                     => \$jacobian           ,
	 offsets                      => \$offsets            ,
	 quiet                        =>  $arguments{'quiet'} ,
	 inversionMethod              => "eigendecomposition" ,
	 productMethod                => "linearSolver"
	);
    $constraint->{'label'        } = "diskGalaxySizeZ0.07";
    $constraint->{'logLikelihood'} = $logLikelihood;
    # Compute the variance in the log-likelihood due to errors in the model.
    my $logLikelihoodVariance = transpose($jacobian) x $modelCovarianceExcluded x $jacobian;
    $constraint->{'logLikelihoodVariance'} = $logLikelihoodVariance->sclr();
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

# Output the results to file if requested.
if ( exists($arguments{'resultFile'}) ) {
    # Identify upper limits.
    my $upperLimits                      = which($sizeData->{'radiusFunction'} < 0.0);
    my $dataRadiusFunction               =       $sizeData->{'radiusFunction'}->flat()->copy();
    $dataRadiusFunction->($upperLimits) .= -$dataRadiusFunction->($upperLimits);
    # For each function, find the peak value, and exclude from the points.
    my $exclude = pdl ones($sizeData->{'radiusFunction'});
    for(my $i=0;$i<$sizeData->{'radiusFunction'}->dim(1);++$i) {
	# Get the index of the peak value in this size function.
	my $maxIndex = maximum_ind($sizeData->{'radiusFunction'}->(:,($i)));
	# Mark this element as excluded.
	$exclude->(($maxIndex),($i)) .= 0;
    }
    # Create excluded copies.
    my $includeIndices              = which($exclude == 1.0);
    my $radiusExcluded              = $sizeData          ->{'radius'        }->flat()->($includeIndices                );
    my $dataRadiusFunctionExcluded  = $dataRadiusFunction                            ->($includeIndices                );
    my $modelRadiusFunctionExcluded = $model             ->{'radiusFunction'}->flat()->($includeIndices                );
    my $dataCovarianceExcluded      = $sizeData          ->{'covariance'    }        ->($includeIndices,$includeIndices);
    my $modelCovarianceExcluded     = $model             ->{'covariance'    }        ->($includeIndices,$includeIndices);
    my $errorExcluded               = sqrt($dataCovarianceExcluded->diagonal(0,1)+$modelCovarianceExcluded->diagonal(0,1));
    # Write out the results.
    my $resultsFile = new PDL::IO::HDF5(">".$arguments{'resultFile'});
    $resultsFile->dataset('x'             )->set($radiusExcluded             );
    $resultsFile->dataset('y'             )->set($modelRadiusFunctionExcluded);
    $resultsFile->dataset('error'         )->set($errorExcluded              );
    $resultsFile->dataset('covariance'    )->set($modelCovarianceExcluded    );
    $resultsFile->dataset('yData'         )->set($dataRadiusFunctionExcluded );
    $resultsFile->dataset('covarianceData')->set($dataCovarianceExcluded     );
}

# Create a plot of the radius function.
if ( exists($arguments{'plotFile'}) ) {
    my @plotFiles;
    for(my $i=0;$i<$sizeData->{'radius'}->dim(1);++$i) {
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	(my $plotFile = $arguments{'plotFile'}) =~ s/\.pdf$/_$i.pdf/;
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
	my $xMinimum = 0.8e3*minimum($sizeData->{'radius'}->(:,($i)));
	my $xMaximum = 1.2e3*maximum($sizeData->{'radius'}->(:,($i)));
	my $yMinimum = 0.0;
	my $yMaximum = 1.2*maximum($sizeData->{'radiusFunction'}->(:,($i))->append($model->{'radiusFunction'}->(:,($i))));
	print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
	print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
	my $plogMassMinimum = sprintf("%5.2f",log10($sizeData->{'massMinimum'}->(($i))));
	my $plogMassMaximum = sprintf("%5.2f",log10($sizeData->{'massMaximum'}->(($i))));
	print $gnuPlot "set title offset 0,-0.5 'Petrosian half-light radius distribution; late-type; \$".$plogMassMinimum." < \\log_{10}(M_\\star/M_\\odot) < ".$plogMassMaximum." \$'\n";
	print $gnuPlot "set xlabel 'Petrosian half-light radius; \$r_{50}\$ [kpc]'\n";
	print $gnuPlot "set ylabel 'Differential fraction; \${\\rm d}F/{\\rm d}\\log_{10}r_{50}\$ [dex\$^{-1}\$]'\n";
	# Extract errors for the observed data.
	my $errorData                               = 
	    sqrt(
		$sizeData->{'covariance'}->(
		    $i*$sizeData->{'radiusCount'}:($i+1)*$sizeData->{'radiusCount'}-1,
		    $i*$sizeData->{'radiusCount'}:($i+1)*$sizeData->{'radiusCount'}-1
		)->diagonal(0,1)
	    );
	# Identify upper limits and construct appropriate inputs for plotting functions.
	my $upperLimits                             = which($sizeData->{'radiusFunction'}->(:,($i)) < 0.0);
	my $dataRadiusFunction                      =       $sizeData->{'radiusFunction'}->(:,($i))->copy();
	$dataRadiusFunction->($upperLimits)        .= -$dataRadiusFunction->($upperLimits);
	my $errorUp                                 = $errorData->copy();
	my $errorDown                               = $errorData->copy();
	$errorUp           ->($upperLimits)        .= 0.0;
	$errorDown         ->($upperLimits)        .= -0.5;
	if ( nelem($upperLimits) > 0 ) {
	    my $shortArrows                             = which($dataRadiusFunction->($upperLimits) < 0.6);
	    $errorDown->($upperLimits)->($shortArrows) .= -$dataRadiusFunction->($upperLimits)->($shortArrows)+0.1;
	}
	&PrettyPlots::Prepare_Dataset(\$plot,
				      1.0e3*$sizeData->{'radius'}->(:,($i)),$dataRadiusFunction,
				      errorUp   => $errorUp,
				      errorDown => $errorDown,
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
				      title     => $sizeData->{'label'}
	    );
	my $errorModel = 
	    sqrt(
		$model->{'covariance'}->(
		    $i*$sizeData->{'radiusCount'}:($i+1)*$sizeData->{'radiusCount'}-1,
		    $i*$sizeData->{'radiusCount'}:($i+1)*$sizeData->{'radiusCount'}-1
		)
		->diagonal(0,1)
	    );
	&PrettyPlots::Prepare_Dataset(\$plot,
				      1.0e3*$sizeData->{'radius'}->(:,($i)),$model->{'radiusFunction'}->(:,($i)),
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

sub assclr {
    # Return the first argument as a scalar.
    return $_[0];
}
