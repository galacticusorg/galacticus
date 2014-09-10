#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = cwd()."/";
}
unshift(@INC,$galacticusPath."perl"); 
use WWW::Curl::Easy;
use WWW::Curl::Form;
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::LinearAlgebra;
use PDL::IO::Misc;
use Data::Dumper;
use DateTime;
use LaTeX::Encode;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require List::ExtraUtils;
require Galacticus::Constraints::Covariances;

# Generate a file containing the projected correlation function along with its covariance matrix.
# Andrew Benson (22-July-2014)

# Global variables.
our $observedSeparation;
our $observedProjectedCorrelationFunction;

# Constants.
my $Pi = 3.1415927;

# Get the parameter file controlling this calculation.
die("Usage: generateCovarianceMatrixProjectedCorrelation.pl <parameterFile> <configFile> <stage>")
    unless ( scalar(@ARGV) == 3 );
my $parameterFile = $ARGV[0];
my $configFile    = $ARGV[1];
my $stage         = $ARGV[2];

# Parse the config file.
my $xml        = new XML::Simple;
my $config     = $xml->XMLin($configFile);
my $sourceLabel = "data";
$sourceLabel = latex_encode($config->{'sourceLabel'})
    if ( exists($config->{'sourceLabel'}) );
my $pbsLabel = "covariance";
$pbsLabel = $config->{'pbsLabel'}
   if ( exists($config->{'pbsLabel'}) );
my $massVariable = "M_\\star";
$massVariable = $config->{'massVariable'}
    if ( exists($config->{'massVariable'}) );

# Create an HDF5 file containing the observed correlation function.
(my $massMinimum, my $massMaximum) = &observedCorrelationFunction($parameterFile,$configFile,$stage);

# Read the parameter file for this covariance calculation.
$parameterFile = $galacticusPath.$parameterFile
    unless ( $parameterFile =~ m/^\// );
my $parameters   = $xml->XMLin($parameterFile);

# Get work directory name.
my $workDirectoryName = `dirname $parameterFile`;
chomp($workDirectoryName);
$workDirectoryName .= "/";

# Compute the covariance matrix.
unless ( $stage == 0 ) {
    # Generate a Pinocchio halo catalog.
    my @pinocchioJobs;
    my @pinocchioHaloCatalogFileNames;
    for(my $iPinocchio=0;$iPinocchio<$config->{'pinocchio'}->{'realizationCount'};++$iPinocchio) {
    	(my $command, my $pinocchioHaloCatalogFileName) = &generateHalosPinocchio($parameterFile,$iPinocchio);
    	push(@pinocchioHaloCatalogFileNames,$pinocchioHaloCatalogFileName);
    	# Submit the job.
    	if ( defined($command) ) {
    	    my $pinocchioJob = 
    	    {
    		batchFile  => $workDirectoryName."/pinocchio/generatePinocchio".$iPinocchio.".pbs",
    		queue      => "batch",
    		nodes      => "nodes=16:ppn=12",
    		wallTime   => "2:00:00",
    		outputFile => $workDirectoryName."/pinocchio/generatePinocchio".$iPinocchio.".log",
    		name       => $pbsLabel."Stage".$stage."GeneratePinocchio".$iPinocchio,
    		commands   => $command
    	    };
    	    push(@pinocchioJobs,$pinocchioJob);	    
    	}
    }
    &Submit_To_PBS(\@pinocchioJobs,1);
    # Compute correlations from Pinocchio mocks.
    my @pinocchioMockJobs;
    for(my $iPinocchio=0;$iPinocchio<$config->{'pinocchio'}->{'realizationCount'};++$iPinocchio) {
	# Skip missing realizations.
	if ( -e $pinocchioHaloCatalogFileNames[$iPinocchio] ) {
	    # Select random origin and rotation for the mock.
	    $parameters->{'parameter'}->{'mockCorrelationFunctionOrigin'        }->{'value'} =          rand()     ." ".        rand()." ".rand();
	    $parameters->{'parameter'}->{'mockCorrelationFunctionRotationVector'}->{'value'} = acos(2.0*rand()-1.0)." ".2.0*$Pi*rand()           ;
	    $parameters->{'parameter'}->{'mockCorrelationFunctionRotationAngle' }->{'value'} =                          2.0*$Pi*rand()           ;
	    # Construct the command to compute correlations from the mock catalog.
	    my $command;
	    $command .= "mpirun -np 1 -hostfile \$PBS_NODEFILE Halo_Model_Mock.exe ".$parameterFile." ".$pinocchioHaloCatalogFileNames[$iPinocchio]." ".$workDirectoryName."pinocchioGalaxies".$iPinocchio.".hdf5\n"
		unless ( -e $workDirectoryName."pinocchioGalaxies".$iPinocchio.".hdf5" );
	    # Generate mass-specific parameter files.
	    for(my $i=0;$i<nelem($massMinimum);++$i) {
		$parameters->{'parameter'}->{'randomSeed'                        }->{'value'} = int(1000000*rand());
		$parameters->{'parameter'}->{'gaussianRandomSeed'                }->{'value'} = int(1000000*rand());
		$parameters->{'parameter'}->{'poissonRandomSeed'                 }->{'value'} = int(1000000*rand());
		$parameters->{'parameter'}->{'mockCorrelationFunctionMassMinimum'}->{'value'} = $massMinimum->(($i))->sclr();
		$parameters->{'parameter'}->{'mockCorrelationFunctionMassMaximum'}->{'value'} = $massMaximum->(($i))->sclr();
		(my $massSpecificParameterFileName = $parameterFile) =~ s/\.xml/$iPinocchio\_$i.xml/;
		unless ( -e $massSpecificParameterFileName ) {
		    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
		    open(oHndl,">".$massSpecificParameterFileName);
		    print oHndl $xmlOutput->XMLout($parameters);
		    close(oHndl);
		}
		# Compute the mock correlation function.
		$command .= "mpirun -np 1 -hostfile \$PBS_NODEFILE Mocks_Correlation_Functions.exe ".$massSpecificParameterFileName." ".$workDirectoryName."pinocchioGalaxies".$iPinocchio.".hdf5 ".$workDirectoryName."pinocchioMock".$iPinocchio."_".$i.".hdf5\n"
		    unless ( -e $workDirectoryName."pinocchioMock".$iPinocchio."_".$i.".hdf5" );
	    }
	    my $pinocchioMockJob = 
	    {
		batchFile  => $workDirectoryName."/pinocchioMockCovariance".$iPinocchio.".pbs",
		queue      => "batch",
		nodes      => "nodes=1:ppn=12",
		wallTime   => "20:00:00",
		outputFile => $workDirectoryName."/pinocchioMockCovariance".$iPinocchio.".log",
		name       => $pbsLabel."Stage".$stage."PinocchioMockCovariance".$iPinocchio,
		commands   => $command
	    };
	    push(@pinocchioMockJobs,$pinocchioMockJob);
	}
    }
    &Submit_To_PBS(\@pinocchioMockJobs,20);
    # Read Pinocchio correlations.
    my $pinocchioData;
    my $separationCount;
    for(my $iPinocchio=0;$iPinocchio<$config->{'pinocchio'}->{'realizationCount'};++$iPinocchio) {
	# Iterate over masses.
	my $complete = 1;
	my @correlation;
	my @correlationSurvey;
	for(my $i=0;$i<nelem($massMinimum);++$i) {
	    # Skip missing realizations.
	    my $pinocchioCorrelationFileName = $workDirectoryName."pinocchioMock".$iPinocchio."_".$i.".hdf5";
	    if ( ! -e $pinocchioCorrelationFileName ) {
		$complete = 0;
		last;
	    }
	    # Read the data.
	    my $pinocchioCorrelationFile = new PDL::IO::HDF5($pinocchioCorrelationFileName);
	    push(@correlation      , $pinocchioCorrelationFile->dataset('projectedCorrelation'      )->get());
	    push(@correlationSurvey, $pinocchioCorrelationFile->dataset('projectedCorrelationSurvey')->get());
	    unless ( defined($separationCount) ) {
		$pinocchioData->{'separation'} =$pinocchioCorrelationFile->dataset('separation')->get();
		$separationCount = nelem($pinocchioData->{'separation'});
	    }
	}
	# Accumulate the data.
	if ( $complete == 1 ) {
	    push(@{$pinocchioData->{'correlation'      }->{'realizations'}},\@correlation      );
	    push(@{$pinocchioData->{'correlationSurvey'}->{'realizations'}},\@correlationSurvey);
	}
    }
    # Construct Pinocchio mean correlations.
    foreach my $correlationType ( 'correlation', 'correlationSurvey' ) {
	my $correlationMean = pdl zeroes($separationCount*nelem($massMinimum));
	for(my $i=0;$i<scalar(@{$pinocchioData->{$correlationType}->{'realizations'}});++$i) {
	    for(my $j=0;$j<nelem($massMinimum);++$j) {
		$correlationMean->($j*$separationCount:($j+1)*$separationCount-1) += $pinocchioData->{$correlationType}->{'realizations'}->[$i]->[$j];
	    }
	}
	$correlationMean /= scalar(@{$pinocchioData->{$correlationType}->{'realizations'}});
	$pinocchioData->{$correlationType}->{'mean'} = $correlationMean;
    }
    # Construct Pinocchio covariances.
    foreach my $correlationType ( 'correlation', 'correlationSurvey' ) {
	my $correlationCovariance = pdl zeroes($separationCount*nelem($massMinimum),$separationCount*nelem($massMinimum));
	for(my $i=0;$i<scalar(@{$pinocchioData->{$correlationType}->{'realizations'}});++$i) {
	    my $correlation = pdl [];
	    for(my $j=0;$j<nelem($massMinimum);++$j) {
		$correlation = $correlation->append($pinocchioData->{$correlationType}->{'realizations'}->[$i]->[$j]);
	    }
	    my $delta = $correlation-$pinocchioData->{$correlationType}->{'mean'};
	    my $outerProduct = outer($delta,$delta);
	    push(@{$pinocchioData->{$correlationType}->{'outerProducts'}},$outerProduct);
	    $correlationCovariance += $outerProduct;
	}
	$correlationCovariance /= scalar(@{$pinocchioData->{$correlationType}->{'realizations'}})-1;
	$pinocchioData->{$correlationType}->{'covariance'} = $correlationCovariance;
    }
    # Generate an N-body halo catalog.
    (my $haloCatalogFileName, my $cosmology) = &generateHalosNBody($parameterFile);
    $parameters->{'parameter'}->{'H_0'                                  }->{'value'} = $cosmology->{'HubbleConstant'                       };
    $parameters->{'parameter'}->{'Omega_DE'                             }->{'value'} = $cosmology->{'OmegaDarkEnergy'                      };
    $parameters->{'parameter'}->{'Omega_Matter'                         }->{'value'} = $cosmology->{'OmegaMatter'                          };
    $parameters->{'parameter'}->{'Omega_b'                              }->{'value'} = $cosmology->{'OmegaBaryon'                          };
    $parameters->{'parameter'}->{'powerSpectrumIndex'                   }->{'value'} = $cosmology->{'spectralIndex'                        };
    $parameters->{'parameter'}->{'sigma_8'                              }->{'value'} = $cosmology->{'sigma8'                               };
    $parameters->{'parameter'}->{'transferFunctionMethod'               }->{'value'} = $cosmology->{'transferFunctionMethod'               };
    $parameters->{'parameter'}->{'virialDensityContrastMethod'          }->{'value'} = $cosmology->{'virialDensityContrastMethod'          };
    $parameters->{'parameter'}->{'virialDensityContrastFoFLinkingLength'}->{'value'} = $cosmology->{'virialDensityContrastFoFLinkingLength'};
    $parameters->{'parameter'}->{'virialDensityContrastFoFDensityRatio' }->{'value'} = $cosmology->{'virialDensityContrastFoFDensityRatio' };
    # Generate N-body galaxy catalogs and compute correlations.
    my @nBodyJobs;
    # Generate galaxy mock catalog.
    for(my $iNBody=0;$iNBody<$config->{'nBody'}->{'realizationCount'};++$iNBody) {
	# Set random seeds for this realization.
	$parameters->{'parameter'}->{'randomSeed'                           }->{'value'} = int(1000000*rand());
	$parameters->{'parameter'}->{'gaussianRandomSeed'                   }->{'value'} = int(1000000*rand());
	$parameters->{'parameter'}->{'poissonRandomSeed'                    }->{'value'} = int(1000000*rand());
	# Select random origin and rotation for the mock.
	$parameters->{'parameter'}->{'mockCorrelationFunctionOrigin'        }->{'value'} =          rand()     ." ".        rand()." ".rand();
	$parameters->{'parameter'}->{'mockCorrelationFunctionRotationVector'}->{'value'} = acos(2.0*rand()-1.0)." ".2.0*$Pi*rand()           ;
	$parameters->{'parameter'}->{'mockCorrelationFunctionRotationAngle' }->{'value'} =                          2.0*$Pi*rand()           ;
	# Output parameter file.
	my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
	my $cosmologySpecificParameterFileName = $workDirectoryName."nBodyMock".$iNBody.".xml";
	open(oHndl,">".$cosmologySpecificParameterFileName);
	print oHndl $xmlOutput->XMLout($parameters);
	close(oHndl);
	# Construct the command to compute correlations from the mock catalog.
	my $command;
	$command .= "mpirun -np 1 -hostfile \$PBS_NODEFILE Halo_Model_Mock.exe ".$cosmologySpecificParameterFileName." ".$haloCatalogFileName." ".$workDirectoryName."nBodyGalaxies".$iNBody.".hdf5\n"
	    unless ( -e $workDirectoryName."nBodyGalaxies".$iNBody.".hdf5" );
	for(my $i=0;$i<nelem($massMinimum);++$i) {
	    # Generate a mass-specific parameter file.
	    $parameters->{'parameter'}->{'randomSeed'                        }->{'value'} = int(1000000*rand());
	    $parameters->{'parameter'}->{'gaussianRandomSeed'                }->{'value'} = int(1000000*rand());
	    $parameters->{'parameter'}->{'poissonRandomSeed'                 }->{'value'} = int(1000000*rand());
	    $parameters->{'parameter'}->{'mockCorrelationFunctionMassMinimum'}->{'value'} = $massMinimum->(($i))->sclr();
	    $parameters->{'parameter'}->{'mockCorrelationFunctionMassMaximum'}->{'value'} = $massMaximum->(($i))->sclr();
	    my $massSpecificParameterFileName = $workDirectoryName."nBodyGalaxies".$iNBody."_".$i.".xml";
	    unless ( -e $massSpecificParameterFileName ) {
		my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
		open(oHndl,">".$massSpecificParameterFileName);
		print oHndl $xmlOutput->XMLout($parameters);
		close(oHndl);
	    }
	    # Compute the mock correlation function.
	    $command .= "mpirun -np 1 -hostfile \$PBS_NODEFILE Mocks_Correlation_Functions.exe ".$massSpecificParameterFileName." ".$workDirectoryName."nBodyGalaxies".$iNBody.".hdf5 ".$workDirectoryName."nBodyMock".$iNBody."_".$i.".hdf5\n"
		unless ( -e $workDirectoryName."nBodyMock".$iNBody."_".$i.".hdf5" );
	}
	# Construct the job.
	my $generateJob = 
	{
	    batchFile  => $workDirectoryName."/nBodyCovariance".$iNBody.".pbs",
	    queue      => "batch",
	    nodes      => "nodes=1:ppn=12",
	    wallTime   => "20:00:00",
	    outputFile => $workDirectoryName."/nBodyCovariance".$iNBody.".log",
	    name       => $pbsLabel."Stage".$stage."NBodyCovariance".$iNBody,
	    commands   => $command
	};
	push(@nBodyJobs,$generateJob);
    }
    &Submit_To_PBS(\@nBodyJobs,20);
    # Read NBody correlations.
    my $nbodyData;
    for(my $iNBody=0;$iNBody<$config->{'nBody'}->{'realizationCount'};++$iNBody) {
	# Iterate over masses.
	my $complete = 1;
	my @correlation;
	my @correlationSurvey;
	for(my $i=0;$i<nelem($massMinimum);++$i) {
	    # Skip missing realizations.
	    my $nbodyCorrelationFileName = $workDirectoryName."nBodyMock".$iNBody."_".$i.".hdf5";
	    if ( ! -e $nbodyCorrelationFileName ) {
		$complete = 0;
		last;
	    }
	    # Read the data.
	    my $nbodyCorrelationFile = new PDL::IO::HDF5($nbodyCorrelationFileName);
	    push(@correlation      , $nbodyCorrelationFile->dataset('projectedCorrelation'      )->get());
	    push(@correlationSurvey, $nbodyCorrelationFile->dataset('projectedCorrelationSurvey')->get());
	}
	# Accumulate the data.
	if ( $complete == 1 ) {
	    push(@{$nbodyData->{'correlation'      }->{'realizations'}},\@correlation      );
	    push(@{$nbodyData->{'correlationSurvey'}->{'realizations'}},\@correlationSurvey);
	}
    }
    # Construct NBody mean correlations.
    foreach my $correlationType ( 'correlation', 'correlationSurvey' ) {
	my $correlationMean = pdl zeroes($separationCount*nelem($massMinimum));
	for(my $i=0;$i<scalar(@{$nbodyData->{$correlationType}->{'realizations'}});++$i) {
	    for(my $j=0;$j<nelem($massMinimum);++$j) {
		$correlationMean->($j*$separationCount:($j+1)*$separationCount-1) += $nbodyData->{$correlationType}->{'realizations'}->[$i]->[$j];
	    }
	}
	$correlationMean /= scalar(@{$nbodyData->{$correlationType}->{'realizations'}});
	$nbodyData->{$correlationType}->{'mean'} = $correlationMean;
    }
    # Construct NBody covariances.
    foreach my $correlationType ( 'correlation', 'correlationSurvey' ) {
	my $correlationCovariance = pdl zeroes($separationCount*nelem($massMinimum),$separationCount*nelem($massMinimum));
	for(my $i=0;$i<scalar(@{$nbodyData->{$correlationType}->{'realizations'}});++$i) {
	    my $correlation = pdl [];
	    for(my $j=0;$j<nelem($massMinimum);++$j) {
		$correlation = $correlation->append($nbodyData->{$correlationType}->{'realizations'}->[$i]->[$j]);
	    }
	    my $delta = $correlation-$nbodyData->{$correlationType}->{'mean'};
	    my $outerProduct = outer($delta,$delta);
	    push(@{$nbodyData->{$correlationType}->{'outerProducts'}},$outerProduct);
	    $correlationCovariance += $outerProduct;
	}
	$correlationCovariance /= scalar(@{$nbodyData->{$correlationType}->{'realizations'}})-1;
	$nbodyData->{$correlationType}->{'covariance'} = $correlationCovariance;
    }
    # Open the covariance HDF5 file.
    my $hdfFile = new PDL::IO::HDF5(">".$parameters->{'parameter'}->{'projectedCorrelationFunctionCovarianceOutputFileName'}->{'value'}); 
    # Get the observed correlation.
    my $correlationSurveyObserved = $hdfFile->dataset('projectedCorrelationFunctionObserved')->get();
    # Apply the "shrinkage" technique of Pope & Szapudi (2008; http://adsabs.harvard.edu/abs/2008MNRAS.389..766P).
    my $covarianceEmpirical;
    my $covariance;
    my $correlation;
    foreach my $correlationType ( 'correlation', 'correlationSurvey' ) {
	# Compute covariance of empirical (i.e. N-body) data.
	my $realizationCount = scalar(@{$nbodyData->{$correlationType}->{'realizations'}});
	my $covarianceEmpirical = pdl zeroes(($separationCount*nelem($massMinimum))**2,($separationCount*nelem($massMinimum))**2);
	for(my $i=0;$i<scalar(@{$nbodyData->{$correlationType}->{'outerProducts'}});++$i) {
	    my $delta = 
		+$nbodyData->{$correlationType}->{'outerProducts'}->[$i]->flat()
		-$nbodyData->{$correlationType}->{'covariance'   }      ->flat()
		*($realizationCount-1)/$realizationCount;
	    $covarianceEmpirical += outer($delta,$delta);
	}
	$covarianceEmpirical *= $realizationCount/($realizationCount-1)**3;
	# Compute the optimal shrinkage intensity.
	my $intensity = sum($covarianceEmpirical)/sum(($pinocchioData->{$correlationType}->{'covariance'}-$nbodyData->{$correlationType}->{'covariance'})**2);
	# Construct optimal covariance estimator.
	$covariance->{$correlationType} =
	    +     $intensity *$pinocchioData->{$correlationType}->{'covariance'}
	    +(1.0-$intensity)*$nbodyData    ->{$correlationType}->{'covariance'};
	# Scale the survey covariance matrix to the observed correlation function.
	if ( $correlationType eq "correlationSurvey" ) {
	    my $scaleFactor = $correlationSurveyObserved->flat()/$pinocchioData->{$correlationType}->{'mean'};
	    $covariance->{$correlationType} *= outer($scaleFactor,$scaleFactor);
	}
	# Compute associated correlation matrix.
	$correlation->{$correlationType} =
	    $covariance->{$correlationType}
	    /outer(
		sqrt($covariance->{$correlationType}->diagonal(0,1)),
		sqrt($covariance->{$correlationType}->diagonal(0,1))
		);
    }
    # Store covariances to file.
    $hdfFile->dataset('covariance'                  )->set(    $covariance       ->{'correlationSurvey'}                                                               );
    $hdfFile->dataset('correlation'                 )->set(    $correlation      ->{'correlationSurvey'}                                                               );
    $hdfFile->dataset('separation'                  )->set(    $pinocchioData                           ->{'separation'}                                               );
    $hdfFile->dataset('projectedCorrelationFunction')->set(    $pinocchioData    ->{'correlationSurvey'}->{'mean'      }->reshape($separationCount,nelem($massMinimum)));  
}

# Open the covariance HDF5 file.
my $hdfFile = new PDL::IO::HDF5(">".$parameters->{'parameter'}->{'projectedCorrelationFunctionCovarianceOutputFileName'}->{'value'});

# Check that separation bins match up.
my $separation = $hdfFile->dataset('separation')->get();
if ( nelem($separation) == nelem($observedSeparation) ) {
    die("generateCovarianceMatrixProjectedCorrelation.pl: masses in data and covariance matrix do not match\n")
	unless ( all(abs($separation-$observedSeparation) < 1.0e-3*$observedSeparation) );
} else {
    die("generateCovarianceMatrix.pl: number of mass bins in data and covariance matrix do not match\n");
}

# Set cosmology scalings for separation and projected correlation function.
$hdfFile->dataset('massMinimum'                 )->attrSet(cosmologyScaling => "luminosity");
$hdfFile->dataset('massMaximum'                 )->attrSet(cosmologyScaling => "luminosity");
$hdfFile->dataset('separation'                  )->attrSet(cosmologyScaling => "angular"   );
$hdfFile->dataset('projectedCorrelationFunction')->attrSet(cosmologyScaling => "none"      );

# Compute the inverse and determinant of the covariance matrix - store to file.
my $covariance             = $hdfFile->dataset('covariance')->get();
my $covarianceZeroDiagonal = $covariance->copy();
$covarianceZeroDiagonal->diagonal(0,1) .= 0.0;
my $inverseCovariance;
my $logDeterminantCovariance;
if ( all($covarianceZeroDiagonal == 0.0) ) {
    $inverseCovariance                 = $covariance->copy();
    $inverseCovariance->diagonal(0,1) .= 1.0/$inverseCovariance->diagonal(0,1);
    $logDeterminantCovariance          = pdl sum(log($covariance->diagonal(0,1)));
} else {
    # Invert the matrix using Cholesky decomposition. Work with a scaled matrix to avoid underflow problems.
    my $scaledCovariance       = $covariance/$covariance->((0),(0));
    $inverseCovariance         = mposinv($scaledCovariance);
    $inverseCovariance        /= $covariance->((0),(0));
    $logDeterminantCovariance  = log(mposdet($scaledCovariance))+nelem($separation)*log($covariance->((0),(0)));
}
$hdfFile->dataset("inverseCovariance"       )->set($inverseCovariance       );
$hdfFile->dataset("logDeterminantCovariance")->set($logDeterminantCovariance);

exit;

sub observedCorrelationFunction {
    # Add the observed projected correlation function to an HDF5 file suitable for use with the covariance matrix calculator.
    my $parameterFile = shift;
    my $configFile    = shift;
    my $stage         = shift;
    # Read the parameter and config files for this covariance calculation.
    my $xml          = new XML::Simple;
    my $parameters   = $xml->XMLin($parameterFile);
    my $config       = $xml->XMLin($configFile   );
    # Remove the old covariance matrix file.
    unlink($parameters->{'parameter'}->{'projectedCorrelationFunctionCovarianceOutputFileName'}->{'value'});
    # Open the covariance HDF5 file.
    my $hdfFile      = new PDL::IO::HDF5(">".$parameters->{'parameter'}->{'projectedCorrelationFunctionCovarianceOutputFileName'}->{'value'});
    # Read the XML data file.
    die("observedCorrelationFunction(): config file must specify observedDataFile")
 	unless ( exists($config->{'observedDataFile'}) );
    my $observed     = $xml->XMLin($config->{'observedDataFile'});
    my @columns      = @{$observed->{'correlationFunction'}->{'columns'}};
    my $i            = -1;
    my $covariance;
    my $correlation;
    my $massMinimum;
    my $massMaximum;
    foreach my $column ( @columns ) {
	++$i;
	unless ( defined($observedProjectedCorrelationFunction) ) {
	    # Get separations.
	    $observedSeparation     = pdl @{$column->{'separation'         }->{'datum'}};	
	    # Convert to linear scaling.
	    if      ( $column->{'separation'}->{'scaling'} eq "linear" ) {
		# Nothing to do.
	    } elsif ( $column->{'separation'}->{'scaling'} eq "log10"  ) {
		$observedSeparation .= 10.0**$observedSeparation;
	    } else {
		die('observedProjectedCorrelationFunction(): unrecognized scaling for separation');
	    }
	    # Convert to "h-free" units.
	    my $H_0              = pdl $parameters->{'parameter'}->{'H_0'}->{'value'};
	    $observedSeparation *= ($H_0/$column->{'separation'}->{'hubble'})**$column->{'separation'}->{'hubbleExponent'};
	    # Initialize correlation function and covariance.
	    $observedProjectedCorrelationFunction = pdl zeroes(nelem($observedSeparation),scalar(@columns));
	    $covariance                           = pdl
		zeroes(
		    nelem($observedSeparation)*scalar(@columns),
		    nelem($observedSeparation)*scalar(@columns)
		)
		if ( $stage == 0 );
	    $correlation         = pdl zeroes(
		nelem($observedSeparation)*scalar(@columns),
		nelem($observedSeparation)*scalar(@columns)
		)
		if ( $stage == 0 );
	    $massMinimum = pdl zeroes(scalar(@columns));
	    $massMaximum = pdl zeroes(scalar(@columns));
	}
	# Get correlation function.
	$observedProjectedCorrelationFunction->(:,($i)) .= pdl @{$column->{'correlationFunction'}->{'datum'}};
	my $errorUp   = pdl @{$column->{'upperError'}->{'datum'}};
	my $errorDown = pdl @{$column->{'lowerError'}->{'datum'}};
	$covariance ->diagonal(0,1)->($i*nelem($observedSeparation):($i+1)*nelem($observedSeparation)-1) .= 0.5*($errorUp+$errorDown)
	    if ( $stage == 0 );
	$correlation->diagonal(0,1)->($i*nelem($observedSeparation):($i+1)*nelem($observedSeparation)-1) .= 1.0
	    if ( $stage == 0 );
	# Convert to linear scaling.
	if      ( $column->{'correlationFunction'}->{'scaling'} eq "linear" ) {
	    # Nothing to do.
	} elsif ( $column->{'correlationFunction'}->{'scaling'} eq "log10"  ) {
	    $observedProjectedCorrelationFunction->(:,($i))  .= 10.0**$observedProjectedCorrelationFunction->(:,($i)) ;
	} else {
	    die('observedProjectedCorrelationFunction(): unrecognized scaling for correlationFunction');
	}
	$massMinimum->(($i)) .= $column->{'mass'}->{'minimum'};
	if ( exists($column->{'mass'}->{'maximum'}) ) {
	    $massMaximum->(($i)) .= $column->{'mass'}->{'maximum'};
	} else {
	    $massMaximum->(($i)) .= 1.0e14;
	}
	# Convert to linear scaling.
	if      ( $column->{'mass'}->{'scaling'} eq "linear" ) {
	    # Nothing to do.
	} elsif ( $column->{'mass'}->{'scaling'} eq "log10"  ) {
	    $massMinimum->(($i)) .= 10.0**$massMinimum->(($i));
	    $massMaximum->(($i)) .= 10.0**$massMaximum->(($i))
		if ( exists($column->{'mass'}->{'maximum'}) );
	} else {
	    die('observedProjectedCorrelationFunction(): unrecognized scaling for mass');
	}
	# Convert to "h-free" units.
	my $H_0              = pdl $parameters->{'parameter'}->{'H_0'}->{'value'};
	$massMinimum->(($i)) *= ($H_0/$column->{'mass'}->{'hubble'})**$column->{'mass'}->{'hubbleExponent'};
	$massMaximum->(($i)) *= ($H_0/$column->{'mass'}->{'hubble'})**$column->{'mass'}->{'hubbleExponent'};
    }

    # Store the correlation function to file.
    $hdfFile->dataset("projectedCorrelationFunctionObserved")->set($observedProjectedCorrelationFunction);
    $hdfFile->dataset("projectedCorrelationFunction"        )->set($observedProjectedCorrelationFunction)
	if ( $stage == 0 );
    $hdfFile->dataset("separationObserved"                  )->set($observedSeparation                  );
    $hdfFile->dataset("separation"                          )->set($observedSeparation                  )
	if ( $stage == 0 );
    $hdfFile->dataset("massMinimum"                         )->set($massMinimum                         );
    $hdfFile->dataset("massMaximum"                         )->set($massMaximum                         );
    $hdfFile->dataset("covariance"                          )->set($covariance                          )
	if ( $stage == 0 );
    $hdfFile->dataset("correlation"                         )->set($correlation                         )
	if ( $stage == 0 );
    # Add a label for the dataset to the file.
    $hdfFile->attrSet(label => $config->{'sourceLabel'})
	if ( exists($config->{'sourceLabel'}) );
    # Set required parameters and output a new parameter file.
    $parameters->{'parameter'}->{'haloModelMockMassMinimum'                }->{'value'} = $massMinimum       ->(( 0))->sclr();
    $parameters->{'parameter'}->{'haloModelMockMassMaximum'                }->{'value'} = $massMaximum       ->((-1))->sclr();
    $parameters->{'parameter'}->{'mockCorrelationFunctionSeparationMinimum'}->{'value'} = $observedSeparation->(( 0))->sclr();
    $parameters->{'parameter'}->{'mockCorrelationFunctionSeparationMaximum'}->{'value'} = $observedSeparation->((-1))->sclr();
    $parameters->{'parameter'}->{'mockCorrelationFunctionSeparationCount'  }->{'value'} = nelem($observedSeparation);
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
    open(oHndl,">".$parameterFile);
    print oHndl $xmlOutput->XMLout($parameters);
    close(oHndl);
    # Return masses used.
    return ($massMinimum, $massMaximum);
}

sub generateHalosPinocchio {
    # Generate a halo catalog from a Pinocchio simulation.
    my $parameterFile = shift;
    my $realization   = shift;
    # Get base directory.
    my $workDirectory = `dirname $parameterFile`;
    chomp($workDirectory);
    $workDirectory .= "/pinocchio/";
    # If the Pinocchio halo file already exists, skip this task.
    my $command;
    unless ( -e $workDirectory."pinocchioHalos".$realization.".hdf5" ) {
	# Read the parameter file for this covariance calculation.
	my $xml          = new XML::Simple;
	my $parameters   = $xml->XMLin($parameterFile);
	# Generate a Pinocchio oututs file.
	my $pinOutFileName = $workDirectory."outputsPinocchio".$realization.".txt";
	open(my $pinOutFile,">".$pinOutFileName);
	print $pinOutFile $parameters->{'parameter'}->{'mockRedshiftMedian' }->{'value'}."\n";
	close($pinOutFile);
	# Generate a Pinocchio parameter file.
	my $littleH0 = $parameters->{'parameter'}->{'H_0'}->{'value'}/100.0;
	my $seed     = int(rand(1000000))+1;
	my $pinParFileName = $workDirectory."parametersPinocchio".$realization.".txt";
	open(my $pinParFile,">".$pinParFileName);
	print $pinParFile "RunFlag            pinMock".$realization."\n";
	print $pinParFile "OutputList         outputsPinocchio".$realization.".txt\n";
	print $pinParFile "Omega0             ".$parameters->{'parameter'}->{'Omega_Matter'}->{'value'}."\n";
	print $pinParFile "OmegaLambda        ".$parameters->{'parameter'}->{'Omega_DE'}->{'value'}."\n";
	print $pinParFile "OmegaBaryon        ".$parameters->{'parameter'}->{'Omega_b'}->{'value'}."\n";
	print $pinParFile "Hubble100          ".$littleH0."\n";
	print $pinParFile "Sigma8             ".$parameters->{'parameter'}->{'sigma_8'}->{'value'}."\n";
	print $pinParFile "PowerSpectrum      ".$parameters->{'parameter'}->{'powerSpectrumIndex'}->{'value'}."\n";
	print $pinParFile "DEw0              -1.0\n";
	print $pinParFile "DEw1               0.0\n";
	print $pinParFile "TabulatedEoSfile   no\n";
	print $pinParFile "BoxSize            ".$parameters->{'parameter'}->{'pinocchioBoxSize'}->{'value'}."\n";
	print $pinParFile "GridSize           ".$parameters->{'parameter'}->{'pinocchioGridSize'}->{'value'}."\n";
	print $pinParFile "RandomSeed         ".$seed."\n";
	print $pinParFile "CatalogInAscii\n";
	print $pinParFile "StartingzForPLC -1\n";
	print $pinParFile "FileWithInputSpectrum           no\n";
	print $pinParFile "InputSpectrum_UnitLength_in_cm  0\n";
	print $pinParFile "WDM_PartMass_in_kev             0.0\n";
	close($pinParFile);
	# Generate Pinnocchio command.
	$command  = "cd ".$workDirectory."\n";
	$command .= "mpirun -np 192 -hostfile \$PBS_NODEFILE pinocchio-3.0.x parametersPinocchio".$realization.".txt\n";
	$command .= "if [ -e pinocchio.pinMock".$realization.".histories.out ]; then\n";
	$command .= "   ".$galacticusPath."constraints/dataAnalysis/scripts/pinocchioToIrate.pl ".$workDirectory." ".$realization." ".$workDirectory."pinocchioHalos".$realization.".hdf5\n";
	$command .= "   if [ -e ".$workDirectory."pinocchioHalos".$realization.".hdf5 ]; then\n";
	$command .= "      rm -f pinocchio.pinMock".$realization.".cosmology.out\n";
	$command .= "      rm -f pinocchio.pinMock".$realization.".histories.out\n";
	$command .= "      rm -f pinocchio.0.0500.pinMock".$realization.".mf.out\n";
	$command .= "      rm -f pinocchio.0.0500.pinMock".$realization.".catalog.out\n";
	$command .= "   else\n";
	$command .= "      echo Pinocchio to IRATE conversion failed\n";
	$command .= "   fi\n";
	$command .= "else\n";
	$command .= "   echo Pinocchio failed\n";
	$command .= "fi\n";
    }
    return $command, $workDirectory."pinocchioHalos".$realization.".hdf5";
}

sub generateHalosNBody {
    # Generate a halo catalog from an N-body simulation.
    my $parameterFile = shift;
    # Get base directory.
    my $workDirectory = `dirname $parameterFile`;
    chomp($workDirectory);
    $workDirectory .= "/nBody/";
    # Define BolshoiP simulation cosmology.
    my $bolshoiCosmology =
    {
	HubbleConstant                        => 70.0,
	OmegaDarkEnergy                       =>  0.69289,
	OmegaMatter                           =>  0.30711,
	OmegaBaryon                           =>  0.048,
	spectralIndex                         =>  0.96,
	sigma8                                =>  0.82,
	transferFunctionMethod                => "CAMB",
	virialDensityContrastMethod           => "friendsOfFriends",
	virialDensityContrastFoFLinkingLength =>  0.17,
	virialDensityContrastFoFDensityRatio  =>  4.688

    };
    # Read the parameter file for this covariance calculation.
    my $xml          = new XML::Simple;
    my $parameters   = $xml->XMLin($parameterFile);
    # Get a list of available redshifts.
    print "Getting Bolshoi-P snapshot redshifts...\n";
    &cosmosimQuery("select zred,snapnum from BolshoiP.Redshifts",$workDirectory."nbodyRedshifts.csv")
	unless ( -e $workDirectory."nbodyRedshifts.csv" );
    # Find redshift most closely matching that required.
    (my $redshift, my $snapshotNumber) = rcols($workDirectory."nbodyRedshifts.csv",1,2,{LINES => "1:", COLSEP => ","});
    my $redshiftTarget     = $parameters->{'parameter'}->{'mockRedshiftMedian' }->{'value'};
    my $haloMassMinimum    = $parameters->{'parameter'}->{'mockHaloMassMinimum'}->{'value'}*($bolshoiCosmology->{'HubbleConstant'}/100.0);
    my $redshiftDifference = abs($redshift-$redshiftTarget);
    my $redshiftSort       = $redshiftDifference->qsorti();
    my $redshiftClosest    = $redshift->($redshiftSort)->(0)->sclr();
    # Get the halos.
    print "Getting Bolshoi-P halos...\n";
    &cosmosimQuery("select x,y,z,vx,vy,vz,mass from BolshoiP.FOF where snapnum=".$snapshotNumber->($redshiftSort)->(0)->sclr()." and mass > ".$haloMassMinimum,$workDirectory."nbodyHalos.csv")
	unless ( -e $workDirectory."nbodyHalos.csv" );
    # Convert to IRATE format.
    unless ( -e $workDirectory."nbodyHalos.hdf5" ) {
	print "Converting Bolshoi-P halos to IRATE format...\n";
	# Read halo data.
	(my $x, my $y, my $z, my $vx, my $vy, my $vz, my $mass) = rcols($workDirectory."nbodyHalos.csv",1,2,3,4,5,6,7,{LINES => "1:", COLSEP => ","});
	# Convert to preferred units.
	$x    /= $bolshoiCosmology->{'HubbleConstant'}/100.0;
	$y    /= $bolshoiCosmology->{'HubbleConstant'}/100.0;
	$z    /= $bolshoiCosmology->{'HubbleConstant'}/100.0;
	$mass /= $bolshoiCosmology->{'HubbleConstant'}/100.0;
	# Construct 3D datasets.
	my $center   = pdl zeroes(3,nelem($x));
	my $velocity = pdl zeroes(3,nelem($x));
	$center  ->((0),:) .= $x;
	$center  ->((1),:) .= $y;
	$center  ->((2),:) .= $z;
	$velocity->((0),:) .= $vx;
	$velocity->((1),:) .= $vy;
	$velocity->((2),:) .= $vz;
	# Output in IRATE format.
	my $outputFile = new PDL::IO::HDF5(">".$workDirectory."nbodyHalos.hdf5");
	# Snapshot group.
	my $snapshot = $outputFile->group('Snapshot00001');
	$snapshot->attrSet(Redshift => pdl $redshiftClosest);
	# Halo catalog.
	my $haloCatalog = $snapshot->group('HaloCatalog');
	$haloCatalog->dataset('Center'  )->set($center  );
	$haloCatalog->dataset('Velocity')->set($velocity);
	$haloCatalog->dataset('Mass'    )->set($mass    );
	$haloCatalog->dataset('Center'  )->attrSet(unitname => "Mpc"   , unitscgs => pdl [3.08568e+24,  0, -1]);
	$haloCatalog->dataset('Velocity')->attrSet(unitname => "km/s"  , unitscgs => pdl [1.00000e+05,  0,  0]);
	$haloCatalog->dataset('Mass'    )->attrSet(unitname => "Msolar", unitscgs => pdl [1.98892e+33,  0,  0]);
	# Cosmology.
	my $cosmology = $outputFile->group('Cosmology');
	$cosmology->attrSet(HubbleParam        => pdl $bolshoiCosmology->{'HubbleConstant '}/100.0);
	$cosmology->attrSet(OmegaBaryon        => pdl $bolshoiCosmology->{'OmegaBaryon'    }      ); 
	$cosmology->attrSet( OmegaLambda       => pdl $bolshoiCosmology->{'OmegaDarkEnergy'}      );   
	$cosmology->attrSet(OmegaMatter        => pdl $bolshoiCosmology->{'OmegaMatter'    }      );    
	$cosmology->attrSet(PowerSpectrumIndex => pdl $bolshoiCosmology->{'spectralIndex'  }      );     
	$cosmology->attrSet(sigma_8            => pdl $bolshoiCosmology->{'sigma8'         }      );
	# Simulation properties.
	my $simulation = $outputFile->group('SimulationProperties');
	$simulation->attrSet(
	    boxSize => pdl 250.0/($bolshoiCosmology->{'HubbleConstant'}/100.0)
	    );
    }
    return ($workDirectory."nbodyHalos.hdf5", $bolshoiCosmology);
}

sub cosmosimQuery {
    # Run a query on the CosmoSim server and download results.
    my $queryString = shift;
    my $fileName    = shift;
    # Parse the Galacticus config file if it is present.
    my $sqlUser;
    my $sqlPassword;
    if ( -e $galacticusPath."/galacticusConfig.xml" ) {
	my $xml    = new XML::Simple();
	my $config = $xml->XMLin($galacticusPath."/galacticusConfig.xml");
	if ( exists($config->{'cosmosimDB'}->{'host'}) ) {
	    my %hosts;
	    if ( exists($config->{'cosmosimDB'}->{'host'}->{'name'}) ) {
		$hosts{'default'} = $config->{'cosmosimDB'}->{'host'};
	    } else {
		%hosts = %{$config->{'cosmosimDB'}->{'host'}};
	    }
	    foreach ( keys(%hosts) ) {
		if ( $_ eq $ENV{'HOSTNAME'} || $_ eq "default" ) {
		    $sqlUser     = $hosts{$_}->{'user'    }
		    if ( exists($hosts{$_}->{'user'    }) );
		    $sqlPassword = $hosts{$_}->{'password'}
		    if ( exists($hosts{$_}->{'password'}) );
		    if ( exists($hosts{$_}->{'passwordFrom'}) ) {
			if ( $hosts{$_}->{'passwordFrom'} eq "input" ) {
			    $sqlPassword = <>;
			    chomp($sqlPassword);
			}
		    }
		}
	    }
	}
    }
    die("generateHalosNBody: CosmoSim database username and password must be defined")
	unless ( defined($sqlUser) && defined($sqlPassword) );
    # Get CUrl objects.
    my $xml      = new XML::Simple    ();
    my $curlPost = new WWW::Curl::Easy();
    my $curlGet  = new WWW::Curl::Easy();
    # Create the query job.
    my $createText;
    my $date       = DateTime->now();
    my $createForm = new WWW::Curl::Form();
    $createForm->formadd("query",$queryString);
    $createForm->formadd("queue","long"      );
    $createForm->formadd("table",$date       );
    $curlPost->setopt(CURLOPT_URL           ,"http://www.cosmosim.org/uws/query");
    $curlPost->setopt(CURLOPT_HTTPPOST      ,$createForm                        );
    $curlPost->setopt(CURLOPT_FOLLOWLOCATION,1                                  );
    $curlPost->setopt(CURLOPT_USERPWD       ,$sqlUser.":".$sqlPassword          );
    $curlPost->setopt(CURLOPT_WRITEDATA     ,\$createText                       );
    print "generateHalosNBody: creating CosmoSim UWS job\n";
    die("generateHalosNBody(): failed to create job")
	unless ( $curlPost->perform() == 0 );
    my $query = $xml->XMLin($createText);
    # Submit the job.
    my $submitText;
    my $submitForm = new WWW::Curl::Form();
    $submitForm->formadd("phase","run");
    $curlPost->setopt(CURLOPT_URL      ,"http://www.cosmosim.org/uws/query/".$query->{'uws:jobId'});
    $curlPost->setopt(CURLOPT_HTTPPOST ,$submitForm                                               );
    $curlPost->setopt(CURLOPT_WRITEDATA,\$submitText                                              );
    print "generateHalosNBody: submitting CosmoSim UWS job\n";
    die("generateHalosNBody(): failed to submit job")
	unless ( $curlPost->perform() == 0 );
    # Check status.
    my $statusText;
    do {
	undef($statusText);
	sleep(10);
	$curlGet->setopt(CURLOPT_URL      ,"http://www.cosmosim.org/uws/query/".$query->{'uws:jobId'}."/phase");
	$curlGet->setopt(CURLOPT_USERPWD  ,$sqlUser.":".$sqlPassword                                          );
	$curlGet->setopt(CURLOPT_WRITEDATA,\$statusText                                                       );
	die("generateHalosNBody(): failed to check status")
	    unless ( $curlGet->perform() == 0 );
	print "generateHalosNBody: CosmoSim UWS job status is '".$statusText."'\n";
	die
	    if ( $statusText eq "ERROR" || $statusText eq "ABORTED" );
    }
    until ( $statusText eq "COMPLETED" );
    # Get results information.
    my $resultsText;
    $curlGet->setopt(CURLOPT_URL      ,"http://www.cosmosim.org/uws/query/".$query->{'uws:jobId'}."/results");
    $curlGet->setopt(CURLOPT_USERPWD  ,$sqlUser.":".$sqlPassword                                            );
    $curlGet->setopt(CURLOPT_WRITEDATA,\$resultsText                                                        );
    die("generateHalosNBody(): failed to get results information")
	unless ( $curlGet->perform() == 0 );
    my $results = $xml->XMLin($resultsText);
    # Download results.
    open(my $resultsFile,">".$fileName);
    $curlGet->setopt(CURLOPT_URL      , $results->{'uws:result'}->{'xlink:href'});
    $curlGet->setopt(CURLOPT_USERPWD  ,$sqlUser.":".$sqlPassword          );
    $curlGet->setopt(CURLOPT_WRITEDATA,$resultsFile                             );
    print "generateHalosNBody: downloading CosmoSim data\n";
    die("generateHalosNBody(): failed to get results information")
	unless ( $curlGet->perform() == 0 );
    close($resultsFile);
    # Delete the job.
    my $deleteText;
    $curlGet->setopt(CURLOPT_URL          , "http://www.cosmosim.org/uws/query/".$query->{'uws:jobId'});
    $curlGet->setopt(CURLOPT_USERPWD      ,$sqlUser.":".$sqlPassword                                  );
    $curlGet->setopt(CURLOPT_CUSTOMREQUEST, "DELETE"                                                  ); 
    $curlGet->setopt(CURLOPT_WRITEDATA    ,\$deleteText                                               );
    print "generateHalosNBody: deleting CosmoSim UWS job\n";
    die("generateHalosNBody(): failed to get results information")
	unless ( $curlGet->perform() == 0 );
}

sub Submit_To_PBS {
    # Submit a job to the PBS queue and wait for it to complete.
    my @pbsStack   = &ExtraUtils::as_array(shift());
    my $jobMaximum =                       shift() ;
    my %pbsJobs;
    # Submit jobs and wait.
    print "Waiting for PBS jobs to finish...\n";
    while ( scalar(keys %pbsJobs) > 0 || scalar(@pbsStack) > 0 ) {
	# Find all PBS jobs that are running.
	my %runningPBSJobs;
	undef(%runningPBSJobs);
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^Job\sId:\s+(\S+)/ ) {$runningPBSJobs{$1} = 1};
	}
	close(pHndl);
	foreach my $jobID ( keys(%pbsJobs) ) {
	    unless ( exists($runningPBSJobs{$jobID}) ) {
		print "PBS job ".$jobID." has finished.\n";
		# Remove the job ID from the list of active PBS jobs.
		delete($pbsJobs{$jobID});
	    }
	}
	# If fewer than maximum number of jobs are in the queue, pop one off the stack.
	if ( scalar(@pbsStack) > 0 && scalar(keys %pbsJobs) < $jobMaximum ) {
	    my $jobDescriptor = pop(@pbsStack);
	    # Assert that job commands must be present.
	    die("Submit_To_PBS: no job commands were supplied")
		unless ( exists($jobDescriptor->{'commands'}) );
	    # Assert that batch file must be present.
	    die("Submit_To_PBS: no batchFile was supplied")
		unless ( exists($jobDescriptor->{'batchFile'}) );
	    # Generate the bacth file.
	    open(oHndl,">".$jobDescriptor->{'batchFile'});
	    print oHndl "#!/bin/bash\n";
	    print oHndl "#PBS -q ".$jobDescriptor->{'queue'}."\n"
		if ( exists($jobDescriptor->{'queue'}) );
	    print oHndl "#PBS -l ".$jobDescriptor->{'nodes'}."\n"
		if ( exists($jobDescriptor->{'nodes'}) );
	    print oHndl "#PBS -j oe\n";
	    print oHndl "#PBS -o ".$jobDescriptor->{'outputFile'}."\n"
		if ( exists($jobDescriptor->{'outputFile'}) );
	    print oHndl "#PBS -l walltime=".$jobDescriptor->{'wallTime'}."\n"
		if ( exists($jobDescriptor->{'wallTime'}) );
	    print oHndl "#PBS -N ".$jobDescriptor->{'name'}."\n"
		if ( exists($jobDescriptor->{'name'}) );
	    print oHndl "#PBS -V\n";
	    print oHndl "cd \$PBS_O_WORKDIR\n";
	    print oHndl "export LD_LIBRARY_PATH=\$HOME/Galacticus/OpenMPI/lib:\$HOME/Galacticus/Tools/lib:\$HOME/Galacticus/Tools/lib64:\$LD_LIBRARY_PATH\n";
	    print oHndl "export PATH=\$HOME/Galacticus/OpenMPI/bin:\$HOME/Galacticus/Tools/bin:\$HOME/perl5/bin:\$PATH\n";
	    print oHndl "export GFORTRAN_ERROR_DUMPCORE=YES\n";
	    print oHndl "export PERL_LOCAL_LIB_ROOT=\"\$HOME/perl5\"\n";
	    print oHndl "export PERL_MB_OPT=\"--install_base \$HOME/perl5\"\n";
	    print oHndl "export PERL_MM_OPT=\"INSTALL_BASE=\$HOME/perl5\"\n";
	    print oHndl "export PERL5LIB=\"\$HOME/perl5/lib/perl5/x86_64-linux-thread-multi:\$HOME/perl5/lib/perl5\"\n";
	    print oHndl "export PYTHONPATH=\"\$HOME/Galacticus/Tools/lib/python:\$HOME/Galacticus/Tools/lib/python2.7:/share/apps/atipa/acms/lib\"\n";
	    print oHndl "ulimit -t unlimited\n";
	    print oHndl "ulimit -c unlimited\n";
	    print oHndl $jobDescriptor->{'commands'};
	    close(oHndl);
	    open(pHndl,"qsub ".$jobDescriptor->{'batchFile'}." |");
	    my $jobID = "";
	    while ( my $line = <pHndl> ) {
		if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
	    }
	    close(pHndl);
	    # Add the job number to the active job hash.
	    unless ( $jobID eq "" ) {
	    	$pbsJobs{$jobID} = 1;
	    }
	    sleep 5;
	} else {
	    # Wait.
	    sleep 60;
	}
    }
}
