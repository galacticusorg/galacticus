# Contains a Perl module which implements halo model calculations of
# galaxy clustering. Calculations in this module generally follow
# notation and conventions of Cooray & Sheth (2002; Physics Reports;
# 372; 1-129).
#
# Andrew Benson (01-September-2010)

package HaloModel;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::GSL::INTEG;
use PDL::GSL::INTERP;
use PDL::GSLSF::ERF;
use Data::Dumper;

# Define Pi.
my $Pi = pdl 3.141592653589793;

# Interpolator;
my $interp;
my $interpXi;

# Wavenumber ranges.
my $logWaveNumberMinimum;
my $logWaveNumberMaximum;

# Shared variables;
my $sigmaJ;
my $separation;
my $virialRadius;
my $growthRate;

# Computes the power spectrum of a selected subset of galaxies. This
# subroutine should be passed a standard data hash reference as used by
# the Galacticus::HDF5 module which has been initialized for a
# particular file and output. Additionally, a PDL containing the indices
# of galaxies which are to be used in computing the power spectrum
# should be provided as the second argument. Additional options can be
# provided as a Perl options hash. Currently supported options are:
#
# space => {redshift|real} - if set to "redshift" the power spectrum
#                            will be computed in redshift space,
#                            otherwise in real space.
sub Compute_Power_Spectrum {
    # Compute the power spectrum of a selection of galaxies, using the halo model.

    # Get the data hash and indices of selected galaxies.
    my $dataBlock = shift;
    my $selected  = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};

    # Check that some galaxies were selected.
    if (nelem($selected) < 1) {die("Compute_Power_Spectrum(): no galaxies were selected")};

    # Determine what type of power spectrum to compute.
    my $redshiftSpace;
    if ( exists($options{'space'}) ) {
	if ( $options{'space'} eq "redshift" ) {
	    $redshiftSpace = 1;
	} else {
	    $redshiftSpace = 0;
	}
    } else {
	$redshiftSpace = 0;
    }

    # Open the file.
    &HDF5::Open_File($dataBlock);

    # Read the linear power spectrum.
    my $waveNumber          = $dataBlock->{'hdf5File'}->group("haloModel")->dataset("wavenumber"   )->get;
    my $linearPowerSpectrum = $dataBlock->{'hdf5File'}->group("haloModel")->dataset("powerSpectrum")->get;

    # Get the linear growth factor for this output.
    my @growthFactor              = $dataBlock->{'hdf5File'}->group("Outputs/Output".$dataBlock->{'output'})->attrGet("linearGrowthFactor"             );
    my @growthFactorLogDerivative = $dataBlock->{'hdf5File'}->group("Outputs/Output".$dataBlock->{'output'})->attrGet("linearGrowthFactorLogDerivative");

    # Scale the power spectrum by the growth factor.
    $linearPowerSpectrum *= $growthFactor[0]**2;

    # Get galaxy data.
    my @properties = ('mergerTreeIndex','nodeIndex','isolatedHostIndex','mergerTreeWeight','nodeBias');
    if ( $redshiftSpace == 1 ) {push(@properties,'nodeVirialVelocity','nodeVirialRadius','basicMass')};
    &HDF5::Get_Dataset($dataBlock,\@properties);
    my $dataSets = $dataBlock->{'dataSets'};

    # Acquire data on profiles and occupancy.
    my $occupancy;
    my $profiles;
    for(my $i=0;$i<nelem($selected);++$i) {
	my $treeIndex = $dataSets->{'mergerTreeIndex'  }->index($selected->index($i));
	my $hostIndex = $dataSets->{'isolatedHostIndex'}->index($selected->index($i));
	# Read Fourier profiles of all relevant dark matter halos.
	unless ( exists($profiles->{$treeIndex}->{$hostIndex}) ) {
	    $profiles->{$treeIndex}->{$hostIndex} = $dataBlock->{'hdf5File'}->group("haloModel/Output".$dataBlock->{'output'}."/mergerTree".$treeIndex)->dataset("fourierProfile".$hostIndex)->get;
	}
	# Compute occupancy of isolated halos.
	++$occupancy->{$treeIndex}->{$hostIndex};
    }

    # Compute mean galaxy number density.
    my $meanDensity = $dataSets->{'mergerTreeWeight'}->index($selected)->sum;

    # Compute redshift space terms if required.
    my $R1;
    my $R2;
    if ( $redshiftSpace == 1 ) {
	# Compute the Hubble parameter at the selected redshift.
	my $expansionFactor = $dataBlock->{'outputs'}->{'expansionFactor'}->index($dataBlock->{'output'}-1);
	my $hubble = $dataBlock->{'parameters'}->{'H_0'}
	*sqrt(
	    $dataBlock->{'parameters'}->{'Omega_Matter'}/($expansionFactor**3)
	    +$dataBlock->{'parameters'}->{'Omega_DE'}
	    +(1.0-$dataBlock->{'parameters'}->{'Omega_Matter'}-$dataBlock->{'parameters'}->{'Omega_DE'})/($expansionFactor**2)
	    );
	# Compute the growth rate (the quantity often approximated as f(Omega)=Omega^0.6 at z=0).
	$growthRate = $growthFactorLogDerivative[0];

	# Construct arrays of wavenumber and power spectrum for interpolation.
	my $logWaveNumber    = log($waveNumber);
	my $logPowerSpectrum = log($linearPowerSpectrum);
	$logWaveNumberMinimum = $logWaveNumber->index(0);
	$logWaveNumberMaximum = $logWaveNumber->index(nelem($logWaveNumber)-1);
	
	# Initialize the interpolator.
	$interp = PDL::GSL::INTERP->init('cspline',$logWaveNumber,$logPowerSpectrum);
	
	# Loop over all halos.
	for(my $i=0;$i<nelem($selected);++$i) {
	    my $weight    = $dataSets->{'mergerTreeWeight' }->index($selected->index($i));
	    my $nodeIndex = $dataSets->{'nodeIndex'        }->index($selected->index($i));
	    my $treeIndex = $dataSets->{'mergerTreeIndex'  }->index($selected->index($i));
	    my $hostIndex = $dataSets->{'isolatedHostIndex'}->index($selected->index($i));
	    if ( $hostIndex == $nodeIndex ) {
		# Compute virial 1-D velocity dispersion. The normalization factor appearing below is taken from Table 2 of
		# (Bryan & Norman; 1998; 495; 80-99) in which they calibrate this relation against N-body simulations.
		my $sigmaVirialNormalization = 0.85;
		my $sigmaVirial1D = sqrt(0.5*$sigmaVirialNormalization)*$dataSets->{'nodeVirialVelocity'}->index($selected->index($i));

		# Get the comoving virial radius of the halo.
		$virialRadius = $dataSets->{'nodeVirialRadius'}->index($selected->index($i))/$expansionFactor;

		# Compute halo-halo velocity dispersion.
		my $absoluteTolerance   = 1.0e-10;
		my $relativeTolerance   = 1.0e-3;
		my $maximumSubdivisions = 1000;		
		$sigmaJ = -1;
		my $sigmaSquaredError;
		my $iError;
		(my $sigmaSquaredMinus1, $sigmaSquaredError, $iError) = gslinteg_qagiu(\&Sigma_Integrand,0.0
										  ,$absoluteTolerance
										  ,$relativeTolerance
										  ,$maximumSubdivisions);
		$sigmaJ = 0;
		(my $sigmaSquared0     , $sigmaSquaredError, $iError) = gslinteg_qagiu(\&Sigma_Integrand,0.0
										  ,$absoluteTolerance
										  ,$relativeTolerance
										  ,$maximumSubdivisions);
		$sigmaJ = 1;
		(my $sigmaSquaredPlus1 , $sigmaSquaredError, $iError) = gslinteg_qagiu(\&Sigma_Integrand,0.0
										  ,$absoluteTolerance
										  ,$relativeTolerance
										  ,$maximumSubdivisions);

		my $sigmaHalo3D = $hubble*$growthRate*$expansionFactor*sqrt($sigmaSquaredMinus1)
		    *sqrt(1.0-$sigmaSquared0**2/$sigmaSquaredPlus1/$sigmaSquaredMinus1);

		# Combine virial and halo velocity dispersions.
		my $sigma = sqrt($sigmaVirial1D**2+($sigmaHalo3D**2)/3.0);

		# Convert velocity dispersion to a comoving distance.
		$sigma /= $hubble*$expansionFactor;
		# Compute the R-factors for 1-halo term.
		my $alpha = $waveNumber*$sigma*sqrt(1.0/2.0);
		$R1->{$treeIndex}->{$hostIndex} = sqrt($Pi)*erf($alpha)/2.0/$alpha;
		$alpha = $waveNumber*$sigma*sqrt(2.0/2.0);
		$R2->{$treeIndex}->{$hostIndex} = sqrt($Pi)*erf($alpha)/2.0/$alpha;
	    }
	}
    }
    
    # Compute the 2-halo effective bias.
    my $Fg;
    my $Fv;
    my $twoHaloFactor;
    for(my $i=0;$i<nelem($selected);++$i) {
	my $weight    = $dataSets->{'mergerTreeWeight' }->index($selected->index($i));
	my $nodeIndex = $dataSets->{'nodeIndex'        }->index($selected->index($i));
	my $treeIndex = $dataSets->{'mergerTreeIndex'  }->index($selected->index($i));
	my $hostIndex = $dataSets->{'isolatedHostIndex'}->index($selected->index($i));
	my $bias      = $dataSets->{'nodeBias'         }->index($selected->index($i));
	if ( $hostIndex == $nodeIndex ) {
	    if ( $redshiftSpace == 1 ) {
		$Fg += $weight*$bias
		    *$profiles->{$treeIndex}->{$hostIndex}
		*$occupancy->{$treeIndex}->{$hostIndex}
		*$R1->{$treeIndex}->{$hostIndex};
		$Fv += $growthRate*$weight*$bias
		    *$profiles->{$treeIndex}->{$hostIndex}
		*$R1->{$treeIndex}->{$hostIndex};
	    } else {
		$twoHaloFactor += $weight*$bias*$profiles->{$treeIndex}->{$hostIndex}*$occupancy->{$treeIndex}->{$hostIndex};
	    }
	}
    }
    if ( $redshiftSpace == 1 ) {
	$twoHaloFactor = ($Fg**2+(2.0/3.0)*$Fg*$Fv+(1.0/5.0)*$Fv**2)/($meanDensity**2);
    } else {
	$twoHaloFactor = ($twoHaloFactor/$meanDensity)**2;
    }

    # Compute 2-halo power spectrum.
    my $twoHaloPowerSpectrum = $linearPowerSpectrum*$twoHaloFactor;
 
    # Compute the 1-halo
    my $oneHaloPowerSpectrum;
    for(my $i=0;$i<nelem($selected);++$i) {
	my $weight    = $dataSets->{'mergerTreeWeight' }->index($selected->index($i));
	my $nodeIndex = $dataSets->{'nodeIndex'        }->index($selected->index($i));
	my $treeIndex = $dataSets->{'mergerTreeIndex'  }->index($selected->index($i));
	my $hostIndex = $dataSets->{'isolatedHostIndex'}->index($selected->index($i));
	my $profileExponent;
	if ( $hostIndex == $nodeIndex ) {
	    my $haloRp = 1.0;
	    if ( $occupancy->{$treeIndex}->{$hostIndex} > 1 ) {
		$profileExponent = 2;
		if ( $redshiftSpace == 1 ) {$haloRp = $R1->{$treeIndex}->{$hostIndex}};
	    } else {
		$profileExponent = 1;
		if ( $redshiftSpace == 1 ) {$haloRp = $R2->{$treeIndex}->{$hostIndex}};
	    }
	    $oneHaloPowerSpectrum += 
		$weight
		*$haloRp
		*($profiles->{$treeIndex}->{$hostIndex}**$profileExponent)
		*$occupancy->{$treeIndex}->{$hostIndex}*($occupancy->{$treeIndex}->{$hostIndex}-1);
	}
    }
    $oneHaloPowerSpectrum /= $meanDensity**2;

    # Compute the total power spectrum.
    my $powerSpectrum = $twoHaloPowerSpectrum+$oneHaloPowerSpectrum;

    # Return the power spectrum.
    return ($waveNumber,$linearPowerSpectrum,$powerSpectrum);
}

# Internal integrand function used in computing halo-halo velocity
# dispersions. See eqn. (142) of Cooray & Sheth (2002; Physics Reports;
# 372; 1-129).
sub Sigma_Integrand {
    # Get the wavenumber.
    my ($myWaveNumber) = @_;

    # Get logarithm of wavenumber.
    my $logWaveNumber   = log($myWaveNumber);

    # Interpolate to get the power spectrum.
    my $myPowerSpectrum = 0.0;    
    $myPowerSpectrum    = exp($interp->eval($logWaveNumber)) if ( $logWaveNumber > $logWaveNumberMinimum && $logWaveNumber < $logWaveNumberMaximum );

    # Compute the window function.
    my $x = $myWaveNumber*$virialRadius;
    my $windowFunction = (3.0/$x**3)*(sin($x)-$x*cos($x));

    # Compute the integrand.
    my $integrand = ($myWaveNumber**(2.0+2.0*$sigmaJ))*$myPowerSpectrum*$windowFunction**2/2.0/($Pi**2);
    return $integrand;
}

# Computes a two-point correlation function from an input power spectrum
# over a specified range of separations. The input wavenumber and power
# spectrum PDLs (first and second arguments) should be generated by the
# Compute_Power_Spectrum subroutine above. The minimum and maximum
# separations to consider are given as the third and fourth arguments,
# while the fifth arguments specifies the number of points per decade in
# separation at which the correlation function should be tabulated.
sub Compute_Correlation_Function {   
    # Get inputs.
    my $waveNumber                = shift;
    my $powerSpectrum             = shift;
    my $separationMinimum         = shift;
    my $separationMaximum         = shift;
    my $separationPointsPerDecade = shift;

    # Set integration accuracy parameters.
    my $absoluteTolerance = 1.0e-3;
    my $maxDivisions      = 1000;

    # Construct arrays of wavenumber and power spectrum for interpolation.
    my $logWaveNumber    = log($waveNumber);
    my $logPowerSpectrum = log($powerSpectrum);

    # Initialize the interpolator.
    $interpXi = PDL::GSL::INTERP->init('cspline',$logWaveNumber,$logPowerSpectrum);

    # Create a list of separations.
    my $separationCount = int(log($separationMaximum/$separationMinimum)/log(10.0)*$separationPointsPerDecade)+1;
    my $separations     = exp(sequence($separationCount)*log($separationMaximum/$separationMinimum)/($separationCount-1)+log($separationMinimum));

    # Loop through separations and compute the correlation function.
    my $correlationFunction = pdl [];
    foreach $separation ( $separations->list ) {
	(my $xi, my $xiError, my $iError) = gslinteg_qawf(\&Correlation_Function_Integrand,$separation,'sin',$waveNumber->index(0),$absoluteTolerance,$maxDivisions);
	print "Compute_Correlation_Function(): failed to compute correlation function\n" unless ( $iError == 0 );
	$correlationFunction = $correlationFunction->append($xi);
    }

    # Return the correlation function.
    return ($separations,$correlationFunction);
}

# Internal integrand function used in Fourier transforming a power
# spectrum into a correlation function.
sub Correlation_Function_Integrand {
    # Get the wavenumber.
    my ($myWaveNumber) = @_;

    # Get logarithm of wavenumber.
    my $logWaveNumber   = log($myWaveNumber);

    # Interpolate to get the power spectrum.
    my $myPowerSpectrum = exp($interpXi->eval($logWaveNumber,{Extrapolate => 1}));

    # Compute the integrand.
    my $integrand       = ($myWaveNumber**2)*$myPowerSpectrum/2.0/($Pi**2)/$separation/$myWaveNumber;
    return $integrand;
}

1;
