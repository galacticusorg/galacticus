# Calculate column density of hydrogen in cm^{-2} to the center of each galaxy.

package ColumnDensity;
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use PDL::Constants qw(PI);
use PDL::GSLSF::PSI;
use PDL::Math;
require Galacticus::HDF5;
require Galacticus::Inclination;

%HDF5::galacticusFunctions = 
    (
     %HDF5::galacticusFunctions,
     "columnDensity(Disk|Spheroid)??\$" => \&ColumnDensity::Get_Column_Density
    );

sub Get_Column_Density {
    my $model       = shift;
    my $dataSetName = $_[0];
    # Define constants.
    my $massHydrogen             = pdl 1.67262158000e-27; # kg
    my $massSolar                = pdl 1.98892000000e+30; # kg
    my $megaParsec               = pdl 3.08568024000e+22; # m
    my $hecto                    = pdl 1.00000000000e+02;
    my $hydrogenByMassPrimordial = pdl 0.76000000000e+00;

    # Determine which components are needed.
    my $includeDisk     = 0;
    my $includeSpheroid = 0;
    $includeDisk     = 1
	if ( $dataSetName eq "columnDensity" || $dataSetName eq "columnDensityDisk"     );
    $includeSpheroid = 1
	if ( $dataSetName eq "columnDensity" || $dataSetName eq "columnDensitySpheroid" );

    # Ensure that we have all of the datasets that we need.
    my @requiredDatasets = ( 'nodeIndex' );
    push(@requiredDatasets,'diskRadius','diskMassGas','inclination')
	if ( $includeDisk     == 1 );
    push(@requiredDatasets,'spheroidMassGas', 'spheroidRadius')
	if ( $includeSpheroid == 1 );
    &HDF5::Get_Dataset($model,\@requiredDatasets);
    my $dataSets = $model->{'dataSets'};
    # Compute spheroid column density if necessary.
    my $sigmaSpheroid = pdl zeroes(nelem($dataSets->{'nodeIndex'}));
    if ( $includeSpheroid == 1 ) {
	# Compute central density of spheroid component.
	my $spheroidMassGas        = $dataSets->{'spheroidMassGas'    };
	my $spheroidRadius    = $dataSets->{'spheroidRadius'};
	my $spheroidDensityCentral = $spheroidMassGas/(2.0*PI*$spheroidRadius**3);
	# Using 0.1*(scale length) for inner spheroid cutoff
	my $spheroidRadiusMinimum              = 0.1*$spheroidRadius;
	my $spheroidRadiusMinimumDimensionless = $spheroidRadiusMinimum/$spheroidRadius;
	# Compute column density to center of spheroid.
	$sigmaSpheroid                        .=
	    0.5
	    *$spheroidDensityCentral
	    *$spheroidRadius
	    *(
		-$spheroidRadius
		*(3.0+2.0*$spheroidRadiusMinimumDimensionless)
		/(2.0*(1.0+$spheroidRadiusMinimumDimensionless)**2)
		+2.0*log(1.0+1.0/$spheroidRadiusMinimumDimensionless)
	    );
	$sigmaSpheroid->where($spheroidRadius == 0.0) .= 0.0;
    }
    # Compute disk column density if necessary.
    my $sigmaDisk = pdl zeroes(nelem($dataSets->{'nodeIndex'}));
    if ( $includeDisk == 1 ) {
	# Specify disk scale height in units of disk scale length.
	my $diskHeightRatio    = pdl 0.1;
	# Compute central density of disk component.
	my $diskMassGas        = $dataSets->{'diskMassGas'    };
	my $diskRadius         = $dataSets->{'diskRadius'};
	my $diskDensityCentral = $diskMassGas/(4.0*PI*$diskRadius**3*$diskHeightRatio);
	# Compute column density to center of disk.
	my $diskInclination    = ($dataSets->{'inclination'})*(PI/180.0);
	my $tangentInclination = abs(tan($diskInclination));
	my $inclinationHeight  = $tangentInclination*$diskHeightRatio;
	my $digamma1           = (gsl_sf_psi(   -$inclinationHeight/4.0))[0];
	my $digamma2           = (gsl_sf_psi(0.5-$inclinationHeight/4.0))[0];	
	$sigmaDisk            .= 
	    (1.0/2.0)
	    *$diskDensityCentral
	    *$diskRadius
	    *sqrt(1.0+1.0/$tangentInclination**2)
	    *$inclinationHeight
	    *($inclinationHeight*($digamma1-$digamma2)-2.0);
	$sigmaDisk->where($diskRadius == 0.0) .= 0.0;
    }
    # Evaluate conversion factor from mass column density to hydrogen column density.
    my $hydrogenFactor = $hydrogenByMassPrimordial*$massSolar/($massHydrogen*($megaParsec*$hecto)**2);
    # Compute the column density.
    $dataSets->{$dataSetName} = ($sigmaDisk+$sigmaSpheroid)*$hydrogenFactor;

}

1;

