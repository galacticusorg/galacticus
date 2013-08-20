#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::LinearAlgebra;
use Data::Dumper;
use Switch;

# Generate a file containing the mass function along with its covariance matrix.
# Andrew Benson (29-June-2012)

# Global variables.
our $observedMass;

# Get the parameter file controlling this calculation.
die("Usage: generateCovarianceMatrix.pl <parameterFile> <configFile>")
    unless ( scalar(@ARGV) == 2 );
my $parameterFile = $ARGV[0];
my $configFile    = $ARGV[1];

# Read the parameter file for this covariance calculation.
my $xml          = new XML::Simple;
my $parameters   = $xml->XMLin($parameterFile);

# Create an HDF5 file containing the observed mass function.
unlink($parameters->{'parameter'}->{'massFunctionCovarianceOutputFileName'}->{'value'});
&observedMassFunction($parameterFile,$configFile);

# Compute the covariance matrix.
system("cd ".$galacticusPath."; make Mass_Function_Covariance.exe");
system("Mass_Function_Covariance.exe ".$parameterFile);

# Open the covariance HDF5 file.
my $hdfFile = new PDL::IO::HDF5(">".$parameters->{'parameter'}->{'massFunctionCovarianceOutputFileName'}->{'value'});

# Check that mass bins match up.
my $mass = $hdfFile->dataset('mass')->get();
if ( nelem($mass) == nelem($observedMass) ) {
    die("generateCovarianceMatrix.pl: masses in data and covariance matrix do not match\n")
	unless ( all(abs($mass-$observedMass) < 1.0e-3*$observedMass) );
} else {
    die("generateCovarianceMatrix.pl: number of mass bins in data and covariance matrix do not match\n");
}

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
    $logDeterminantCovariance  = log(mposdet($scaledCovariance))+nelem($mass)*log($covariance->((0),(0)));
}
$hdfFile->dataset("inverseCovariance"       )->set($inverseCovariance       );
$hdfFile->dataset("logDeterminantCovariance")->set($logDeterminantCovariance);

# Add a label for the dataset to the file.
$hdfFile->attrSet(label => $parameters->{'sourceLabel'})
    if ( exists($parameters->{'sourceLabel'}) );

exit;

sub observedMassFunction {
    # Add the observed mass function to an HDF5 file suitable for use with the covariance matrix calculator.
    my $parameterFile = shift;
    my $configFile    = shift;

    # Read the parameter and config files for this covariance calculation.
    my $xml          = new XML::Simple;
    my $parameters   = $xml->XMLin($parameterFile);
    my $config       = $xml->XMLin($configFile   );

    # Open the covariance HDF5 file.
    my $hdfFile      = new PDL::IO::HDF5(">".$parameters->{'parameter'}->{'massFunctionCovarianceOutputFileName'}->{'value'});

    # Read the XML data file.
    die("observedMassFunction(): config file must specify observedDataFile")
 	unless ( exists($config->{'observedDataFile'}) );
    my $observed     = $xml->XMLin($config->{'observedDataFile'});
    $observedMass    = pdl @{$observed->{'massFunction'}->{'columns'}->{'mass'        }->{'datum'}};
    my $massFunction = pdl @{$observed->{'massFunction'}->{'columns'}->{'massFunction'}->{'datum'}};
    
    # Convert to linear scaling.
    if      ( $observed->{'massFunction'}->{'columns'}->{'mass'}->{'scaling'} eq "linear" ) {
	# Nothing to do.
    } elsif ( $observed->{'massFunction'}->{'columns'}->{'mass'}->{'scaling'} eq "log10"  ) {
	$observedMass .= 10.0**$observedMass;
    } else {
	die('observedMassFunction(): unrecognized scaling for mass');
    }
    
    if      ( $observed->{'massFunction'}->{'columns'}->{'massFunction'}->{'scaling'} eq "linear" ) {
	# Nothing to do.
    } elsif ( $observed->{'massFunction'}->{'columns'}->{'massFunction'}->{'scaling'} eq "log10"  ) {
	$massFunction .= 10.0**$massFunction;
    } else {
	die('observedMassFunction(): unrecognized scaling for massFunction');
    }

    # Convert to "h-free" units.
    my $H_0        = pdl $parameters->{'parameter'}->{'H_0'}->{'value'};
    $observedMass *= ($H_0/$observed->{'massFunction'}->{'columns'}->{'mass' }->{'hubble'})**$observed->{'massFunction'}->{'columns'}->{'mass' }->{'hubbleExponent'};
    $massFunction *= ($H_0/$observed->{'massFunction'}->{'columns'}->{'massFunction'}->{'hubble'})**$observed->{'massFunction'}->{'columns'}->{'massFunction'}->{'hubbleExponent'};

    # Convert mass function to per log(M), i.e. per Np ("Neper"; http://en.wikipedia.org/wiki/Neper).
    if      ( $observed->{'massFunction'}->{'columns'}->{'massFunction'}->{'units'} eq "Mpc^-3 dex^-1" ) {
	# Convert from per log10(M) to per log(M).
	$massFunction /= log(10.0);
    } elsif ( $observed->{'massFunction'}->{'columns'}->{'massFunction'}->{'units'} eq "Mpc^-3 Np^-1"  ) {
	# Nothing to do.
    } else {
	die("observedMassFunction(): unrecognized units for massFunction");
    }
    
    # Store the mass function to file.
    $hdfFile->dataset("massFunctionObserved")->set($massFunction);
    $hdfFile->dataset("massFunctionObserved")->attrSet(hubbleExponent => $observed->{'massFunction'}->{'columns'}->{'massFunction'}->{'hubbleExponent'});

}

exit;
