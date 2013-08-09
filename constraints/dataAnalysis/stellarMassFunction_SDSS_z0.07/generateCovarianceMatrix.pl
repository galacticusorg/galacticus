#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
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

# Generate a file containing the SDSS stellar mass function from Li & White (2009) along with its covariance matrix.
# Andrew Benson (29-June-2012)

# Global variables.
our $stellarMass;

# Get the parameter file controlling this calculation.
die("Usage: generateCovarianceMatrix.pl <parameterFile>")
    unless ( scalar(@ARGV) == 1 );
my $parameterFile = $ARGV[0];

# Read the parameter file for this covariance calculation.
my $xml          = new XML::Simple;
my $parameters   = $xml->XMLin($parameterFile);

# Create an HDF5 file containing the observed mass function.
unlink($parameters->{'parameter'}->{'massFunctionCovarianceOutputFileName'}->{'value'});
&observedMassFunction($parameterFile);

# Compute the covariance matrix.
system("cd ".$galacticusPath."; make Mass_Function_Covariance.exe");
system("Mass_Function_Covariance.exe ".$parameterFile);

# Open the covariance HDF5 file.
my $hdfFile = new PDL::IO::HDF5(">".$parameters->{'parameter'}->{'massFunctionCovarianceOutputFileName'}->{'value'});

# Check that mass bins match up.
my $mass = $hdfFile->dataset('mass')->get();
if ( nelem($mass) == nelem($stellarMass) ) {
    die("generateCovarianceMatrix.pl: masses in data and covariance matrix do not match\n")
	unless ( all(abs($mass-$stellarMass) < 1.0e-3*$stellarMass) );
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
$hdfFile->attrSet(label => "Li & White (2009)");

exit;

sub observedMassFunction {
    # Add the observed mass function to an HDF5 file suitable for use with the covariance matrix calculator.
    my $parameterFile = shift;

    # Read the parameter file for this covariance calculation.
    my $xml          = new XML::Simple;
    my $parameters   = $xml->XMLin($parameterFile);

    # Open the covariance HDF5 file.
    my $hdfFile = new PDL::IO::HDF5(">".$parameters->{'parameter'}->{'massFunctionCovarianceOutputFileName'}->{'value'});

    # Read the XML data file.
    my $observed     = $xml->XMLin($galacticusPath."data/observations/massFunctionsStellar/Stellar_Mass_Function_Li_White_2009.xml");
    $stellarMass     = pdl @{$observed->{'stellarMassFunction'}->{'columns'}->{'stellarMass' }->{'datum'}};
    my $massFunction = pdl @{$observed->{'stellarMassFunction'}->{'columns'}->{'massFunction'}->{'datum'}};

    # Convert to non-logarithmic.
    $stellarMass  .= 10.0**$stellarMass ;
    $massFunction .= 10.0**$massFunction;

    # Convert to "h-free" units.
    my $H_0        = pdl $parameters->{'parameter'}->{'H_0'}->{'value'};
    $stellarMass  *= ($H_0/$observed->{'stellarMassFunction'}->{'columns'}->{'stellarMass' }->{'hubble'})**$observed->{'stellarMassFunction'}->{'columns'}->{'stellarMass' }->{'hubbleExponent'};
    $massFunction *= ($H_0/$observed->{'stellarMassFunction'}->{'columns'}->{'massFunction'}->{'hubble'})**$observed->{'stellarMassFunction'}->{'columns'}->{'massFunction'}->{'hubbleExponent'};

    # Convert from per log10(M) to per log(M).
    $massFunction /= log(10.0);
    
    # Store the mass function to file.
    $hdfFile->dataset("massFunctionObserved")->set($massFunction);
    $hdfFile->dataset("mass"                )->attrSet(hubbleExponent => $observed->{'stellarMass'        }->{'columns'}->{'massFunction'}->{'hubbleExponent'});
    $hdfFile->dataset("massFunctionObserved")->attrSet(hubbleExponent => $observed->{'stellarMassFunction'}->{'columns'}->{'massFunction'}->{'hubbleExponent'});

}

exit;
