#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::MatrixOps;
use Scalar::Util qw(reftype);
use Data::Dumper;

# Generates a configuration for cosmological parameters based on the WMAP-9 covariance matrix and suitable for use with the
# Galacticus MCMC constraints architecture. Parameters are expressed as linear combiantions of six independent normal deviates
# (labelled "cosmology0" through "cosmology5"). The output of this script can be included into the base parameter file for a
# Galacticus MCMC run, and actual cosmological parameters set by referencing those created in this file.
# Andrew Benson (25-May-2012)

# Mapping of WMAP-9 parameter names to more descriptive parameter names.
my %parameterNameMapping = 
    (
     H_0     => "wmap9HubbleConstant",
     sigma_8 => "wmap9Sigma8",
     omega_M => "wmap9OmegaMatter",
     omega_B => "wmap9OmegaBaryon",
     n_s     => "wmap9PowerSpectrumIndex",
     tau     => "wmap9ReionizationSuppressionOpticalDepth"
    );

# Read the parameters and their covariances.
my $xml = new XML::Simple;
my $data = $xml->XMLin($ENV{'GALACTICUS_DATA_PATH'}."/static/cosmology/Cosmological_Parameters_WMAP-9.xml");
my $parameterCount = 0;
my %parameterMap;
foreach my $parameter ( @{$data->{'parameter'}} ) {
    $parameterMap{$parameter->{'label'}} = $parameterCount;
    ++$parameterCount;
}

# Create the covariance matrix and means vector.
my $mean       = pdl zeroes($parameterCount);
my $covariance = pdl zeroes($parameterCount,$parameterCount);
foreach my $parameterA ( @{$data->{'parameter'}} ) {
    my $indexA = $parameterMap{$parameterA->{'label'}};
    $mean(($indexA)) .= $parameterA->{'mean'};
    my @parametersB;
    if (reftype($parameterA->{'parameter'}) eq 'ARRAY') {
	@parametersB = @{$parameterA->{'parameter'}};
    } else {
	@parametersB = ( $parameterA->{'parameter'} );
    }
    foreach my $parameterB ( @parametersB ) {
	my $indexB = $parameterMap{$parameterB->{'label'}};
	$covariance(($indexA),($indexB)) .= $parameterB->{'covariance'};
	$covariance(($indexB),($indexA)) .= $parameterB->{'covariance'};
    }
}

# Perform a Cholesky decomposition on the covariance matrix.
my $choleskyDecomposed = mchol($covariance);

# Construct the config data.
my $parameterConfig;
foreach my $parameter ( keys(%parameterMap) ) {
    # Find the index of this parameter in the vectors.
    my $index = $parameterMap{$parameter};
    # Begin by setting the parameter equal to its mean value.
    my $value = "=".$mean->(($index))->string();
    # Loop over all independent normal deviates.
    for(my $i=0;$i<$parameterCount;++$i) {
	# Extract the coefficient for this random deviate and (if it is non-zero) add it multiplied by the appropriate random
	# deviate.
	my $coefficient = $choleskyDecomposed->(($index),($i));
	$value .= "+[cosmology".$i."]*".$coefficient
	    if ( $coefficient != 0.0 );
    }
    # Extract the parameter name, mapping to the Galacticus name if necessary.
    my $parameterName = $parameter;
    $parameterName = $parameterNameMapping{$parameter}
         if ( exists($parameterNameMapping{$parameter}) );
    # If the parameter is a density parameter then it has been multiplied by h^2 in the WMAP-9 analysis. Undo this.
    $value = "(".$value.")/([wmap9OmegaMatterHubbleConstant]/100.0)**2"
	if ( $parameter =~ m/omega/ );
    # Store the parameter config in our array.
    $parameterConfig->{$parameterName} = {"value" => $value};
}

# Output the configuration.
my $xmlOut = new XML::Simple (RootName=>"parameters");
print $xmlOut->XMLout($parameterConfig);

exit;
