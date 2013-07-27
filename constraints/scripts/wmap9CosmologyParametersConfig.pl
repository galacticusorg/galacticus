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
unshift(@INC, $galacticusPath."perl"); 
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::MatrixOps;
use UNIVERSAL 'isa';
use Data::Dumper;

# Generates a configuration for cosmological parameters based on the WMAP-9 covariance matrix and suitable for use with the
# Galacticus+BIE constraints architecture. Parameters are expressed as linear combiantions of six independent normal deviates
# (labelled "cosmology0" through "cosmology5"). The output of this script can be dropped into the "parameters" section of a
# Galacticus+BIE config file.
# Andrew Benson (25-May-2012)

# Mapping of WMAP-9 parameter names to Galacticus input parameter names.
my %parameterNameMapping = 
    (
     omega_M => "Omega_Matter",
     omega_B => "Omega_b",
     n_s     => "powerSpectrumIndex",
     tau     => "reionizationSuppressionOpticalDepth"
    );

# Read the parameters and their covariances.
my $xml = new XML::Simple;
my $data = $xml->XMLin($galacticusPath."/data/cosmology/Cosmological_Parameters_WMAP-9.xml");
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
    if (isa($parameterA->{'parameter'},'ARRAY')) {
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
    my $value = $mean->(($index))->string();
    # Loop over all independent normal deviates.
    for(my $i=0;$i<$parameterCount;++$i) {
	# Extract the coefficient for this random deviate and (if it is non-zero) add it multiplied by the appropriate random
	# deviate.
	my $coefficient = $choleskyDecomposed->(($index),($i));
	$value .= "+\%cosmology".$i."*".$coefficient
	    if ( $coefficient != 0.0 );
    }
    # Extract the parameter name, mapping to the Galacticus name if necessary.
    my $parameterName = $parameter;
    $parameterName = $parameterNameMapping{$parameter}
         if ( exists($parameterNameMapping{$parameter}) );
    # If the parameter is a density parameter then it has been multiplied by h^2 in the WMAP-9 analysis. Undo this.
    $value = "(".$value.")/(%H_0/100.0)**2"
	if ( $parameter =~ m/omega/ );
    # Store the parameter config in our array.
    push(
	@{$parameterConfig->{'parameter'}},
	{
	    name   => $parameterName,
	    define => $value
	}
	);
}

# Output the configuration.
my $xmlOut = new XML::Simple (NoAttr=>1, RootName=>"parameters");
print $xmlOut->XMLout($parameterConfig);

exit;
