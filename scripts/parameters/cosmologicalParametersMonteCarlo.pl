#!/usr/bin/env perl
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::MatrixOps;
use UNIVERSAL 'isa';
use Data::Dumper;

# Generate sets of cosmological parameters drawn at random from the WMAP-7 constraints using the full covariance matrix.
# Andrew Benson (15-September-2010)

# Read the parameters and their covariances.
$xml = new XML::Simple;
$data = $xml->XMLin($galacticusPath."data/Cosmological_Parameters_WMAP-7.xml");
$parameterCount = 0;
foreach $parameter ( @{$data->{'parameter'}} ) {
    $parameterMap{$parameter->{'label'}} = $parameterCount;
    ++$parameterCount;
}

# Create the covariance matrix and means vector.
$mean       = pdl zeroes($parameterCount);
$covariance = pdl zeroes($parameterCount,$parameterCount);
foreach $parameterA ( @{$data->{'parameter'}} ) {
    $indexA = $parameterMap{$parameterA->{'label'}};
    $mean(($indexA)) .= $parameterA->{'mean'};
    if (isa($parameterA->{'parameter'},'ARRAY')) {
	@parametersB = @{$parameterA->{'parameter'}};
    } else {
	@parametersB = ( $parameterA->{'parameter'} );
    }
    foreach $parameterB ( @parametersB ) {
	$indexB = $parameterMap{$parameterB->{'label'}};
	$covariance(($indexA),($indexB)) .= $parameterB->{'covariance'};
	$covariance(($indexB),($indexA)) .= $parameterB->{'covariance'};
    }
}

# Perform a Cholesky decomposition on the covariance matrix.
$choleskyDecomposed = mchol($covariance);

# Generate Gaussian random numbers.
$deviates = grandom($parameterCount);

# Generate a set of parameters.
$parameters = $mean + ($deviates x $choleskyDecomposed);

# Compute required parameters.
$Omega_M            = $parameters(($parameterMap{'omega_M'}))/($parameters(($parameterMap{'H_0'}))/100.0)**2;
$Omega_DE           = 1.0-$Omega_M;
$Omega_b            = $parameters(($parameterMap{'omega_b'}))/($parameters(($parameterMap{'H_0'}))/100.0)**2;
$sigma_8            = $parameters(($parameterMap{'sigma_8'}));
$H_0                = $parameters(($parameterMap{'H_0'}));
$powerSpectrumIndex = $parameters(($parameterMap{'n_s'}));

# Construct data for XML output.
${$parameterData->{'parameter'}}[++$#{$parameterData->{'parameter'}}] = { name => "Omega_M"           , value => $Omega_M           ->list};
${$parameterData->{'parameter'}}[++$#{$parameterData->{'parameter'}}] = { name => "Omega_DE"          , value => $Omega_DE          ->list};
${$parameterData->{'parameter'}}[++$#{$parameterData->{'parameter'}}] = { name => "Omega_b"           , value => $Omega_b           ->list};
${$parameterData->{'parameter'}}[++$#{$parameterData->{'parameter'}}] = { name => "sigma_8"           , value => $sigma_8           ->list};
${$parameterData->{'parameter'}}[++$#{$parameterData->{'parameter'}}] = { name => "H_0"               , value => $H_0               ->list};
${$parameterData->{'parameter'}}[++$#{$parameterData->{'parameter'}}] = { name => "powerSpectrumIndex", value => $powerSpectrumIndex->list};

# Output data as XML.
$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
print $xmlOutput->XMLout($parameterData);

exit;
