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
use UNIVERSAL qw(isa);

# Generate sets of cosmological parameters drawn at random from the WMAP-9 constraints using the full covariance matrix.
# Andrew Benson (15-September-2010)

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

# Generate Gaussian random numbers.
my $deviates = grandom($parameterCount);

# Generate a set of parameters.
my $parameters = $mean + ($deviates x $choleskyDecomposed);

# Compute required parameters.
my $Omega_M            = $parameters(($parameterMap{'omega_M'}))/($parameters(($parameterMap{'H_0'}))/100.0)**2;
my $Omega_DE           = 1.0-$Omega_M;
my $Omega_b            = $parameters(($parameterMap{'omega_b'}))/($parameters(($parameterMap{'H_0'}))/100.0)**2;
my $sigma_8            = $parameters(($parameterMap{'sigma_8'}));
my $H_0                = $parameters(($parameterMap{'H_0'}));
my $powerSpectrumIndex = $parameters(($parameterMap{'n_s'}));

# Construct data for XML output.
my $parameterData =
{
    cosmologyParameters =>
    {
	value => "simple",
	"Omega_M"            => {value => $Omega_M           ->sclr()},
	"Omega_DE"           => {value => $Omega_DE          ->sclr()},
	"Omega_b"            => {value => $Omega_b           ->sclr()},
	"sigma_8"            => {value => $sigma_8           ->sclr()},
	"H_0"                => {value => $H_0               ->sclr()},
	"powerSpectrumIndex" => {value => $powerSpectrumIndex->sclr()}
    }
};

# Output data as XML.
my $xmlOutput = new XML::Simple (RootName=>"parameters");
print $xmlOutput->XMLout($parameterData);

exit;
