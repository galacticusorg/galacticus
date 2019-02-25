#!/usr/bin/env perl
use strict;
use warnings;
use lib './perl';
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::MatrixOps;
use Galacticus::Constraints::Parameters;
use Galacticus::Options;
use Data::Dumper;
use List::ExtraUtils;

# Fit a multivariate normal distribution to MCMC chain samples and
# output a suitable XML definition file to define priors for the
# parameters corresponding to this distribution. This is accomplished
# by definine standard normal distributed meta-parameters and defining
# the actual parameters in terms of these meta-parameters.
# Andrew Benson (07-February-2017)

# Get command line arguments.
die("Usage: fitMultivariateNormal.pl <configFile> [options...]")
    unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %options = 
    (
     sampleCount => 0                 , # Use all chain states by default.
     metaName    => "metaParameter"   ,
     outputFile  => "multiVariate.xml"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Build any parameter name map.
my %parameterNameMap;
if ( exists($options{'parameterName'}) ) {
    foreach ( &List::ExtraUtils::as_array($options{'parameterName'}) ) {
	if ( $_ =~ m/([^\/]+)\/([^\/]+)/ ) {
	    my $inputName   = $1;
	    my $outputName  = $2;
	    $parameterNameMap{$inputName} = $outputName;
	} else {
	    die("fitMultivariateNormal.pl: parameterName option has incorrect syntax");
	}
    }
}

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFile);

# Read a sample matrix of all chain states.
(my $sampleMatrix, my $sampleLikelihood) = &Galacticus::Constraints::Parameters::Sample_Matrix($config,\%options);

# Determine which, if any, parameters are being excluded.
my $j = -1;
my $activeParameters = pdl [];
my $modelParameters = exists($config->{'posteriorSampleSimulationMethod'}->{'modelParameterMethod'}) ? $config->{'posteriorSampleSimulationMethod'}->{'modelParameterMethod'} : $config->{'modelParameterMethod'};
foreach my $parameter ( @{$modelParameters} ) {
    next 
	unless ( exists($parameter->{'distributionFunction1DPrior'}) );
    ++$j;
    $activeParameters = $activeParameters->append($j)
	unless ( exists($options{'exclude'}) && grep {$parameter->{'name'}->{'value'} eq $_} &List::ExtraUtils::as_array($options{'exclude'}) );
}
my $sampleMatrixExcluded = $sampleMatrix->($activeParameters,:);

# Compute the mean parameters.
my $parameterMean = $sampleMatrixExcluded->xchg(0,1)->average();

# Compute the covariance of parameters.
my $parameterCovariance = pdl zeroes(nelem($parameterMean),nelem($parameterMean));
for(my $i=0;$i<nelem($parameterMean);++$i) {
    for(my $j=$i;$j<nelem($parameterMean);++$j) {
	my $covariance = 
	    sum(
		+($sampleMatrixExcluded->(($i),:)-$parameterMean->(($i)))
		*($sampleMatrixExcluded->(($j),:)-$parameterMean->(($j)))
	    )
	    /$sampleMatrixExcluded->dim(1);
	$parameterCovariance->(($i),($j)) .= $covariance;
	$parameterCovariance->(($j),($i)) .= $covariance;
    }
}

# Find Cholesky decomposition of the covariance matrix.
my $cholesky = mchol($parameterCovariance);
# Construct a parameter definition file to describe this multivariate normal distribution.
my $i = -1;
my $parametersOutput;
foreach my $parameter ( @{$modelParameters} ) {
    next 
	unless ( exists($parameter->{'distributionFunction1DPrior'}) );
    next
 	if ( exists($options{'exclude'}) && grep {$parameter->{'name'}->{'value'} eq $_} &List::ExtraUtils::as_array($options{'exclude'}) );
    ++$i;    
    my $definition = $parameterMean->(($i))->string();
    for(my $j=0;$j<nelem($parameterMean);++$j) {
	my $coefficient = $cholesky->(($i),($j));
	$definition .= "+\%[".$options{'metaName'}.$j."]*".$coefficient
	    if ( $coefficient != 0.0 );
    }
    my $name  = exists($parameterNameMap{$parameter->{'name'}->{'value'}}) ? $parameterNameMap{$parameter->{'name'}->{'value'}} : $parameter->{'name' }->{'value'};
    my %parameterOutput = 
	(
	 value      => "derived",
         name       => {value => $name      },
	 definition => {value => $definition}
	);
    push(@{$parametersOutput->{'parameter'}},\%parameterOutput);
}
# Append meta-parameters.
for(my $i=0;$i<nelem($parameterMean);++$i) {
    my %parameterOutput =
        (
	 value => "active",
         name                            => 
	 {
	     value => $options{'metaName'}.$i
	 },
	 distributionFunction1DPerturber => 
	 {
	     value  => "cauchy",
	     scale  => {value => 1.0e-5},
	     median => {value => 0.0e+0}
	 },
	 operatorUnaryMapper => 
	 {
	     value => "identity"
	 },
	 distributionFunction1DPrior => 
	 {
	     value      => "normal",
	     mean       => {value => +0.0},
	     variance   => {value => +1.0},
	     limitLower => {value => -5.0},
	     limitUpper => {value => +5.0}	     
	 }	 
        );
    push(@{$parametersOutput->{'modelParameterMethod'}},\%parameterOutput);
}

# Output the parameter definition file.
my $xml = new XML::Simple();
open(my $outputFile,">".$options{'outputFile'});
print $outputFile $xml->XMLout($parametersOutput, RootName => "parameters", KeyAttr => ["value"]);
close($outputFile);

exit 0;
