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
    foreach ( @{$options{'parameterName'}} ) {
	if ( $_ =~ m/([^:]+):([^:]+):([^:]+)/ ) {
	    my $inputName   = $1;
	    my $outputName  = $2;
	    my $outputLabel = $3;
	    $parameterNameMap{$inputName} = {name => $outputName, label => $outputLabel};
	} else {
	    die("fitMultivariateNormal.pl: parameterName option has incorrect syntax");
	}
    }
}

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFile);

# Read a sample matrix of all chain states.
my $sampleMatrix = &Galacticus::Constraints::Parameters::Sample_Matrix($config,\%options);

# Compute the mean parameters.
my $parameterMean = $sampleMatrix->xchg(0,1)->average();

# Compute the covariance of parameters.
my $parameterCovariance = pdl zeroes(nelem($parameterMean),nelem($parameterMean));
for(my $i=0;$i<nelem($parameterMean);++$i) {
    for(my $j=$i;$j<nelem($parameterMean);++$j) {
	my $covariance = 
	    sum(
		+($sampleMatrix->(($i),:)-$parameterMean->(($i)))
		*($sampleMatrix->(($j),:)-$parameterMean->(($j)))
	    )
	    /$sampleMatrix->dim(1);
	$parameterCovariance->(($i),($j)) .= $covariance;
	$parameterCovariance->(($j),($i)) .= $covariance;
    }
}

# Find Cholesky decomposition of the covariance matrix.
my $cholesky = mchol($parameterCovariance);
# Construct a parameter definition file to describe this multivariate normal distribution.
my $i = -1;
my $parametersOutput;
foreach my $parameter ( @{$config->{'parameters'}->{'parameter'}} ) {
    next 
	unless ( exists($parameter->{'prior'}) );
    ++$i;
    my $definition      = $parameterMean->(($i))->string();
    for(my $j=0;$j<nelem($parameterMean);++$j) {
	my $coefficient = $cholesky->(($i),($j));
	$definition .= "+\%[".$options{'metaName'}.$j."]*".$coefficient
	    if ( $coefficient != 0.0 );
    }
    my $name  = exists($parameterNameMap{$parameter->{'name'}}) ? $parameterNameMap{$parameter->{'name'}}->{'name' } : $parameter->{'name' };
    my $label = exists($parameterNameMap{$parameter->{'name'}}) ? $parameterNameMap{$parameter->{'name'}}->{'label'} : $parameter->{'label'};
    my %parameterOutput = 
	(
	 name   => $name      ,
	 label  => $label     ,
	 define => $definition
	);
    push(@{$parametersOutput->{'parameter'}},\%parameterOutput);
}
# Append meta-parameters.
for(my $i=0;$i<nelem($parameterMean);++$i) {
    my %parameterOutput =
        (
         name   => $options{'metaName'}.$i,
	 label  => "\\hbox{".$options{'metaName'}.$i."}",
	 random => {
	     scale => 1.0e-5,
	     type => "Cauchy",
	     median => 0.0
	 },
	 mapping => {
	     type => "linear"
	 },
	 prior => {
	     distribution => {
		 type     => "normal",
		 mean     => +0.0,
		 variance => +1.0,
		 minimum  => -5.0,
		 maximum  => +5.0
	     }
	 }	 
        );
    push(@{$parametersOutput->{'parameter'}},\%parameterOutput);
}

# Output the parameter definition file.
my $xml = new XML::Simple();
open(my $outputFile,">".$options{'outputFile'});
print $outputFile $xml->XMLout($parametersOutput, RootName => "parameters", NoAttr => 1);
close($outputFile);

exit 0;
