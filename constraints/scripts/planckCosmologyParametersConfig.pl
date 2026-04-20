#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::IO::Misc;
use PDL::MatrixOps;
use Data::Dumper;
use File::Temp qw(tempdir);

# Generates a configuration for cosmological parameters based on the Planck covariance matrix and suitable for use with the
# Galacticus MCMC constraints architecture. Parameters are expressed as linear combiantions of six independent normal deviates
# (labelled "cosmology0" through "cosmology5"). The output of this script can be included into the base parameter file for a
# Galacticus MCMC run, and actual cosmological parameters set by referencing those created in this file.
# Andrew Benson (16-October-2014)
 
# Get a temporary directory.
my $dir = tempdir( CLEANUP => 1 );

# Download Planck MCMC chains data set.
system("wget \"http://pla.esac.esa.int/pla-sl/data-action?COSMOLOGY.COSMOLOGY_OID=1200\" -O ".$dir."/COM_CosmoParams_base-plikHM-TT-lowTEB_R2.00.tar.gz")
    unless ( -e $dir."/COM_CosmoParams_base-plikHM-TT-lowTEB_R2.00.tar.gz" );

# Unpack Planck MCMC chains data set.
system("cd ".$dir."; tar xvfz COM_CosmoParams_base-plikHM-TT-lowTEB_R2.00.tar.gz")
    unless ( -e $dir."/base/plikHM_TT_lowTEB_lensing" );

# Specify Planck directory.
my $planckDirectoryName = $dir."/base/plikHM_TT_lowTEB_lensing/";

# Specify parameter names of interest.
my %parameterMap = 
    (
     "omegabh2"  => "planckCosmologyOmegaBaryon",
     "omegamh2*" => "planckCosmologyOmegaMatter",
     "tau"       => "planckCosmologyReionizationSuppressionOpticalDepth",
     "ns"        => "planckCosmologyPowerSpectrumIndex",
     "H0*"       => "planckCosmologyHubbleConstant",
     "sigma8*"   => "planckCosmologySigma8"
    );

# Read Planck parameter names.
my $parameters;
open(my $parameterNameFile,$planckDirectoryName."/base_plikHM_TT_lowTEB_lensing_post_BAO_H070p6_JLA.paramnames");
my $i = -1;
while ( my $line = <$parameterNameFile> ) {
    ++$i;
    $line =~ s/^\s*//;
    my @columns = split(/\s+/,$line);
    if ( exists($parameterMap{$columns[0]}) ) {
	$parameters->{$columns[0]}->{'column'} = $i+2;
	$parameters->{$columns[0]}->{'chain' } = pdl [];
    }
}
close($parameterNameFile);

# Check all required parameters were found.
foreach ( keys(%parameterMap) ) {
    die("planckCosmologyParametersConfig.pl: parameter '".$_."' not found")
	unless ( exists($parameters->{$_}) );
}

# Read Planck chains.
my @columnsToRead = map {$parameters->{$_}->{'column'}} sort(keys(%parameterMap));
opendir(my $planckDirectory,$planckDirectoryName);
while ( my $fileName = readdir($planckDirectory) ) {
    if ( $fileName =~ m/base_plikHM_TT_lowTEB_lensing_post_BAO_H070p6_JLA_\d+\.txt/ ) {
	my @columnData = rcols($planckDirectoryName."/".$fileName,@columnsToRead);
	# Append to chains.
	my $i = -1;
	foreach ( sort(keys(%parameterMap)) ) {
	    ++$i;
	    $parameters->{$_}->{'chain'} = $parameters->{$_}->{'chain'}->append($columnData[$i]);
	}
    }
}
closedir($planckDirectory);

# Find the means of the parameters.
my $parameterCount = scalar(keys(%parameterMap));
my $mean           = pdl zeroes($parameterCount);
$i                 = -1;
foreach ( sort(keys(%parameterMap)) ) {
    ++$i;
    $mean->(($i)) .= average($parameters->{$_}->{'chain'});
}

# Find the covariances of the parameters.
my $covariance = pdl zeroes($parameterCount,$parameterCount);
$i             = -1;
foreach my $iKey ( sort(keys(%parameterMap)) ) {
    ++$i;
    my$j = -1;
    foreach my $jKey ( sort(keys(%parameterMap)) ) {
	++$j;
	$covariance->(($i),($j)) .= 
	    average(
		($parameters->{$iKey}->{'chain'}-$mean->(($i)))
		*
		($parameters->{$jKey}->{'chain'}-$mean->(($j)))
	    );
    }
}

# Perform a Cholesky decomposition on the covariance matrix.
my $choleskyDecomposed = mchol($covariance);

# Construct the config data.
my $parameterConfig;
my $index = -1;
foreach my $parameter ( sort(keys(%parameterMap)) ) {
    # Find the index of this parameter in the vectors.
    ++$index;
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
    $parameterName = $parameterMap{$parameter}
          unless ( $parameterMap{$parameter} eq "" );
    # If the parameter is a density parameter then it has been multiplied by h^2 in the Planck analysis. Undo this.
    $value = "(".$value.")/(%[planckCosmologyHubbleConstant]/100.0)**2"
 	if ( $parameter =~ m/omega/ );
    # Add limits.
    $value = "min(".$value.",1.0)"
	if ( $parameter eq "omegamh2*" );
    $value = "max(".$value.",0.0)"
	if ( $parameter eq "tau" );    
    # Store the parameter config in our array.
    $parameterConfig->{$parameterName} = {"value" => $value};
}

# Output the configuration.
my $xmlOut = new XML::Simple (RootName=>"parameters");
print $xmlOut->XMLout($parameterConfig);

exit;
