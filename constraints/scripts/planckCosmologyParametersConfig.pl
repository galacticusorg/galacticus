#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::LinearAlgebra;
use PDL::IO::Misc;
use PDL::MatrixOps;
use UNIVERSAL 'isa';
use Data::Dumper;

# Generates a configuration for cosmological parameters based on the Planck covariance matrix and suitable for use with the
# Galacticus+MCMC constraints architecture. Parameters are expressed as linear combiantions of six independent normal deviates
# (labelled "cosmology0" through "cosmology5"). The output of this script can be dropped into the "parameters" section of a
# Galacticus+MCMC config file.
# Andrew Benson (16-October-2014)

# Download Planck MCMC chains data set.
system("wget \"http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_CosmoParams_base_planck_lowl_lowLike_highL_post_lensing_R1.10.tar.gz\" -O ".$galacticusPath."aux/COM_CosmoParams_base_planck_lowl_lowLike_highL_post_lensing_R1.10.tar.gz")
    unless ( -e $galacticusPath."aux/COM_CosmoParams_base_planck_lowl_lowLike_highL_post_lensing_R1.10.tar.gz" );

# Unpack Planck MCMC chains data set.
system("cd ".$galacticusPath."aux; tar xvfz COM_CosmoParams_base_planck_lowl_lowLike_highL_post_lensing_R1.10.tar.gz")
    unless ( -e $galacticusPath."aux/base_planck_lowl_lowLike_highL_post_lensing" );

# Specify Planck directory.
my $planckDirectoryName = $galacticusPath."aux/base_planck_lowl_lowLike_highL_post_lensing/base/planck_lowl_lowLike_highL";

# Specify parameter names of interest.
my %parameterMap = 
    (
     "omegabh2"  => "Omega_b",
     "omegamh2*" => "Omega_Matter",
     "tau"       => "reionizationSuppressionOpticalDepth",
     "ns"        => "powerSpectrumIndex",
     "H0*"       => "H_0",
     "sigma8*"   => "sigma_8"
    );

# Read Planck parameter names.
my $parameters;
open(my $parameterNameFile,$planckDirectoryName."/base_planck_lowl_lowLike_highL_post_lensing.paramnames");
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
    if ( $fileName =~ m/base_planck_lowl_lowLike_highL_post_lensing_\d+\.txt/ ) {
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
    $parameterName = $parameterMap{$parameter}
          unless ( $parameterMap{$parameter} eq "" );
    # If the parameter is a density parameter then it has been multiplied by h^2 in the Planck analysis. Undo this.
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
