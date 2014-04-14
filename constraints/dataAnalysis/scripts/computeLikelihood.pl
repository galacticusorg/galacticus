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
use Data::Dumper;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
require System::Redirect;

# Compute the likelihood for a conditional mass function model to fit mass function data.
# Andrew Benson (06-July-2012)

# Parse the constraint config file for parameters.
die("Usage: computeLikelihood.pl <configFile> <mpiRank> <param1> [<param2>......]")
    unless ( scalar(@ARGV) > 2 );

# Parse the config file.
my $configFile = $ARGV[0];
my $xml        = new XML::Simple;
my $config     = $xml->XMLin($configFile, KeyAttr => 0);

# Get the MPI rank.
my $mpiRank = sprintf("%4.4d",$ARGV[1]);

# Define parameter names.
my @parameters = (
		  "conditionalMassFunctionBehrooziAlphaSatellite",
		  "conditionalMassFunctionBehrooziLog10M1",
		  "conditionalMassFunctionBehrooziLog10Mstar0",
		  "conditionalMassFunctionBehrooziBeta",
		  "conditionalMassFunctionBehrooziDelta",
		  "conditionalMassFunctionBehrooziGamma",
		  "conditionalMassFunctionBehrooziSigmaLogMstar",
		  "conditionalMassFunctionBehrooziBCut",
		  "conditionalMassFunctionBehrooziBSatellite",
		  "conditionalMassFunctionBehrooziBetaCut",
		  "conditionalMassFunctionBehrooziBetaSatellite"
		  );

my $parameterCount = scalar(@parameters);

# Read the base parameters.
my $newParameters = $xml->XMLin($config->{'parameterFile'});

# Create an array of new parameters.
die("computeLikelihood.pl: number of supplied arguments does not match number of parameters")
    unless ( scalar(@ARGV) == $parameterCount+2 );
for(my $i=0;$i<$parameterCount;++$i) {
    $newParameters->{'parameter'}->{$parameters[$i]}->{'value'} = $ARGV[$i+2];
}

# Set the output file name.
$newParameters->{'parameter'}->{'conditionalMassFunctionOutputFileName'}->{'value'} = $config->{'workDirectory'}."massFunction_".$mpiRank.".hdf5";

# Clean up.
unlink(
       $config->{'workDirectory'}."glcLikelihood_".$mpiRank.".xml",
       $config->{'workDirectory'}."massFunction_".$mpiRank.".hdf5"
       );

# Output the parameter file.
my $xmlOut = new XML::Simple(RootName=>"parameters", NoAttr => 1);
open(oHndl,">".$config->{'workDirectory'}."glcLikelihood_".$mpiRank.".xml");
print oHndl $xmlOut->XMLout($newParameters);
close(oHndl);

# Generate the mass function.
SystemRedirect::tofile("Conditional_Mass_Function.exe ".$config->{'workDirectory'}."glcLikelihood_".$mpiRank.".xml","/dev/null");
system("mkdir -p ".$config->{'workDirectory'}."failures; mv -f ".$config->{'workDirectory'}."glcLikelihood_".$mpiRank.".xml ".$config->{'workDirectory'}."failures/")
    unless ( $? == 0 );

# Read mass function.
my $massFunctionFile  = new PDL::IO::HDF5($config->{'workDirectory'}."massFunction_".$mpiRank.".hdf5");
my $massFunctionModel = $massFunctionFile->dataset('massFunction')->get();

# Read the covariance matrix and observed mass function.
my $covarianceFile       = new PDL::IO::HDF5($config->{'constraintFile'});
my $massFunctionObserved = $covarianceFile->dataset('massFunctionObserved'    )->get();
my $logDeterminant       = $covarianceFile->dataset('logDeterminantCovariance')->get();
my $inverseCovariance    = $covarianceFile->dataset('inverseCovariance'       )->get();

# Compute the likelihood.
my $difference       = $massFunctionModel-$massFunctionObserved;
my $vCv              = $difference x $inverseCovariance x transpose($difference);
my $logLikelihood    = -0.5*$vCv((0),(0))-0.5*nelem($massFunctionObserved)*log(2.0*3.1415927)-0.5*$logDeterminant;

# Output the likelihood.
my $likelihoodFile = $config->{'workDirectory'}."glcLikelihood.dat_".$mpiRank;
open(oHndl,">".$likelihoodFile);
print oHndl $logLikelihood."\n";
close(oHndl);

exit;
