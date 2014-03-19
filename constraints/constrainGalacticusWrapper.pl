#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use XML::Simple;
use UNIVERSAL;
use Data::Dumper;
use Sys::CPU;
use List::Util qw(min max);
require Galacticus::Constraints::Parameters;

# A wrapper script to be called by the Galacticus likelihood function of the Galacticus constraints infrastructure to drive a
# Galacticus model and compute likelihood.
# Andrew Benson (04-October-2011)

# Parse the constraint config file for parameters.
die("Usage: constrainGalacticusWrapper.pl <configFile> <mpiRank> <likelihoodFile> <temperature> <param1> [<param2>......]") 
    unless ( scalar(@ARGV) > 4 );
my $configFile = $ARGV[0];
my $config     = &Parameters::Parse_Config($configFile);
my @parameters;
if ( UNIVERSAL::isa($config->{'parameters'}->{'parameter'},"ARRAY") ) {
    @parameters = @{$config->{'parameters'}->{'parameter'}};
} else {
    push(@parameters,$config->{'parameters'}->{'parameter'});
}
# Count active parameters.
my $parameterCount = 0;
for(my $i=0;$i<scalar(@parameters);++$i) {
    ++$parameterCount 
	if ( exists($parameters[$i]->{'prior'}) );
}

# Get the MPI rank.
my $mpiRank        = $ARGV[1];

# Get the name for the likelihood file.
my $likelihoodFile = $ARGV[2];

# Get the temperature.
my $temperature    = $ARGV[3];

# Convert command line arguments to a parameter structure.
die("constrainGalacticusWrapper.pl: number of supplied arguments does not match number of parameters") 
    unless ( scalar(@ARGV) == $parameterCount+4 );

my $j = -1;
my %parameterValues;
for(my $i=0;$i<scalar(@parameters);++$i) {
    if ( exists($parameters[$i]->{'prior'}) ) {
	++$j;
	$parameterValues{$parameters[$i]->{'name'}} = $ARGV[$j+4];
    }
}

# Set the values of any parameters that are defined in terms of other parameters.
my $failCount = 1;
while ( $failCount > 0 ) {
    $failCount = 0;
    for(my $i=0;$i<scalar(@parameters);++$i) {
	if ( exists($parameters[$i]->{'define'}) ) {
	    die ("constrainGalacticusWrapper.pl: cannot specify a prior for a defined parameter")
		if ( exists($parameters[$i]->{'prior'}) );
	    # Attempt to replace named parameters in the definition with their values.
	    while ( $parameters[$i]->{'define'} =~ m/\%([a-zA-Z0-9_]+)/ ) {
		my $parameterName = $1;
		if ( exists($parameterValues{$parameterName}) ) {
		    $parameters[$i]->{'define'} =~ s/\%$parameterName/$parameterValues{$parameterName}/g;
		} else {
		    ++$failCount;
		    last;
		}
		$parameterValues{$parameters[$i]->{'name'}} = eval($parameters[$i]->{'define'})
		    unless ( $parameters[$i]->{'define'} =~ m/\%([a-zA-Z0-9_]+)/ );
	    }
	}
    }
}

# Create an array of new parameters.
my $newParameters;
for(my $i=0;$i<scalar(@parameters);++$i) {
    push(
	 @{$newParameters->{'parameter'}},
	 {
	     name  =>                  $parameters[$i]->{'name'} ,
	     value => $parameterValues{$parameters[$i]->{'name'}}
	 }
	 );
}

# Find the scratch directory.
my $scratchDirectory = $config->{'likelihood'}->{'workDirectory'}."/mcmc";
$scratchDirectory = $config->{'likelihood'}->{'scratchDirectory'}
    if ( exists($config->{'likelihood'}->{'scratchDirectory'}) );

# Expand any environment variable names in the scratch directory.
while ( $scratchDirectory =~ m/\$([_A-Z]+)/ ) {
    my $environmentVariableName  = $1;
    my $environmentVariableValue = $ENV{$environmentVariableName};
    $scratchDirectory =~ s/\$$environmentVariableName/$environmentVariableValue/g;
}

# Ensure scratch and work directories exist.
system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'}."/mcmc")
    unless ( -e $config->{'likelihood'}->{'workDirectory'}."/mcmc" );
system("mkdir -p ".$scratchDirectory)
    unless ( -e $scratchDirectory );

# Report.
if ( exists($config->{'likelihood'}->{'report'}) ) {
    if ( $config->{'likelihood'}->{'report'} eq "yes" ) {
	print "Report from constrainGalacticusWrapper.pl:\n";
	print "  MPI rank is : ".$mpiRank."\n";
	print "  Output in   : ".$scratchDirectory."/newParameters_".$mpiRank.".xml\n";
	print "  Parameters  : \n";
	print Dumper($newParameters);
    }
}

# Output the parameter file.
my $xmlOut = new XML::Simple(RootName=>"parameters", NoAttr => 1);
open(oHndl,">".$scratchDirectory."/newParameters_".$mpiRank.".xml");
print oHndl $xmlOut->XMLout($newParameters);
close(oHndl);

# Find the set of base parameters to use.
my $baseParameters = "parameters.xml";
$baseParameters = $config->{'likelihood'}->{'baseParameters'}
   if ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Run the Galacticus constraint script.
my $runCommand = "./constraints/constrainGalacticus.pl ".$configFile." ".$scratchDirectory." ".$config->{'likelihood'}->{'compilation'}." ".$scratchDirectory."/newParameters_".$mpiRank.".xml ".$config->{'likelihood'}->{'workDirectory'}." --baseParameters ".$baseParameters." --make no --timing no --output ".$likelihoodFile." --galacticusFile constrainGalacticus_".$mpiRank.".hdf5 --suffix ".$mpiRank." --temperature ".$temperature." --failCount ".$config->{'likelihood'}->{'workDirectory'}."/failureCount.txt --failArchive ".$config->{'likelihood'}->{'workDirectory'}."/failure";
if ( exists($config->{'likelihood'}->{'cpulimit'}) ) {
    my $cpuLimit = $config->{'likelihood'}->{'cpulimit'};
    if ( exists($config->{'likelihood'}->{'threads'}) ) { 
	$cpuLimit *= $config->{'likelihood'}->{'threads'};
    } else {
	$cpuLimit *= Sys::CPU::cpu_count();
    }
    $runCommand .= " --cpulimit ".$cpuLimit;
}
if ( exists($config->{'likelihood'}->{'threads'}) ) {
    $runCommand .= " --threads ".$config->{'likelihood'}->{'threads'};
}
if ( exists($config->{'likelihood'}->{'saveState'}) ) {
    $runCommand .= " --saveState ".$config->{'likelihood'}->{'saveState'};
}
if ( exists($config->{'likelihood'}->{'storeResults'}) ) {
    $runCommand .= " --storeResults ".$config->{'likelihood'}->{'workDirectory'}."/mcmc/"
	if ( $config->{'likelihood'}->{'storeResults'} eq "yes" );
}
$runCommand .= " --cleanUp "   .$config->{'likelihood'}->{'cleanUp'   } if ( exists($config->{'likelihood'}->{'cleanUp'   }) );
$runCommand .= " --randomize " .$config->{'likelihood'}->{'randomize' } if ( exists($config->{'likelihood'}->{'randomize' }) );
system($runCommand);

exit;
