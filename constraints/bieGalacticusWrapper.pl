#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
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

# A wrapper script to be called by the GalacticusLikelihoodFunction of the Bayesian Inference Engine to drive a Galacticus model
# and compute likelihood.
# Andrew Benson (04-October-2011)

# Parse the constraint config file for parameters.
die("Usage: bieGalacticusWrapper.pl <configFile> <mpiRank> <param1> [<param2>......]") unless ( scalar(@ARGV) > 2 );
my $configFile = $ARGV[0];
my $config = &Parameters::Parse_Config($configFile);
my @parameters;
if ( UNIVERSAL::isa($config->{'parameters'}->{'parameter'},"ARRAY") ) {
    @parameters = @{$config->{'parameters'}->{'parameter'}};
} else {
    push(@parameters,$config->{'parameters'}->{'parameter'});
}
# Count active parameters.
my $parameterCount = 0;
for(my $i=0;$i<scalar(@parameters);++$i) {
    ++$parameterCount if ( exists($parameters[$i]->{'prior'}) );
}

# Get the MPI rank.
my $mpiRank = sprintf("%4.4d",$ARGV[1]);

# Convert command line arguments to a parameter structure.
die("bieGalacticusWrapper.pl: number of supplied arguments does not match number of parameters") unless ( scalar(@ARGV) == $parameterCount+2 );
my $j = -1;
my %parameterValues;
for(my $i=0;$i<scalar(@parameters);++$i) {
    if ( exists($parameters[$i]->{'prior'}) ) {
	++$j;
	$parameterValues{$parameters[$i]->{'name'}} = $ARGV[$j+2];
	# Unmap from logarithmic space if necessary.
	if ( $parameters[$i]->{'prior'}->{'distribution'} eq "uniform" ) {
	    if ( exists($parameters[$i]->{'prior'}->{'mapping'}) ) {
		if ( $parameters[$i]->{'prior'}->{'mapping'} eq "logarithmic" ) {
		    $parameterValues{$parameters[$i]->{'name'}} = exp($parameterValues{$parameters[$i]->{'name'}});
		}
	    }
	}
    }
}

# Set the values of any parameters that are defined in terms of other parameters.
my $failCount = 1;
while ( $failCount > 0 ) {
    $failCount = 0;
    for(my $i=0;$i<scalar(@parameters);++$i) {
	if ( exists($parameters[$i]->{'define'}) ) {
	    die ("bieGalacticusWrapper.pl: cannot specify a prior for a defined parameter")
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
	     name  => $parameters[$i  ]->{'name'},
	     value => $parameterValues{$parameters[$i]->{'name'}}
	 }
	 );
}

# Find the temperature.
my $temperature = 1.0;
$temperature = $config->{'temperature'}
    if ( exists($config->{'temperature'}) );

# Find the scratch directory.
my $scratchDirectory = $config->{'workDirectory'}."/mcmc";
$scratchDirectory = $config->{'scratchDirectory'}
    if ( exists($config->{'scratchDirectory'}) );

# Expand any environment variable names in the scratch directory.
while ( $scratchDirectory =~ m/\$([_A-Z]+)/ ) {
    my $environmentVariableName  = $1;
    my $environmentVariableValue = $ENV{$environmentVariableName};
    $scratchDirectory =~ s/\$$environmentVariableName/$environmentVariableValue/g;
}

# Ensure scratch and work directories exist.
system("mkdir -p ".$config->{'workDirectory'}."/mcmc")
    unless ( -e $config->{'workDirectory'}."/mcmc" );
system("mkdir -p ".$scratchDirectory)
    unless ( -e $scratchDirectory );

# Report.
if ( exists($config->{'report'}) ) {
    if ( $config->{'report'} eq "yes" ) {
	print "Report from bieGalacticusWrapper.pl:\n";
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
$baseParameters = $config->{'baseParameters'} if ( exists($config->{'baseParameters'}) );

# Run the Galacticus constraint script.
my $runCommand = "./constraints/constrainGalacticus.pl ".$scratchDirectory." ".$config->{'compilation'}." ".$scratchDirectory."/newParameters_".$mpiRank.".xml ".$config->{'workDirectory'}." --baseParameters ".$baseParameters." --make no --timing no --output ".$scratchDirectory."/glcLikelihood.dat_".$mpiRank." --galacticusFile constrainGalacticus_".$mpiRank.".hdf5 --suffix ".$mpiRank." --temperature ".$temperature." --failCount ".$config->{'workDirectory'}."/failureCount.txt --failArchive ".$config->{'workDirectory'}."/failure";
if ( exists($config->{'cpulimit'}) ) {
    my $cpuLimit = $config->{'cpulimit'};
    if ( exists($config->{'galacticusThreads'}) ) { 
	$cpuLimit *= $config->{'galacticusThreads'};
    } else {
	$cpuLimit *= Sys::CPU::cpu_count();
    }
    $runCommand .= " --cpulimit ".$cpuLimit;
}
if ( exists($config->{'galacticusThreads'}) ) {
    $runCommand .= " --galacticusThreads ".$config->{'galacticusThreads'};
}
if ( exists($config->{'galacticusSaveState'}) ) {
    $runCommand .= " --galacticusSaveState ".$config->{'galacticusSaveState'};
}
if ( exists($config->{'storeResults'}) ) {
    $runCommand .= " --storeResults ".$config->{'workDirectory'}."/mcmc/"
	if ( $config->{'storeResults'} eq "yes" );
}
$runCommand .= " --reuseTrees ".$config->{'reuseTrees'} if ( exists($config->{'reuseTrees'}) );
$runCommand .= " --cleanUp ".$config->{'cleanUp'} if ( exists($config->{'cleanUp'}) );
$runCommand .= " --randomize ".$config->{'randomize'} if ( exists($config->{'randomize'}) );
system($runCommand);

exit;
