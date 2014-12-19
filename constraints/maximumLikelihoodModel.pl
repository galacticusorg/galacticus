#!/usr/bin/env perl
use strict;
use warnings;
use lib './perl';
use XML::Simple;
use Data::Dumper;
use System::Redirect;
use File::NFSLock;
use Fcntl qw(:DEFAULT :flock);
require Galacticus::Constraints::Parameters;
require Galacticus::Options;

# Run the current maximum likelihood model from a constraint run.
# Andrew Benson (04-December-2012)

# Get command line arguments.
die("Usage: maximumLikelihoodModel.pl <configFile> [options...]")
    unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = (
		 runModel => "yes"
		 );
&Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Parameters::Parse_Config($configFile);

# Validate the config file.
die("maximumLikelihoodModel.pl: workDirectory must be specified in config file" ) 
    unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("maximumLikelihoodModel.pl: compilation must be specified in config file"   )
    unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("maximumLikelihoodModel.pl: baseParameters must be specified in config file") 
    unless ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Determine the work directory.
my $workDirectory  = $config->{'likelihood'}->{'workDirectory'};

# Determine the MCMC directory.
my $mcmcDirectory  = $workDirectory."/mcmc/";

# Determine the maximum likelihood model directory.
my $maximumLikelihoodDirectory  = $mcmcDirectory."/maximumLikelihoodModel/";

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'likelihood'}->{'compilation'},$config->{'likelihood'}->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Set an output file name.
system("mkdir -p ".$workDirectory."/maximumLikelihoodModel");
$parameters->{'parameter'}->{'galacticusOutputFileName'}->{'value'} = $maximumLikelihoodDirectory."/galacticus.hdf5";

# Set a random number seed.
$parameters->{'parameter'}->{'randomSeed'}->{'value'} = int(rand(10000))+1
    if ( exists($config->{'likelihood'}->{'randomize'}) && $config->{'likelihood'}->{'randomize'} eq "yes" );

# Determine number of chains.
my $chainCount = 0;
while () {
    ++$chainCount;
    my $chainFileName = sprintf("%s/chains_%4.4i.log",$mcmcDirectory,$chainCount);
    last
	unless ( -e $chainFileName );
}

# Parse the chains to find the maximum likelihood model.
my $maximumLikelihood = -1e30;
my @maximumLikelihoodParameters;
for(my $i=0;$i<$chainCount;++$i) {
    open(iHndl,$mcmcDirectory."/chains_".sprintf("%4.4i",$i).".log");
    while ( my $line = <iHndl> ) {
	unless ( $line =~ m/^\"/ ) {
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    my @columns = split(/\s+/,$line);
	    if ( $columns[4] > $maximumLikelihood ) {
		$maximumLikelihood           = $columns[4];
		@maximumLikelihoodParameters = @columns[5..$#columns];
	    }
	}
    }
    close(iHndl);
}
print "Maximum likelihood: ".$maximumLikelihood."\n";
print "  Model parameters:\n".join("\t",@maximumLikelihoodParameters)."\n";

# Convert these values into a parameter array.
my $newParameters = &Parameters::Convert_Parameters_To_Galacticus($config,@maximumLikelihoodParameters);

# Apply to parameters.
$parameters->{'parameter'}->{$_->{'name'}}->{'value'} = $_->{'value'}
    foreach ( @{$newParameters->{'parameter'}} );

# Write the modified parameters to file.
system("mkdir -p ".$maximumLikelihoodDirectory);
&Parameters::Output($parameters,$maximumLikelihoodDirectory."/parameters.xml");

# Finish if we're not to run the model.
exit
    if ( $arguments{'runModel'} eq "no" );

# Run the Galacticus model.
system("make Galacticus.exe")
    unless ( -e "Galacticus.exe" );
die("maximumLikelihoodModel.pl: failed to build Galacticus.exe")
    unless ( $? == 0 );
my $glcCommand;
$glcCommand .= "./Galacticus.exe ".$maximumLikelihoodDirectory."/parameters.xml";
my $logFile = $maximumLikelihoodDirectory."/galacticus.log";
&SystemRedirect::tofile($glcCommand,$logFile);
die("maximumLikelihoodModel.pl: Galacticus model failed")
    unless ( $? == 0 );

# Perform processing of the model, accumulating likelihood as we go.
my $logLikelihood = 0.0;
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $xml = new XML::Simple;
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Run the analysis code.
    my $analysisCode = $constraintDefinition->{'analysis'};
    my $plotFile = $constraintDefinition->{'label'};
    $plotFile =~ s/\./_/g;
    system($analysisCode." ".$maximumLikelihoodDirectory."/galacticus.hdf5 --outputFile ".$maximumLikelihoodDirectory."/likelihood".ucfirst($constraintDefinition->{'label'}).".xml --modelDiscrepancies ".$mcmcDirectory."/modelDiscrepancy --plotFile ".$maximumLikelihoodDirectory."/".$plotFile.".pdf --resultsFile ".$maximumLikelihoodDirectory."/".$constraintDefinition->{'label'}.":results.hdf5");
    die("maximumLikelihoodModel.pl: analysis code failed")
	unless ( $? == 0 );
    # Read the likelihood.
    my $likelihood = $xml->XMLin($maximumLikelihoodDirectory."/likelihood".ucfirst($constraintDefinition->{'label'}).".xml");
    die("maximumLikelihoodModel.pl: likelihood calculation failed")
	if ( $likelihood->{'logLikelihood'} eq "nan" );
    # Extract the likelihood and weight it.
    my $thisLogLikelihood = $likelihood->{'logLikelihood'};
    $thisLogLikelihood *= $constraintDefinition->{'weight'}
        if ( defined($constraintDefinition->{'weight'}) );
    # Accumulate the likelihood.
    $logLikelihood += $thisLogLikelihood;
}

# Write the final likelihood to file.
open(oHndl,">".$maximumLikelihoodDirectory."/likelihood.xml");
print oHndl $logLikelihood."\n";
close(oHndl);

exit;
