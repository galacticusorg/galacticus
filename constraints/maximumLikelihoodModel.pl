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

# Run the current maximum likelihood model from a Galacticus+BIE constraint run.
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
die("maximumLikelihoodModel.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'workDirectory' }) );
die("maximumLikelihoodModel.pl: compilation must be specified in config file"   ) unless ( exists($config->{'compilation'   }) );
die("maximumLikelihoodModel.pl: baseParameters must be specified in config file") unless ( exists($config->{'baseParameters'}) );

# Determine the work directory.
my $workDirectory  = $config->{'workDirectory'};

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'compilation'},$config->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Set an output file name.
system("mkdir -p ".$workDirectory."/maximumLikelihoodModel");
$parameters->{'parameter'}->{'galacticusOutputFileName'}->{'value'} = $workDirectory."/maximumLikelihoodModel/galacticus.hdf5";

# Set a random number seed.
$parameters->{'parameter'}->{'randomSeed'}->{'value'} = int(rand(10000))+1
    if ( exists($config->{'randomize'}) && $config->{'randomize'} eq "yes" );

# Parse the statefile to find the maximum likelihood model.
my $maximumLikelihood = -1e30;
my @maximumLikelihoodParameters;
open(iHndl,$workDirectory."/mcmc/galacticusBIE.statelog");
while ( my $line = <iHndl> ) {
    unless ( $line =~ m/^\"/ ) {
	$line =~ s/^\s*//;
	$line =~ s/\s*$//;
	my @columns = split(/\s+/,$line);
	if ( $columns[2] > $maximumLikelihood ) {
	    $maximumLikelihood = $columns[2];
	    @maximumLikelihoodParameters = @columns[5..$#columns];
	}
    }
}
close(iHndl);
print $maximumLikelihood."\n";

# Convert these values into a parameter array.
my $newParameters = &Parameters::Convert_BIE_Parameters_To_Galacticus($config,@maximumLikelihoodParameters);

# Apply to parameters.
$parameters->{'parameter'}->{$_->{'name'}}->{'value'} = $_->{'value'}
    foreach ( @{$newParameters->{'parameter'}} );

# Write the modified parameters to file.
&Parameters::Output($parameters,$workDirectory."/maximumLikelihoodModel/parameters.xml");

# Finish if we're not to run the model.
exit
    if ( $arguments{'runModel'} eq "no" );

# Run the Galacticus model.
system("make Galacticus.exe")
    unless ( -e "Galacticus.exe" );
die("maximumLikelihoodModel.pl: failed to build Galacticus.exe")
    unless ( $? == 0 );
my $glcCommand;
$glcCommand .= "./Galacticus.exe ".$workDirectory."/maximumLikelihoodModel/parameters.xml";
my $logFile = $workDirectory."/maximumLikelihoodModel/galacticus.log";
SystemRedirect::tofile($glcCommand,$logFile);
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
    system($analysisCode." ".$workDirectory."/maximumLikelihoodModel/galacticus.hdf5 --outputFile ".$workDirectory."/maximumLikelihoodModel/likelihood".ucfirst($constraintDefinition->{'label'}).".xml --modelDiscrepancies ".$workDirectory."/modelDiscrepancy --plotFile ".$workDirectory."/maximumLikelihoodModel/".$plotFile.".pdf --accuracyFile ".$workDirectory."/maximumLikelihoodModel/".$constraintDefinition->{'label'}.":accuracy.xml");
    die("maximumLikelihoodModel.pl: analysis code failed")
	unless ( $? == 0 );
    # Read the likelihood.
    my $likelihood = $xml->XMLin($workDirectory."/maximumLikelihoodModel/likelihood".ucfirst($constraintDefinition->{'label'}).".xml");
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
open(oHndl,">".$workDirectory."/maximumLikelihoodModel/likelihood.xml");
print oHndl $logLikelihood."\n";
close(oHndl);

exit;
