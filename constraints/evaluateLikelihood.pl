#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use XML::Simple;
use Data::Dumper;
use System::Redirect;
use File::NFSLock;
use Fcntl qw(:DEFAULT :flock);
use PDL;
use PDL::IO::HDF5;
use File::Slurp;
use Galacticus::Constraints::Parameters;
use Galacticus::Options;

# Run the current maximum likelihood model from a constraint run.
# Andrew Benson (04-December-2012)

# Get command line arguments.
die("Usage: evaluateLikelihood.pl <configFile> <modelFile>")
    unless ( scalar(@ARGV) >= 2 );
my $configFileName = $ARGV[0];
my $modelFileName  = $ARGV[1];
# Create a hash of named arguments.
my %arguments;
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFileName);

# Validate the config file.
die("maximumLikelihoodModel.pl: workDirectory must be specified in config file" ) 
    unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("maximumLikelihoodModel.pl: compilation must be specified in config file"   )
    unless ( exists($config->{'likelihood'}->{'compilation'   }) );

# Determine the work directory.
my $workDirectory  = $config->{'likelihood'}->{'workDirectory'};

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($config->{'likelihood'}->{'compilation'},$config->{'likelihood'}->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Perform processing of the model, accumulating likelihood as we go.
my $logLikelihood = 0.0;
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $xml = new XML::Simple;
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Run the analysis code.
    my $analysisCode = $constraintDefinition->{'analysis'};
    my $options = "";
    $options = " ".$constraintDefinition->{'analysisArguments'}
        if ( exists($constraintDefinition->{'analysisArguments'}) );
    $options .= " --modelDiscrepancies ".$workDirectory."/modelDiscrepancy"
	if ( -e $workDirectory."/modelDiscrepancy" );
    system($analysisCode." ".$modelFileName." --outputFile /dev/shm/likelihood.tmp.".ucfirst($constraintDefinition->{'label'}).".xml ".$options);
    die("maximumLikelihoodModel.pl: analysis code failed")
	unless ( $? == 0 );
    # Read the likelihood.
    my $likelihood = $xml->XMLin("/dev/shm/likelihood.tmp.".ucfirst($constraintDefinition->{'label'}).".xml");
    die("maximumLikelihoodModel.pl: likelihood calculation failed")
	if ( $likelihood->{'logLikelihood'} eq "nan" );
    print $constraintDefinition->{'label'}."\t".$likelihood->{'logLikelihood'}."\n";
    # Extract the likelihood and weight it.
    my $thisLogLikelihood = $likelihood->{'logLikelihood'};
    $thisLogLikelihood *= $constraintDefinition->{'weight'}
        if ( defined($constraintDefinition->{'weight'}) );
    # Accumulate the likelihood.
    $logLikelihood += $thisLogLikelihood;
}

# Write the final likelihood.
print "Likelihood = ".$logLikelihood."\n";

exit 0;
