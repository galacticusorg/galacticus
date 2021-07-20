#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use System::Redirect;
use File::NFSLock;
use Fcntl qw(:DEFAULT :flock);
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use File::Slurp;
use Galacticus::Constraints::Parameters;
use Galacticus::Options;
use List::ExtraUtils;
use Scalar::Util 'reftype';

# Run the current maximum likelihood model from a constraint run.
# Andrew Benson (04-December-2012)

# Get command line arguments.
die("Usage: maximumLikelihoodModel.pl <parameterFile> [options...]")
    unless ( scalar(@ARGV) >= 1 );
my $parameterFile = $ARGV[0];
# Create a hash of named options.
my %options = (
    runModel   => "yes",
    directory  => "maximumLikelihoodModel",
    chain      => "all",
    launcher   => "local",
    label      => "maximumLikelihood",
    mpi        => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::parseConfig($parameterFile);
 
# Extract MCMC directory.
my $logFileRoot = $config->{'posteriorSampleSimulation'}->{'logFileRoot'}->{'value'};
(my $mcmcDirectory  = $logFileRoot) =~ s/\/[^\/]+$//;    

# Find number of steps taken.
my $stepCount = &Galacticus::Constraints::Parameters::stepCount($config,\%options);

# Add a parameter with the current step count.
my @derivedParameters = &List::ExtraUtils::as_array($options{'derivedParameter'})
    if ( exists($options{'derivedParameter'}) );
push(@derivedParameters,"posteriorSimulationStep=".$stepCount);
delete($options{'derivedParameter'});
@{$options{'derivedParameter'}} = @derivedParameters;

# Get maximum likelihood parameter vector.
(my $maximumPosteriorParameters, my $maximumPosterior) = &Galacticus::Constraints::Parameters::maximumPosteriorParameterVector($config,\%options);
print "Maximum posterior: ".$maximumPosterior."\n";
print " Model parameters:\n".join("\t",$maximumPosteriorParameters->list())."\n";

# Get a list of likelihood models to evaluate.
my @models     = &Galacticus::Constraints::Parameters::modelList($config,\%options);
my @jobs           ;
my $modelCount = -1;
foreach my $model ( @models ) {
    # Determine the maximum likelihood model directory.
    my $maximumLikelihoodDirectory  = ($options{'directory'} =~ m/^\// ? $options{'directory'} : $mcmcDirectory."/".$options{'directory'}).(scalar(@models) == 1 ? "" : "_".(++$modelCount))."/";
    # Get the parameters.
    my $parameters = &Galacticus::Constraints::Parameters::parameterVectorApply($config,$model,$maximumPosteriorParameters,\%options);
    # Set an output file name.
    $parameters->{'galacticusOutputFileName'}->{'value'} = $maximumLikelihoodDirectory."/galacticus.hdf5";
    # Write the modified parameters to file.
    &Galacticus::Constraints::Parameters::parametersOutput($parameters,$maximumLikelihoodDirectory."/parameters.xml");    
    # Finish if we're not to run the model.
    next
     	if ( $options{'runModel'} eq "no" );
    # Stage the Galacticus model.
    push(@jobs,&Galacticus::Constraints::Parameters::stageModel($maximumLikelihoodDirectory."/parameters.xml",$options{'label'},\%options));
}

# Launch all models.
&Galacticus::Constraints::Parameters::launchModel(\@jobs,\%options);

exit 0;
