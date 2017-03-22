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
die("Usage: maximumLikelihoodModel.pl <configFile> [options...]")
    unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = (
    runModel   => "yes",
    directory  => "maximumLikelihoodModel",
    fixedTrees => "config",
    chain      => "all"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFile);

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
my $logFileRoot = $config->{'simulation'}->{'logFileRoot'};
(my $mcmcDirectory  = $logFileRoot) =~ s/\/[^\/]+$//;

# Determine the maximum likelihood model directory.
my $maximumLikelihoodDirectory  = $mcmcDirectory."/".$arguments{'directory'}."/";

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($config->{'likelihood'}->{'compilation'},$config->{'likelihood'}->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Set an output file name.
system("mkdir -p ".$maximumLikelihoodDirectory);
$parameters->{'galacticusOutputFileName'}->{'value'} = $maximumLikelihoodDirectory."/galacticus.hdf5";

# Set a random number seed.
$parameters->{'randomSeed'}->{'value'} = int(rand(10000))+1
    if ( exists($config->{'likelihood'}->{'randomize'}) && $config->{'likelihood'}->{'randomize'} eq "yes" );

# Determine number of chains.
my $chainCount = 0;
while () {
    ++$chainCount;
    my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,$chainCount);
    last
	unless ( -e $chainFileName );
}

# Parse the chains to find the maximum likelihood model.
my $maximumLikelihood = -1e30;
my @maximumLikelihoodParameters;
for(my $i=0;$i<$chainCount;++$i) {
    next
	unless
	(
	 $arguments{'chain'} eq "all"
	 ||
	 $arguments{'chain'} == $i
	);
    open(iHndl,sprintf("%s_%4.4i.log",$logFileRoot,$i));
    while ( my $line = <iHndl> ) {
	unless ( $line =~ m/^\"/ ) {
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    my @columns = split(/\s+/,$line);
	    if ( $columns[4] > $maximumLikelihood ) {
		$maximumLikelihood           = $columns[4];
		@maximumLikelihoodParameters = @columns[6..$#columns];
	    }
	}
    }
    close(iHndl);
}
print "Maximum likelihood: ".$maximumLikelihood."\n";
print "  Model parameters:\n".join("\t",@maximumLikelihoodParameters)."\n";

# Convert these values into a parameter array.
my $newParameters = &Galacticus::Constraints::Parameters::Convert_Parameters_To_Galacticus($config,@maximumLikelihoodParameters);

# Apply to parameters.
for my $newParameterName ( keys(%{$newParameters}) ) {
    my $parameter = $parameters;
    foreach ( split(/\-\>/,$newParameterName) ) {
	$parameter->{$_}->{'value'} = undef()
	    unless ( exists($parameter->{$_}) );
	$parameter = $parameter->{$_};
    }
    $parameter->{'value'} = $newParameters->{$newParameterName};
}

# Apply any parameters from command line.
foreach my $argument ( keys(%arguments) ) {
    if ( $argument =~ m/^parameter:(.*)/ ) {
	my @parametersSet = split(":",$1);
	my $parameter     = $parameters;
	foreach my $parameterSet ( @parametersSet ) {
	    if ( $parameterSet =~ m/([^\{]*)(\{{0,1}([^\}]*)\}{0,1})/ ) {
		my $parameterName  = $1;
		my $parameterValue = $3;
		$parameter->{$parameterName}->{'value'} = $parameterValue
		    if ( defined($2) );
		$parameter = $parameter->{$parameterName};
	    } else {
		die("malformed parameter definition");
	    }
	}
	$parameter->{'value'} = $arguments{$argument};
    }
}

# If fixed sets of trees are to be used, modify parameters to use them.
my $useFixedTrees = $arguments{'fixedTrees'} eq "config" ? $config->{'likelihood'}->{'useFixedTrees'} : $arguments{'fixedTrees'};
if ( $useFixedTrees eq "yes" ) {
    # Get a lock on the tree file.
    my $fixedTreeDirectory;
    if ( exists($config->{'likelihood'}->{'fixedTreesInScratch'}) && $config->{'likelihood'}->{'fixedTreesInScratch'} eq "yes" ) {
	$fixedTreeDirectory = $config->{'likelihood'}->{'scratchDirectory'}."/";
    } else {
	$fixedTreeDirectory = $config->{'likelihood'}->{'workDirectory'   }."/trees/";
    }
    my $fixedTreeFile      = $fixedTreeDirectory."fixedTrees".$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".hdf5";
    # Modify parameters to use the tree file.
    $parameters->{'mergerTreeConstructMethod'}->{'value'} = "read";
    $parameters->{'mergerTreeReadFileName'   }->{'value'} = $fixedTreeFile;
}

# Write the modified parameters to file.
system("mkdir -p ".$maximumLikelihoodDirectory);
&Galacticus::Constraints::Parameters::Output($parameters,$maximumLikelihoodDirectory."/parameters.xml");

# Finish if we're not to run the model.
exit
    if ( $arguments{'runModel'} eq "no" );

# Run the Galacticus model.
system("make Galacticus.exe")
    unless ( -e "Galacticus.exe" );
die("maximumLikelihoodModel.pl: failed to build Galacticus.exe")
    unless ( $? == 0 );
my $glcCommand;
$glcCommand .= "time ./Galacticus.exe ".$maximumLikelihoodDirectory."/parameters.xml";
my $logFile = $maximumLikelihoodDirectory."/galacticus.log";
&System::Redirect::tofile($glcCommand,$logFile);
die("maximumLikelihoodModel.pl: Galacticus model failed")
    unless ( $? == 0 );

# Store raw XML parameter file in the model.
&storeXML($maximumLikelihoodDirectory."/galacticus.hdf5",$maximumLikelihoodDirectory."/parameters.xml");

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
    my $options = "";
    $options = " ".$constraintDefinition->{'analysisArguments'}
        if ( exists($constraintDefinition->{'analysisArguments'}) );
    $options .= " --modelDiscrepancies ".$workDirectory."/modelDiscrepancy"
	if ( -e $workDirectory."/modelDiscrepancy" );
    system($analysisCode." ".$maximumLikelihoodDirectory."/galacticus.hdf5 --outputFile ".$maximumLikelihoodDirectory."/likelihood".ucfirst($constraintDefinition->{'label'}).".xml --plotFile ".$maximumLikelihoodDirectory."/".$plotFile.".pdf --resultsFile ".$maximumLikelihoodDirectory."/".$constraintDefinition->{'label'}.":results.hdf5".$options);
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

sub storeXML {
    my $galacticusFileName = shift();
    my $parametersFileName = shift();
    my $galacticusModel = new PDL::IO::HDF5(">".$galacticusFileName);
    my $parametersGroup = $galacticusModel->group('Parameters');
    my $parametersRaw   = read_file($parametersFileName);
    $parametersGroup->attrSet('rawXML' => $parametersRaw);    
}
