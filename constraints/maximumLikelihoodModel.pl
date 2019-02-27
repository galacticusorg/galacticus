#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
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
use List::ExtraUtils;
use Scalar::Util 'reftype';

# Run the current maximum likelihood model from a constraint run.
# Andrew Benson (04-December-2012)

# Get command line arguments.
die("Usage: maximumLikelihoodModel.pl <parameterFile> [options...]")
    unless ( scalar(@ARGV) >= 1 );
my $parameterFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments = (
    runModel   => "yes",
    directory  => "maximumLikelihoodModel",
    fixedTrees => "config",
    chain      => "all",
    mpi        => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);
# Parse the constraint config file.
my $xml = new XML::Simple();
my $config = $xml->XMLin($parameterFile);
# Determine the MCMC directory.
my $logFileRoot = $config->{'posteriorSampleSimulationMethod'}->{'logFileRoot'}->{'value'};
(my $mcmcDirectory  = $logFileRoot) =~ s/\/[^\/]+$//;    
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
my @parameterDefinitions;
my @modelParameters = &List::ExtraUtils::as_array($config->{'posteriorSampleSimulationMethod'}->{'modelParameterMethod'});
my $activeParameter = -1;
foreach my $modelParameter ( @modelParameters ) {
    if ( $modelParameter->{'value'} eq "active" ) {
	++$activeParameter;
	push(@parameterDefinitions,$modelParameter->{'name'}->{'value'}."=".$maximumLikelihoodParameters[$activeParameter]);
    } elsif ( $modelParameter->{'value'} eq "derived" ) {
	push(@parameterDefinitions,$modelParameter->{'name'}->{'value'}."==".$modelParameter->{'definition'}->{'value'});
    }
}
my $newParameters = &Galacticus::Constraints::Parameters::Convert_Parameters_To_Galacticus(\@parameterDefinitions);
# Begin building a stack of likelihood models to evaluate.
my @allParameterNames = map {$_->{'name'}->{'value'}} &List::ExtraUtils::as_array($config->{'posteriorSampleSimulationMethod'}->{'modelParameterMethod'});
my @modelStack = ( 
    {
	posteriorSampleLikelihoodMethod => $config->{'posteriorSampleLikelihoodMethod'},
	parameters                      => \@allParameterNames
    } 
    );
my $singleModel = $config->{'posteriorSampleLikelihoodMethod'}->{'value'} eq "galaxyPopulation";
# Iterate over likelihood models.
my $modelCount = 0;
while ( scalar(@modelStack) > 0 ) {
    # Pop a model off of the stack.
    my $model = shift(@modelStack);
    # Add any sub-models to the stack.
    if ( exists($model->{'posteriorSampleLikelihoodMethod'}->{'posteriorSampleLikelihoodMethod'}) ) {
	my @likelihoods           = &List::ExtraUtils::as_array($model->{'posteriorSampleLikelihoodMethod'}->{'posteriorSampleLikelihoodMethod'});
	my @parameterMapsActive   = &List::ExtraUtils::as_array($model->{'posteriorSampleLikelihoodMethod'}->{'parameterMap'                   })
	    if ( exists($model->{'posteriorSampleLikelihoodMethod'}->{'parameterMap'        }) );
	my @parameterMapsInactive = &List::ExtraUtils::as_array($model->{'posteriorSampleLikelihoodMethod'}->{'parameterInactiveMap'           })
	    if ( exists($model->{'posteriorSampleLikelihoodMethod'}->{'parameterInactiveMap'}) );
	for(my $i=0;$i<scalar(@likelihoods);++$i) {
	    my @parameterNames;
	    push(@parameterNames,split(" ",$parameterMapsActive  [$i]->{'value'}))
		if ( @parameterMapsActive   );
	    push(@parameterNames,split(" ",$parameterMapsInactive[$i]->{'value'}))
		if ( @parameterMapsInactive );
	    my $childModel;
	    $childModel->{'posteriorSampleLikelihoodMethod'} = $likelihoods[$i];
	    $childModel->{'parameters'} = \@parameterNames;
	    push(
		@modelStack,
		$childModel
		);
	}
    }
    # Ignore non-Galacticus likelihood models.
    next
	unless ( $model->{'posteriorSampleLikelihoodMethod'}->{'value'} eq "galaxyPopulation" );
    # Determine the maximum likelihood model directory.
    my $maximumLikelihoodDirectory  = ($arguments{'directory'} =~ m/^\// ? $arguments{'directory'} : $mcmcDirectory."/".$arguments{'directory'});
    ++$modelCount;
    $maximumLikelihoodDirectory .= "_".$modelCount
	unless ( $singleModel );
    $maximumLikelihoodDirectory .= "/";
    # Get the parameters structure.
    my $parameters = $xml->XMLin($model->{'posteriorSampleLikelihoodMethod'}->{'baseParametersFileName'}->{'value'});
    # Apply any command line parameters.
    &Galacticus::Constraints::Parameters::Apply_Command_Line_Parameters($parameters,\%arguments);
    # Set an output file name.
    system("mkdir -p ".$maximumLikelihoodDirectory);
    $parameters->{'galacticusOutputFileName'}->{'value'} = $maximumLikelihoodDirectory."/galacticus.hdf5";
    # Apply to parameters.
    for my $newParameterName ( keys(%{$newParameters}) ) {
	next
	    unless ( grep {$_ eq $newParameterName} @{$model->{'parameters'}} );
	my $parameter = $parameters;
	my $valueIndex;
	if ( $newParameterName =~ m/^(.*)\{(\d+)\}$/ ) {
	    $newParameterName = $1;
	    $valueIndex = $2;
	}
	foreach ( split(/::/,$newParameterName) ) {
	    # Check if the parameter name contains an array reference.
	    if ( $_ =~ m/^(.*)\[(\d+)\]$/ ) {
		# Parameter name contains array reference. Step through to the relevant parameter in the list. If the parameter is
		# not an array, allow this only if the array index given is zero.
		if ( reftype($parameter->{$1}) eq "ARRAY" ) {
		    $parameter->{$1}->[$2]->{'value'} = undef()
			unless ( scalar(@{$parameter->{$1}}) > $2 );
		    $parameter = $parameter->{$1}->[$2];
		} else {
		    die('maximumLikelihoodModel.pl: attempt to access non-existant array')
			unless ( $2 == 0 );
		    $parameter->{$1}->{'value'} = undef()
			unless ( exists($parameter->{$1}) );
		    $parameter = $parameter->{$1};
		}
	    } else {
		# Parameter does not contain an array reference - so simply step through to the named parameter.
		$parameter->{$_}->{'value'} = undef()
		    unless ( exists($parameter->{$_}) );
		$parameter = $parameter->{$_};
	    }
	}
	# Test if the parameter name contains a value index.
	if ( defined($valueIndex) ) {
	    # A value index is given - set the relevant entry.
	    my @values = split(" ",$parameter->{'value'});
	    $values[$valueIndex] = $newParameters->{$newParameterName."{".$valueIndex."}"};
	    $parameter->{'value'} = join(" ",@values);
	} else {
	    # No value index is given - simply set the value of the parameter.
	    $parameter->{'value'} = $newParameters->{$newParameterName};
	}    
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
    # Write the modified parameters to file.
    system("mkdir -p ".$maximumLikelihoodDirectory);
    &Galacticus::Constraints::Parameters::Output($parameters,$maximumLikelihoodDirectory."/parameters.xml");    
    # Finish if we're not to run the model.
    next
	if ( $arguments{'runModel'} eq "no" );
    # Run the Galacticus model.
    my $executable = exists($model->{'posteriorSampleLikelihoodMethod'}->{'executable'}) ? $model->{'posteriorSampleLikelihoodMethod'}->{'executable'}->{'value'} : "Galacticus.exe";
    $executable = $arguments{'executable'}
        if ( exists($arguments{'executable'}) );
    my $glcCommand = ($arguments{'mpi'} > 0 ? "mpirun -np ".$arguments{'mpi'}." " : "").$executable." ".$maximumLikelihoodDirectory."/parameters.xml";
    my $logFile = $maximumLikelihoodDirectory."/galacticus.log";
    &System::Redirect::tofile($glcCommand,$logFile);
    die("maximumLikelihoodModel.pl: Galacticus model failed")
	unless ( $? == 0 );
    # Store raw XML parameter file in the model.
    &storeXML($maximumLikelihoodDirectory."/galacticus".($arguments{'mpi'} > 0 ? ":MPI".sprintf("%4.4i",$arguments{'mpi'}) : "").".hdf5",$maximumLikelihoodDirectory."/parameters.xml");
}
exit;

sub storeXML {
    my $galacticusFileName = shift();
    my $parametersFileName = shift();
    my $galacticusModel    = new PDL::IO::HDF5(">".$galacticusFileName);
    my $parametersGroup    = $galacticusModel->group('Parameters');
    my $parametersRaw      = read_file($parametersFileName);
    $parametersGroup->attrSet('rawXML' => $parametersRaw);    
}
