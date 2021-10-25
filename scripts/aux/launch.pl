#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Data::Dumper;
use Clone qw(clone);
use Digest::MD5 qw(md5_hex);
use Scalar::Util qw(reftype);
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::Local;
use Galacticus::Launch::PBS;
use Galacticus::Launch::MonolithicPBS;
use Galacticus::Launch::Slurm;
use List::ExtraUtils;
use Scalar::Util qw(looks_like_number reftype);

# Script to launch sets of Galacticus models, iterating over sets of parameters and performing analysis
# on the results. Supports launching on a variety of platforms via modules.
# Andrew Benson (11-June-2010; major re-write 01-February-2014)

# Get command line arguments.
die("Usage: launch.pl <runFile>")
    unless ( scalar(@ARGV) >= 1 );
my $launchFileName = $ARGV[0];

# Extract options.
my %arguments = 
    (
     instance => "1:1"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%arguments);

# Parse the launch script.
my $launchScript          = &Parse_Launch_Script            ($launchFileName);

# Read in any configuration options.
$launchScript->{'config'} = &Galacticus::Options::LoadConfig(               );

# Override launch method if specified on the command line.
$launchScript->{'launchMethod'} = $arguments{'launchMethod'}
    if ( exists($arguments{'launchMethod'}) );

# Check for an instance number for this launch.
if ( $arguments{"instance"} =~ m/(\d+):(\d+)/ ) {
    $launchScript->{'thisInstance' } = $1;
    $launchScript->{'instanceCount'} = $2;
    print " -> launching instance ".$launchScript->{'thisInstance'}." of ".$launchScript->{'instanceCount'}."\n"
	if ( $launchScript->{'verbosity'} > 0 );
} else {
    die("launch.pl: 'instance' argument syntax error");
}

# Validate the launch method.
die("launch.pl: unrecognized launch method")
    unless ( exists($Galacticus::Launch::Hooks::moduleHooks{$launchScript->{'launchMethod'}}) );
&{$Galacticus::Launch::Hooks::moduleHooks{$launchScript->{'launchMethod'}}->{'validate'}}($launchScript);

# Construct models.
my @jobs = &Construct_Models($launchScript);

# Launch models.
&{$Galacticus::Launch::Hooks::moduleHooks{$launchScript->{'launchMethod'}}->{'launch'}}
		      (\@jobs,$launchScript,\%arguments);

exit;

sub Construct_Models {
    # Constructs model directories and parameter files. Returns an array of model jobs.
    # Get arguments script.
    my $launchScript = shift();
    # Set initial value of random seed.
    my $randomSeed   = 219;
    # Initialize a model counter.
    $launchScript->{'modelCounter'} = -1;
    # Initialize array of model jobs.
    my @jobs;
    # Loop through all model sets.
    my $iModelSet    = -1;
    foreach my $parameterSet ( @{$launchScript->{'parameters'}} ) {
	# Increment model set counter.
	++$iModelSet;	
	# Set base model name.
	my $modelBaseName = "galacticus";
	$modelBaseName    = $parameterSet->{'label'}
	   if ( exists($parameterSet->{'label'}) );
	# Create an array of hashes giving the parameters for this parameter set.
	my @parameterSets = &unfoldParameters($parameterSet);
	# Loop over all models and run them.
	my $iModel     = 0;
	my $mergeGroup = 0;
	foreach my $parameterData ( @parameterSets ) {
	    ++$mergeGroup;
	    # Split the job across multiple workers if needed.
	    for(my $workerInstance=1;$workerInstance<=$launchScript->{'splitModels'};++$workerInstance) {
		# Increment model counters.
		++$iModel;
		++$launchScript->{'modelCounter'};
		# Increment random seed.
		++$randomSeed;
		# Evaluate on which instance this model should be launched.
		my $runOnInstance = ($launchScript->{'modelCounter'} % $launchScript->{'instanceCount'}) + 1;	    
		# Set worker instance if needed.
		if ( $launchScript->{'splitModels'} > 1 ) {
		    $parameterData->{'treeEvolveWorkerNumber'}->{'value'} = $workerInstance;
		    $parameterData->{'treeEvolveWorkerCount' }->{'value'} = $launchScript->{'splitModels'};
		}
		# Specify the output directory.
		my $modelLabel    = exists($parameterData->{'label'}) ? $parameterData->{'label'}->{'value'} : $iModelSet.":".$iModel;
		my $galacticusOutputDirectory = $launchScript->{'modelRootDirectory'}
		."/".$modelBaseName."_".$modelLabel;
		# Change output directory name to an md5 hash if so requested.
		my $descriptor;
		if ( $launchScript->{'md5Names'} eq "yes" ) {
		    my $md5 = md5_hex(Dumper($parameterData));
		    $galacticusOutputDirectory =  $launchScript->{'modelRootDirectory'}
		    ."/".$modelBaseName."_".$md5;
		}
		# If the output directory does not exist, then create it.
		unless (
		    -e $galacticusOutputDirectory 
		    ||
		    $runOnInstance != $launchScript->{'thisInstance'}
		    ) {
		    system("mkdir -p ".$galacticusOutputDirectory);
		    # Output the MD5 descriptor to this directory.
		    if ( defined($descriptor) ) {
			open(oHndl,">".$galacticusOutputDirectory."/md5Descriptor.txt");
			print oHndl $descriptor."\n";
			close(oHndl);
		    }
		    # Specify the output file.
		    my $galacticusOutputFile = $galacticusOutputDirectory."/galacticus.hdf5";		
		    # Read the default set of parameters.
		    my $parameters;
		    unless ( $launchScript->{'baseParameters'} eq "" ) {
			my $xml  = new XML::Simple;
			$parameters = $xml->XMLin($launchScript->{'baseParameters'});
		    }
		    # Set the output file name.
		    $parameters->{'galacticusOutputFileName'}->{'value'} 
		        = &{$Galacticus::Launch::Hooks::moduleHooks{$launchScript->{'launchMethod'}}->{'outputFileName'}}
		           ($galacticusOutputFile,$launchScript);
		    # Set the random seed.
		    $parameters->{'randomSeed'}->{'value'} = $randomSeed 
			unless ( exists($parameters->{'randomSeed'}) );
		    # Set a state restore file.
		    if ( $launchScript->{'useStateFile'} eq "yes" ) {
			(my $stateFile = $parameters->{'galacticusOutputFileName'}->{'value'}) =~ s/\.hdf5//;
			$parameters->{'stateFileRoot'}->{'value'} = $stateFile;
		    }
		    # Transfer parameters for this model to the active parameter set.
		    foreach my $parameter ( keys(%{$parameterData}) ) {
			$parameters->{$parameter} = $parameterData->{$parameter};
		    }
		    # Output the parameters as an XML file.
		    my $xmlOutput = new XML::Simple (RootName=>"parameters");
		    open(outHndl,">".$galacticusOutputDirectory."/parameters.xml");
		    print outHndl $xmlOutput->XMLout($parameters);
		    close(outHndl);
		    undef($parameters);
		    # Generate analysis code for this job.
		    my $analysisCode;
		    $analysisCode = $launchScript->{'analysisScript'}." ".$galacticusOutputDirectory."\n"
			if ( $launchScript->{'doAnalysis'} eq "yes" );
		    # Push the job onto the stack.
		    push(
			@jobs,
			{
			    label        => $modelLabel                    ,
			    directory    => $galacticusOutputDirectory     ,
			    analysis     => $analysisCode                  ,
			    mergeGroup   => $mergeGroup                    ,
			    modelCounter => $launchScript->{'modelCounter'}
			}
			);
		}	    
	    }
	}
    }
    return @jobs;
}

sub unfoldParameters {
    # Create an array of hashes which give the parameter values for each model.
    my $parameterSet = shift();
    # Walk through the parameter structure, identifying any parameters which are arrays.
    my @parametersIn  = (
	{
	    parameter => $parameterSet,
	    parent    => undef(),
	    name      => undef()
	}
	);
    my @parametersIntermediate;
    while ( scalar(@parametersIn) > 0 ) {
	# Get a parameter set from the input list.
	my $parameters = shift(@parametersIn);
	# Walk through all parameters.
	my @parameterStack = ( $parameters );
	my $cloned         = 0;
	while ( scalar(@parameterStack) > 0 ) {
	    my $parameter = shift(@parameterStack);
	    if ( ! reftype($parameter->{'parameter'}) ) {
		# Simple parameter lacking subparameters - ignore.
	    } elsif ( reftype($parameter->{'parameter'}) eq "ARRAY" ) {
		# Duplicate the parameter structure, making a copy for each array element and push back onto the stack, unless
		# explicitly forbidden, in which case just push the parameter onto the stack.
		my $i             = -1;
		my $clonesAdded   =  0;
		my $parameterCopy = clone($parameter->{'parameter'});
		foreach my $element ( @{$parameterCopy} ) {
		    ++$i;
		    if ( exists($element->{'iterable'}) && $element->{'iterable'} eq "no" ) {
			# Element is non-iterable, push onto the stack.
			push(
			    @parameterStack,
			    {
				parameter => $parameter->{'parameter'}->[$i],
				parent    => $parameter->{'parent'   }      ,
				name      => $parameter->{'name'     } 
			    }
			    );
		    } else {
			# Element is iterable - clone parameters.
			$parameter->{'parent'}->{$parameter->{'name'}} = $element;
			push(@parametersIn,clone($parameters));			
			$clonesAdded = 1;
		    }
		}
		if ( $clonesAdded ) {
		    $cloned = 1;
		    last;
		}
	    } elsif ( reftype($parameter->{'parameter'}) eq "HASH" ) {
		# Contains sub-parameters. Push them to the parameter stack.
		push(
		    @parameterStack,
		    {
			parameter => $parameter->{'parameter'}->{$_},
			parent    => $parameter->{'parameter'}      ,
			name      =>                             $_
		    }
		    )
		    foreach ( keys(%{$parameter->{'parameter'}}) );
	    }
	}
	push(@parametersIntermediate,$parameters->{'parameter'})
	    unless ( $cloned );
    }
    # Parameter sets all now have single values of each parameter. Search for cases where a parameter's position in the parameter
    # hierarchy should be moved.
    my @parametersOut;
    while ( scalar(@parametersIntermediate) > 0 ) {
	my $parameters = shift(@parametersIntermediate);
	my $modified   = 0;
  	# Walk through all parameters.
	my @parameterStack = 
	    (
	     {
		 parameter => $parameters,
		 parent    => undef(),
		 name      => undef()
	     }
	    );
	# Build a stack that we can traverse depth-first.
	my @parameterStackNew;
	while ( scalar(@parameterStack) > 0 ) {
	    my $parameter = shift(@parameterStack);
	    if ( reftype($parameter->{'parameter'}) && reftype($parameter->{'parameter'}) eq "HASH" ) {
		# Contains sub-parameters. Push them to the parameter stack.
		push(@parameterStackNew,$parameter);
		push(
		    @parameterStack,
		    {
			parameter => $parameter->{'parameter'}->{$_},
			parent    => $parameter                     ,
			name      =>                             $_
		    }
		    )
		    foreach ( keys(%{$parameter->{'parameter'}}) );
	    }
	}
      	# Step through the stack handling any parameterLevel attributes.
	while ( scalar(@parameterStackNew) > 0 ) {
	    my $parameter = pop(@parameterStackNew);
	    if ( exists($parameter->{'parameter'}->{'parameterLevel'}) ) {
		if ( looks_like_number($parameter->{'parameter'}->{'parameterLevel'}) && $parameter->{'parameter'}->{'parameterLevel'} < 0 ) {
		    # Move this parameter up.
		    my $parameterTarget = $parameter;
		    for(;$parameter->{'parameter'}->{'parameterLevel'}<=0;++$parameter->{'parameter'}->{'parameterLevel'}) {
			$parameterTarget = $parameterTarget->{'parent'};
		    }
		    delete($parameter->{'parameter'}->{'parameterLevel'});
		    delete($parameter->{'parent'}->{'parameter'}->{$parameter->{'name'}});
		    delete($parameterTarget->{'parameter'}->{$parameter->{'name'}});
		    $parameterTarget->{'parameter'}->{$parameter->{'name'}} = $parameter->{'parameter'};
		} elsif ( $parameter->{'parameter'}->{'parameterLevel'} eq "top" ) {
		    # Move this parameter to the top of the parameter hierarchy.
		    delete($parameter->{'parameter'}->{'parameterLevel'});
		    $parameters->{$parameter->{'name'}} = $parameter->{'parameter'};
		    delete($parameter->{'parent'}->{'parameter'}->{$parameter->{'name'}});
		} else {
		    die("unknown parameterLevel");
		}
		# Mark the parameters as modified and exit the parameter walk.
		$modified = 1;
		last;
	    }
	}
	# Transfer parameters to the output list unless modified, in which case push them back onto our stack for further
	# processing.
	if ( $modified ) {
	    push(@parametersIntermediate,$parameters);
	} else {
	    push(@parametersOut         ,$parameters);
	}
    }
    # Return the list of parameters.
    return @parametersOut;
}

sub Parse_Launch_Script {
    # Parse the launch script.
    my $launchFileName = shift;
    # Load XML.
    my $xml          = new XML::Simple;
    my $launchScript = $xml->XMLin($launchFileName, KeyAttr => "", ForceArray => [ "modify", "value", "parameter", "parameters", "requirement" ]);
    # Assign defaults.
    my %defaults = 
	(
	 verbosity          => 0         ,
	 md5Names           => "no"      ,
	 useStateFile       => "no"      ,
	 compressModels     => "no"      ,
	 launchMethod       => "local"   ,
	 modelRootDirectory => "./models",
	 baseParameters     => ""        ,
	 doAnalysis         => "no"      ,
	 splitModels        => 1
	);
    foreach ( keys(%defaults) ) {
	$launchScript->{$_} = $defaults{$_}
	unless ( exists($launchScript->{$_}) );
    }
    # Return the script.
    return $launchScript;
}
