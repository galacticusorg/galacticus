#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use Data::Dumper;
use Clone qw(clone);
use Digest::MD5 qw(md5_hex);
use Scalar::Util qw(reftype);
require Galacticus::Options;
require Galacticus::Launch::Hooks;
require Galacticus::Launch::Local;
require Galacticus::Launch::PBS;
require Galacticus::Launch::MonolithicPBS;
require Galacticus::Launch::Slurm;
require List::ExtraUtils;

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
&Options::Parse_Options(\@ARGV,\%arguments);

# Parse the launch script.
my $launchScript          = &Parse_Launch_Script    ($launchFileName);

# Read in any configuration options.
$launchScript->{'config'} = &Parse_Galacticus_Config(               );

# Check for an instance number for this launch.
if ( $arguments{"instance"} =~ m/(\d+):(\d+)/ ) {
    $launchScript->{'thisInstance' }  = $1;
    $launchScript->{'instanceCount'} = $2;
    print " -> launching instance ".$launchScript->{'thisInstance'}." of ".$launchScript->{'instanceCount'}."\n"
	if ( $launchScript->{'verbosity'} > 0 );
} else {
    die("launch.pl: 'instance' argument syntax error");
}

# Validate the launch method.
die("launch.pl: unrecognized launch method")
    unless ( exists($Hooks::moduleHooks{$launchScript->{'launchMethod'}}) );
&{$Hooks::moduleHooks{$launchScript->{'launchMethod'}}->{'validate'}}($launchScript);

# Construct models.
my @jobs = &Construct_Models($launchScript);
## AJB HACK
exit;

# Launch models.
&{$Hooks::moduleHooks{$launchScript->{'launchMethod'}}->{'launch'}}
		      (\@jobs,$launchScript);

exit;

sub Construct_Models {
    # Constructs model directories and parameter files. Returns an array of model jobs.
    # Get arguments script.
    my $launchScript  = shift;
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
	my @parameterHashes = &Create_Parameter_Hashes($parameterSet);
	# Loop over all models and run them.
	my $iModel     = 0;
	my $mergeGroup = 0;
	foreach my $parameterData ( @parameterHashes ) {
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
		my $modelLabel    = exists($parameterData->{'label'}) ? $parameterData->{'label'}->[0]->{'value'} : $iModelSet.":".$iModel;
		my $galacticusOutputDirectory =  $launchScript->{'modelRootDirectory'}
		."/".$modelBaseName."_".$modelLabel;
		# Change output directory name to an md5 hash if so requested.
		my $descriptor;
		if ( $launchScript->{'md5Names'} eq "yes" ) {
		    foreach my $parameter ( keys(%{$parameterData}) ) {
			$descriptor .= $parameter.":".$parameterData->{$parameter}->{'value'}.":";
		    }
		    my $md5 = md5_hex($descriptor);
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
			my $data = $xml->XMLin($launchScript->{'baseParameters'});
			foreach my $parameterName ( keys(%{$data}) ) {
			    $parameters->{$parameterName}->[0] = $data->{$parameterName}
			        if ( reftype($data->{$parameterName}) );
			}
		    }
		    # Set the output file name.
		    $parameters->{'galacticusOutputFileName'}->[0]->{'value'} 
		        = &{$Hooks::moduleHooks{$launchScript->{'launchMethod'}}->{'outputFileName'}}
		           ($galacticusOutputFile,$launchScript);
		    # Set the random seed.
		    $parameters->{'randomSeed'}->[0]->{'value'} = $randomSeed 
			unless ( exists($parameters->{'randomSeed'}) );
		    # Set a state restore file.
		    if ( $launchScript->{'useStateFile'} eq "yes" ) {
			(my $stateFile = $parameters->{'galacticusOutputFileName'}->[0]->{'value'}) =~ s/\.hdf5//;
			$parameters->{'stateFileRoot'}->[0]->{'value'} = $stateFile;
		    }
		    # Transfer parameters for this model from the array of model parameter hashes to the
		    # active hash.
		    foreach my $parameter ( keys(%{$parameterData}) ) {
			$parameters->{$parameter} = $parameterData->{$parameter};
		    }
		    # Transfer values from the active hash to an array suitable for XML output.
		    my $data;
		    foreach my $name ( sort(keys(%{$parameters})) ) {
			$data->{$name} = $parameters->{$name};
			foreach ( @{$data->{$name}} ) {
			    $_->{'value'} =~ s/\%\%galacticusOutputPath\%\%/$galacticusOutputDirectory/g
				unless ( $name eq "label" );
			}
		    }
		    # Output the parameters as an XML file.
		    my $xmlOutput = new XML::Simple (RootName=>"parameters");
		    open(outHndl,">".$galacticusOutputDirectory."/parameters.xml");
		    print outHndl $xmlOutput->XMLout($data);
		    close(outHndl);
		    undef($parameters);
		    # Generate analysis code for this job.
		    my $analysisCode;
		    if ( $launchScript->{'doAnalysis'} eq "yes" ) {
			if ( $launchScript->{'analysisScript'} =~ m/\.xml$/ ) {
			    $analysisCode = $galacticusPath."scripts/analysis/Galacticus_Compute_Fit.pl ".$galacticusOutputDirectory."/galacticus.hdf5 ".$galacticusOutputDirectory." ".$launchScript->{'analysisScript'}."\n";
			} else {
			    $analysisCode = $launchScript->{'analysisScript'}." ".$galacticusOutputDirectory."\n";
			}
		    }
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

sub Create_Parameter_Hashes {
    # Create an array of hashes which give the parameter values for each model.
    # Get the input parameters structure.
    my $parameterSet = shift;
    # Convert to a more convenient hash structure.
    my $hash = &Parameters_To_Hash($parameterSet);
    # Populate an array of hashes with this initial hash.
    my @toProcessHashes = ( $hash );
    my @processedHashes;
    # Flatten the hash.
    while ( scalar(@toProcessHashes) > 0 ) {
	# Shift the first hash off the array.
	my $hash = shift(@toProcessHashes);
	# Record of whether this hash is flat. Assume it is initially.
	my $isFlat = 1;
	# Create a stack of parameters.
	my @stack = map {$hash->{$_}} keys(%{$hash});
	# Iterate over parameters.
	while ( scalar(@stack) > 0 ) {
	    # Pop a node array from the stack.
	    my $nodeArray = pop(@stack);
	    # Iterate over nodes.
	    foreach my $node ( @{$nodeArray} ) {		
		# Push any subparameters onto the stack.
		push(@stack,map {$node->{'subParameters'}->{$_}} keys(%{$node->{'subParameters'}}))
		    if ( exists($node->{'subParameters'}) );
		# Check for non-flat structure.
		if ( scalar(@{$node->{'values'}}) > 1 ) {
		    # Parameter has multiple values. Iterate over them, generating a new hash for each value, and
		    # push these new hashes back onto the stack.
		    my @values = @{$node->{'values'}};
		    foreach my $value ( @values ) {
			@{$node->{'values'}} = ( $value );
			my $newHash = clone($hash);
			push(
			    @toProcessHashes,
			    $newHash
			);
		    }
		    $isFlat = 0;
		} elsif ( exists(${$node->{'values'}}[0]->{'subtree'}) ) {
		    # Parameter has only one value, but it has sub-structure. Promote that substructure and push
		    # the hash back onto the stack.
		    foreach my $subName ( keys(%{${$node->{'values'}}[0]->{'subtree'}}) ) {		
			# Look for parameter definitions which modify existing values, and apply them.			
			my $subValuesCount = -1;
			foreach my $subValues ( @{${$node->{'values'}}[0]->{'subtree'}->{$subName}}  ) {
			    ++$subValuesCount;
			    foreach my $subValue ( @{$subValues->{'values'}} ) {
				if ( $subValue->{'value'} =~ m/\%\%modify\%\%([^\%]+)\%\%([^\%]+)\%\%/ ) {
				    my $find        = $1;
				    my $replaceBase = $2;
				    my @newValues = @{$hash->{$subName}->{'values'}};
				    foreach my $newValue ( @newValues ) {
					if ( my @captures = $newValue->{'value'} =~ m/$find/ ) {
					    my $replace = $replaceBase;
					    for(my $i=scalar(@captures)-1;$i>=0;--$i) {
						my $j = $i+1;
						$replace =~ s/\\$j/$captures[$i]/g;
					    }
					    $newValue->{'value'} =~ s/$find/$replace/;
					}
				    }
				    @{$subValues} = @newValues;
				}
			    }
			    if ( $subName eq "label" && exists($hash->{$subName}) ) {
				$hash->{$subName}->[$subValuesCount]->{'values'}->[0]->{'value'} .= "_".$subValues->{'values'}->[0]->{'value'};
			    } else {
				@{$hash->{$subName}->[$subValuesCount]->{'values'}} = @{$subValues->{'values'}};
			    }
			    $hash->{$subName}->[$subValuesCount]->{'subParameters'} = $subValues->{'subParameters'}
			        if ( exists($subValues->{'subParameters'}) );
			}
		    }
		    delete(${$node->{'values'}}[0]->{'subtree'});
		    push(
			@toProcessHashes,
			$hash
		    );
		    $isFlat = 0;	
		}
	    }
	    # Exit parameter iteration if the hash was found to be not flat.
	    last
		if ( $isFlat == 0 );
	}
	# If the hash is flat, then push it onto the processed hashes array.
	push(
	    @processedHashes,
	    $hash
	    ) if ( $isFlat == 1 );	
    }
    # Hashes are now flattened. Convert to simple form.
    foreach my $hash ( @processedHashes ) {
    	my @stack = map {$hash->{$_}} keys(%{$hash});	
    	while ( scalar(@stack) > 0 ) {
	    # Pop a node array from the stack.
	    my $nodeArray = pop(@stack);
	    # Iterate over nodes.
	    foreach my $node ( @{$nodeArray} ) {
		my $value = $node->{'values'}->[0]->{'value'};
		$value =~ s/^\s*//;
		$value =~ s/\s*$//;
		$node->{'value'} = $value;
		delete($node->{'values'});
		# Handle sub-parameters:
		#  --> Transfer them out of the "subParameters" element so that they will be in the correct location in the output XML;
		#  --> Push them onto the stack for conversion to simple form.
		if ( exists($node->{'subParameters'}) ) {
		    foreach ( keys(%{$node->{'subParameters'}}) ) {		    
			$node->{$_} = $node->{'subParameters'}->{$_};		   
			push(@stack,$node->{$_});
		    }
		    delete($node->{'subParameters'});
		}
	    }
    	}
    }
    # Return the result.
    return @processedHashes;
}

sub Parameters_To_Hash {
    # Convert an input parameter structure (as read from a Galacticus parameters XML file) into a more convenient internal hash.
    my $parameters = shift;
    my $hash;
    foreach my $name ( keys(%{$parameters}) ) {
	next
	    unless ( reftype($parameters->{$name}) );
	# Iterate over all subelements of this node.
	my $elementCounter = -1;
	foreach my $element ( &ExtraUtils::as_array($parameters->{$name}) ) {
	    my $subParameters;
	    ++$elementCounter;
	    foreach my $subElementName ( keys(%{$element}) ) { 
		if ( $subElementName eq "value" ) {
		    foreach my $value ( &ExtraUtils::as_array($element->{'value'}) ) {
			if ( UNIVERSAL::isa($value,"HASH") ) {
			    # This value of the parameter contains a subtree. Get a hash representation by calling ourself recursively.
			    push(
				@{$hash->{$name}->[$elementCounter]->{'values'}},
				{
				    value   => $value->{'content'},
				    subtree => &Parameters_To_Hash($value)
				}
				);
			} else {
			    # This value has no subtree, so just store the value.
			    push(
				@{$hash->{$name}->[$elementCounter]->{'values'}},
				{
				    value   => $value
				}
				);
			}
		    }
		} elsif ( $subElementName eq "modify" ) {
		    foreach my $modify ( &ExtraUtils::as_array($element->{'modify'}) ) {
			if ( exists($modify->{'parameter'}) ) {
			    # This modify of the parameter contains a subtree. Get a hash representation by calling ourself recursively.
			    push(
				@{$hash->{$name}->[$elementCounter]->{'values'}},
				{
				    modify  => $modify->{'content'},
				    subtree => &Parameters_To_Hash($modify)
				}
				);
			} else {
			    # This modify has no subtree, so just store the value.
			    push(
				@{$hash->{$name}->[$elementCounter]->{'values'}},
				{
				    value   => "\%\%modify\%\%".$modify->{'find'}."\%\%".$modify->{'replace'}."\%\%"
				}
				);
			}
		    }
		} else {
		    # For any other element, accumulate to set of subparameters.
		    $subParameters->{$subElementName} = $element->{$subElementName};
		}
	    }
	    # Process and attach any subparameters.
	    $hash->{$name}->[$elementCounter]->{'subParameters'} = &Parameters_To_Hash($subParameters)
		if ( $subParameters );	
 	}
    }
    return $hash;
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
	 verbosity          => 0                                                                  ,
	 md5Names           => "no"                                                               ,
	 useStateFile       => "no"                                                               ,
	 compressModels     => "no"                                                               ,
	 launchMethod       => "local"                                                            ,
	 modelRootDirectory => "./models"                                                         ,
	 baseParameters     => ""                                                                 ,
	 doAnalysis         => "no"                                                               ,
	 splitModels        => 1                                                                  ,
	 analysisScript     => $galacticusPath."data/analyses/Galacticus_Compute_Fit_Analyses.xml"
	);
    foreach ( keys(%defaults) ) {
	$launchScript->{$_} = $defaults{$_}
	unless ( exists($launchScript->{$_}) );
    }
    # Return the script.
    return $launchScript;
}

sub Parse_Galacticus_Config {
    # Parse any local configuration.
    my $config;
    if ( -e $galacticusPath."galacticusConfig.xml" ) {
	# Load XML.
	my $xml          = new XML::Simple;
	$config = $xml->XMLin($galacticusPath."galacticusConfig.xml");
    }
    # Return the config.
    return $config;
}
