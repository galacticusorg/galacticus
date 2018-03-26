# Contains a Perl module which implements various useful functionality for constructing parameter files for
# Galacticus when fitting to constraints.

package Galacticus::Constraints::Parameters;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use XML::LibXML;
use XML::Simple;
use XML::Twig;
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
use Clone qw(clone);
use List::Util;
use List::ExtraUtils;
use Galacticus::Launch::PBS;
use Scalar::Util 'reftype';
use Storable qw(dclone);

sub Parse_Config {
    # Get the config file name.
    my $configFile = shift;
    # Get any options.
    my %options;
    (%options) = @_
	if ( scalar(@_) > 0 );
    # Get the config.
    my $config;
    if ( exists($options{'useStored'}) && $options{'useStored'} == 1 && -e $configFile.".store" ) {
	$config = retrieve($configFile.".store");
    } else {
	# Parse the content.
	my $xml    = new XML::Simple;
	my $parser = XML::LibXML->new();
	my $dom    = $parser->load_xml(location => $configFile);
	$parser->process_xincludes($dom);
	my $twig = XML::Twig->new (comments => 'drop', pretty_print => 'indented');
	$twig->parse($dom->serialize());
	my $contentCommentless = $twig->sprint();
	$config = $xml->XMLin($contentCommentless, KeyAttr => 0);
	if ( UNIVERSAL::isa($config->{'parameters'},"ARRAY") ) {
	    my @parameters;
	    push(@parameters,&List::ExtraUtils::as_array($_->{'parameter'}))
		foreach (  &List::ExtraUtils::as_array($config->{'parameters'}) );
	    delete($config->{'parameters'});
	    @{$config->{'parameters'}->{'parameter'}} = @parameters;
	}
    }
    return $config;
}

sub Output {
    # Output a set of parameters to file.
    my $parameters = shift();
    my $fileName   = shift();
    # Test for parameters with no value element - ensure that they are forced to be arrays. This is done in a kludgey way by
    # inserting an empty "value" element. A better approach would be to avoid using XML::Simple, and instead use a better XML
    # parser which gives control over whether structure is serialized as elements or attributes.
    my $parametersCopy = dclone($parameters);
    my @nodes = map {$parametersCopy->{$_}} keys(%{$parametersCopy});
    while ( scalar(@nodes) > 0 ) {
    	my $node = pop(@nodes);
    	if ( UNIVERSAL::isa($node,"ARRAY") ) {
    	    unshift(@nodes,map {\$_} @{$node});
    	} elsif ( UNIVERSAL::isa($node,"HASH") ) {
    	    unless ( grep {$_ eq "value"} keys(%{$node}) ) {
		$node->{'value'} = "";
	    }
    		foreach ( keys(%{$node}) ) {
    		    unshift(@nodes,\${$node}{$_});
    	    }
    	}
    }
    # Create an XML output object.
    my $xmlOutput = new XML::Simple (RootName=>"parameters");
    # Output the parameters to file.
    open(my $outputFile,">".$fileName);
    print $outputFile $xmlOutput->XMLout($parametersCopy);
    close $outputFile;
}

sub Compilation {
    # Process a compilation of constraints, adjusting parameters as necessary and return a parameter hash.
    my $compilationFileName    = shift();
    my $baseParametersFileName = shift();
    my $adjustMasses           = shift();
    my $adjustOutputs          = shift();
    # Create an XML worker object.
    my $xml = new XML::Simple();
    # Retrieve the compilation file.
    my $compilationFilePath = ($compilationFileName =~ m/\//) ? $compilationFileName : "constraints/compilations/".$compilationFileName;
    my $compilation =
	-e $compilationFilePath.".store" 
	?
	retrieve($compilationFilePath.".store")
	:
	$xml->XMLin($compilationFilePath, KeyAttr => "");
    # Parse the base set of parameters.
    my $parameters = 
	-e $baseParametersFileName.".store" 
	?
	retrieve($baseParametersFileName.".store") 
	:
	$xml->XMLin($baseParametersFileName);
    # Specify default values for parameters which can be adjusted by constraint definitions.
    my %outputLuminosities;
    my %optionsOn;
    my %outputRedshifts;
    if ( exists($parameters->{'outputRedshifts'}) ) {
	foreach my $redshift ( split(" ",$parameters->{'outputRedshifts'}->{'value'}) ) {
	    $outputRedshifts{&redshiftPrecision($redshift)} = 1;
	}
    }
    my $haloMassResolution;
    if ( exists($parameters->{'mergerTreeMassResolutionMethod'}) ) {
	if ( $parameters->{'mergerTreeMassResolutionMethod'}->{'value'} eq "fixed" ) {
	    $haloMassResolution = $parameters->{'mergerTreeMassResolutionMethod'}->{'massResolution'       }->{'value'};
	} elsif ( $parameters->{'mergerTreeMassResolutionMethod'}->{'value'} eq "scaled" ) {
	    $haloMassResolution = $parameters->{'mergerTreeMassResolutionMethod'}->{'massResolutionMinimum'}->{'value'};
	}
    }
    my $haloMassMinimum = 
	exists($parameters->{'mergerTreeBuildHaloMassMinimum'}) 
	?
	$parameters->{'mergerTreeBuildHaloMassMinimum'}->{'value'} 
        :
	undef();
    my $haloMassMaximum = 
	exists($parameters->{'mergerTreeBuildHaloMassMaximum'}) 
	?
	$parameters->{'mergerTreeBuildHaloMassMaximum'}->{'value'} 
        :
	undef();
    # Scan through all constraints.
    my @constraints = exists($compilation->{'constraint'}) ? &List::ExtraUtils::as_array($compilation->{'constraint'}) : ();
    foreach my $constraint ( @constraints ) {
	# Check we have a definition file.
	die("Compilation(): compilation must specify a definition for each constraint")
	    unless ( defined($constraint->{'definition'}) );
	# Parse the definition file.
	my $constraintDefinition =
	    -e $constraint->{'definition'}.".store"
	    ?
	    retrieve($constraint->{'definition'}.".store") 
	    :
	    $xml->XMLin($constraint->{'definition'},KeyAttr => "");
	# Extract any required output redshift.
	if ( defined($constraintDefinition->{'outputRedshift'}) ) {
	    foreach my $redshift ( &List::ExtraUtils::as_array($constraintDefinition->{'outputRedshift'}) ) {
		# Add each redshift (at a specific precision) to the list of output redshifts.
		$outputRedshifts{&redshiftPrecision($redshift)} = 1;
	    }
	}
	# Extract any filter definitions.
	if ( defined($constraintDefinition->{'luminosity'}) ) {
	    foreach my $luminosity ( &List::ExtraUtils::as_array($constraintDefinition->{'luminosity'}) ) {
		# Encode the luminosity into a unique key.
		my $key = join(":",($luminosity->{'filter'},sprintf("%6.4f",$luminosity->{'redshift'}),$luminosity->{'frame'}));
		$outputLuminosities{$key} = 1;
	    }
	}
	# Ensure that the halo mass resolution is sufficient.
	if ( defined($constraintDefinition->{'haloMassResolution'}) ) {
	    $haloMassResolution = $constraintDefinition->{'haloMassResolution'}
	        unless ( defined($haloMassResolution) && $constraintDefinition->{'haloMassResolution'} > $haloMassResolution );
	}
	# Ensure that the minimum halo mass is sufficient.
	if ( defined($constraintDefinition->{'haloMassMinimum'}) ) {
	    $haloMassMinimum = $constraintDefinition->{'haloMassMinimum'}
	        unless ( defined($haloMassMinimum   ) && $constraintDefinition->{'haloMassMinimum'   } > $haloMassMinimum    );
	}
	# Ensure that the maximum halo mass is sufficient.
	if ( defined($constraintDefinition->{'haloMassMaximum'}) ) {
	    $haloMassMaximum = $constraintDefinition->{'haloMassMaximum'}
	        unless ( defined($haloMassMaximum   ) && $constraintDefinition->{'haloMassMaximum'   } < $haloMassMaximum    );
	}
	# Accumulate any options that must be switched on.
	if ( defined($constraintDefinition->{'optionOn'}) ) {
	    for my $option ( &List::ExtraUtils::as_array($constraintDefinition->{'optionOn'}) ) {$optionsOn{$option} = 1};
	}    
	# Accumulate any parameters that must be set.
	if ( defined($constraintDefinition->{'parameter'}) ) {
	    for my $parameter ( &List::ExtraUtils::as_array($constraintDefinition->{'parameter'}) ) {
		my $name  = $parameter->{'name' };
		my $value = $parameter->{'value'};
		$value =~ s/^\s*(.*?)\s*$/$1/;
		my @values     = split(/\s+/,$value);
		my @currentValues;
		@currentValues = split(/\s+/,$parameters->{$name}->{'value'}) if ( exists($parameters->{$name}) );
		my $accumulation = "overwrite";
		$accumulation = $parameter->{'accumulation'} if ( exists($parameter->{'accumulation'}) );
		if ( $accumulation eq "overwrite" ) {
		    @currentValues = @values;
		} elsif ( $accumulation eq "combine" ) {
		    push(@currentValues,@values);
		} elsif ( $accumulation eq "unique" ) {
		    push(@currentValues,@values);
		    @currentValues = uniq(@currentValues);
		}
		$parameters->{$name}->{'value'} = join(" ",@currentValues);
	    }
	}
    }
    # Construct a sorted list of redshifts, suitable for output as a parameter. Ensure that we always have at least one output
    # redshift.
    if ( $adjustOutputs eq "T" ) {
	my @outputRedshiftList = keys(%outputRedshifts);
	if ( scalar(@outputRedshiftList) == 0 ) {
	    push(@outputRedshiftList,&redshiftPrecision(0.0));
	} else {
	    @outputRedshiftList = sort(@outputRedshiftList);
	}
	# Modify the set of output redshifts.
	$parameters->{'outputRedshifts'}->{'value'} = join(" ",@outputRedshiftList);
    }
    # Modify the minimum and maximum halo masses.
    $haloMassResolution = 5.00e09 unless ( defined($haloMassResolution) );
    $haloMassMinimum    = 1.00e10 unless ( defined($haloMassMinimum   ) );
    $haloMassMaximum    = 1.01e10 unless ( defined($haloMassMaximum   ) );
    if ( $adjustMasses eq "T" ) {
	$parameters->{'mergerTreeBuildHaloMassMinimum'}->{'value'} = $haloMassMinimum;
	$parameters->{'mergerTreeBuildHaloMassMaximum'}->{'value'} = $haloMassMaximum;
	if ( exists($parameters->{'mergerTreeMassResolutionMethod'}) ) {
	    if ( $parameters->{'mergerTreeMassResolutionMethod'}->{'value'} eq "fixed" ) {
		$parameters->{'mergerTreeMassResolutionMethod'}->{'massResolution'       }->{'value'} = $haloMassResolution;
	    } elsif ( $parameters->{'mergerTreeMassResolutionMethod'}->{'value'} eq "scaled" ) {
		$parameters->{'mergerTreeMassResolutionMethod'}->{'massResolutionMinimum'}->{'value'} = $haloMassResolution;
	    }
	} else {
	    $parameters->{'mergerTreeMassResolutionMethod'}                    ->{'value'} = "fixed";
	    $parameters->{'mergerTreeMassResolutionMethod'}->{'massResolution'}->{'value'} = $haloMassResolution;
	}
    }
    # Set required options on.
    $parameters->{$_}->{'value'} = "true" 
	foreach ( keys(%optionsOn) );
    # Construct luminosity requirements.
    $parameters->{'luminosityFilter'  }->{'value'} = join(" ",map {(split(/:/,$_))[0]} keys(%outputLuminosities))
	if ( scalar(keys(%outputLuminosities)) > 0 );
    $parameters->{'luminosityRedshift'}->{'value'} = join(" ",map {(split(/:/,$_))[1]} keys(%outputLuminosities))
	if ( scalar(keys(%outputLuminosities)) > 0 );
    $parameters->{'luminosityType'    }->{'value'} = join(" ",map {(split(/:/,$_))[2]} keys(%outputLuminosities))
	if ( scalar(keys(%outputLuminosities)) > 0 );
    # Return the parameter hash.
    return (\@constraints,$parameters);
}

sub Active_Parameter_Names {
    # Return a list of active parameter names.
    my $config = shift;
    my @names;
    # Extract parameters from config file.
    my @parameters = &List::ExtraUtils::as_array($config->{'parameters'}->{'parameter'});
    # Extract names.
    my $parameterCount = 0;
    for(my $i=0;$i<scalar(@parameters);++$i) {
	push(@names,$parameters[$i]->{'name'})
	    if ( exists($parameters[$i]->{'prior'}) );
    }
    return @names;
}

sub Convert_Parameters_To_Galacticus {
    my @parameterDefinitions = &List::ExtraUtils::as_array(shift());
    # Build an array of parameter names and values.
    my @parameters;
    my %parameterValues;
    foreach my $parameterDefinition ( @parameterDefinitions ) {
	my $parameter;
	if      ( $parameterDefinition =~ m/^([^=]+)==([^=]+)$/ ) {
	    $parameter->{'name'      } = $1;
	    $parameter->{'definition'} = $2;
	} elsif ( $parameterDefinition =~ m/([^=]+)=([^=]+)$/  ) {
	    $parameter->{'name'      } = $1;
	    $parameterValues{$parameter->{'name'}} = $2;
	} else {
	    die('unable to parse parameter definition');
	}
	push(@parameters,$parameter);
    }
    # Set the values of any parameters that are defined in terms of other parameters.
    my $failCount  = 1;
    my $iterations = 0;
    while ( $failCount > 0 ) {
	$failCount = 0;
	++$iterations;
	die("Convert_Parameters_To_Galacticus: Failed to resolve parameter definitions")
	    if ( $iterations > 100000 );
	for(my $i=0;$i<scalar(@parameters);++$i) {
	    if ( exists($parameters[$i]->{'definition'}) ) {
		# Attempt to replace named parameters in the definition with their values.
		while ( $parameters[$i]->{'definition'} =~ m/\%\[([a-zA-Z0-9_\.:]+)\]/ ) {
		    my $parameterName = $1;
		    if ( exists($parameterValues{$parameterName}) ) {
			$parameters[$i]->{'definition'} =~ s/\%\[$parameterName\]/$parameterValues{$parameterName}/g;
		    } else {
			++$failCount;
			last;
		    }
		}
		$parameterValues{$parameters[$i]->{'name'}} = eval($parameters[$i]->{'definition'})
		    unless ( $parameters[$i]->{'definition'} =~ m/\%\[([a-zA-Z0-9_\.:]+)\]/ );
	    }
	}
    }
    # Create a hash of new parameters.
    my $newParameters;
    for(my $i=0;$i<scalar(@parameters);++$i) {
	$newParameters->{$parameters[$i]->{'name'}} = $parameterValues{$parameters[$i]->{'name'}};
    }
    # Return the parameter has.
    return $newParameters;
}

sub Chains_Count {
    # Return a count of the number of chains used.
    my $config    =   shift() ;
    my %arguments = %{shift()};
    # Find the work directory.
    my $workDirectory = $config->{'likelihood'}->{'workDirectory'};
    # Determine number of chains.
    my $logFileRoot = $config->{'simulation'}->{'logFileRoot'};
    my $chainCount  = 0;
    while () {
	++$chainCount;
	my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,$chainCount);
	last
	unless ( -e $chainFileName );
    }
    return $chainCount;
}

sub Parameters_Count {
    # Return a count of the number of parameters used.
    my $config    =   shift() ;
    my %arguments = %{shift()};
    # Find the work directory.
    my $workDirectory = $config->{'likelihood'}->{'workDirectory'};
    # Read the first chain.
    my $logFileRoot = $config->{'simulation'}->{'logFileRoot'};
    my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,0);
    my $parameterCount;
    open(iHndl,$chainFileName);
    while ( my $line = <iHndl> ) {
	unless ( $line =~ m/^\"/ ) {
	    my @columns = split(" ",$line);
	    $parameterCount = scalar(@columns)-6;
	}
    }
    die("Galacticus::Constraints::Parameters::Parameters_Count(): could not determine number of parameters")
	unless ( defined($parameterCount) );
    return $parameterCount;
}

sub Sample_Matrix {
    # Generate a matrix of parameters sampled from the posterior distribution.
    my $config    =   shift() ;
    my %arguments = %{shift()};
    # Find the work directory.
    my $workDirectory = $config->{'likelihood'}->{'workDirectory'};
    # Determine number of chains.
    my $logFileRoot = $config->{'simulation'}->{'logFileRoot'};    
    my $chainCount  = &Chains_Count($config,\%arguments);
    # Build a list of outlier chains.
    my @outlierChains = exists($arguments{'outliers'}) ? split(/,/,$arguments{'outliers'}) : ();
    # Parse the chains to find all viable parameter sets.
    my $iSample = -1;
    my $viable = pdl [];
    my $parameterCount;
    for(my $i=0;$i<$chainCount;++$i) {
	# Skip outlier chains.
	next
	    if ( grep {$_ eq $i} @outlierChains );
	# Skip non-selected chain.
	next
	    if ( defined($arguments{'selectChain'}) && $i != $arguments{'selectChain'} );
	# Parse the chain file.
	my @suffixes = ( "" );
	push(@suffixes,"Previous")
	    if ( exists($arguments{'includePrevious'}) && $arguments{'includePrevious'} eq "yes" );
	foreach my $suffix ( @suffixes ) {
	    my $chainFileName = sprintf("%s%s_%4.4i.log",$logFileRoot,$suffix,$i);
	    next
		unless ( -e $chainFileName );
	    open(iHndl,$chainFileName);
	    while ( my $line = <iHndl> ) {
		unless ( $line =~ m/^\"/ ) {
		    ++$iSample;
		    my @columns = split(" ",$line);
		    my $accept = 1;
		    # Skip unconverged states unless explicitly allowed.
		    $accept = 0
			if ( $columns[3] eq "F" && ( ! exists($arguments{'useUnconverged'}) || $arguments{'useUnconverged'} eq "no" ) );
		    # Skip chains before the given start point.
		    $accept = 0
			if ( exists($arguments{'sampleFrom'}) && $columns[0] < $arguments{'sampleFrom'} );
		    $viable = $viable->append($iSample)
			if ( $accept == 1 );
		    $parameterCount = scalar(@columns)-6;
		}
	    } 
	    close(iHndl);
	}
    }
    die("Galacticus::Constraints::Parameters::Sample_Matrix(): no viable samples found")
	if ( nelem($viable) == 0 );
    # Sample parameters.
    my $randomSample = exists($arguments{'sampleCount'}) && $arguments{'sampleCount'} > 0 ? 1 : 0;
    $arguments{'sampleCount'} = nelem($viable)
	if ( ! exists($arguments{'sampleCount'}) || $arguments{'sampleCount'} <= 0 );
    my $sampleIndex;
    if ( $randomSample == 1 ) {
	$sampleIndex = pdl long(nelem($viable)*random($arguments{'sampleCount'}));
    } else {
	$sampleIndex = pdl sequence(nelem($viable));
    }
    $sampleIndex .= $sampleIndex->qsort();
    # Build the matrix.
    my $sampleMatrix = pdl zeroes($parameterCount,nelem($sampleIndex));
    $iSample    = -1;
    my $jSample =  0;
    for(my $i=0;$i<$chainCount;++$i) {
	# Skip outlier chains.
	next
	    if ( grep {$_ eq $i} @outlierChains );
	# Skip non-selected chain.
	next
	    if ( defined($arguments{'selectChain'}) && $i != $arguments{'selectChain'} );
	# Parse the chain file.
	my @suffixes = ( "" );
	push(@suffixes,"Previous")
	    if ( exists($arguments{'includePrevious'}) && $arguments{'includePrevious'} eq "yes" );
	foreach my $suffix ( @suffixes ) {
	    my $chainFileName = sprintf("%s%s_%4.4i.log",$logFileRoot,$suffix,$i);
	    next
		unless ( -e $chainFileName );
	    open(iHndl,$chainFileName);
	    while ( my $line = <iHndl> ) {
		unless ( $line =~ m/^\"/ ) {
		    last
			if ( $jSample >= nelem($sampleIndex) );
		    ++$iSample;
		    my @columns = split(" ",$line);
		    my $accept = 1;
		    # Skip unconverged states unless explicitly allowed.
		    $accept = 0
			if ( $columns[3] eq "F" && ( ! exists($arguments{'useUnconverged'}) || $arguments{'useUnconverged'} eq "no" ) );
		    # Skip chains before the given start point.
		    $accept = 0
			if ( exists($arguments{'sampleFrom'}) && $columns[0] < $arguments{'sampleFrom'} );
		    if ( $accept && $iSample == $viable->($sampleIndex)->(($jSample)) ) {
			$sampleMatrix->(:,($jSample)) .= pdl @columns[6..$#columns];
			++$jSample;			
		    }
		}
	    } 
	    close(iHndl);
	}
    }
    # If parameters are to be mapped, do so.
    if ( exists($arguments{'parametersMapped'}) && $arguments{'parametersMapped'} eq "yes" ) {
	my $i = -1;
	foreach my $parameter ( @{$config->{'parameters'}->{'parameter'}} ) {
	    if ( exists($parameter->{'prior'}) ) {
		++$i;
		if ( $parameter->{'mapping'}->{'type'} eq "logarithmic" ) {
		    $sampleMatrix->(($i),:) .= log($sampleMatrix->(($i),:));		    
		}
	    }
	}
    }
    # Return the matrix.
    return $sampleMatrix;
}

sub Sample_Models {
    # Generate a sample of models from the posterior distribution.
    my $config    =   shift() ;
    my %arguments = %{shift()};
    (my %options) = @_
	if ( scalar(@_) > 0 );
    # Find the work directory.
    my $workDirectory = $config->{'likelihood'}->{'workDirectory'};
    # Get a hash of the parameter values.
    my $compilationFile;
    if ( exists($arguments{'compilationOverride'}) ) {
	$compilationFile = $arguments{'compilationOverride'};
    } else {
	$compilationFile = $config->{'likelihood'}->{'compilation'};
    }
    (my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($compilationFile,$config->{'likelihood'}->{'baseParameters'});
    my @constraints = @{$constraintsRef};
    # Get a matrix of sampled states.
    my $sampleMatrix = &Sample_Matrix($config,\%arguments);
    # Run model for each sample.
    my $sampleDirectory = exists($arguments{'sampleDirectory'}) ? $arguments{'sampleDirectory'}."/" : $workDirectory."/posteriorSample/";
    my @pbsStack;
    for (my $i=0;$i<$sampleMatrix->dim(1);++$i) {
	# Create an output directory.
	my $modelDirectory = $sampleDirectory.$i."/";
	system("mkdir -p ".$modelDirectory);
	my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
	# Check if the model has already been run.
	unless ( -e $galacticusFileName ) {
	    # Convert these values into a parameter array.
	    my $currentConfig = clone($config);
	    my $newParameters = &Convert_Parameters_To_Galacticus($currentConfig, $sampleMatrix->(:,($i))->list());
	    # Increment the random number seed.
	    $parameters->{'randomSeed'}->{'value'} += $config->{'likelihood'}->{'threads'};
	    # Clone parameters.
	    my $currentParameters = clone($parameters);
	    # Apply to parameters.
	    &Apply_Parameters($currentParameters,$newParameters);
	    # Apply any parameters from command line.
	    &Apply_Command_Line_Parameters($currentParameters,\%arguments);
	    # Apply any other parameter modifications requested by the caller.
	    &{$options{'parametersModifier'}}($currentParameters)
		if ( exists($options{'parametersModifier'}) );
	    # Specify the output file name.
	    $currentParameters->{'galacticusOutputFileName'}->{'value'} = $galacticusFileName;
	    # Write the modified parameters to file.
	    &Output($currentParameters,$modelDirectory."parameters.xml");
	    # Construct the tasks to perform.
	    my $command;
	    if ( exists($config->{'likelihood'}->{'environment'}) ) {
		foreach ( &List::ExtraUtils::as_array($config->{'likelihood'}->{'environment'}) ) {
		    $command .= "export ".$_."\n";
		}
	    }
	    $command .= "ulimit -t unlimited\n";
	    $command .= "ulimit -c unlimited\n";
	    $command .= "export OMP_NUM_THREADS=".$config->{'likelihood'}->{'threads'}."\n";
	    $command .= "mpirun --bynode -np 1 ".(exists($config->{'likelihood'}->{'executable'}) ? $config->{'likelihood'}->{'executable'} : "Galacticus.exe")." ".$modelDirectory."parameters.xml\n";
	    unless ( exists($arguments{'analyze'}) && $arguments{'analyze'} eq "no" ) {
		foreach my $constraint ( @constraints ) {
		    # Parse the definition file.
		    my $xml = new XML::Simple;
		    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});	    
		    # Insert code to run the analysis code.
		    my $analysisCode = $constraintDefinition->{'analysis'};
		    (my $plotLabel = $constraintDefinition->{'label'}) =~ s/\./_/g;
		    $command .= $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5 --plotFile ".$modelDirectory."/".$plotLabel.".pdf --outputFile ".$modelDirectory."/".$constraintDefinition->{'label'}."Likelihood.xml --modelDiscrepancies ".$workDirectory."/modelDiscrepancy\n";
		}
	    }
	    # Create a PBS job.
	    my %job =
		(
		 launchFile => $modelDirectory."/launch.pbs",
		 label      => $config->{'likelihood'}->{'name'}."_ppc".$i,
		 logFile    => $modelDirectory."/launch.log",
		 command    => $command
		);
	    foreach ( 'ppn', 'walltime', 'memory' ) {
		$job{$_} = $arguments{$_}
		if ( exists($arguments{$_}) );
	    }
	    # Queue the calculation.
	    push(@pbsStack,\%job);
	}
    }
    # Send jobs to PBS.
    &Galacticus::Launch::PBS::SubmitJobs(\%arguments,@pbsStack)
     	if ( scalar(@pbsStack) > 0 );
    # Return the number of models sampled.
    return ($sampleMatrix->dim(1), $sampleDirectory);
}

sub Maximum_Likelihood_Vector {
    my $config    =   shift() ;
    my %arguments = %{shift()};
    # Determine the MCMC directory.
    my $logFileRoot = $config->{'simulation'}->{'logFileRoot'};
    (my $mcmcDirectory  = $logFileRoot) =~ s/\/[^\/]+$//;    
    # Determine number of chains.
    my $chainCount = &Chains_Count($config,\%arguments);
    # Parse the chains to find the maximum likelihood model.
    my $maximumLikelihood = -1e30;
    my @maximumLikelihoodParameters;
    my @chainFiles;
    for(my $i=0;$i<$chainCount;++$i) {
	next
	    unless
	    (
	     ! exists($arguments{'chain'})
	     ||
	     $arguments{'chain'} eq "all"
	     ||
	     $arguments{'chain'} == $i
	    );
	my $chainFileName         = sprintf("%s_%4.4i.log"        ,$logFileRoot,$i);
	my $chainFilePreviousName = sprintf("%sPrevious_%4.4i.log",$logFileRoot,$i);
	push(@chainFiles,$chainFileName        );
	push(@chainFiles,$chainFilePreviousName)
	    if ( exists($arguments{'includePrevious'}) && $arguments{'includePrevious'} eq "yes" && -e $chainFilePreviousName );
    }
    foreach my $chainFile ( @chainFiles ) {
	open(iHndl,$chainFile);
	while ( my $line = <iHndl> ) {
	    unless ( $line =~ m/^\"/ ) {
		$line =~ s/^\s*//;
		$line =~ s/\s*$//;
		my @columns = split(/\s+/,$line);
		# Determine if state is accepted.
		my $accept = 1;
		# Skip unconverged states unless explicitly allowed.
		$accept = 0
		    if ( $columns[3] eq "F" && ( ! exists($arguments{'useUnconverged'}) || $arguments{'useUnconverged'} eq "no" ) );
		if ( $accept == 1 && $columns[4] > $maximumLikelihood ) {
		    $maximumLikelihood           = $columns[4];
		    @maximumLikelihoodParameters = @columns[6..$#columns];
		}
	    }
	}
	close(iHndl);
    }
    # Convert parameters to a PDL.
    my $maximumLikelihoodVector = pdl @maximumLikelihoodParameters;
    # Return the vector.
    return $maximumLikelihoodVector;
}

sub Maximum_Likelihood_Parameters {
    my $config    =   shift() ;
    my %arguments = %{shift()};
    # Get a vector of maximum likelihood parameters.
    my $maximumLikelihoodVector = &Maximum_Likelihood_Vector($config,\%arguments);	
    # Get the parameters.
    my $parameters = &Convert_Parameter_Vector_To_Galacticus($config,$maximumLikelihoodVector);
    return $parameters;
}

sub Convert_Parameter_Vector_To_Galacticus {
    # Convert a PDL vector of parameter values into a Galacticus parameter structure.
    my $config          = shift();
    my $parameterVector = shift();
    # Get a hash of the parameter values.
    (my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($config->{'likelihood'}->{'compilation'},$config->{'likelihood'}->{'baseParameters'});
    # Convert vector into a parameter array.
    my $newParameters = &Galacticus::Constraints::Parameters::Convert_Parameters_To_Galacticus($config,$parameterVector->list());
    # Apply to parameters.
    &Apply_Parameters($parameters,$newParameters);
    return $parameters;
}

sub Apply_Parameters {
    # Apply a set of new parameters from a constraint procedure to an existing parameters data structure.
    my $parameters    = shift();
    my $newParameters = shift();
    for my $newParameterName ( keys(%{$newParameters}) ) {
	my $parameter = $parameters;
	my $valueIndex;
	if ( $newParameterName =~ m/^(.*)\{(\d+)\}$/ ) {
	    $newParameterName = $1;
	    $valueIndex = $2;
	}
	foreach ( split(/\-\>/,$newParameterName) ) {
	    # Check if the parameter name contains an array reference.
	    if ( $_ =~ m/^(.*)\[(\d+)\]$/ ) {
		# Parameter name contains array reference. Step through to the relevant parameter in the list. If the parameter is
		# not an array, allow this only if the array index given is zero.
		if ( reftype($parameter->{$1}) eq "ARRAY" ) {
		    $parameter->{$1}->[$2]->{'value'} = undef()
			unless ( scalar(@{$parameter->{$1}}) > $2 );
		    $parameter = $parameter->{$1}->[$2];
		} else {
		    die('constrainGalacticus.pl: attempt to access non-existant array')
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
}

sub Apply_Command_Line_Parameters {
    my $parameters =   shift() ;
    my %arguments  = %{shift()};
    # Apply any parameters from command line.
    foreach my $argument ( keys(%arguments) ) {
	if ( $argument =~ m/^parameter:(.*)/ ) {
	    my @parametersSet = split(":",$1);
	    my $parameter     = $parameters;
	    my $parameterPrevious;
	    my $parameterNameCurrent;
	    foreach my $parameterSet ( @parametersSet ) {
		if ( $parameterSet =~ m/([^\{]*)(\{{0,1}([^\}]*)\}{0,1})/ ) {
		    my $parameterName  = $1;
		    my $parameterValue = $3;
		    $parameter->{$parameterName}->{'value'} = $parameterValue
			if ( defined($2) );
		    $parameterNameCurrent = $parameterName;
		    $parameterPrevious    = $parameter;
		    $parameter            = $parameter->{$parameterName};
		} else {
		    die("malformed parameter definition");
		}
	    }
	    if ( $arguments{$argument} eq "DELETE" ) {
		delete($parameterPrevious->{$parameterNameCurrent});
	    } else {
		$parameter->{'value'} = $arguments{$argument};
	    }
	}
    }
}

sub redshiftPrecision {
    # Return a redshift formatted to given precision.
    my $redshift = shift();
    return sprintf("%6.4f",$redshift);
}

1;
