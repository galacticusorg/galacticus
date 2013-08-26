# Contains a Perl module which implements various useful functionality for constructing parameter files for
# Galacticus when fitting to constraints.

package Parameters;
use strict;
use warnings;
use XML::Simple;
use XML::SAX;
use XML::SAX::Writer;
use XML::Filter::XInclude;
use XML::Twig;
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
use Clone qw(clone);
use List::Util;

sub Parse_Config {
    # Get the config file name.
    my $configFile = shift;
    # Parse the content.
    my $xml    = new XML::Simple;
    my $content;
    my $parser = XML::SAX::ParserFactory->parser(
	Handler => XML::Filter::XInclude->new(
	    Handler => XML::SAX::Writer->new(Output => \$content)
	)
	);
    $parser->parse_uri($configFile);
    my $twig = XML::Twig->new (comments => 'drop', pretty_print => 'indented');
    $twig->parse($content);
    my $contentCommentless = $twig->sprint();
    my $config = $xml->XMLin($contentCommentless, KeyAttr => 0);
    if ( UNIVERSAL::isa($config->{'parameters'},"ARRAY") ) {
	my @parameters;
	push(@parameters,@{$_->{'parameter'}})
	    foreach ( @{$config->{'parameters'}} );
	delete($config->{'parameters'});
	@{$config->{'parameters'}->{'parameter'}} = @parameters;
    }
    return $config;
}

sub Output {
    # Output a set of parameters to file.
    my $parameters = shift;
    my $fileName   = shift;
    # Create an XML output object.
    my $xmlOutput  = new XML::Simple (RootName=>"parameters", NoAttr => 1);
    # Transfer parameters to a suitable array structure.
    my $outputParameters;
    push(@{$outputParameters->{'parameter'}},{name => $_, value => $parameters->{'parameter'}->{$_}->{'value'}})
	foreach ( keys(%{$parameters->{'parameter'}}) );
    # Output the parameters to file.
    open(pHndl,">".$fileName);
    print pHndl $xmlOutput->XMLout($outputParameters);
    close pHndl;
}

sub Compilation {
    # Process a compilation of constraints, adjusting parameters as necessary and return a parameter hash.
    my $compilationFileName    = shift;
    my $baseParametersFileName = shift;
    # Create an XML worker object.
    my $xml = new XML::Simple;
    # Parse the constraint compilation file.
    my $compilation = $xml->XMLin("constraints/compilations/".$compilationFileName);
    # Specify default values for parameters which can be adjusted by constraint definitions.
    my %outputRedshifts;
    my %outputLuminosities;
    my %optionsOn;
    my $haloMassResolution;
    my $haloMassMinimum;
    my $haloMassMaximum;
    # Parse the base set of parameters.
    my $parameters = $xml->XMLin($baseParametersFileName);
    # Scan through all constraints.
    my @constraints;
    if ( ref($compilation->{'constraint'}) eq "ARRAY" ) {
	@constraints = @{$compilation->{'constraint'}};
    } else {
	push(@constraints,$compilation->{'constraint'});
    }
    foreach my $constraint ( @constraints ) {
	# Check we have a definition file.
	die("Compilation(): compilation must specify a definition for each constraint")
	    unless ( defined($constraint->{'definition'}) );
	# Parse the definition file.
	my $constraintDefinition = $xml->XMLin($constraint->{'definition'},KeyAttr => "");
	# Extract any required output redshift.
	if ( defined($constraintDefinition->{'outputRedshift'}) ) {
	    my @redshifts;
	    if ( ref($constraintDefinition->{'outputRedshift'}) eq "ARRAY" ) {
		@redshifts = @{$constraintDefinition->{'outputRedshift'}};
	    } else {
		push(@redshifts,$constraintDefinition->{'outputRedshift'});
	    }
	    foreach my $redshift ( @redshifts ) {
		# Add each redshift (at a specific precision) to the list of output redshifts.
		my $precisionRedshift = sprintf("%6.4f",$redshift);
		$outputRedshifts{$precisionRedshift} = 1;
	    }
	}
	# Extract any filter definitions.
	if ( defined($constraintDefinition->{'luminosity'}) ) {
	    my @luminosities;
	    if ( ref($constraintDefinition->{'luminosity'}) eq "ARRAY" ) {
		@luminosities = @{$constraintDefinition->{'luminosity'}};
	    } else {
		push(@luminosities,$constraintDefinition->{'luminosity'});
	    }
	    foreach my $luminosity ( @luminosities ) {
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
	    unless ( defined($haloMassMinimum) && $constraintDefinition->{'haloMassMinimum'} > $haloMassMinimum );
	}
	# Ensure that the maximum halo mass is sufficient.
	if ( defined($constraintDefinition->{'haloMassMaximum'}) ) {
	    $haloMassMaximum = $constraintDefinition->{'haloMassMaximum'}
	    unless ( defined($haloMassMaximum) && $constraintDefinition->{'haloMassMaximum'} < $haloMassMaximum );
	}
	# Accumulate any options that must be switched on.
	if ( defined($constraintDefinition->{'optionOn'}) ) {
	    my @options;
	    if ( ref($constraintDefinition->{'optionOn'}) eq "ARRAY" ) {
		@options = @{$constraintDefinition->{'optionOn'}};
	    } else {
		push(@options,$constraintDefinition->{'optionOn'});
	    }
	    for my $option ( @options ) {$optionsOn{$option} = 1};
	}    
	# Accumulate any parameters that must be set.
	if ( defined($constraintDefinition->{'parameter'}) ) {
	    my @extraParameters;
	    if ( ref($constraintDefinition->{'parameter'}) eq "ARRAY" ) {
		@extraParameters = @{$constraintDefinition->{'parameter'}};
	    } else {
		push(@extraParameters,$constraintDefinition->{'parameter'});
	    }
	    for my $parameter ( @extraParameters ) {
		my $name  = $parameter->{'name' };
		my $value = $parameter->{'value'};
		$value =~ s/^\s*(.*?)\s*$/$1/;
		my @values     = split(/\s+/,$value);
		my @currentValues;
		@currentValues = split(/\s+/,$parameters->{'parameter'}->{$name}->{'value'}) if ( exists($parameters->{'parameter'}->{$name}) );
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
		$parameters->{'parameter'}->{$name}->{'value'} = join(" ",@currentValues);
	    }
	}
    }
    # Construct a sorted list of redshifts, suitable for output as a parameter. Ensure that we always have at least one output
    # redshift.
    my @outputRedshiftList = keys(%outputRedshifts);
    if ( scalar(@outputRedshiftList) == 0 ) {
	push(@outputRedshiftList,"0.0000");
    } else {
	@outputRedshiftList = sort(@outputRedshiftList);
    }
    # Modify the set of output redshifts.
    $parameters->{'parameter'}->{'outputRedshifts'}->{'value'} = join(" ",@outputRedshiftList);
    # Modify the minimum and maximum halo masses.
    $haloMassResolution = 5.00e09 unless ( defined($haloMassResolution) );
    $haloMassMinimum    = 1.00e10 unless ( defined($haloMassMinimum   ) );
    $haloMassMaximum    = 1.01e10 unless ( defined($haloMassMaximum   ) );
    $parameters->{'parameter'}->{'mergerTreeBuildMassResolutionFixed'}->{'value'} = $haloMassResolution;
    $parameters->{'parameter'}->{'mergerTreeBuildHaloMassMinimum'    }->{'value'} = $haloMassMinimum   ;
    $parameters->{'parameter'}->{'mergerTreeBuildHaloMassMaximum'    }->{'value'} = $haloMassMaximum   ;
    # Set required options on.
    $parameters->{'parameter'}->{$_}->{'value'} = "true" 
	foreach ( keys(%optionsOn) );
    # Construct luminosity requirements.
    $parameters->{'parameter'}->{'luminosityFilter'  }->{'value'} = join(" ",map {(split(/:/,$_))[0]} keys(%outputLuminosities))
	if ( scalar(keys(%outputLuminosities)) > 0 );
    $parameters->{'parameter'}->{'luminosityRedshift'}->{'value'} = join(" ",map {(split(/:/,$_))[1]} keys(%outputLuminosities))
	if ( scalar(keys(%outputLuminosities)) > 0 );
    $parameters->{'parameter'}->{'luminosityType'    }->{'value'} = join(" ",map {(split(/:/,$_))[2]} keys(%outputLuminosities))
	if ( scalar(keys(%outputLuminosities)) > 0 );
    # Return the parameter hash.
    return (\@constraints,$parameters);
}

sub Convert_BIE_Parameters_To_Galacticus {
    my $config = shift;
    my @values = @_;

    # Extract parameters from config file.
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
    die("Convert_BIE_Parameters_To_Galacticus: number of supplied values does not match number of parameters")
	unless ( scalar(@values) == $parameterCount );
    # Map values to parameters, undoing any logarithmic mapping.
    my $j = -1;
    my %parameterValues;
    for(my $i=0;$i<scalar(@parameters);++$i) {
	if ( exists($parameters[$i]->{'prior'}) ) {
	    ++$j;
	    $parameterValues{$parameters[$i]->{'name'}} = $values[$j];
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
    return $newParameters;
}

sub Sample_Models {
    # Generate a sample of models from the posterior distribution.
    my $config    = shift;
    my %arguments = %{$_[0]};
    # Find the work directory.
    my $workDirectory = $config->{'workDirectory'};
    # Get a hash of the parameter values.
    my $compilationFile;
    if ( exists($arguments{'compilationOverride'}) ) {
	$compilationFile = $arguments{'compilationOverride'};
    } else {
	$compilationFile = $config->{'compilation'};
    }
    (my $constraintsRef, my $parameters) = &Parameters::Compilation($compilationFile,$config->{'baseParameters'});
    my @constraints = @{$constraintsRef};
    # Parse the statefile to find all parameter values sampled by the chains.
    my @chainParameters;
    open(iHndl,$workDirectory."/mcmc/galacticusBIE.statelog");
    while ( my $line = <iHndl> ) {
	unless ( $line =~ m/^\"/ ) {
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    my @columns = split(/\s+/,$line);
	    push(@{$chainParameters[++$#chainParameters]},@columns[5..$#columns]);
	}
    }
    close(iHndl);   
    # Select viable parameter sets.
    my @outlierChains = split(/,/,$arguments{'outliers'});
    my @chainParametersViable;
    my $chainCount = $config->{'threads'};
    for(my $i=0;$i<scalar(@chainParameters);++$i) {
	my $accept = 1;
	$accept = 0
	    if ( $arguments{'sampleFrom'} > 0 && $i < scalar(@chainParameters)-$arguments{'sampleFrom'} );
	my $chainNumber = $i % $chainCount;
	foreach ( @outlierChains ) {
	    $accept = 0
		if ( $chainNumber == $_ );
	}
	push(@{$chainParametersViable[++$#chainParametersViable]},@{$chainParameters[$i]})
	    if ( $accept == 1 );
    }
    # Sample parameters.
    $arguments{'sampleCount'} = scalar(@chainParametersViable)
	if ( $arguments{'sampleCount'} < 0 );
    my $sampleIndex = pdl long(scalar(@chainParametersViable)*random($arguments{'sampleCount'}));
    # Run model for each sample.
    my $sampleDirectory = $workDirectory."/posteriorSampleModels/";
    $sampleDirectory = $arguments{'sampleDirectory'}."/"
       if ( exists($arguments{'sampleDirectory'}) );
    my @pbsStack;
    for (my $i=0;$i<nelem($sampleIndex);++$i) {
	# Create an output directory.
	my $modelDirectory = $sampleDirectory.$i."/";
	system("mkdir -p ".$modelDirectory);
	my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
	# Check if the model has already been run.
	unless ( -e $galacticusFileName ) {
	    # Convert these values into a parameter array.
	    my $j = $sampleIndex->(($i))->sclr();
	    my $currentConfig = clone($config);
	    my $newParameters = &Convert_BIE_Parameters_To_Galacticus($currentConfig,@{$chainParametersViable[$j]});    
	    # Increment the random number seed.
	    $parameters->{'parameter'}->{'randomSeed'}->{'value'} += $config->{'threadsPerNode'};
	    # Clone parameters.
	    my $currentParameters = clone($parameters);
	    # Apply to parameters.
	    $currentParameters->{'parameter'}->{$_->{'name'}}->{'value'} = $_->{'value'}
	       foreach ( @{$newParameters->{'parameter'}} );    
	    # Apply any parameter overrides from the command line.
	    foreach ( keys(%arguments) ) {
		if ( $_ =~ m/^parameterOverride:(.+)/ ) {
		    my $parameterName = $1;
		    $currentParameters->{'parameter'}->{$parameterName}->{'value'} = $arguments{$_};
		}
	    }
	    # Specify the output file name.
	    $currentParameters->{'parameter'}->{'galacticusOutputFileName'}->{'value'} = $galacticusFileName;
	    # Write the modified parameters to file.
	    &Output($currentParameters,$modelDirectory."parameters.xml");
	    # Create a batch script for PBS.
	    my $batchScriptFileName = $modelDirectory."/launch.pbs";
	    open(oHndl,">".$batchScriptFileName);
	    print oHndl "#!/bin/bash\n";
	    print oHndl "#PBS -N ".$config->{'name'}."_ppc".$i."\n";
	    print oHndl "#PBS -l walltime=".$config->{'walltimeLimit'}."\n"
		if ( exists($config->{'walltimeLimit'}) );
	    print oHndl "#PBS -l mem=".$config->{'memoryLimit'}."\n"
		if ( exists($config->{'memoryLimit'}) );
	    my $threadsPerNode = 1;
	    $threadsPerNode = $config->{'threadsPerNode'}
	    if ( exists($config->{'threadsPerNode'}) );
	    print oHndl "#PBS -l nodes=1:ppn=".$threadsPerNode."\n";
	    print oHndl "#PBS -j oe\n";
	    print oHndl "#PBS -o ".$modelDirectory."/launch.log\n";
	    print oHndl "#PBS -V\n";
	    print oHndl "cd \$PBS_O_WORKDIR\n";
	    if ( exists($config->{'environment'}) ) {
		my @environment;
		if ( UNIVERSAL::isa($config->{'environment'},"ARRAY") ) {
		    push(@environment,@{$config->{'environment'}});
		} else {
		    push(@environment,  $config->{'environment'} );
		}
		foreach ( @environment ) {
		    print oHndl "export ".$_."\n";
		}
	    }
	    print oHndl "ulimit -t unlimited\n";
	    print oHndl "ulimit -c unlimited\n";
	    print oHndl "export OMP_NUM_THREADS=".$config->{'threadsPerNode'}."\n";
	    print oHndl "mpirun --bynode -np 1 Galacticus.exe ".$modelDirectory."parameters.xml\n";
	    foreach my $constraint ( @constraints ) {
		# Parse the definition file.
		my $xml = new XML::Simple;
		my $constraintDefinition = $xml->XMLin($constraint->{'definition'});	    
		# Insert code to run the analysis code.
		my $analysisCode = $constraintDefinition->{'analysis'};
		print oHndl $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".xml --modelDiscrepancies ".$workDirectory."/modelDiscrepancy\n";
	    }
	    close(oHndl);
	    # Queue the calculation.
	    push(
		@pbsStack,
		$batchScriptFileName
		);
	}
    }
    # Send jobs to PBS.
    my $jobMaximum = 10;
    $jobMaximum = $arguments{'pbsJobMaximum'}
        if ( exists($arguments{'pbsJobMaximum'}) );
    &PBS_Submit($jobMaximum,@pbsStack)
	if ( scalar(@pbsStack) > 0 );
    # Return the number of models sampled.
    return nelem($sampleIndex);
}

sub PBS_Submit {
    # Submit jobs to PBS and wait for them to finish.
    my $jobMaximum = shift;
    my @pbsStack   = @_;
    my %pbsJobs;
    # Submit jobs and wait.
    print "Waiting for PBS jobs to finish...\n";
    while ( scalar(keys %pbsJobs) > 0 || scalar(@pbsStack) > 0 ) {
	# Find all PBS jobs that are running.
	my %runningPBSJobs;
	undef(%runningPBSJobs);
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^Job\sId:\s+(\S+)/ ) {$runningPBSJobs{$1} = 1};
	}
	close(pHndl);
	foreach my $jobID ( keys(%pbsJobs) ) {
	    unless ( exists($runningPBSJobs{$jobID}) ) {
		print "PBS job ".$jobID." has finished.\n";
		# Remove the job ID from the list of active PBS jobs.
		delete($pbsJobs{$jobID});
	    }
	}
	# If fewer than the maximum number of jobs are in the queue, pop one off the stack.
	if ( scalar(@pbsStack) > 0 && scalar(keys %pbsJobs) < $jobMaximum ) {
	    my $batchScript = pop(@pbsStack);
	    # Submit the PBS job.
	    open(pHndl,"qsub ".$batchScript."|");
	    my $jobID = "";
	    while ( my $line = <pHndl> ) {
	    	if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
	    }
	    close(pHndl);	    
	    # Add the job number to the active job hash.
	    unless ( $jobID eq "" ) {
	    	$pbsJobs{$jobID} = 1;
	    }
	    sleep 1;
	} else {
	    # Wait.
	    sleep 5;
	}
    }
}

1;
