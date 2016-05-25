# Contains a Perl module which implements various useful functionality for constructing parameter files for
# Galacticus when fitting to constraints.

package Parameters;
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
use XML::LibXML;
use XML::Simple;
use XML::Twig;
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
use Clone qw(clone);
use List::Util;
use Storable;
require List::ExtraUtils;
require Galacticus::Launch::PBS;

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
	    push(@parameters,&ExtraUtils::as_array($_->{'parameter'}))
		foreach (  &ExtraUtils::as_array($config->{'parameters'}) );
	    delete($config->{'parameters'});
	    @{$config->{'parameters'}->{'parameter'}} = @parameters;
	}
    }
    return $config;
}

sub Output {
    # Output a set of parameters to file.
    my $parameters = shift;
    my $fileName   = shift;
    # Create an XML output object.
    my $xmlOutput  = new XML::Simple (RootName=>"parameters");
    # Output the parameters to file.
    open(pHndl,">".$fileName);
    print pHndl $xmlOutput->XMLout($parameters);
    close pHndl;
}

sub Compilation {
    # Process a compilation of constraints, adjusting parameters as necessary and return a parameter hash.
    my $compilationFileName    = shift;
    my $baseParametersFileName = shift;
    # Create an XML worker object.
    my $xml = new XML::Simple;
    # Retrieve the compilation file.
    my $compilationFilePath = "constraints/compilations/".$compilationFileName;
    my $compilation;
    if ( -e $compilationFilePath.".store" ) {
	$compilation = retrieve($compilationFilePath.".store");
    } else {
	# Create an XML worker object.
	# Parse the constraint compilation file.
	$compilation = $xml->XMLin($compilationFilePath);
    }
    # Specify default values for parameters which can be adjusted by constraint definitions.
    my %outputRedshifts;
    my %outputLuminosities;
    my %optionsOn;
    my $haloMassResolution;
    my $haloMassMinimum;
    my $haloMassMaximum;
    # Parse the base set of parameters.
    my $parameters;
    if ( -e $baseParametersFileName.".store" ) {
	$parameters = retrieve($baseParametersFileName.".store");
    } else {
	$parameters = $xml->XMLin($baseParametersFileName);
    }
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
	my $constraintDefinition;
	if ( -e $constraint->{'definition'}.".store" ) {
	    $constraintDefinition = retrieve($constraint->{'definition'}.".store");
	} else {
	    $constraintDefinition = $xml->XMLin($constraint->{'definition'},KeyAttr => "");
	}
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
    my @outputRedshiftList = keys(%outputRedshifts);
    if ( scalar(@outputRedshiftList) == 0 ) {
	push(@outputRedshiftList,"0.0000");
    } else {
	@outputRedshiftList = sort(@outputRedshiftList);
    }
    # Modify the set of output redshifts.
    $parameters->{'outputRedshifts'}->{'value'} = join(" ",@outputRedshiftList);
    # Modify the minimum and maximum halo masses.
    $haloMassResolution = 5.00e09 unless ( defined($haloMassResolution) );
    $haloMassMinimum    = 1.00e10 unless ( defined($haloMassMinimum   ) );
    $haloMassMaximum    = 1.01e10 unless ( defined($haloMassMaximum   ) );
    $parameters->{'mergerTreeBuildMassResolutionFixed'}->{'value'} = $haloMassResolution;
    $parameters->{'mergerTreeBuildHaloMassMinimum'    }->{'value'} = $haloMassMinimum   ;
    $parameters->{'mergerTreeBuildHaloMassMaximum'    }->{'value'} = $haloMassMaximum   ;
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

sub Convert_Parameters_To_Galacticus {
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
    die("Convert_Parameters_To_Galacticus: number of supplied values does not match number of parameters")
	unless ( scalar(@values) == $parameterCount );

    # Map values to parameters, undoing any logarithmic mapping.
    my $j = -1;
    my %parameterValues;
    for(my $i=0;$i<scalar(@parameters);++$i) {
	if ( exists($parameters[$i]->{'prior'}) ) {
	    ++$j;
	    $parameterValues{$parameters[$i]->{'name'}} = $values[$j];	  
	}
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
	    if ( exists($parameters[$i]->{'define'}) ) {
		die ("Convert_Parameters_To_Galacticus: cannot specify a prior for a defined parameter")
		    if ( exists($parameters[$i]->{'prior'}) );
		# Attempt to replace named parameters in the definition with their values.
		while ( $parameters[$i]->{'define'} =~ m/\%\[([a-zA-Z0-9_\.\-\>]+)\]/ ) {
		    my $parameterName = $1;
		    if ( exists($parameterValues{$parameterName}) ) {
			$parameters[$i]->{'define'} =~ s/\%\[$parameterName\]/$parameterValues{$parameterName}/g;
		    } else {
			++$failCount;
			last;
		    }
		}
		$parameterValues{$parameters[$i]->{'name'}} = eval($parameters[$i]->{'define'})
		    unless ( $parameters[$i]->{'define'} =~ m/\%\[([a-zA-Z0-9_\.\-\>]+)\]/ );
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

sub Sample_Models {
    # Generate a sample of models from the posterior distribution.
    my $config    =   shift() ;
    my %arguments = %{shift()};
    # Find the work directory.
    my $workDirectory = $config->{'likelihood'}->{'workDirectory'};
    # Get a hash of the parameter values.
    my $compilationFile;
    if ( exists($arguments{'compilationOverride'}) ) {
	$compilationFile = $arguments{'compilationOverride'};
    } else {
	$compilationFile = $config->{'likelihood'}->{'compilation'};
    }
    (my $constraintsRef, my $parameters) = &Parameters::Compilation($compilationFile,$config->{'likelihood'}->{'baseParameters'});
    my @constraints = @{$constraintsRef};
    # Determine number of chains.
    my $logFileRoot = $config->{'simulation'}->{'logFileRoot'};
    my $chainCount  = 0;
    while () {
	++$chainCount;
	my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,$chainCount);
	last
	unless ( -e $chainFileName );
    }
    print "Found ".$chainCount." chains\n";
    # Build a list of outlier chains.
    my @outlierChains = split(/,/,$arguments{'outliers'});
    print scalar(@outlierChains) > 0 ? "Outlier chains are: ".join(", ",@outlierChains)."\n" : "No outlier chains\n";
    # Parse the chains to find all parameter values sampled by the chains.
    my @chainParameters;
    for(my $i=0;$i<$chainCount;++$i) {
	# Skip outlier chains.
	next
	    if ( grep {$_ eq $i} @outlierChains );
	# Parse the chain file.
	my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,$i);
	my $step = 0;
	open(iHndl,$chainFileName);
	while ( my $line = <iHndl> ) {
	    unless ( $line =~ m/^\"/ ) {
		++$step;
		my @columns = split(" ",$line);
		my $accept = 1;
		# Skip unconverged states unless explicitly allowed.
		$accept = 0
		    if ( $columns[3] eq "F" && ( ! exists($arguments{'useUnconverged'}) || $arguments{'useUnconverged'} eq "no" ) );
		# Skip chains before the given start point.
		$accept = 0
		    if ( exists($arguments{'sampleFrom'}) && $step < $arguments{'sampleFrom'} );
		push(@{$chainParameters[++$#chainParameters]},@columns[6..$#columns])
		    if ( $accept );
	    }
	}
    }
    close(iHndl);
    print "Found ".scalar(@chainParameters)." chain states\n";
    # Sample parameters.
    $arguments{'sampleCount'} = scalar(@chainParameters)
	if ( $arguments{'sampleCount'} < 0 );
    my $sampleIndex = pdl long(scalar(@chainParameters)*random($arguments{'sampleCount'}));
    # Run model for each sample.
    my $sampleDirectory = exists($arguments{'sampleDirectory'}) ? $arguments{'sampleDirectory'}."/" : $workDirectory."/posteriorSample/";
    my @pbsStack;
    for (my $i=0;$i<nelem($sampleIndex);++$i) {
	# Create an output directory.
	my $modelDirectory = $sampleDirectory.$i."/";
	system("mkdir -p ".$modelDirectory);
	my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
	# Check if the model has already been run.
	unless ( -e $galacticusFileName ) {
	    # Convert these values into a parameter array.
	    my $j             = $sampleIndex->(($i))->sclr();
	    my $currentConfig = clone($config);
	    my $newParameters = &Convert_Parameters_To_Galacticus($currentConfig,@{$chainParameters[$j]});    
	    # Increment the random number seed.
	    $parameters->{'randomSeed'}->{'value'} += $config->{'likelihood'}->{'threads'};
	    # Clone parameters.
	    my $currentParameters = clone($parameters);
	    # Apply to parameters.
	    for my $newParameterName ( keys(%{$newParameters}) ) {
		my $parameter = $currentParameters;
		foreach ( split(/\-\>/,$newParameterName) ) {
		    $parameter->{$_}->{'value'} = undef()
			unless ( exists($parameter->{$_}) );
		    $parameter = $parameter->{$_};
		}
		$parameter->{'value'} = $newParameters->{$newParameterName};
	    }
	    # Apply any parameters from command line.
	    foreach my $argument ( keys(%arguments) ) {
		if ( $argument =~ m/^parameterOverride:(.*)/ ) {
		    my @parametersSet = split(":",$1);
		    my $parameter     = $currentParameters;
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
	    # Specify the output file name.
	    $currentParameters->{'galacticusOutputFileName'}->{'value'} = $galacticusFileName;
	    # Write the modified parameters to file.
	    &Output($currentParameters,$modelDirectory."parameters.xml");
	    # Construct the tasks to perform.
	    my $command;
	    if ( exists($config->{'likelihood'}->{'environment'}) ) {
		my @environment;
		if ( UNIVERSAL::isa($config->{'likelihood'}->{'environment'},"ARRAY") ) {
		    push(@environment,@{$config->{'likelihood'}->{'environment'}});
		} else {
		    push(@environment,  $config->{'likelihood'}->{'environment'} );
		}
		foreach ( @environment ) {
		    $command .= "export ".$_."\n";
		}
	    }
	    $command .= "ulimit -t unlimited\n";
	    $command .= "ulimit -c unlimited\n";
	    $command .= "export OMP_NUM_THREADS=".$config->{'likelihood'}->{'threads'}."\n";
	    $command .= "mpirun --bynode -np 1 Galacticus.exe ".$modelDirectory."parameters.xml\n";
	    foreach my $constraint ( @constraints ) {
		# Parse the definition file.
		my $xml = new XML::Simple;
		my $constraintDefinition = $xml->XMLin($constraint->{'definition'});	    
		# Insert code to run the analysis code.
		my $analysisCode = $constraintDefinition->{'analysis'};
		(my $plotLabel = $constraintDefinition->{'label'}) =~ s/\./_/g;
		$command .= $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5 --plotFile ".$modelDirectory."/".$plotLabel.".pdf --outputFile ".$modelDirectory."/".$constraintDefinition->{'label'}."Likelihood.xml --modelDiscrepancies ".$workDirectory."/modelDiscrepancy\n";
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
    &PBS::SubmitJobs(\%arguments,@pbsStack)
     	if ( scalar(@pbsStack) > 0 );
    # Return the number of models sampled.
    return (nelem($sampleIndex), $sampleDirectory);
}

sub Maximum_Likelihood_Parameters {
    my $config    =   shift() ;
    my %arguments = %{shift()};
    # Get a hash of the parameter values.
    (my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'likelihood'}->{'compilation'},$config->{'likelihood'}->{'baseParameters'});
   # Determine the MCMC directory.
    my $logFileRoot = $config->{'simulation'}->{'logFileRoot'};
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
    # Convert these values into a parameter array.
    my $newParameters = &Parameters::Convert_Parameters_To_Galacticus($config,@maximumLikelihoodParameters);
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
    return $parameters;
}

sub Apply_Command_Line_Parameters {
    my $parameters =   shift() ;
    my %arguments  = %{shift()};
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
}

1;
