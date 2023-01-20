# Contains a Perl module which implements various useful functionality for constructing parameter files for
# Galacticus when fitting to constraints.

package Galacticus::Constraints::Parameters;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::LibXML;
use XML::Simple;
use XML::Twig;
use Data::Dumper;
use Clone qw(clone);
use List::Util;
use List::ExtraUtils;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Local;
use Galacticus::Launch::Slurm;
use Scalar::Util qw(reftype looks_like_number);
use Storable qw(dclone);
use PDL;
use PDL::NiceSlice;

sub parseConfig {
    # Parse a config file, handling xi:xinclude elements.
    my $configFileName = shift();
    my $xml            = new XML::Simple();
    my $parser         = XML::LibXML->new();
    my $dom            = $parser->load_xml(location => $configFileName);
    $parser->process_xincludes($dom);
    my $config = $xml->XMLin($dom->serialize());
    return $config;
}

sub logFileRoot {
    # Return the root of the log (chain) file names.
    my $config  =   shift() ;
    my %options = %{shift()};
    return exists($options{'chainRoot'}) ? $options{'chainRoot'} : $config->{'posteriorSampleSimulation'}->{'logFileRoot'}->{'value'};
}

sub stepCount {
    # Return a count of the number of steps taken.
    my $config  =   shift() ;
    my %options = %{shift()};
    # Read the first chain.
    my $logFileRoot   = &logFileRoot($config,\%options);
    my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,0);
    my $stepCount;
    open(my $chainFile,$chainFileName);
    while ( my $line = <$chainFile> ) {
	unless ( $line =~ m/^\"/ ) {
	    my @columns = split(" ",$line);
	    $stepCount = $columns[0];
	}
    }
    die("Galacticus::Constraints::Parameters::parametersCount(): could not determine number of steps")
	unless ( defined($stepCount) );
    return $stepCount;
}

sub chainCount {
    # Return a count of the number of chains used.
    my $config  =   shift() ;
    my %options = %{shift()};
    # Determine number of chains.
    my $logFileRoot = &logFileRoot($config,\%options);
    my $chainCount  = 0;
    while () {
	++$chainCount;
	my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,$chainCount);
	last
	    unless ( -e $chainFileName );
    }
    return $chainCount;
}

sub parameterCount {
    # Return a count of the number of parameters used.
    my $config  =   shift() ;
    my %options = %{shift()};
    # Read the first chain.
    my $logFileRoot   = &logFileRoot($config,\%options);
    my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,0);
    my $parameterCount;
    open(my $chainFile,$chainFileName);
    while ( my $line = <$chainFile> ) {
	unless ( $line =~ m/^\"/ ) {
	    my @columns = split(" ",$line);
	    $parameterCount = scalar(@columns)-6;
	}
    }
    die("Galacticus::Constraints::Parameters::parametersCount(): could not determine number of parameters")
	unless ( defined($parameterCount) );
    $parameterCount /= 2
	if ( $config->{'posteriorSampleSimulation'}->{'value'} eq "particleSwarm" );
    return $parameterCount;
}

sub maximumPosteriorParameterVector {
    my $config  =   shift() ;
    my %options = %{shift()};
    # Determine the MCMC directory.
    my $logFileRoot = &logFileRoot($config,\%options);
    (my $mcmcDirectory  = $logFileRoot) =~ s/\/[^\/]+$//;    
    # Determine number of chains.
    my $chainCount = &chainCount($config,\%options);
    # Parse the chains to find the maximum likelihood model.
    my $maximumLikelihood;
    my @maximumLikelihoodParameters;
    my @chainFiles;
    for(my $i=0;$i<$chainCount;++$i) {
	next
	    unless
	    (
	     ! exists($options{'chain'})
	     ||
	     $options{'chain'} eq "all"
	     ||
	     $options{'chain'} == $i
	    );
	my $chainFileName         = sprintf("%s_%4.4i.log"        ,$logFileRoot,$i);
	my $chainFilePreviousName = sprintf("%sPrevious_%4.4i.log",$logFileRoot,$i);
	push(@chainFiles,$chainFileName        );
	push(@chainFiles,$chainFilePreviousName)
	    if ( exists($options{'includePrevious'}) && $options{'includePrevious'} eq "yes" && -e $chainFilePreviousName );
    }
    foreach my $chainFile ( @chainFiles ) {
	my $step = 0;
	open(iHndl,$chainFile);
	while ( my $line = <iHndl> ) {
	    ++$step;
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    my @columns = split(/\s+/,$line);
	    # Determine if state is accepted.
	    my $accept = 1;
	    $accept = 0
		if ( exists($options{'burnCount'}) && $step <= $options{'burnCount'} );
	    if ( $accept == 1 && (! defined($maximumLikelihood) || $columns[4] > $maximumLikelihood ) ) {
		$maximumLikelihood           = $columns[4];
		@maximumLikelihoodParameters = @columns[6..$#columns];
	    }
	}
	close(iHndl);
    }
    # Check that parameters were found.
    die('found no acceptable states')
	unless ( @maximumLikelihoodParameters );
    # Convert parameters to a PDL.
    my $maximumLikelihoodVector;
    if ( $config->{'posteriorSampleSimulation'}->{'value'} eq "particleSwarm" ) {
	my $indexMaximum = scalar(@maximumLikelihoodParameters)/2-1;
	$maximumLikelihoodVector = pdl @maximumLikelihoodParameters[0..$indexMaximum];
    } else {
	$maximumLikelihoodVector = pdl @maximumLikelihoodParameters                  ;
    }
    # Return the vector.
    return ($maximumLikelihoodVector, $maximumLikelihood);
}

sub parameterVector {
    my $config  =   shift() ;
    my $chain   =   shift() ;
    my $state   =   shift() ;
    my %options = %{shift()};
    # Determine the MCMC directory.
    my $logFileRoot    = &logFileRoot($config,\%options);
    (my $mcmcDirectory = $logFileRoot) =~ s/\/[^\/]+$//;    
    # Determine number of chains.
    my $chainCount = &chainCount($config,\%options);
    die('Galacticus::Constraints::Parameters::parameterVector(): chain index out of range')
	if ( $chain >= $chainCount );
    # Parse the state option.
    unless ( looks_like_number($state) ) {
	if ( $state eq "initial" ) {
	    $state =  1;
	} elsif ( $state eq "final" ) {
	    $state = -1;
	} else {
	    die('Galacticus::Constraints::Parameters::parameterVector(): invalid state');
	}
    }
    # Parse the chain.
    my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,$chain);
    my $stateCount    = 0;
    my @state;
    open(my $chainFile,$chainFileName);
    while ( my $line = <$chainFile> ) {
	++$stateCount;
	$line =~ s/^\s*//;
	$line =~ s/\s*$//;
	my @columns = split(/\s+/,$line);
	# Determine whether to use this state.
	if (
	    $state == -1
	    ||
	    $state == $stateCount
	    ) {
	    @state = @columns[6..$#columns];
	}
    }
    close($chainFile);
    # Convert parameters to a PDL.
    my $parameterVector = pdl @state;
    # Return the vector.
    return $parameterVector;
}

sub parameterMatrix {
    my $config  =   shift() ;
    my $chain   =   shift() ;
    my %options = %{shift()};
    # Determine the MCMC directory.
    my $logFileRoot    = &logFileRoot($config,\%options);
    (my $mcmcDirectory = $logFileRoot) =~ s/\/[^\/]+$//;    
    # Determine number of chains.
    my $chainCount = &chainCount($config,\%options);
    die('Galacticus::Constraints::Parameters::parameterMatrix(): chain index out of range')
	if ( $chain ne "all" && $chain >= $chainCount );
    # Parse the chain.
    my $chainMinimum;
    my $chainMaximum;
    if ( $chain eq "all" ) {
	$chainMinimum = 0;
	$chainMaximum = $chainCount-1;
    } else {
	$chainMinimum = $chain;
	$chainMaximum = $chain;
    }
    my $countChainsActive = $chainMaximum-$chainMinimum+1;
    my $parameterMatrix = pdl zeros(&parameterCount($config,\%options),&stepCount($config,\%options)*$countChainsActive);
    my $logLikelihood   = pdl zeros(                                   &stepCount($config,\%options)*$countChainsActive);
    my $stateCount      = -1;
    for(my $chainIndex=$chainMinimum;$chainIndex <= $chainMaximum;++$chainIndex) {
	my $chainFileName = sprintf("%s_%4.4i.log",$logFileRoot,$chainIndex);
	open(my $chainFile,$chainFileName);
	while ( my $line = <$chainFile> ) {
	    ++$stateCount;
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    my @columns = split(/\s+/,$line);
	    $logLikelihood->(($stateCount)) .= $columns[4];
	    my @state = @columns[6..$#columns];
	    for(my $i=0;$i<scalar(@state);++$i) {
		$parameterMatrix->(($i),($stateCount)) .= $state[$i];
	    }
	}
	close($chainFile);
    }
    if ( exists($options{'includeLikelihood'}) && $options{'includeLikelihood'} ) {
	return ($parameterMatrix, $logLikelihood);
    } else {
	return  $parameterMatrix                 ;
    }
}

sub parameterFind {
    # Locate a named parameter in a parameters structure.
    my $parameters    = shift();
    my $parameterName = shift();
    my $valueIndex             ;
    if ( $parameterName =~ m/^(.*)\{(\d+)\}$/ ) {
	$parameterName = $1;
	$valueIndex    = $2;
    }
    my $parameter = $parameters;
    foreach ( split(/::/,$parameterName) ) {
	# Check if the parameter name contains an array reference.
	if ( $_ =~ m/^(.*)\[(\d+)\]$/ ) {
	    # Parameter name contains array reference. Step through to the relevant parameter in the list. If the parameter is
	    # not an array, allow this only if the array index given is zero.
	    return ( undef(), undef() )
		unless ( defined($parameter->{$1}) );
	    if ( reftype($parameter->{$1}) eq "ARRAY" ) {
		$parameter->{$1}->[$2]->{'value'} = undef()
		    unless ( scalar(@{$parameter->{$1}}) > $2 );
		$parameter = $parameter->{$1}->[$2];
	    } else {
		die('Galacticus::Constraints::Parameters::parameterFind(): attempt to access non-existant array')
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
    return ( $parameter, $valueIndex );
}

sub parameterValueGet {
    # Get a value from a parameter.
    my $parameter  = shift();
    my $valueIndex = shift();
    my $value               ;
    if ( defined($valueIndex) ) {
	my @values = split(" ",$parameter->{'value'});
	$value  = $values[$valueIndex];
    } else {
	$value  = $parameter->{'value'};
    }
    return $value;
}

sub parameterValueSet {
    # Set a value in a parameter;
    my $parameter  = shift();
    my $valueIndex = shift();
    my $value      = shift();
    if ( defined($valueIndex) ) {
	my @values = split(" ",$parameter->{'value'});
	$values[$valueIndex] = $value ;
	$parameter->{'value'} = join(" ",@values);
    } else {
	$parameter->{'value'} = $value;
    }
}

sub applyCommandLineParameters {
    # Apply any parameters from command line.
    my $parameters =   shift() ;
    my %options    = %{shift()};
    foreach my $argument ( keys(%options) ) {
	if ( $argument =~ m/^parameter:(.*)/ ) {
	    my $parameterName = $1;
	    (my $parameter, my $valueIndex) = &parameterFind($parameters,$parameterName);
	    &parameterValueSet($parameter,$valueIndex,$options{$argument});
	}
    }
}

sub parameterVectorApply {
    # Convert a PDL vector of parameter values into a parameter structure.
    my $config          =   shift() ;
    my $model           =   shift() ;
    my $parameterVector =   shift() ;
    my %options         = %{shift()};
    # Get a hash of the parameter values.
    my $xml            = new XML::Simple();
    my $parser         = XML::LibXML->new();
    my $fileName       = exists($options{'baseParameters'}) ? $options{'baseParameters'} : $model->{'posteriorSampleLikelihood'}->{'baseParametersFileName'}->{'value'};
    my $dom            = $parser->load_xml(location => $fileName);
    $parser->process_xincludes($dom);
    my $parameters     = $xml->XMLin($dom->serialize());
    # Apply vector of parameter values to parameters structure.
    my $i               = -1;  
    foreach my $modelParameter ( &List::ExtraUtils::as_array($config->{'posteriorSampleSimulation'}->{'modelParameter'}) ) {
	# Determine the value to set.
	my $parameterValue;
	if      ( $modelParameter->{'value'} eq "active" ) {
	    ++$i;
	    $parameterValue = $parameterVector->(($i))->sclr();
	} elsif ( $modelParameter->{'value'} eq "derived") {
	    $parameterValue = $modelParameter->{'definition'}->{'value'};
	} else {
	    die('Galacticus::Constraints::Parameters::parameterVectorApply(): unknown parameter type');
	}
	# Find the parameter and set its value.
	(my $parameter, my $valueIndex) = &parameterFind($parameters,$modelParameter->{'name'}->{'value'});
	&parameterValueSet($parameter,$valueIndex,$parameterValue)
	    if ( grep {$modelParameter->{'name'}->{'value'} eq $_} @{$model->{'parameters'}} );
    }
    # Find any derived parameters specified as options.
    my @derivedParameters;
    if ( exists($options{'derivedParameter'}) ) {	
	foreach my $derivedParameter ( &List::ExtraUtils::as_array($options{'derivedParameter'}) ) {
	    if ( $derivedParameter =~ m/(.*)=(.*)/ ) {
		my $parameter = {
		    value      => "derived",
		    name       => {value => $1},
		    definition => {value => $2}
		};
		push(@derivedParameters,$parameter);
		(my $parameterDependent, my $valueIndexDependent) = &parameterFind($parameters,$parameter->{'name'}->{'value'});
		my $valueDependent = &parameterValueGet($parameterDependent,$valueIndexDependent);
		&parameterValueSet($parameterDependent,$valueIndexDependent,$parameter->{'definition'}->{'value'});
	    } else {
		die("Galacticus::Constraints::Parameters::parameterVectorApply(): Can not parse derived parameter definition '".$derivedParameter."'");
	    }
	}
    }
    # Resolve any derived parameters.
    my $dependenciesResolved = 0;
    while ( ! $dependenciesResolved ) {
	my $progress = 0;
	$dependenciesResolved = 1;
	foreach my $modelParameter ( &List::ExtraUtils::as_array($config->{'posteriorSampleSimulation'}->{'modelParameter'}), @derivedParameters ) {
	    # Skip parameters not applicable to this model.
	    next
		unless (
		    ( grep {$modelParameter->{'name'}->{'value'} eq $_                     }  @{$model->{'parameters'}} )
		    ||  
		    ( grep {$modelParameter->{'name'}->{'value'} eq $_->{'name'}->{'value'}}  @derivedParameters        )
		);
	    # Skip non-defined parameters.
	    next
		unless ( $modelParameter->{'value'} eq "derived" );
	    # Find the parameter.
	    (my $parameterDerived, my $valueIndexDerived) = &parameterFind($parameters,$modelParameter->{'name'}->{'value'});
	    my $definition = &parameterValueGet($parameterDerived,$valueIndexDerived);
	    # Look for names of dependencies.
	    if ( $definition =~ m/\%\[([^\%]+)\]/ ) {
		# Get the dependent parameter name.
		my $parameterDependentName = $1;
		# First replace any constants (which are known to the libmatheval library, see: https://www.gnu.org/software/libmatheval/manual/libmatheval.html#evaluator_005fcreate).
		my $ln10    = log(10.0);
		$definition =~ s/(^|\W)ln10(\W|$)/$1$ln10$2/g;
                # Extract the value of the dependent parameter.  
		(my $parameterDependent, my $valueIndexDependent) = &parameterFind($parameters,$parameterDependentName);
		my $valueDependent = &parameterValueGet($parameterDependent,$valueIndexDependent);
		# Check if the dependent parameter is resolved.
		if ( $valueDependent !~ m/\%\[/ ) {
		    # Replace the dependent parameter name with the value in the target.
		    $definition =~ s/\%\[$parameterDependentName\]/$valueDependent/g;
		    $progress   = 1;
		}
		# Check if target can now be evaluated.
		if ( $definition  !~ m/\%\[/ ) {
		    $definition = eval($definition);
		} else {
		    $dependenciesResolved = 0;
		}
		# Store the definition.
		&parameterValueSet($parameterDerived,$valueIndexDerived,$definition);
	    }
	}
	die('Galacticus::Constraints::Parameters::parameterVectorApply(): unable to resolve parameter dependencies')
	    unless ( $progress || $dependenciesResolved );
    }
    # Apply any command line parameters.
    &applyCommandLineParameters($parameters,\%options);
    return $parameters;
}

sub parametersOutput {
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
    system("mkdir -p `dirname ".$fileName."`");
    open(my $outputFile,">".$fileName);
    print $outputFile $xmlOutput->XMLout($parametersCopy);
    close $outputFile;
}

sub modelList {
    # Construct a list of all "galaxyPopulation" models.
    my $config  =   shift() ;
    my %options = %{shift()};
    # Begin building a stack of likelihood models to evaluate.
    my @galaxyPopulationModels;
    my $instance               = 0;
    my @allParameterNames      = map {$_->{'name'}->{'value'}} &List::ExtraUtils::as_array($config->{'posteriorSampleSimulation'}->{'modelParameter'});
    my @modelStack             = ( 
	{
	    posteriorSampleLikelihood => $config->{'posteriorSampleLikelihood'},
	    parameters                => \@allParameterNames
	} 
	);
    # Iterate over likelihood models.
    while ( scalar(@modelStack) > 0 ) {
	# Pop a model off of the stack.
	my $model = shift(@modelStack);
	if ( $model->{'posteriorSampleLikelihood'}->{'value'} eq "galaxyPopulation" ) {
	    ++$instance;
	    push(@galaxyPopulationModels,$model)
		if ( ! exists($options{'modelInstance'}) || $options{'modelInstance'} eq "all" || $options{'modelInstance'} == $instance );
	}
	# Add any sub-models to the stack.
	if ( exists($model->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}) ) {
	    my @likelihoods           = &List::ExtraUtils::as_array($model->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'});
	    my @parameterMapsActive   = &List::ExtraUtils::as_array($model->{'posteriorSampleLikelihood'}->{'parameterMap'             })
		if ( exists($model->{'posteriorSampleLikelihood'}->{'parameterMap'        }) );
	    my @parameterMapsInactive = &List::ExtraUtils::as_array($model->{'posteriorSampleLikelihood'}->{'parameterInactiveMap'     })
		if ( exists($model->{'posteriorSampleLikelihood'}->{'parameterInactiveMap'}) );
	    for(my $i=0;$i<scalar(@likelihoods);++$i) {
		my @parameterNames;
		push(@parameterNames,split(" ",$parameterMapsActive  [$i]->{'value'}))
		    if ( @parameterMapsActive   );
		push(@parameterNames,split(" ",$parameterMapsInactive[$i]->{'value'}))
		    if ( @parameterMapsInactive );
		my $childModel;
		$childModel->{'posteriorSampleLikelihood'} = $likelihoods[$i];
		$childModel->{'parameters'               } = \@parameterNames;
		push(
		    @modelStack,
		    $childModel
		    );
	    }
	}
    }
    return @galaxyPopulationModels;
}

sub stageModel {
    # Stage a model.
    my $parameterFileName   =   shift() ;
    my $label               =   shift() ;
    my %options             = %{shift()};
    (my $modelDirectoryName = $parameterFileName) =~ s/[^\/]+$//;
    my $job;
    $job->{$_} = $options{$_}
        foreach ( keys(%options) );
    $job->{'command'   } = "Galacticus.exe ".$parameterFileName."\n";
    $job->{'launchFile'} = $modelDirectoryName.$label.".sh" ;
    $job->{'logFile'   } = $modelDirectoryName.$label.".log";
    $job->{'label'     } = $label;
    return $job;
}

sub launchModel {
    # Launch a set of models.
    my @jobs    = @{shift()};
    my %options = %{shift()};
    &{$Galacticus::Launch::Hooks::moduleHooks{$options{'launcher'}}->{'jobArrayLaunch'}}(\%options,@jobs);
}

sub parameterNames {
    # Return a list of active parameter names.
    my $config = shift();
    my @names;
    # Extract parameters from config file.
    my @parameters = &List::ExtraUtils::as_array($config->{'posteriorSampleSimulation'}->{'modelParameter'});
    # Extract names.
    foreach my $parameter ( @parameters ) {
	push(@names,$parameter->{'name'}->{'value'})
	    if ( $parameter->{'value'} eq "active" );
    }
    return @names;
}

sub parameterMappings {
    # Return a list of active parameter mappings.
    my $config = shift();
    my @mappings;
    # Extract parameters from config file.
    my @parameters = &List::ExtraUtils::as_array($config->{'posteriorSampleSimulation'}->{'modelParameter'});
    # Extract names.
    foreach my $parameter ( @parameters ) {
	push(@mappings,$parameter->{'operatorUnaryMapper'}->{'value'})
	    if ( $parameter->{'value'} eq "active" );
    }
    return @mappings;
}

sub step {
    # Heaviside step function, required for some derived parameter evaluations.
    my $x = shift();
    return $x >= 0.0 ? 1.0 : 0.0;
}

# sub Sample_Models {
#     # Generate a sample of models from the posterior distribution.
#     my $config  =   shift() ;
#     my %options = %{shift()};
#     # Find the work directory.
#     my $workDirectory = $config->{'likelihood'}->{'workDirectory'};
#     # Get a hash of the parameter values.
#     my $xml        = new XML::Simple();
#     my $parameters = $xml->XMLin($config->{'posteriorSampleLikelihood'}->{'baseParameters'}->{'value'});
#     # Get a matrix of sampled states.
#     my $sampleMatrix = &Sample_Matrix($config,\%options);
#     # Run model for each sample.
#     my $sampleDirectory = exists($options{'sampleDirectory'}) ? $options{'sampleDirectory'}."/" : $workDirectory."/posteriorSample/";
#     my @pbsStack;
#     for (my $i=0;$i<$sampleMatrix->dim(1);++$i) {
# 	# Create an output directory.
# 	my $modelDirectory = $sampleDirectory.$i."/";
# 	system("mkdir -p ".$modelDirectory);
# 	my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
# 	# Check if the model has already been run.
# 	unless ( -e $galacticusFileName ) {
# 	    # Convert these values into a parameter array.
# 	    my $currentConfig = clone($config);
# 	    my $newParameters = &Convert_Parameters_To_Galacticus($currentConfig, $sampleMatrix->(:,($i))->list());
# 	    # Increment the random number seed.
# 	    $parameters->{'randomSeed'}->{'value'} += $config->{'likelihood'}->{'threads'};
# 	    # Clone parameters.
# 	    my $currentParameters = clone($parameters);
# 	    # Apply to parameters.
# 	    &Apply_Parameters($currentParameters,$newParameters);
# 	    # Apply any parameters from command line.
# 	    &applyCommandLineParameters($currentParameters,\%options);
# 	    # Apply any other parameter modifications requested by the caller.
# 	    &{$options{'parametersModifier'}}($currentParameters)
# 		if ( exists($options{'parametersModifier'}) );
# 	    # Specify the output file name.
# 	    $currentParameters->{'outputFileName'}->{'value'} = $galacticusFileName;
# 	    # Write the modified parameters to file.
# 	    &Output($currentParameters,$modelDirectory."parameters.xml");
# 	    # Construct the tasks to perform.
# 	    my $command;
# 	    if ( exists($config->{'likelihood'}->{'environment'}) ) {
# 		foreach ( &List::ExtraUtils::as_array($config->{'likelihood'}->{'environment'}) ) {
# 		    $command .= "export ".$_."\n";
# 		}
# 	    }
# 	    $command .= "ulimit -t unlimited\n";
# 	    $command .= "ulimit -c unlimited\n";
# 	    $command .= "export OMP_NUM_THREADS=".$config->{'likelihood'}->{'threads'}."\n";
# 	    $command .= "mpirun --bynode -np 1 ".(exists($config->{'likelihood'}->{'executable'}) ? $config->{'likelihood'}->{'executable'} : "Galacticus.exe")." ".$modelDirectory."parameters.xml\n";	    
# 	    # Create a PBS job.
# 	    my %job =
# 		(
# 		 launchFile => $modelDirectory."/launch.pbs",
# 		 label      => $config->{'likelihood'}->{'name'}."_ppc".$i,
# 		 logFile    => $modelDirectory."/launch.log",
# 		 command    => $command
# 		);
# 	    foreach ( 'ppn', 'walltime', 'memory' ) {
# 		$job{$_} = $options{$_}
# 		if ( exists($options{$_}) );
# 	    }
# 	    # Queue the calculation.
# 	    push(@pbsStack,\%job);
# 	}
#     }
#     # Send jobs to PBS.
#     &Galacticus::Launch::PBS::SubmitJobs(\%options,@pbsStack)
#      	if ( scalar(@pbsStack) > 0 );
#     # Return the number of models sampled.
#     return ($sampleMatrix->dim(1), $sampleDirectory);
# }

1;
