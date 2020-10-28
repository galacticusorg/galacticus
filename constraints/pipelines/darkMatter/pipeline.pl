#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use Cwd;
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Galacticus::Constraints::Parameters;

# A constraint pipeline which constrains various dark matter physics models in Galacticus.
# Andrew Benson (22-September-2020)

# Get command line options.
my %options;
$options{'outputDirectory'} = getcwd()."/pipeline";
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Specify the pipeline path.
my $pipelinePath = $ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/";

# Make the output directory.
mkdir($options{'outputDirectory'});

# Get an XML parser.
my $xml = new XML::Simple;

# Define tasks in the pipeline.
my @tasks =
    (
     {
	 label       => "haloMassFunction"                            ,
	 config      => "haloMassFunctionConfig.xml"                  ,
	 base        => "haloMassFunctionBase.xml"                    ,
	 ppn         => 16                                            ,
	 nodes       =>  1                                            ,
	 postProcess => $pipelinePath."haloMassFunctionPostProcess.pl",
	 parameterMap =>
	 {
	     "haloMassFunctionMethod::haloMassFunctionMethod::haloMassFunctionConditioned::a"               => "haloMassFunctionMethod::a"            ,
	     "haloMassFunctionMethod::haloMassFunctionMethod::haloMassFunctionConditioned::p"               => "haloMassFunctionMethod::p"            ,
	     "haloMassFunctionMethod::haloMassFunctionMethod::haloMassFunctionConditioned::normalization"   => "haloMassFunctionMethod::normalization",
	     "haloMassFunctionMethod::haloMassFunctionMethod::haloMassFunctionUnconditioned::a"             => "haloMassFunctionMethod::a"            ,
	     "haloMassFunctionMethod::haloMassFunctionMethod::haloMassFunctionUnconditioned::p"             => "haloMassFunctionMethod::p"            ,
	     "haloMassFunctionMethod::haloMassFunctionMethod::haloMassFunctionUnconditioned::normalization" => "haloMassFunctionMethod::normalization",
	 }
     },
     {
	 label       => "final"    ,
	 base        => "final.xml"
     }
    );

# Initialize the list of parameters whose values have been determined.
my %parametersDetermined;

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Iterate through tasks.
foreach my $task ( @tasks ) {
    # Parse the config and base parameter files for this task.
    my $config = $xml->XMLin($pipelinePath.$task->{'config'})
	if ( exists($task->{'config'}) );
    my $base   = $xml->XMLin($pipelinePath.$task->{'base'  });
    
    # Copy in current best parameters.
    &applyParameters($base,$task,\%parametersDetermined);

    # Change output paths.
    $base->{'galacticusOutputFileName'}->{'value'} = $options{'outputDirectory'}."/".$task->{'label'}.".hdf5" ;
    if ( defined($config) ) {
	$config->{'galacticusOutputFileName'       }                 ->{'value'} = $options{'outputDirectory'}."/".$task->{'label'}.".hdf5" ;
	$config->{'posteriorSampleSimulationMethod'}->{'logFileRoot'}->{'value'} = $options{'outputDirectory'}."/".$task->{'label'}."Chains";
    }
    
    # Write out modified config and base parameter files.
    my $baseFileName = $options{'outputDirectory'}."/".$task->{'base'};
    open(my $baseFile,">".$baseFileName);
    print $baseFile $xml->XMLout($base,RootName => "parameters");
    close($baseFile);
    if ( defined($config) ) {
	my $configFileName = $options{'outputDirectory'}."/".$task->{'config'};
	open(my $configFile,">".$configFileName);
	print $configFile $xml->XMLout($config,RootName => "parameters");
	close($configFile);

	# Generate the job.
	if ( defined($config) ) {
	    my $job;
	    $job->{'command'   } = "Galacticus.exe ".$configFileName;
	    $job->{'launchFile'} = $options{'outputDirectory'}."/".$task->{'label'}.".sh" ;
	    $job->{'logFile'   } = $options{'outputDirectory'}."/".$task->{'label'}.".log";
	    $job->{'label'     } = "darkMatterPipeline".ucfirst($task->{'label'});
	    $job->{'ppn'       } = $task->{'ppn'  };
	    $job->{'nodes'     } = $task->{'nodes'};
	    $job->{'mpi'       } = "yes";
	    my @jobs = ( $job );
	    
	    # Run the job.
	    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs)
	    	unless ( -e $options{'outputDirectory'}."/".$task->{'label'}."Chains_0000.log" );
	    
	    # Extract results.
	    my %posteriorOptions =
		(
		);
	    (my $parametersMaximumLikelihood) = &Galacticus::Constraints::Parameters::maximumPosteriorParameterVector($config,\%posteriorOptions);
	    my @parameterNames                = &Galacticus::Constraints::Parameters::parameterNames                 ($config                   );
	    for(my $i=0;$i<scalar(@parameterNames);++$i) {
		$parametersDetermined{$parameterNames[$i]} = $parametersMaximumLikelihood->(($i));
	    }
	    
	    # Output the determined parameters.
	    open(my $resultsFile,">".$options{'outputDirectory'}."/results.txt");
	    foreach my $parameterName ( sort(keys(%parametersDetermined)) ) {
		print $resultsFile $parameterName."\t".$parametersDetermined{$parameterName}."\n";
	    }
	    close($resultsFile);

	    # Copy in new best parameters.
	    &applyParameters($base,$task,\%parametersDetermined);
	    open(my $baseFile,">".$baseFileName);
	    print $baseFile $xml->XMLout($base,RootName => "parameters");
	    close($baseFile);

	    # Postprocess the results.
	    if ( exists($task->{'postProcess'}) ) {
		system($task->{'postProcess'}." ".$options{'outputDirectory'});
	    }
	    
	}
	
    }

}

exit 0;

sub applyParameters {
    # Apply determined parameters to the base parameters.
    my $base = shift();
    my $task = shift();    
    my %parametersDetermined = %{shift()};
    if ( exists($task->{'parameterMap'}) ) {
	# Apply parameters according to the task's map.
	foreach my $parameterName ( keys(%{$task->{'parameterMap'}}) ) {
	    if ( exists($parametersDetermined{$task->{'parameterMap'}->{$parameterName}}) ) {
		(my $parameter, my $valueIndex) = &Galacticus::Constraints::Parameters::parameterFind($base,$parameterName);
		&Galacticus::Constraints::Parameters::parameterValueSet($parameter,$valueIndex,sclr($parametersDetermined{$task->{'parameterMap'}->{$parameterName}}));
	    }
	}
    } else {
	# Apply all parameters directly.
	foreach my $parameterName ( keys(%parametersDetermined) ) {
	    (my $parameter, my $valueIndex) = &Galacticus::Constraints::Parameters::parameterFind($base,$parameterName);
	    &Galacticus::Constraints::Parameters::parameterValueSet($parameter,$valueIndex,sclr($parametersDetermined{$parameterName}));
	}
    }
}
