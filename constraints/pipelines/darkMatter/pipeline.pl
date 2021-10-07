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
use DateTime;
use Scalar::Util qw(reftype);
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
	 label        => "haloMassFunction"                            ,
	 config       => "haloMassFunctionConfig.xml"                  ,
	 base         => [ "haloMassFunctionBase.xml" ]                ,
	 ppn          => 16                                            ,
	 nodes        =>  1                                            ,
	 postProcess  => $pipelinePath."haloMassFunctionPostProcess.pl",
	 parameterMap =>
	 {
	     "haloMassFunction::haloMassFunction::haloMassFunctionConditioned::a"               => "haloMassFunction::a"            ,
	     "haloMassFunction::haloMassFunction::haloMassFunctionConditioned::p"               => "haloMassFunction::p"            ,
	     "haloMassFunction::haloMassFunction::haloMassFunctionConditioned::normalization"   => "haloMassFunction::normalization",
	     "haloMassFunction::haloMassFunction::haloMassFunctionUnconditioned::a"             => "haloMassFunction::a"            ,
	     "haloMassFunction::haloMassFunction::haloMassFunctionUnconditioned::p"             => "haloMassFunction::p"            ,
	     "haloMassFunction::haloMassFunction::haloMassFunctionUnconditioned::normalization" => "haloMassFunction::normalization",
	 }
     },
     {
	 label        => "progenitorMassFunction"                            ,
	 config       => "progenitorMassFunctionConfig.xml"                  ,
	 base         =>
	     [
	      "progenitorMassFunctionBaseHugeMDPL.xml"        ,
	      "progenitorMassFunctionBaseBigMDPL.xml"         ,
	      "progenitorMassFunctionBaseMDPL2.xml"           ,
	      "progenitorMassFunctionBaseSMDPL.xml"           ,
	      "progenitorMassFunctionBaseVSMDPL.xml"          ,
	      "progenitorMassFunctionBaseCaterpillar_LX12.xml",
	      "progenitorMassFunctionBaseCaterpillar_LX13.xml",
	      "progenitorMassFunctionBaseCaterpillar_LX14.xml"
	     ],
	 suffix       =>
	     [
	      "HugeMDPL"        ,
	      "BigMDPL"         ,
	      "MDPL2"           ,
	      "SMDPL"           ,
	      "VSMDPL"          ,
	      "Caterpillar_LX12",
	      "Caterpillar_LX13",
	      "Caterpillar_LX14"
	     ],
	 ppn          => 16                                                  ,
	 nodes        =>  1                                                  ,
	 postProcess  => $pipelinePath."progenitorMassFunctionPostProcess.pl"
     },
     {
	 label        => "spinConcentration"          ,
	 config       => "spinConcentrationConfig.xml",
	 base         =>
	     [
	      "spinConcentrationBaseHugeMDPL.xml",
	      "spinConcentrationBaseBigMDPL.xml" ,
	      "spinConcentrationBaseMDPL2.xml"   ,
	      "spinConcentrationBaseSMDPL.xml"   ,
	      "spinConcentrationBaseVSMDPL.xml"
	     ],
	 suffix       =>
	     [
	      "HugeMDPL",
	      "BigMDPL" ,
	      "MDPL2"   ,
	      "SMDPL"   ,
	      "VSMDPL"
	     ],
	 ppn          => 16                                             ,
	 nodes        =>  1                                             ,
	 postProcess  => $pipelinePath."spinConcentrationPostProcess.pl"
     },
     {
	 label       => "final"    ,
	 base        => [ "final.xml" ]
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
    my $config =      $xml->XMLin($pipelinePath.$task->{'config'})
	if ( exists($task->{'config'}) );
    my @bases  = map {$xml->XMLin($pipelinePath.$_)} @{$task->{'base'}};
    
    # Copy in current best parameters.
    &applyParameters($_,$task,\%parametersDetermined)
	foreach ( @bases );

    # Change output paths.
    my $i = -1;
    foreach my $base ( @bases ) {
	++$i;
	my $suffix = exists($task->{'suffix'}) ? $task->{'suffix'}->[$i] : "";
	$base->{'galacticusOutputFileName'}->{'value'} = $options{'outputDirectory'}."/".$task->{'label'}.$suffix.".hdf5";
    }
    if ( defined($config) ) {
	$config->{'galacticusOutputFileName' }                 ->{'value'} = $options{'outputDirectory'}."/".$task->{'label'}.".hdf5" ;
	$config->{'posteriorSampleSimulation'}->{'logFileRoot'}->{'value'} = $options{'outputDirectory'}."/".$task->{'label'}."Chains";
	my $i = -1;
	# Set base file names in the config file.
	foreach my $base ( @bases ) {
	    ++$i;
	    my $baseFileName = $options{'outputDirectory'}."/".$task->{'base'}->[$i];
	    if      ( exists($config->{'posteriorSampleLikelihood'}->{'baseParametersFileName'         }) ) {
		die("config file has insufficient base parameter file names")
		    if ( $i > 0 );
		$config->{'posteriorSampleLikelihood'}->{'baseParametersFileName'}->{'value'} = $baseFileName;
	    } elsif ( exists($config->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}) ) {
		if ( reftype($config->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}) eq "ARRAY" ) {
		    die("config file has insufficient base parameter file names")
			if ( $i >= scalar(@{$config->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}}) );
		    $config->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}->[$i]->{'baseParametersFileName'}->{'value'} = $baseFileName;
		} else {
		    die("config file has insufficient base parameter file names")
			if ( $i > 0 );
		    $config->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}      ->{'baseParametersFileName'}->{'value'} = $baseFileName;
		}
	    }
	}
    }
    
    # Write out modified config and base parameter files.
    $i = -1;
    foreach my $base ( @bases ) {
	++$i;
	my $baseFileName = $options{'outputDirectory'}."/".$task->{'base'}->[$i];
	open(my $baseFile,">".$baseFileName);
	print $baseFile $xml->XMLout($base,RootName => "parameters");
	close($baseFile);
    }
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
	    unless ( -e $options{'outputDirectory'}."/".$task->{'label'}."Chains_0000.log" ) {
		my $timeBegin = DateTime->now();
		print "Running constraint for '".$task->{'label'}."'\n";
		print "\tbegin at ".$timeBegin->datetime()."\n";
		&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs);
		my $timeEnd   = DateTime->now();
		print "\t done at ".$timeEnd  ->datetime()."\n";
	    }
	    
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
	    $i = -1;
	    foreach my $base ( @bases ) {
		++$i;
		my $baseFileName = $options{'outputDirectory'}."/".$task->{'base'}->[$i];
		&applyParameters($base,$task,\%parametersDetermined);
		open(my $baseFile,">".$baseFileName);
		print $baseFile $xml->XMLout($base,RootName => "parameters");
		close($baseFile);
	    }

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
	    next
		unless ( exists($parametersDetermined{$task->{'parameterMap'}->{$parameterName}}) );
	    (my $parameter, my $valueIndex) = &Galacticus::Constraints::Parameters::parameterFind($base,$parameterName);
	    &Galacticus::Constraints::Parameters::parameterValueSet($parameter,$valueIndex,sclr($parametersDetermined{$task->{'parameterMap'}->{$parameterName}}))
		if ( defined($parameter) );
	}
    } else {
	# Apply all parameters directly.
	foreach my $parameterName ( keys(%parametersDetermined) ) {
	    (my $parameter, my $valueIndex) = &Galacticus::Constraints::Parameters::parameterFind($base,$parameterName);
	    &Galacticus::Constraints::Parameters::parameterValueSet($parameter,$valueIndex,sclr($parametersDetermined{$parameterName}))
		if ( defined($parameter) );
	}
    }
}
