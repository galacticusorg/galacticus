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
use List::ExtraUtils;
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
$options{'outputDirectory'       } = getcwd()."/pipeline";
$options{'updateResults'         } = "yes";
$options{'maximum'               } = "posterior";
$options{'writeParameters'       } = "yes";
$options{'callgrind'             } = "no";
$options{'generateContent'       } = "yes";
$options{'haloMassFunction:nodes'} = 16;
$options{'haloMassFunction:ppn'  } = 32;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Specify the pipeline path.
my $pipelinePath = $ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/";

# Make the output directory.
mkdir($options{'outputDirectory'});

# Write our own run command to file.
open(my $command,">",$options{'outputDirectory'}."/command.txt");
print $command join(" ",$0,@ARGV)."\n";
close($command);

# Get an XML parser.
my $xml = new XML::Simple();

# Define tasks in the pipeline.
my @tasks =
    (
     {
	 label       => "haloMassFunction"                ,
	 ppn         => $options{'haloMassFunction:ppn'  },
	 nodes       => $options{'haloMassFunction:nodes'},
	 postprocess => 1
     },
     # {
     # 	 label        => "progenitorMassFunction"                            ,
     # 	 config       => "progenitorMassFunctionConfig.xml"                  ,
     # 	 base         =>
     # 	     [
     # 	      "progenitorMassFunctionBaseHugeMDPL.xml"        ,
     # 	      "progenitorMassFunctionBaseBigMDPL.xml"         ,
     # 	      "progenitorMassFunctionBaseMDPL2.xml"           ,
     # 	      "progenitorMassFunctionBaseSMDPL.xml"           ,
     # 	      "progenitorMassFunctionBaseVSMDPL.xml"          ,
     # 	      "progenitorMassFunctionBaseCaterpillar_LX12.xml",
     # 	      "progenitorMassFunctionBaseCaterpillar_LX13.xml",
     # 	      "progenitorMassFunctionBaseCaterpillar_LX14.xml"
     # 	     ],
     # 	 ppn          => 16                                                  ,
     # 	 nodes        =>  1                                                  ,
     # 	 postProcess  => $pipelinePath."progenitorMassFunctionPostProcess.pl"
     # },
     # {
     # 	 label        => "spinConcentration"          ,
     # 	 config       => "spinConcentrationConfig.xml",
     # 	 base         =>
     # 	     [
     # 	      "spinConcentrationBaseHugeMDPL.xml",
     # 	      "spinConcentrationBaseBigMDPL.xml" ,
     # 	      "spinConcentrationBaseMDPL2.xml"   ,
     # 	      "spinConcentrationBaseSMDPL.xml"   ,
     # 	      "spinConcentrationBaseVSMDPL.xml"
     # 	     ],
     # 	 ppn          => 16                                             ,
     # 	 nodes        =>  1                                             ,
     # 	 postProcess  => $pipelinePath."spinConcentrationPostProcess.pl"
     # },
     # {
     # 	 label       => "final"    ,
     # 	 base        => [ "final.xml" ]
     # }
    );

# Initialize the list of parameters whose values have been determined.
my %parametersDetermined;

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Iterate through tasks.
foreach my $task ( @tasks ) {
    print "Begin task: ".$task->{'label'}."\n";
    # Generate all required files for this task.
    if ( $options{'generateContent'} eq "yes" ) {
	print "  Generating content...\n";
	system($pipelinePath.$task->{'label'}."GenerateContent.pl --pipelinePath ".$pipelinePath.&Galacticus::Options::Serialize_Options(\%options));
	print "  ...done\n";
    }

    # Construct the post-processing command and write to file.
    my $postProcessCommand = $pipelinePath.$task->{'label'}."PostProcess.pl --pipelinePath ".$pipelinePath.&Galacticus::Options::Serialize_Options(\%options);
    open(my $postProcessCommandFile,">",$options{'outputDirectory'}."/postProcessCommand.txt");
    print $postProcessCommandFile $postProcessCommand."\n";
    close($postProcessCommandFile);
    
    # Parse the config and base parameter files for this task.
    print "  Processing base parameter files...\n";
    my $configFileName = $options{'outputDirectory'}."/".$task->{'label'}."Config.xml";
    my $config = $xml->XMLin($configFileName);
    my $parser = XML::LibXML->new();
    foreach my $likelihoodModel ( &List::ExtraUtils::as_array($config->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}) ) {
	my $likelihoodModel_ = $likelihoodModel;
	while ( exists($likelihoodModel_->{'posteriorSampleLikelihood'}) ) {
	    $likelihoodModel_ = $likelihoodModel_->{'posteriorSampleLikelihood'};
	}
	my @redshifts = split(" ",$likelihoodModel_->{'redshifts'}->{'value'});
	foreach my $redshift ( @redshifts ) {
	    (my $baseFileName = $likelihoodModel_->{'baseParametersFileName'}->{'value'}) =~ s/_z[\d\.]+\d/_z$redshift/;
	    my $dom          = $parser->load_xml(location => $baseFileName);
	    $parser->process_xincludes($dom);
	    push(@{$likelihoodModel_->{'baseParameters'}},$xml->XMLin($dom->serialize()));
	}
    }
    
    # Copy in current best parameters.
    &applyParameters($config,\%parametersDetermined);
    print "  ...done\n";
    
    # Write out modified config and base parameter files.
    if ( $options{'writeParameters'} eq "yes" ) {
	print "  Writing updated parameter files...\n";
	&writeParameters($config);
	print "  ...done\n";
    }

    # Create the MCMC job.
    my $job;
    $job->{'command'   } = ($options{'callgrind'} eq "yes" ? "--output-filename ".$options{'outputDirectory'}."/".$task->{'label'}.".vlog valgrind --tool=callgrind " : "")."Galacticus.exe ".$configFileName;
    $job->{'launchFile'} = $options{'outputDirectory'}."/".$task->{'label'}.".sh" ;
    $job->{'logFile'   } = $options{'outputDirectory'}."/".$task->{'label'}.".log";
    $job->{'label'     } = "darkMatterPipeline".ucfirst($task->{'label'});
    $job->{'ppn'       } = $task->{'ppn'  };
    $job->{'nodes'     } = $task->{'nodes'};
    $job->{'mpi'       } = "yes";
    my @jobs = ( $job );
    # Run the job.
    unless ( -e $config->{'posteriorSampleSimulation'}->{'logFileRoot'}->{'value'}."_0000.log" ) {
	my $timeBegin = DateTime->now();
	print "  Running MCMC for '".$task->{'label'}."'\n";
	print "  \tbegin at ".$timeBegin->datetime()."\n";
	&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs);
	my $timeEnd   = DateTime->now();
	print "  \t done at ".$timeEnd  ->datetime()."\n";
	print "  ...done\n";
    }

    # Get results.
    if ( $options{'updateResults'} eq "yes" ) {
	# Extract results.
	print "  Extracting maximum ".$options{'maximum'}." parameter vector...\n";
	my %posteriorOptions =
	    (
	    );
	my $parametersMaximumLikelihood;
	if ( $options{'maximum'} eq "likelihood" ) {
	    ($parametersMaximumLikelihood) = &Galacticus::Constraints::Parameters::maximumLikelihoodParameterVector($config,\%posteriorOptions);
	} else {
	    ($parametersMaximumLikelihood) = &Galacticus::Constraints::Parameters::maximumPosteriorParameterVector ($config,\%posteriorOptions);
	}
	my @parameterNames                 = &Galacticus::Constraints::Parameters::parameterNames                  ($config                   );
	for(my $i=0;$i<scalar(@parameterNames);++$i) {
	    $parametersDetermined{$parameterNames[$i]} = $parametersMaximumLikelihood->(($i));
	}
	print "  ...done\n";

	# Output the determined parameters.
	print "  Outputing maximum likelihood parameters...\n";
	open(my $resultsFile,">".$options{'outputDirectory'}."/results.txt");
	foreach my $parameterName ( sort(keys(%parametersDetermined)) ) {
	    print $resultsFile $parameterName."\t".$parametersDetermined{$parameterName}."\n";
	}
	close($resultsFile);
	print "  ...done\n";
    } else {
	# Not updating results - attempt to read the prior results.
	die("No results.txt file exists")
	    unless ( -e $options{'outputDirectory'}."/results.txt" );
	print "  Reading prior maximum likelihood parameters...\n";
	open(my $resultsFile,"<".$options{'outputDirectory'}."/results.txt");
	while ( my $line = <$resultsFile> ) {
	    if ( $line =~ m/([a-zA-Z0-9:\/]+)\s+([\+\-\d\.e]+)/ ) {
		$parametersDetermined{$1} = pdl $2;
	    } else {
		die("Can not parse results.txt file");
	    }
	}
	close($resultsFile);
	print "  ...done\n";
    }

    # Copy in new best parameters.
    if ( $options{'writeParameters'} eq "yes" ) {
	print "  Updating parameter files...\n";
	&applyParameters($config,\%parametersDetermined);
	&writeParameters($config                       );
	print "  ...done\n";
    }

    # Postprocess the results.
    if ( $task->{'postprocess'} ) {
	print "  Postprocessing...\n";
	system($postProcessCommand);
	print "  ...done\n";
    }
}

exit 0;

sub applyParameters {
    # Apply determined parameters to the base parameters.
    my $config               =   shift() ;    
    my %parametersDetermined = %{shift()};
    my $i                    = -1;
    foreach my $likelihoodModel ( &List::ExtraUtils::as_array($config->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}) ) {
	++$i;
	my @parametersMapped;
	if ( exists($config->{'posteriorSampleLikelihood'}->{'parameterMap'}) ) {
	    my @parameterMap = &List::ExtraUtils::as_array($config->{'posteriorSampleLikelihood'}->{'parameterMap'});
	    @parametersMapped = split(" ",$parameterMap[$i]->{'value'});
	}
	my $likelihoodModel_ = $likelihoodModel;
	while ( exists($likelihoodModel_->{'posteriorSampleLikelihood'}) ) {
	    $likelihoodModel_ = $likelihoodModel_->{'posteriorSampleLikelihood'};
	}
	my @redshifts = split(" ",$likelihoodModel_->{'redshifts'}->{'value'});
	my $j = -1;
	foreach my $redshift ( @redshifts ) {
	    ++$j;
	    foreach my $parameterName ( keys(%parametersDetermined) ) {
		next
		    unless ( ! @parametersMapped || grep {$_ eq $parameterName} @parametersMapped );
		(my $parameter, my $valueIndex) = &Galacticus::Constraints::Parameters::parameterFind($likelihoodModel_->{'baseParameters'}->[$j],$parameterName);
		&Galacticus::Constraints::Parameters::parameterValueSet($parameter,$valueIndex,sclr($parametersDetermined{$parameterName}))
		    if ( defined($parameter) );
	    }
	}
    }
}

sub writeParameters {
    # Write updated parameter sets back to file.
    my $config = shift();
    foreach my $likelihoodModel ( &List::ExtraUtils::as_array($config->{'posteriorSampleLikelihood'}->{'posteriorSampleLikelihood'}) ) {
	my $likelihoodModel_ = $likelihoodModel;
	while ( exists($likelihoodModel_->{'posteriorSampleLikelihood'}) ) {
	    $likelihoodModel_ = $likelihoodModel_->{'posteriorSampleLikelihood'};
	}
	my @redshifts = split(" ",$likelihoodModel_->{'redshifts'}->{'value'});
	my $j = -1;
	foreach my $redshift ( @redshifts ) {
	    ++$j;
	    (my $baseFileName = $likelihoodModel_->{'baseParametersFileName'}->{'value'}) =~ s/_z[\d\.]+\d/_z$redshift/;
	    open(my $baseFile,">".$baseFileName);
	    print $baseFile $xml->XMLout($likelihoodModel_->{'baseParameters'}->[$j],RootName => "parameters");
	    close($baseFile);
	}
    }
}
