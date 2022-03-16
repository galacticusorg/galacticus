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
	 base         =>
	     [
	      "haloMassFunctionBase_VSMDPL_z0.000.xml"                 ,
	      "haloMassFunctionBase_VSMDPL_z0.490.xml"                 ,
	      "haloMassFunctionBase_VSMDPL_z0.990.xml"                 ,
	      "haloMassFunctionBase_VSMDPL_z2.030.xml"                 ,
	      "haloMassFunctionBase_VSMDPL_z3.040.xml"                 ,
	      "haloMassFunctionBase_SMDPL_z0.000.xml"                  ,
	      "haloMassFunctionBase_SMDPL_z0.505.xml"                  ,
	      "haloMassFunctionBase_SMDPL_z1.000.xml"                  ,
	      "haloMassFunctionBase_SMDPL_z2.021.xml"                  ,
	      "haloMassFunctionBase_SMDPL_z3.032.xml"                  ,
	      "haloMassFunctionBase_MDPL2_z0.000.xml"                  ,
	      "haloMassFunctionBase_MDPL2_z0.490.xml"                  ,
	      "haloMassFunctionBase_MDPL2_z0.987.xml"                  ,
	      "haloMassFunctionBase_MDPL2_z2.028.xml"                  ,
	      "haloMassFunctionBase_MDPL2_z3.127.xml"                  ,
	      "haloMassFunctionBase_BigMDPL_z0.000.xml"                ,
	      "haloMassFunctionBase_BigMDPL_z0.492.xml"                ,
	      "haloMassFunctionBase_BigMDPL_z1.000.xml"                ,
	      "haloMassFunctionBase_BigMDPL_z2.145.xml"                ,
	      "haloMassFunctionBase_HugeMDPL_z0.000.xml"               ,
	      "haloMassFunctionBase_HugeMDPL_z0.490.xml"               ,
	      "haloMassFunctionBase_HugeMDPL_z0.987.xml"               ,
	      "haloMassFunctionBase_HugeMDPL_z2.028.xml"               ,
              "haloMassFunctionBase_MilkyWay_Halo014_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo247_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo327_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo414_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo460_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo530_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo569_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo628_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo749_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo8247_z0.000.xml"      ,
              "haloMassFunctionBase_MilkyWay_Halo852_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo925_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo939_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo9829_z0.000.xml"      ,
              "haloMassFunctionBase_MilkyWay_Halo023_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo268_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo349_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo415_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo469_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo558_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo570_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo641_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo797_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo825_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo878_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo926_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo967_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo990_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo119_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo288_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo374_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo416_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo490_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo567_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo573_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo675_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo800_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo829_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo881_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo937_z0.000.xml"       ,
              "haloMassFunctionBase_MilkyWay_Halo9749_z0.000.xml"      ,
              "haloMassFunctionBase_MilkyWay_WDM1_Halo416_z0.000.xml"  ,
              "haloMassFunctionBase_MilkyWay_WDM5_Halo416_z0.000.xml"  ,
              "haloMassFunctionBase_MilkyWay_WDM10_Halo416_z0.000.xml" ,
	     ]                                                         ,
	 suffix       =>
	     [
	      "VSMDPL_z0.000"                                          ,
	      "VSMDPL_z0.490"                                          ,
	      "VSMDPL_z0.990"                                          ,
	      "VSMDPL_z2.030"                                          ,
	      "VSMDPL_z3.040"                                          ,
	      "SMDPL_z0.000"                                           ,
	      "SMDPL_z0.505"                                           ,
	      "SMDPL_z1.000"                                           ,
	      "SMDPL_z2.021"                                           ,
	      "SMDPL_z3.032"                                           ,
	      "MDPL2_z0.000"                                           ,
	      "MDPL2_z0.490"                                           ,
	      "MDPL2_z0.987"                                           ,
	      "MDPL2_z2.028"                                           ,
	      "MDPL2_z3.127"                                           ,
	      "BigMDPL_z0.000"                                         ,
	      "BigMDPL_z0.492"                                         ,
	      "BigMDPL_z1.000"                                         ,
	      "BigMDPL_z2.145"                                         ,
	      "HugeMDPL_z0.000"                                        ,
	      "HugeMDPL_z0.490"                                        ,
	      "HugeMDPL_z0.987"                                        ,
	      "HugeMDPL_z2.028"                                        ,
              "MilkyWay_Halo014_z0.000"                                ,
              "MilkyWay_Halo247_z0.000"                                ,
              "MilkyWay_Halo327_z0.000"                                ,
              "MilkyWay_Halo414_z0.000"                                ,
              "MilkyWay_Halo460_z0.000"                                ,
              "MilkyWay_Halo530_z0.000"                                ,
              "MilkyWay_Halo569_z0.000"                                ,
              "MilkyWay_Halo628_z0.000"                                ,
              "MilkyWay_Halo749_z0.000"                                ,
              "MilkyWay_Halo8247_z0.000"                               ,
              "MilkyWay_Halo852_z0.000"                                ,
              "MilkyWay_Halo925_z0.000"                                ,
              "MilkyWay_Halo939_z0.000"                                ,
              "MilkyWay_Halo9829_z0.000"                               ,
              "MilkyWay_Halo023_z0.000"                                ,
              "MilkyWay_Halo268_z0.000"                                ,
              "MilkyWay_Halo349_z0.000"                                ,
              "MilkyWay_Halo415_z0.000"                                ,
              "MilkyWay_Halo469_z0.000"                                ,
              "MilkyWay_Halo558_z0.000"                                ,
              "MilkyWay_Halo570_z0.000"                                ,
              "MilkyWay_Halo641_z0.000"                                ,
              "MilkyWay_Halo797_z0.000"                                ,
              "MilkyWay_Halo825_z0.000"                                ,
              "MilkyWay_Halo878_z0.000"                                ,
              "MilkyWay_Halo926_z0.000"                                ,
              "MilkyWay_Halo967_z0.000"                                ,
              "MilkyWay_Halo990_z0.000"                                ,
              "MilkyWay_Halo119_z0.000"                                ,
              "MilkyWay_Halo288_z0.000"                                ,
              "MilkyWay_Halo374_z0.000"                                ,
              "MilkyWay_Halo416_z0.000"                                ,
              "MilkyWay_Halo490_z0.000"                                ,
              "MilkyWay_Halo567_z0.000"                                ,
              "MilkyWay_Halo573_z0.000"                                ,
              "MilkyWay_Halo675_z0.000"                                ,
              "MilkyWay_Halo800_z0.000"                                ,
              "MilkyWay_Halo829_z0.000"                                ,
              "MilkyWay_Halo881_z0.000"                                ,
              "MilkyWay_Halo937_z0.000"                                ,
              "MilkyWay_Halo9749_z0.000"                               ,
              "MilkyWay_WDM1_Halo416_z0.000"                           ,
              "MilkyWay_WDM5_Halo416_z0.000"                           ,
              "MilkyWay_WDM10_Halo416_z0.000"                          ,
	     ]                                                         ,
	 ppn          => 16                                            ,
	 nodes        =>  2                                            ,
	 postProcess  => $pipelinePath."haloMassFunctionPostProcess.pl"
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
    my $parser = XML::LibXML->new();
    my @bases;
    foreach my $baseFileName ( @{$task->{'base'}} ) {
	my $dom = $parser->load_xml(location => $pipelinePath.$baseFileName);
	$parser->process_xincludes($dom);
	push(@bases,$xml->XMLin($dom->serialize()));
    }
   
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
