#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use File::Copy;
use Data::Dumper;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Galacticus::Constraints::Parameters;

# Construct halo mass function data from a variety of cosmological N-body simulations.
# Andrew Benson (14-October-2020)

# Get command line options.
my %options =
    (
     submitSleepDuration =>  5,
     waitSleepDuration   => 30,
     pbsJobMaximum       => 64
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Validate required parameters are present.
die('simulationDataPath is required but is not present')
    unless ( exists($options{'simulationDataPath'}) );

# Define simulations to process.
my @simulations =
(
 # {
 #     label               => "VSMDPL",
 #     description         => "Halo mass function for non-backsplash z=0 halos from the VSMDPL simulation.",
 #     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
 #     simulationURL       => "https://www.cosmosim.org/cms/simulations/vsmdpl/",
 #     hubbleConstant      => 0.6777,
 #     massParticle        => 6.2e6
 # },
 # {
 #     label               => "SMDPL",
 #     description         => "Halo mass function for non-backsplash z=0 halos from the SMDPL simulation.",
 #     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
 #     simulationURL       => "https://www.cosmosim.org/cms/simulations/smdpl/",
 #     hubbleConstant      => 0.6777,
 #     massParticle        => 9.63e7
 # },
 {
     label               => "MDPL2",
     description         => "Halo mass function for non-backsplash z=0 halos from the MDPL2 simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/mdpl2/",
     hubbleConstant      => 0.6777,
     massParticle        => 1.51e9
 },
 # {
 #     label               => "BigMDPL",
 #     description         => "Halo mass function for non-backsplash z=0 halos from the BigMDPL simulation.",
 #     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
 #     simulationURL       => "https://www.cosmosim.org/cms/simulations/bigmdpl/",
 #     hubbleConstant      => 0.6777,
 #     massParticle        => 2.359e10
 # },
 {
     label               => "HugeMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the HugeMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/hugemdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 7.9e10
 }
);

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Get an XML object.
my $xml = new XML::Simple();

# Iterate over simulations to identify always-isolated halos.
my @jobs;
foreach my $simulation ( @simulations ) {
    # Convert particle mass to Solar masses.
    $simulation->{'massParticle'} /= $simulation->{'hubbleConstant'};
    # Determine minimum and maximum masses for the mass function.
    $simulation->{'massMinimum'} = 10.0**(int(log($simulation->{'massParticle'})/log(10.0))+2);
    $simulation->{'massMaximum'} = 1.0e16;
    # Construct the simulation path.
    $simulation->{'path'} = $options{'simulationDataPath'};
    $simulation->{'path'} .= "/"
	unless ( $simulation->{'path'} =~ m/\/$/ );
    $simulation->{'path'} .= $simulation->{'label'}."/";
    # Identify always isolated halos at z=0 and export these to IRATE format files.
    ## Parse the base parameters.
    my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/haloMassFunctionIdentifyAlwaysIsolated.xml");
    ## Iterate over subvolumes.
    for(my $i=0;$i<10;++$i) {
	for(my $j=0;$j<10;++$j) {
	    for(my $k=0;$k<10;++$k) {
		# Skip if the file exists.
		next
		    if ( -e $simulation->{'path'}."alwaysIsolated_z0.000_subVolume".$i."_".$j."_".$k.".hdf5" );
		# Modify file names.
		$parameters->{'nbodyImporter'}                              ->{'fileName'}->{'value'} = $simulation->{'path'}."tree_"                          .$i."_".$j."_".$k.".dat" ;
		$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[3]->{'fileName'}->{'value'} = $simulation->{'path'}."alwaysIsolated_z0.000_subVolume".$i."_".$j."_".$k.".hdf5";
		# Write parmeter file.
		my $parameterFileName = $simulation->{'path'}."identifyAlwaysIsolated_".$i."_".$j."_".$k.".xml";
		open(my $outputFile,">",$parameterFileName);
		print $outputFile $xml->XMLout($parameters, RootName => "parameters");
		close($outputFile);
		# Generate a job.
		my $job;
		$job->{'command'   } =
		    "Galacticus.exe ".$parameterFileName;
		$job->{'launchFile'} = $simulation->{'path'}."identifyAlwaysIsolated_".$i."_".$j."_".$k.".sh" ;
		$job->{'logFile'   } = $simulation->{'path'}."identifyAlwaysIsolated_".$i."_".$j."_".$k.".log";
		$job->{'label'     } =                 "identifyAlwaysIsolated_".$i."_".$j."_".$k             ;
		$job->{'ppn'       } = 1;
		$job->{'nodes'     } = 1;
		$job->{'mpi'       } = "yes";
		push(@jobs,$job);
	    }
	}
    }
}
&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs)
    if ( scalar(@jobs) > 0 );

# Iterate over simulations to construct the mass functions.
my @massFunctionJobs;
foreach my $simulation ( @simulations ) {
    ## Parse the base parameters.
    my $massFunctionParameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/haloMassFunctionCompute.xml");
    ## Iterate over subvolumes.
    my @nbodyImporters;
    for(my $i=0;$i<10;++$i) {
	for(my $j=0;$j<10;++$j) {
	    for(my $k=0;$k<10;++$k) {
		# Add an importer for this subvolume.
		push(
		    @nbodyImporters,
		    {
			value    => "IRATE"                                                                              ,
			fileName => {value => $simulation->{'path'}."alwaysIsolated_z0.000_subVolume".$i."_".$j."_".$k.".hdf5"},
			snapshot => {value => "1"}
		    }
		    );
	    }
	}
    }
    ## Compute the mass function.
    unless ( -e $simulation->{'path'}."haloMassFunction_z0.000.hdf5" ) {
	## Modify parameters.
	@{$massFunctionParameters->{'nbodyImporter'     }->{'nbodyImporter'}} = @nbodyImporters;
	$massFunctionParameters  ->{'outputFileName'}->{'value'}  = $simulation->{'path'}."haloMassFunction_z0.000.hdf5";
$massFunctionParameters  ->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'values'             }->{'value'} = $simulation->{'massParticle'       };
	$massFunctionParameters  ->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'description'        }->{'value'} = $simulation->{'description'        };
	$massFunctionParameters  ->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'simulationReference'}->{'value'} = $simulation->{'simulationReference'};
	$massFunctionParameters  ->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'simulationURL'      }->{'value'} = $simulation->{'simulationURL'      };
	$massFunctionParameters  ->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'massMinimum'        }->{'value'} = $simulation->{'massMinimum'        };
	$massFunctionParameters  ->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'massMaximum'        }->{'value'} = $simulation->{'massMaximum'        };
	## Write the parameter file.
	my $parameterFileName = $simulation->{'path'}."haloMassFunction_z0.000.xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($massFunctionParameters, RootName => "parameters");
	close($outputFile);
	## Construct the job.
	my $job;
	$job->{'command'   } =
	    "Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $simulation->{'path'}."haloMassFunction_z0.000.sh" ;
	$job->{'logFile'   } = $simulation->{'path'}."haloMassFunction_z0.000.log";
	$job->{'label'     } =                       "haloMassFunction_z0.000"    ;
	$job->{'ppn'       } = 16;
	$job->{'nodes'     } = 1;
	$job->{'mpi'       } = "no";
	push(@massFunctionJobs,$job);
    }
}
&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@massFunctionJobs)
    if ( scalar(@massFunctionJobs) > 0 );

# Move the resulting mass functions.
foreach my $simulation ( @simulations ) {
    copy($simulation->{'path'}."haloMassFunction_z0.000:MPI0000.hdf5",$ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$simulation->{'label'}."_z0.000.hdf5");
}

exit 0;
