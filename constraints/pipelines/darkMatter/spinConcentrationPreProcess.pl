#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use File::Copy;
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::IO::HDF5;
use List::Util;
use List::ExtraUtils;
use Data::Dumper;
use Storable qw(dclone);
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Galacticus::Constraints::Parameters;

# Construct distribution functions of spins and concentrations from a variety of cosmological N-body simulations.
# Andrew Benson (26-April-2021)

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
 {
     label                              => "VSMDPL",
     description                        => "distribution function for non-backsplash z=0 parent halos from the VSMDPL simulation.",
     simulationReference                => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL                      => "https://www.cosmosim.org/cms/simulations/vsmdpl/",
     hubbleConstant                     => 0.6777,
     massParticle                       => 6.2e6,
     snapshot                           => 150,
     redshift                           => 0.0,
     timeRecent                         => 0.0,
     particleCountMinimum               => 300,
     energyEstimateParticleCountMaximum => 1000.0,
     builder                            => \&cosmoSimBuilder
 },
 {
     label                              => "SMDPL",
     description                        => "function for non-backsplash z=0 parent halos from the SMDPL simulation.",
     simulationReference                => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL                      => "https://www.cosmosim.org/cms/simulations/smdpl/",
     hubbleConstant                     => 0.6777,
     massParticle                       => 9.63e7,
     snapshot                           => 116,
     redshift                           => 0.0,
     timeRecent                         => 0.0,
     particleCountMinimum               => 300,
     energyEstimateParticleCountMaximum => 1000.0,
     builder                            => \&cosmoSimBuilder
 },
 {
     label                              => "MDPL2",
     description                        => "function for non-backsplash z=0 parent halos from the MDPL2 simulation.",
     simulationReference                => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL                      => "https://www.cosmosim.org/cms/simulations/mdpl2/",
     hubbleConstant                     => 0.6777,
     massParticle                       => 1.51e9,
     snapshot                           => 125,
     redshift                           => 0.0,
     timeRecent                         => 0.0,
     particleCountMinimum               => 300,
     energyEstimateParticleCountMaximum => 1000.0,
     builder                            => \&cosmoSimBuilder
 },
 {
     label                              => "BigMDPL",
     description                        => "function for non-backsplash z=0 parent halos from the BigMDPL simulation.",
     simulationReference                => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL                      => "https://www.cosmosim.org/cms/simulations/bigmdpl/",
     hubbleConstant                     => 0.6777,
     massParticle                       => 2.359e10,
     snapshot                           => 79,
     redshift                           => 0.0,
     timeRecent                         => 0.0,
     particleCountMinimum               => 300,
     energyEstimateParticleCountMaximum => 1000.0,
     builder                            => \&cosmoSimBuilder
 },
 {
     label                              => "HugeMDPL",
     description                        => "function for non-backsplash z=0 parent halos from the HugeMDPL simulation.",
     simulationReference                => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL                      => "https://www.cosmosim.org/cms/simulations/hugemdpl/",
     hubbleConstant                     => 0.6777,
     massParticle                       => 7.9e10,
     snapshot                           => 102,
     redshift                           => 0.0,
     timeRecent                         => 0.0,
     particleCountMinimum               => 300,
     energyEstimateParticleCountMaximum => 1000.0,
     builder                            => \&cosmoSimBuilder
 }
);

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Get an XML object.
my $xml = new XML::Simple();

# Iterate over simulations and build jobs.
my $jobs;
foreach my $simulation ( @simulations ) {
    &{$simulation->{'builder'}}($simulation,\%options);
}

# Launch job sequences.
foreach my $jobsList ( @{$jobs} ) {
    next
	unless ( defined($jobsList) );
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@{$jobsList})
	if ( scalar(@{$jobsList}) > 0 );
}

exit 0;

sub cosmoSimBuilder {
    my $simulation =   shift() ;
    my %options    = %{shift()};
    # Convert particle mass to Solar masses.
    $simulation->{'massParticle'} /= $simulation->{'hubbleConstant'};
    # Determine minimum and maximum halo masses to consider.
    $simulation->{'massHaloMinimum'} = 1.0e03*$simulation->{'massParticle'};
    $simulation->{'massHaloMaximum'} = 1.0e16;
    # Construct the simulation path.
    $simulation->{'path'} = $options{'simulationDataPath'};
    $simulation->{'path'} .= "/"
	unless ( $simulation->{'path'} =~ m/\/$/ );
    $simulation->{'path'} .= "CosmoSim/".$simulation->{'label'}."/";
    # Identify always isolated halos and export these to IRATE format files.
    ## Parse the base parameters.
    my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/spinConcentrationIdentifyAlwaysIsolated.xml");
    ## Set the snapshot to select.
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[3]->{'selectedValues'}->{'value'} = $simulation->{'snapshot'};
    ## Iterate over subvolumes.
    for(my $i=0;$i<10;++$i) {
	for(my $j=0;$j<10;++$j) {
	    for(my $k=0;$k<10;++$k) {
		# Skip if the file exists.
		next
		    if ( -e $simulation->{'path'}."alwaysIsolated_spinConcentration_subVolume".$i."_".$j."_".$k.".hdf5" );
		# Modify file names.
		$parameters->{'nbodyImporter'}                              ->{'fileName'}->{'value'} = $simulation->{'path'}."tree_"                                     .$i."_".$j."_".$k.".dat" ;
		$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[6]->{'fileName'}->{'value'} = $simulation->{'path'}."alwaysIsolated_spinConcentration_subVolume".$i."_".$j."_".$k.".hdf5";
		# Write parmeter file.
		my $parameterFileName = $simulation->{'path'}."identifyAlwaysIsolated_spinConcentration_".$i."_".$j."_".$k.".xml";
		open(my $outputFile,">",$parameterFileName);
		print $outputFile $xml->XMLout($parameters, RootName => "parameters");
		close($outputFile);
		# Generate a job.
		my $job;
		$job->{'command'   } =
		    "Galacticus.exe ".$parameterFileName;
		$job->{'launchFile'} = $simulation->{'path'}."identifyAlwaysIsolated_spinConcentration_".$i."_".$j."_".$k.".sh" ;
		$job->{'logFile'   } = $simulation->{'path'}."identifyAlwaysIsolated_spinConcentration_".$i."_".$j."_".$k.".log";
		$job->{'label'     } =                       "identifyAlwaysIsolated_spinConcentration_".$i."_".$j."_".$k       ;
		$job->{'ppn'       } = 1;
		$job->{'nodes'     } = 1;
		$job->{'mpi'       } = "yes";
		push(@{$jobs->[0]},$job);
	    }
	}
    }
    ## Parse the base parameters.
    my $spinDistributionFunctionParameters          = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/spinDistributionFunctionCompute.xml"         );
    my $concentrationDistributionFunctionParameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/concentrationDistributionFunctionCompute.xml");
    ## Iterate over subvolumes.
    my @spinNbodyImporters;
    my @concentrationNbodyImporters;
    for(my $i=0;$i<10;++$i) {
	for(my $j=0;$j<10;++$j) {
	    for(my $k=0;$k<10;++$k) {
		# Add an importer for this subvolume.
		push(
		    @spinNbodyImporters,
		    {
			value      => "IRATE"                                                                                               ,
			fileName   => {value => $simulation->{'path'}."alwaysIsolated_spinConcentration_subVolume".$i."_".$j."_".$k.".hdf5"},
			properties => {value => "massVirial spin"},
			snapshot   => {value => "1"}
		    }
		    );
		push(
		    @concentrationNbodyImporters,
		    {
			value      => "IRATE"                                                                                               ,
			fileName   => {value => $simulation->{'path'}."alwaysIsolated_spinConcentration_subVolume".$i."_".$j."_".$k.".hdf5"},
			properties => {value => "massVirial radiusVirial radiusScale"},
			snapshot   => {value => "1"}
		    }
		    );
	    }
	}
    }
    ## Compute the spin distribution function.
    unless ( -e $simulation->{'path'}."spinDistributionFunctions:MPI0000.hdf5" ) {
	## Modify parameters.
	@{$spinDistributionFunctionParameters->{'nbodyImporter'     }->{'nbodyImporter'}} = @spinNbodyImporters;
	$spinDistributionFunctionParameters  ->{'galacticusOutputFileName'}                                                       ->{'value'} =                  $simulation->{'path'                              }."spinDistributionFunctions.hdf5";
	$spinDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[0]->{'values'             }->{'value'} = sprintf("%12.6e",$simulation->{'massParticle'                      })." ".
	                                                                                                                            sprintf("%12.6e",$simulation->{'redshift'                          })." ".
	                                                                                                                            sprintf("%12.6e",$simulation->{'timeRecent'                        })." ".
	                                                                                                                            sprintf("%12.6e",$simulation->{'massHaloMinimum'                   })." ".
	                                                                                                                            sprintf("%12.6e",$simulation->{'massHaloMaximum'                   })." ".
	                                                                                                                            sprintf("%8i"   ,$simulation->{'particleCountMinimum'              })." ".
	                                                                                                                            sprintf("%12.6e",$simulation->{'energyEstimateParticleCountMaximum'})." ".
	                                                                                                                                             $simulation->{'label'                             }                                 ;
	$spinDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'description'        }->{'value'} = "Spin ".         $simulation->{'description'                       }                                 ;
	$spinDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'simulationReference'}->{'value'} =                  $simulation->{'simulationReference'               }                                 ;
	$spinDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'simulationURL'      }->{'value'} =                  $simulation->{'simulationURL'                     }                                 ;
	$spinDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massMinimum'        }->{'value'} =                  $simulation->{'massHaloMinimum'                   }                                 ;
	$spinDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massMaximum'        }->{'value'} =                  $simulation->{'massHaloMaximum'                   }                                 ;
	## Write the parameter file.
	my $parameterFileName = $simulation->{'path'}."spinDistributionFunctions.xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($spinDistributionFunctionParameters, RootName => "parameters");
	close($outputFile);
	## Construct the job.
	my $job;
	$job->{'command'   } =
	    "Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $simulation->{'path'}."spinDistributionFunctions.sh" ;
	$job->{'logFile'   } = $simulation->{'path'}."spinDistributionFunctions.log";
	$job->{'label'     } =                       "spinDistributionFunctions"    ;
	$job->{'ppn'       } = 16;
	$job->{'nodes'     } = 1;
	$job->{'mpi'       } = "no";
	$job->{'onCompletion'} =
	{
	    function  => \&copyFile,
	    arguments => [ $simulation->{'path'}."spinDistributionFunctions:MPI0000.hdf5", $ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/spinDistributionFunctions_".$simulation->{'label'}.".hdf5" ]
	};
	push(@{$jobs->[1]},$job);
    }
    ## Compute the concentration distribution functions.
    my $logMassMinimum = int(log10($simulation->{'massHaloMinimum'})*2.0+0.5)/2.0;
    while ( $logMassMinimum < 15.0 ) {
	my $logMassMaximum = $logMassMinimum+0.5;
	my $suffix = "_logM".$logMassMinimum."_".$logMassMaximum;
	unless ( -e $simulation->{'path'}."concentrationDistributionFunction".$suffix.":MPI0000.hdf5" ) {
	    ## Modify parameters.
	    @{$concentrationDistributionFunctionParameters->{'nbodyImporter'     }->{'nbodyImporter'}} = @concentrationNbodyImporters;
	    $concentrationDistributionFunctionParameters  ->{'galacticusOutputFileName'}                                                       ->{'value'} =                  $simulation->{'path'               }."concentrationDistributionFunctions".$suffix.".hdf5";
	    $concentrationDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[0]->{'values'             }->{'value'} = sprintf("%12.6e",$simulation->{'massParticle'       })." ".
	                                                                                                                                         sprintf("%12.6e",$simulation->{'redshift'           })." ".
	                                                                                                                                         sprintf("%12.6e",$simulation->{'timeRecent'         })." ".
	                                                                                                                                         sprintf("%12.6e",10.0**$logMassMinimum              )." ".
	                                                                                                                                         sprintf("%12.6e",10.0**$logMassMaximum              )." ".
	                                                                                                                                                          $simulation->{'label'              };
	    $concentrationDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'description'        }->{'value'} = "Concentration ".$simulation->{'description'        };
	    $concentrationDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'simulationReference'}->{'value'} = $simulation->{'simulationReference'};
	    $concentrationDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'simulationURL'      }->{'value'} = $simulation->{'simulationURL'      };
	    $concentrationDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massMinimum'        }->{'value'} = 10.0**$logMassMinimum;
	    $concentrationDistributionFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massMaximum'        }->{'value'} = 10.0**$logMassMaximum;
	    ## Write the parameter file.
	    my $parameterFileName = $simulation->{'path'}."concentrationDistributionFunctions".$suffix.".xml";
	    open(my $outputFile,">",$parameterFileName);
	    print $outputFile $xml->XMLout($concentrationDistributionFunctionParameters, RootName => "parameters");
	    close($outputFile);
	    ## Construct the job.
	    my $job;
	    $job->{'command'   } =
		"Galacticus.exe ".$parameterFileName;
	    $job->{'launchFile'} = $simulation->{'path'}."concentrationDistributionFunctions".$suffix.".sh" ;
	    $job->{'logFile'   } = $simulation->{'path'}."concentrationDistributionFunctions".$suffix.".log";
	    $job->{'label'     } =                       "concentrationDistributionFunctions".$suffix       ;
	    $job->{'ppn'       } = 16;
	    $job->{'nodes'     } = 1;
	    $job->{'mpi'       } = "no";
	    $job->{'onCompletion'} =
	    {
		function  => \&copyFile,
		arguments => [ $simulation->{'path'}."concentrationDistributionFunctions".$suffix.":MPI0000.hdf5", $ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/concentrationDistributionFunctions_".$simulation->{'label'}.$suffix.".hdf5" ]
	    };
	    push(@{$jobs->[1]},$job);
	}
	$logMassMinimum += 0.5;
    }
}

sub copyFile {
    # Perform a file copy.
    my $from   = shift();
    my $to     = shift();
    my $jobID  = shift();
    my $status = shift();
    if ( $status == 0 ) {
	copy($from,$to);
    } else {
	print "Warning: job '".$jobID."' failed\n";
    }
}
