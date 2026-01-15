#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use File::Copy;
use Data::Dumper;
use DateTime;
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Galacticus::Constraints::Parameters;
use Galacticus::Constraints::Simulations qw(iterate);

# Analyze a variety of cosmological N-body simulations to extract statistics of interest.
# Andrew Benson (14-October-2020)

# Get command line options.
my %options =
    (
     submitSleepDuration =>  5,
     waitSleepDuration   => 30,
     pbsJobMaximum       => 64,
     slurmJobMaximum     => 64,
     ompThreads          => "max",
     analyses            => "haloMassFunction,subhaloStatistics"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Validate required parameters are present.
foreach ( "simulationDataPath", "pipelinePath" ) {
    die('"'.$_.'" is required but is not present')
	unless ( exists($options{$_}) );
}

# Ensure paths are correctly suffixed.
foreach my $path ( 'simulationDataPath', 'pipelinePath' ) {
    $options{$path} .= "/"
	unless ( $options{$path} =~ m/\/$/ );
}

# Determine which steps are to be included bases on the analyses requested.
my %activeAnalyses = map {$_ => 1} split(/,/,$options{'analyses'});
my %activeSteps;
if ( exists($activeAnalyses{'haloMassFunction' }) ) {
    push(@{$activeSteps{"identifyAlwaysIsolated"}},"haloMassFunction" );
    push(@{$activeSteps{"extractHalos"          }},"haloMassFunction" );
    push(@{$activeSteps{"massFunctions"         }},"haloMassFunction" );
}
if ( exists($activeAnalyses{'subhaloStatistics'}) ) {
    push(@{$activeSteps{"identifyAlwaysIsolated"}},"subhaloStatistics");
    push(@{$activeSteps{"extractSubhalos"       }},"subhaloStatistics");
    push(@{$activeSteps{"subhaloFunctions"      }},"subhaloStatistics");
}

# Establish a set of processor functions for handling the analysis of individual simulation suites.
my %suites =
    (
     MDPL     =>
     {
	 analyses =>
	     [
	      "haloMassFunction"
	     ]
     },
     Symphony =>
     {
	 analyses =>
	     [
	      "haloMassFunction" ,
	      "subhaloStatistics"
	     ],
	 steps =>
	 {
	     identifyAlwaysIsolated =>
	     {
		 processParameters => \&symphonyProcessIdentifyAlwaysIsolated
	     },
	     extractHalos           =>
	     {
		 preprocess        =>
		     [
		      \&symphonyPreProcessExtractHalosLocate        ,
		      \&symphonyPreProcessExtractHalosUncontaminated,
		     ],
		 processParameters => \&symphonyProcessExtractHalos,
		 postprocess       =>
		     [
		      \&symphonyPostprocessSelectInSphere      ,
		      \&symphonyPostprocessSelectInICs         ,
		      \&symphonyPostprocessAnalyze             ,
		      \&symphonyPostprocessSetVolume
		     ]
	     },
	     extractSubhalos        =>
	     {
		 preprocess        =>
		     [
		      \&symphonyPreProcessExtractHalosLocate        ,
		      \&symphonyPreProcessExtractHalosUncontaminated,
		     ],
		 processParameters => \&symphonyProcessExtractSubhalos
	     },
	     massFunctions          =>
	     {
		 postprocess       =>
		     [
		      \&symphonyPostProcessMassFunction
		     ]
	     },
	     subhaloFunctions       =>
	     {
		 processParameters => \&symphonyProcessSubhaloFunctions
	     }
	 }
     }
    );

# The COZMIC suite uses the same processing pipeline as the Symphony suite.
$suites{'COZMIC'} = $suites{'Symphony'};

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Determine number of OpenMP threads to use.
my $ompThreads = $options{'ompThreads'} eq "max" ? $queueConfig->{'ppn'} : $options{'ompThreads'};

# Parse the simulations definition file.
my $xml = new XML::Simple();
my $simulations = $xml->XMLin(
    $options{'pipelinePath'}."simulations.xml",
    ForceArray => [ "suite"         , "group"         , "simulation"         , "resolution"           ],
    KeyAttr    => {  suite => "name",  group => "name",  simulation => "name",  resolution  => "name" }
    );

# Get the list of entries to process.
my @entries = &iterate($simulations,\%options, stopAfter => "realization");

# Define steps to perform (if required).
my @stepIDs = (
    "identifyAlwaysIsolated",
    "extractHalos"          ,
    "extractSubhalos"       ,
    "massFunctions"         ,
    "subhaloFunctions"
  );

# Perform steps.
foreach my $stepID ( @stepIDs ) {
    # Skip inactive steps.
    next
	unless ( exists($activeSteps{$stepID}) );
    # Perform any preprocessing.
    {
	my $workDone  =  1;
	my $iteration = -1;
	while ( $workDone ) {
	    my @preprocessingJobs;
	    ++$iteration;
	    $workDone = 0;
	    foreach my $entry ( @entries ) {
		# If a custom preprocess function is defined, call it.
		if (
		    exists(  $suites{$entry->{'suite'}->{'name'}}->{'steps'}->{$stepID}->{'preprocess'} )
		     &&
		    scalar(@{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{$stepID}->{'preprocess'}}) > $iteration
		    ) {
		    $workDone = 1;
		    &{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{$stepID}->{'preprocess'}->[$iteration]}($entry,\@preprocessingJobs,\%options);
		}
	    }
	    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@preprocessingJobs)
		if ( scalar(@preprocessingJobs) > 0 );
	}
    }
    # Call the function to perform the step.
    &{\&{"step".ucfirst($stepID)}}(\@entries,\%suites,\%activeSteps,$queueManager,\%options);
    # Perform any postprocessing.
    {
	my $workDone  =  1;
	my $iteration = -1;
	while ( $workDone ) {
	    my @postprocessingJobs;
	    ++$iteration;
	    $workDone = 0;
	    foreach my $entry ( @entries ) {
		# If a custom postprocess function is defined, call it.
		if (
		    exists(  $suites{$entry->{'suite'}->{'name'}}->{'steps'}->{$stepID}->{'postprocess'} )
		     &&
		    scalar(@{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{$stepID}->{'postprocess'}}) > $iteration
		    ) {
		    $workDone = 1;
		    &{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{$stepID}->{'postprocess'}->[$iteration]}($entry,\@postprocessingJobs,\%options);
		}
	    }
	    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@postprocessingJobs)
		if ( scalar(@postprocessingJobs) > 0 );
	}
    }
}

# Store the resulting data to the datasets repo.
# Iterate over all active steps
foreach my $stepID ( keys(%activeAnalyses) ) {
    # Determine how to store.
    if ( $stepID eq "haloMassFunction" ) {
	# For halo mass functions we simply copy the mass function file.
	# Iterate over all entries.
	foreach my $entry ( @entries ) {
	    # Iterate over expansion factors.
	    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
		copy($entry->{'path'}."haloMassFunction_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5",$ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'resolution'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$entry->{'realization'}."_".$epoch->{'redshiftLabel'}.".hdf5");
	    }
	}
    } elsif ( $stepID eq "subhaloStatistics" ) {
	# For subhalo statistics, we combine statistics across all available realizations.
	# Initialize data structure to accumulate the statistics.
	my $subhaloStatistics;
	# Iterate over all entries.
	foreach my $entry ( @entries ) {
	    # Iterate over expansion factors.
	    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
		# Build a label unique for the set of realizations to which this entry belongs.
		my $label = $entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'resolution'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$epoch->{'redshiftLabel'};
		# Store relevant metadata.
		(my $refShort = $entry->{'group'}->{'metaData'}->{'reference'}) =~ s/^([^\)]+)\s+\((\d+).*\)/($1 $2)/;
		$subhaloStatistics->{$label}->{'redshift'    } = pdl $epoch                         ->{'redshift' }              ;
		$subhaloStatistics->{$label}->{'reference'   } =     $entry->{'group'}->{'metaData'}->{'reference'}              ;
		$subhaloStatistics->{$label}->{'referenceURL'} =     $entry->{'group'}->{'metaData'}->{'url'      }              ;
		$subhaloStatistics->{$label}->{'label'       } =     $entry->{'suite'}              ->{'name'     }." ".$refShort;
		# Commence reading data for this realization.
		my $entryData  = new PDL::IO::HDF5($entry->{'path'}."subhaloFunctions_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5");
		my $simulation = $entryData->group('simulation0001');
		# Count the number of realizations in this set.
		$subhaloStatistics->{$label}->{'countRealizations'} += 1;
		# Accumulate subhalo mass function.
		{
		    my $subhaloMassFunction = $simulation         ->group  ('subhaloMassFunction')       ;
		    my $count               = $subhaloMassFunction->dataset('count'              )->get();
		    my $massFunction        = $subhaloMassFunction->dataset('massFunction'       )->get();
		    my $massRatio           = $subhaloMassFunction->dataset('massRatio'          )->get();
		    (my $massHost)          = $subhaloMassFunction->attrGet('massHost'           )       ;
		    if ( $subhaloStatistics->{$label}->{'countRealizations'} == 1 ) {
			$subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'count'       } = pdl zeros($massFunction);
			$subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massFunction'} = pdl zeros($massFunction);
			$subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massRatio'   } =           $massRatio    ;
			$subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'    } = pdl      [             ];
		    } else {
			die("mass ratios do not align")
			    unless( all($massRatio == $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massRatio'}) );
		    }
		    $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'count'       } += float($count);
		    $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massFunction'} += $massFunction;
		    $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'    } =  $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'}->append($massHost);
		}
		# Accumulate subhalo radial function.
		{
		    my $subhaloRadiusFunction = $simulation           ->group  ('subhaloRadiusFunction')       ;
		    (my $massMinimum)         = $subhaloRadiusFunction->attrGet('massMinimum'          )       ;
		    my $count                 = $subhaloRadiusFunction->dataset('count'                )->get();
		    my $radialDistribution    = $subhaloRadiusFunction->dataset('radialDistribution'   )->get();
		    my $radiusRatio           = $subhaloRadiusFunction->dataset('radiusRatio'          )->get();
		    if ( $subhaloStatistics->{$label}->{'countRealizations'} == 1 ) {
			$subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'count'             } = pdl zeros($radialDistribution);
			$subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'radialDistribution'} = pdl zeros($radialDistribution);
			$subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'radiusRatio'       } =           $radiusRatio        ;
			$subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'massMinimum'       } =           $massMinimum        ;
		    } else {
			die("radius ratios do not align")
			    unless( all($radiusRatio == $subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'radiusRatio'}) );
		    }
		    $subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'count'             } += float($count);
		    $subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'radialDistribution'} += $radialDistribution;
		}
		# Accumulate subhalo velocity maximum function.
		{
		    my $subhaloVelocityMaximumFunction = $simulation                    ->group  ('subhaloVelocityMaximumMeanFunction')       ;
		    my $count                          = $subhaloVelocityMaximumFunction->dataset('count'                             )->get();
		    my $velocityMaximumMean            = $subhaloVelocityMaximumFunction->dataset('velocityMaximumMean'               )->get();
		    my $velocityMaximumMeanError       = $subhaloVelocityMaximumFunction->dataset('velocityMaximumMeanError'          )->get();
		    my $mass                           = $subhaloVelocityMaximumFunction->dataset('mass'                              )->get();
		    if ( $subhaloStatistics->{$label}->{'countRealizations'} == 1 ) {
			$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'count'                      } = pdl zeros($velocityMaximumMean);
			$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'countRealizations'          } = pdl zeros($velocityMaximumMean);
			$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMean'        } = pdl zeros($velocityMaximumMean);
			$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanSquared' } = pdl zeros($velocityMaximumMean);
			$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanVariance'} = pdl zeros($velocityMaximumMean);
			$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'mass'                       } =           $mass                ;
		    } else {
			die("mass ratios do not align")
			    unless( all($mass == $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'mass'}) );
		    }
		    my $nonZero = which($velocityMaximumMean > 0.0);
		    $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'count'                      }             += float($count                   )   ;
		    $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'countRealizations'          }->($nonZero) += 1
			if ( nelem($nonZero) > 0 );
		    $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMean'        }             +=       $velocityMaximumMean         ;
		    $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanSquared' }             +=       $velocityMaximumMean      **2;
		    $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanVariance'}             +=       $velocityMaximumMeanError **2;
		}
	    }
	}
	# Generate the output.
	foreach my $label ( keys(%{$subhaloStatistics}) ) {
	    # Take averages over realizations.
	    $subhaloStatistics->{$label}->{'subhaloMassFunction'           }->{'count'                      }             /= $subhaloStatistics->{$label}                                    ->{'countRealizations'}            ;
	    $subhaloStatistics->{$label}->{'subhaloMassFunction'           }->{'massFunction'               }             /= $subhaloStatistics->{$label}                                    ->{'countRealizations'}            ;
	    $subhaloStatistics->{$label}->{'subhaloRadiusFunction'         }->{'count'                      }             /= $subhaloStatistics->{$label}                                    ->{'countRealizations'}            ;
	    $subhaloStatistics->{$label}->{'subhaloRadiusFunction'         }->{'radialDistribution'         }             /= $subhaloStatistics->{$label}                                    ->{'countRealizations'}            ;
	    my $nonZero = which($subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'countRealizations'} > 0);
	    if ( nelem($nonZero) > 0 ) {
		$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMean'        }->($nonZero) /= $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'countRealizations'}->($nonZero);
		$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanSquared' }->($nonZero) /= $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'countRealizations'}->($nonZero);
		$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanVariance'}->($nonZero) /= $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'countRealizations'}->($nonZero);
	    }
	    # For the mean velocity, compute the mean and variance by treating individual realizations as components of a mixture
	    # distribution (https://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians/16609#16609).
	    $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanError'   }
	    = sqrt(
		+$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanVariance'}
		+$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanSquared' }
		-$subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMean'        }**2
		);
	    # Divide by √N to get the error on the mean.
	    my $nonZeroCount = which($subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'count'} > 0);
	    $subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanError'}->($nonZeroCount) /= sqrt($subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'count'})->($nonZeroCount)
		if ( nelem($nonZeroCount) > 0 );
	    # Output host halo masses.
	    my $hostHaloFile = new PDL::IO::HDF5(">".$ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/hostHaloMasses_".$label.".hdf5");
	    ## Workaround for failure to write length 1 datasets. Simply replicate the mass - we use these only to set a
	    ## distribution of halo masses, so this makes no practical difference.
	    $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'} = $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'}->append($subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'})
		if ( nelem($subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'}) == 1 );
	    my $treeWeights = pdl ones($subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'});
	    $hostHaloFile->dataset('treeRootMass')->set($subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massHost'},unlimited => 1);
	    $hostHaloFile->dataset('treeWeight'  )->set($treeWeights                                                       ,unlimited => 1);
	    # Output the data.
	    my $storeFile = new PDL::IO::HDF5(">".$ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/subhaloDistributions_".$label.".hdf5");
	    # Metadata.
	    $storeFile->attrSet(label        => $subhaloStatistics->{$label}->{'label'       });
	    $storeFile->attrSet(redshift     => $subhaloStatistics->{$label}->{'redshift'    });
	    $storeFile->attrSet(reference    => $subhaloStatistics->{$label}->{'reference'   });
	    $storeFile->attrSet(referenceURL => $subhaloStatistics->{$label}->{'referenceURL'});
	    $storeFile->attrSet(provenance   => "Computed from Symphony Rockstar `tree_?_?_?.dat` files at ".DateTime->now().".");
	    # Subhalo mass function.
	    my $subhaloMassFunction = $storeFile->group("massFunction");
	    $subhaloMassFunction->attrSet(selection => "Subhalos within the host virial radius at z=".sprintf("%5.3f",$subhaloStatistics->{$label}->{'redshift'}).".");
	    $subhaloMassFunction->dataset('massRatio'        )->set(     $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'massRatio' }                                                     ,unlimited => 1);
	    $subhaloMassFunction->dataset('massFunction'     )->set(     $subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'count'     }                                                     ,unlimited => 1);
	    $subhaloMassFunction->dataset('massFunctionError')->set(sqrt($subhaloStatistics->{$label}->{'subhaloMassFunction'}->{'count'     }/$subhaloStatistics->{$label}->{'countRealizations'}),unlimited => 1);
	    # Subhalo radius function.
	    my $subhaloRadiusFunction = $storeFile->group("radialDistribution");
	    $subhaloRadiusFunction->attrSet(selection => "Subhalos within the host virial radius at z=".sprintf("%5.3f",$subhaloStatistics->{$label}->{'redshift'}).".");
	    $subhaloRadiusFunction->dataset('radiusFractional'       )->set(     $subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'radiusRatio'}                                                     ,unlimited => 1);
	    $subhaloRadiusFunction->dataset('radialDistribution'     )->set(     $subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'count'      }                                                     ,unlimited => 1);
	    $subhaloRadiusFunction->dataset('radialDistributionError')->set(sqrt($subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'count'      }/$subhaloStatistics->{$label}->{'countRealizations'}),unlimited => 1);
	    $subhaloRadiusFunction->attrSet(massMinimum => $subhaloStatistics->{$label}->{'subhaloRadiusFunction'}->{'massMinimum'});
	    # Subhalo velocity maximum function.
	    my $subhaloVelocityMaximumFunction = $storeFile->group("velocityMaximum");
	    $subhaloVelocityMaximumFunction->attrSet(selection => "Subhalos within the host virial radius at z=".sprintf("%5.3f",$subhaloStatistics->{$label}->{'redshift'}).".");
	    $subhaloVelocityMaximumFunction->dataset('mass'                    )->set($subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'mass'                    },unlimited => 1);
	    $subhaloVelocityMaximumFunction->dataset('velocityMaximumMean'     )->set($subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMean'     },unlimited => 1);
	    $subhaloVelocityMaximumFunction->dataset('velocityMaximumMeanError')->set($subhaloStatistics->{$label}->{'subhaloVelocityMaximumFunction'}->{'velocityMaximumMeanError'},unlimited => 1);
	}
    }
}

exit 0;

sub stepIdentifyAlwaysIsolated {
    # Identify always isolated halos.
    my @entries      = @{shift()};
    my %suites       = %{shift()};
    my %activeSteps  = %{shift()};
    my $queueManager =   shift() ;
    my %options      = %{shift()};
    my @jobsIdentify;
    foreach my $entry ( @entries ) {
	# Skip this entry if this suite does not require any analysis for which this step is required.
	next
	    unless ( grep {exists($activeAnalyses{$_})} @{$suites{$entry->{'suite'}->{'name'}}->{'analyses'}} );
	# Build redshift labels.
	@{$entry->{'resolution'}->{'epochs'}} = map {my $redshift = 1.0/$_-1.0; {expansionFactor => $_, redshift => $redshift, redshiftLabel => sprintf("z%5.3f",$redshift)}} @{$entry->{'resolution'}->{'expansionFactors'}};
	# Determine minimum and maximum masses for the mass function.
	$entry->{'massMinimum'} = 10.0**(int(log($entry->{'resolution'}->{'massParticle'})/log(10.0))+2);
	$entry->{'massMaximum'} = 1.0e16;
	# Construct the simulation path.
	$entry->{'path'} = $options{'simulationDataPath'}.$entry->{'suite'}->{'name'}."/".$entry->{'group'}->{'name'}."/".$entry->{'resolution'}->{'name'}."/".$entry->{'simulation'}->{'name'}."/".$entry->{'realization'}."/";
	# Iterate over subvolumes.
	for(my $i=0;$i<$entry->{'resolution'}->{'subvolumes'};++$i) {
	    for(my $j=0;$j<$entry->{'resolution'}->{'subvolumes'};++$j) {
		for(my $k=0;$k<$entry->{'resolution'}->{'subvolumes'};++$k) {
		    # Skip any missing tree - not all simulations have the full complement of subvolumes.
		    next
			unless ( -e $entry->{'path'}."tree_".$i."_".$j."_".$k.".dat" );
		    # Parse the base parameters.
		    my $parameters = $xml->XMLin($options{'pipelinePath'}."identifyAlwaysIsolated.xml");
		    # Modify file names.
		    $parameters->{'outputFileName'}                                      ->{'value'} = $entry->{'path'}."identifyAlwaysIsolatedGLC_".$i."_".$j."_".$k.".hdf5";
		    $parameters->{'nbodyImporter' }                        ->{'fileName'}->{'value'} = $entry->{'path'}."tree_"                     .$i."_".$j."_".$k.".dat" ;
		    $parameters->{'nbodyOperator' }->{'nbodyOperator'}->[4]->{'fileName'}->{'value'} = $entry->{'path'}."alwaysIsolated_subVolume"  .$i."_".$j."_".$k.".hdf5";
		    # Modify cosmological parameters.
		    $parameters->{'cosmologyParameters'}->{$_}->{'value'} = $entry->{'suite'}->{'cosmology'}->{$_}
		        foreach ( 'HubbleConstant', 'OmegaMatter', 'OmegaDarkEnergy', 'OmegaBaryon' );
		    # If a custom parameter processing function is defined, call it.
		    &{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'identifyAlwaysIsolated'}->{'processParameters'}}($entry,undef(),$parameters,\%options)
			if ( exists($suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'identifyAlwaysIsolated'}->{'processParameters'}) );
		    # Write parmeter file.
		    my $parameterFileName = $entry->{'path'}."identifyAlwaysIsolated_".$i."_".$j."_".$k.".xml";
		    open(my $outputFile,">",$parameterFileName);
		    print $outputFile $xml->XMLout($parameters, RootName => "parameters");
		    close($outputFile);
		    # Skip if the file exists.
		    next
			if ( -e $entry->{'path'}."alwaysIsolated_subVolume".$i."_".$j."_".$k.".hdf5" );
		    # Generate a job.
		    my $job;
		    $job->{'command'   } =
			$ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$parameterFileName;
		    $job->{'launchFile'} = $entry->{'path'}."identifyAlwaysIsolated_".$i."_".$j."_".$k.".sh" ;
		    $job->{'logFile'   } = $entry->{'path'}."identifyAlwaysIsolated_".$i."_".$j."_".$k.".log";
		    $job->{'label'     } =                  "identifyAlwaysIsolated_".$i."_".$j."_".$k       ;
		    $job->{'ppn'       } = $ompThreads;
		    $job->{'ompThreads'} = $ompThreads;
		    $job->{'nodes'     } = 1 ;
		    $job->{'mem'       } = $entry->{'resolution'}->{'name'} eq "resolutionX64" ? "64G" : "8G";
		    $job->{'walltime'  } = "8:00:00";
		    $job->{'mpi'       } = "no";
		    push(@jobsIdentify,$job);
		}
	    }
	}
    }
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobsIdentify)
	if ( scalar(@jobsIdentify) > 0 );
}

sub stepExtractHalos {
    # Extract halo snapshots from simulations.
    my @entries      = @{shift()};
    my %suites       = %{shift()};
    my %activeSteps  = %{shift()};
    my $queueManager =   shift() ;
    my %options      = %{shift()};
    my @jobsExtract;
    foreach my $entry ( @entries ) {
	# Iterate over subvolumes.
	for(my $i=0;$i<$entry->{'resolution'}->{'subvolumes'};++$i) {
	    for(my $j=0;$j<$entry->{'resolution'}->{'subvolumes'};++$j) {
		for(my $k=0;$k<$entry->{'resolution'}->{'subvolumes'};++$k) {
		    # Skip any missing tree - not all simulations have the full complement of subvolumes.
		    next
			unless ( -e $entry->{'path'}."tree_".$i."_".$j."_".$k.".dat" );
		    # Iterate over expansion factors.
		    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
			my $expansionFactorLow  = (1.0-5.0e-4)*$epoch->{'expansionFactor'};
			my $expansionFactorHigh = (1.0+5.0e-4)*$epoch->{'expansionFactor'};
			# Parse the base parameters.
			my $parameters = $xml->XMLin($options{'pipelinePath'}."extractHalosSnapshot.xml");
			# Modify file names.
			$parameters->{'outputFileName'}                                           ->{'value'} = $entry->{'path'}."alwaysIsolated_subVolumeGLC"                     .$i."_".$j."_".$k.".hdf5";
			$parameters->{'nbodyImporter' }                        ->{'fileName'     }->{'value'} = $entry->{'path'}."alwaysIsolated_subVolume"                        .$i."_".$j."_".$k.".hdf5";
			$parameters->{'nbodyImporter' }                        ->{'properties'   }->{'value'} = "particleID isFlyby expansionFactor massVirial hostedRootID";
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]->{'propertyNames'}->{'value'} = "isFlyby expansionFactor";
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]->{'rangeLow'     }->{'value'} = "0 ".$expansionFactorLow ;
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]->{'rangeHigh'    }->{'value'} = "0 ".$expansionFactorHigh;
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[1]->{'propertyNames'}->{'value'} = "isFlyby expansionFactor";
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[2]->{'fileName'     }->{'value'} = $entry->{'path'}."nonFlyby_".$epoch->{'redshiftLabel'}."_subVolume".$i."_".$j."_".$k.".hdf5";
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[2]->{'redshift'     }->{'value'} =                              $epoch->{'redshift'     }                                      ;
			# If a custom parameter processing function is defined, call it.
			&{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'extractHalos'}->{'processParameters'}}($entry,$epoch->{'expansionFactor'},$parameters,\%options)
			    if ( exists($suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'extractHalos'}->{'processParameters'}) );
			# Write parameter file.
			my $parameterFileName = $entry->{'path'}."identifyNonFlyby_".$epoch->{'redshiftLabel'}."_".$i."_".$j."_".$k.".xml";
			open(my $outputFile,">",$parameterFileName);
			print $outputFile $xml->XMLout($parameters, RootName => "parameters");
			close($outputFile);
			# Skip if the file exists.
			next
			    if ( -e $entry->{'path'}."nonFlyby_".$epoch->{'redshiftLabel'}."_subVolume".$i."_".$j."_".$k.".hdf5" );
			# Generate a job.
			my $job;
			$job->{'command'   } =
			    $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$parameterFileName;
			$job->{'launchFile'} = $entry->{'path'}."identifyNonFlyby_".$epoch->{'redshiftLabel'}."_".$i."_".$j."_".$k.".sh" ;
			$job->{'logFile'   } = $entry->{'path'}."identifyNonFlyby_".$epoch->{'redshiftLabel'}."_".$i."_".$j."_".$k.".log";
			$job->{'label'     } =                  "identifyNonFlyby_".$epoch->{'redshiftLabel'}."_".$i."_".$j."_".$k       ;
			$job->{'ppn'       } = 1;
			$job->{'ompThreads'} = 1;
			$job->{'nodes'     } = 1;
			my $mem = "8G";
			$mem = "64G"
			    if ( $entry->{'resolution'}->{'name'} eq "resolutionX64" );
			$job->{'mem'       } = $mem;
			$job->{'walltime'  } = "8:00:00";
			$job->{'mpi'       } = "no";
			push(@jobsExtract,$job);
		    }
		}
	    }
	}
    }
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobsExtract)
	if ( scalar(@jobsExtract) > 0 );
}

sub stepExtractSubhalos {
    # Extract subhalo snapshots from simulations.
    my @entries      = @{shift()};
    my %suites       = %{shift()};
    my %activeSteps  = %{shift()};
    my $queueManager =   shift() ;
    my %options      = %{shift()};
    my @jobsExtract;
    foreach my $entry ( @entries ) {
	# Iterate over subvolumes.
	for(my $i=0;$i<$entry->{'resolution'}->{'subvolumes'};++$i) {
	    for(my $j=0;$j<$entry->{'resolution'}->{'subvolumes'};++$j) {
		for(my $k=0;$k<$entry->{'resolution'}->{'subvolumes'};++$k) {
		    # Skip any missing tree - not all simulations have the full complement of subvolumes.
		    next
			unless ( -e $entry->{'path'}."tree_".$i."_".$j."_".$k.".dat" );
		    # Iterate over expansion factors.
		    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
			my $expansionFactorLow  = (1.0-5.0e-4)*$epoch->{'expansionFactor'};
			my $expansionFactorHigh = (1.0+5.0e-4)*$epoch->{'expansionFactor'};
			# Parse the base parameters.
			my $parameters = $xml->XMLin($options{'pipelinePath'}."extractSubhalosSnapshot.xml");
			# Modify file names.
			$parameters->{'outputFileName'}                                           ->{'value'} = $entry->{'path'}."subhalos_subVolumeGLC"                     .$i."_".$j."_".$k.".hdf5";
			$parameters->{'nbodyImporter' }                        ->{'fileName'     }->{'value'} = $entry->{'path'}."alwaysIsolated_subVolume"                  .$i."_".$j."_".$k.".hdf5";
			$parameters->{'nbodyImporter' }                        ->{'properties'   }->{'value'} = "particleID isFlyby expansionFactor massVirial velocityMaximum";
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]->{'propertyNames'}->{'value'} = "isFlyby expansionFactor";
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]->{'rangeLow'     }->{'value'} = "1 ".$expansionFactorLow ;
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]->{'rangeHigh'    }->{'value'} = "1 ".$expansionFactorHigh;
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[1]->{'propertyNames'}->{'value'} = "isFlyby expansionFactor";
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[2]->{'fileName'     }->{'value'} = $entry->{'path'}."subhalos_".$epoch->{'redshiftLabel'}."_subVolume".$i."_".$j."_".$k.".hdf5";
			$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[2]->{'redshift'     }->{'value'} =                              $epoch->{'redshift'     }                                      ;
			# If a custom parameter processing function is defined, call it.
			&{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'extractSubhalos'}->{'processParameters'}}($entry,$epoch->{'expansionFactor'},$parameters,\%options)
			    if ( exists($suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'extractSubhalos'}->{'processParameters'}) );
			# Write parameter file.
			my $parameterFileName = $entry->{'path'}."identifySubhalos_".$epoch->{'redshiftLabel'}."_".$i."_".$j."_".$k.".xml";
			open(my $outputFile,">",$parameterFileName);
			print $outputFile $xml->XMLout($parameters, RootName => "parameters");
			close($outputFile);
			# Skip if the file exists.
			next
			    if ( -e $entry->{'path'}."subhalos_".$epoch->{'redshiftLabel'}."_subVolume".$i."_".$j."_".$k.".hdf5" );
			# Generate a job.
			my $job;
			$job->{'command'   } =
			    $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$parameterFileName;
			$job->{'launchFile'} = $entry->{'path'}."identifySubhalos_".$epoch->{'redshiftLabel'}."_".$i."_".$j."_".$k.".sh" ;
			$job->{'logFile'   } = $entry->{'path'}."identifySubhalos_".$epoch->{'redshiftLabel'}."_".$i."_".$j."_".$k.".log";
			$job->{'label'     } =                  "identifySubhalos_".$epoch->{'redshiftLabel'}."_".$i."_".$j."_".$k       ;
			$job->{'ppn'       } = 1;
			$job->{'ompThreads'} = 1;
			$job->{'nodes'     } = 1;
			my $mem = "8G";
			$mem = "64G"
			    if ( $entry->{'resolution'}->{'name'} eq "resolutionX64" );
			$job->{'mem'       } = $mem;
			$job->{'walltime'  } = "8:00:00";
			$job->{'mpi'       } = "no";
			push(@jobsExtract,$job);
		    }
		}
	    }
	}
    }
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobsExtract)
	if ( scalar(@jobsExtract) > 0 );
}

sub stepMassFunctions {
    # Construct halo mass functions from simulations.
    my @entries      = @{shift()};
    my %suites       = %{shift()};
    my %activeSteps  = %{shift()};
    my $queueManager =   shift() ;
    my %options      = %{shift()};
    my @massFunctionJobs;
    foreach my $entry ( @entries ) {
	# Iterate over expansion factors.
	foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	    # Iterate over subvolumes.
	    my @nbodyImporters;
	    for(my $i=0;$i<$entry->{'resolution'}->{'subvolumes'};++$i) {
		for(my $j=0;$j<$entry->{'resolution'}->{'subvolumes'};++$j) {
		    for(my $k=0;$k<$entry->{'resolution'}->{'subvolumes'};++$k) {
			# Skip any missing tree - not all simulations have the full complement of subvolumes.
			next
			    unless ( -e $entry->{'path'}."tree_".$i."_".$j."_".$k.".dat" );
			# Add an importer for this subvolume.
			push(
			    @nbodyImporters,
			    {
				value      => "IRATE"                                                                                                ,
				fileName   => {value => $entry->{'path'}."nonFlyby_".$epoch->{'redshiftLabel'}."_subVolume".$i."_".$j."_".$k.".hdf5"},
				snapshot   => {value => "1"},
				properties => {value => "massVirial"}
			    }
			    );
		    }
		}
	    }
	    # Compute the mass function.
	    unless ( -e $entry->{'path'}."haloMassFunction_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5" ) {
		# Parse the base parameters.
		my $massFunctionParameters = $xml->XMLin($options{'pipelinePath'}."haloMassFunctionCompute.xml");
		# Modify parameters.
		@{$massFunctionParameters->{'nbodyImporter' }->{'nbodyImporter'}}                                          = @nbodyImporters;
		$massFunctionParameters  ->{'outputFileName'}                                                  ->{'value'} = $entry->{'path'       }."haloMassFunction_".$epoch->{'redshiftLabel'}.".hdf5";
		$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[0]->{'values'             }->{'value'} = $entry->{'resolution' }              ->{'massParticle'};
		$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'simulationReference'}->{'value'} = $entry->{'group'      }->{'metaData'}->{'reference'   };
		$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'simulationURL'      }->{'value'} = $entry->{'group'      }->{'metaData'}->{'url'         };
		$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'massMinimum'        }->{'value'} = $entry->{'massMinimum'}                                ;
		$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'massMaximum'        }->{'value'} = $entry->{'massMaximum'}                                ;
		$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'description'        }->{'value'} = "Halo mass function of non-flyby halos for the ".$entry->{'suite'}->{'name'}." ".$entry->{'group'}->{'name'}." ".$entry->{'simulation'}->{'name'}." ".$epoch->{'redshiftLabel'}." simulation";
		# If a custom parameter processing function is defined, call it.
		&{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'massFunctions'}->{'processParameters'}}($entry,$epoch->{'expansionFactor'},$massFunctionParameters,\%options)
		    if ( exists($suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'massFunctions'}->{'processParameters'}) );
		# Write the parameter file.
		my $parameterFileName = $entry->{'path'}."haloMassFunction_".$epoch->{'redshiftLabel'}.".xml";
		open(my $outputFile,">",$parameterFileName);
		print $outputFile $xml->XMLout($massFunctionParameters, RootName => "parameters");
		close($outputFile);
		# Construct the job.
		my $job;
		$job->{'command'   } =
		    $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$parameterFileName;
		$job->{'launchFile'} = $entry->{'path'}."haloMassFunction_".$epoch->{'redshiftLabel'}.".sh" ;
		$job->{'logFile'   } = $entry->{'path'}."haloMassFunction_".$epoch->{'redshiftLabel'}.".log";
		$job->{'label'     } =                  "haloMassFunction_".$epoch->{'redshiftLabel'}       ;
		$job->{'ppn'       } = $ompThreads;
		$job->{'ompThreads'} = $ompThreads;
		$job->{'nodes'     } = 1;
		$job->{'mem'       } = "8G";
		$job->{'walltime'  } = "8:00:00";
		$job->{'mpi'       } = "no";
		push(@massFunctionJobs,$job);
	    }
	}
    }
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@massFunctionJobs)
	if ( scalar(@massFunctionJobs) > 0 );
}


sub stepSubhaloFunctions {
    # Construct subhalo mass, radial, and Vmax functions from simulations.
    my @entries      = @{shift()};
    my %suites       = %{shift()};
    my %activeSteps  = %{shift()};
    my $queueManager =   shift() ;
    my %options      = %{shift()};
    my @subhaloFunctionJobs;
    foreach my $entry ( @entries ) {
	# Iterate over expansion factors.
	foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	    # Iterate over subvolumes.
	    my @nbodyImporters;
	    for(my $i=0;$i<$entry->{'resolution'}->{'subvolumes'};++$i) {
		for(my $j=0;$j<$entry->{'resolution'}->{'subvolumes'};++$j) {
		    for(my $k=0;$k<$entry->{'resolution'}->{'subvolumes'};++$k) {
			# Skip any missing tree - not all simulations have the full complement of subvolumes.
			next
			    unless ( -e $entry->{'path'}."tree_".$i."_".$j."_".$k.".dat" );
			# Add an importer for this subvolume.
			push(
			    @nbodyImporters,
			    {
				value      => "IRATE"                                                                                                ,
				fileName   => {value => $entry->{'path'}."subhalos_".$epoch->{'redshiftLabel'}."_subVolume".$i."_".$j."_".$k.".hdf5"},
				snapshot   => {value => "1"},
				properties => {value => "position massVirial velocityMaximum"}
			    }
			    );
		    }
		}
	    }
	    # Compute the subhalo functions.
	    unless ( -e $entry->{'path'}."subhaloFunctions_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5" ) {
		# Parse the base parameters.
		my $subhaloFunctionParameters = $xml->XMLin($options{'pipelinePath'}."subhaloFunctionsCompute.xml");
		# Modify parameters.
		@{$subhaloFunctionParameters->{'nbodyImporter' }->{'nbodyImporter'}}                                          = @nbodyImporters;
		$subhaloFunctionParameters  ->{'outputFileName'}                                                  ->{'value'} = $entry->{'path'       }."subhaloFunctions_".$epoch->{'redshiftLabel'}.".hdf5";
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[0]->{'values'             }->{'value'} = $entry->{'resolution' }              ->{'massParticle'};
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[2]->{'simulationReference'}->{'value'} = $entry->{'group'      }->{'metaData'}->{'reference'   };
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[2]->{'simulationURL'      }->{'value'} = $entry->{'group'      }->{'metaData'}->{'url'         };
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[2]->{'description'        }->{'value'} = "Subhalo mass function of non-flyby halos for the ".$entry->{'suite'}->{'name'}." ".$entry->{'group'}->{'name'}." ".$entry->{'simulation'}->{'name'}." ".$epoch->{'redshiftLabel'}." simulation";
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[3]->{'simulationReference'}->{'value'} = $entry->{'group'      }->{'metaData'}->{'reference'   };
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[3]->{'simulationURL'      }->{'value'} = $entry->{'group'      }->{'metaData'}->{'url'         };
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[3]->{'description'        }->{'value'} = "Subhalo radial distribution function of non-flyby halos for the ".$entry->{'suite'}->{'name'}." ".$entry->{'group'}->{'name'}." ".$entry->{'simulation'}->{'name'}." ".$epoch->{'redshiftLabel'}." simulation";
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[4]->{'simulationReference'}->{'value'} = $entry->{'group'      }->{'metaData'}->{'reference'   };
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[4]->{'simulationURL'      }->{'value'} = $entry->{'group'      }->{'metaData'}->{'url'         };
		$subhaloFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[4]->{'description'        }->{'value'} = "Subhalo \$V_\\mathrm{max}\$ function of non-flyby halos for the ".$entry->{'suite'}->{'name'}." ".$entry->{'group'}->{'name'}." ".$entry->{'simulation'}->{'name'}." ".$epoch->{'redshiftLabel'}." simulation";
		# If a custom parameter processing function is defined, call it.
		&{$suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'subhaloFunctions'}->{'processParameters'}}($entry,$epoch->{'expansionFactor'},$subhaloFunctionParameters,\%options)
		    if ( exists($suites{$entry->{'suite'}->{'name'}}->{'steps'}->{'subhaloFunctions'}->{'processParameters'}) );
		# Write the parameter file.
		my $parameterFileName = $entry->{'path'}."subhaloFunctions_".$epoch->{'redshiftLabel'}.".xml";
		open(my $outputFile,">",$parameterFileName);
		print $outputFile $xml->XMLout($subhaloFunctionParameters, RootName => "parameters");
		close($outputFile);
		# Construct the job.
		my $job;
		$job->{'command'   } =
		    $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$parameterFileName;
		$job->{'launchFile'} = $entry->{'path'}."subhaloFunctions_".$epoch->{'redshiftLabel'}.".sh" ;
		$job->{'logFile'   } = $entry->{'path'}."subhaloFunctions_".$epoch->{'redshiftLabel'}.".log";
		$job->{'label'     } =                  "subhaloFunctions_".$epoch->{'redshiftLabel'}       ;
		$job->{'ppn'       } = $ompThreads;
		$job->{'ompThreads'} = $ompThreads;
		$job->{'nodes'     } = 1;
		$job->{'mem'       } = "8G";
		$job->{'walltime'  } = "8:00:00";
		$job->{'mpi'       } = "no";
		push(@subhaloFunctionJobs,$job);
	    }
	}
    }
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@subhaloFunctionJobs)
	if ( scalar(@subhaloFunctionJobs) > 0 );
}

## Symphony(/COZMIC)-specific processors.

sub symphonyProcessIdentifyAlwaysIsolated {
    # Set the appropriate cosmology.
    my $entry           =   shift() ;
    my $expansionFactor =   shift() ;
    my $parameters      =   shift() ;
    my %options         = %{shift()};
    # Find the host halo ID for this simulation.
    die("can not find host halo ID for ".$entry->{'suite'}->{'name'}."; ".$entry->{'group'}->{'name'}."; ".$entry->{'resolution'}->{'name'}."; ".$entry->{'simulation'}->{'name'}."; ".$entry->{'realization'})
	unless ( exists($entry->{'simulation'}->{'hostHaloIDs'}->{$entry->{'realization'}}) );
    $entry->{'hostHaloID'} = $entry->{'simulation'}->{'hostHaloIDs'}->{$entry->{'realization'}};
    # Add read of additional columns.
    my @propertiesImport = split(" ",$parameters->{'nbodyImporter'}->{'readColumns'}->{'value'});
    foreach my $property ( "X", "Y", "Z", "Rvir", "rs", "Vmax" ) {
	$parameters->{'nbodyImporter'}->{'readColumns'}->{'value'} .= " ".$property
	    unless ( grep {$_ eq $property} @propertiesImport );
    }
    # Remove any delete of these properties.
    my @propertiesDeleted = split(" ",$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'propertyNames'}->{'value'});
    my @propertiesToDelete;
    foreach my $property ( @propertiesDeleted ) {
	push(@propertiesToDelete,$property)
	    unless
	    (
	     $property eq "position"
	     ||
	     $property eq "radiusScale"
	     ||
	     $property eq "radiusVirial"
	     ||
	     $property eq "velocityMaximum"
	    );
    }
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'propertyNames'}->{'value'} = join(" ",@propertiesToDelete);
    # Add physical to comoving conversion.
    splice(
	@{$parameters->{'nbodyOperator'}->{'nbodyOperator'}},
	2,0,
	{
	    value => "physicalToComoving"
	}
	);
}

sub symphonyPreProcessExtractHalosLocate {
    # Identify the primary progenitor.
    my $entry   =   shift() ;
    my $jobs    =   shift() ;
    my %options = %{shift()};
    # Find the host halo ID for this realization.
    die("can not find host halo ID for ".$entry->{'suite'}->{'name'}."; ".$entry->{'group'}->{'name'}."; ".$entry->{'resolution'}->{'name'}."; ".$entry->{'simulation'}->{'name'}."; ".$entry->{'realization'})
	unless ( exists($entry->{'simulation'}->{'hostHaloIDs'}->{$entry->{'realization'}}) );
    my $hostHaloID = $entry->{'simulation'}->{'hostHaloIDs'}->{$entry->{'realization'}};
    # Iterate over expansion factors.
    my $job;
    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	# Look for an existing record of the primary halo.
	my $primaryHaloFileName    = $entry->{'path'}."primaryHalo_".$epoch->{'redshiftLabel'}.".xml";
	my $xml                    = new XML::Simple();
	unless ( -e $primaryHaloFileName ) {
	    # Generate a job to extract the primary halo.
	    unless ( defined($job) ) {
		$job->{'command'   } = "";
		$job->{'launchFile'} = $entry->{'path'}."zoomInExtract_".$entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$entry->{'realization'}.".sh" ;
		$job->{'logFile'   } = $entry->{'path'}."zoomInExtract_".$entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$entry->{'realization'}.".log";
		$job->{'label'     } =                  "zoomInExtract_".$entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$entry->{'realization'}       ;
		$job->{'ppn'       } = 1;
		$job->{'ompThreads'} = 1;
		$job->{'nodes'     } = 1;
		my $mem = "16G";
		$mem = "32G"
		    if ( $entry->{'group'}->{'name'} eq "Group" );
		$mem = "64G"
		    if ( $entry->{'resolution'}->{'name'} eq "resolutionX64" );
		$job->{'mem'       } = $mem;
		$job->{'walltime'  } = "8:00:00";
		$job->{'mpi'       } = "no";
	    }
	    $job->{'command'} .= $options{'pipelinePath'}."zoomInExtract.py ".$entry->{'path'}." ".$primaryHaloFileName." ".$epoch->{'expansionFactor'}." ".$entry->{'suite'}->{'cosmology'}->{'hubbleConstant'}." ".$entry->{'resolution'}->{'massParticle'}." ".$entry->{'resolution'}->{'haloMassFunction'}->{'massHostLogMin'}->{'value'}." ".$entry->{'resolution'}->{'haloMassFunction'}->{'massHostLogMax'}->{'value'}." ".$hostHaloID."\n";
	}
    }   
    push(@{$jobs},$job)
	if ( defined($job) );
}

sub symphonyPreProcessExtractHalosUncontaminated {
    # Find the uncontaminated region.
    my $entry   =   shift() ;
    my $jobs    =   shift() ;
    my %options = %{shift()};
    # Iterate over expansion factors.
    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	# Look for an existing record of the primary halo.
	my $primaryHaloFileName  = $entry->{'path'}."primaryHalo_".$epoch->{'redshiftLabel'}.".xml";
	# Find properties of the central halo.
	my $xml                  = new XML::Simple();
	my $primaryHaloData      = $xml->XMLin($primaryHaloFileName);
	my $expansionFactorLabel = sprintf("sphericalOrigin:a%5.3f",$epoch->{'expansionFactor'});
	$entry->{$expansionFactorLabel} =
	    $primaryHaloData->{'x'}." ".
	    $primaryHaloData->{'y'}." ".
	    $primaryHaloData->{'z'}    ;
	$entry->{'massCentral'} = pdl $primaryHaloData->{'mc'};
	# Generate job to find the uncontaminated region.
	my $uncontaminatedFileName = $entry->{'path'}."uncontaminated_".$epoch->{'redshiftLabel'}.".hdf5";
	unless ( -e $uncontaminatedFileName ) {
	    my $parametersUncontaminated = $xml->XMLin($options{'pipelinePath'}."zoomInSelectUncontaminated.xml");
	    # Modify file names.
	    (my $snapshot) = map {$entry->{'resolution'}->{'expansionFactors'}->[$_] == $epoch->{'expansionFactor'} ? $entry->{'resolution'}->{'snapshots'}->[$_] : ()} 0..$#{$entry->{'resolution'}->{'expansionFactors'}};
	    $parametersUncontaminated->{'nbodyImporter'}                        ->{'fileName'      }->{'value'} = $entry->{'path'}."snapshots/snapshot_".$snapshot;
	    $parametersUncontaminated->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'point'         }->{'value'} = $entry->{$expansionFactorLabel};
	    $parametersUncontaminated->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'fileName'      }->{'value'} = $uncontaminatedFileName;
	    $parametersUncontaminated                                           ->{'outputFileName'}->{'value'} = $entry->{'path'}."uncontaminatedExtract_".$epoch->{'redshiftLabel'}.".hdf5";
	    # Write parameter file.
	    my $parameterFileName = $entry->{'path'}."uncontaminated_".$epoch->{'redshiftLabel'}.".xml";
	    open(my $outputFile,">",$parameterFileName);
	    print $outputFile $xml->XMLout($parametersUncontaminated, RootName => "parameters");
	    close($outputFile);
	    # Generate a job to extract the uncontaminated region around the primary halo.
	    my $job;
	    $job->{'command'   } = $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$entry->{'path'}."uncontaminated_".$epoch->{'redshiftLabel'}.".xml";
	    $job->{'launchFile'} = $entry->{'path'}."uncontaminatedExtract_".$entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$entry->{'realization'}."_".$epoch->{'redshiftLabel'}.".sh" ;
	    $job->{'logFile'   } = $entry->{'path'}."uncontaminatedExtract_".$entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$entry->{'realization'}."_".$epoch->{'redshiftLabel'}.".log";
	    $job->{'label'     } =                  "uncontaminatedExtract_".$entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$entry->{'realization'}."_".$epoch->{'redshiftLabel'}       ;
	    $job->{'ppn'       } = $ompThreads;
	    $job->{'ompThreads'} = $ompThreads;
	    $job->{'nodes'     } = 1;
	    $job->{'mem'       } = $entry->{'resolution'}->{'name'} eq "resolutionX64" ? "196G" : "32G";
	    $job->{'walltime'  } = "8:00:00";
	    $job->{'mpi'       } = "no";
	    push(@{$jobs},$job)
	}
    }
}

sub symphonyProcessExtractHalos {
    # Determine the central point, and extent of the high-resolution region in a zoom in simulation.
    my $entry           =   shift() ;
    my $expansionFactor =   shift() ;
    my $parameters      =   shift() ;
    my %options         = %{shift()};
    print "Processing halo extraction for ".$entry->{'suite'}->{'name'}." : ".$entry->{'group'}->{'name'}." : ".$entry->{'simulation'}->{'name'}." : ".$entry->{'realization'}."\n";
    # Set path names.
    my $redshift               =  1.0/$expansionFactor-1.0;
    my $epoch->{'redshiftLabel'}          = sprintf("z%5.3f",$redshift);
    my $expansionFactorLabel   = sprintf("sphericalOrigin:a%5.3f",$expansionFactor);
    my $primaryHaloFileName    = $entry->{'path'}."primaryHalo_"   .$epoch->{'redshiftLabel'}.".xml";
    my $uncontaminatedFileName = $entry->{'path'}."uncontaminated_".$epoch->{'redshiftLabel'}.".hdf5";
    # Find properties of the central halo.
    my $xml                    = new XML::Simple();
    my $primaryHaloData        = $xml->XMLin($primaryHaloFileName);
    # Find uncontaminated region.
    my $uncontaminatedFile      = new PDL::IO::HDF5($uncontaminatedFileName);
    my $uncontaminatedParticles = $uncontaminatedFile     ->group  ('Snapshot00001/HaloCatalog');
    (my $radiusUncontaminated)  = $uncontaminatedParticles->attrGet('radiusUncontaminated'     );
    $primaryHaloData->{'ur'} = $radiusUncontaminated->sclr();
    open(my $primaryHaloDataFile,">",$primaryHaloFileName);
    print $primaryHaloDataFile $xml->XMLout($primaryHaloData, RootName => "primaryHalo");
    close($primaryHaloDataFile);
    # Add read of virial and scale radii columns.
    $parameters->{'nbodyImporter'}                        ->{'properties'   }->{'value'} .= " radiusScale radiusVirial";
    # Add read of (x,y,z) coordinate columns, and subsequent delete.
    $parameters->{'nbodyImporter'}                        ->{'properties'   }->{'value'} .= " position"         ;
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'propertyNames'}->{'value'} .= " distanceFromPoint";
    # Add calculation of distance from primary halo.
    $entry->{$epoch->{'redshiftLabel'}}->{'sphericalRadiusMinimum'} = 0.0;
    $entry->{$epoch->{'redshiftLabel'}}->{'sphericalRadiusMaximum'} = $radiusUncontaminated->sclr();
    splice(
	@{$parameters->{'nbodyOperator'}->{'nbodyOperator'}},
	1,0,
	{
	    value         => "distanceFromPoint",
	    point         => {value => $entry->{$expansionFactorLabel}}
	},
	{
	    value         =>           "filterProperties"                                              ,
	    propertyNames => {value => "distanceFromPoint"                                            },
	    rangeLow      => {value => $entry->{$epoch->{'redshiftLabel'}}->{'sphericalRadiusMinimum'}},
	    rangeHigh     => {value => $entry->{$epoch->{'redshiftLabel'}}->{'sphericalRadiusMaximum'}},
	},
	{
	    value         => "setBoxSize",
	    boxSize       => {value => $primaryHaloData->{'l'}}
	}
	);
}

sub symphonyProcessExtractSubhalos {
    # Determine the central point, and extent of the host virial radius in a zoom in simulation.
    my $entry           =   shift() ;
    my $expansionFactor =   shift() ;
    my $parameters      =   shift() ;
    my %options         = %{shift()};
    print "Processing subhalo extraction for ".$entry->{'suite'}->{'name'}." : ".$entry->{'group'}->{'name'}." : ".$entry->{'simulation'}->{'name'}." : ".$entry->{'realization'}."\n";
    # Set path names.
    my $redshift                 = 1.0/$expansionFactor-1.0;
    my $epoch->{'redshiftLabel'} = sprintf("z%5.3f",$redshift);
    my $expansionFactorLabel     = sprintf("sphericalOrigin:a%5.3f",$expansionFactor);
    my $primaryHaloFileName      = $entry->{'path'}."primaryHalo_".$epoch->{'redshiftLabel'}.".xml";
    # Find properties of the central halo.
    my $xml                    = new XML::Simple();
    my $primaryHaloData        = $xml->XMLin($primaryHaloFileName);
    # Add read of (x,y,z) coordinate columns, and subsequent delete.
    $parameters->{'nbodyImporter'}                        ->{'properties'   }->{'value'} .= " position"         ;
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'propertyNames'}->{'value'} .= " distanceFromPoint";
    # Add calculation of distance from primary halo.
    splice(
	@{$parameters->{'nbodyOperator'}->{'nbodyOperator'}},
	1,0,
	{
	    value         => "distanceFromPoint",
	    point         => {value => $entry->{$expansionFactorLabel}}
	},
	{
	    value         =>           "filterProperties"      ,
	    propertyNames => {value => "distanceFromPoint"    },
	    rangeLow      => {value => 0.0                    },
	    rangeHigh     => {value => $primaryHaloData->{'r'}},
	}
	);
}

sub symphonyProcessSubhaloFunctions {
    # Determine the host halo mass and mass ratio range.
    my $entry           =   shift() ;
    my $expansionFactor =   shift() ;
    my $parameters      =   shift() ;
    my %options         = %{shift()};
    # Set path names.
    my $redshift                 = 1.0/$expansionFactor-1.0;
    my $redshiftLabel            = sprintf("z%5.3f",$redshift);
    my $primaryHaloFileName      = $entry->{'path'}."primaryHalo_".$redshiftLabel.".xml";

    # Find all primary halos.
    my $massPrimaryMaximum = 0.0;
    {
	my $xml          = new XML::Simple();
	my @realizations = split(" ",$entry->{'simulation'}->{'realizations'}->{'value'});
	foreach my $realization ( @realizations ) {
	    (my $path = $entry->{'path'}) =~ s/\/$entry->{'realization'}\//\/$realization\//;
	    my $primaryHaloFileName = $path."primaryHalo_".$redshiftLabel.".xml";
	    my $primaryHaloData        = $xml->XMLin($primaryHaloFileName);
	    $massPrimaryMaximum = $primaryHaloData->{'m'}
	        if ( $primaryHaloData->{'m'} > $massPrimaryMaximum );
	}
    }
    # Find properties of the central halo.
    my $xml                    = new XML::Simple();
    my $primaryHaloData        = $xml->XMLin($primaryHaloFileName);
    # Set particle count limits. Limits on analysis masses are chosen such thatg bins are sligned with integer log₁₀ masses, with
    # the lowest mass bin center chosen to be at least half of a bin width above the lowest allowed mass to ensure that only
    # sufficiently well-resolved subhalos are included.
    my $massCountParticlesMinimum      =  300;
    my $structureCountParticlesMinimum = 2000;
    # Set the host position.
    my $expansionFactorLabel = sprintf("sphericalOrigin:a%5.3f",$expansionFactor);
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'point'             }->{'value'} = $entry->{$expansionFactorLabel};
    # Subhalo mass function.
    ## Set the host halo mass.
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'massHost'          }->{'value'} = $primaryHaloData->{'m'};
    ## Determine a suitable mass ratio range.
    my $massFunctionCountPerDecade                                                            = 3;
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'massRatioMinimum'  }->{'value'} = sclr(10.0**((floor(log($massCountParticlesMinimum*$entry->{'resolution'}->{'massParticle'}/$massPrimaryMaximum)/log(10.0)*$massFunctionCountPerDecade)+1.5)/$massFunctionCountPerDecade));
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'massRatioMaximum'  }->{'value'} = 1.0;
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'massCountPerDecade'}->{'value'} = $massFunctionCountPerDecade;

    # Subhalo radial function.
    ## Set the host virial radius.
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[3]->{'radiusVirialHost'  }->{'value'} = $primaryHaloData->{'r'};
    ## Set the minimum subhalo mass.
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[3]->{'massMinimum'       }->{'value'} = $massCountParticlesMinimum*$entry->{'resolution'}->{'massParticle'};
    # Subhalo Vmax function.
    ## Determine a suitable mass range.
    my $velocityFunctionCountPerDecade                                                        = 5;
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'massMinimum'       }->{'value'} = sclr(10.0**((floor(log($structureCountParticlesMinimum*$entry->{'resolution'}->{'massParticle'})/log(10.0)*$velocityFunctionCountPerDecade)+1.5)/$velocityFunctionCountPerDecade));
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'massMaximum'       }->{'value'} = sclr(10.0**( floor(log($massPrimaryMaximum                                                     )/log(10.0)                                )+1.0                                 ));
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'massCountPerDecade'}->{'value'} = $velocityFunctionCountPerDecade;
}

sub symphonyPostprocessSelectInSphere {
    # Select all particles within the sphere of interest.
    my $entry   =   shift() ;
    my $jobs    =   shift() ;
    my %options = %{shift()};
    # Generate a parameter file to extract particles within the spherical volume of interest.
    my $xml        = new XML::Simple();
    my $parameters = $xml->XMLin($options{'pipelinePath'}."zoomInSelectInSphere.xml");
    # Iterate over expansion factors.
    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	my $expansionFactorLabel = sprintf("sphericalOrigin:a%5.3f",$epoch->{'expansionFactor'});
	(my $snapshot) = map {$entry->{'resolution'}->{'expansionFactors'}->[$_] == $epoch->{'expansionFactor'} ? $entry->{'resolution'}->{'snapshots'}->[$_] : ()} 0..$#{$entry->{'resolution'}->{'expansionFactors'}};
	$parameters->{'outputFileName'}                                       ->{'value'} = $entry->{'path'}."selectInSphereGLC_".$epoch->{'redshiftLabel'}.".hdf5";
	$parameters->{'nbodyImporter' }                        ->{'fileName' }->{'value'} = $entry->{'path'}."snapshots/snapshot_".$snapshot;
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[4]->{'fileName' }->{'value'} = $entry->{'path'}."selectedParticles_".$epoch->{'redshiftLabel'}.".hdf5";
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[4]->{'redshift' }->{'value'} =                                       $epoch->{'redshift'     }        ;
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]->{'point'    }->{'value'} = $entry->{$expansionFactorLabel}                                        ;
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[1]->{'rangeLow' }->{'value'} = $entry->{$epoch->{'redshiftLabel'}}->{'sphericalRadiusMinimum'} ;
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[1]->{'rangeHigh'}->{'value'} = $entry->{$epoch->{'redshiftLabel'}}->{'sphericalRadiusMaximum'} ;
	my $parameterFileName = $entry->{'path'}."selectInSphere_".$epoch->{'redshiftLabel'}.".xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($parameters, RootName => "parameters");
	close($outputFile);
	my $job;
	$job->{'command'   } =
	    $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $entry->{'path'}."selectInSphere_".$epoch->{'redshiftLabel'}.".sh" ;
	$job->{'logFile'   } = $entry->{'path'}."selectInSphere_".$epoch->{'redshiftLabel'}.".log";
	$job->{'label'     } =                  "selectInSphere_".$epoch->{'redshiftLabel'}       ;
	$job->{'ppn'       } = $ompThreads;
	$job->{'ompThreads'} = $ompThreads;
	$job->{'nodes'     } =  1;
	my $mem = "16G";
	$mem = "32G"
	    if ( $entry->{'group'     }->{'name'} eq "Group"         );
	$mem = "196G"
	    if ( $entry->{'resolution'}->{'name'} eq "resolutionX64" );
	$job->{'mem'       } = $mem;
	$job->{'walltime'  } = "8:00:00";
	$job->{'mpi'       } = "no";
	push(@{$jobs},$job)
	    unless ( -e $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'fileName' }->{'value'} );
    }
}

sub symphonyPostprocessSelectInICs {
    # Identify selected particles in the initial conditions.
    my $entry   =   shift() ;
    my $jobs    =   shift() ;
    my %options = %{shift()};
    # Find the redshift of the ICs.
    my $musicFileName = $entry->{'path'}."music.conf";
    die("missing file: ".$musicFileName)
	unless ( -e $musicFileName );
    open(my $musicFile,$musicFileName);
    while ( my $line = <$musicFile> ) {
	if ( $line =~ m/^zstart\s*=\s*([\d\.]+)/ ) {
	    $entry->{'redshiftICs'} = $1;
	}
    }
    close($musicFile);
    # Generate a parameter file to extract selected particles in the ICs.
    my $xml        = new XML::Simple();
    my $parameters = $xml->XMLin($options{'pipelinePath'}."zoomInSelectInICs.xml");
    # Iterate over expansion factors.
    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	$parameters->{'nbodyImporter'}                        ->{'fileName'           }->{'value'} =                $entry->{'path'       }."ic/ic_gadget_dist"                                        ;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'idSelectionFileName'}->{'value'} =                $entry->{'path'       }."selectedParticles_".$epoch->{'redshiftLabel'}.    ".hdf5" ;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'fileName'           }->{'value'} =                $entry->{'path'       }."selectedParticles_".$epoch->{'redshiftLabel'}."_ICs.hdf5" ;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'redshift'           }->{'value'} = sprintf("%.3f",$entry->{'redshiftICs'}                                                           );
	$parameters                                           ->{'outputFileName'     }->{'value'} =                $entry->{'path'       }."selectInICs_"      .$epoch->{'redshiftLabel'}.    ".hdf5" ;
	my $parameterFileName = $entry->{'path'}."selectInICs_".$epoch->{'redshiftLabel'}.".xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($parameters, RootName => "parameters");
	close($outputFile);
	my $memory = "16G";
	$memory = "32G"
	    if ( $entry->{'group'}->{'name'} eq "Group" );
	$memory = "32G"
	    if ( $entry->{'group'}->{'name'} eq "MilkyWay" && $entry->{'resolution'}->{'name'} eq "resolutionX8"  );
	$memory = "196G"
	    if ( $entry->{'group'}->{'name'} eq "MilkyWay" && $entry->{'resolution'}->{'name'} eq "resolutionX64" );
	my $job;
	$job->{'command'   } =
	    $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $entry->{'path'}."selectInICs_".$epoch->{'redshiftLabel'}.".sh" ;
	$job->{'logFile'   } = $entry->{'path'}."selectInICs_".$epoch->{'redshiftLabel'}.".log";
	$job->{'label'     } =                  "selectInICs_".$epoch->{'redshiftLabel'}       ;
	$job->{'ppn'       } = $ompThreads;
	$job->{'ompThreads'} = $ompThreads;
	$job->{'nodes'     } =  1;
	$job->{'mem'       } = $memory;
	$job->{'walltime'  } = "8:00:00";
	$job->{'mpi'       } = "no";
	push(@{$jobs},$job)
	    unless ( -e $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'fileName' }->{'value'} );
    }
}

sub symphonyPostprocessAnalyze {
    # Measure mass and overdensity of the selected region.
    my $entry   =   shift() ;
    my $jobs    =   shift() ;
    my %options = %{shift()};
    # Generate a parameter file to measure mass and overdensity.
    my $xml        = new XML::Simple();
    my $parameters = $xml->XMLin($options{'pipelinePath'}."zoomInAnalyze.xml");
    # Iterate over expansion factors.
    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	$parameters->{'outputFileName'}                                                         ->{'value'} =                $entry->{'path'       }."environment_"      .$epoch->{'redshiftLabel'}.".hdf5"     ;
	$parameters->{'nbodyImporter' }->{'nbodyImporter'}->[0]                   ->{'fileName'}->{'value'} =                $entry->{'path'       }."selectedParticles_".$epoch->{'redshiftLabel'}."_ICs.hdf5" ;
	$parameters->{'nbodyImporter' }->{'nbodyImporter'}->[1]                   ->{'fileName'}->{'value'} =                $entry->{'path'       }."ic/ic_gadget_dist"                             ;
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]                   ->{'values'  }->{'value'} = sprintf("%.3f",$entry->{'redshiftICs'}                                                );
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[5]->{'nbodyOperator'}->{'fileName'}->{'value'} =                $entry->{'path'       }."allParticles_"     .$epoch->{'redshiftLabel'}."_ICs.hdf5" ;
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[5]->{'nbodyOperator'}->{'redshift'}->{'value'} = sprintf("%.3f",$entry->{'redshiftICs'}                                                );
	my $parameterFileName = $entry->{'path'}."analyze_".$epoch->{'redshiftLabel'}.".xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($parameters, RootName => "parameters");
	close($outputFile);
	my $job;
	$job->{'command'   } =
	    $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $entry->{'path'}."analyze_".$epoch->{'redshiftLabel'}.".sh" ;
	$job->{'logFile'   } = $entry->{'path'}."analyze_".$epoch->{'redshiftLabel'}.".log";
	$job->{'label'     } =                  "analyze_".$epoch->{'redshiftLabel'}       ;
	$job->{'ppn'       } = $ompThreads;
	$job->{'ompThreads'} = $ompThreads;
	$job->{'nodes'     } =  1;
	my $mem = "48G";
	$mem = "64G"
	    if ( $entry->{'resolution'}->{'name'} eq "resolutionX8"  );
	$mem = "256G"
	    if ( $entry->{'resolution'}->{'name'} eq "resolutionX64" );
	$job->{'mem'       } = $mem;
	$job->{'walltime'  } = "8:00:00";
	$job->{'mpi'       } = "no";
	push(@{$jobs},$job)
	    unless ( -e $entry->{'path'}."environment_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5" );
    }
}

sub symphonyPostprocessSetVolume {
    # Set the effective Lagrangian volume for the selected region.
    my $entry   =   shift() ;
    my $jobs    =   shift() ;
    my %options = %{shift()};
    # Iterate over expansion factors.
    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	# Extract the mass of the region.
	my $analysisFile          = new PDL::IO::HDF5($entry->{'path'}."selectedParticles_".$epoch->{'redshiftLabel'}.".hdf5");
	die("error: 'massTotal' attribute is missing from file '".$entry->{'path'}."selectedParticles_".$epoch->{'redshiftLabel'}.".hdf5'")
	    unless ( grep {$_ eq "massTotal"} $analysisFile->group('Snapshot00001')->group('HaloCatalog')->attrs() );
	(my $mass)                = $analysisFile->group('Snapshot00001')->group('HaloCatalog')->attrGet('massTotal');
	# Find properties of the central halo.
	my $primaryHaloFileName   = $entry->{'path'}."primaryHalo_".$epoch->{'redshiftLabel'}.".xml";
	my $xml                   = new XML::Simple();
	my $primaryHaloData       = $xml->XMLin($primaryHaloFileName);
	my $massPrimary           = $primaryHaloData->{'m'};
	# Compute the mean density of the universe.
	my $gravitationalConstant = pdl 4.3011827419096073e-9;
	my $densityMean           = 3.0*$entry->{'suite'}->{'cosmology'}->{'OmegaMatter'}*$entry->{'suite'}->{'cosmology'}->{'HubbleConstant'}**2/8.0/PI/$gravitationalConstant;
	# Compute the volume of a cube containing the mass of our selected region.
	die("Primary halo mass exceeds region mass")
	    if ( $massPrimary >= $mass );
	my $boxSize                                = (    ($mass-$massPrimary)/$densityMean       )**(1.0/3.0);
	$entry->{$epoch->{'redshiftLabel'}}->{'radiusRegion'} = (3.0* $mass              /$densityMean/4.0/PI)**(1.0/3.0);
	$entry->{$epoch->{'redshiftLabel'}}->{'massRegion'  } =       $mass                                              ;
	$entry->{$epoch->{'redshiftLabel'}}->{'massPrimary' } =             $massPrimary                                 ;
	# Set the box size.
	my $halosFile = new PDL::IO::HDF5(">".$entry->{'path'}."nonFlyby_".$epoch->{'redshiftLabel'}."_subVolume0_0_0.hdf5");
	$halosFile->group('Snapshot00001'       )->group('HaloCatalog')->attrSet(boxSize => $boxSize);
	$halosFile->group('SimulationProperties')                      ->attrSet(boxSize => $boxSize);
    }
}

sub symphonyPostProcessMassFunction {
    # Store the region mass and overdensity in the mass function file.
    my $entry   =   shift() ;
    my $jobs    =   shift() ;
    my %options = %{shift()};
    # Iterate over expansion factors.
    foreach my $epoch ( @{$entry->{'resolution'}->{'epochs'}} ) {
	# Extract overdensity from the analysis file.
	my $analysisFile         = new PDL::IO::HDF5($entry->{'path'}."environment_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5");
	die("error: 'simulation0002' group is missing from file '".$entry->{'path'}."environment_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5'")
	    unless ( grep {$_ eq "simulation0002"} $analysisFile->groups() );
	foreach my $attributeName ( 'massTotal', 'convexHullOverdensity' ) {
	    die("error: '".$attributeName."' attribute is missing from file '".$entry->{'path'}."environment_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5'")
		unless ( grep {$_ eq $attributeName} $analysisFile->group('simulation0002')->attrs() );
	}
	(my $mass)               = $analysisFile->group('simulation0002')->attrGet('massTotal'            );
	(my $overdensity)        = $analysisFile->group('simulation0002')->attrGet('convexHullOverdensity');
	# Store these to the halo mass function.
	my $haloMassFunctionFile = new PDL::IO::HDF5(">".$entry->{'path'}."haloMassFunction_".$epoch->{'redshiftLabel'}.":MPI0000.hdf5");
	my $simulationGroup      = $haloMassFunctionFile->group('simulation0001')                                                       ;
	$simulationGroup->attrSet(massPrimary            => $entry->{$epoch->{'redshiftLabel'}}->{'massPrimary' });
	$simulationGroup->attrSet(massRegion             => $entry->{$epoch->{'redshiftLabel'}}->{'massRegion'  });
	$simulationGroup->attrSet(radiusRegion           => $entry->{$epoch->{'redshiftLabel'}}->{'radiusRegion'});
	$simulationGroup->attrSet(overdensityEnvironment => $overdensity                              );
	$simulationGroup->attrSet(massEnvironment        => $mass                                     );
    }
}
