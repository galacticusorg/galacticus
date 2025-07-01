#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use File::Copy;
use File::Find;
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::IO::HDF5;
use List::Util;
use List::ExtraUtils;
use List::MoreUtils;
use Data::Dumper;
use Storable qw(dclone);
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Galacticus::Constraints::Parameters;

# Construct progenitor halo mass function data from a variety of cosmological N-body simulations.
# Andrew Benson (21-October-2020)

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
     label               => "VSMDPL",
     description         => "Progenitor halo mass function for non-backsplash z=0 parent halos from the VSMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/vsmdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 6.2e6,
     snapshots           => "150 132 119 100 77 51",
     builder             => \&cosmoSimBuilder
 },
 {
     label               => "SMDPL",
     description         => "Progenitor halo mass function for non-backsplash z=0 parent halos from the SMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/smdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 9.63e7,
     snapshots           => "116 70 48 38 26 12",
     builder             => \&cosmoSimBuilder
 },
 {
     label               => "MDPL2",
     description         => "Progenitor halo mass function for non-backsplash z=0 parent halos from the MDPL2 simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/mdpl2/",
     hubbleConstant      => 0.6777,
     massParticle        => 1.51e9,
     snapshots           => "125 124 120 107 94 75 52 26",
     builder             => \&cosmoSimBuilder
 },
 {
     label               => "BigMDPL",
     description         => "Progenitor halo mass function for non-backsplash z=0 parent halos from the BigMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/bigmdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 2.359e10,
     snapshots           => "79 33 11 9 5 2",
     builder             => \&cosmoSimBuilder
 },
 {
     label               => "HugeMDPL",
     description         => "Progenitor halo mass function for non-backsplash z=0 parent halos from the HugeMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/hugemdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 7.9e10,
     snapshots           => "102 84 71 52 29 3",
     builder             => \&cosmoSimBuilder
 },
# {
#     label               => "Caterpillar",
#     description         => "Progenitor halo mass function for non-backsplash z=0 parent halos from the Caterpillar simulations.",
#     simulationReference => "Griffen et al.; 2016; ApJ; 818; 10",
#     simulationURL       => "https://www.caterpillarproject.org/",
#     hubbleConstant      => 0.6711,
#     builder             => \&caterpillarBuilder
# }
 {
    label               => "Symphony",
    description         => "Progenitor halo mass function for non-backsplash z=0 parent halos from the Symphony ZoomIn simulations.",
    simulationReference => "Nadler et al.; 2022;",
    simulationURL       => "https://web.stanford.edu/group/gfc/symphony/build/html/index.html",
    hubbleConstant      => 0.7,
    builder             => \&symphonyZoomInBuilder
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
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@{$jobsList})
	if ( scalar(@{$jobsList}) > 0 );
}

exit 0;

sub cosmoSimBuilder {
    my $simulation =   shift() ;
    my %options    = %{shift()};
    # Convert particle mass to Solar masses.
    $simulation->{'massParticle'} /= $simulation->{'hubbleConstant'};
    # Determine minimum and maximum parent halo masses for the mass function.
    $simulation->{'massParentMinimum'} = 10.0**(int(log($simulation->{'massParticle'})/log(10.0))+3);
    $simulation->{'massParentMaximum'} = 1.0e16;
    # Determine minimum and maximum progenitor mass ratios for the mass function.
    $simulation->{'massRatioProgenitorMinimum'} = 10.0**(int(log($simulation->{'massParticle'})/log(10.0))+2-16);
    $simulation->{'massRatioProgenitorMaximum'} = 10.0;
    # Determine parent and progenitor snapshots.
    my @snapshots = split(" ",$simulation->{'snapshots'});
    $simulation->{'snapshotParents'     } = &List::Util::max (                                                           @snapshots);
    $simulation->{'snapshotsProgenitors'} =              join(" ",map {$_ == $simulation->{'snapshotParents'} ? () : $_} @snapshots);
    # Construct the simulation path.
    $simulation->{'path'} = $options{'simulationDataPath'};
    $simulation->{'path'} .= "/"
	unless ( $simulation->{'path'} =~ m/\/$/ );
    $simulation->{'path'} .= "CosmoSim/".$simulation->{'label'}."/";
    # Identify non fly by halos and export these to IRATE format files.
    ## Parse the base parameters.
    my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/progenitorMassFunctionIdentifyNonFlyby.xml");
    ## Set the snapshots to select.
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[3]->{'selectedValues'}->{'value'} = $simulation->{'snapshots'};
    ## Remove hostedRootID selection.
    splice(@{$parameters->{'nbodyOperator'}->{'nbodyOperator'}},4,1);
    ## Iterate over subvolumes.
    for(my $i=0;$i<10;++$i) {
	for(my $j=0;$j<10;++$j) {
	    for(my $k=0;$k<10;++$k) {
		# Skip if the file exists.
		next
		    if ( -e $simulation->{'path'}."nonFlyby_progenitors_subVolume".$i."_".$j."_".$k.".hdf5" );
		# Modify file names.
		$parameters->{'nbodyImporter'}                        ->{'fileName'}->{'value'} = $simulation->{'path'}."tree_"                               .$i."_".$j."_".$k.".dat" ;
		$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[5]->{'fileName'}->{'value'} = $simulation->{'path'}."nonFlyby_progenitors_subVolume".$i."_".$j."_".$k.".hdf5";
		# Write parmeter file.
		my $parameterFileName = $simulation->{'path'}."identifyNonFlyby_progenitors_".$i."_".$j."_".$k.".xml";
		open(my $outputFile,">",$parameterFileName);
		print $outputFile $xml->XMLout($parameters, RootName => "parameters");
		close($outputFile);
		# Generate a job.
		my $job;
		$job->{'command'   } =
		    "./Galacticus.exe ".$parameterFileName;
		$job->{'launchFile'} = $simulation->{'path'}."identifyNonFlyby_progenitors_".$i."_".$j."_".$k.".sh" ;
		$job->{'logFile'   } = $simulation->{'path'}."identifyNonFlyby_progenitors_".$i."_".$j."_".$k.".log";
		$job->{'label'     } =                       "identifyNonFlyby_progenitors_".$i."_".$j."_".$k       ;
		$job->{'ppn'       } = 1;
		$job->{'nodes'     } = 1;
		$job->{'mpi'       } = "yes";
		push(@{$jobs->[0]},$job);
	    }
	}
    }
    # Iterate over simulations to construct the progenitor mass functions.
    foreach my $simulation ( @simulations ) {
	## Parse the base parameters.
	my $massFunctionParameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/progenitorMassFunctionCompute.xml");
	## Iterate over subvolumes.
	my @nbodyImporters;
	for(my $i=0;$i<10;++$i) {
	    for(my $j=0;$j<10;++$j) {
		for(my $k=0;$k<10;++$k) {
		    # Add an importer for this subvolume.
		    push(
			@nbodyImporters,
			{
			    value      => "IRATE"                                                                              ,
			    fileName   => {value => $simulation->{'path'}."nonFlyby_progenitors_subVolume".$i."_".$j."_".$k.".hdf5"},
			    properties => {value => "massVirial expansionFactor hostedRootID snapshotID"},
			    snapshot   => {value => "1"}
			}
			);
		}
	    }
	}
	## Compute the progenitor mass function.
	unless ( -e $simulation->{'path'}."progenitorsMassFunctions.hdf5" ) {
	    ## Modify parameters.
	    @{$massFunctionParameters->{'nbodyImporter'           }->{'nbodyImporter'}} = @nbodyImporters;
	    $massFunctionParameters  ->{'outputFileName'}                                                        ->{'value'} = $simulation->{'path'                      }."progenitorMassFunctions.hdf5";
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[0]->{'values'                    }->{'value'} = $simulation->{'massParticle'              }                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'description'               }->{'value'} = $simulation->{'description'               }                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'simulationReference'       }->{'value'} = $simulation->{'simulationReference'       }                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'simulationURL'             }->{'value'} = $simulation->{'simulationURL'             }                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'massParentMinimum'         }->{'value'} = $simulation->{'massParentMinimum'         }                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'massParentMaximum'         }->{'value'} = $simulation->{'massParentMaximum'         }                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'massRatioProgenitorMinimum'}->{'value'} = $simulation->{'massRatioProgenitorMinimum'}                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'massRatioProgenitorMaximum'}->{'value'} = $simulation->{'massRatioProgenitorMaximum'}                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'snapshotParents'           }->{'value'} = $simulation->{'snapshotParents'           }                               ;
	    $massFunctionParameters  ->{'nbodyOperator'           }->{'nbodyOperator'}->[1]->{'snapshotsProgenitors'      }->{'value'} = $simulation->{'snapshotsProgenitors'      }                               ;
	    ## Write the parameter file.
	    my $parameterFileName = $simulation->{'path'}."progenitorMassFunctions.xml";
	    open(my $outputFile,">",$parameterFileName);
	    print $outputFile $xml->XMLout($massFunctionParameters, RootName => "parameters");
	    close($outputFile);
	    ## Construct the job.
	    my $job;
	    $job->{'command'   } =
		"./Galacticus.exe ".$parameterFileName;
	    $job->{'launchFile'} = $simulation->{'path'}."progenitorMassFunctions.sh" ;
	    $job->{'logFile'   } = $simulation->{'path'}."progenitorMassFunctions.log";
	    $job->{'label'     } =                       "progenitorMassFunctions"    ;
	    $job->{'ppn'       } = 16;
	    $job->{'nodes'     } = 1;
	    $job->{'mpi'       } = "no";
	    $job->{'onCompletion'} =
	    {
		function  => \&copyFile,
		arguments => [ $simulation->{'path'}."progenitorMassFunctions:MPI0000.hdf5", $ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/progenitorMassFunctions_".$simulation->{'label'}.".hdf5" ]
	    };
	    push(@{$jobs->[1]},$job);
	}
    }
}

sub caterpillarBuilder {
    # Build jobs for the Caterpillar simulations.
    my $simulation =   shift() ;
    my %options    = %{shift()};
    # Construct the simulation path.
    $simulation->{'path'} = $options{'simulationDataPath'};
    $simulation->{'path'} .= "/"
	unless ( $simulation->{'path'} =~ m/\/$/ );
    $simulation->{'path'} .= "Caterpillar/";
    # Set snapshots to process for each resolution, corresponding to z ~ 0.0, 0.02, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0.
    my $snapshots = "319 314 295 232 189 145 111 70";
    # Read the parent halo file.
    (my $parentId, my $LX, my $zoomId, my $min2, my $mass, my $rvir, my $badFlag, my $badSubF) = rcols($simulation->{'path'}."parent_zoom_index.txt",1,3,5,6,10,11,12,13,{COLSEP => qr/\s*\|\s*/, LINES => "1:"});
    # Find the resolution levels available for each parent along with the highest resolution level achieved.
    my $parents;
    for(my $i=0;$i<nelem($parentId);++$i) {
	my $parentId_ = $parentId->(($i))->sclr();
	$parents->{$parentId_}->{'ID'} = $parentId_;
	# Store virial mass.
	$parents->{$parentId_}->{'mass'}->{$LX->(($i))->sclr()} = $mass->(($i))->sclr();
	# Append to the available levels.
	push(@{$parents->{$parentId_}->{'levels'}},$LX->(($i))->sclr());
	# Find the final snapshot for this case.
	my $snapshotFinal;
	my $parentPath      = $simulation->{'path'}."H".$parentId_."_LX".$LX->(($i))->sclr();
	my $gotHostedRootId = 0;
	opendir(my $parentDirectory,$parentPath);
	while ( my $treeFileName = readdir($parentDirectory) ) {
	    next
		unless ( $treeFileName =~ m/^tree_\d_\d_\d\.dat/ );
	    my $gotSnapshotFinal = 0;
	    open(my $treeFile,$parentPath."/".$treeFileName);
	    while ( my $line = <$treeFile> ) {
		next
		    if ( $line =~ m/^#/ || $line =~ m/^\d+$/ );
		my @columns = split(" ",$line);
		if ( scalar(@columns) >= 29 &&  $columns[0] == 1.0000 && $columns[30] == $zoomId->(($i))->sclr() ) {
		    $parents->{$parentId_}->{$LX->(($i))->sclr()}->{'hostedRootId'} = $columns[1];
		    $gotHostedRootId = 1;
		}
		if ( $columns[0] == 1.0000 && ! $gotSnapshotFinal ) {
		    $snapshotFinal = $columns[31];
		    push(@{$parents->{$parentId_}->{'treeFiles'}->{$LX->(($i))->sclr()}},$treeFileName);
		    $gotSnapshotFinal = 1;
		}
		last
		    if ( $gotSnapshotFinal && $gotHostedRootId );
	    }
	    close($treeFile);
	}
	closedir($parentDirectory);
	$parents->{$parentId_}->{'snapshotFinal'}->{$LX->(($i))->sclr()} = $snapshotFinal;
	die("unable to find hostedRootId for parent ".$parentId->(($i))->sclr()." LX".$LX->(($i))->sclr()." zoomID=".$zoomId->(($i))->sclr())
	    unless ( $gotHostedRootId );	
    }
    my $levelMaxMin   = 20;
    my $levelMaxCount     ;
    foreach my $parent ( &List::ExtraUtils::hashList($parents) ) {
	$parent->{'levelMax'} = &List::Util::max(@{$parent->{'levels'}});
	++$levelMaxCount->{$parent->{'levelMax'}};
	$levelMaxMin = $parent->{'levelMax'}
	if ( $parent->{'levelMax'} < $levelMaxMin );
    }
    # Initialize lists of parent halo masses.
    my $massParent;
    # For each parent, construct a job to identify non-flyby progenitors.
    ## Parse the base parameters.
    my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/progenitorMassFunctionIdentifyNonFlyby.xml");
    ## Adjust comsological model.
    $parameters->{'cosmologyParameters'}->{'HubbleConstant' }->{'value'} = 67.11;
    $parameters->{'cosmologyParameters'}->{'OmegaMatter'    }->{'value'} =  0.32;
    $parameters->{'cosmologyParameters'}->{'OmegaDarkEnergy'}->{'value'} =  0.68;
    $parameters->{'cosmologyParameters'}->{'OmegaBaryon'    }->{'value'} =  0.05;
    ## Iterate over parents.
    foreach my $parent ( &List::ExtraUtils::hashList($parents) ) {
	# Construct the file name.
	my $parentDirectoryName = "H".$parent->{'ID'}."_LX".$parent->{'levelMax'}."/";
	my $outputFileName = $simulation->{'path'}.$parentDirectoryName."identifyNonFlyby_progenitors.hdf5";
	# Skip if the file already exists.
	next
	    if ( -e $outputFileName );
	# Clone parameters.
	my $parameters_ = dclone($parameters);
	# Set output file.
	$parameters_->{'nbodyOperator'}->{'nbodyOperator'}->[6]->{'fileName'}->{'value'} = $outputFileName;
	# Set the snapshots to select.
	$parameters_->{'nbodyOperator'}->{'nbodyOperator'}->[3]->{'selectedValues'}->{'value'} = $snapshots;
	# Set the hostedRootID to select.
	$parameters_->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'selectedValues'}->{'value'} = $parent->{$parent->{'levelMax'}}->{'hostedRootId'};
	# Add an operator to shift snapshot numbers. All Caterpillar simulations should have 320 snapshots, except for some
	# early-run simulations which had 256 (64 extra early snapshots were added to later-run simulations). There are fewer
	# snapshots in the Rockstar trees because Rockstar ignores snapshots with no halos and renumbers. So, here we want to
	# adjust the snapshot numbers such that a=1 is always snapshot 319.
	my $snapshotOffset = 319-$parent->{'snapshotFinal'}->{$parent->{'levelMax'}};
	my $shiftOperator =
	{
	    value        =>           "shiftProperty" ,
	    propertyName => {value => "snapshotID"   },
	    shiftBy      => {value => $snapshotOffset}
	};
	unshift(@{$parameters_->{'nbodyOperator'}->{'nbodyOperator'}},$shiftOperator);
	# Add all tree files.
	$parameters_->{'nbodyImporter'} =
	{
	    value => "merge"
	};
	@{$parameters_->{'nbodyImporter'}->{'nbodyImporter'}} =
	    map 
	{
	    {
		value       => "rockstar"                                                    ,
		fileName    => {value => $simulation->{'path'}.$parentDirectoryName.$_      },
		readColumns => {value => "id scale desc_id pid upid mmp Mvir scale Snap_num"}
	    }
	}
	@{$parent->{'treeFiles'}->{$parent->{'levelMax'}}};
	# Write parmeter file.
	my $parameterFileName = $simulation->{'path'}.$parentDirectoryName."identifyNonFlyby_progenitors.xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($parameters_, RootName => "parameters");
	close($outputFile);
	# Generate a job.
	my $job;
	$job->{'command'   } =
	    "./Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $simulation->{'path'}.$parentDirectoryName."identifyNonFlyby_progenitors.sh" ;
	$job->{'logFile'   } = $simulation->{'path'}.$parentDirectoryName."identifyNonFlyby_progenitors.log";
	$job->{'label'     } = "Caterpillar_H".$parent->{'ID'}."_LX".$parent->{'levelMax'};
	$job->{'ppn'       } = 1;
	$job->{'nodes'     } = 1;
	$job->{'mem'       } = "16gb";
	$job->{'mpi'       } = "yes";
	push(@{$jobs->[0]},$job);
	# Add the parent mass to the relevant list.
	push(@{$massParent->{$parent->{'levelMax'}}},$parent->{'mass'}->{$parent->{'levelMax'}});
    }
    # Write files of parent halo masses.
    foreach my $resolution ( keys(%{$massParent}) ) {
	my $masses     = pdl @{$massParent->{$resolution}};
	my $weights    = pdl ones(nelem($masses));
	my $massesFile = new PDL::IO::HDF5(">".$simulation->{'path'}."masses_LX".$resolution.".hdf5");
	$massesFile->dataset('treeRootMass')->set($masses );
	$massesFile->dataset('treeWeight'  )->set($weights);
    }
    # Generate jobs to construct mass functions.
    ## Parse the base parameters.
    my $massFunctionParameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/progenitorMassFunctionCompute.xml");
    ## Iterate over resolutions.
    foreach my $resolution ( keys(%{$massParent}) ) {
	# Set output file name.
	my $outputFileName = $simulation->{'path'}."progenitorMassFunctions_LX".$resolution.".hdf5";
	next
	    if ( -e $simulation->{'path'}."progenitorMassFunctions_LX".$resolution.":MPI0000.hdf5" );
	# Add an importer for each parent.
	@{$massFunctionParameters->{'nbodyImporter'}->{'nbodyImporter'}} =
	    map
	    {
		$_->{'levelMax'} == $resolution
		    ?
		    
		{
		    value      => "IRATE",
		    fileName   => {value => $simulation->{'path'}."H".$_->{'ID'}."_LX".$_->{'levelMax'}."/identifyNonFlyby_progenitors.hdf5"},
		    properties => {value => "massVirial expansionFactor hostedRootID snapshotID"},
		    snapshot   => {value => "1"}
		}
		:
		    ()
	    }
	&List::ExtraUtils::hashList($parents);
	# Determine minimum and maximum progenitor mass ratios for the mass function.
	my $massParticle;
	if      ( $resolution == 15 ) {
	    $massParticle = 3.7317e3;
	} elsif ( $resolution == 14 ) {
	    $massParticle = 2.9854e4;
	} elsif ( $resolution == 13 ) {
	    $massParticle = 2.3883e5;
	} elsif ( $resolution == 12 ) {
	    $massParticle = 1.9106e6;
	} elsif ( $resolution == 11 ) {
	    $massParticle = 1.5285e7;
	} elsif ( $resolution == 10 ) {
	    $massParticle = 1.2228e8;
	} else {
	    die("unknown resolution");
	}
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[0]->{'values'                    }->{'value'} = $massParticle;
	# Determine minimum and maximum parent halo masses for the mass function.
	my $massParentMinimum = 1.0e11;
	my $massParentMaximum = 1.0e13;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massParentMinimum'         }->{'value'} = $massParentMinimum;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massParentMaximum'         }->{'value'} = $massParentMaximum;
	# Define minimum and maximum mass ratios
	my $massRatioProgenitorMinimum = 10.0**(int(log($massParticle)/log(10.0))+2-log10($massParentMaximum));
	my $massRatioProgenitorMaximum = 10.0;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massRatioProgenitorMinimum'}->{'value'} = $massRatioProgenitorMinimum;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massRatioProgenitorMaximum'}->{'value'} = $massRatioProgenitorMaximum;
	# Set snapshots.
	my @snapshotList         = split(" ",$snapshots);
	my $snapshotParents      = &List::Util::max (                                           @snapshotList);
	my $snapshotsProgenitors =              join(" ",map {$_ == $snapshotParents ? () : $_} @snapshotList);
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'snapshotParents'           }->{'value'} = $snapshotParents;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'snapshotsProgenitors'      }->{'value'} = $snapshotsProgenitors;
	# Set other task properties.
	$massFunctionParameters  ->{'outputFileName'}                                                              ->{'value'} = $outputFileName;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'description'               }->{'value'} = $simulation->{'description'        }                                                 ;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'simulationReference'       }->{'value'} = $simulation->{'simulationReference'}                                                 ;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'simulationURL'             }->{'value'} = $simulation->{'simulationURL'      }                                                 ;
	## Write the parameter file.
	my $parameterFileName = $simulation->{'path'}."progenitorMassFunctions_LX".$resolution.".xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($massFunctionParameters, RootName => "parameters");
	close($outputFile);
	## Construct the job.
	my $job;
	$job->{'command'   } =
	    "./Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $simulation->{'path'}."progenitorMassFunctions_LX".$resolution.".sh" ;
	$job->{'logFile'   } = $simulation->{'path'}."progenitorMassFunctions_LX".$resolution.".log";
	$job->{'label'     } =                       "progenitorMassFunctions_LX".$resolution       ;
	$job->{'ppn'       } = 16;
	$job->{'nodes'     } = 1;
	$job->{'mpi'       } = "no";
	$job->{'onCompletion'} =
	{
	    function  => \&copyFile,
	    arguments => [ $simulation->{'path'}."progenitorMassFunctions_LX".$resolution.":MPI0000.hdf5", $ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/progenitorMassFunctions_".$simulation->{'label'}."_LX".$resolution.".hdf5" ]
	};
	push(@{$jobs->[1]},$job);
    }   
}

sub symphonyZoomInBuilder {
    # Build jobs for the Caterpillar simulations.
    my $simulation =   shift() ;
    my %options    = %{shift()};
    # Construct the simulation path.
    $simulation->{'path'} = $options{'simulationDataPath'};
    $simulation->{'path'} .= "/"
        unless ( $simulation->{'path'} =~ m/\/$/ );
    $simulation->{'path'} .= "Symphony_ZoomIns/";

    # Set snapshots to process for each resolution, corresponding to z ~ 0.0, 0.02, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0.
    # expansion factors:  1.0000,   0.66503,   0.50239,   0.32987,   0.20064 (grabs from snapshot #)
    my $snapshots = "235  203  181 148 109";
    
    # Initialize lists of parent halo masses.
    my $massParent;
    # For each parent, construct a job to identify non flyby progenitors.
    ## Parse the base parameters.
    my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/progenitorMassFunctionIdentifyNonFlyby.xml");
    ## Adjust comsological model.
    $parameters->{'cosmologyParameters'}->{'HubbleConstant' }->{'value'} = 70.00;
    $parameters->{'cosmologyParameters'}->{'OmegaMatter'    }->{'value'} =  0.286;
    $parameters->{'cosmologyParameters'}->{'OmegaDarkEnergy'}->{'value'} =  0.714;
    $parameters->{'cosmologyParameters'}->{'OmegaBaryon'    }->{'value'} =  0.047;
    # List available simulations
    my @models;
    find( (sub {my @path = split(/\//,$File::Find::dir);push(@models,$_ eq "tree_0_0_0.dat" ? {simulation => $path[-2], realization => $path[-1]}  : ())}), ($options{'simulationDataPath'}."/Symphony_ZoomIns") );
    my $xml         = new XML::Simple();
    my $hostHaloIDs = $xml->XMLin("constraints/pipelines/darkMatter/symphonyZoomInHostHaloIDs.xml");
    foreach my $model ( @models ) {
	$model->{'hostHaloID'} = $hostHaloIDs->{$model->{'simulation'}}->{$model->{'realization'}};
	print "Model: ".$model->{'simulation'}."; Realization: ".$model->{'realization'}."; HostHaloID: ".$model->{'hostHaloID'}."\n";	
	my $parentDirectoryName = $model->{'simulation'}."/".$model->{'realization'}."/";
	my $outputFileName = $simulation->{'path'}.$parentDirectoryName."identifyNonFlyby_progenitors.hdf5";
	# Skip if the file already exists.
	next
	    if ( -e $outputFileName );
	# Read the tree file.
	open(my $treeFile,$simulation->{'path'}.$parentDirectoryName."tree_0_0_0.dat");
	while ( my $line = <$treeFile> ) {
	    # Skip header lines.
	    next
		if ( $line =~ m/^#/ || $line =~ m/^\d+\s*$/ );
	    # Split the line into columns.
	    my @columns = split(" ",$line);
	    # Halo ID is in the 2nd column, halo mass in the 11th.
	    # Check for a match to the host halo ID, and record the mass.
	    if ( $columns[1] == $model->{'hostHaloID'} ) {
		$model->{'hostHaloMass'} = $columns[10]/$simulation->{'hubbleConstant'};
		# We can exit reading the file now as we've found the host halo of interest.
		last;
	    }
	}
	close($treeFile);
	# Check that we got the host halo mass.
	print $model->{'simulation'}."\t".$model->{'realization'}."\t".$model->{'hostHaloMass'}."\n";
	die("unable to find host halo mass")
	    unless ( exists($model->{'hostHaloMass'}) );
	# Clone parameters.
	my $parameters_ = dclone($parameters);
	# Set output file.
	$parameters_->{'nbodyOperator'}->{'nbodyOperator'}->[6]->{'fileName'}->{'value'} = $outputFileName;
	# Set the snapshots to select.
	$parameters_->{'nbodyOperator'}->{'nbodyOperator'}->[3]->{'selectedValues'}->{'value'} = $snapshots;
	# Set the hostedRootID to select.
	$parameters_->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'selectedValues'}->{'value'} = $model->{'hostHaloID'};
	
	# Add all tree files.
	$parameters_->{'nbodyImporter'} =
	{
	    value       => "rockstar"                                                    ,
	    fileName    => {value => $simulation->{'path'}.$parentDirectoryName."tree_0_0_0.dat"   },
	    readColumns => {value => "id scale desc_id pid upid mmp Mvir Rvir X Y Z scale Snap_num"}
	};
	# Write parameter file.
	my $parameterFileName = $simulation->{'path'}.$parentDirectoryName."identifyNonFlyby_progenitors.xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($parameters_, RootName => "parameters");
	close($outputFile);
	# Generate a job.
	my $job;
	$job->{'command'   } =
	    "./Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $simulation->{'path'}.$parentDirectoryName."identifyNonFlyby_progenitors.sh" ;
	$job->{'logFile'   } = $simulation->{'path'}.$parentDirectoryName."identifyNonFlyby_progenitors.log";
	$job->{'label'     } = "Symphony".$model->{'simulation'}." ".$model->{'realization'};
	$job->{'ppn'       } = 1;
	$job->{'nodes'     } = 1;
	$job->{'mem'       } = "16gb";
	$job->{'mpi'       } = "yes";
	push(@{$jobs->[0]},$job);
    }

    ## NOTE: We want to group the models such that all realizations of "MilkyWay" simulations for example as in one "group". We'll
    ##       want to combine the progenitor mass functions of all realizations in a group. In the following, we first make a list
    ##       of unique group names, and then create a dictionary of group names in which each entry is a list of the available
    ##       realizations for that group.    
    # Create a list of model groups (e.g. "MilkyWay", "MilkyWay_WDM3", etc.).
    my @modelGroups = &List::MoreUtils::uniq(sort(map {$_->{'simulation'}} @models));
    # Create a dictionary of lists of available realizations in each model group.
    my %modelRealizations;
    foreach my $modelGroup ( @modelGroups ) {
	@{$modelRealizations{$modelGroup}} = map {$_->{'simulation'} eq $modelGroup ? $_->{'realization'} : ()} @models;
    }
    # Write files of parent halo masses.
    foreach my $modelGroup ( @modelGroups ) {
	my $masses     = pdl [ map {$_->{'simulation'} eq $modelGroup ? $_->{'hostHaloMass'} : ()} @models ];
	my $weights    = pdl ones(nelem($masses));
	my $massesFile = new PDL::IO::HDF5(">".$simulation->{'path'}."hostHaloMasses_".$modelGroup.".hdf5");
	$massesFile->dataset('treeRootMass')->set($masses );
	$massesFile->dataset('treeWeight'  )->set($weights);
    }

    # Generate jobs to construct mass functions.
    ## Parse the base parameters.
    my $massFunctionParameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/progenitorMassFunctionCompute.xml");

    ## Iterate over the "modelGroup". We want to generate one job for each model group - that job will combine all realizations in that model group.
    foreach my $modelGroup ( @modelGroups ) {
	# Set output file name.
	my $outputFileName = $simulation->{'path'}."progenitorMassFunctions_".$modelGroup.".hdf5";
	next
            if ( -e $simulation->{'path'}."progenitorMassFunctions_".$modelGroup.":MPI0000.hdf5" );
	# Add an importer for each parent.
	@{$massFunctionParameters->{'nbodyImporter'}->{'nbodyImporter'}} =
            map
	{{
    	    value      => "IRATE",
    	    fileName   => {value => $simulation->{'path'}.$modelGroup."/".$_."/identifyNonFlyby_progenitors.hdf5"},
    	    properties => {value => "massVirial expansionFactor hostedRootID snapshotID"},
    	    snapshot   => {value => "1"}
    	    }} @{$modelRealizations{$modelGroup}};
	
	# Determine minimum and maximum progenitor mass ratios for the mass function.
	my $massParticle;
	if      ( $modelGroup =~ m/LMC/ ) {
	    $massParticle = 3.5247625e4;
	} elsif ( $modelGroup =~ m/MilkyWay/ ) {
	    if ( $modelGroup =~ m/hires/ ){
            	$massParticle = 3.5247625e4;
	    } else {
    	    	$massParticle = 2.81981e5;
	    }
	} else {
	    die("unknown resolution");
	}
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[0]->{'values'                    }->{'value'} = $massParticle;
	# Determine minimum and maximum parent halo masses for the mass function.
	my $massParentMinimum;
	my $massParentMaximum;
	if ($modelGroup =~ m/MilkyWay/){
	    $massParentMinimum = 1.0e11;
	    $massParentMaximum = 1.0e13;
	}
	elsif ($modelGroup =~ m/LMC/) {
	    $massParentMinimum = 1.0e10;
	    $massParentMaximum = 1.0e12;
	}
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massParentMinimum'         }->{'value'} = $massParentMinimum;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massParentMaximum'         }->{'value'} = $massParentMaximum;
	# Define minimum and maximum mass ratios
	my $massRatioProgenitorMinimum = 10.0**(int(log($massParticle)/log(10.0))+2-log10($massParentMaximum));
	my $massRatioProgenitorMaximum = 10.0;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massRatioProgenitorMinimum'}->{'value'} = $massRatioProgenitorMinimum;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'massRatioProgenitorMaximum'}->{'value'} = $massRatioProgenitorMaximum;
	# Set snapshots.
	my @snapshotList         = split(" ",$snapshots);
	my $snapshotParents      = &List::Util::max (                                           @snapshotList);
	my $snapshotsProgenitors =              join(" ",map {$_ == $snapshotParents ? () : $_} @snapshotList);
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'snapshotParents'           }->{'value'} = $snapshotParents;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'snapshotsProgenitors'      }->{'value'} = $snapshotsProgenitors;
	# Set other task properties.
	$massFunctionParameters  ->{'outputFileName'}                                                              ->{'value'} = $outputFileName;
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'description'               }->{'value'} = $simulation->{'description'        };
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'simulationReference'       }->{'value'} = $simulation->{'simulationReference'};
	$massFunctionParameters  ->{'nbodyOperator'     }->{'nbodyOperator'}->[1]->{'simulationURL'             }->{'value'} = $simulation->{'simulationURL'      };
	## Write the parameter file.
	my $parameterFileName = $simulation->{'path'}."progenitorMassFunctions_".$modelGroup.".xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($massFunctionParameters, RootName => "parameters");
	close($outputFile);
	## Construct the job.
	my $job;
	$job->{'command'   } =
	    "./Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $simulation->{'path'}."progenitorMassFunctions_".$modelGroup.".sh" ;
	$job->{'logFile'   } = $simulation->{'path'}."progenitorMassFunctions_".$modelGroup.".log";
	$job->{'label'     } =                       "progenitorMassFunctions_".$modelGroup    ;
	$job->{'ppn'       } = 16;
	$job->{'nodes'     } = 1;
	$job->{'mpi'       } = "no";
	$job->{'onCompletion'} =
	{
	    function  => \&copyFile,
	    arguments => [ $simulation->{'path'}."progenitorMassFunctions_".$modelGroup.":MPI0000.hdf5", $ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/progenitorMassFunctions_".$simulation->{'label'}."_".$modelGroup.".hdf5" ]
	};
	push(@{$jobs->[1]},$job);
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
