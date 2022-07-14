#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use Data::Dumper;
use Galacticus::Options;

# Run a set of merger trees with forests split and not split. Compare the results which should be identical.
# Andrew Benson (19-May-2016)

# Test is currently being skipped until MPI implementation is completed (split forest processing is currently broken).
print "SKIPPED: split forest processing is currently broken until MPI implementation is completed\n";
exit 0;

# Global status for all tests.
my $statusGlobal = 0;
# Define a set of merger tree forest files that we want to process.
my @forestFiles =
    (
     # milli-Millennium forest number 12000000 - a modestly sized forest which can be split and has sufficiently complex structure
     # that it has uncovered bugs in split forest handling in the past.
     {
     	 label       => "milliMillennium-forest12000000-noMerging",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree12000000.hdf5",
     	 parameters  => 
     	 {
     	     mergerTreeConstructor => {satelliteMergingTimescalesSubresolution => {value => "infinite"}}
     	 }
     },
     {
     	 label       => "milliMillennium-forest12000000-merging",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree12000000.hdf5",
     	 skipOrphans => 1,     
     	 parameters  => 
     	 {
     	     virialOrbit                                 => {value => "fixed"            },
     	     mergerTreeConstructor => {
		 satelliteMergingTimescalesSubresolution => {value => "boylanKolchin2008"},
		 presetMergerTimes                       => {value => "true"             },
		 presetMergerNodes                       => {value => "false"            }
	     }
     	 }
     },
     {
     	 label       => "milliMillennium-forest12000000-mergingPresetTargets",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree12000000.hdf5",
     	 skipOrphans => 1,     
     	 parameters  => 
     	 {
     	     virialOrbit                                 => {value => "fixed"            },
     	     mergerTreeConstructor => {
		 satelliteMergingTimescalesSubresolution => {value => "boylanKolchin2008"},
		 presetMergerTimes                       => {value => "true"             },
		 presetMergerNodes                       => {value => "true"             }
	     }
     	 }
     },
     # A milli-Millennium tree which fails unless we (correctly) nullify the merge target pointer of a node which never merges but
     # which originated in a tree in the forest where the original branch it should merge with terminated due to an inter-tree
     # transfer.
     {
     	 label       => "milliMillennium-forest3000055000000-mergingPresetTargets",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree3000055000000.hdf5",
     	 parameters  => 
     	 {
     	     virialOrbit                                 => {value => "fixed"            } ,
     	     mergerTreeConstructor => {
		 satelliteMergingTimescalesSubresolution => {value => "boylanKolchin2008"},
		 presetMergerTimes                       => {value => "true"             },
		 presetMergerNodes                       => {value => "true"             }		 
	     }
	 }
     },
     # A milli-Millennium tree which fails if the "final time in tree" used in limiting node evolution is computed across all
     # trees in a forest rather than just in the local tree.
      {
      	 label       => "milliMillennium-forest2000024000000-noMerging",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree2000024000000.hdf5",
     	 parameters  => 
     	  {
     	      mergerTreeConstructor => {
		  satelliteMergingTimescalesSubresolution => {value => "infinite"},
		  presetMergerNodes                       => {value => "true"    }
	      }
     	  }
     },
     # A milli-Millennium tree in which a mergee is transferred as part of an inter-tree event but the destination node is a
     # clone. Requires that promotion of a cloned node to a satellite parent correctly moves any satellites in the clone into the
     # host halo of the parent.
     {
     	 label       => "milliMillennium-forest1000019000000-noMerging",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree1000019000000.hdf5",
     	 parameters  => 
     	  {
     	      mergerTreeConstructor => {
		  satelliteMergingTimescalesSubresolution => {value => "infinite"},
		  presetMergerNodes                       => {value => "true"    }
	      }
     	  }
     },
     # A milli-Millennium tree in which a satellite is orphanized with a merge target that is a satellite. Fails unless we
     # explicitly check for that merge target being a satellite and move the new host to the merge target's parent.
     {
     	 label       => "milliMillennium-forest1000020000000-noMerging",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree1000020000000.hdf5",
     	 parameters  => 
     	  {
     	      mergerTreeConstructor => {satelliteMergingTimescalesSubresolution => {value => "infinite"}},
     	      mergerTreeConstructor => {presetMergerNodes                       => {value => "true"    }}
     	  }
     },
     # A milli-Millennium tree which fails unless orphaned mergees are reassigned to their merge target branch in cases where their
     # parent undergoes a node merger.
     {
     	 label       => "milliMillennium-forest4000027000000-noMerging",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree4000027000000.hdf5",
     	 parameters  => 
     	  {
     	      mergerTreeConstructor => {
		  satelliteMergingTimescalesSubresolution => {value => "infinite"},
		  presetMergerNodes                       => {value => "true"    }
	      }
     	  }
     },
     # A milli-Millennium tree which fails (due to an orphan galaxy merging in the split case, but not merging in the unsplit case)
     # unless orphaned galaxies are ignored.
     {
     	 label       => "milliMillennium-forest2000035000000-merging",
     	 fileName    => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree2000035000000.hdf5",
     	 skipOrphans => 1,     
     	 parameters  => 
	 {
	     virialOrbit                                 => {value => "fixed"            },
	     mergerTreeConstructor => {
		 satelliteMergingTimescalesSubresolution => {value => "boylanKolchin2008"},
		 presetMergerTimes                       => {value => "true"             },
		 presetMergerNodes                       => {value => "false"            }
	     }
	 }
     }
    );

# Find the location of the Millennium Database data.
my $skipMillennium           = 0;
my $millenniumDatabaseConfig = &Galacticus::Options::Config('millenniumDB');
if ( defined($millenniumDatabaseConfig) ) {
    if ( exists($millenniumDatabaseConfig->{'path'}) ) {
	print $millenniumDatabaseConfig->{'path'}."\n";
	unless ( -e $millenniumDatabaseConfig->{'path'}."/milliMillennium/milliMillennium.hdf5" ) {
	    print "SKIPPED: milli-Millennium data not available\n";
	    $skipMillennium = 1;
	}	    
    } else {
	print "SKIPPED: Millennium database data not available\n";
	$skipMillennium = 1;
    }
} else {
    print "SKIPPED: Millennium database location undefined\n";
    $skipMillennium = 1;
}
push(
    @forestFiles,
    {
 	label       => "milliMillennium",
 	fileName    => $millenniumDatabaseConfig->{'path'}."/milliMillennium/milliMillennium.hdf5",
	skipOrphans => 1,     
	parameters  => 
	{
	    mergerTreeConstructor => {satelliteMergingTimescalesSubresolution => {value => "boylanKolchin2008"}},
	    virialOrbit                                                       => {value => "fixed"            } ,
	    mergerTreeConstructor => {presetMergerTimes                       => {value => "true"             }},
	    mergerTreeConstructor => {presetMergerNodes                       => {value => "true"             }}
	}
    }
    )
    unless ( $skipMillennium );
    
# Define the set of properties that we want to compare.
my @properties = ( "nodeIndex", "positionPositionX", "positionPositionY", "positionPositionZ", "positionVelocityX", "positionVelocityY", "positionVelocityZ", "satelliteNodeIndex", "satelliteMergeTime", "satelliteBoundMass", "nodeIsIsolated", "basicTimeLastIsolated", "parentIndex", "basicMass", "satelliteStatus" );

# Define the two types of model we want to run.
my @types = ( 'split', 'unsplit' );

# Create directories needed.
system("mkdir -p outputs/test-splitForests");

# Locate a scratch directory.
my $scratchConfig = &Galacticus::Options::Config('scratch');

# Iterate over sets of merger tree forest files.
foreach my $forestFile ( @forestFiles ) {
    # Write initial report.
    my $status = 0;
    print "\n\nRunning split forest tests for forest: ".$forestFile->{'label'}."\n";
    # Determine output directory for these forests.
    my $outputDirectoryName = "outputs/test-splitForests/".$forestFile->{'label'}."/";
    system("mkdir -p ".$outputDirectoryName);
    # Iterate over models, running them, reading in their data, and generating sort indices into the node indices.
    my $data;
    foreach ( @types ) {
	unless ( -e $outputDirectoryName.$_.'.hdf5' ) {
	    print "Running ".$_." model...\n";
	    # Read and modify parameter file.
	    my $xml        = new XML::Simple(RootName => "parameters");
	    my $parameters = $xml->XMLin("parameters/test-splitForests-".$_.".xml");
	    $parameters->{'outputFileName'}                 ->{'value'} = "testSuite/".$outputDirectoryName.$_.".hdf5";
	    $parameters->{'mergerTreeConstructor'   }->{'fileNames'  }->{'value'} = $forestFile->{'fileName'};
	    $parameters->{'task'                    }->{'suspendPath'}->{'value'} = defined($scratchConfig) ? $scratchConfig->{'path'} : ".";
	    my @stack  = ( { node => $parameters, content => $forestFile->{'parameters'} } );
	    while ( scalar(@stack) > 0 ) {
		my $entry = pop(@stack);
		foreach ( keys(%{$entry->{'content'}}) ) {
		    if ( $_ eq "value" ) {
			$entry->{'node'}->{'value'} = $entry->{'content'}->{'value'};
		    } else {
			$entry->{'node'}->{$_}->{'value'} = ""
			    unless ( exists($entry->{'node'}->{$_}) );
			push(@stack,{node => $entry->{'node'}->{$_}, content => $entry->{'content'}->{$_}});
		    }
		}
	    }
	    my $parameterFileName = $outputDirectoryName.$_.".xml";
	    open(my $parameterFile,">".$parameterFileName);
	    print $parameterFile $xml->XMLout($parameters);
	    close($parameterFile);
	    # Run the model.
	    system("cd ..; ./Galacticus.exe testSuite/".$parameterFileName);
	    unless ( $? == 0 ) { 
	    	print "FAILED {".$forestFile->{'label'}."}: model '".$_."' did not complete\n";
	    	$status       = 1;
	    	$statusGlobal = 1;
	    }
	}
	if ( $status == 0 ) {
	    $data->{$_}->{'file' } = new PDL::IO::HDF5($outputDirectoryName.$_.'.hdf5');
	    $data->{$_}->{'nodes'} = $data->{$_}->{'file'}->group("Outputs/Output1/nodeData");
	    foreach my $property ( @properties ) {
		$data->{$_}->{'properties'}->{$property} = $data->{$_}->{'nodes'}->dataset($property)->get();
	    }
	    my $selection;
	    if ( exists($forestFile->{'skipOrphans'}) && $forestFile->{'skipOrphans'} ) {
		$selection = which($data->{$_}->{'properties'}->{'satelliteStatus'} != 2);
	    } else {
		$selection = pdl sequence(nelem($data->{$_}->{'properties'}->{'nodeIndex'}));
	    }
	    $data->{$_}->{'selectionSortIndex'} = $data->{$_}->{'properties'}->{'nodeIndex'}->($selection)->qsorti();
	    $data->{$_}->{'rank'              } = $selection->($data->{$_}->{'selectionSortIndex'});
	    # Get tree indices.
	    foreach my $treeDatasetName ( "mergerTreeIndex", "mergerTreeStartIndex", "mergerTreeCount" ) {
		$data->{$_}->{$treeDatasetName} = $data->{$_}->{'file'}->dataset("Outputs/Output1/".$treeDatasetName)->get();
	    }
	}
    }
    next
	unless ( $status == 0 );
    # Test for equal numbers of nodes.
    my $orphansText = exists($forestFile->{'skipOrphans'}) && $forestFile->{'skipOrphans'} ? " (ignoring orphaned galaxies) " : "";
    print "Testing for equal numbers of nodes".$orphansText."...\n";
    unless ( nelem($data->{'split'}->{'rank'}) == nelem($data->{'unsplit'}->{'rank'}) ) {
	for(my $i=0;$i<nelem($data->{'split'}->{'rank'});++$i) {
	    print "In split but not unsplit: ".$data->{'split'}->{'properties'}->{'nodeIndex'}->($data->{'split'}->{'rank'})->(($i))."\n"
		unless ( any($data->{'unsplit'}->{'properties'}->{'nodeIndex'}->($data->{'unsplit'}->{'rank'}) == $data->{'split'}->{'properties'}->{'nodeIndex'}->($data->{'split'}->{'rank'})->(($i))) );
	}
	for(my $i=0;$i<nelem($data->{'unsplit'}->{'rank'});++$i) {
	    print "In unsplit but not split: ".$data->{'unsplit'}->{'properties'}->{'nodeIndex'}->($data->{'unsplit'}->{'rank'})->(($i))."\n"
		unless ( any($data->{'split'}->{'properties'}->{'nodeIndex'}->($data->{'split'}->{'rank'}) == $data->{'unsplit'}->{'properties'}->{'nodeIndex'}->($data->{'unsplit'}->{'rank'})->(($i))) );
	}
	print "FAILED {".$forestFile->{'label'}."}: number of nodes differ\n";
	$status       = 2;
	$statusGlobal = 1;
    }
    next
	unless ( $status == 0 );
    print "...done\n";

    # Test for equal properties.
    foreach my $property ( @properties ) {
	print "Testing for equal node ".$property.$orphansText."...\n";
	for(my $i=0;$i<nelem($data->{'split'}->{'rank'});++$i) {
	    my $js = $data->{  'split'}->{'rank'}->(($i));
	    my $ju = $data->{'unsplit'}->{'rank'}->(($i));
	    # Test for equality - allow some tolerance for floating point values.
	    my $areEqual;
	    if ( $data->{'split'}->{'properties'}->{$property}->type() eq "double" ) {
		$areEqual =
		    (
		         ($data->{'split'}->{'properties'}->{$property}->(($js)) -  $data->{'unsplit'}->{'properties'}->{$property}->(($ju)))
		     <
		     +0.5e+0
		     *1.0e-6
		     *abs($data->{'split'}->{'properties'}->{$property}->(($js)) +  $data->{'unsplit'}->{'properties'}->{$property}->(($ju)))
		    )
		    ||
		    (
		     $data->{'split'  }->{'properties'}->{$property}->(($js)) == 0.0
		     &&
		     $data->{'unsplit'}->{'properties'}->{$property}->(($ju)) == 0.0
		    )
	    } else {
		$areEqual =
		        ($data->{'split'}->{'properties'}->{$property}->(($js)) == $data->{'unsplit'}->{'properties'}->{$property}->(($ju)))
	    }
	    unless ( $areEqual ) {
		# Identify to which tree the problem node belongs.
		my $nodeCount = 0;
		my $iTree = 0;
		while ( $data->{'split'}->{'mergerTreeStartIndex'}->(($iTree))+$data->{'split'}->{'mergerTreeCount'}->(($iTree))-1 < $js ) {
		    ++$iTree;
		}
		# Report on the problem.
		print $i."\t".$data->{'split'}->{'properties'}->{'nodeIndex'}->(($js))."\t".$data->{'unsplit'}->{'properties'}->{'nodeIndex'}->(($ju))."\t:\t".$data->{'split'}->{'properties'}->{$property}->(($js))."\t".$data->{'unsplit'}->{'properties'}->{$property}->(($ju))."\t:\t".$data->{'split'}->{'properties'}->{'nodeIsIsolated'}->(($js))."\t".$data->{'unsplit'}->{'properties'}->{'nodeIsIsolated'}->(($ju))." : ".$data->{'split'}->{'mergerTreeIndex'}->(($iTree))."\n";
		print "FAILED {".$forestFile->{'label'}."}: '".$property."' mismatch\n";
		$status       = 3;
		$statusGlobal = 1;
	    }
	}
	print "...done\n";
    }
    next
	unless ( $status == 0 );
    # No failures found, report success.
    print "SUCCESS {".$forestFile->{'label'}."}: split forests match unsplit forests\n";
}
# Final report.
if ( $statusGlobal == 0 ) {
    print "SUCCESS: all split forest tests were successful\n";
} else {
    print "FAILED: some split forest tests failed - see preceeding reports\n";
}
exit 0;
