#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/../perl';
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use Galacticus::Options;

# Run a set of merger trees with forests split and not split. Compare the results which should be identical.
# Andrew Benson (19-May-2016)

# Global status for all tests.
my $statusGlobal = 0;
# Define a set of merger tree forest files that we want to process.
my @forestFiles =
    (
     # milli-Millennium forest number 12000000 - a modestly sized forest which can be split and has sufficiently complex structure
     # that it has uncovered bugs in split forest handling in the past.
     {
     	 label      => "milliMillennium-forest12000000-noMerging",
     	 fileName   => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree12000000.hdf5",
     	 parameters => 
     	  {
     	      mergerTreeReadSubresolutionMergingMethod => {value => "infinite"}
     	  }
     },
     {
     	 label      => "milliMillennium-forest12000000-merging",
     	 fileName   => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree12000000.hdf5",
     	 parameters => 
     	 {
     	     mergerTreeReadSubresolutionMergingMethod => {value => "boylanKolchin2008"},
     	     virialOrbitMethod                        => {value => "fixed"            },
     	     mergerTreeReadPresetMergerTimes          => {value => "true"             },
     	     mergerTreeReadPresetMergerNodes          => {value => "false"            }
     	 }
     },
     {
     	 label      => "milliMillennium-forest12000000-mergingPresetTargets",
     	 fileName   => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree12000000.hdf5",
     	 parameters => 
     	 {
     	     mergerTreeReadSubresolutionMergingMethod => {value => "boylanKolchin2008"},
     	     virialOrbitMethod                        => {value => "fixed"            },
     	     mergerTreeReadPresetMergerTimes          => {value => "true"             },
     	     mergerTreeReadPresetMergerNodes          => {value => "true"             }
     	 }
     },
     # A milli-Millennium tree which fails unless we (correctly) nullify the merge target pointer of a node which never merges but
     # which originated in a tree in the forest where the original branch it should merge with terminated due to an inter-tree
     # transfer.
     {
     	 label      => "milliMillennium-forest3000055000000-mergingPresetTargets",
     	 fileName   => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree3000055000000.hdf5",
     	 parameters => 
     	 {
     	     mergerTreeReadSubresolutionMergingMethod => {value => "boylanKolchin2008"},
     	     virialOrbitMethod                        => {value => "fixed"            },
     	     mergerTreeReadPresetMergerTimes          => {value => "true"             },
     	     mergerTreeReadPresetMergerNodes          => {value => "true"             }
     	 }
     },
     # A milli-Millennium tree which fails if the "final time in tree" used in limiting node evolution is computed across all
     # trees in a forest rather than just in the local tree.
     {
     	 label      => "milliMillennium-forest2000024000000-noMerging",
     	 fileName   => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree2000024000000.hdf5",
     	 parameters => 
     	  {
     	      mergerTreeReadSubresolutionMergingMethod => {value => "infinite"},
     	      mergerTreeReadPresetMergerNodes          => {value => "true"    }
     	  }
     },
     # A milli-Millennium tree in which a mergee is transferred as part of an inter-tree event but the destination node is a
     # clone. Requires that promotion of a cloned node to a satellite parent correctly moves any satellites in the clone into the
     # host halo of the parent.
     {
     	 label      => "milliMillennium-forest1000019000000-noMerging",
     	 fileName   => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree1000019000000.hdf5",
     	 parameters => 
     	  {
     	      mergerTreeReadSubresolutionMergingMethod => {value => "infinite"},
     	      mergerTreeReadPresetMergerNodes          => {value => "true"    }
     	  }
     },
     # A milli-Millennium tree in which a satellite is orphanized with a merge target that is a satellite. Fails unless we
     # explicitly check for that merge target being a satellite and move the new host to the merge target's parent.
     {
	 label      => "milliMillennium-forest1000020000000-noMerging",
	 fileName   => "testSuite/data/mergerTrees/splitForests-milliMillennium-tree1000020000000.hdf5",
	 parameters => 
	  {
	      mergerTreeReadSubresolutionMergingMethod => {value => "infinite"},
	      mergerTreeReadPresetMergerNodes          => {value => "true"    }
	  }
     }
    );

# Find the location of the Millennium Database data.
my $millenniumDatabaseConfig = &Galacticus::Options::Config('millenniumDB');
if ( defined($millenniumDatabaseConfig) ) {
    if ( exists($millenniumDatabaseConfig->{'path'}) ) {
	print $millenniumDatabaseConfig->{'path'}."\n";
	unless ( -e $millenniumDatabaseConfig->{'path'}."/milliMillennium/milliMillennium.hdf5" ) {
	    print "SKIPPED: milli-Millennium data not available\n";
	    exit;
	}	    
    } else {
	print "SKIPPED: Millennium database data not available\n";
	exit;
    }
} else {
    print "SKIPPED: Millennium database location undefined\n";
    exit;
}
## AJB HACK
# push(
#     @forestFiles,
#     {
# 	label    => "milliMillennium",
# 	fileName => $millenniumDatabaseConfig->{'path'}."/milliMillennium/milliMillennium.hdf5"
#     );
    
# Define the set of properties that we want to compare.
my @properties = ( "nodeIndex", "positionPositionX", "positionPositionY", "positionPositionZ", "positionVelocityX", "positionVelocityY", "positionVelocityZ", "satelliteNodeIndex", "satelliteMergeTime", "satelliteBoundMass", "nodeIsIsolated", "basicTimeLastIsolated", "parentIndex", "basicMass" );

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
	    $parameters->{'galacticusOutputFileName'}->{'value'} = "testSuite/".$outputDirectoryName.$_.".hdf5";
	    $parameters->{'mergerTreeReadFileName'  }->{'value'} = $forestFile->{'fileName'};
	    $parameters->{'treeEvolveSuspendPath'   }->{'value'} = defined($scratchConfig) ? $scratchConfig->{'path'} : ".";
	    $parameters->{$_} = $forestFile->{'parameters'}->{$_}
	       foreach ( keys(%{$forestFile->{'parameters'}}) );
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
	    $data->{$_}->{'rank'} = $data->{$_}->{'properties'}->{'nodeIndex'}->qsorti();
	    # Get tree indices.
	    foreach my $treeDatasetName ( "mergerTreeIndex", "mergerTreeStartIndex", "mergerTreeCount" ) {
		$data->{$_}->{$treeDatasetName} = $data->{$_}->{'file'}->dataset("Outputs/Output1/".$treeDatasetName)->get();
	    }
	}
    }
    next
	unless ( $status == 0 );
    # Test for equal numbers of nodes.
    print "Testing for equal numbers of nodes...\n";
    unless ( nelem($data->{'split'}->{'rank'}) == nelem($data->{'unsplit'}->{'rank'}) ) {
	for(my $i=0;$i<nelem($data->{'split'}->{'properties'}->{'nodeIndex'});++$i) {
	    print "In split but not unsplit: ".$data->{'split'}->{'properties'}->{'nodeIndex'}->(($i))."\n"
		unless ( any($data->{'unsplit'}->{'properties'}->{'nodeIndex'} == $data->{'split'}->{'properties'}->{'nodeIndex'}->(($i))) );
	}
	for(my $i=0;$i<nelem($data->{'unsplit'}->{'properties'}->{'nodeIndex'});++$i) {
	    print "In unsplit but not split: ".$data->{'unsplit'}->{'properties'}->{'nodeIndex'}->(($i))."\n"
		unless ( any($data->{'split'}->{'properties'}->{'nodeIndex'} == $data->{'unsplit'}->{'properties'}->{'nodeIndex'}->(($i))) );
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
	print "Testing for equal node ".$property."...\n";
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
		while ( $data->{'split'}->{'mergerTreeStartIndex'}->(($iTree))+$data->{'split'}->{'mergerTreeCount'}->(($iTree)) < $js ) {
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
