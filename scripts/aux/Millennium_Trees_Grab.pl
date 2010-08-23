#!/usr/bin/env perl
use Data::Dumper;

# Grab tree data from the Millennium databases and store as an HDF5 file that Galacticus can read.
# Andrew Benson (18-Mar-2010)

# Create a hash of named arguments.
$iArg = -1;
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Specify user and password.
$sqlUser     = $arguments{"user"};
$sqlPassword = $arguments{"password"};

# Specify any selection.
$selection   = $arguments{"select"};

# Specify the output file.
if ( exists($argument{"output"}) ) {
    $outputFile = $argument{"output"};
} else {
    $outputFile = "Millennium_Trees.hdf5";
}

# Build the SQL query.
# Basic query.
$sqlQuery = "http://www.g-vo.org/MyMillennium3?action=doQuery&SQL=select node.treeId, node.haloId, node.descendantId, node.snapNum, node.redshift, node.m_tophat, node.x, node.y, node.z, node.velX, node.velY, node.velZ, node.spinX, node.spinY, node.spinZ, node.halfmassRadius from millimil..MPAHalo node, millimil..MPAHalo root where root.haloId = node.treeId";
# Append any required selection.
$sqlQuery .= " and ".$selection unless ( $selection eq "" );
# Add an order by statement.
$sqlQuery .= " order by node.treeId";

# Initalize the tree ID to a negative number.
$lastTreeID = -1;

# Initialize merger tree counter.
$mergerTree = 0;

# Specify cosmological parameters for Millennium Simulation.
$H0 = 73.0;

# Specify volume weight for Millennium trees.
$boxLength    = 500.0/($H0/100.0); # Length of the Millennium simulation cube.
$volumeWeight = 1.0/$boxLength**3;

# Specify the mass unit used in the Millennium Simulation.
$massUnit = 1.0e10/($H0/100.0);

# Specify the length unit used in the Millennium Simulation.
$lengthUnit = 1.0/($H0/100.0);

# Open a pipe to retrieve the data.
$getCommand = "wget";
$getCommand .= " --http-user="  .$sqlUser     unless ( $sqlUser     eq "" );
$getCommand .= " --http-passwd=".$sqlPassword unless ( $sqlPassword eq "" );
$getCommand .= " \"".$sqlQuery."\" -O - |";
$totalNodes = 0;
open(getHndl,$getCommand);
while ( $line = <getHndl> ) {
    if ( $line =~ m/^\d/ ) {
	@columns = split(/,/,$line);
	$treeID                 = $columns[ 0];
	$nodeIndex              = $columns[ 1];
	$redshifts[$columns[3]] = $columns[ 4];
	$nodeMass               = $columns[ 5]*$massUnit;
	$nodeX                  = $columns[ 6]*$lengthUnit;
	$nodeY                  = $columns[ 7]*$lengthUnit;
	$nodeZ                  = $columns[ 8]*$lengthUnit;
	$nodeVx                 = $columns[ 9];
	$nodeVy                 = $columns[10];
	$nodeVz                 = $columns[11];
	$nodeSpinX              = $columns[12]*$lengthUnit;
	$nodeSpinY              = $columns[13]*$lengthUnit;
	$nodeSpinZ              = $columns[14]*$lengthUnit;
	$nodeHalfMassRadius     = $columns[15]*$lengthUnit;

	%nodeHash =  (
	    nodeIndex          => $nodeIndex,
	    parentNode         => $columns[2],
	    level              => $columns[3],
	    nodeRedshift       => $columns[4],
	    nodeMass           => $nodeMass,
	    nodeX              => $nodeX,
	    nodeY              => $nodeY,
	    nodeZ              => $nodeZ,
	    nodeVx             => $nodeVz,
	    nodeVy             => $nodeVy,
	    nodeVz             => $nodeVz,
	    nodeSpinX          => $nodeSpinX,
	    nodeSpinY          => $nodeSpinY,
	    nodeSpinZ          => $nodeSpinY,
	    nodeHalfMassRadius => $nodeHalfMassRadius
	    );

	if ( $treeID != $lastTreeID ) {
	    # New tree was found:
	    # Millennium has some nodes with zero mass, which we have removed. Remove any nodes which point to those nodes.
	    $reCheck = 1;
	    while ( $reCheck == 1 ) {
		$reCheck = 0;
		foreach $index ( keys(%treeData) ) {
		    $parentIndex = $treeData{$index}->{'parentNode'};
		    unless ( exists($treeData{$parentIndex}) || $parentIndex < 0 ) {
			delete($treeData{$index});
			$reCheck = 1;
		    }
		}
	    }
	    # Process the tree to build child and sibling pointers.
	    foreach $index ( keys(%treeData) ) {
		undef(@childNodes);
		foreach $index1 ( keys(%treeData) ) {
		    if ( $treeData{$index1}->{'parentNode'} == $index ) {$childNodes[++$#childNodes] = $index1};
		}
		@sortedChildNodes = sort {$treeData{$a}->{'nodeMass'} cmp $treeData{$b}->{'nodeMass'}} @childNodes;
		if ( $#sortedChildNodes >= 0 ) {$treeData{$index}->{'childNode'} = $sortedChildNodes[0]};
		for($iChild=0;$iChild<$#sortedChildNodes;++$iChild) {
		    $index1 = $sortedChildNodes[$iChild  ];
		    $index2 = $sortedChildNodes[$iChild+1];
		    $treeData{$index1}->{'siblingNode'} = $index2;
		}
	    }
	    # Get a list of nodes sorted by level order.
	    @sortedIndices = sort {$treeData{$a}->{'level'} cmp $treeData{$b}->{'level'}} keys(%treeData);
	    $totalNodes += $#sortedIndices;
	    if ( $#sortedIndices >= 0 ) {
		# Output tree index and volume weight.
		open(dataHndl,">data.tmp");
		print dataHndl $lastTreeID."\n";
		close(dataHndl);
		open(configFile,">config.tmp");
		print configFile "PATH /mergerTrees/mergerTree".$mergerTree."/treeIndex\n";
		print configFile "INPUT-CLASS TEXTIN\n";
		print configFile "INPUT-SIZE 64\n";
		print configFile "OUTPUT-CLASS IN\n";
		print configFile "OUTPUT-SIZE 64\n";
		print configFile "RANK 1\n";
		print configFile "DIMENSION-SIZES 1\n";
		print configFile "OUTPUT-ARCHITECTURE NATIVE\n";
		close(configFile);
		system("h5import data.tmp -c config.tmp -o ".$outputFile);
		unlink("data.tmp","config.tmp");
		open(dataHndl,">data.tmp");
		print dataHndl $volumeWeight."\n";
		close(dataHndl);
		open(configFile,">config.tmp");
		print configFile "PATH /mergerTrees/mergerTree".$mergerTree."/volumeWeight\n";
		print configFile "INPUT-CLASS TEXTFP\n";
		print configFile "INPUT-SIZE 64\n";
		print configFile "OUTPUT-CLASS FP\n";
		print configFile "OUTPUT-SIZE 64\n";
		print configFile "RANK 1\n";
		print configFile "DIMENSION-SIZES 1\n";
		print configFile "OUTPUT-ARCHITECTURE NATIVE\n";
		close(configFile);
		system("h5import data.tmp -c config.tmp -o ".$outputFile);
		unlink("data.tmp","config.tmp");
		# Output the nodes.
		++$mergerTree;
		foreach $dataSet ( "nodeIndex", "childNode", "parentNode", "siblingNode", "nodeMass", "nodeRedshift",
				   "nodeX", "nodeY", "nodeZ", "nodeVx", "nodeVy", "nodeVz", "nodeSpinX", "nodeSpinY",
				   "nodeSpinZ", "nodeHalfMassRadius" ) {
		    open(dataHndl,">data.tmp");
		    foreach $index ( @sortedIndices ) {
			if ( exists(${$treeData{$index}}{$dataSet}) ) {
			    print dataHndl $treeData{$index}->{$dataSet}."\n";
			} else {
			    print dataHndl "-1\n";
			}
		    }
		    close(dataHndl);
		    open(configFile,">config.tmp");
		    print configFile "PATH /mergerTrees/mergerTree".$mergerTree."/".$dataSet."\n";
		    if ( $dataSet eq "nodeIndex" || $dataSet =~ m/Node$/ ) {
			print configFile "INPUT-CLASS TEXTIN\n";
			print configFile "INPUT-SIZE 64\n";
			print configFile "OUTPUT-CLASS IN\n";
			print configFile "OUTPUT-SIZE 64\n";
		    } else {
			print configFile "INPUT-CLASS TEXTFP\n";
			print configFile "INPUT-SIZE 64\n";
			print configFile "OUTPUT-CLASS FP\n";
			print configFile "OUTPUT-SIZE 64\n";
		    }
		    print configFile "RANK 1\n";
		    $nodeCount = $#sortedIndices+1;
		    print configFile "DIMENSION-SIZES ".$nodeCount."\n";
		    print configFile "OUTPUT-ARCHITECTURE NATIVE\n";
		    $chunkSize = 10;
		    if ( $nodeCount > $chunkSize ) {
			print configFile "CHUNKED-DIMENSION-SIZES ".$chunkSize."\n";
			print configFile "COMPRESSION-TYPE GZIP\n";
			print configFile "COMPRESSION-PARAM 1\n";
		    }
		    close(configFile);
		    system("h5import data.tmp -c config.tmp -o ".$outputFile);
		    unlink("data.tmp","config.tmp");
		}
	    }
	    undef(%treeData);
	}

	if ( $nodeMass > 0.0 ) {%{$treeData{$nodeIndex}} = %nodeHash};

	# Store the index of the current tree.
	$lastTreeID = $treeID;
    }
}
close(getHndl);
print "Total number of nodes added = ".$totalNodes."\n";

exit;
