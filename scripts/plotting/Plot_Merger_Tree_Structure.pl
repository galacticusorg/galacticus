#!/usr/bin/env perl
use POSIX;
use GraphViz;
use Math::SigFigs;

# A simple example of how to use the GraphViz package to make plots of merger trees from Galacticus output. This
# version plots a merger tree with each node labelled by the log10 of its mass and the node index (in []). Edges
# connect parents to their children and are labelled with the time (in Gyr) at which the merging event occurs.
#
# Andrew Benson (28-Apr-2010)

# Check that syntax was correct.
if ( $#ARGV != 2 ) {die("Usage: Plot_Merger_Tree_Structure.pl <galacticusFile> <treeNumber> <outputFile>")};

# Get name of Galacticus output file.
$fileName   = $ARGV[0];

# Get number of tree to plot.
$treeNumber = $ARGV[1];

# Get name of output file.
$outputFile = $ARGV[2];

# Create name of group holding the merger tree structure.
$structureGroup = $fileName."/mergerTreeStructures/mergerTree".$treeNumber."/";

# Create a GraphViz object.
my $graphVizObject = GraphViz->new(
    width    =>  8.5,
    height   => 11.0,
    directed => 1,
    node     => {fontname  =>'arial', fontsize => '120', shape => 'circle'},
    edge     => {fontname  =>'arial', fontsize => '120'}
    );

# Create a set of nodes - we need the node index, mass and time from the Galacticus output file.
open(indexPipe,"h5ls -S -d ".$structureGroup."nodeIndex |");
open( massPipe,"h5ls -S -d ".$structureGroup."nodeMass  |");
open( timePipe,"h5ls -S -d ".$structureGroup."nodeTime  |");
# Read a line from the index pipe.
$massRoot = -1.0;
$timeRoot = -1.0;
while ( $index = <indexPipe> ) {
    $index =~ s/[\s\n]//g;
    # Get lines from mass and time pipes also.
    $mass = <massPipe>;
    $mass =~ s/[\s\n]//g;
    $time = <timePipe>;
    $time =~ s/[\s\n]//g;
    # Check if this is a data line.
    if ( $index =~ m/^\s*\d+\s*$/ ) {
	# If root mass has not yet been set, set it now.
	if ( $massRoot < 0.0 ) {$massRoot = $mass};
	# If root time has not yet been set, set it now.
	if ( $timeRoot < 0.0 ) {$timeRoot = $time};
	# Generate node label.
	$label            = FormatSigFigs(POSIX::log10($mass),4)." [".$index."]";
	# Convert mass to a plottable value.
	$mass             = 10.0*(POSIX::log10($mass/$massRoot)+3.0)/3.0;
	# Convert time to a plottable value.
	$realTime[$index] = $time;
	$times   [$index] = int(1000.0*(POSIX::log10($time/$timeRoot)+3.0)/3.0);
	# Create a node.
	$graphVizObject->add_node($index, height => $mass, label => $label);
    }
}
close(indexPipe);
close( massPipe);
close( timePipe);

# Create links to parent nodes.
open(indexPipe ,"h5ls -S -d ".$structureGroup."nodeIndex       |");
open(parentPipe,"h5ls -S -d ".$structureGroup."parentNodeIndex |");
while ( $index = <indexPipe> ) {
    $index =~ s/[\s\n]//g;
    $parentIndex = <parentPipe>;
    $parentIndex =~ s/[\s\n]//g;
    @columns = split(/\s+/,$line);
    # Check if this is a data line.
    if ( $index =~ m/^\s*\d+\s*$/ ) {
	if ( $parentIndex > 0 ) {
	    $dist = $times[$parentIndex]-$times[$index];
	    $graphVizObject->add_edge($parentIndex => $index, minlen => $dist, label => FormatSigFigs($realTime[$index],4));
	}
    }
}
close(indexPipe);
close(parentPipe);

open(outfile,">treeTemp.ps");
print outfile $graphVizObject->as_ps;
close(outfile);
system("ps2pdf treeTemp.ps ".$outputFile);
unlink("treeTemp.ps");

exit;
