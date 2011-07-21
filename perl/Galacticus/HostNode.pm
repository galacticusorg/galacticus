# Contains a Perl module which implements host halo property calculations for Galacticus.

package HostNode;
use PDL;
use Galacticus::HDF5;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "hostNodeMass" => "HostNode::Get_Host_Node_Mass"
    );

my $status = 1;
$status;

sub Get_Host_Node_Mass {
    $dataSet     = shift;
    $dataSetName = $_[0];
    &HDF5::Get_Dataset($dataSet,[
			   "nodeMass",
			   "nodeIndex",
			   "parentNode",
			   "nodeIsIsolated"
		       ]
	);
    $dataSets = $dataSet->{'dataSets'};

    # Create a copy of the node mass data.
    $hostNodeMass = $dataSets->{"nodeMass"}->copy();

    # Identify isolated nodes.
    $isolatedNodes = which($dataSets->{"nodeIsIsolated"} == 1);

    # Loop over all isolated nodes.
    for($i=0;$i<nelem($isolatedNodes);++$i) {
	# Find satellite nodes in the current isolated node.
	$satelliteNodes = which(
	    $dataSets->{"nodeIsIsolated"} == 0 &
	    $dataSets->{"parentNode"} == $dataSets->{"nodeIndex"}->index($isolatedNodes->index($i))
	    );
	# Set the host node mass of these satellites to the node mass of their host.
	$hostNodeMass->index($satelliteNodes) .= $dataSets->{"nodeMass"}->index($isolatedNodes->index($i));
    }

    # Transfer to the output data structure.
    $dataSets->{"hostNodeMass"} = $hostNodeMass;
}
