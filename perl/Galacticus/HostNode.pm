# Contains a Perl module which implements host halo property calculations for Galacticus.

package HostNode;
use PDL;
require Galacticus::HDF5;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "hostNodeMass" => \&HostNode::Get_Host_Node_Mass
    );

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

    # Build hash of isolated node identifiers.
    my %hostMass;
    for($i=0;$i<nelem($isolatedNodes);++$i) {
	my $key = $dataSets->{"mergerTreeIndex"}->index($isolatedNodes->index($i)).":".$dataSets->{"nodeIndex"}->index($isolatedNodes->index($i));
	$hostMass{$key} = $dataSets->{"nodeMass"}->index($isolatedNodes->index($i));
    }

    # Loop over all nodes.
    for($i=0;$i<nelem($hostNodeMass);++$i) {
	unless ( $dataSets->{"nodeIsIsolated"}->index($i) == 1 ) {
	    my $key = $dataSets->{"mergerTreeIndex"}->index($i).":".$dataSets->{"parentNode"}->index($i);
	    die("Galacticus::HostNode - host node not found")
		unless ( exists($hostMass{$key}) );
	    $hostNodeMass->index($i) .= $hostMass{$key};
	}
    }

    # Transfer to the output data structure.
    $dataSet->{'dataSets'}->{"hostNodeMass"} = $hostNodeMass;
}

1;
