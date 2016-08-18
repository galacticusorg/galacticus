# Contains a Perl module which implements host halo property calculations for Galacticus.

package Galacticus::HostNode;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use Galacticus::HDF5;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
    "^hostNodeMass[\\d\\.]*\$" => \&Galacticus::HostNode::Get_Host_Node_Mass
    );

sub Get_Host_Node_Mass {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Determine which mass is required.
    my$nodeMassDataSet;
    if ( $dataSetName =~ m/^hostNodeMass([\d\.]*)$/ ) {
	$nodeMassDataSet = "nodeMass".$1;
	# For the usual, total virial mass, switch the dataset to "basicMass".
	$nodeMassDataSet =~ s/^nodeMass$/basicMass/;
    } else {
	die("Get_Host_Node_Mass: dataset name is not recognized");
    }

    &Galacticus::HDF5::Get_Dataset($model,[
			   $nodeMassDataSet,
			   "nodeIndex",
			   "parentIndex",
			   "nodeIsIsolated"
		       ]
	);
    my $dataSets = $model->{'dataSets'};

    # Create a copy of the node mass data.
    my $hostNodeMass = $dataSets->{$nodeMassDataSet}->copy();
    
    # Find parent node masses taking advantage of the fact that trees are stored depth first.
    my $indexStart = 0;
    for(my $i=0;$i<nelem($dataSets->{$nodeMassDataSet});++$i) {
    	if ($dataSets->{"nodeIsIsolated"}->(($i)) == 1) {
    	    $hostNodeMass->($indexStart:$i) .= $dataSets->{$nodeMassDataSet}->(($i));
    	    $indexStart = $i+1;
    	}
    }
    
    # Transfer to the output data structure.
    $model->{'dataSets'}->{$dataSetName} = $hostNodeMass;
}

1;
