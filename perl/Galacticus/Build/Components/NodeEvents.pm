# Contains a Perl module which builds the nodeEvent class.

package NodeEvents;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
use Text::Template 'fill_in_string';
require Galacticus::Build::Components::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     nodeEvents => 
     {
	 types     =>
	     [
	      \&Build_Node_Event_Class
	     ],
         functions =>
	     [
	      \&Node_Event_Task_Interface
	     ]
     }
    );

sub Build_Node_Event_Class {
    # Build the "nodeEvent" class.
    my $build = shift();
    # Add data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "integer",
	     type       => "kind=kind_int8",
	     attributes => [ "public" ],
	     variables  => [ "ID" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "treeNode",
	     attributes => [ "pointer", "public" ],
	     variables  => [ "node" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "public" ],
	     variables  => [ "time" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "nodeEvent",
	     attributes => [ "public", "pointer" ],
	     variables  => [ "next" ]
	 },
	 {
	     intrinsic  => "procedure",
	     type       => "nodeEventTask",
	     attributes => [ "public", "pointer" ],
	     variables  => [ "task" ]
	 }
	);
    # Create the tree node class.
    $build->{'types'}->{'nodeEvent'} = {
	name           => "nodeEvent",
	comment        => "Type for events attached to nodes.",
	isPublic       => 1,
	dataContent    => \@dataContent
    };
}

sub Node_Event_Task_Interface {
    # Define an interface for "nodeEvent" task functions.
    my $build = shift();
    # Define the interface.
    $build->{'interfaces'}->{'nodeEventTask'} =
    {
	name      => "nodeEventTask"                  ,
	comment   => "Interface for node event tasks.",
	intrinsic => "logical"                        ,
	data      =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeEvent",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "thisEvent" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "pointer", "intent(inout)" ],
		 variables  => [ "thisNode" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "deadlockStatus" ]
	     }
	    ]		
    };
}

1;
