# Contains a Perl module which builds the treeNode class.

package TreeNodes;
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
require Galacticus::Build::Components::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodes => 
     {
	 types     =>
	     [
	      \&Build_Tree_Node_Class
	     ],
	 functions =>
	     [
	      \&Insert_Interrupt_Interface
	     ]
     }
    );

sub Build_Tree_Node_Class {
    # Build the "treeNode" class.
    my $build = shift();
    # Define bound functions.
    my @typeBoundFunctions = 
	(
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "type"                                                                                                            ,
	     function    => "Tree_Node_Type"                                                                                                  ,
	     description => "Return the type of this node."                                                                                   ,
	     returnType  => "\\textcolor{red}{\\textless type(varying\\_string)\\textgreater}"                                                ,
	     arguments   => ""

	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "index"                                                                                                           ,
	     function    => "Tree_Node_Index"                                                                                                 ,
	     description => "Return the index of this node."                                                                                  ,
	     returnType  => "\\textcolor{red}{\\textless integer(kind\\_int8)\\textgreater}"                                                  ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "indexSet"                                                                                                        ,
	     function    => "Tree_Node_Index_Set"                                                                                             ,
	     description => "Set the index of this node."                                                                                     ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\textcolor{red}{\\textless integer(kind\\_int8)\\textgreater} index\\argin"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "timeStep"                                                                                                        ,
	     function    => "Tree_Node_Time_Step"                                                                                             ,
	     description => "Return the time-step last used by this node."                                                                    ,
	     returnType  => "\\doublezero"                                                                                                    ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "timeStepSet"                                                                                                     ,
	     function    => "Tree_Node_Time_Step_Set"                                                                                         ,
	     description => "Set the time-step used by this node."                                                                            ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\doublezero\ index\\argin"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "uniqueID"                                                                                                        ,
	     function    => "Tree_Node_Unique_ID"                                                                                             ,
	     description => "Return the unique identifier for this node."                                                                     ,
	     returnType  => "\\textcolor{red}{\\textless integer(kind\\_int8)\\textgreater}"                                                  ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "uniqueIDSet"                                                                                                     ,
	     function    => "Tree_Node_Unique_ID_Set"                                                                                         ,
	     description => "Set the unique identifier for this node."                                                                        ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\textcolor{red}{\\textless integer(kind\\_int8)\\textgreater} uniqueID\\argin"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "initialize"                                                                                                      ,
	     function    => "treeNodeInitialize"                                                                                              ,
	     description => "Initialize this node (assigns a unique identifier, creates generic components)."                                 ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\textcolor{red}{\\textless integer(kind\\_int8)\\textgreater} index\\argin"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "destroy"                                                                                                         ,
	     function    => "treeNodeDestroy"                                                                                                 ,
	     description => "Destroy this node."                                                                                              ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "componentBuilder"                                                                                                ,
	     function    => "Tree_Node_Component_Builder"                                                                                     ,
	     description => "Build components in this node given an XML description of their properties."                                     ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\textcolor{red}{\\textless *type(node)\\textgreater} nodeDefinition\\argin"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "removeFromHost"                                                                                                  ,
	     function    => "Tree_Node_Remove_From_Host"                                                                                      ,
	     description => "Remove this node from the satellite population of its host halo."                                                ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "removeFromMergee"                                                                                                ,
	     function    => "Tree_Node_Remove_From_Mergee"                                                                                    ,
	     description => "Remove this node from the list of mergees associated with its merge target."                                     ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "isPrimaryProgenitor"                                                                                             ,
	     function    => "Tree_Node_Is_Primary_Progenitor"                                                                                 ,
	     description => "Return true if this node is the primary progenitor of its descendent, false otherwise."                          ,
	     returnType  => "\\logicalzero"                                                                                                   ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     function    => "Tree_Node_Is_Primary_Progenitor_Of_Index"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     function    => "Tree_Node_Is_Primary_Progenitor_Of_Node" 
	 },
	 {
	     type        => "generic"                                                                                                         ,
	     name        => "isPrimaryProgenitorOf"                                                                                           ,
	     function    => ["Tree_Node_Is_Primary_Progenitor_Of_Index","Tree_Node_Is_Primary_Progenitor_Of_Node"]                            ,
	     description => "Return true is this node is the primary progenitor of the specified (by index or pointer) node, false otherwise.",
	     returnType  => "\\logicalzero"                                                                                                   ,
	     arguments   => "\\textcolor{red}{\\textless integer(kind\\_int8)\\textgreater} targetNodeIndex\\argin|\\textcolor{red}{\\textless *type(treeNode)\\textgreater} targetNode\\argin"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     function    => "Tree_Node_Is_Progenitor_Of_Index"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     function    => "Tree_Node_Is_Progenitor_Of_Node" 
	 },
	 {
	     type        => "generic"                                                                                                         ,
	     name        => "isProgenitorOf"                                                                                                  ,
	     function    => ["Tree_Node_Is_Progenitor_Of_Index","Tree_Node_Is_Progenitor_Of_Node"]                                            ,
	     description => "Return true is this node is a progenitor of the specified (by index or pointer) node, false otherwise."          ,
	     returnType  => "\\logicalzero"                                                                                                   ,
	     arguments   => "\\textcolor{red}{\\textless integer(kind\\_int8)\\textgreater} targetNodeIndex\\argin|\\textcolor{red}{\\textless *type(treeNode)\\textgreater} targetNode\\argin"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "isOnMainBranch"                                                                                                  ,
	     function    => "Tree_Node_Is_On_Main_Branch"                                                                                     ,
	     description => "Return true if this node is on the main branch of its tree, false otherwise."                                    ,
	     returnType  => "\\logicalzero"                                                                                                   ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "isSatellite"                                                                                                     ,
	     function    => "Tree_Node_Is_Satellite"                                                                                          ,
	     description => "Return true if this node is a satellite, false otherwise."                                                       ,
	     returnType  => "\\logicalzero"                                                                                                   ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "lastSatellite"                                                                                                   ,
	     function    => "Tree_Node_Get_Last_Satellite"                                                                                    ,
	     description => "Return a pointer to the last satellite in the list of satellites beloning to this node."                         ,
	     returnType  => "\\textcolor{red}{\\textless *type(treeNode)\\textgreater}"                                                       ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "earliestProgenitor"                                                                                              ,
	     function    => "Tree_Node_Get_Earliest_Progenitor"                                                                               ,
	     description => "Return a pointer to the earliest progenitor (along the main branch) of this node."                               ,
	     returnType  => "\\textcolor{red}{\\textless *type(treeNode)\\textgreater}"                                                       ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "mergesWith"                                                                                                      ,
	     function    => "Tree_Node_Merges_With_Node"                                                                                      ,
	     description => "Return a pointer to the node with which this node will merge."                                                   ,
	     returnType  => "\\textcolor{red}{\\textless *type(treeNode)\\textgreater}"                                                       ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "walkBranch"                                                                                                      ,
	     function    => "treeNodeWalkBranch"                                                                                              ,
	     description => "Return a pointer to the next node when performing a walk of a single branch of the tree, excluding satellites."  ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\textcolor{red}{\\textless *type(treeNode)\\textgreater} startNode\\arginout"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "walkBranchWithSatellites"                                                                                        ,
	     function    => "treeNodeWalkBranchWithSatellites"                                                                                ,
	     description => "Return a pointer to the next node when performing a walk of a single branch of the tree, including satellites."  ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\textcolor{red}{\\textless *type(treeNode)\\textgreater} startNode\\arginout"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "walkTree"                                                                                                        ,
	     function    => "treeNodeWalkTree"                                                                                                ,
	     description => "Return a pointer to the next node when performing a walk of the entire tree, excluding satellites."              ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => ""                                                       
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "walkTreeUnderConstruction"                                                                                       ,
	     function    => "treeNodeWalkTreeUnderConstruction"                                                                               ,
	     description => "Return a pointer to the next node when performing a walk of a tree under construction."                          ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\textcolor{red}{\\textless *type(treeNode)\\textgreater} nextNode\\arginout"
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "walkTreeWithSatellites"                                                                                          ,
	     function    => "treeNodeWalkTreeWithSatellites"                                                                                  ,
	     description => "Return a pointer to the next node when performing a walk of the entire tree, including satellites."              ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "destroyBranch"                                                                                                   ,
	     function    => "treeNodeDestroyBranch"                                                                                           ,
	     description => "Destroy a branch of a merger tree rooted at this node."                                                          ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "createEvent"                                                                                                     ,
	     function    => "Tree_Node_Create_Event"                                                                                          ,
	     description => "Create a {\\normalfont \\ttfamily nodeEvent} object in this node."                                               ,
	     returnType  => "\\textcolor{red}{\\textless *type(nodeEvent)\\textgreater}"                                                      ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                                       ,
	     name        => "removePairedEvent"                                                                                               ,
	     function    => "Tree_Node_Remove_Paired_Event"                                                                                   ,
	     description => "Remove a paired {\\normalfont \\ttfamily nodeEvent} from this node."                                             ,
	     returnType  => "\\void"                                                                                                          ,
	     arguments   => "\\textcolor{red}{\\textless type(nodeEvent)\\textgreater} event\\argin"
	 }
	);
    # Add data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "integer",
	     type       => "kind=kind_int8",
	     variables  => [ "indexValue", "uniqueIdValue" ]
	 },
	 {
	     intrinsic  => "double precision",
	     variables  => [ "timeStepValue" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "treeNode",
	     attributes => [ "pointer", "public" ],
	     variables  => [ "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode" ]
	 },
	 {
	     intrinsic  => "logical",
	     attributes => [ "public" ],
	     variables  => [ "isPhysicallyPlausible" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "nodeEvent",
	     attributes => [ "public", "pointer" ],
	     variables  => [ "event" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "mergerTree",
	     attributes => [ "public", "pointer" ],
	     variables  => [ "hostTree" ]
	 }
	);
    foreach ( @{$build->{'componentClassList'}} ) {
	push(
	    @dataContent,
	    {
		intrinsic  => "class",
		type       => "nodeComponent".ucfirst($_),
		attributes => [ "allocatable", "dimension(:)" ],
		variables  => [ "component".ucfirst($_) ],
		comment    => "A generic ".$_." object."
	    }
	    );
    }
    # Create the tree node class.
    $build->{'types'}->{'treeNode'} = {
	name           => "treeNode",
	comment        => "A class for \\glspl{node} in merger trees.",
	isPublic       => 1,
	boundFunctions => \@typeBoundFunctions,
	dataContent    => \@dataContent
    };
}

sub Insert_Interrupt_Interface {
    # Insert the interrupt procedure interface.
    my $build = shift();
    # Define the interface.
    $build->{'interfaces'}->{'interruptTask'} =
    {
	name      => "interruptTask"                                        ,
	comment   => "Interface for differential evolution interrupt tasks.",
	intrinsic => "void"                                                 ,
	data      =>
	    [
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "pointer", "intent(inout)" ],
		 variables  => [ "node" ]
	     }
	    ]		
    };    
}

1;
