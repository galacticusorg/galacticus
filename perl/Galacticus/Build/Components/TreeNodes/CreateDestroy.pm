# Contains a Perl module which handles creation and destruction of the treeNode class.

package CreateDestroy;
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
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;
require Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodesCreateDestroy =>
     {
	 functions =>
	     [
	      \&Tree_Node_Creation    ,
	      \&Tree_Node_Finalization
	     ]
     }
    );

sub Tree_Node_Creation {
    # Generate a function to create tree nodes.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeInitialize",
	description => "Initialize a {\\normalfont \\ttfamily treeNode} object.",
	modules     =>
	    [
	     "Galacticus_Error"
	     ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 type       => "kind=kind_int8",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "index" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "mergerTree",
		 attributes => [ "intent(in   )", "optional", "target" ],
		 variables  => [ "hostTree" ]
	     }
	    ]
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
! Ensure pointers are nullified.
{join("",map {"nullify (self%".$_.")\n"} ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode", "event" ))}
{join("",map {"allocate(self%component".$_->{'name'}."(1))\n"} &ExtraUtils::hashList($build->{'componentClasses'}))}
select type (self)
type is (treeNode)
{join("",map {"   self%component".$_->{'name'}."(1)%hostNode => self\n"} &ExtraUtils::hashList($build->{'componentClasses'}))}
end select
! Assign a host tree if supplied.
if (present(hostTree)) self%hostTree => hostTree
! Assign index if supplied.
if (present(index)) call self%indexSet(index)
! Assign a unique ID.
!$omp critical(UniqueID_Assign)
uniqueIDCount=uniqueIDCount+1
if (uniqueIDCount <= 0) call Galacticus_Error_Report('treeNodeInitialize','ran out of unique ID numbers')
self%uniqueIdValue=uniqueIDCount
!$omp end critical(UniqueID_Assign)
! Assign a timestep.
self%timeStepValue=-1.0d0
CODE
    # Add the function to the functions list.
    push(
	@{$build->{'functions'}},
	$function
	);
}

sub Tree_Node_Finalization {
    # Generate a function to finalize tree nodes.
    my $build = shift();
    # Generate the function.
    my $function =
    {
	type        => "void",
	name        => "treeNodeDestroy",
	description => "Destroy a {\\normalfont \\ttfamily treeNode} object.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeEvent",
		 attributes => [ "pointer" ],
		 variables  => [ "thisEvent", "pairEvent", "lastEvent", "nextEvent" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "pairMatched" ]
	     }
	    ],
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
! Destroy all components.
{join(" ",map {"call self%".$_->{'name'}."Destroy()\n"} &ExtraUtils::hashList($build->{'componentClasses'}))}
! Remove any events attached to the node, along with their paired event in other nodes.
thisEvent => self%event
do while (associated(thisEvent))
    ! Locate the paired event and remove it.
    pairEvent => thisEvent%node%event
    lastEvent => thisEvent%node%event
    ! Iterate over all events.
    pairMatched=.false.
    do while (associated(pairEvent).and..not.pairMatched)
       ! Match the paired event ID with the current event ID.
       if (pairEvent%ID == thisEvent%ID) then
          pairMatched=.true.
          if (associated(pairEvent,thisEvent%node%event)) then
             thisEvent%node  %event => pairEvent%next
             lastEvent       => thisEvent%node %event
          else
             lastEvent%next  => pairEvent%next
          end if
          nextEvent => pairEvent%next
          deallocate(pairEvent)
          pairEvent => nextEvent
       else
          lastEvent => pairEvent
          pairEvent => pairEvent%next
       end if
    end do
    if (.not.pairMatched) call Galacticus_Error_Report('treeNodeDestroy','unable to find paired event')
    nextEvent => thisEvent%next
    deallocate(thisEvent)
    thisEvent => nextEvent
end do
CODE
    # Add the function to the functions list.
    push(
	@{$build->{'functions'}},
	$function
	);
}

1;
