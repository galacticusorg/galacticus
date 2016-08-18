# Contains a Perl module which handles creation and destruction of the treeNode class.

package Galacticus::Build::Components::TreeNodes::CreateDestroy;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodesCreateDestroy =>
     {
	 functions =>
	     [
	      \&Tree_Node_Creation    ,
	      \&Tree_Node_Builder     ,
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
{join("",map {"allocate(self%component".$_->{'name'}."(1))\n"} &List::ExtraUtils::hashList($build->{'componentClasses'}))}
select type (self)
type is (treeNode)
{join("",map {"   self%component".$_->{'name'}."(1)%hostNode => self\n"} &List::ExtraUtils::hashList($build->{'componentClasses'}))}
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

sub Tree_Node_Builder {
    # Generate a function to build tree nodes from an XML definition.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeComponentBuilder",
	description => "Build components in a {\\normalfont \\ttfamily treeNode} object given an XML definition.",
	modules     =>
	    [
	     "FoX_DOM",
	     "Hashes"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "node",
		 attributes => [ "intent(in   )", "pointer" ],
		 variables  => [ "nodeDefinition" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "node",
		 attributes => [ "pointer" ],
		 variables  => [ "componentDefinition" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "nodeList",
		 attributes => [ "pointer" ],
		 variables  => [ "componentList" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i", "j", "componentCount" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "integerScalarHash",
		 variables  => [ "componentIndex" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=128",
		 variables  => [ "nodeName" ]
	     }
	    ]
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
select type (self)
type is (treeNode)
   call componentIndex%initialize()
   !$omp critical (FoX_DOM_Access)
CODE
    foreach $code::component ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    componentList => getChildNodes(nodeDefinition)
    componentCount=0
    do i=0,getLength(componentList)-1
      componentDefinition => item(componentList,i)
      if (getNodeName(componentDefinition) == '{$component->{'name'}}') componentCount=componentCount+1
    end do
    if (componentCount > 0) then
      if (allocated(self%component{ucfirst($component->{'name'})})) deallocate(self%component{ucfirst($component->{'name'})})
      allocate(self%component{ucfirst($component->{'name'})}(componentCount),source=default{ucfirst($component->{'name'})}Component)
      call componentIndex%set('{$component->{'name'}}',0)
    end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   componentCount=getLength(componentList)
   !$omp end critical (FoX_DOM_Access)
   do i=0,componentCount-1
CODE
    foreach $code::component ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
     !$omp critical (FoX_DOM_Access)
     componentDefinition => item(componentList,i)
     nodeName=getNodeName(componentDefinition)
     !$omp end critical (FoX_DOM_Access)
     if (trim(nodeName) == '{$component->{'name'}}') then
       j=componentIndex%value('{$component->{'name'}}')
       j=j+1
       self%component{ucfirst($component->{'name'})}(j)%hostNode => self
       call self%component{ucfirst($component->{'name'})}(j)%builder(componentDefinition)
       call componentIndex%set('{$component->{'name'}}',j)
     end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  end do
  call componentIndex%destroy()
end select
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
{join(" ",map {"call self%".$_->{'name'}."Destroy()\n"} &List::ExtraUtils::hashList($build->{'componentClasses'}))}
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
