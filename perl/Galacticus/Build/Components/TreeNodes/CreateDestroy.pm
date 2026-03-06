# Contains a Perl module which handles creation and destruction of the treeNode class.

package Galacticus::Build::Components::TreeNodes::CreateDestroy;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;
use Data::Dumper;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodesCreateDestroy =>
     {
	 functions =>
	     [
	      \&Tree_Node_Creation         ,
	      \&Tree_Node_Builder          ,
	      \&Tree_Node_Finalization
	      ],
	 classIteratedFunctions =>
	     [
	      \&Tree_Node_Class_Creation   ,
	      \&Tree_Node_Class_Destruction
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
	     "Error"
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
{join("",map {"allocate(self%component".$_."(1))\n"} @{$build->{'componentClassListActive'}})}
select type (self)
type is (treeNode)
{join("",map {"   self%component".$_."(1)%hostNode => self\n"} @{$build->{'componentClassListActive'}})}
end select
! Assign a host tree if supplied.
if (present(hostTree)) then
   self%hostTree => hostTree
else
   self%hostTree => null()
end if
! Assign index if supplied.
if (present(index)) then
 call self%indexSet(index)
else
 call self%indexSet(-1_kind_int8)
end if
! Assign a unique ID.
!$omp critical(UniqueID_Assign)
uniqueIDCount=uniqueIDCount+1
if (uniqueIDCount <= 0) call Error_Report('ran out of unique ID numbers'//\{introspection:location\})
self%uniqueIdValue=uniqueIDCount
!$omp end critical(UniqueID_Assign)
! Assign a timestep and subsampling weight.
self%timeStepValue         =-1.0d0
self%subsamplingWeightValue= 1.0d0
! Initialize physical state.
self%isSolvable            =.true.
self%isPhysicallyPlausible =.true.
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "initialize"
	}
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
	     "FoX_DOM, only : node, nodeList, getChildNodes, getLength, getNodeName, item",
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
		 type       => "integerHash",
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
	next
	    unless ( grep {$code::component->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
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
	next
	    unless ( grep {$code::component->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
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
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "componentBuilder"
	}
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
	modules     =>
	    [
	     "Error"
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
{join(" ",map {"call self%".$_."Destroy()\n"} @{$build->{'componentClassListActive'}})}
! Destroy any formation node.
if (associated(self%formationNode)) then
   call self%formationNode%destroy()
   deallocate(self%formationNode)
   nullify(self%formationNode)
end if
! Remove any events attached to the node, along with their paired event in other nodes.
thisEvent => self%event
do while (associated(thisEvent))
   ! If a paired node is given, remove any paired event from it.
   if (associated(thisEvent%node)) then
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
      if (.not.pairMatched) call Error_Report('unable to find paired event'//\{introspection:location\})
   end if
   nextEvent => thisEvent%next
   deallocate(thisEvent)
   thisEvent => nextEvent
end do
CODE
     # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "destroy"
	}
	);
}

sub Tree_Node_Class_Creation {
    # Generate a function to create component classes.
    my $build    = shift();
    $code::class = shift();
    return
	unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
    my $function =
    {
	type        => "void",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."Create",
	description => "Create the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component of {\\normalfont \\ttfamily self}.",
	modules     =>
	    [
	     "ISO_Varying_String",
	     "Display",
	     "Error",
	     "String_Handling"
  	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "target", "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "template" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };
    my @nonNullComponents = map {$_->{'name'} eq "null" ? () : $_->{'name'}} @{$build->{'componentClasses'}->{$code::class->{'name'}}->{'members'}};
    $code::nonNullComponents = "char(10)//".join("//char(10)//",map {"'   ".$_."'"} sort(@nonNullComponents));
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
if (displayVerbosity() >= verbosityLevelInfo) then
  block
    type(varying_string) :: message
    message='Creating {$class->{'name'}} in node '
    message=message//self%index()
    call displayMessage(message,verbosityLevelInfo)
  end block
end if
if (present(template)) then
   allocate(self%component{ucfirst($class->{'name'})}(1),source=template)
else
   select type (default{ucfirst($class->{'name'})}Component)
   type is (nodeComponent{ucfirst($class->{'name'})}Null)
      block
         type(varying_string) :: message
         message=         'creation of the {$class->{'name'}} component requested, but that component is null'//char(10)
         message=message//'please select a non-null {$class->{'name'}} component - available options are:'
         message=message//{$nonNullComponents}
         call displayMessage(message,verbosityLevelSilent)
         call Error_Report('refusing to create null instance'//\{introspection:location\})
      end block
   class default
      allocate(self%component{ucfirst($class->{'name'})}(1),source=default{ucfirst($class->{'name'})}Component)
   end select
end if
select type (self)
type is (treeNode)
  do i=1,size(self%component{ucfirst($class->{'name'})})
    self%component{ucfirst($class->{'name'})}(i)%hostNode => self
    call self%component{ucfirst($class->{'name'})}(i)%initialize()
  end do
end select
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => $code::class->{'name'}."Create", 
	}
	);
}


sub Tree_Node_Class_Destruction {
    # Generate a function to destroy component classes.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "void",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."Destroy",
	description => "Destroy the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component of {\\normalfont \\ttfamily self}",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%destroy()
  end do
  deallocate (self%component{ucfirst($class->{'name'})})
end if
CODE
    } else {
	$function->{'modules'} = [ "Error" ];
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
call Error_Report("component '{$class->{'name'}}' is not active"//\{introspection:location\})
CODE
    }
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => $code::class->{'name'}."Destroy", 
	}
	);
}


1;
