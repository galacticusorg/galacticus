# Contains a Perl module which provides various utility functions for tree nodes.

package Galacticus::Build::Components::TreeNodes::Utils;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodeUtils =>
     {
	 functions =>
	     [
	      \&Tree_Node_Copy             ,
	      \&Tree_Node_Move             ,
	      \&Tree_Node_Mass_Distribution,
	      \&Tree_Node_Mass_Baryonic
	     ]
     }
    );

sub Tree_Node_Copy {
    # Generate a function to copy one tree node to another.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeCopyNodeTo",
	description => "Make a copy of {\\normalfont \\ttfamily self} in {\\normalfont \\ttfamily targetNode}.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "targetNode" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "skipFormationNode", "skipEvent" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "skipFormationNodeActual", "skipEventActual" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };    
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
skipFormationNodeActual=.false.
skipEventActual        =.false.
if (present(skipFormationNode)) skipFormationNodeActual=skipFormationNode
if (present(skipEvent        )) skipEventActual        =skipEvent
{join("",map {"targetNode%".$_." =  self%".$_."\n"} ( "uniqueIdValue", "indexValue", "timeStepValue", "subsamplingWeightValue" ))}
{join("",map {"targetNode%".$_." => self%".$_."\n"} ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "event", "formationNode", "hostTree" ))}
if (skipFormationNodeActual) targetNode%formationNode => null()
if (skipEventActual        ) targetNode%event         => null()
CODE
    # Iterate over all component classes
    $function->{'content'} .=
	join("",map {"targetNode%component".ucfirst($_)." = self%component".ucfirst($_)."\n"} @{$build->{'componentClassListActive'}});
    # Update target node pointers.
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
select type (targetNode)
type is (treeNode)
CODE
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   do i=1,size(self%component{ucfirst($class->{'name'})})
     targetNode%component{ucfirst($class->{'name'})}(i)%hostNode => targetNode
   end do
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end select
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "copyNodeTo", 
	    returnType  => "\\void", 
	    arguments   => "\\textcolor{red}{\\textless class(treeNode)\\textgreater} targetNode\\arginout, \\logicalzero\\ [skipFormationNode]\\argin"
	}
	);
}

sub Tree_Node_Move {
    # Generate a function to move components of one tree node to another.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeComponentsMove",
	description => "Move components from {\\normalfont \\ttfamily self} to {\\normalfont \\ttfamily targetNode}.",
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
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "targetNode" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(targetNode%component{ucfirst($class->{'name'})})) then
  do i=1,size(targetNode%component{ucfirst($class->{'name'})})
    call targetNode%component{ucfirst($class->{'name'})}(i)%destroy()
  end do
  deallocate(targetNode%component{ucfirst($class->{'name'})})
end if
if (allocated(self      %component{ucfirst($class->{'name'})})) then
   call Move_Alloc(self%component{ucfirst($class->{'name'})},targetNode%component{ucfirst($class->{'name'})})
   do i=1,size(targetNode%component{ucfirst($class->{'name'})})
     targetNode%component{ucfirst($class->{'name'})}(i)%hostNode => targetNode
   end do
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "moveComponentsTo"
	}
	);
}

sub Tree_Node_Mass_Distribution {
    # Generate a function to construct and return the mass distribution associated with a node.
    my $build = shift();
    my $function =
    {
	type        => "class(massDistributionClass), pointer => massDistribution_",
	name        => "treeNodeMassDistribution",
	description => "Construct and return the mass distribution associated with {\\normalfont \\ttfamily self}.",
	modules     =>
	    [
	     "Mass_Distributions        , only : massDistributionClass       , massDistributionComposite, massDistributionList   , massDistributionZero, kinematicsDistributionClass, kinematicsDistributionIsothermal"                             ,
	     "Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType  , enumerationWeightByType, componentTypeAll    , componentTypeDarkMatterOnly, massTypeAll                      , massTypeDark, weightByMass"
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
		 type       => "enumerationComponentTypeType",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "componentType" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "enumerationMassTypeType",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "massType" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "enumerationWeightByType",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "weightBy" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "weightIndex" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "enumerationComponentTypeType",
		 variables  => [ "componentType_", "componentType__" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "enumerationMassTypeType",
		 variables  => [ "massType_", "massType__" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "enumerationWeightByType",
		 variables  => [ "weightBy_" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "weightIndex_" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "massDistributionClass",
		 attributes => [ "pointer" ],
		 variables  => [ "massDistributionComponent" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i", "iMassDistribution", "iMassDistributionAll", "iMassDistributionConstruct", "iEmpty" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "construct", "isDarkMatterOnly" ]
	     },
	     {
		 intrinsic  => "integer",
		 type       => "kind_int8",
		 variables  => [ "uniqueID", "uniqueIDParent" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "massDistributionList",
		 attributes => [ "pointer" ],
		 variables  => [ "massDistributionList_", "massDistributionListCopy_", "next_", "nextCopy_" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "kinematicsDistributionClass",
		 attributes => [ "pointer" ],
		 variables  => [ "kinematicsDistribution_" ]
	     }
	    ]
    };
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
! Set defaults.
if (present(componentType)) then
 componentType_=componentType
else
 componentType_=componentTypeAll
end if
if (present(massType)) then
 massType_=massType
else
 massType_=massTypeAll
end if
if (present(weightBy)) then
 weightBy_=weightBy
else
 weightBy_=weightByMass
end if
if (present(weightIndex)) then
 weightIndex_=weightIndex
else
 weightIndex_=-1
end if
! Search for a match to our ID.
iMassDistribution   =0
iMassDistributionAll=0
iEmpty              =0
uniqueID            =self       %uniqueID()
if (associated(self%parent)) then
 uniqueIDParent  =self%parent%uniqueID()
else
 uniqueIDParent  =-1_kind_int8
end if
do i=1,massDistributionsCount
 if   (                                                        &
    &   massDistributions__(i)%uniqueID      == uniqueID       &
    &  .and.                                                   &
    &   massDistributions__(i)%componentType == componentType_ &
    &  .and.                                                   &
    &   massDistributions__(i)%massType      == massType_      &
    &  .and.                                                   &
    &   massDistributions__(i)%weightBy      == weightBy_      &
    &  .and.                                                   &
    &   massDistributions__(i)%weightIndex   == weightIndex_   &
    & ) then
  iMassDistribution=i
  exit
 else if (massDistributions__(i)%uniqueID < 0_kind_int8 .and. iEmpty == 0) then
  iEmpty           =i
 end if
 if   (                                                          &
    &   massDistributions__(i)%uniqueID      == uniqueID         &
    &  .and.                                                     &
    &   massDistributions__(i)%componentType == componentTypeAll &
    &  .and.                                                     &
    &   massDistributions__(i)%massType      == massTypeAll      &
    &  .and.                                                     &
    &   massDistributions__(i)%weightBy      == weightBy_        &
    &  .and.                                                     &
    &   massDistributions__(i)%weightIndex   == weightIndex_     &
    & ) then
  iMassDistributionAll=i
 end if
end do
! If we found no match, we need to create the distribution.
construct=.false.
if (iMassDistribution == 0) then
 isDarkMatterOnly=.false.
 if      (componentType_       == componentTypeDarkMatterOnly) then
  isDarkMatterOnly=.true.
  construct       =.true.
  componentType__ =componentTypeDarkMatterOnly
  massType__      =massType_
 else if (iMassDistributionAll == 0                          ) then
  construct       =.true.
  componentType__ =componentTypeAll
  massType__      =massTypeAll
 end if
 ! If no existing all/all mass distribution matched.....
 if (construct) then
  if      (iEmpty     /= 0) then  
   ! If we have an empty slot, use that.
   iMassDistributionConstruct=iEmpty
  else
   ! Simply use the next slot, unless it is occupied by a parent node massDistribution (unless we have no choice because we have run out of slots).
   do i=1,massDistributionsCount+1
    massDistributionsLast=mod(massDistributionsLast,massDistributionsCount)+1
    if (massDistributions__(massDistributionsLast)%uniqueID == uniqueIDParent) cycle
    exit
   end do
   iMassDistributionConstruct=    massDistributionsLast
   !![
   <objectDestructor name="massDistributions__(iMassDistributionConstruct)%massDistribution_"/>
   !!]
   massDistributions__(massDistributionsLast)%uniqueID=-huge(kind_int8)
  end if
  if (isDarkMatterOnly) then
   iMassDistribution   =iMassDistributionConstruct
  else
   iMassDistributionAll=iMassDistributionConstruct
  end if
  if (.not.associated(massDistributions__(iMassDistributionConstruct)%massDistribution_)) then
   massDistributionList_     => null()
   massDistributionListCopy_ => null()
   next_                     => null()
   nextCopy_                 => null()
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (allocated(self%component{ucfirst($class->{'name'})})) then
     do i=1,size(self%component{ucfirst($class->{'name'})})
       massDistributionComponent => self%component{ucfirst($class->{'name'})}(i)%massDistribution(componentType__,massType__,weightBy_,weightIndex_)
       if (associated(massDistributionComponent)) then
          if (associated(massDistributionList_)) then
            allocate(next_    %next)
            allocate(nextCopy_%next)
            next_     => next_    %next
            nextCopy_ => nextCopy_%next
          else
            allocate(massDistributionList_    )
            allocate(massDistributionListCopy_)
            next_     => massDistributionList_
            nextCopy_ => massDistributionListCopy_
          end if
          next_    %massDistribution_ => massDistributionComponent
          nextCopy_%massDistribution_ => massDistributionComponent
          next_    %next              => null()
          nextCopy_%next              => null()
       end if
     end do
  end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   allocate(massDistributionComposite :: massDistributions__(iMassDistributionConstruct)%massDistribution_)
   select type (massDistribution__ => massDistributions__(iMassDistributionConstruct)%massDistribution_)
   type is (massDistributionComposite)
     !![
     <referenceConstruct isResult="yes" object="massDistribution__" constructor="massDistributionComposite(massDistributionList_)"/>
     !!]
   end select
   massDistributions__(iMassDistributionConstruct)%uniqueID     =self%uniqueID        ()
   massDistributions__(iMassDistributionConstruct)%componentType=     componentType__
   massDistributions__(iMassDistributionConstruct)%massType     =     massType__
   massDistributions__(iMassDistributionConstruct)%weightBy     =     weightBy_
   massDistributions__(iMassDistributionConstruct)%weightIndex  =     weightIndex_
   next_ => massDistributionListCopy_
   do while (associated(next_))
    !![
    <objectDestructor name="next_%massDistribution_" nullify="no"/>
    !!]
    nextCopy_ => next_%next
    deallocate(next_)
    next_ => nextCopy_
   end do
   nullify(massDistributionList_)
  end if
 end if
 if (isDarkMatterOnly) then
  ! We already have the relevant mass distribution constructed in the iMassDistribution slot - nothing more to do here.
 else if (componentType_ == componentTypeAll .and. massType_ == massTypeAll) then
  ! The all/all mass distribution was required - we have just created it, so return it.
  iMassDistribution=iMassDistributionAll
 else
  ! Some other mass distribution was required - get it as a subset of the all/all mass distribution.
  iEmpty=0
  do i=1,massDistributionsCount
   if (massDistributions__(i)%uniqueID < 0_kind_int8 .and. iEmpty == 0) iEmpty=i
  end do
  if (iEmpty /= 0) then
   ! If we have an empty slot, use that.
   iMassDistribution=iEmpty
  else
   ! Simply use the next slot, unless it is occupied by a parent node massDistribution (unless we have no choice because we have run out of slots).
   do i=1,massDistributionsCount+1
    massDistributionsLast=mod(massDistributionsLast,massDistributionsCount)+1
    massDistributionsLast=mod(massDistributionsLast,massDistributionsCount)+1
    ! But never replace the all/all distribution.
    if   (                                                                              &
       &   massDistributions__(massDistributionsLast)%uniqueID      == uniqueID         &
       &  .and.                                                                         &
       &   massDistributions__(massDistributionsLast)%componentType == componentTypeAll &
       &  .and.                                                                         &
       &   massDistributions__(massDistributionsLast)%massType      == massTypeAll      &
       &  .and.                                                                         &
       &   massDistributions__(massDistributionsLast)%weightBy      == weightBy_        &
       &  .and.                                                                         &
       &   massDistributions__(massDistributionsLast)%weightIndex   == weightIndex_     &
       & ) cycle
    if   (                                                                              &
       &   massDistributions__(massDistributionsLast)%uniqueID      == uniqueIDParent   &
       & ) cycle
    exit
   end do
   iMassDistribution=massDistributionsLast
   !![
   <objectDestructor name="massDistributions__(iMassDistribution)%massDistribution_"/>
   !!]
   massDistributions__(massDistributionsLast)%uniqueID=-huge(kind_int8)
  end if
  massDistributions__(iMassDistribution)%massDistribution_ => massDistributions__(iMassDistributionAll)%massDistribution_%subset        (componentType_,massType_)
  massDistributions__(iMassDistribution)%uniqueID          =  self                                                       %uniqueID      (                        )
  massDistributions__(iMassDistribution)%componentType     =                                                              componentType_
  massDistributions__(iMassDistribution)%massType          =                                                              massType_
  massDistributions__(iMassDistribution)%weightBy          =                                                              weightBy_
  massDistributions__(iMassDistribution)%weightIndex       =                                                              weightIndex_
  if (.not.associated(massDistributions__(iMassDistribution)%massDistribution_)) then
   allocate(massDistributionZero :: massDistributions__(iMassDistribution)%massDistribution_)
   select type (massDistributions___ => massDistributions__(iMassDistribution)%massDistribution_)
   type is (massDistributionZero)
    !![
    <referenceConstruct object="massDistributions___" constructor="massDistributionZero(dimensionless=.false.)"/>
    !!]
    allocate(kinematicsDistributionIsothermal :: kinematicsDistribution_)
    select type (kinematicsDistribution_)
    type is (kinematicsDistributionIsothermal)
     !![
     <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionIsothermal(velocityDispersion_=0.0d0)"/>
     !!]
    end select
    call massDistributions___%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
   end select
  end if
 end if
end if
!![
<referenceAcquire target="massDistribution_" source="massDistributions__(iMassDistribution)%massDistribution_"/>
!!]
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "massDistribution"
	}
	);
}

sub Tree_Node_Mass_Baryonic {
    # Generate a function to return the total baryonic mass associated with a node.
    my $build    = shift();
    my $function =
    {
	type        => "double precision",
	name        => "treeNodeMassBaryonic",
	description => "Return the total baryonic mass associated with {\\normalfont \\ttfamily self}.",
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
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
treeNodeMassBaryonic=0.0d0
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (allocated(self%component{ucfirst($class->{'name'})})) then
     do i=1,size(self%component{ucfirst($class->{'name'})})
       treeNodeMassbaryonic=treeNodeMassBaryonic+self%component{ucfirst($class->{'name'})}(i)%massBaryonic()
     end do
  end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "massBaryonic"
	}
	);
}

1;
