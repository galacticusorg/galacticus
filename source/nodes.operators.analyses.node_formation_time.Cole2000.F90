!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

  !!{
  Implements a node operator class that computes the formation time for each node using the definition of \cite{cole_hierarchical_2000}.
  !!}

  !![
  <nodeOperator name="nodeOperatorNodeFormationTimeCole2000">
   <description>A node operator class that computes the formation time for each node using the definition of \cite{cole_hierarchical_2000}.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorNodeFormationTimeCole2000
     !!{
     A node operator class that computes the formation time for each node using the definition of \cite{cole_hierarchical_2000}.
     !!}
     private
     integer          :: nodeFormationTimeID
     logical          :: reformationOnPromotionOnly
     double precision :: massFactorReformation
   contains
     !![
     <methods>
       <method method="reform" description="Implements a halo reformation event."/>
     </methods>
     !!]
     procedure :: nodeInitialize        => nodeFormationTimeCole2000NodeInitialize
     procedure :: nodePromote           => nodeFormationTimeCole2000NodePromote
     procedure :: differentialEvolution => nodeFormationTimeCole2000DifferentialEvolution
     procedure :: reform                => nodeFormationTimeCole2000Reform
  end type nodeOperatorNodeFormationTimeCole2000
  
  interface nodeOperatorNodeFormationTimeCole2000
     !!{
     Constructors for the \refClass{nodeOperatorNodeFormationTimeCole2000} node operator class.
     !!}
     module procedure nodeFormationTimeCole2000ConstructorParameters
     module procedure nodeFormationTimeCole2000ConstructorInternal
  end interface nodeOperatorNodeFormationTimeCole2000

  ! Submodule-scope pointer to self for use in callbacks.
  class(nodeOperatorNodeFormationTimeCole2000), pointer :: self_
  !$omp threadprivate(self_)
  
contains

  function nodeFormationTimeCole2000ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorNodeFormationTimeCole2000} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorNodeFormationTimeCole2000)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    logical                                                                :: reformationOnPromotionOnly
    double precision                                                       :: massFactorReformation

    !![
    <inputParameter>
      <name>massFactorReformation</name>
      <defaultValue>2.0d0</defaultValue>
      <description>The factor by which halo mass must have increased to trigger a new formation event.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reformationOnPromotionOnly</name>
      <defaultValue>.false.</defaultValue>
      <description>
	Specifies whether halo reformation should occur only at node promotion events, or at the precise time that the halo mass
	has increased sufficiently in mass.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorNodeFormationTimeCole2000(massFactorReformation,reformationOnPromotionOnly)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nodeFormationTimeCole2000ConstructorParameters

  function nodeFormationTimeCole2000ConstructorInternal(massFactorReformation,reformationOnPromotionOnly) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorNodeFormationTimeCole2000} node operator class.
    !!}
    implicit none
    type            (nodeOperatorNodeFormationTimeCole2000)               :: self
    logical                                               , intent(in   ) :: reformationOnPromotionOnly
    double precision                                      , intent(in   ) :: massFactorReformation
    !![
    <constructorAssign variables="massFactorReformation, reformationOnPromotionOnly"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="nodeFormationTime" id="self%nodeFormationTimeID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function nodeFormationTimeCole2000ConstructorInternal

  subroutine nodeFormationTimeCole2000NodeInitialize(self,node)
    !!{
    Initialize node formation times.
    !!}
    implicit none
    class(nodeOperatorNodeFormationTimeCole2000), intent(inout), target  :: self
    type (treeNode                             ), intent(inout), target  :: node

    ! Initialize the formation time in any leaf node.
    if (.not.associated(node%firstChild)) call self%reform(node)
    return
  end subroutine nodeFormationTimeCole2000NodeInitialize
 
  subroutine nodeFormationTimeCole2000NodePromote(self,node)
    !!{
    Check for a reformation event during a halo promotion.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeFormationTimeCole2000), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    class(nodeComponentBasic                   ), pointer       :: basicFormation, basicParent

    if (self%reformationOnPromotionOnly) return
    basicParent    => node%parent       %basic()
    basicFormation => node%formationNode%basic()
    if (basicParent%mass() > self%massFactorReformation*basicFormation%mass())  call self%reform(node)
    return
  end subroutine nodeFormationTimeCole2000NodePromote

  subroutine nodeFormationTimeCole2000DifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Check for a reformation event in a halo.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, nodeComponentBasic
    implicit none
    class    (nodeOperatorNodeFormationTimeCole2000), intent(inout), target  :: self
    type     (treeNode                             ), intent(inout), target  :: node
    logical                                         , intent(inout)          :: interrupt
    procedure(interruptTask                        ), intent(inout), pointer :: functionInterrupt
    integer                                         , intent(in   )          :: propertyType
    class    (nodeComponentBasic                   )               , pointer :: basicFormation   , basic

    ! Return immediately if inactive variables are requested or if reformation happens only on halo formation.
    if (propertyInactive(propertyType) ) return
    if (self%reformationOnPromotionOnly) return
    ! Check if the halo has grown sufficiently in mass to trigger a new formation event.
    basic          => node              %basic()
    basicFormation => node%formationNode%basic()
    if (basic%mass() > self%massFactorReformation*basicFormation%mass()) then
       interrupt         =  .true.
       self_             => self
       functionInterrupt => reformOnInterrupt
       return
    end if
    return
  end subroutine nodeFormationTimeCole2000DifferentialEvolution

  subroutine reformOnInterrupt(node,timeEnd)
    !!{
    Wrapper function to perform node reformation during interrupt of differential evolution.
    !!}
    type            (treeNode), intent(inout), target   :: node
    double precision          , intent(in   ), optional :: timeEnd
    !$GLC attributes unused :: timeEnd

    call self_%reform(node)
    return
  end subroutine reformOnInterrupt

  subroutine nodeFormationTimeCole2000Reform(self,node)
    !!{
    Creates a halo formation time component for {\normalfont \ttfamily node}. This function is also used to ``reform'' the halo, since it
    simply resets the formation time and mass to the current values.
    !!}
    use :: Events_Halo_Formation, only : Event_Halo_Formation
    use :: Galacticus_Nodes     , only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeFormationTimeCole2000), intent(inout), target  :: self
    type (treeNode                             ), intent(inout), target  :: node
    class(nodeComponentBasic                   )               , pointer :: basic
    
    ! Trigger a halo formation event.
    call Event_Halo_Formation(node)
    ! Make a copy of the formation node, and decouple it from the tree, using the parentNode pointer to point to the node of which
    ! it is the formation node.
    if (associated(node%formationNode)) then
       call node%formationNode%destroy()
       deallocate(node%formationNode)
    end if
    allocate(node%formationNode)
    call node%copyNodeTo(node%formationNode,skipFormationNode=.true.)
    node %formationNode%parent         => node
    node %formationNode%firstChild     => null()
    node %formationNode%sibling        => null()
    node %formationNode%firstSatellite => null()
    node %formationNode%firstMergee    => null()
    node %formationNode%mergeTarget    => null()
    node %formationNode%siblingMergee  => null()
    node %formationNode%formationNode  => null()
    basic                              => node%basic()
    call basic%floatRank0MetaPropertySet(self%nodeFormationTimeID,basic%time())
    return
  end subroutine nodeFormationTimeCole2000Reform
