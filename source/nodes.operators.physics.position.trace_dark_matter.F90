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
  Implements a node operator class that sets the positions of subhalos to trace the dark matter component of their host halo.
  !!}

  use :: Dark_Matter_Halo_Scales       , only : darkMatterHaloScaleClass
  use :: Satellite_Oprhan_Distributions, only : satelliteOrphanDistributionTraceDarkMatter
  !![
  <nodeOperator name="nodeOperatorPositionTraceDarkMatter">
   <description>A node operator class that sets the positions of subhalos to trace the dark matter component of their host halo.</description>
   <deepCopy>
    <functionClass variables="satelliteOrphanDistribution_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="satelliteOrphanDistribution_"/>
   </stateStorable>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorPositionTraceDarkMatter
     !!{
     A node operator class that sets the positions of subhalos to trace the dark matter component of their host halo.
     !!}
     private
     class(darkMatterHaloScaleClass                  ), pointer :: darkMatterHaloScale_         => null()
     type (satelliteOrphanDistributionTraceDarkMatter), pointer :: satelliteOrphanDistribution_ => null()
   contains
     !![
     <methods>
       <method method="assignPosition" description="Assign a position to a node such that it traces the dark matter of its host."/>
     </methods>
     !!]
     procedure :: nodeInitialize => positionTraceDarkMatterNodeInitialize
     procedure :: nodesMerge     => positionTraceDarkMatterNodesMerge
     procedure :: autoHook       => positionTraceDarkMatterAutoHook
     procedure :: assignPosition => positionTraceDarkMatterAssignPosition
  end type nodeOperatorPositionTraceDarkMatter
  
  interface nodeOperatorPositionTraceDarkMatter
     !!{
     Constructors for the \refClass{nodeOperatorPositionTraceDarkMatter} node operator class.
     !!}
     module procedure positionTraceDarkMatterConstructorParameters
     module procedure positionTraceDarkMatterConstructorInternal
  end interface nodeOperatorPositionTraceDarkMatter
  
contains
  
  function positionTraceDarkMatterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorPositionTraceDarkMatter} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorPositionTraceDarkMatter)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
     
    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodeOperatorPositionTraceDarkMatter(darkMatterHaloScale_)
    !![
    <objectDestructor name="darkMatterHaloScale_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function positionTraceDarkMatterConstructorParameters

  function positionTraceDarkMatterConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorPositionTraceDarkMatter} node operator class.
    !!}
    implicit none
    type (nodeOperatorPositionTraceDarkMatter)                        :: self
    class(darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    allocate(self%satelliteOrphanDistribution_)
    !![
    <referenceConstruct owner="self" isResult="yes" object="satelliteOrphanDistribution_" constructor="satelliteOrphanDistributionTraceDarkMatter(darkMatterHaloScale_)"/>
    !!]
    return
  end function positionTraceDarkMatterConstructorInternal

  subroutine positionTraceDarkMatterAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAtLevel, satelliteHostChangeEvent
    implicit none
    class(nodeOperatorPositionTraceDarkMatter), intent(inout) :: self

    call satelliteHostChangeEvent%attach(self,positionTraceDarkMatterSatelliteHostChange,openMPThreadBindingAtLevel,label='positionTraceDarkMatter')
    return
  end subroutine positionTraceDarkMatterAutoHook

  subroutine positionTraceDarkMatterDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorPositionTraceDarkMatter} node operator class.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAtLevel, satelliteHostChangeEvent
    implicit none
    type(nodeOperatorPositionTraceDarkMatter), intent(inout) :: self

    if (satelliteHostChangeEvent%isAttached(self,positionTraceDarkMatterSatelliteHostChange)) call satelliteHostChangeEvent%detach(self,positionTraceDarkMatterSatelliteHostChange)
    !![
    <objectDestructor name="self%darkMatterHaloScale_"        />
    <objectDestructor name="self%satelliteOrphanDistribution_"/>
    !!]
    return
  end subroutine positionTraceDarkMatterDestructor

  subroutine positionTraceDarkMatterNodeInitialize(self,node)
    !!{
    Assign the position of a subhalo during tree initialization.
    !!}
    implicit none
    class(nodeOperatorPositionTraceDarkMatter), intent(inout), target :: self
    type (treeNode                           ), intent(inout), target :: node
    
    call self%assignPosition(node,checkSatelliteStatus=.true.)
    return
  end subroutine positionTraceDarkMatterNodeInitialize

  subroutine positionTraceDarkMatterSatelliteHostChange(self,node)
    !!{
    Update the position of a subhalo when it changes host.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node
    
    select type (self)
    class is (nodeOperatorPositionTraceDarkMatter)
       call self%assignPosition(node,checkSatelliteStatus=.false.)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine positionTraceDarkMatterSatelliteHostChange
  
  subroutine positionTraceDarkMatterNodesMerge(self,node)
    !!{
    Update the position of the node after it merges with a new host.
    !!}
    implicit none
    class(nodeOperatorPositionTraceDarkMatter), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    
    call self%assignPosition(node,checkSatelliteStatus=.false.)
    return
  end subroutine positionTraceDarkMatterNodesMerge

  subroutine positionTraceDarkMatterAssignPosition(self,node,checkSatelliteStatus)
    !!{
    Assign a position to a node such that it traces the dark matter of the host halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition
    implicit none
    class  (nodeOperatorPositionTraceDarkMatter), intent(inout) :: self
    type   (treeNode                           ), intent(inout) :: node
    logical                                     , intent(in   ) :: checkSatelliteStatus
    class  (nodeComponentPosition              ), pointer       :: position
    
    position => node%position(autoCreate=.true.)
    if (.not.checkSatelliteStatus .or. node%isSatellite()) then
       call position%positionSet(self%satelliteOrphanDistribution_%position(node))
    else
       call position%positionSet([0.0d0,0.0d0,0.0d0])
    end if
    return
  end subroutine positionTraceDarkMatterAssignPosition
  
