!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  Implements a node operator class that computes the time of the most recent major merger between galaxies.
  !!}

  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass

  !![
  <nodeOperator name="nodeOperatorGalaxyMajorMergerTime">
   <description>A node operator class that computes the time of the most recent major merger between galaxies.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGalaxyMajorMergerTime
     !!{
     A node operator class that computes the time of the most recent major merger between galaxies.
     !!}
     private
     class  (mergerMassMovementsClass), pointer :: mergerMassMovements_    => null()
     integer                                    :: galaxyMajorMergerTimeID
  contains
     final     ::                   galaxyMajorMergerTimeDestructor
     procedure :: nodeInitialize => galaxyMajorMergerTimeNodeInitialize
     procedure :: autoHook       => galaxyMajorMergerTimeAutoHook
  end type nodeOperatorGalaxyMajorMergerTime
  
  interface nodeOperatorGalaxyMajorMergerTime
     !!{
     Constructors for the {\normalfont \ttfamily galaxyMajorMergerTime} node operator class.
     !!}
     module procedure galaxyMajorMergerTimeConstructorParameters
     module procedure galaxyMajorMergerTimeConstructorInternal
  end interface nodeOperatorGalaxyMajorMergerTime
  
contains

  function galaxyMajorMergerTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily galaxyMajorMergerTime} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorGalaxyMajorMergerTime)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(mergerMassMovementsClass         ), pointer       :: mergerMassMovements_
    
    !![
    <objectBuilder class="mergerMassMovements" name="mergerMassMovements_" source="parameters"/>
    !!]
    self=nodeOperatorGalaxyMajorMergerTime(mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerMassMovements_"/>
    !!]
    return
  end function galaxyMajorMergerTimeConstructorParameters

  function galaxyMajorMergerTimeConstructorInternal(mergerMassMovements_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily galaxyMajorMergerTime} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type (nodeOperatorGalaxyMajorMergerTime)                        :: self
    class(mergerMassMovementsClass         ), intent(in   ), target :: mergerMassMovements_
    !![
    <constructorAssign variables="*mergerMassMovements_"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="galaxyMajorMergerTime" id="self%galaxyMajorMergerTimeID" isCreator="yes"/>
    !!]
    return
  end function galaxyMajorMergerTimeConstructorInternal

  subroutine galaxyMajorMergerTimeAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel
    implicit none
    class(nodeOperatorGalaxyMajorMergerTime), intent(inout) :: self

    call satelliteMergerEvent%attach(self,satelliteMerger ,openMPThreadBindingAtLevel,label='galaxyMajorMergerTime')
    return
  end subroutine galaxyMajorMergerTimeAutoHook

  subroutine galaxyMajorMergerTimeDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily galaxyMajorMergerTime} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorGalaxyMajorMergerTime), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger )) call satelliteMergerEvent%detach(self,satelliteMerger )
    !![
    <objectDestructor name="self%mergerMassMovements_"/>
    !!]
    return
  end subroutine galaxyMajorMergerTimeDestructor

  subroutine galaxyMajorMergerTimeNodeInitialize(self,node)
    !!{
    Initialize galaxyMajorMergerTime level data.
    !!}
    use :: Galacticus_Nodes   , only : nodeComponentBasic
    implicit none
    class(nodeOperatorGalaxyMajorMergerTime), intent(inout)          :: self
    type (treeNode                         ), intent(inout), target  :: node
    class(nodeComponentBasic               )               , pointer :: basic

    basic => node%basic()
    call basic%metaPropertySet(self%galaxyMajorMergerTimeID,-1.0d0)
    return
  end subroutine galaxyMajorMergerTimeNodeInitialize

  subroutine satelliteMerger(self,node)
    !!{
    Record times of galaxy-galaxy major mergers.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (*                 ), intent(inout)          :: self
    type   (treeNode          ), intent(inout), target  :: node
    type   (treeNode          )               , pointer :: nodeHost
    class  (nodeComponentBasic)               , pointer :: basicHost
    integer                                             :: destinationGasSatellite, destinationGasHost       , &
         &                                                 destinationStarsHost   , destinationStarsSatellite
    logical                                             :: mergerIsMajor
 
    select type (self)
    class is (nodeOperatorGalaxyMajorMergerTime)
    call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
     if (mergerIsMajor) then
       ! Find the node to merge with.
       nodeHost  => node    %mergesWith()
       basicHost => nodeHost%basic     ()
       ! Record the merger time.
       call basicHost%metaPropertySet(self%galaxyMajorMergerTimeID,basicHost%time())
    end if
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger
