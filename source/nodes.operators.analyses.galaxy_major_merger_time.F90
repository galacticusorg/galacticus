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
  Implements a node operator class that computes the time of the most recent major merger between galaxies.
  !!}

  use, intrinsic :: ISO_C_Binding                   , only : c_size_t
  use            :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass

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
     integer(c_size_t                )          :: countTimesMaximum
  contains
     final     ::            galaxyMajorMergerTimeDestructor
     procedure :: autoHook=> galaxyMajorMergerTimeAutoHook
  end type nodeOperatorGalaxyMajorMergerTime
  
  interface nodeOperatorGalaxyMajorMergerTime
     !!{
     Constructors for the \refClass{nodeOperatorGalaxyMajorMergerTime} node operator class.
     !!}
     module procedure galaxyMajorMergerTimeConstructorParameters
     module procedure galaxyMajorMergerTimeConstructorInternal
  end interface nodeOperatorGalaxyMajorMergerTime
  
contains

  function galaxyMajorMergerTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorGalaxyMajorMergerTime} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorGalaxyMajorMergerTime)                :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    class  (mergerMassMovementsClass         ), pointer       :: mergerMassMovements_
    integer(c_size_t                         )                :: countTimesMaximum
    
    !![
    <inputParameter>
      <name>countTimesMaximum</name>
      <defaultValue>huge(0_c_size_t)</defaultValue>
      <description>The maximum number of major merger times to accumulate for each node. Defaults to the maximum possible.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerMassMovements" name="mergerMassMovements_" source="parameters"/>
    !!]
    self=nodeOperatorGalaxyMajorMergerTime(countTimesMaximum,mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerMassMovements_"/>
    !!]
    return
  end function galaxyMajorMergerTimeConstructorParameters

  function galaxyMajorMergerTimeConstructorInternal(countTimesMaximum,mergerMassMovements_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorGalaxyMajorMergerTime} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type   (nodeOperatorGalaxyMajorMergerTime)                        :: self
    class  (mergerMassMovementsClass         ), intent(in   ), target :: mergerMassMovements_
    integer(c_size_t                         ), intent(in   )         :: countTimesMaximum
    !![
    <constructorAssign variables="countTimesMaximum, *mergerMassMovements_"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="galaxyMajorMergerTime" id="self%galaxyMajorMergerTimeID" rank="1" isCreator="yes"/>
    !!]
    return
  end function galaxyMajorMergerTimeConstructorInternal

  subroutine galaxyMajorMergerTimeAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel, dependencyDirectionAfter, dependencyRegEx
    implicit none
    class(nodeOperatorGalaxyMajorMergerTime), intent(inout) :: self
    type (dependencyRegEx                  ), dimension(1)  :: dependenciesSatelliteMerger 

    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='galaxyMajorMergerTime',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine galaxyMajorMergerTimeAutoHook
  
  subroutine galaxyMajorMergerTimeDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorGalaxyMajorMergerTime} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorGalaxyMajorMergerTime), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    !![
    <objectDestructor name="self%mergerMassMovements_"/>
    !!]
    return
  end subroutine galaxyMajorMergerTimeDestructor

  subroutine satelliteMerger(self,node)
    !!{
    Record times of galaxy-galaxy major mergers.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: Satellite_Merging_Mass_Movements, only : enumerationDestinationMergerType
    implicit none
    class           (*                               ), intent(inout)              :: self
    type            (treeNode                        ), intent(inout), target      :: node
    type            (treeNode                        )               , pointer     :: nodeHost
    class           (nodeComponentBasic              )               , pointer     :: basicHost
    double precision                                  , dimension(:) , allocatable :: majorMergerTimesCurrent, majorMergerTimesNew
    type            (enumerationDestinationMergerType)                             :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                                            destinationGasHost     , destinationStarsHost
    logical                                                                        :: mergerIsMajor
 
    select type (self)
    class is (nodeOperatorGalaxyMajorMergerTime)
    call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
     if (mergerIsMajor) then
       ! Find the node to merge with.
       nodeHost  => node    %mergesWith()
       basicHost => nodeHost%basic     ()
       ! Append the merger time.
       majorMergerTimesCurrent=basicHost%floatRank1MetaPropertyGet(self%galaxyMajorMergerTimeID)
       if (size(majorMergerTimesCurrent) < self%countTimesMaximum) then
          allocate(majorMergerTimesNew(size(majorMergerTimesCurrent)+1_c_size_t))
          majorMergerTimesNew(1_c_size_t:size(majorMergerTimesCurrent)           )=majorMergerTimesCurrent(          :                             )
       else
          allocate(majorMergerTimesNew(size(majorMergerTimesCurrent)           ))
          majorMergerTimesNew(1_c_size_t:size(majorMergerTimesCurrent)-1_c_size_t)=majorMergerTimesCurrent(2_c_size_t:size(majorMergerTimesCurrent))
       end if
       majorMergerTimesNew(size(majorMergerTimesNew))=basicHost%time()
       call basicHost%floatRank1MetaPropertySet(self%galaxyMajorMergerTimeID,majorMergerTimesNew)
    end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger
