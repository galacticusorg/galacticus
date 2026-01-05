!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implements a node operator class that records the times of gas-mass-based major mergers between galaxies.
  !!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t

  !![
  <nodeOperator name="nodeOperatorGalaxyGasMajorMergerTime">
   <description>A node operator class that records the times of gas-mass-based major mergers between galaxies.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGalaxyGasMajorMergerTime
     !!{
     A node operator class that records the times of gas-mass-based major mergers between galaxies.
     !!}
     private
     integer                    :: galaxyGasMajorMergerTimeID
     integer         (c_size_t) :: countTimesMaximum
     double precision           :: ratioGasMajorMerger
  contains
     final     ::             galaxyGasMajorMergerTimeDestructor
     procedure :: autoHook => galaxyGasMajorMergerTimeAutoHook
  end type nodeOperatorGalaxyGasMajorMergerTime
  
  interface nodeOperatorGalaxyGasMajorMergerTime
     !!{
     Constructors for the \refClass{nodeOperatorGalaxyGasMajorMergerTime} node operator class.
     !!}
     module procedure galaxyGasMajorMergerTimeConstructorParameters
     module procedure galaxyGasMajorMergerTimeConstructorInternal
  end interface nodeOperatorGalaxyGasMajorMergerTime
  
contains

  function galaxyGasMajorMergerTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorGalaxyGasMajorMergerTime} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorGalaxyGasMajorMergerTime)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    integer         (c_size_t                            )                :: countTimesMaximum
    double precision                                                      :: ratioGasMajorMerger
    
    !![
    <inputParameter>
      <name>countTimesMaximum</name>
      <defaultValue>huge(0_c_size_t)</defaultValue>
      <description>The maximum number of major merger times to accumulate for each node. Defaults to the maximum possible.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ratioGasMajorMerger</name>
      <defaultValue>0.25d0</defaultValue>
      <description>The gas mass ratio threshold defining a major merger.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorGalaxyGasMajorMergerTime(countTimesMaximum,ratioGasMajorMerger)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyGasMajorMergerTimeConstructorParameters

  function galaxyGasMajorMergerTimeConstructorInternal(countTimesMaximum,ratioGasMajorMerger) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorGalaxyGasMajorMergerTime} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type            (nodeOperatorGalaxyGasMajorMergerTime)                        :: self
    integer         (c_size_t                            ), intent(in   )         :: countTimesMaximum
    double precision                                      , intent(in   )         :: ratioGasMajorMerger
    !![
    <constructorAssign variables="countTimesMaximum, ratioGasMajorMerger"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="galaxyGasMajorMergerTime" id="self%galaxyGasMajorMergerTimeID" rank="1" isCreator="yes"/>
    !!]
    return
  end function galaxyGasMajorMergerTimeConstructorInternal

  subroutine galaxyGasMajorMergerTimeAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel, dependencyDirectionAfter, dependencyRegEx
    implicit none
    class(nodeOperatorGalaxyGasMajorMergerTime), intent(inout) :: self
    type (dependencyRegEx                     ), dimension(1)  :: dependenciesSatelliteMerger 

    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='preAnalysis:galaxyGasMajorMergerTime',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine galaxyGasMajorMergerTimeAutoHook
  
  subroutine galaxyGasMajorMergerTimeDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorGalaxyGasMajorMergerTime} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorGalaxyGasMajorMergerTime), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    return
  end subroutine galaxyGasMajorMergerTimeDestructor

  subroutine satelliteMerger(self,node)
    !!{
    Record times of galaxy-galaxy major mergers.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentBasic              , nodeComponentDisk, nodeComponentSpheroid
    use :: Satellite_Merging_Mass_Movements, only : enumerationDestinationMergerType
    implicit none
    class           (*                    ), intent(inout)              :: self
    type            (treeNode             ), intent(inout), target      :: node
    type            (treeNode             )               , pointer     :: nodeHost
    class           (nodeComponentBasic   )               , pointer     :: basicHost
    class           (nodeComponentDisk    )               , pointer     :: disk                   , diskHost
    class           (nodeComponentSpheroid)               , pointer     :: spheroid               , spheroidHost
    double precision                       , dimension(:) , allocatable :: majorMergerTimesCurrent, majorMergerTimesNew
    logical                                                             :: mergerIsMajor
    double precision                                                    :: massGas                , massGasHost
    
    select type (self)
    class is (nodeOperatorGalaxyGasMajorMergerTime)
       nodeHost      =>  node    %mergesWith()
       disk          =>  node    %disk      ()
       spheroid      =>  node    %spheroid  ()
       diskHost      =>  nodeHost%disk      ()
       spheroidHost  =>  nodeHost%spheroid  ()
       massGas       =  +disk    %massGas()+spheroid    %massGas()
       massGasHost   =  +diskHost%massGas()+spheroidHost%massGas()
       mergerIsMajor =   massGas > self%ratioGasMajorMerger*massGasHost
       if (mergerIsMajor) then
       ! Find the node to merge with.
       basicHost => nodeHost%basic()
       ! Append the merger time.
       majorMergerTimesCurrent=basicHost%floatRank1MetaPropertyGet(self%galaxyGasMajorMergerTimeID)
       if (size(majorMergerTimesCurrent) < self%countTimesMaximum) then
          allocate(majorMergerTimesNew(size(majorMergerTimesCurrent)+1_c_size_t))
          majorMergerTimesNew(1_c_size_t:size(majorMergerTimesCurrent)           )=majorMergerTimesCurrent(          :                             )
       else
          allocate(majorMergerTimesNew(size(majorMergerTimesCurrent)           ))
          majorMergerTimesNew(1_c_size_t:size(majorMergerTimesCurrent)-1_c_size_t)=majorMergerTimesCurrent(2_c_size_t:size(majorMergerTimesCurrent))
       end if
       majorMergerTimesNew(size(majorMergerTimesNew))=basicHost%time()
       call basicHost%floatRank1MetaPropertySet(self%galaxyGasMajorMergerTimeID,majorMergerTimesNew)
    end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger
