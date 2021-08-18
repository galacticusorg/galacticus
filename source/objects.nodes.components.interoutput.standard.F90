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
Contains a module which implements the standard indices component.
!!}

module Node_Component_Inter_Output_Standard
  !!{
  Implements the standard indices component.
  !!}
  use :: Output_Times                    , only : outputTimesClass
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass
  use :: Star_Formation_Rates_Disks      , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids  , only : starFormationRateSpheroidsClass
  implicit none
  private
  public :: Node_Component_Inter_Output_Standard_Rate_Compute      , Node_Component_Inter_Output_Standard_Reset    , &
       &    Node_Component_Interoutput_Standard_Thread_Uninitialize, Node_Component_Inter_Output_Standard_Scale_Set, &
       &    Node_Component_Interoutput_Standard_Thread_Initialize

  !![
  <component>
   <class>interOutput</class>
   <name>standard</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>diskStarFormationRate</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <classDefault>0.0d0</classDefault>
      <output unitsInSI="massSolar/gigaYear" comment="Disk star formation rate averaged over time between current and previous outputs."/>
    </property>
    <property>
      <name>spheroidStarFormationRate</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar/gigaYear" comment="Spheroid star formation rate averaged over time between current and previous outputs."/>
      <classDefault>0.0d0</classDefault>
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(outputTimesClass               ), pointer :: outputTimes_
  class(mergerMassMovementsClass       ), pointer :: mergerMassMovements_
  class(starFormationRateDisksClass    ), pointer :: starFormationRateDisks_
  class(starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_
  !$omp threadprivate(outputTimes_,mergerMassMovements_,starFormationRateDisks_,starFormationRateSpheroids_)

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Interoutput_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Interoutput_Standard_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node standard interoutput module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent       , openMPThreadBindingAtLevel, dependencyRegEx, dependencyDirectionAfter
    use :: Galacticus_Nodes, only : defaultInteroutputComponent
    use :: Input_Parameters, only : inputParameter             , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies

    if (defaultInteroutputComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent%attach(defaultInteroutputComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentInteroutputStandard',dependencies=dependencies)
       !![
       <objectBuilder class="outputTimes"                name="outputTimes_"                source="parameters_"/>
       <objectBuilder class="mergerMassMovements"        name="mergerMassMovements_"        source="parameters_"/>
       <objectBuilder class="starFormationRateDisks"     name="starFormationRateDisks_"     source="parameters_"/>
       <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Interoutput_Standard_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Interoutput_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Interoutput_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the tree node standard interoutput module.
    !!}
    use :: Events_Hooks    , only : satelliteMergerEvent
    use :: Galacticus_Nodes, only : defaultInteroutputComponent
    implicit none

    if (defaultInteroutputComponent%standardIsActive()) then
       call satelliteMergerEvent%detach(defaultInteroutputComponent,satelliteMerger)
       !![
       <objectDestructor name="outputTimes_"               />
       <objectDestructor name="mergerMassMovements_"       />
       <objectDestructor name="starFormationRateDisks_"    />
       <objectDestructor name="starFormationRateSpheroids_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Interoutput_Standard_Thread_Uninitialize

  !![
  <scaleSetTask>
   <unitName>Node_Component_Inter_Output_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Inter_Output_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentInterOutput   , nodeComponentInterOutputStandard, nodeComponentSpheroid, &
          &                         treeNode         , defaultInterOutputComponent
    implicit none
    type            (treeNode                ), intent(inout), pointer :: node
    class           (nodeComponentInterOutput)               , pointer :: interOutput
    class           (nodeComponentDisk       )               , pointer :: disk
    class           (nodeComponentSpheroid   )               , pointer :: spheroid
    double precision                          , parameter              :: massMinimum    =1.0d0
    double precision                          , parameter              :: timeScale      =1.0d0
    double precision                                                   :: mass

    ! Check if we are the default method.
    if (.not.defaultInterOutputComponent%standardIsActive()) return
    ! Get the interoutput component.
    interOutput => node%interOutput()
    ! Check if component is of standard class.
    select type (interOutput)
    class is (nodeComponentInterOutputStandard)
       ! Get disk and spheroid components.
       disk     => node%disk    ()
       spheroid => node%spheroid()
       ! Set scale for masses.
       mass   = disk%massGas    ()+spheroid%massGas    () &
            &  +disk%massStellar()+spheroid%massStellar()
       call interOutput%    diskStarFormationRateScale(max(mass,massMinimum)/timeScale)
       call interOutput%spheroidStarFormationRateScale(max(mass,massMinimum)/timeScale)
    end select
    return
  end subroutine Node_Component_Inter_Output_Standard_Scale_Set

  !![
  <rateComputeTask>
   <unitName>Node_Component_Inter_Output_Standard_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Inter_Output_Standard_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the exponential disk node mass rate of change.
    !!}
    use :: Galacticus_Nodes, only : defaultInteroutputComponent, interruptTask        , nodeComponentBasic  , nodeComponentDisk, &
          &                         nodeComponentInterOutput   , nodeComponentSpheroid, propertyTypeInactive, treeNode
    implicit none
    type            (treeNode                    ), intent(inout)          :: node
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    integer                                       , intent(in   )          :: propertyType
    class           (nodeComponentInterOutput    )               , pointer :: interOutput
    class           (nodeComponentDisk           )               , pointer :: disk
    class           (nodeComponentSpheroid       )               , pointer :: spheroid
    class           (nodeComponentBasic          )               , pointer :: basic
    double precision                                                       :: diskStarFormationRate, spheroidStarFormationRate, &
         &                                                                    timeCurrent          , timeOutputNext           , &
         &                                                                    timeOutputPrevious

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if the standard inter-output component is not active.
    if (.not.defaultInteroutputComponent%standardIsActive()) return
    ! Get the disk and check that it is of our class.
    interOutput => node%interOutput()
    ! Get disk and spheroid star formation rates.
    disk                      => node                       %disk    (    )
    spheroid                  => node                       %spheroid(    )
    diskStarFormationRate     =  starFormationRateDisks_    %rate    (node)
    spheroidStarFormationRate =  starFormationRateSpheroids_%rate    (node)
    ! Find the time interval between previous and next outputs.
    basic              => node %basic()
    timeCurrent        =  basic%time ()
    timeOutputPrevious =  outputTimes_%timePrevious(timeCurrent)
    timeOutputNext     =  outputTimes_%timeNext    (timeCurrent)
    ! Return if there is no next output.
    if (timeOutputNext     < 0.0d0) return
    ! Set previous time to zero if there is no previous output.
    if (timeOutputPrevious < 0.0d0) timeOutputPrevious=0.0d0
    ! Accumulate rates.
    call interOutput%    diskStarFormationRateRate(    diskStarFormationRate/(timeOutputNext-timeOutputPrevious),interrupt,interruptProcedure)
    call interOutput%spheroidStarFormationRateRate(spheroidStarFormationRate/(timeOutputNext-timeOutputPrevious),interrupt,interruptProcedure)
    return
  end subroutine Node_Component_Inter_Output_Standard_Rate_Compute

  !![
  <mergerTreeExtraOutputTask>
   <unitName>Node_Component_Inter_Output_Standard_Reset</unitName>
  </mergerTreeExtraOutputTask>
  !!]
  subroutine Node_Component_Inter_Output_Standard_Reset(node,iOutput,treeIndex,nodePassesFilter,treeLock)
    !!{
    Reset interoutput accumulated quantities.
    !!}
    use            :: Galacticus_Nodes, only : nodeComponentInterOutput, nodeComponentInterOutputStandard, treeNode, defaultInterOutputComponent
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Kind_Numbers    , only : kind_int8
    use            :: Locks           , only : ompLock
    implicit none
    type   (treeNode                ), intent(inout), pointer :: node
    integer(c_size_t                ), intent(in   )          :: iOutput
    integer(kind=kind_int8          ), intent(in   )          :: treeIndex
    logical                          , intent(in   )          :: nodePassesFilter
    type   (ompLock                 ), intent(inout)          :: treeLock
    class  (nodeComponentInterOutput)               , pointer :: interOutput
    !$GLC attributes unused :: iOutput, nodePassesFilter, treeIndex, treeLock
    
    ! Check if we are the default method.
    if (.not.defaultInterOutputComponent%standardIsActive()) return
    ! Get the interoutput component and check it is of our class.
    interOutput => node%interOutput()
    select type (interOutput)
       class is (nodeComponentInterOutputStandard)
       call interOutput%    diskStarFormationRateSet(0.0d0)
       call interOutput%spheroidStarFormationRateSet(0.0d0)
    end select
    return
  end subroutine Node_Component_Inter_Output_Standard_Reset

  subroutine satelliteMerger(self,node)
    !!{
    Remove any inter-output quantities associated with {\normalfont \ttfamily node} and add them to the merge target.
    !!}
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentInterOutput, nodeComponentInterOutputStandard, treeNode
    use :: Satellite_Merging_Mass_Movements    , only : destinationMergerDisk   , destinationMergerSpheroid       , destinationMergerUnmoved
    implicit none
    class  (*                       ), intent(inout) :: self
    type   (treeNode                ), intent(inout) :: node
    type   (treeNode                ), pointer       :: nodeHost
    class  (nodeComponentInterOutput), pointer       :: interOutputHost        , interOutput
    integer                                          :: destinationGasSatellite, destinationGasHost       , &
         &                                              destinationStarsHost   , destinationStarsSatellite
    logical                                          :: mergerIsMajor
    !$GLC attributes unused :: self

    ! Get the inter-output component.
    interOutput => node%interOutput()
    ! Ensure that it is of the standard class.
    select type (interOutput)
    class is (nodeComponentInterOutputStandard)
       ! Find the node to merge with.
       nodeHost        => node    %mergesWith ()
       interOutputHost => nodeHost%interOutput()
       ! Get mass movement descriptors.
       call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       ! Move the star formation rates from secondary to primary.
       select case (destinationStarsSatellite)
       case (destinationMergerDisk    )
          call interOutputHost%    diskStarFormationRateSet(                                             &
               &                                             interOutputHost%    diskStarFormationRate() &
               &                                            +interOutput    %    diskStarFormationRate() &
               &                                           )
          call interOutputHost%spheroidStarFormationRateSet(                                             &
               &                                             interOutputHost%spheroidStarFormationRate() &
               &                                            +interOutput    %spheroidStarFormationRate() &
               &                                           )
       case (destinationMergerSpheroid)
       case default
          call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       ! Zero rates in the secondary,
       call interOutput%    diskStarFormationRateSet(                                             &
            &                                         0.0d0                                       &
            &                                       )
       call interOutput%spheroidStarFormationRateSet(                                             &
            &                                         0.0d0                                       &
            &                                       )
       ! Move star formation rates within the host if necessary.
       select case (destinationStarsHost)
       case (destinationMergerDisk)
          call interOutputHost%    diskStarFormationRateSet(                                             &
               &                                             interOutputHost%    diskStarFormationRate() &
               &                                            +interOutputHost%spheroidStarFormationRate() &
               &                                           )
          call interOutputHost%spheroidStarFormationRateSet(                                             &
               &                                             0.0d0                                       &
               &                                           )
       case (destinationMergerSpheroid)
          call interOutputHost%spheroidStarFormationRateSet(                                             &
               &                                             interOutputHost%spheroidStarFormationRate() &
               &                                            +interOutputHost%    diskStarFormationRate() &
               &                                           )
          call interOutputHost%    diskStarFormationRateSet(                                             &
               &                                             0.0d0                                       &
               &                                           )
       case (destinationMergerUnmoved)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
    end select
    return
  end subroutine satelliteMerger

end module Node_Component_Inter_Output_Standard
