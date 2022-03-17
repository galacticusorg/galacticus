!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements an extension of the standard hot halo node component which tracks the metals arriving from
outflows.
!!}

module Node_Component_Hot_Halo_Outflow_Tracking
  !!{
  Implements an extension of the standard hot halo node component which tracks the metals arriving from
  outflows.
  !!}
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  implicit none
  private
  public :: Node_Component_Hot_Halo_Outflow_Tracking_Rate_Compute     , Node_Component_Hot_Halo_Outflow_Tracking_Scale_Set          , &
       &    Node_Component_Hot_Halo_Outflow_Tracking_Thread_Initialize, Node_Component_Hot_Halo_Outflow_Tracking_Thread_Uninitialize, &
       &    Node_Component_Hot_Halo_Outflow_Tracking_State_Store      , Node_Component_Hot_Halo_Outflow_Tracking_State_Restore

  !![
  <component>
   <class>hotHalo</class>
   <name>outflowTracking</name>
   <extends>
    <class>hotHalo</class>
    <name>standard</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>trackedOutflowMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass in the hot phase of the hot halo arrived via direct outflow."/>
    </property>
    <property>
      <name>trackedOutflowAbundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the hot phase of the hot halo arrived via direct outflow."/>
    </property>
   </properties>
   <bindings>
    <binding method="massRemovalRate" function="Node_Component_Hot_Halo_Outflow_Tracking_Mass_Removal_Rate" bindsTo="component" />
   </bindings>
   <functions>objects.nodes.components.hot_halo.outflow_tracking.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
  !$omp threadprivate(darkMatterHaloScale_)

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Hot_Halo_Outflow_Tracking_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Outflow_Tracking_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node very simple disk profile module.
    !!}
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultHotHaloComponent%outflowTrackingIsActive()) then
       !![
       <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Hot_Halo_Outflow_Tracking_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Hot_Halo_Outflow_Tracking_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Outflow_Tracking_Thread_Uninitialize()
    !!{
    Uninitializes the tree node very simple disk profile module.
    !!}
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none

    if (defaultHotHaloComponent%outflowTrackingIsActive()) then
       !![
       <objectDestructor name="darkMatterHaloScale_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Hot_Halo_Outflow_Tracking_Thread_Uninitialize

  !![
  <rateComputeTask>
   <unitName>Node_Component_Hot_Halo_Outflow_Tracking_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Hot_Halo_Outflow_Tracking_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the hot halo node mass rate of change.
    !!}
    use :: Abundances_Structure                 , only : abundances              , operator(*)
    use :: Galacticus_Nodes                     , only : defaultHotHaloComponent , interruptTask, nodeComponentHotHalo, nodeComponentHotHaloOutflowTracking, &
          &                                              propertyInactive        , treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : hotHaloOutflowReturnRate
    implicit none
    type            (treeNode                    ), intent(inout)          :: node
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    integer                                       , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo        )               , pointer :: hotHalo
    double precision                                                       :: massReturnRate
    type            (abundances                  )                         :: abundancesReturnRate
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%outflowTrackingIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Act only if this hot halo is of our class.
    select type (hotHalo)
    class is (nodeComponentHotHaloOutflowTracking)
       ! Add the rate of abundances return from the outflowed component.
       massReturnRate      =hotHaloOutflowReturnRate*hotHalo%outflowedMass      ()/darkMatterHaloScale_%timescaleDynamical(node)
       abundancesReturnRate=hotHaloOutflowReturnRate*hotHalo%outflowedAbundances()/darkMatterHaloScale_%timescaleDynamical(node)
       call hotHalo%trackedOutflowMassRate      (      massReturnRate)
       call hotHalo%trackedOutflowAbundancesRate(abundancesReturnRate)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Outflow_Tracking_Rate_Compute

  !![
  <scaleSetTask>
   <unitName>Node_Component_Hot_Halo_Outflow_Tracking_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Hot_Halo_Outflow_Tracking_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure, only : unitAbundances
    use :: Galacticus_Nodes    , only : nodeComponentBasic     , nodeComponentHotHalo, nodeComponentHotHaloOutflowTracking, treeNode, &
         &                              defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                      , parameter              :: scaleMassRelative=1.0d-3
    double precision                                               :: massVirial

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%outflowTrackingIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of the standard class.
    select type (hotHalo)
    class is (nodeComponentHotHaloOutflowTracking)
       ! Get the basic component.
       basic => node%basic()
       ! Get virial properties.
       massVirial=basic%mass()
       ! Set a scale for the tracked abundances.
       call hotHalo%trackedOutflowMassScale      (               massVirial*scaleMassRelative)
       call hotHalo%trackedOutflowAbundancesScale(unitAbundances*massVirial*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Outflow_Tracking_Scale_Set

  !![
  <stateStoreTask>
   <unitName>Node_Component_Hot_Halo_Outflow_Tracking_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Hot_Halo_Outflow_Tracking_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentHotHalo -> outflowTracking',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Outflow_Tracking_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Hot_Halo_Outflow_Tracking_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Hot_Halo_Outflow_Tracking_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentHotHalo -> outflowTracking',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Outflow_Tracking_State_Restore

end module Node_Component_Hot_Halo_Outflow_Tracking
