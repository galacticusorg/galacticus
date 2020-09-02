!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements the standard mass flow statistics component.

module Node_Component_Mass_Flow_Statistics_Standard
  !% Implements the standard mass flow statistics component.
  use :: Cooling_Rates, only : coolingRateClass
  implicit none
  private
  public :: Node_Component_Mass_Flow_Statistics_Standard_Merger_Tree_Init , Node_Component_Mass_Flow_Statistics_Standard_Scale_Set    , &
       &    Node_Component_Mass_Flow_Statistics_Standard_Extra_Output     , Node_Component_Mass_Flow_Statistics_Standard_Rate_Compute , &
       &    Node_Component_Mass_Flow_Statistics_Standard_Thread_Initialize, Node_Component_Mass_Flow_Statistics_Standard_Thread_Uninit, &
       &    Node_Component_Mass_Flow_Statistics_Standard_Initialize

  !# <component>
  !#  <class>massFlowStatistics</class>
  !#  <name>standard</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>cooledMass</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Cumulative mass of gas cooled directly onto galaxy."/>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Objects used by this component.
  class(coolingRateClass), pointer :: coolingRate_
  !$omp threadprivate(coolingRate_)

  ! Options controlling module behavior.
  logical :: massFlowStatisticsResetOnOutput

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Initialize(parameters_)
    !% Initializes the standard mass flow statistics component.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    !# <inputParameter>
    !#   <name>massFlowStatisticsResetOnOutput</name>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>Specifies whether or not mass flow statistics should be reset to zero at each output.</description>
    !#   <source>parameters_</source>
    !# </inputParameter>
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Thread_Initialize(parameters_)
    !% Initializes the tree node standard mass flow statistics module.
    use :: Galacticus_Nodes, only : defaultMassFlowStatisticsComponent
    use :: Input_Parameters, only : inputParameter                    , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultMassFlowStatisticsComponent%standardIsActive()) then
       !# <objectBuilder class="coolingRate" name="coolingRate_" source="parameters_"/>
    end if
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Thread_Uninit</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Thread_Uninit()
    !% Uninitializes the tree node standard mass flow statistics module.
    use :: Galacticus_Nodes, only : defaultMassFlowStatisticsComponent
    implicit none

    if (defaultMassFlowStatisticsComponent%standardIsActive()) then
       !# <objectDestructor name="coolingRate_"/>
    end if
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Thread_Uninit

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Merger_Tree_Init</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Merger_Tree_Init(node)
    !% Initialize the mass flow statistics component by creating components in nodes and computing formation times.
    use :: Galacticus_Nodes, only : defaultMassFlowStatisticsComponent, nodeComponentMassFlowStatistics, nodeComponentMassFlowStatisticsStandard, treeNode
    implicit none
    type (treeNode                       ), pointer, intent(inout) :: node
    class(nodeComponentMassFlowStatistics), pointer                :: massFlowStatistics

    ! Return immediately if this class is not active.
    if (.not.defaultMassFlowStatisticsComponent%standardIsActive()) return
    ! Create a mass flow statistics component and initialize it.
    massFlowStatistics => node%massFlowStatistics(autoCreate=.true.)
    select type (massFlowStatistics)
    class is (nodeComponentMassFlowStatisticsStandard)
       call massFlowStatistics%cooledMassSet(0.0d0)
    end select
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Merger_Tree_Init

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Rate_Compute(node,odeConverged,interrupt,interruptProcedure,propertyType)
    !% Compute rates of change of properties in the standard implementation of the basic component.
    use :: Galacticus_Nodes, only : defaultMassFlowStatisticsComponent, nodeComponentMassFlowStatistics, nodeComponentMassFlowStatisticsStandard, propertyTypeInactive, &
          &                         treeNode
    implicit none
    type     (treeNode                       ), pointer, intent(inout) :: node
    logical                                            , intent(in   ) :: odeConverged
    logical                                   ,          intent(inout) :: interrupt
    procedure(                               ), pointer, intent(inout) :: interruptProcedure
    integer                                   , intent(in   )          :: propertyType
    class    (nodeComponentMassFlowStatistics), pointer                :: massFlowStatistics
    !$GLC attributes unused :: interrupt, interruptProcedure, odeConverged

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultMassFlowStatisticsComponent%standardIsActive()) return
    ! Get the massFlowStatistics component.
    massFlowStatistics => node%massFlowStatistics()
    ! Ensure that it is of the standard class.
    select type (massFlowStatistics)
    class is (nodeComponentMassFlowStatisticsStandard)
       ! Cooled mass rate simply equals the cooling rate.
       call massFlowStatistics%cooledMassRate(coolingRate_%rate(node))
    end select
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Scale_Set(node)
    !% Set scales for properties in the standard implementation of the massFlowStatistics component.
    use :: Galacticus_Nodes, only : nodeComponentBasic                , nodeComponentMassFlowStatistics, nodeComponentMassFlowStatisticsStandard, treeNode, &
         &                          defaultMassFlowStatisticsComponent
    implicit none
    type            (treeNode                       ), pointer, intent(inout) :: node
    double precision                                 , parameter              :: scaleMassRelative =1.0d-6
    class           (nodeComponentMassFlowStatistics), pointer                :: massFlowStatistics
    class           (nodeComponentBasic             ), pointer                :: basic

    ! Check if we are the default method.
    if (.not.defaultMassFlowStatisticsComponent%standardIsActive()) return
    ! Get the massFlowStatistics component.
    massFlowStatistics => node%massFlowStatistics()
    ! Ensure that it is of the standard class.
    select type (massFlowStatistics)
    class is (nodeComponentMassFlowStatisticsStandard)
       ! Set scale for cooled mass.
       basic => node%basic()
       call massFlowStatistics%cooledMassScale(scaleMassRelative*basic%mass())
    end select
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Scale_Set

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Extra_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Extra_Output(node,iOutput,treeIndex,nodePassesFilter)
    !% Reset mass flow statistics at output time.
    use            :: Galacticus_Nodes, only : nodeComponentMassFlowStatistics, nodeComponentMassFlowStatisticsStandard, treeNode, defaultMassFlowStatisticsComponent
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Kind_Numbers    , only : kind_int8
    implicit none
    type            (treeNode                       ), intent(inout), pointer :: node
    integer         (kind=kind_int8                 ), intent(in   )          :: treeIndex
    integer         (kind=c_size_t                  ), intent(in   )          :: iOutput
    logical                                          , intent(in   )          :: nodePassesFilter
    class           (nodeComponentMassFlowStatistics),                pointer :: massFlowStatistics
    !$GLC attributes unused :: iOutput, nodePassesFilter, treeIndex

    ! Check if we are the default method.
    if (.not.defaultMassFlowStatisticsComponent%standardIsActive()) return
    ! Return immediately if we are not to reset mass flow statistics at output time.
    if (.not.massFlowStatisticsResetOnOutput) return
    ! Get the massFlowStatistics component.
    massFlowStatistics => node%massFlowStatistics()
    ! Ensure that it is of the standard class.
    select type (massFlowStatistics)
    class is (nodeComponentMassFlowStatisticsStandard)
       call massFlowStatistics%cooledMassSet(0.0d0)
    end select
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Extra_Output

end module Node_Component_Mass_Flow_Statistics_Standard
