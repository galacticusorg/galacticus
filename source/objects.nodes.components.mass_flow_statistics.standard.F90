!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Mass_Flow_Statistics_Standard_Merger_Tree_Init, Node_Component_Mass_Flow_Statistics_Standard_Scale_Set   , &
       &    Node_Component_Mass_Flow_Statistics_Standard_Extra_Output    , Node_Component_Mass_Flow_Statistics_Standard_Rate_Compute

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

  ! Options controlling module behavior.
  logical :: massFlowStatisticsResetOnOutput

  ! Record of whether this module has been initialized.
  logical :: moduleInitialized=.false.

contains

  subroutine Node_Component_Mass_Flow_Statistics_Standard_Initialize()
    !% Initializes the standard mass flow statistics component.
    use Input_Parameters
    implicit none

    ! Test whether module is already initialize.
    !$omp critical (Node_Component_Mass_Flow_Statistics_Standard_Initialize)
    if (.not.moduleInitialized) then
       !@ <inputParameter>
       !@   <name>massFlowStatisticsResetOnOutput</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not mass flow statistics should be reset to zero at each output.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('massFlowStatisticsResetOnOutput',massFlowStatisticsResetOnOutput,defaultValue=.true.)
       ! Record that the module is now initialized.
       moduleInitialized=.true.     
    end if
    !$omp end critical (Node_Component_Mass_Flow_Statistics_Standard_Initialize)
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Initialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Merger_Tree_Init</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Merger_Tree_Init(node)
    !% Initialize the mass flow statistics component by creating components in nodes and computing formation times.
    use Dark_Matter_Halo_Formation_Times
    implicit none
    type (treeNode                       ), pointer, intent(inout) :: node
    class(nodeComponentMassFlowStatistics), pointer                :: massFlowStatistics
    
    ! Return immediately if this class is not active.
    if (.not.defaultMassFlowStatisticsComponent%standardIsActive()) return
    
    ! Ensure that the module is initialized.
    call Node_Component_Mass_Flow_Statistics_Standard_Initialize()

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
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Rate_Compute(node,odeConverged,interrupt,interruptProcedure)
    !% Compute rates of change of properties in the standard implementation of the basic component.
    use Cooling_Rates
    implicit none
    type     (treeNode                       ), pointer, intent(inout) :: node
    logical                                            , intent(in   ) :: odeConverged
    logical                                   ,          intent(inout) :: interrupt
    procedure(                               ), pointer, intent(inout) :: interruptProcedure
    class    (nodeComponentMassFlowStatistics), pointer                :: massFlowStatistics
    !GCC$ attributes unused :: interrupt, interruptProcedure, odeConverged
    
    ! Get the massFlowStatistics component.
    massFlowStatistics => node%massFlowStatistics()
    ! Ensure that it is of the standard class.
    select type (massFlowStatistics)
    class is (nodeComponentMassFlowStatisticsStandard)
       ! Cooled mass rate simply equals the cooling rate.
       call massFlowStatistics%cooledMassRate(Cooling_Rate(node))
    end select
    return
  end subroutine Node_Component_Mass_Flow_Statistics_Standard_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Mass_Flow_Statistics_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Mass_Flow_Statistics_Standard_Scale_Set(node)
    !% Set scales for properties in the standard implementation of the massFlowStatistics component.
    implicit none
    type            (treeNode                       ), pointer, intent(inout) :: node
    double precision                                 , parameter              :: scaleMassRelative =1.0d-6
    class           (nodeComponentMassFlowStatistics), pointer                :: massFlowStatistics
    class           (nodeComponentBasic             ), pointer                :: basic

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
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Nodes
    use Kind_Numbers
    implicit none
    type            (treeNode                       ), intent(inout), pointer :: node
    integer         (kind=kind_int8                 ), intent(in   )          :: treeIndex
    integer         (kind=c_size_t                  ), intent(in   )          :: iOutput
    logical                                          , intent(in   )          :: nodePassesFilter
    class           (nodeComponentMassFlowStatistics),                pointer :: massFlowStatistics
    !GCC$ attributes unused :: iOutput, nodePassesFilter, treeIndex
    
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
