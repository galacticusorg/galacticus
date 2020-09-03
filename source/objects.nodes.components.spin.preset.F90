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

!% Contains a module which implements the preset spin component.

module Node_Component_Spin_Preset
  !% Implements the preset spin component.
  implicit none
  private
  public :: Node_Component_Spin_Preset_Initialize         , Node_Component_Spin_Preset_Rate_Compute     , &
       &    Node_Component_Spin_Preset_Scale_Set          , Node_Component_Spin_Preset_Thread_Initialize, &
       &    Node_Component_Spin_Preset_Thread_Uninitialize

  !# <component>
  !#  <class>spin</class>
  !#  <name>preset</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>spin</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="0.0d0" comment="Spin parameter of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>spinGrowthRate</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Spin_Preset_Initialize</unitName>
  !#  <sortName>spin</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Spin_Preset_Initialize(node)
    !% Initialize the spin of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentBasic  , nodeComponentSpin, nodeComponentSpinPreset, treeNode, &
         &                          defaultSpinComponent
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    class           (nodeComponentSpin )               , pointer :: spinParent , spin
    class           (nodeComponentBasic)               , pointer :: basicParent, basic
    double precision                                             :: deltaTime

    ! Check if we are the default method.
    if (.not.defaultSpinComponent%presetIsActive()) return
    ! Ensure that the spin component is of the preset class.
    spin => node%spin()
    select type (spin)
    class is (nodeComponentSpinPreset)
       ! Check if this node is the primary progenitor.
       if (node%isPrimaryProgenitor()) then
          ! It is, so compute the spin growth rate.
          basic   => node       %basic()
          basicParent => node%parent%basic()
          spinParent  => node%parent%spin ()
          deltaTime=basicParent%time()-basic%time()
          if (deltaTime > 0.0d0) then
             call spin%spinGrowthRateSet((spinParent%spin()-spin%spin())/deltaTime)
          else
             call spin%spinGrowthRateSet(0.0d0)
          end if
       else
          ! It is not, so set spin growth rate to zero.
          call spin%spinGrowthRateSet(0.0d0)
       end if
    end select
    return
  end subroutine Node_Component_Spin_Preset_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Spin_Preset_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Spin_Preset_Thread_Initialize(parameters_)
    !% Initializes the tree node preset spin module.
    use :: Events_Hooks    , only : nodePromotionEvent  , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultSpinComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_
    
    if (defaultSpinComponent%presetIsActive()) &
         & call nodePromotionEvent%attach(defaultSpinComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSpinPreset')
    return
  end subroutine Node_Component_Spin_Preset_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Spin_Preset_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Spin_Preset_Thread_Uninitialize()
    !% Uninitializes the tree node preset spin module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (defaultSpinComponent%presetIsActive()) &
         & call nodePromotionEvent%detach(defaultSpinComponent,nodePromotion)
    return
  end subroutine Node_Component_Spin_Preset_Thread_Uninitialize
  
  !# <rateComputeTask>
  !#  <unitName>Node_Component_Spin_Preset_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Spin_Preset_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Compute rates of change of properties in the preset implementation of the spin component.
    use :: Galacticus_Nodes, only : defaultSpinComponent, nodeComponentSpin, nodeComponentSpinPreset, propertyTypeInactive, &
          &                         treeNode
    implicit none
    type     (treeNode         ), intent(inout), pointer :: node
    logical                     , intent(inout)          :: interrupt
    procedure(                 ), intent(inout), pointer :: interruptProcedure
    integer                     , intent(in   )          :: propertyType
    class    (nodeComponentSpin)               , pointer :: spin
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%presetIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the preset class.
    select type (spin)
    class is (nodeComponentSpinPreset)
       ! Rate of change is set equal to the precomputed growth rate.
       call spin%spinRate(spin%spinGrowthRate())
    end select
    return
  end subroutine Node_Component_Spin_Preset_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Spin_Preset_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Spin_Preset_Scale_Set(node)
    !% Set scales for properties in the preset implementation of the spin component.
    use :: Galacticus_Nodes, only : nodeComponentSpin, nodeComponentSpinPreset, treeNode, defaultSpinComponent
    implicit none
    type (treeNode         ), intent(inout), pointer :: node
    class(nodeComponentSpin)               , pointer :: spin

    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%presetIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the preset class.
    select type (spin)
    class is (nodeComponentSpinPreset)
       ! Set scale for spin.
       call spin%spinScale(spin%spin())
    end select
    return
  end subroutine Node_Component_Spin_Preset_Scale_Set

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the spin of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentSpin, treeNode
    implicit none
    class(*                 ), intent(inout)          :: self
    type (treeNode          ), intent(inout), target  :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentSpin )               , pointer :: spinParent , spin
    class(nodeComponentBasic)               , pointer :: basicParent, basic
    !$GLC attributes unused :: self

    spin        => node      %spin  ()
    nodeParent  => node      %parent
    basic       => node      %basic ()
    basicParent => nodeParent%basic ()
    if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust the spin growth rate to that of the parent node.
    spinParent => nodeParent%spin()
    call spin%spinSet          (spinParent%spin          ())
    call spin%spinGrowthRateSet(spinParent%spinGrowthRate())
    return
  end subroutine nodePromotion

end module Node_Component_Spin_Preset
