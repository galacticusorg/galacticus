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

!% Contains a module which implements the preset spin component.

module Node_Component_Spin_Preset
  !% Implements the preset spin component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Spin_Preset_Initialize  , Node_Component_Spin_Preset_Promote     , &
       &    Node_Component_Spin_Preset_Scale_Set   , Node_Component_Spin_Preset_Rate_Compute

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
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    class           (nodeComponentSpin )               , pointer :: spinParent , spin
    class           (nodeComponentBasic)               , pointer :: basicParent, basic
    double precision                                             :: deltaTime

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

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Spin_Preset_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Spin_Preset_Rate_Compute(node,odeConverged,interrupt,interruptProcedure)
    !% Compute rates of change of properties in the preset implementation of the spin component.
    implicit none
    type     (treeNode         ), intent(inout), pointer :: node
    logical                     , intent(in   )          :: odeConverged
    logical                     , intent(inout)          :: interrupt
    procedure(                 ), intent(inout), pointer :: interruptProcedure
    class    (nodeComponentSpin)               , pointer :: spin
    !GCC$ attributes unused :: interrupt, interruptProcedure, odeConverged
    
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
    implicit none
    type            (treeNode         ), intent(inout), pointer :: node
    double precision                   , parameter              :: timeScale        =1.0d-3
    double precision                   , parameter              :: scaleMassRelative=1.0d-6
    class           (nodeComponentSpin)               , pointer :: spin

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

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Spin_Preset_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Spin_Preset_Promote(node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the spin of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentSpin )               , pointer :: spinParent , spin
    class(nodeComponentBasic)               , pointer :: basicParent, basic

    ! Ensure that the spin component is of the preset class.
    spin => node%spin()
    select type (spin)
    class is (nodeComponentSpinPreset)
       nodeParent => node%parent
       basic      => node%basic()
       basicParent => nodeParent%basic()
       if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('Node_Component_Spin_Random_Promote','node has not been evolved to its parent')
       ! Adjust the spin growth rate to that of the parent node.
       spinParent => nodeParent%spin()
       call spin%spinSet          (spinParent%spin          ())
       call spin%spinGrowthRateSet(spinParent%spinGrowthRate())
    end select
    return
  end subroutine Node_Component_Spin_Preset_Promote

end module Node_Component_Spin_Preset
