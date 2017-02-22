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

!% Contains a module which implements the preset 3-D spin component.

module Node_Component_Spin_Preset3D
  !% Implements the preset spin component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Spin_Preset3D_Initialize  , Node_Component_Spin_Preset3D_Promote     , &
       &    Node_Component_Spin_Preset3D_Scale_Set   , Node_Component_Spin_Preset3D_Rate_Compute

  !# <component>
  !#  <class>spin</class>
  !#  <name>preset3D</name>
  !#  <extends>
  !#   <class>spin</class>
  !#   <name>preset</name>
  !#  </extends>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>spinVector</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output labels="[X,Y,Z]" unitsInSI="0.0d0" comment="Spin vector of the node."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>spinVectorGrowthRate</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Spin_Preset3D_Initialize</unitName>
  !#  <sortName>spin</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Spin_Preset3D_Initialize(node)
    !% Initialize the spin of {\normalfont \ttfamily node}.
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    class           (nodeComponentSpin )               , pointer :: parentSpinComponent , spin
    class           (nodeComponentBasic)               , pointer :: parentBasicComponent, basic
    double precision                                             :: deltaTime

    ! Ensure that the spin component is of the preset class.
    spin => node%spin()
    select type (spin)
    class is (nodeComponentSpinPreset3D)
       ! Check if this node is the primary progenitor.
       if (node%isPrimaryProgenitor()) then
          ! It is, so compute the spin vector growth rate.
          basic   => node       %basic()
          parentBasicComponent => node%parent%basic()
          parentSpinComponent  => node%parent%spin ()
          deltaTime=parentBasicComponent%time()-basic%time()
          if (deltaTime > 0.0d0) then
             call spin%spinVectorGrowthRateSet((parentSpinComponent%spinVector()-spin%spinVector())/deltaTime)
          else
             call spin%spinVectorGrowthRateSet([0.0d0,0.0d0,0.0d0])
          end if
       else
          ! It is not, so set spin growth rate to zero.
          call spin%spinVectorGrowthRateSet([0.0d0,0.0d0,0.0d0])
       end if
    end select
    return
  end subroutine Node_Component_Spin_Preset3D_Initialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Spin_Preset3D_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Spin_Preset3D_Rate_Compute(node,odeConverged,interrupt,interruptProcedure)
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
    class is (nodeComponentSpinPreset3D)
       ! Rate of change is set equal to the precomputed growth rate.
       call spin%spinVectorRate(spin%spinVectorGrowthRate())
    end select
    return
  end subroutine Node_Component_Spin_Preset3D_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Spin_Preset3D_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Spin_Preset3D_Scale_Set(node)
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
    class is (nodeComponentSpinPreset3D)
       ! Set scale for spin.
       call spin%spinVectorScale([1.0d0,1.0d0,1.0d0]*spin%spin())
    end select
    return
  end subroutine Node_Component_Spin_Preset3D_Scale_Set

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Spin_Preset3D_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Spin_Preset3D_Promote(node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the spin of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentSpin )               , pointer :: spinParent , spin
    class(nodeComponentBasic)               , pointer :: basicParent, basic

    ! Ensure that the spin component is of the preset3D class.
    spin => node%spin()
    select type (spin)
    class is (nodeComponentSpinPreset3D)
       nodeParent  => node      %parent
       basic       => node      %basic ()
       basicParent => nodeParent%basic ()
       if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('Node_Component_Spin_Preset3D_Promote','node has not been evolved to its parent')
       ! Adjust the spin growth rate to that of the parent node.
       spinParent => nodeParent%spin()
       call spin%spinVectorSet          (spinParent%spinVector          ())
       call spin%spinVectorGrowthRateSet(spinParent%spinVectorGrowthRate())
    end select
    return
  end subroutine Node_Component_Spin_Preset3D_Promote
  
end module Node_Component_Spin_Preset3D
