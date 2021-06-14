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

!% Contains a module which implements the a vector spin component for tree nodes.

module Node_Component_Spin_Vector
  !% Implements the vector spin component.
  implicit none
  private
  public :: Node_Component_Spin_Vector_Thread_Initialize, Node_Component_Spin_Vector_Thread_Uninitialize, &
       &    Node_Component_Spin_Vector_Scale_Set
  
  !# <component>
  !#  <class>spin</class>
  !#  <name>vector</name>
  !#  <extends>
  !#   <class>spin</class>
  !#   <name>scalar</name>
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

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Spin_Vector_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Spin_Vector_Thread_Initialize(parameters_)
    !% Initializes the tree node vector spin module.
    use :: Events_Hooks    , only : nodePromotionEvent  , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultSpinComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_
    
    if (defaultSpinComponent%vectorIsActive()) &
         & call nodePromotionEvent%attach(defaultSpinComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSpinVector')
    return
  end subroutine Node_Component_Spin_Vector_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Spin_Vector_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Spin_Vector_Thread_Uninitialize()
    !% Uninitializes the tree node preset spin module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (defaultSpinComponent%vectorIsActive()) &
         & call nodePromotionEvent%detach(defaultSpinComponent,nodePromotion)
    return
  end subroutine Node_Component_Spin_Vector_Thread_Uninitialize
 
  !# <scaleSetTask>
  !#  <unitName>Node_Component_Spin_Vector_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Spin_Vector_Scale_Set(node)
    !% Set scales for properties in the preset implementation of the spin component.
    use :: Galacticus_Nodes, only : nodeComponentSpin, nodeComponentSpinVector, treeNode, defaultSpinComponent
    implicit none
    type (treeNode         ), intent(inout), pointer :: node
    class(nodeComponentSpin)               , pointer :: spin

    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%vectorIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the preset class.
    select type (spin)
    class is (nodeComponentSpinVector)
       ! Set scale for spin.
       call spin%spinVectorScale([1.0d0,1.0d0,1.0d0]*spin%spin())
    end select
    return
  end subroutine Node_Component_Spin_Vector_Scale_Set

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
    call spin%spinVectorSet          (spinParent%spinVector          ())
    call spin%spinVectorGrowthRateSet(spinParent%spinVectorGrowthRate())
    return
  end subroutine nodePromotion

end module Node_Component_Spin_Vector
