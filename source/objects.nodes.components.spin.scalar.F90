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
Contains a module implementing a scalar spin component for tree nodes.
!!}

module Node_Component_Spin_Scalar
  !!{
  Implement a scalar spin component for tree nodes.
  !!}
  implicit none
  private
  public :: Node_Component_Spin_Scalar_Thread_Initialize, Node_Component_Spin_Scalar_Thread_Uninitialize, &
       &    Node_Component_Spin_Scalar_Scale_Set

  !![
  <component>
   <class>spin</class>
   <name>scalar</name>
   <isDefault>true</isDefault>
   <properties>
    <property>
      <name>spin</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="0.0d0" comment="Spin parameter of the node."/>
    </property>
    <property>
      <name>spinGrowthRate</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
  </component>
  !!]

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Spin_Scalar_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Spin_Scalar_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node scalar spin module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent  , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultSpinComponent
    use :: Input_Parameters, only : inputParameter      , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultSpinComponent%scalarIsActive()) &
         & call nodePromotionEvent%attach(defaultSpinComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSpinScalar')
    return
  end subroutine Node_Component_Spin_Scalar_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Spin_Scalar_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Spin_Scalar_Thread_Uninitialize()
    !!{
    Uninitializes the tree node scalar spin module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (defaultSpinComponent%scalarIsActive() .and. nodePromotionEvent%isAttached(defaultSpinComponent,nodePromotion)) &
         & call nodePromotionEvent%detach(defaultSpinComponent,nodePromotion)
    return
  end subroutine Node_Component_Spin_Scalar_Thread_Uninitialize

  !![
  <scaleSetTask>
   <unitName>Node_Component_Spin_Scalar_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Spin_Scalar_Scale_Set(node)
    !!{
    Set scales for properties in the scalar implementation of the spin component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpin, nodeComponentSpinScalar, treeNode, defaultSpinComponent
    implicit none
    type            (treeNode         ), intent(inout), pointer :: node
    class           (nodeComponentSpin)               , pointer :: spin
    double precision                   , parameter              :: spinMinimum=1.0d-6

    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%scalarIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the scalar class.
    select type (spin)
    class is (nodeComponentSpinScalar)
       ! Set scale for spin.
       call spin%spinScale(max(spin%spin(),spinMinimum))
    end select
    return
  end subroutine Node_Component_Spin_Scalar_Scale_Set

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the spin of {\normalfont \ttfamily node}
    to be that of its parent.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentSpin, treeNode
    implicit none
    class(*                 ), intent(inout) :: self
    type (treeNode          ), intent(inout) :: node
    type (treeNode          ), pointer       :: nodeParent
    class(nodeComponentSpin ), pointer       :: spinParent , spin
    class(nodeComponentBasic), pointer       :: basicParent, basic
    !$GLC attributes unused :: self

    nodeParent  => node      %parent
    spin        => node      %spin  ()
    spinParent  => nodeParent%spin  ()
    basic       => node      %basic ()
    basicParent => nodeParent%basic ()
    call spin%spinSet          (spinParent%spin          ())
    call spin%spinGrowthRateSet(spinParent%spinGrowthRate())
    return
  end subroutine nodePromotion

end module Node_Component_Spin_Scalar
