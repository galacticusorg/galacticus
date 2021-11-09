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
Contains a module which implements the a vector angular momentum component for halos.
!!}

module Node_Component_Halo_Angular_Momentum_Vector
  !!{
  Implements the vector spin component.
  !!}
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  implicit none
  private
  public :: Node_Component_Halo_Angular_Momentum_Vector_Thread_Initialize, Node_Component_Halo_Angular_Momentum_Vector_Thread_Uninitialize, &
       &    Node_Component_Halo_Angular_Momentum_Vector_Scale_Set
  
  !![
  <component>
   <class>spin</class>
   <name>vector</name>
   <extends>
    <class>spin</class>
    <name>scalar</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>angularMomentumVector</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output labels="[X,Y,Z]" unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum vector of the DMO halo."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>angularMomentumVectorGrowthRate</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_
  !$omp threadprivate(darkMatterProfileDMO_)

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_Thread_Initialize(parameters_)
    !!{
    Initializes the halo vector angular momentum module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent  , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultSpinComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_
    
    if (defaultSpinComponent%vectorIsActive()) then
       !![
       <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters_"/>
       !!]
       call nodePromotionEvent%attach(defaultSpinComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSpinVector')
    end if
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_Thread_Uninitialize()
    !!{
    Uninitializes the tree node preset spin module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (defaultSpinComponent%vectorIsActive()) then
       !![
       <objectDestructor name="darkMatterProfileDMO_"/>
       !!]
       if (nodePromotionEvent%isAttached(defaultSpinComponent,nodePromotion)) call nodePromotionEvent%detach(defaultSpinComponent,nodePromotion)
    end if
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_Thread_Uninitialize
 
  !![
  <scaleSetTask>
   <unitName>Node_Component_Halo_Angular_Momentum_Vector_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Halo_Angular_Momentum_Vector_Scale_Set(node)
    !!{
    Set scales for properties in the preset implementation of the spin component.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentSpin                      , nodeComponentSpinVector, treeNode, defaultSpinComponent
    implicit none
    type            (treeNode         ), intent(inout), pointer   :: node
    class           (nodeComponentSpin)               , pointer   :: spin
    double precision                                  , parameter :: spinMinimum=1.0d-6

    ! Return immediately if this class is not in use.
    if (.not.defaultSpinComponent%vectorIsActive()) return
    ! Get the spin component.
    spin => node%spin()
    ! Ensure that it is of the preset class.
    select type (spin)
    class is (nodeComponentSpinVector)
       ! Set scale for spin.
       call spin%angularMomentumVectorScale([1.0d0,1.0d0,1.0d0]*max(spin%angularMomentum(),spinMinimum*Dark_Matter_Halo_Angular_Momentum_Scale(node,darkMatterProfileDMO_)))
    end select
    return
  end subroutine Node_Component_Halo_Angular_Momentum_Vector_Scale_Set

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the angular momentum of {\normalfont \ttfamily node}
    to be that of its parent.
    !!}
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
    ! Adjust the angular momentum growth rate to that of the parent node.
    spinParent => nodeParent%spin()
    call spin%angularMomentumVectorSet          (spinParent%angularMomentumVector          ())
    call spin%angularMomentumVectorGrowthRateSet(spinParent%angularMomentumVectorGrowthRate())
    return
  end subroutine nodePromotion

end module Node_Component_Halo_Angular_Momentum_Vector
