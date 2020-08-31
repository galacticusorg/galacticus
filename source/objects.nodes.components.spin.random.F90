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

!% Contains a module of spin tree node methods.

module Node_Component_Spin_Random
  !% Implement random spin tree node method.
  use :: Halo_Spin_Distributions, only : haloSpinDistributionClass
  implicit none
  private
  public :: Node_Component_Spin_Random_Initialize       , Node_Component_Spin_Random_Initialize_Spins   , &
       &    Node_Component_Spin_Random_Thread_Initialize, Node_Component_Spin_Random_Thread_Uninitialize

  !# <component>
  !#  <class>spin</class>
  !#  <name>random</name>
  !#  <isDefault>true</isDefault>
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
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <getFunction bindsTo="component">SpinRandomSpinGrowthRate</getFunction>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.spin.random.bound_functions.inc</functions>
  !# </component>

  ! Objects used by this component.
  class(haloSpinDistributionClass), pointer:: haloSpinDistribution_
  !$omp threadprivate(haloSpinDistribution_)

  ! The factor by which the mass of a node must increase before its spin parameter is re-chosen.
  double precision :: randomSpinResetMassFactor

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Spin_Random_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Spin_Random_Initialize(parameters_)
    !% Initializes the random spin component module.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    !# <inputParameter>
    !#   <name>randomSpinResetMassFactor</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>2.0d0</defaultValue>
    !#   <description>The factor by which a node must increase in mass before its spin parameter is reset.</description>
    !#   <source>parameters_</source>
    !#   <type>real</type>
    !# </inputParameter>
    return
  end subroutine Node_Component_Spin_Random_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Spin_Random_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Spin_Random_Thread_Initialize(parameters_)
    !% Initializes the tree node random spin module.
    use :: Events_Hooks    , only : nodePromotionEvent  , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultSpinComponent
    use :: Input_Parameters, only : inputParameter      , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultSpinComponent%randomIsActive()) then
       !# <objectBuilder class="haloSpinDistribution" name="haloSpinDistribution_" source="parameters_"/>
       call nodePromotionEvent%attach(defaultSpinComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSpinRandom')
    end if
    return
  end subroutine Node_Component_Spin_Random_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Spin_Random_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Spin_Random_Thread_Uninitialize()
    !% Uninitializes the tree node random spin module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none

    if (defaultSpinComponent%randomIsActive()) then
       !# <objectDestructor name="haloSpinDistribution_"/>
       call nodePromotionEvent%detach(defaultSpinComponent,nodePromotion)
    end if
    return
  end subroutine Node_Component_Spin_Random_Thread_Uninitialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Spin_Random_Initialize_Spins</unitName>
  !#  <sortName>spin</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Spin_Random_Initialize_Spins(node)
    !% Initialize the spin of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : defaultSpinComponent, nodeComponentBasic, nodeComponentSpin, treeNode
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    type            (treeNode          )               , pointer :: nodeRelated
    class           (nodeComponentSpin )               , pointer :: spinRelated    , spin
    class           (nodeComponentBasic)               , pointer :: basicRelated
    double precision                                             :: previousSetMass, previousSetSpin

    ! Check if we are the default method.
    if (defaultSpinComponent%randomIsActive()) then
       ! Get the basic component.
       spin => node%spin()
       ! Ensure that the spin has not yet been assigned for this node.
       select type (spin)
       type is (nodeComponentSpin)
           ! Walk the tree back along primary children to the earliest such progenitor.
          nodeRelated => node
          do while (associated(nodeRelated%firstChild))
             nodeRelated => nodeRelated%firstChild
          end do
          ! Walk forward through the branch, assigning spins. If the mass of the halo exceeds that of the halo for which we last
          ! selected a spin by a given factor, then select a new spin from the distribution. Otherwise, use the previously
          ! assigned spin.
          spinRelated     => nodeRelated          %spin  (autoCreate =.true.)
          basicRelated    => nodeRelated          %basic (                  )
          previousSetSpin =  haloSpinDistribution_%sample(nodeRelated       )
          previousSetMass =  basicRelated         %mass  (                  )
          call spinRelated%spinSet(previousSetSpin)
          do while (nodeRelated%isPrimaryProgenitor())
             nodeRelated           => nodeRelated%parent
             basicRelated => nodeRelated%basic()
             if (basicRelated%mass() > randomSpinResetMassFactor*previousSetMass) then
                previousSetSpin=haloSpinDistribution_%sample(nodeRelated)
                previousSetMass=basicRelated         %mass  (           )
             end if
             spinRelated => nodeRelated%spin(autoCreate=.true.)
             call spinRelated%spinSet(previousSetSpin)
          end do
       end select
    end if
    return
  end subroutine Node_Component_Spin_Random_Initialize_Spins
  
  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the spin of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentSpin, nodeComponentSpinRandom, treeNode
    implicit none
    class(*                 ), intent(inout) :: self
    type (treeNode          ), intent(inout) :: node
    type (treeNode          ), pointer       :: nodeParent
    class(nodeComponentSpin ), pointer       :: spinParent , spin
    class(nodeComponentBasic), pointer       :: basicParent, basic
    !$GLC attributes unused :: self

    ! Ensure that the spin component is of the random class.
    spin        => node     %spin  ()
    nodeParent  => node     %parent
    basic       => node     %basic ()
    basicParent => nodeParent%basic()
    if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust the spin to that of the parent node.
    spinParent => nodeParent%spin()
    call spin%spinSet(spinParent%spin())
    return
  end subroutine nodePromotion

end module Node_Component_Spin_Random
