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

!% Contains a module of halo formation time methods.

module Node_Component_Formation_Times_Mass_Fraction
  !% Implement tracking of halo formation times.
  implicit none
  private
  public :: Node_Component_Formation_Times_Mass_Fraction_Tree_Initialize, Node_Component_Formation_Times_Mass_Fraction_Initialize   , &
       &    Node_Component_Formation_Times_Mass_Fraction_Thread_Init    , Node_Component_Formation_Times_Mass_Fraction_Thread_Uninit

  !# <component>
  !#  <class>formationTime</class>
  !#  <name>massFraction</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>formationTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="gigaYear" comment="The time at which a fixed fraction of the node''s mass was assembled."/>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Fractional mass of primary progenitor used to define formation time.
  double precision :: formationTimeMassFraction

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Formation_Times_Mass_Fraction_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Formation_Times_Mass_Fraction_Initialize(parameters_)
    !% Initializes the tree node formation time tracking module.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    !# <inputParameter>
    !#   <name>formationTimeMassFraction</name>
    !#   <defaultValue>0.5d0</defaultValue>
    !#   <source>parameters_</source>
    !#   <description>Fractional mass of primary progenitor used to define formation time.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    return
  end subroutine Node_Component_Formation_Times_Mass_Fraction_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Formation_Times_Mass_Fraction_Thread_Init</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Formation_Times_Mass_Fraction_Thread_Init(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent       , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultFormationTimeComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultFormationTimeComponent%massFractionIsActive()) &
         call nodePromotionEvent%attach(defaultFormationTimeComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentFormationTimeMassFraction')
    return
  end subroutine Node_Component_Formation_Times_Mass_Fraction_Thread_Init

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Formation_Times_Mass_Fraction_Thread_Uninit</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Formation_Times_Mass_Fraction_Thread_Uninit()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultFormationTimeComponent
    implicit none

    if (defaultFormationTimeComponent%massFractionIsActive()) &
         & call nodePromotionEvent%detach(defaultFormationTimeComponent,nodePromotion)
    return
  end subroutine Node_Component_Formation_Times_Mass_Fraction_Thread_Uninit

  subroutine nodePromotion(self,node)
    !% Handle node promotion for formation times.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentFormationTime, nodeComponentFormationTimeMassFraction, treeNode
    implicit none
    class(*                         ), intent(inout)          :: self
    type (treeNode                  ), intent(inout), target  :: node
    class(nodeComponentFormationTime)               , pointer :: formationTime, formationTimeParent
    type (treeNode                  )               , pointer :: nodeParent
    !$GLC attributes unused :: self
    
    formationTime       => node      %formationTime()
    nodeParent          => node      %parent
    formationTimeParent => nodeParent%formationTime()
    ! Adjust the formation time to that of the parent node.
    call formationTime%formationTimeSet(formationTimeParent%formationTime())
    return
  end subroutine nodePromotion

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Formation_Times_Mass_Fraction_Tree_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Formation_Times_Mass_Fraction_Tree_Initialize(node)
    !% Initialize the formation node pointer for any childless node.
    use :: Galacticus_Nodes, only : defaultFormationTimeComponent, nodeComponentBasic, nodeComponentFormationTime, treeNode
    implicit none
    type            (treeNode                  ), intent(inout), pointer :: node
    type            (treeNode                  )               , pointer :: nodeProgenitor
    class           (nodeComponentBasic        )               , pointer :: basic         , basicProgenitor
    class           (nodeComponentFormationTime)               , pointer :: formationTime
    double precision                                                     :: timeFormation , massFormation

    ! Return immediately if this implementation is not active.
    if (.not.defaultFormationTimeComponent%massFractionIsActive()) return
    ! Compute the formation time.
    formationTime => node%formationTime(autoCreate=.true.)
    ! Determine the formation mass.
    basic         =>  node %basic()
    massFormation =  +basic%mass ()             &
         &           *formationTimeMassFraction
    ! Step through progenitors to find the formation time.
    timeFormation  =  -1.0d0 ! Assume formation term is indeterminate by default.
    nodeProgenitor => node
    do while (associated(nodeProgenitor%firstChild))
       basic           => nodeProgenitor           %basic()
       basicProgenitor => nodeProgenitor%firstChild%basic()
       if     (                                         &
            &   basic          %mass() >  massFormation &
            &  .and.                                    &
            &   basicProgenitor%mass() <= massFormation &
            & ) then
          ! Interpolate to find the formation time.
          timeFormation=+(basic%time         ()-basicProgenitor%time()) &
               &        /(basic%mass         ()-basicProgenitor%mass()) &
               &        *(      massFormation  -basicProgenitor%mass()) &
               &                               +basicProgenitor%time()
          exit
       end if
       nodeProgenitor => nodeProgenitor%firstChild
    end do
    call formationTime%formationTimeSet(timeFormation)
    return
  end subroutine Node_Component_Formation_Times_Mass_Fraction_Tree_Initialize

end module Node_Component_Formation_Times_Mass_Fraction
