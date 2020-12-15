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

!% Contains a module with the standard implementation of basic tree node methods.

module Node_Component_Basic_Non_Evolving
  !% A non-evolving implementation of basic tree node methods.
  implicit none
  private
  public :: Node_Component_Basic_Non_Evolving_Rate_Compute     , Node_Component_Basic_Non_Evolving_Scale_Set          , &
       &    Node_Component_Basic_Non_Evolving_Thread_Initialize, Node_Component_Basic_Non_Evolving_Thread_Uninitialize

  !# <component>
  !#  <class>basic</class>
  !#  <name>nonEvolving</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>mass</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="massSolar" comment="Total mass of the node, assuming univeral baryon fraction."/>
  !#   </property>
  !#   <property>
  !#     <name>time</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>timeLastIsolated</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <getFunction bindsTo="component">BasicNonEvolvingTimeLastIsolated</getFunction>
  !#     <output unitsInSI="gigaYear" comment="Time at which node was last an isolated halo."/>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.basic.non_evolving.bound_functions.inc</functions>
  !# </component>

contains

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Basic_Non_Evolving_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Basic_Non_Evolving_Thread_Initialize(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent   , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultBasicComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultBasicComponent%nonEvolvingIsActive()) &
         call nodePromotionEvent%attach(defaultBasicComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentBasicNonEvolving')
    return
  end subroutine Node_Component_Basic_Non_Evolving_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Basic_Non_Evolving_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Basic_Non_Evolving_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none

    if (defaultBasicComponent%nonEvolvingIsActive()) &
         & call nodePromotionEvent%detach(defaultBasicComponent,nodePromotion)
    return
  end subroutine Node_Component_Basic_Non_Evolving_Thread_Uninitialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Basic_Non_Evolving_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Basic_Non_Evolving_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Compute rates of change of properties in the standard implementation of the basic component.
    use :: Galacticus_Nodes, only : defaultBasicComponent, nodeComponentBasic, nodeComponentBasicNonEvolving, propertyTypeInactive, &
          &                         treeNode
    implicit none
    type     (treeNode          ), intent(inout)          :: node
    logical                      , intent(inout)          :: interrupt
    procedure(                  ), intent(inout), pointer :: interruptProcedure
    integer                      , intent(in   )          :: propertyType
    class    (nodeComponentBasic)               , pointer :: basicComponent
    !$GLC attributes unused :: interrupt, interruptProcedure


    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultBasicComponent%nonEvolvingIsActive()) return
    ! Get the basic component.
    basicComponent => node%basic()
    ! Ensure that it is of the non-evolving class.
    select type (basicComponent)
    class is (nodeComponentBasicNonEvolving)
       ! Time rate of change is unity, by definition.
       call basicComponent%timeRate(1.0d0)
    end select
    return
  end subroutine Node_Component_Basic_Non_Evolving_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Basic_Non_Evolving_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Basic_Non_Evolving_Scale_Set(node)
    !% Set scales for properties in the standard implementation of the basic component.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentBasicNonEvolving, treeNode
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    double precision                    , parameter              :: timeScale     =1.0d-3
    class           (nodeComponentBasic)               , pointer :: basicComponent

    ! Get the basic component.
    basicComponent => node%basic()
    ! Ensure that it is of the standard class.
    select type (basicComponent)
    class is (nodeComponentBasicNonEvolving)
       ! Set scale for time.
       call basicComponent%timeScale(timeScale)
    end select
    return
  end subroutine Node_Component_Basic_Non_Evolving_Scale_Set

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the mass of {\normalfont \ttfamily node}
    !% to be that of its parent.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , nodeComponentBasicNonEvolving, treeNode
    implicit none
    class(*                 ), intent(inout)          :: self
    type (treeNode          ), intent(inout), target  :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentBasic)               , pointer :: basicParent, basic
    !$GLC attributes unused :: self

    basic       => node      %basic ()
    nodeParent  => node      %parent
    basicParent => nodeParent%basic ()
    ! Ensure the two halos exist at the same time.
    if (basic%time() /= basicParent%time()) call Galacticus_Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust the mass to that of the parent node.
    call basic%massSet(basicParent%mass())
    return
  end subroutine nodePromotion

end module Node_Component_Basic_Non_Evolving
