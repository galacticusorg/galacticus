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

!% Contains a module which extends the standard implementation of basic component to track the maximum progenitor mass.

module Node_Component_Basic_Standard_Tracking
  !% Extends the standard implementation of basic component to track the maximum progenitor mass.
  implicit none
  private
  public :: Node_Component_Basic_Standard_Tree_Tracking_Initialize    , Node_Component_Basic_Standard_Tracking_Thread_Initialize, &
       &    Node_Component_Basic_Standard_Tracking_Thread_Uninitialize

  !# <component>
  !#  <class>basic</class>
  !#  <name>standardTracking</name>
  !#  <extends>
  !#   <class>basic</class>
  !#   <name>standard</name>
  !#  </extends>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>massMaximum</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Basic_Standard_Tree_Tracking_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Basic_Standard_Tree_Tracking_Initialize(node)
    !% Set the mass accretion rate for {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentBasicStandardTracking, treeNode
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    type            (treeNode          )               , pointer :: nodeProgenitor
    class           (nodeComponentBasic)               , pointer :: basic         , basicProgenitor
    double precision                                             :: massMaximum

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the standard tracking class.
    select type (basic)
       class is (nodeComponentBasicStandardTracking)
          ! Find the maximum progenitor mass.
       nodeProgenitor => node
       massMaximum    =  0.0d0
       do while (associated(nodeProgenitor))
          basicProgenitor => nodeProgenitor%basic     ()
          massMaximum     =  max(massMaximum,basicProgenitor%mass())
          nodeProgenitor  => nodeProgenitor%firstChild
       end do
       call basic%massMaximumSet(massMaximum)
    end select
    return
  end subroutine Node_Component_Basic_Standard_Tree_Tracking_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Basic_Standard_Tracking_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Basic_Standard_Tracking_Thread_Initialize(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent   , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultBasicComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultBasicComponent%standardTrackingIsActive()) &
         call nodePromotionEvent%attach(defaultBasicComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentBasicStandardTracking')
    return
  end subroutine Node_Component_Basic_Standard_Tracking_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Basic_Standard_Tracking_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Basic_Standard_Tracking_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none

    if (defaultBasicComponent%standardTrackingIsActive()) &
         & call nodePromotionEvent%detach(defaultBasicComponent,nodePromotion)
    return
  end subroutine Node_Component_Basic_Standard_Tracking_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the maximum mass of {\normalfont \ttfamily
    !% node} to be that of its parent.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic     , treeNode
    implicit none
    class(*                 ), intent(inout)          :: self
    type (treeNode          ), intent(inout), target  :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentBasic)               , pointer :: basicParent, basic
    !$GLC attributes unused :: self

    basic       => node      %basic ()
    nodeParent  => node      %parent
    basicParent => nodeParent%basic ()
    ! Adjust the maximum mass to that of the parent node.
    call basic%massMaximumSet(basicParent%massMaximum())
    return
  end subroutine nodePromotion

end module Node_Component_Basic_Standard_Tracking
