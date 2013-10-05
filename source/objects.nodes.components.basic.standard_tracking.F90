!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Basic_Standard_Tree_Tracking_Initialize, Node_Component_Basic_Standard_Tracking_Promote

  !# <component>
  !#  <class>basic</class>
  !#  <name>standardTracking</name>
  !#  <extends>
  !#   <class>basic</class>
  !#   <name>standard</name>
  !#  </extends>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>massMaximum</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Basic_Standard_Tree_Tracking_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Basic_Standard_Tree_Tracking_Initialize(thisNode)
    !% Set the mass accretion rate for {\tt thisNode}.
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    type            (treeNode          )               , pointer :: progenitorNode
    class           (nodeComponentBasic)               , pointer :: thisBasic     , progenitorBasic
    double precision                                             :: massMaximum

    ! Get the basic component.
    thisBasic => thisNode%basic()
    ! Ensure that it is of the standard tracking class.
    select type (thisBasic)
       class is (nodeComponentBasicStandardTracking)
          ! Find the maximum progenitor mass.
       progenitorNode => thisNode
       massMaximum    =  0.0d0
       do while (associated(progenitorNode))
          progenitorBasic => progenitorNode%basic     ()
          massMaximum     =  max(massMaximum,progenitorBasic%mass())
          progenitorNode  => progenitorNode%firstChild
       end do
       call thisBasic%massMaximumSet(massMaximum)
    end select
    return
  end subroutine Node_Component_Basic_Standard_Tree_Tracking_Initialize

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Basic_Standard_Tracking_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Basic_Standard_Tracking_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the maximum mass of {\tt
    !% thisNode} to be that of its parent.
    use Galacticus_Error
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode
    type (treeNode          )               , pointer :: parentNode
    class(nodeComponentBasic)               , pointer :: parentBasic, thisBasic

    ! Get the basic component.
    thisBasic => thisNode%basic()
    ! Ensure that it is of the standard tracking class.
    select type (thisBasic)
       class is (nodeComponentBasicStandardTracking)
          ! Get the parent node and its basic component.
       parentNode  => thisNode  %parent
       parentBasic => parentNode%basic()
       ! Adjust the maximum mass to that of the parent node.
       call thisBasic%massMaximumSet(parentBasic%massMaximum())
    end select
    return
  end subroutine Node_Component_Basic_Standard_Tracking_Promote

end module Node_Component_Basic_Standard_Tracking
