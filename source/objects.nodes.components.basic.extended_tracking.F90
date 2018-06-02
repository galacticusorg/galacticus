!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which extends the extended implementation of basic component to track the maximum progenitor mass.

module Node_Component_Basic_Extended_Tracking
  !% Extends the extended implementation of basic component to track the maximum progenitor mass.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Basic_Extended_Tree_Tracking_Initialize, Node_Component_Basic_Extended_Tracking_Promote

  !# <component>
  !#  <class>basic</class>
  !#  <name>extendedTracking</name>
  !#  <extends>
  !#   <class>basic</class>
  !#   <name>standardExtended</name>
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
  !#  <unitName>Node_Component_Basic_Extended_Tree_Tracking_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Basic_Extended_Tree_Tracking_Initialize(node)
    !% Set the mass accretion rate for {\normalfont \ttfamily node}.
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    type            (treeNode          )               , pointer :: nodeProgenitor
    class           (nodeComponentBasic)               , pointer :: basic         , basicProgenitor
    double precision                                             :: massMaximum

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the extended tracking class.
    select type (basic)
    class is (nodeComponentBasicExtendedTracking)
       ! Find the maximum progenitor mass.
       nodeProgenitor => node
       massMaximum    =  0.0d0
       do while (associated(nodeProgenitor))
          basicProgenitor => nodeProgenitor%basic()
          massMaximum     =  max(massMaximum,basicProgenitor%mass())
          nodeProgenitor  => nodeProgenitor%firstChild
       end do
       call basic%massMaximumSet(massMaximum)
    end select
    return
  end subroutine Node_Component_Basic_Extended_Tree_Tracking_Initialize

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Basic_Extended_Tracking_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Basic_Extended_Tracking_Promote(node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the maximum
    !% mass of {\normalfont \ttfamily node} to be that of its parent.
    use Galacticus_Error
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    type (treeNode          )               , pointer :: nodeParent
    class(nodeComponentBasic)               , pointer :: basicParent, basic

    ! Get the basic component.
    basic => node%basic()
    ! Ensure that it is of the extended tracking class.
    select type (basic)
    class is (nodeComponentBasicExtendedTracking)
       ! Get the parent node and its basic component.
       nodeParent  => node      %parent
       basicParent => nodeParent%basic ()
       ! Adjust the maximum mass to that of the parent node.
       call basic%massMaximumSet(basicParent%massMaximum())
    end select
    return
  end subroutine Node_Component_Basic_Extended_Tracking_Promote

end module Node_Component_Basic_Extended_Tracking
