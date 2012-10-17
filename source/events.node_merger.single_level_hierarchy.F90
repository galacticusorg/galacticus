!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements merger of nodes utilizing a single level substructure hierarchy.

module Events_Node_Mergers_SLH
  !% Implements merger of nodes utilizing a single level substructure hierarchy.
  implicit none
  private
  public :: Events_Node_Merger_Initialize_SLH
  
contains

  !# <nodeMergersMethod>
  !#  <unitName>Events_Node_Merger_Initialize_SLH</unitName>
  !# </nodeMergersMethod>
  subroutine Events_Node_Merger_Initialize_SLH(nodeMergersMethod,Events_Node_Merger_Do)
    !% Determine if use of this method is requested and set procedure pointer appropriately if it is.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: nodeMergersMethod
    procedure(),          pointer, intent(inout) :: Events_Node_Merger_Do

    if (nodeMergersMethod == 'singleLevelHierarchy') Events_Node_Merger_Do => Events_Node_Merger_Do_SLH

    return
  end subroutine Events_Node_Merger_Initialize_SLH
  
  subroutine Events_Node_Merger_Do_SLH(thisNode)
    !% Processes a node merging event, utilizing a single level substructure hierarchy.
    use Tree_Nodes
    use Satellite_Promotion
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(treeNode),                pointer :: parentNode,childNode,thisSatellite,nextSatellite

    ! Get the parent node.
    parentNode => thisNode%parentNode

    ! Uncouple thisNode from the children of its parent.
    childNode => parentNode%childNode
    do while (.not.associated(childNode%siblingNode,thisNode))
       childNode => childNode%siblingNode
    end do
    childNode%siblingNode => thisNode%siblingNode

    ! Unset the sibling pointer for this node.
    thisNode%siblingNode => null()

    ! Add it to the list of satellite nodes associated with its parent.
    if (associated(parentNode%satelliteNode)) then
       call parentNode%lastSatellite(thisSatellite)
       thisSatellite%siblingNode => thisNode
    else
       parentNode%satelliteNode => thisNode
    end if

    ! Move any of its own satellites to become satellites of the parent and set their parent node pointers appropriately.
    if (associated(thisNode%satelliteNode)) then
       thisSatellite => thisNode%satelliteNode
       do while (associated(thisSatellite))
          ! Find next sibling satellite.
          nextSatellite => thisSatellite%siblingNode
          ! Move the satellite to the new parent.
          call Satellite_Move_To_New_Host(thisSatellite,parentNode)
          ! Move to the next sibling satellite.
          thisSatellite => nextSatellite
       end do
       thisNode%satelliteNode => null()
    end if
    return
  end subroutine Events_Node_Merger_Do_SLH

end module Events_Node_Mergers_SLH
