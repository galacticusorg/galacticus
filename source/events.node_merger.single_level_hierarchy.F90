!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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

    if (nodeMergersMethod == 'single level hierarchy') Events_Node_Merger_Do => Events_Node_Merger_Do_SLH

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
