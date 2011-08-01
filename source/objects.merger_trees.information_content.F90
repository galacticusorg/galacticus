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


!% Contains a module which computes the cladistic information content (CIC) of a merger tree.

module Merger_Trees_Information_Content
  !% Implements computation of the cladistic information content (CIC) of a merger tree.
  implicit none
  private
  public :: Merger_Tree_Cladistic_Information_Content

contains

  double precision function Merger_Tree_Cladistic_Information_Content(thisTree)
    !% Compute the cladistic information content (CIC; \citealt{thorley_information_1998}) of a merger tree in bits.
    use Factorials
    use Merger_Trees
    use Tree_Nodes
    implicit none
    type(mergerTree), intent(in) :: thisTree
    type(treeNode),   pointer    :: thisNode,childNode
    integer                      :: leafCount,childCount
    double precision             :: logPermittedBifurcations,logPossibleBifurcations  

    ! Walk the tree, counting the number of leaves and accumulated the log of the number of permitted bifurcations.
    leafCount               =0
    logPermittedBifurcations=0.0d0
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       ! Check for leaf nodes.
       if (.not.associated(thisNode%childNode)) then
          ! Increment leaf node counter.
          leafCount=leafCount+1
       else
          ! Count child nodes.
          childCount=0
          childNode => thisNode%childNode
          do while (associated(childNode))
             childCount=childCount+1
             childNode => childNode%siblingNode
          end do
          ! Increment the number of permitted bifurcations based on this number of children.
          logPermittedBifurcations=logPermittedBifurcations+Logarithmic_Double_Factorial(2*childCount-3)
       end if
       call thisNode%walkTree()
    end do

    ! Compute logarithm of the possible bifurcations.
    logPossibleBifurcations=Logarithmic_Double_Factorial(2*leafCount-3)
    
    ! Compute the CIC in bits.
    Merger_Tree_Cladistic_Information_Content=(logPossibleBifurcations-logPermittedBifurcations)/dlog(2.0d0)
    return
  end function Merger_Tree_Cladistic_Information_Content
  
end module Merger_Trees_Information_Content
