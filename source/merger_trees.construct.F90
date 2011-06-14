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


!% Contains a module which constructs/destructs merger trees.

module Merger_Tree_Construction
  !% Constructs/destructs merger trees.
  use ISO_Varying_String
  use Merger_Trees
  private
  public :: Merger_Tree_Create

  ! Flag to indicate if this module has been initialized.  
  logical              :: mergerTreeConstructInitialized=.false.

  ! Name of tree construction method used.
  type(varying_string) :: mergerTreeConstructMethod
  ! Pointer to the subroutine that tabulates the transfer function and template interface for that subroutine.
  procedure(Merger_Tree_Construct_Template), pointer :: Merger_Tree_Construct => null()
  abstract interface
     subroutine Merger_Tree_Construct_Template(thisTree,skipTree)
       import mergerTree
       type(mergerTree), intent(inout) :: thisTree
       logical,          intent(in)    :: skipTree
     end subroutine Merger_Tree_Construct_Template
  end interface
  
contains

  function Merger_Tree_Create(skipTree) result(thisTree)
    !% Creates a merger tree.
    use Input_Parameters
    use Galacticus_Error
    use Memory_Management
    !# <include directive="mergerTreeConstructMethod" type="moduleUse">
    include 'merger_trees.construct.modules.inc'
    !# </include>
    implicit none
    type(mergerTree), pointer    :: thisTree
    logical,          intent(in) :: skipTree

    !$omp critical(Merger_Tree_Construct_Initialization) 
    ! Initialize if necessary.
    if (.not.mergerTreeConstructInitialized) then
       !@ <inputParameter>
       !@   <name>mergerTreeConstructMethod</name>
       !@   <defaultValue>build</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Selects the method to be used constructing merger trees.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeConstructMethod',mergerTreeConstructMethod,defaultValue='build')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="mergerTreeConstructMethod" type="code" action="subroutine">
       !#  <subroutineArgs>mergerTreeConstructMethod,Merger_Tree_Construct</subroutineArgs>
       include 'merger_trees.construct.inc'
       !# </include>
       if (.not.associated(Merger_Tree_Construct)) call Galacticus_Error_Report('Merger_Tree_Create','method '&
            &//char(mergerTreeConstructMethod)//' is unrecognized')
       mergerTreeConstructInitialized=.true.
    end if
    !$omp end critical(Merger_Tree_Construct_Initialization)

    ! Create the object.
    allocate(thisTree)
    call Memory_Usage_Record(sizeof(thisTree))
    thisTree%baseNode => null()

    ! Flag that the tree is uninitialized. Some construction methods may opt to fully initialize the tree, in which case they will
    ! reset this to true.
    thisTree%initialized=.false.

    ! Call the routine to construct the merger tree.
    call Merger_Tree_Construct(thisTree,skipTree)

    ! Deallocate the tree if no nodes were created.
    if (.not.associated(thisTree%baseNode)) then
       call Memory_Usage_Record(sizeof(thisTree),addRemove=-1)
       deallocate(thisTree)
       thisTree => null()
    end if
    return
  end function Merger_Tree_Create

end module Merger_Tree_Construction
