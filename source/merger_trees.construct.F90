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

!% Contains a module which constructs/destructs merger trees.

module Merger_Tree_Construction
  !% Constructs/destructs merger trees.
  use ISO_Varying_String
  use Merger_Trees
  implicit none
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

    ! Initialize if necessary.
    if (.not.mergerTreeConstructInitialized) then
       !$omp critical(Merger_Tree_Construct_Initialization) 
       if (.not.mergerTreeConstructInitialized) then
          !@ <inputParameter>
          !@   <name>mergerTreeConstructMethod</name>
          !@   <defaultValue>build</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Selects the method to be used for constructing merger trees.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeConstructMethod',mergerTreeConstructMethod,defaultValue='build')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="mergerTreeConstructMethod" type="functionCall" functionType="void">
          !#  <functionArgs>mergerTreeConstructMethod,Merger_Tree_Construct</functionArgs>
          include 'merger_trees.construct.inc'
          !# </include>
          if (.not.associated(Merger_Tree_Construct)) call Galacticus_Error_Report('Merger_Tree_Create','method '&
               &//char(mergerTreeConstructMethod)//' is unrecognized')
          mergerTreeConstructInitialized=.true.
       end if
       !$omp end critical(Merger_Tree_Construct_Initialization)
    end if

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
