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

!% Contains a module which provides filtering of output.

module Galacticus_Merger_Tree_Output_Filters
  !% Provides filtering of output.
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output_Filter, Galacticus_Merger_Tree_Output_Filter_Initialize

  ! Flag indicating whether filters have been initialized.
  logical                                         :: filtersInitialized=.false.

  ! Count of the number of filters being applied.
  integer                                         :: filterCount       =0

  ! List of filters to be used.
  type(varying_string), allocatable, dimension(:) :: mergerTreeOutputFilters

contains

  subroutine Galacticus_Merger_Tree_Output_Filter_Initialize()
    !% Initialize the output filter subsystem.
    use Input_Parameters
    use Memory_Management
    !# <include directive="mergerTreeOutputFilterInitialize" type="moduleUse">
    include 'galacticus.output.merger_tree.filters.initialize.modules.inc'
    !# </include>
    implicit none

    ! Initialize filters if necessary.
    if (.not.filtersInitialized) then
       ! Determine how many filters are to be applied.
       filterCount=Get_Input_Parameter_Array_Size('mergerTreeOutputFilters')
       ! Allocate filter array and read filter names.
       if (filterCount > 0) then
          allocate(mergerTreeOutputFilters(filterCount))
          call Memory_Usage_Record(sizeof(mergerTreeOutputFilters))
          !@ <inputParameter>
          !@   <name>mergerTreeOutputFilters</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    A list of filters that should be applied when deciding which galaxies to output.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1..*</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeOutputFilters',mergerTreeOutputFilters)

          !# <include directive="mergerTreeOutputFilterInitialize" type="functionCall" functionType="void">
          !#  <functionArgs>mergerTreeOutputFilters</functionArgs>
          include 'galacticus.output.merger_tree.filters.initialize.inc'
          !# </include>

       end if

       ! Flag that filters are now initialized.
       filtersInitialized=.true.
    end if

    return
  end subroutine Galacticus_Merger_Tree_Output_Filter_Initialize

  logical function Galacticus_Merger_Tree_Output_Filter(thisNode)
    !% Return true if {\tt thisNode} should be included in the output. Always arbitrary filters to block output of {\tt thisNode}.
    use Galacticus_Nodes
    !# <include directive="mergerTreeOutputFilter" type="moduleUse">
    include 'galacticus.output.merger_tree.filters.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Ensure the filter subsystem is initialized.
    call Galacticus_Merger_Tree_Output_Filter_Initialize()

    ! Assume galaxy will be output by default.
    Galacticus_Merger_Tree_Output_Filter=.true.
    ! Return immediately if no filters were defined.
    if (filterCount == 0) return
    !# <include directive="mergerTreeOutputFilter" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode,Galacticus_Merger_Tree_Output_Filter</functionArgs>
    include 'galacticus.output.merger_tree.filters.inc'
    !# </include>
    return
  end function Galacticus_Merger_Tree_Output_Filter

end module Galacticus_Merger_Tree_Output_Filters
