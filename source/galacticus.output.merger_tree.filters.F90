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


!% Contains a module which provides filtering of output.

module Galacticus_Merger_Tree_Output_Filters
  !% Provides filtering of output.
  use ISO_Varying_String
  private
  public :: Galacticus_Merger_Tree_Output_Filter

  ! Flag indicating whether filters have been initialized.
  logical                                         :: filtersInitialized=.false.

  ! Count of the number of filters being applied.
  integer                                         :: filterCount       =0

  ! List of filters to be used.
  type(varying_string), allocatable, dimension(:) :: mergerTreeOutputFilters

contains

  logical function Galacticus_Merger_Tree_Output_Filter(thisNode)
    !% Return true if {\tt thisNode} should be included in the output. Always arbitrary filters to block output of {\tt thisNode}.
    use Input_Parameters
    use Tree_Nodes
    use Memory_Management
    !# <include directive="mergerTreeOutputFilter" type="moduleUse">
    include 'galacticus.output.merger_tree.filters.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

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
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeOutputFilters',mergerTreeOutputFilters)
       end if
       ! Flag that filters are now initialized.
       filtersInitialized=.true.
    end if

    ! Assume galaxy will be output by default.
    Galacticus_Merger_Tree_Output_Filter=.true.
    ! Return immediately if no filters were defined.
    if (filterCount == 0) return
    !# <include directive="mergerTreeOutputFilter" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,mergerTreeOutputFilters,Galacticus_Merger_Tree_Output_Filter</subroutineArgs>
    include 'galacticus.output.merger_tree.filters.inc'
    !# </include>
    return
  end function Galacticus_Merger_Tree_Output_Filter

end module Galacticus_Merger_Tree_Output_Filters
