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

!% Contains a module which implements calculations of merger tree building mass resolutions.

module Merger_Trees_Build_Mass_Resolution
  !% Implements calculations of merger tree branching probabilities.
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Build_Mass_Resolution

  ! Flag to indicate if this module has been initialized.  
  logical                                               :: moduleInitialized=.false.

  ! Pointer to the function that returns the mass resolution.
  procedure(Merger_Tree_Build_Mass_Resolution), pointer :: Merger_Tree_Build_Mass_Resolution_Get => null()
 
contains

  subroutine Merger_Tree_Build_Mass_Resolution_Initialize
    !% Initializes the merger tree mass resolution module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="mergerTreesBuildMassResolution" type="moduleUse">
    include 'merger_trees.construct.build.mass_resolution.modules.inc'
    !# </include>
    implicit none
    type(varying_string) :: mergerTreesBuildMassResolutionMethod

    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Merger_Tree_Build_Mass_Resolution_Initialization) 
       if (.not.moduleInitialized) then
          ! Get the mass resolution method parameter.
          !@ <inputParameter>
          !@   <name>mergerTreesBuildMassResolutionMethod</name>
          !@   <defaultValue>fixed</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the mass resolution to use when building merger trees.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreesBuildMassResolutionMethod',mergerTreesBuildMassResolutionMethod,defaultValue='fixed')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="mergerTreesBuildMassResolution" type="functionCall" functionType="void">
          !#  <functionArgs>mergerTreesBuildMassResolutionMethod,Merger_Tree_Build_Mass_Resolution_Get</functionArgs>
          include 'merger_trees.construct.build.mass_resolution.inc'
          !# </include>
          if (.not.associated(Merger_Tree_Build_Mass_Resolution_Get)) call Galacticus_Error_Report('Merger_Tree_Build_Mass_Resolution_Initialize','method '//char(mergerTreesBuildMassResolutionMethod)//' is unrecognized')

          moduleInitialized=.true.
       end if
       !$omp end critical(Merger_Tree_Build_Mass_Resolution_Initialization)
    end if
    return
  end subroutine Merger_Tree_Build_Mass_Resolution_Initialize

  double precision function Merger_Tree_Build_Mass_Resolution(thisTree)
    !% Return the mass resolution to use when building {\tt thisTree}.
    use Merger_Trees
    implicit none
    type(mergerTree), intent(in   ) :: thisTree

    ! Initialize if necessary.
    call Merger_Tree_Build_Mass_Resolution_Initialize

    ! Call the function to compute the mass resolution
    Merger_Tree_Build_Mass_Resolution=Merger_Tree_Build_Mass_Resolution_Get(thisTree)
    return
  end function Merger_Tree_Build_Mass_Resolution
  
end module Merger_Trees_Build_Mass_Resolution
