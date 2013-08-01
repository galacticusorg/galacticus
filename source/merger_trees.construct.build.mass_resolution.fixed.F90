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

!% Contains a module which implements a fixed mass resolution for building merger trees.

module Merger_Trees_Build_Mass_Resolution_Fixed
  !% Implements a fixed mass resolution for building merger trees.
  implicit none
  private
  public :: Merger_Trees_Build_Mass_Resolution_Fixed_Initialize

  ! The fixed mass resolution.
  double precision :: mergerTreeBuildMassResolutionFixed

contains
  
  !# <mergerTreesBuildMassResolution>
  !#  <unitName>Merger_Trees_Build_Mass_Resolution_Fixed_Initialize</unitName>
  !# </mergerTreesBuildMassResolution>
  subroutine Merger_Trees_Build_Mass_Resolution_Fixed_Initialize(mergerTreesBuildMassResolutionMethod,Merger_Tree_Build_Mass_Resolution_Get)
    !% Initialize the modified Press-Schechter branching routines.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type     (varying_string  ),          intent(in   ) :: mergerTreesBuildMassResolutionMethod
    procedure(double precision), pointer, intent(inout) :: Merger_Tree_Build_Mass_Resolution_Get
    
    if (mergerTreesBuildMassResolutionMethod == 'fixed') then
       Merger_Tree_Build_Mass_Resolution_Get => Merger_Tree_Build_Mass_Resolution_Fixed
       !@ <inputParameter>
       !@   <name>mergerTreeBuildMassResolutionFixed</name>
       !@   <defaultValue>$5\times 10^9$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum mass (in units of $M_\odot$) of halos to be resolved in merger trees that are built.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildMassResolutionFixed',mergerTreeBuildMassResolutionFixed,defaultValue=5.0d9)
    end if
    return
  end subroutine Merger_Trees_Build_Mass_Resolution_Fixed_Initialize

  double precision function Merger_Tree_Build_Mass_Resolution_Fixed(thisTree)
    !% Returns a fixed mass resolution to use when building merger trees.
    use Merger_Trees
    implicit none
    type(mergerTree), intent(in   ) :: thisTree

    Merger_Tree_Build_Mass_Resolution_Fixed=mergerTreeBuildMassResolutionFixed
    return
  end function Merger_Tree_Build_Mass_Resolution_Fixed
  
end module Merger_Trees_Build_Mass_Resolution_Fixed
