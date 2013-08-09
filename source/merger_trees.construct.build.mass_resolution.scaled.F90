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

!% Contains a module which implements a scaled mass resolution for building merger trees.

module Merger_Trees_Build_Mass_Resolution_Scaled
  !% Implements a scaled mass resolution for building merger trees.
  implicit none
  private
  public :: Merger_Trees_Build_Mass_Resolution_Scaled_Initialize

  ! The scaled mass resolution.
  double precision :: mergerTreeBuildMassResolutionScaledMinimum,mergerTreeBuildMassResolutionScaledFraction

contains
  
  !# <mergerTreesBuildMassResolution>
  !#  <unitName>Merger_Trees_Build_Mass_Resolution_Scaled_Initialize</unitName>
  !# </mergerTreesBuildMassResolution>
  subroutine Merger_Trees_Build_Mass_Resolution_Scaled_Initialize(mergerTreesBuildMassResolutionMethod,Merger_Tree_Build_Mass_Resolution_Get)
    !% Initialize the modified Press-Schechter branching routines.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type     (varying_string  ),          intent(in   ) :: mergerTreesBuildMassResolutionMethod
    procedure(double precision), pointer, intent(inout) :: Merger_Tree_Build_Mass_Resolution_Get
    
    if (mergerTreesBuildMassResolutionMethod == 'scaled') then
       Merger_Tree_Build_Mass_Resolution_Get => Merger_Tree_Build_Mass_Resolution_Scaled
       !@ <inputParameter>
       !@   <name>mergerTreeBuildMassResolutionScaledMinimum</name>
       !@   <defaultValue>$5\times 10^9$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum mass (in units of $M_\odot$) of halos to be resolved in merger trees that are built using the ``scaled'' mass resolution algorithm.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildMassResolutionScaledMinimum',mergerTreeBuildMassResolutionScaledMinimum,defaultValue=5.0d9)
       !@ <inputParameter>
       !@   <name>mergerTreeBuildMassResolutionScaledFraction</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The fraction of the tree's root node mass to be used for the mass resolution in merger trees that are built using the ``scaled'' mass resolution algorithm.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeBuildMassResolutionScaledFraction',mergerTreeBuildMassResolutionScaledFraction,defaultValue=1.0d-3)
    end if
    return
  end subroutine Merger_Trees_Build_Mass_Resolution_Scaled_Initialize

  double precision function Merger_Tree_Build_Mass_Resolution_Scaled(thisTree)
    !% Returns a scaled mass resolution to use when building merger trees.
    use Galacticus_Nodes
    use Galacticus_Nodes
    implicit none
    type (mergerTree        ), intent(in   ) :: thisTree
    class(nodeComponentBasic), pointer       :: baseNodeBasic

    ! Get the basic component of the tree's base node.
    baseNodeBasic => thisTree%baseNode%basic()
    ! Compute the mass resolution.
    Merger_Tree_Build_Mass_Resolution_Scaled=                                    &
         & max(                                                                  &
         &     mergerTreeBuildMassResolutionScaledMinimum                      , &
         &     mergerTreeBuildMassResolutionScaledFraction*baseNodeBasic%mass()  &
         &    )
    return
  end function Merger_Tree_Build_Mass_Resolution_Scaled
  
end module Merger_Trees_Build_Mass_Resolution_Scaled
