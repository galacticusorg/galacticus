!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which dumps the structure of entire merger trees.

module Merger_Tree_Dump_Structure
  !% Dumps the structure of entire merger trees.
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Structure_Dump

  ! Flag indicating if module is initialized.
  logical                          :: moduleInitialized=.false.

  ! Flag indicating if output is required.
  logical                          :: mergerTreeStructureDump

  ! Directory to which tree structure should be dumped.
  type            (varying_string) :: mergerTreeStructureDumpDirectory

  ! Halo mass range to dump.
  double precision                 :: mergerTreeStructureDumpMassMinimum,mergerTreeStructureDumpMassMaximum

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Structure_Dump</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Structure_Dump(thisTree)
    !% Output the structure of {\tt thisTree}.
    use Galacticus_Nodes
    use Input_Parameters
    use Merger_Trees_Dump
    implicit none
    type (mergerTree        ), intent(in   ),  target :: thisTree
    type (mergerTree        ), pointer                :: currentTree
    class(nodeComponentBasic), pointer                :: baseNodeBasic

    ! Check if module is initialized.
    if (.not.moduleInitialized) then
       !$omp critical(structureDumpModuleInitialize)
       if (.not.moduleInitialized) then
          ! Get parameter specifying if output is required.
          !@ <inputParameter>
          !@   <name>mergerTreeStructureDump</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to output the structure of merger trees prior to evolution.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeStructureDump',mergerTreeStructureDump,defaultValue=.false.)
          !@ <inputParameter>
          !@   <name>mergerTreeStructureDumpDirectory</name>
          !@   <defaultValue>{\tt .}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies the directory to which merger tree structure should be dumped.
          !@   </description>
          !@   <type>text</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeStructureDumpDirectory',mergerTreeStructureDumpDirectory,defaultValue=".")
          !@ <inputParameter>
          !@   <name>mergerTreeStructureDumpMassMinimum</name>
          !@   <defaultValue>$0M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies the minimum root mass for which merger tree structure should be dumped.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeStructureDumpMassMinimum',mergerTreeStructureDumpMassMinimum,defaultValue=0.0d0)
          !@ <inputParameter>
          !@   <name>mergerTreeStructureDumpMassMaximum</name>
          !@   <defaultValue>$10^{30}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies the minimum root mass for which merger tree structure should be dumped.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeStructureDumpMassMaximum',mergerTreeStructureDumpMassMaximum,defaultValue=1.0d30)
          ! Flag that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(structureDumpModuleInitialize)
    end if

    ! Output the tree structure history.
    if (mergerTreeStructureDump) then
       ! Iterate over trees.
       currentTree => thisTree
       do while (associated(currentTree))
          ! Dump the tree.
          baseNodeBasic => currentTree%baseNode%basic()
          if     (                                                                           &
               &   baseNodeBasic%mass() >= mergerTreeStructureDumpMassMinimum                &
               &  .and.                                                                      &
               &   baseNodeBasic%mass() < mergerTreeStructureDumpMassMaximum                 &
               & )                                                                           &
               & call Merger_Tree_Dump(                                                      &
               &                       currentTree%index                                   , &
               &                       currentTree%baseNode                                , &
               &                       scaleNodesByLogMass=.true.                          , &
               &                       edgeLengthsToTimes =.true.                          , &
               &                       path               =mergerTreeStructureDumpDirectory  &
               &                      )
          ! Move to the next tree.
          currentTree => currentTree%nextTree
       end do
    end if

    return
  end subroutine Merger_Tree_Structure_Dump

end module Merger_Tree_Dump_Structure
