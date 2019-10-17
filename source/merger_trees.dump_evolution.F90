!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module which dumps the evolution of merger trees to XML.

module Merger_Trees_Dump_Evolution
  !% Dumps the structure of entire merger trees.
  use :: ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Dump_Evolution, Merger_Tree_Dump_Evolution_Close

  ! Flag indicating if module is initialized.
  logical                 :: moduleInitialized=.false.

  ! Flag indicating if output is required.
  logical                 :: mergerTreeEvolutionDump

  ! File to which tree evolution should be dumped.
  type   (varying_string) :: mergerTreeEvolutionDumpFileName

  ! File unit for output.
  integer                 :: mergerTreeEvolutionDumpFileUnit

contains

  !# <postEvolveTask>
  !# <unitName>Merger_Tree_Dump_Evolution</unitName>
  !# </postEvolveTask>
  subroutine Merger_Tree_Dump_Evolution(thisNode)
    !% Trim histories attached to the disk.
    use :: Galacticus_Nodes, only : treeNode
    use :: Input_Parameters, only : globalParameters, inputParameter
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (.not.moduleInitialized) then
       !$omp critical (Merger_Tree_Dump_Evolution)
       if (.not.moduleInitialized) then
          !# <inputParameter>
          !#   <name>mergerTreeEvolutionDump</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not to output the evolution of merger trees.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>mergerTreeEvolutionDumpFileName</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>var_str('mergerTreeEvolution.xml')</defaultValue>
          !#   <description>Specifies the file to which merger tree evolution should be dumped.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>text</type>
          !# </inputParameter>
          ! Open a file for output.
          if (mergerTreeEvolutionDump)  then
             open(newUnit=mergerTreeEvolutionDumpFileUnit,file=char(mergerTreeEvolutionDumpFileName),status='unknown',form='formatted')
             write (mergerTreeEvolutionDumpFileUnit,'(a)') '<evolution>'
          end if
          ! Record that the module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Dump_Evolution)
    end if
    ! Dump the node.
    if (mergerTreeEvolutionDump.and.thisNode%isOnMainBranch()) call thisNode%serializeXML(mergerTreeEvolutionDumpFileUnit)
    return
  end subroutine Merger_Tree_Dump_Evolution

  !# <hdfPreCloseTask>
  !# <unitName>Merger_Tree_Dump_Evolution_Close</unitName>
  !# </hdfPreCloseTask>
  subroutine Merger_Tree_Dump_Evolution_Close()
    !% Close the merger tree evolution dump file.
    implicit none

    if (mergerTreeEvolutionDump) then
       write (mergerTreeEvolutionDumpFileUnit,'(a)') '</evolution>'
       close(mergerTreeEvolutionDumpFileUnit)
    end if
    return
  end subroutine Merger_Tree_Dump_Evolution_Close

end module Merger_Trees_Dump_Evolution
