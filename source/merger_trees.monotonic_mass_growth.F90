!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which enforces monotonic mass growth along merger tree branches.

module Merger_Trees_Monotonic_Mass_Growth
  !% Enforces monotonic mass growth along merger tree branches.
  implicit none
  private
  public :: Merger_Tree_Monotonic_Mass_Growth
  
  ! Flag indicating if module is initialized.
  logical          :: monotonicGrowthModuleInitialized=.false.
 
  ! Flag indicating if pruning is required.
  logical          :: mergerTreeEnforceMonotonicGrowth
contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Monotonic_Mass_Growth</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Monotonic_Mass_Growth(thisTree)
    !% Enforce monotonic mass growth along branches of {\tt thisTree}.
    use Merger_Trees
    use Galacticus_Nodes
    use Input_Parameters
    implicit none
    type (mergerTree        ), intent(in) :: thisTree
    type (treeNode          ), pointer    :: thisNode,progenitorNode
    class(nodeComponentBasic), pointer    :: thisBasicComponent,progenitorBasicComponent
    logical                               :: didModifyTree
    double precision                      :: progenitorMass

    ! Check if module is initialized.
    if (.not.monotonicGrowthModuleInitialized) then
       !$omp critical (Merger_Tree_Monotonic_Mass_Growth_Initialize)
       if (.not.monotonicGrowthModuleInitialized) then
          ! Get parameter specifying if monotonic growth should be enforced.
          !@ <inputParameter>
          !@   <name>mergerTreeEnforceMonotonicGrowth</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not enforce monotonic mass growth along the branches of merger trees.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeEnforceMonotonicGrowth',mergerTreeEnforceMonotonicGrowth,defaultValue=.false.)
          ! Flag that module is initialized.
          monotonicGrowthModuleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Monotonic_Mass_Growth_Initialize)
    end if

    ! Enforce monotonic growth if necessary.
    if (mergerTreeEnforceMonotonicGrowth) then
       didModifyTree=.true.
       do while (didModifyTree)
          didModifyTree=.false.
          ! Get root node of the tree.       
          thisNode => thisTree%baseNode
          ! Walk the tree.
          do while (associated(thisNode))
             ! Find nodes that have children.
             if (associated(thisNode%firstChild)) then
                ! Find the mass of all progenitor nodes.
                progenitorMass =  0.0d0
                progenitorNode => thisNode%firstChild
                do while (associated(progenitorNode))
                   progenitorBasicComponent => progenitorNode%basic()
                   progenitorMass =  progenitorMass+progenitorBasicComponent%mass()
                   progenitorNode => progenitorNode%sibling
                end do
                ! Find nodes which are less massive than the sum of their progenitors.
                thisBasicComponent => thisNode%basic()
                if (thisBasicComponent%mass() < progenitorMass) then
                   call thisBasicComponent%massSet(progenitorMass)
                   didModifyTree=.true.
                end if
             end if
             ! Walk to the next node in the tree.             
             call thisNode%walkTree(thisNode)
          end do
        end do
     end if
    return
  end subroutine Merger_Tree_Monotonic_Mass_Growth
  
end module Merger_Trees_Monotonic_Mass_Growth
