!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a pruning-by-mass operator on merger trees.

  !# <mergerTreeOperator name="mergerTreeOperatorPruneByMass">
  !#  <description>Provides a pruning-by-mass operator on merger trees.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneByMass
     !% A pruning-by-mass merger tree operator class.
     private
     double precision :: massThreshold
     logical          :: preservePrimaryProgenitor
   contains
     final     ::            pruneByMassDestructor
     procedure :: operate => pruneByMassOperate
  end type mergerTreeOperatorPruneByMass

  interface mergerTreeOperatorPruneByMass
     !% Constructors for the pruning-by-mass merger tree operator class.
     module procedure pruneByMassConstructorParameters
     module procedure pruneByMassConstructorInternal
  end interface mergerTreeOperatorPruneByMass

contains

  function pruneByMassConstructorParameters(parameters)
    !% Constructor for the prune-by-mass merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorPruneByMass)                :: pruneByMassConstructorParameters
    type(inputParameters              ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>massThreshold</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <variable>pruneByMassConstructorParameters%massThreshold</variable>
    !#   <description>Threshold mass below which merger tree branches should be pruned.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>preservePrimaryProgenitor</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <variable>pruneByMassConstructorParameters%preservePrimaryProgenitor</variable>
    !#   <description>If true, primary progenitor status is preserved even if the primary progenitor is pruned from the tree.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    return
  end function pruneByMassConstructorParameters

  function pruneByMassConstructorInternal(massThreshold,preservePrimaryProgenitor)
    !% Internal constructor for the prune-by-mass merger tree operator class.
    implicit none
    type            (mergerTreeOperatorPruneByMass)                :: pruneByMassConstructorInternal
    double precision                               , intent(in   ) :: massThreshold
    logical                                        , intent(in   ) :: preservePrimaryProgenitor
    
    pruneByMassConstructorInternal%massThreshold            =massThreshold
    pruneByMassConstructorInternal%preservePrimaryProgenitor=preservePrimaryProgenitor
    return
  end function pruneByMassConstructorInternal

  elemental subroutine pruneByMassDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorPruneByMass), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    ! Nothing to do.
    return
  end subroutine pruneByMassDestructor

  subroutine pruneByMassOperate(self,tree)
    !% Perform a prune-by-mass operation on a merger tree.
    use Merger_Trees_Pruning_Utilities
    implicit none
    class  (mergerTreeOperatorPruneByMass), intent(inout)         :: self
    type   (mergerTree                   ), intent(inout), target :: tree
    type   (treeNode                     ), pointer               :: nodeNext     , nodePrevious, &
         &                                                           node         , nodeWork
    class  (nodeComponentBasic           ), pointer               :: basic        , basicPrevious
    type   (mergerTree                   ), pointer               :: currentTree
    logical                                                       :: didPruning

    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       didPruning=.true.
       do while (didPruning)
          didPruning=.false.
          ! Get root node of the tree.
          node  => currentTree%baseNode
          basic => node       %basic   ()
          if (basic%mass() < self%massThreshold) then
             ! Entire tree is below threshold. Destroy all but this base node. (Leaving just
             ! the base node makes the tree inert - i.e. it can not do anything.)
             call Merger_Tree_Prune_Clean_Branch (node)
             nodeWork => node%firstChild
             do while (associated(nodeWork))
                nodeNext => nodeWork%sibling
                call Merger_Tree_Prune_Clean_Branch(nodeWork)
                call nodeWork%destroyBranch()
                deallocate(nodeWork)
                nodeWork => nodeNext
             end do
             nullify(node%firstChild)
             nodeWork => node%firstSatellite
             do while (associated(nodeWork))
                nodeNext => nodeWork%sibling
                call Merger_Tree_Prune_Clean_Branch(nodeWork)
                call nodeWork%destroyBranch()
                deallocate(nodeWork)
                nodeWork => nodeNext
             end do
             nullify(node%firstSatellite)
          else
             ! Walk the tree, pruning branches.
             do while (associated(node))
                basic => node%basic()
                ! Record the parent node to which we will return.
                nodePrevious => node%parent
                if (basic%mass() < self%massThreshold) then
                   didPruning=.true.
                   ! Decouple from other nodes.
                   basicPrevious => nodePrevious%basic()
                   call Merger_Tree_Prune_Unlink_Parent(node,nodePrevious,basicPrevious%mass() < self%massThreshold,self%preservePrimaryProgenitor)
                   ! Clean the branch.
                   call Merger_Tree_Prune_Clean_Branch (node)
                   ! Destroy the branch.
                   call node%destroyBranch()
                   deallocate(node)
                   ! Return to parent node.
                   node => nodePrevious
                end if
                node => node%walkTree()
             end do
          end if
       end do
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    ! Uniqueify nodes.
    call Merger_Tree_Prune_Uniqueify_IDs(tree)
    return
  end subroutine pruneByMassOperate
