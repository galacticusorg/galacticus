!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements a prune-by-lightcone operator on merger trees.

  use Geometry_Lightcones
  
  !# <mergerTreeOperator name="mergerTreeOperatorPruneLightcone" defaultThreadPrivate="yes">
  !#  <description>Provides a pruning-by-lightcone operator on merger trees. Trees which have no nodes which lie within the lightcone are completely pruned away.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneLightcone
     !% A pruning-by-mass merger tree operator class.
     private
     class(geometryLightconeClass), pointer :: geometryLightcone_
   contains
     final     ::            pruneLightconeDestructor
     procedure :: operate => pruneLightconeOperate
  end type mergerTreeOperatorPruneLightcone

  interface mergerTreeOperatorPruneLightcone
     !% Constructors for the pruning-by-lightcone merger tree operator class.
     module procedure pruneLightconeConstructorParameters
     module procedure pruneLightconeConstructorInternal
  end interface mergerTreeOperatorPruneLightcone

contains

  function pruneLightconeConstructorParameters(parameters)
    !% Constructor for the prune-by-lightcone merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorPruneLightcone)                :: pruneLightconeConstructorParameters
    type(inputParameters                 ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <objectBuilder class="geometryLightcone" name="pruneLightconeConstructorParameters%geometryLightcone_" source="parameters"/>
    return
  end function pruneLightconeConstructorParameters

  function pruneLightconeConstructorInternal(geometryLightcone_)
    !% Internal constructor for the prune-by-lightcone merger tree operator class.
    implicit none
    type (mergerTreeOperatorPruneLightcone)                        :: pruneLightconeConstructorInternal
    class(geometryLightconeClass          ), intent(in   ), target :: geometryLightcone_
    !# <constructorAssign variables="*geometryLightcone_"/>

    return
  end function pruneLightconeConstructorInternal

  subroutine pruneLightconeDestructor(self)
    !% Destructor for the lightcone merger tree operator function class.
    implicit none
    type(mergerTreeOperatorPruneLightcone), intent(inout) :: self
    
    !# <objectDestructor name="self%geometryLightcone_"/>
    return
  end subroutine pruneLightconeDestructor

  subroutine pruneLightconeOperate(self,tree)
    !% Perform a prune-by-lightcone operation on a merger tree.
    use Merger_Trees_Pruning_Utilities
    implicit none
    class  (mergerTreeOperatorPruneLightcone), intent(inout)         :: self
    type   (mergerTree                      ), intent(inout), target :: tree
    type   (treeNode                        ), pointer               :: node       , nodeNext
    type   (mergerTree                      ), pointer               :: treeCurrent, treeNext

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Get root node of the tree.
       node  => treeCurrent%baseNode
       do while (associated(node))
          if (self%geometryLightcone_%isInLightcone(node,atPresentEpoch=.false.)) return
          node => node%walkTree()
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Entire forest is outside lightcone. Destroy all but the base node in the first tree. (Leaving just the base node makes the
    ! tree inert - i.e. it can not do anything.) Destroy any additional trees in the forest.
    node => tree%baseNode%firstChild
    do while (associated(node))
       nodeNext => node%sibling
       call Merger_Tree_Prune_Clean_Branch(node)
       call node%destroyBranch()
       deallocate(node)
       node => nodeNext
    end do
    treeCurrent          => tree%nextTree
    tree       %nextTree => null()
    do while (associated(treeCurrent))
       treeNext => treeCurrent%nextTree
       ! Destroy the tree.
       call treeCurrent%destroy()
       deallocate(treeCurrent)
       ! Move to the next tree.
       treeCurrent => treeNext
    end do
    return
  end subroutine pruneLightconeOperate
