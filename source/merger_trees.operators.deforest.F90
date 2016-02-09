!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a deforestation operator on merger trees (i.e. removes all but the most massive tree in a forest).

  !# <mergerTreeOperator name="mergerTreeOperatorDeforest">
  !#  <description>Provides a deforestation operator for merger trees. Given a forest, this operator will destroy all but the first tree in the forest.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorDeforest
     !% A deforestation merger tree operator class.
     private
   contains
     final     ::            deforestDestructor
     procedure :: operate => deforestOperate
  end type mergerTreeOperatorDeforest

  interface mergerTreeOperatorDeforest
     !% Constructors for the deforestation merger tree operator class.
     module procedure deforestConstructorParameters
     module procedure deforestConstructorInternal
  end interface mergerTreeOperatorDeforest

contains

  function deforestConstructorParameters(parameters)
    !% Constructor for the deforestation merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorDeforest)                :: deforestConstructorParameters
    type(inputParameters           ), intent(in   ) :: parameters
    
    return
  end function deforestConstructorParameters

  function deforestConstructorInternal()
    !% Internal constructor for the deforestation merger tree operator class.
    implicit none
    type(mergerTreeOperatorDeforest) :: deforestConstructorInternal
    
    return
  end function deforestConstructorInternal

  elemental subroutine deforestDestructor(self)
    !% Destructor for the deforestation merger tree operator function class.
    implicit none
    type(mergerTreeOperatorDeforest), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine deforestDestructor

  subroutine deforestOperate(self,tree)
    !% Perform a deforestation operation on a merger tree.
    implicit none
    class           (mergerTreeOperatorDeforest), intent(inout)         :: self
    type            (mergerTree                ), intent(inout), target :: tree
    type            (treeNode                  ), pointer               :: baseNode       , nodeNext
    type            (mergerTree                ), pointer               :: currentTree
    class           (nodeComponentBasic        ), pointer               :: basic
    double precision                                                    :: massRootMaximum
    integer                                                             :: treeIndex      , massRootMaximumIndex

    ! Iterate over trees to find the most massive.
    currentTree     => tree
    massRootMaximum =  0.0d0
    treeIndex       =  0
    do while (associated(currentTree))
       treeIndex=treeIndex+1
       basic => currentTree%baseNode%basic()
       if (basic%mass() > massRootMaximum) then
          massRootMaximum     =basic    %mass()
          massRootMaximumIndex=treeIndex
       end if
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
     ! Iterate over trees.
    currentTree => tree
    treeIndex   =  0
    do while (associated(currentTree))
       treeIndex=treeIndex+1
       if (treeIndex /= massRootMaximumIndex) then
          ! Get root node of the tree.
          baseNode => currentTree%baseNode%firstChild
          do while (associated(baseNode))
             nodeNext => baseNode%sibling
             call currentTree%destroyBranch(baseNode)
             baseNode => nodeNext
          end do
       end if
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine deforestOperate
  
