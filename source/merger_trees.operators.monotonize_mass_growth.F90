!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Contains a module which implements a merger tree operator which makes mass growth along
  !% branch monotonically increasing.
  
  !# <mergerTreeOperator name="mergerTreeOperatorMonotonizeMassGrowth">
  !#  <description>
  !#   A merger tree operator which makes mass growth along merger tree branches monotonic.
  !# </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorMonotonizeMassGrowth
     !% A merger tree operator class makes mass growth along branch monotonically increasing.
     private
   contains
     final     ::             monotonizeMassGrowthDestructor
     procedure :: operate  => monotonizeMassGrowthOperate
  end type mergerTreeOperatorMonotonizeMassGrowth
  
  interface mergerTreeOperatorMonotonizeMassGrowth
     !% Constructors for the mass growth monotonizing merger tree operator class.
     module procedure monotonizeMassGrowthConstructorParameters
     module procedure monotonizeMassGrowthConstructorInternal
  end interface mergerTreeOperatorMonotonizeMassGrowth

contains

  function monotonizeMassGrowthConstructorParameters(parameters)
    !% Constructor for the mass growth monotonizing merger tree operator class which takes a
    !% parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorMonotonizeMassGrowth)                :: monotonizeMassGrowthConstructorParameters
    type(inputParameters                       ), intent(in   ) :: parameters
        
    return
  end function monotonizeMassGrowthConstructorParameters

  function monotonizeMassGrowthConstructorInternal()
    !% Internal constructor for the mass growth monotonizing merger tree operator class.
    implicit none
    type(mergerTreeOperatorMonotonizeMassGrowth) :: monotonizeMassGrowthConstructorInternal

    return
  end function monotonizeMassGrowthConstructorInternal

  elemental subroutine monotonizeMassGrowthDestructor(self)
    !% Destructor for the mass growth monotonizing merger tree operator function class.
    implicit none
    type(mergerTreeOperatorMonotonizeMassGrowth), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine monotonizeMassGrowthDestructor

  subroutine monotonizeMassGrowthOperate(self,tree)
    !% Perform a mass growth monotonizing operation on a merger tree.
    implicit none
    class           (mergerTreeOperatorMonotonizeMassGrowth), intent(inout)         :: self
    type            (mergerTree                            ), intent(inout), target :: tree
    type            (treeNode                              ), pointer               :: nodeProgenitor          , node
    class           (nodeComponentBasic                    ), pointer               :: basicProgenitor, basic
    type            (mergerTree                            ), pointer               :: treeCurrent
    logical                                                                         :: didModifyTree
    double precision                                                                :: massProgenitor

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       didModifyTree=.true.
       do while (didModifyTree)
          didModifyTree=.false.
          ! Get root node of the tree.
          node => treeCurrent%baseNode
          ! Walk the tree.
          do while (associated(node))
             ! Find nodes that have children.
             if (associated(node%firstChild)) then
                ! Find the mass of all progenitor nodes.
                massProgenitor =  0.0d0
                nodeProgenitor => node%firstChild
                do while (associated(nodeProgenitor))
                   basicProgenitor =>  nodeProgenitor %basic  ()
                   massProgenitor  =  +massProgenitor            &
                        &             +basicProgenitor%mass   ()
                   nodeProgenitor  =>  nodeProgenitor %sibling
                end do
                ! Find nodes which are less massive than the sum of their progenitors.
                basic => node%basic()
                if (basic%mass() < massProgenitor) then
                   call basic%massSet(massProgenitor)
                   didModifyTree=.true.
                end if
             end if
             ! Walk to the next node in the tree.
             node => node%walkTree()
          end do
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine monotonizeMassGrowthOperate
  
