!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module which implements \href{https://en.wikipedia.org/wiki/Binary_search_tree}{binary search tree} data structures.
!!}

module Binary_Search_Trees
  !!{
  Implements \href{https://en.wikipedia.org/wiki/Binary_search_tree}{binary search tree} data structures.
  !!}
  implicit none
  private
  public :: binaryTree, binaryTreeNode

  type :: binaryTree
     !!{
     A type for implementing \href{https://en.wikipedia.org/wiki/Binary_search_tree}{binary search tree} data structures.
     !!}
     type(binaryTreeNode), pointer :: root => null()
   contains
     !![
     <methods>
       <method method ="insert"  description="Insert a new node into the binary search tree."/>
       <method method ="bracket" description="Bracket a value in the binary search tree."    />
     </methods>
     !!]
     final     ::            binaryTreeDestructor
     procedure :: insert  => binaryTreeInsert
     procedure :: bracket => binaryTreeBracket
  end type binaryTree
  
  type :: binaryTreeNode
     !!{
     A type for implementing \href{https://en.wikipedia.org/wiki/Binary_search_tree}{binary search tree} data structures.
     !!}
     type            (binaryTreeNode), pointer :: left   => null(), right => null(), &
          &                                       parent => null()
     double precision                          :: key             , value
   contains
     !![
     <methods>
       <method method ="successor"   description="Return the successor node in the tree."                    />
       <method method ="predecessor" description="Return the predescessor node in the tree."                 />
       <method method ="treeMinimum" description="Return the minimum node in the sub-tree of the given node."/>
       <method method ="treeMaximum" description="Return the maximum node in the sub-tree of the given node."/>
     </methods>
     !!]
     final     ::                binaryTreeNodeDestructor
     procedure :: successor   => binaryTreeNodeSuccessor
     procedure :: predecessor => binaryTreeNodePredecessor
     procedure :: treeMinimum => binaryTreeNodeSubTreeMinimum
     procedure :: treeMaximum => binaryTreeNodeSubTreeMaximum
  end type binaryTreeNode

  interface binaryTreeNode
     !!{
     Constructor for binary search tree nodes.
     !!}
     module procedure :: binaryTreeNodeConstructor
  end interface binaryTreeNode
  

contains

  subroutine binaryTreeDestructor(self)
    !!{
    Destructor for binary search trees.
    !!}
    implicit none
    type(binaryTree), intent(inout) :: self

    if (associated(self%root)) deallocate(self%root)
    return
  end subroutine binaryTreeDestructor

  subroutine binaryTreeInsert(self,key,value)
    !!{
    Insert a new node to a binary search tree.
    !!}
    implicit none
    class           (binaryTree    ), intent(inout)          :: self
    double precision                , intent(in   )          :: key    , value
    type            (binaryTreeNode)               , pointer :: current, previous, &
         &                                                      new

    new      => binaryTreeNode(key,value)
    previous => null()
    current  => self%root
    do while (associated(current))
       previous => current
       if (new%key < current%key) then
          current => current%left
       else
          current => current%right
       end if
    end do
    new%parent => previous
    if (.not.associated(previous)) then
       self    %root  => new
    else if (new%key < previous%key) then
       previous%left  => new
    else
       previous%right => new
    end if
    return
  end subroutine binaryTreeInsert
  
  subroutine binaryTreeBracket(self,key,left,right)
    !!{
    Find nodes which bracket a key.
    !!}
    implicit none
    class           (binaryTree    ), intent(in   )          :: self
    double precision                , intent(in   )          :: key
    type            (binaryTreeNode), intent(  out), pointer :: left   , right
    type            (binaryTreeNode)               , pointer :: current

    left    => null()
    right   => null()
    current => self%root
    do while (associated(current))
       if (key == current%key) then
          left    => current
          right   => current
          exit
       else if (key < current%key) then
          right   => current
          current => current%left
       else
          left    => current
          current => current%right
       end if
    end do
    return
  end subroutine binaryTreeBracket
  
  function binaryTreeNodeConstructor(key,value) result(self)
    !!{
    Constructor for binary search tree nodes.
    !!}
    implicit none
    type            (binaryTreeNode), pointer       :: self
    double precision                , intent(in   ) :: key, value

    allocate(self)
    self%key   =  key
    self%value =  value
    self%left  => null()
    self%right => null()
    return
  end function binaryTreeNodeConstructor

  subroutine binaryTreeNodeDestructor(self)
    !!{
    Destructor for binary search trees.
    !!}
    implicit none
    type(binaryTreeNode), intent(inout) :: self

    if (associated(self%left )) deallocate(self%left )
    if (associated(self%right)) deallocate(self%right)
    return
  end subroutine binaryTreeNodeDestructor

  function binaryTreeNodeSuccessor(self) result(successor)
    !!{
    Return a pointer to the successor node of the given node.
    !!}
    implicit none
    class(binaryTreeNode), intent(inout), target  :: self
    type (binaryTreeNode)               , pointer :: successor, current
    
    if (associated(self%right)) then
       successor => self%right%treeMinimum()
    else
       current   => self
       successor => self%parent
       do while (associated(successor) .and. associated(current,successor%right))
          current   => successor
          successor => successor%parent
       end do
    end if
    return
  end function binaryTreeNodeSuccessor
  
  function binaryTreeNodePredecessor(self) result(predecessor)
    !!{
    Return a pointer to the predecessor node of the given node.
    !!}
    implicit none
    class(binaryTreeNode), intent(inout), target  :: self
    type (binaryTreeNode)               , pointer :: predecessor, current
    
    if (associated(self%left)) then
       predecessor => self%left%treeMaximum()
    else
       current     => self
       predecessor => self%parent
       do while (associated(predecessor) .and. associated(current,predecessor%left))
          current     => predecessor
          predecessor => predecessor%parent
       end do
    end if
    return
  end function binaryTreeNodePredecessor
  
  function binaryTreeNodeSubTreeMinimum(self) result(treeMinimum)
    !!{
    Return a pointer to the minimum node in the sub-tree of the given node.
    !!}
    implicit none
    class(binaryTreeNode), intent(inout), target  :: self
    type (binaryTreeNode)               , pointer :: treeMinimum
    
    treeMinimum => self
    do while (associated(treeMinimum%left))
       treeMinimum => treeMinimum%left
    end do
    return
  end function binaryTreeNodeSubTreeMinimum
  
  function binaryTreeNodeSubTreeMaximum(self) result(treeMaximum)
    !!{
    Return a pointer to the maximum node in the sub-tree of the given node.
    !!}
    implicit none
    class(binaryTreeNode), intent(inout), target  :: self
    type (binaryTreeNode)               , pointer :: treeMaximum
    
    treeMaximum => self
    do while (associated(treeMaximum%right))
       treeMaximum => treeMaximum%right
    end do
    return
  end function binaryTreeNodeSubTreeMaximum
  
end module Binary_Search_Trees
